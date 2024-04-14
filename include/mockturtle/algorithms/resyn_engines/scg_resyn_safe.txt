/* mockturtle: C++ logic network library
 * Copyscght (C) 2018-2022  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the scghts to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyscght notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file scg_resyn.hpp
  \brief Resynthesis by extraction of functional cuts

  \author Andrea Costamagna
*/

#pragma once

#include "../../utils/index_list.hpp"
#include "../../utils/node_map.hpp"
#include "../../utils/spfd_utils.hpp"
#include "../../utils/stopwatch.hpp"
#include "../../utils/tech_library.hpp"
#include <mockturtle/algorithms/emap2.hpp>

#include <fmt/format.h>
#include <kitty/kitty.hpp>

#include <algorithm>
#include <mutex>
#include <optional>
#include <thread>
#include <type_traits>
#include <vector>
#include <fstream>


namespace mockturtle
{

namespace scopt
{

bool VERBOSE{ false };

enum support_selection_t
{
  GREEDY,
  NGREEDY,
  PIVOT,
};

std::mt19937 RIGRNG( 5 );

struct scg_resyn_static_params
{
  using base_type = scg_resyn_static_params;

  /*! \brief Whether to copy truth tables. */
  static constexpr bool copy_tts{ false };

  /*! \brief Reserved capacity for divisor truth tables (number of divisors). */
  static constexpr uint32_t reserve{ 200u };

  /*! \brief Whether to preserve depth. */
  static constexpr bool preserve_depth{ false };

  /*! \brief Whether the divisors have uniform costs (size and depth, whenever relevant). */
  static constexpr bool uniform_div_cost{ true };

  static constexpr uint32_t max_support_size{ 6u };
  static constexpr uint32_t fraction_of_10{ 10 };

  static constexpr int max_fanin_size = -1;
  static constexpr bool accept_worse{ false };
  static constexpr bool on_the_fly{ false };

  static constexpr uint32_t nBest{ 2 };

  static constexpr support_selection_t support_selection{ GREEDY };

  using truth_table_storage_type = void;
  using node_type = void;
};

template<class TT>
struct scg_resyn_static_params_default : public scg_resyn_static_params
{
  using truth_table_storage_type = std::vector<TT>;
  using node_type = uint32_t;
};

template<class Ntk, support_selection_t SUP_SEL, uint32_t SUPP_SIZE, int K = -1, int NRELAX = 0>
struct scg_resyn_static_params_for_sim_resub : public scg_resyn_static_params
{
  using truth_table_storage_type = incomplete_node_map<kitty::partial_truth_table, Ntk>;
  using node_type = typename Ntk::node;
  static constexpr support_selection_t support_selection = SUP_SEL;
  static constexpr uint32_t max_support_size = SUPP_SIZE;
  static constexpr int max_fanin_size = K;
  static constexpr bool accept_worse = NRELAX > 0;
};

template<class Ntk, support_selection_t SUP_SEL, uint32_t NumVars, uint32_t SUPP_SIZE, int K = -1, int NRELAX = 0>
struct scg_resyn_static_params_for_sim_resub_static : public scg_resyn_static_params
{
  using truth_table_storage_type = incomplete_node_map<kitty::static_truth_table<NumVars>, Ntk>;
  using node_type = typename Ntk::node;
  static constexpr support_selection_t support_selection = SUP_SEL;
  static constexpr uint32_t max_support_size = SUPP_SIZE;
  static constexpr int max_fanin_size = K;
  static constexpr bool accept_worse = NRELAX > 0;
};

struct scg_resyn_stats
{
  /*! \brief Time for finding 0-resub and collecting unate literals. */
  stopwatch<>::duration time_0resub{ 0 };

  /*! \brief Time for finding 0-resub and collecting unate literals. */
  stopwatch<>::duration time_supp{ 0 };

  /*! \brief Time for finding resub. */
  stopwatch<>::duration time_resub{ 0 };

  /*! \brief Time for sorting unate literals and unate pairs. */
  stopwatch<>::duration time_sort{ 0 };

  /*! \brief Time for collecting unate pairs. */
  stopwatch<>::duration time_collect_pairs{ 0 };

  /*! \brief Time for dividing the target and recursive call. */
  stopwatch<>::duration time_divide{ 0 };

  void report() const
  {
    fmt::print( "[i]         <xag_resyn_decompose>\n" );
    fmt::print( "[i]             0-resub      : {:>5.2f} secs\n", to_seconds( time_0resub ) );
    fmt::print( "[i]             k-resub      : {:>5.2f} secs\n", to_seconds( time_resub ) );
    fmt::print( "[i]             sort         : {:>5.2f} secs\n", to_seconds( time_sort ) );
    fmt::print( "[i]             collect pairs: {:>5.2f} secs\n", to_seconds( time_collect_pairs ) );
    fmt::print( "[i]             dividing     : {:>5.2f} secs\n", to_seconds( time_divide ) );
  }
};

/*! \brief Logic resynthesis engine for AIGs or XAGs.
 *
 * The algorithm is based on ABC's implementation in `giaResub.c` by Alan Mishchenko.
 *
 * Divisors are classified as positive unate (not overlapping with target offset),
 * negative unate (not overlapping with target onset), or binate (overlapping with
 * both onset and offset). Furthermore, pairs of binate divisors are combined with
 * an AND operation and considering all possible input polarities and again classified
 * as positive unate, negative unate or binate. Simple solutions of zero cost
 * (one unate divisor), one node (two unate divisors), two nodes (one unate divisor +
 * one unate pair), and three nodes (two unate pairs) are exhaustively examined.
 * When no simple solutions can be found, the algorithm heuristically chooses an unate
 * divisor or an unate pair to divide the target function with and recursively calls
 * itself to decompose the remainder function.
   \verbatim embed:rst

   Example

   .. code-block:: c++

      using TT = kitty::static_truth_table<6>;
      const std::vector<aig_network::node> divisors = ...;
      const node_map<TT, aig_network> tts = ...;
      const TT target = ..., care = ...;
      xag_resyn_stats st;
      xag_resyn_decompose<TT, node_map<TT, aig_network>, false, false, aig_network::node> resyn( st );
      auto result = resyn( target, care, divisors.begin(), divisors.end(), tts );
   \endverbatim
 */

template<class TT, class static_params = scg_resyn_static_params_default<TT>, support_selection_t SUP_SEL = GREEDY>
class scg_resyn_decompose
{

public:
  using stats = scg_resyn_stats;
  using index_list_t = large_lig_index_list;
  using truth_table_t = TT;
  using truth_tableK_t = kitty::static_truth_table<static_params::max_support_size>;

private:
  struct scored_div
  {
    scored_div( uint32_t l, uint32_t s )
        : div( l ), score( s )
    {}

    bool operator==( scored_div const& other ) const
    {
      return div == other.lit;
    }

    bool operator<( scored_div const& other ) const
    {
      return score < other.score;
    }

    bool operator>( scored_div const& other ) const
    {
      return score > other.score;
    }

    uint32_t div;
    uint32_t score;
  };

  struct fscored_div
  {
    fscored_div( uint32_t l, double s )
        : div( l ), score( s )
    {}

    bool operator==( fscored_div const& other ) const
    {
      return div == other.lit;
    }

    bool operator<( fscored_div const& other ) const
    {
      return score < other.score;
    }

    bool operator>( fscored_div const& other ) const
    {
      return score > other.score;
    }

    uint32_t div;
    double score;
  };

public:
  explicit scg_resyn_decompose( std::vector<gate> const& gates, stats& st ) noexcept
      : st( st ), _lps(false),  _gates(gates), _tech_lib( gates, _tps ), _lib( _resyn, _lps )
  {
    exact_library_params lps;
    lps.compute_dc_classes = true;

    static_assert( std::is_same_v<typename static_params::base_type, scg_resyn_static_params>, "Invalid static_params type" );
    static_assert( !( static_params::uniform_div_cost && static_params::preserve_depth ), "If depth is to be preserved, divisor depth cost must be provided (usually not uniform)" );
    divisors.reserve( static_params::reserve );
    RIGRNG.seed( 5 );

    auto [buf_area, buf_delay, buf_id] = _tech_lib.get_buffer_info();
    auto [inv_area, inv_delay, inv_id] = _tech_lib.get_inverter_info();

    _area_th = std::min( buf_area, inv_area );


    kitty::static_truth_table<4u> tt;
    int i{0};
    std::string line;
    std::ifstream fTts ("sky130.tts");
    if (fTts.is_open())
    {
      while ( std::getline (fTts,line) )
      {
        kitty::create_from_binary_string( tt, line );
        _pClassMap[tt._bits]=i++;
      }
      fTts.close();
    }

    std::ifstream fAreas ("sky130.area");
    if (fAreas.is_open())
    {
      while ( std::getline (fAreas,line) )
      {
        _areas.push_back( std::stof( line ) );
      }
      fAreas.close();
    }


    std::ifstream fLists ("sky130.list");
    if (fLists.is_open())
    {
      while ( std::getline (fLists,line) )
      {
        std::vector<uint32_t> list;
        uint32_t number;
        std::istringstream line_stream(line);
        while (line_stream >> number)
        {
            list.push_back(number);
        }
        _idlists.push_back( list );
      }
      fLists.close();
    }

  }

  /*! \brief Perform XAG resynthesis.
   *
   * `tts[*begin]` must be of type `TT`.
   * Moreover, if `static_params::copy_tts = false`, `*begin` must be of type `static_params::node_type`.
   *
   * \param target Truth table of the target function.
   * \param care Truth table of the care set.
   * \param begin Begin iterator to divisor nodes.
   * \param end End iterator to divisor nodes.
   * \param tts A data structure (e.g. std::vector<TT>) that stores the truth tables of the divisor functions.
   * \param max_size Maximum number of nodes allowed in the dependency circuit.
   */
  template<class iterator_type,
           bool enabled = static_params::uniform_div_cost && !static_params::preserve_depth, typename = std::enable_if_t<enabled>>
  std::optional<index_list_t> operator()( TT const& target, TT const& care, iterator_type begin, iterator_type end, typename static_params::truth_table_storage_type const& tts, double max_size )
  {
    static_assert( static_params::copy_tts || std::is_same_v<typename std::iterator_traits<iterator_type>::value_type, typename static_params::node_type>, "iterator_type does not dereference to static_params::node_type" );

    ptts = &tts;
    on_off_sets[0] = ~target & care;
    on_off_sets[1] = target & care;

    _uSPFD.init( target, care );

    divisors.clear(); /* clear previous data and reserve 1 dummy node for constant */
    scored_divs.clear();

    while ( begin != end )
    {
      if constexpr ( static_params::copy_tts )
      {
        divisors.emplace_back( ( *ptts )[*begin] );
        scored_divs.emplace_back( divisors.size(), _uSPFD.evaluate( get_div( divisors.size() - 1 ) ) );
      }
      else
      {
        divisors.emplace_back( *begin );
        scored_divs.emplace_back( divisors.size(), _uSPFD.evaluate( get_div( divisors.size() - 1 ) ) );
      }
      ++begin;
    }

    call_with_stopwatch( st.time_sort, [&]() {
      std::sort( scored_divs.begin(), scored_divs.end() );
    } );

    return compute_function( max_size );
  }

  template<class iterator_type,
           bool enabled = static_params::uniform_div_cost && !static_params::preserve_depth, typename = std::enable_if_t<enabled>>
  std::optional<index_list_t> operator()( TT const& target, TT const& care, iterator_type begin, iterator_type end, std::vector<double> const& idelays, typename static_params::truth_table_storage_type const& tts, double max_size )
  {
    static_assert( static_params::copy_tts || std::is_same_v<typename std::iterator_traits<iterator_type>::value_type, typename static_params::node_type>, "iterator_type does not dereference to static_params::node_type" );

    ptts = &tts;
    on_off_sets[0] = ~target & care;
    on_off_sets[1] = target & care;

    _uSPFD.init( target, care );

    divisors.clear();//resize( 1 ); /* clear previous data and reserve 1 dummy node for constant */
    scored_divs.clear();

    while ( begin != end )
    {
      if constexpr ( static_params::copy_tts )
      {
        divisors.emplace_back( ( *ptts )[*begin] );
        scored_divs.emplace_back( divisors.size(), _uSPFD.evaluate( get_div( divisors.size() - 1 ) ) );
      }
      else
      {
        divisors.emplace_back( *begin );
        scored_divs.emplace_back( divisors.size(), _uSPFD.evaluate( get_div( divisors.size() - 1 ) ) );
      }
      ++begin;
    }

    call_with_stopwatch( st.time_sort, [&]() {
      std::sort( scored_divs.begin(), scored_divs.end() );
    } );
    _idelays=idelays;
    return compute_function( max_size );
  }

  template<class iterator_type, class Fn,
           bool enabled = !static_params::uniform_div_cost && !static_params::preserve_depth, typename = std::enable_if_t<enabled>>
  std::optional<index_list_t> operator()( TT const& target, TT const& care, iterator_type begin, iterator_type end, typename static_params::truth_table_storage_type const& tts, Fn&& size_cost, uint32_t max_size = std::numeric_limits<uint32_t>::max() )
  {}

  template<class iterator_type, class Fn,
           bool enabled = !static_params::uniform_div_cost && static_params::preserve_depth, typename = std::enable_if_t<enabled>>
  std::optional<index_list_t> operator()( TT const& target, TT const& care, iterator_type begin, iterator_type end, typename static_params::truth_table_storage_type const& tts, Fn&& size_cost, Fn&& depth_cost, uint32_t max_size = std::numeric_limits<uint32_t>::max(), uint32_t max_depth = std::numeric_limits<uint32_t>::max() )
  {}

private:
  std::optional<index_list_t> compute_function( double num_inserts )
  {
    index_list.clear();
    index_list.reset_area();
    index_list.add_inputs( divisors.size() );
    //std::cout << to_index_list_string(index_list) << std::endl;

    auto const lit = compute_function_rec( num_inserts );
    if ( lit )
    {
      index_list.add_output( *lit );
      return index_list;
    }
    return std::nullopt;
  }

  std::optional<uint32_t> compute_function_rec( double num_inserts )
  {

    /* try 0-resub and collect unate literals */
    auto const res0 = call_with_stopwatch( st.time_0resub, [&]() {
      return try_0resub( num_inserts );
    } );
    if ( res0 )
    {
      return *res0;
    }

    if ( num_inserts <= _area_th )
    {
      return std::nullopt;
    }

    auto const supp = call_with_stopwatch( st.time_supp, [&]() {
      return find_support();
    } );

    /* try n-resub */
    auto const resn = call_with_stopwatch( st.time_resub, [&]() {
      return try_nresub( num_inserts );
    } );
    if ( resn )
    {
      return *resn;
    }

    return std::nullopt;
  }

  /* See if there is a constant or divisor covering all on-set bits or all off-set bits.
     1. Check constant-resub
     2. Collect unate literals
     3. Find 0-resub (both positive unate and negative unate) and collect binate (neither pos nor neg unate) divisors
   */
  std::optional<uint32_t> try_0resub( double max_area )
  {
    auto [buf_area, buf_delay, buf_id] = _tech_lib.get_buffer_info();
    auto [inv_area, inv_delay, inv_id] = _tech_lib.get_inverter_info();

    num_bits[0] = kitty::count_ones( on_off_sets[0] ); /* off-set */
    num_bits[1] = kitty::count_ones( on_off_sets[1] ); /* on-set */
    if ( num_bits[0] == 0 )
    {
      return 1;
    }
    if ( num_bits[1] == 0 )
    {
      return 0;
    }

    for ( auto v = 0u; v < divisors.size(); ++v )
    {
      bool unateness[4] = { false, false, false, false };
      /* check intersection with off-set */
      if ( kitty::intersection_is_empty<TT, 1, 1>( get_div( v ), on_off_sets[0] ) )
      {
        unateness[0] = true;
      }
      else if ( kitty::intersection_is_empty<TT, 0, 1>( get_div( v ), on_off_sets[0] ) )
      {
        unateness[1] = true;
      }

      /* check intersection with on-set */
      if ( kitty::intersection_is_empty<TT, 1, 1>( get_div( v ), on_off_sets[1] ) )
      {
        unateness[2] = true;
      }
      else if ( kitty::intersection_is_empty<TT, 0, 1>( get_div( v ), on_off_sets[1] ) )
      {
        unateness[3] = true;
      }

      /* 0-resub */
      if ( unateness[0] && unateness[3] )
      {
        //if( ( buf_area < max_area ) )
        {
        //  assert( v < index_list.num_pis() );
        //printf("lit out %d\n", (v+1) << 1);
          return (v+1) << 1;
        //  return index_list.add_function( std::vector{ v << 1 }, _gates[buf_id].function, buf_area, buf_id );
        }
      }
      if ( unateness[1] && unateness[2] )
      {
        if( (inv_area < max_area) )
        {
          //assert( v < index_list.num_pis() );
          //printf("lit out not(%d)\n", (v+1) << 1);

          return index_list.add_function( std::vector{ (v+1) << 1 }, _gates[inv_id].function, inv_area, inv_id );
        }
      }
    }
    return std::nullopt;
  }

  /* See if we cna define a new function of the other divisors
   */
  std::optional<uint32_t> try_nresub( double max_inserts )
  {
    auto supp = find_support();

    if ( supp )
    {
      if( ((*supp).size()>4) || static_params::on_the_fly )
      {
        auto const [func, care] = extract_functionality_from_signatures( *supp );
        return _map_on_the_fly( *supp, func, care, max_inserts );
      }
      else
      {
        auto const [func, care] = extract_functionality_from_signatures4( *supp );
        return _map_with_database( *supp, func, care, max_inserts );
      }
    }
    /* resynthesis */
    return std::nullopt;
  }

  std::tuple<kitty::dynamic_truth_table, kitty::dynamic_truth_table> extract_functionality_from_signatures( std::vector<uint32_t> const& supp )
  {
    assert( supp.size() <= static_params::max_support_size );

    std::vector<kitty::dynamic_truth_table> xs;
    for ( uint32_t i{ 0 }; i < supp.size(); ++i )
    {
      xs.emplace_back( supp.size() );
      kitty::create_nth_var( xs[i], i );
    }

    kitty::dynamic_truth_table func_s( supp.size() );
    kitty::dynamic_truth_table care_s = func_s.construct();
    auto temp = _uSPFD.care.construct();
    auto temp_s = func_s.construct();

    for ( uint32_t m{ 0u }; m < ( 1u << supp.size() ); ++m )
    {
      temp = temp | ~temp;
      temp_s = temp_s | ~temp_s;

      for ( uint32_t l{ 0u }; l < supp.size(); ++l )
      {
        if ( ( m >> l ) & 0x1 == 0x1 )
        {
          temp &= get_div( supp[l] );
          temp_s &= xs[l];
        }
        else
        {
          temp &= ~get_div( supp[l] );
          temp_s &= ~xs[l];
        }
      }

      if ( kitty::count_ones( temp & _uSPFD.care ) > 0 ) // care value
      {
        care_s |= temp_s;

        if ( kitty::count_ones( temp & _uSPFD.func[1] ) > 0 )
        {
          func_s |= temp_s;
        }
      }
    }
    auto rnd_tt = func_s.construct();
    kitty::create_random( rnd_tt, _seed++ );

    func_s |= ( rnd_tt & ~care_s );
    // kitty::print_binary( func_s );printf("\n");

    return std::make_tuple( func_s, care_s );
  }

  std::tuple<kitty::static_truth_table<4u>, kitty::static_truth_table<4u>> extract_functionality_from_signatures4( std::vector<uint32_t> const& supp )
  {
    assert( supp.size() <= 4u );

    std::vector<kitty::static_truth_table<4u>> xs;
    for ( uint32_t i{ 0 }; i < supp.size(); ++i )
    {
      xs.emplace_back();
      kitty::create_nth_var( xs[i], i );
    }

    kitty::static_truth_table<4u> func_s;
    kitty::static_truth_table<4u> care_s = func_s.construct();
    kitty::static_truth_table<4u> temp_s = func_s.construct();
    auto temp = _uSPFD.care.construct();

    for ( uint32_t m{ 0u }; m < ( 1u << supp.size() ); ++m )
    {
      temp = temp | ~temp;
      temp_s = temp_s | ~temp_s;

      for ( uint32_t l{ 0u }; l < supp.size(); ++l )
      {
        if ( ( m >> l ) & 0x1 == 0x1 )
        {
          temp &= get_div( supp[l] );
          temp_s &= xs[l];
        }
        else
        {
          temp &= ~get_div( supp[l] );
          temp_s &= ~xs[l];
        }
      }

      if ( kitty::count_ones( temp & _uSPFD.care ) > 0 ) // care value
      {
        care_s |= temp_s;

        if ( kitty::count_ones( temp & _uSPFD.func[1] ) > 0 )
        {
          func_s |= temp_s;
        }
      }
    }
    auto rnd_tt = func_s.construct();
    kitty::create_random( rnd_tt, _seed++ );

    func_s |= ( rnd_tt & ~care_s );
    // kitty::print_binary( func_s );printf("\n");

    return std::make_tuple( func_s, care_s );
  }


#pragma region synthesis

#pragma region synthesize aig
  aig_network::signal synthesize_aig_inplace( aig_network& aig, std::vector<signal<aig_network>>& pis, kitty::dynamic_truth_table tt, kitty::dynamic_truth_table mk )
  {
    auto fout = synthesize_aig_rec( aig, pis, tt, mk );
    aig.create_po( fout );
    return fout;
  }

  aig_network::signal synthesize_aig_rec( aig_network& aig, std::vector<signal<aig_network>> pis, kitty::dynamic_truth_table const& tt, kitty::dynamic_truth_table const& mk )
  {
    if ( kitty::is_const0( ( tt & mk ) ) )
      return aig.get_constant(false);
    if ( kitty::equal( tt & mk, mk ) )
      return aig.get_constant(true);
    if ( pis.size() == 1u )
      return kitty::is_normal( tt ) ? pis[0] : !pis[0];

    uint32_t idx = pis.size() - 1;
    uint32_t impurity;
    uint32_t best_impurity=std::numeric_limits<uint32_t>::max();
    for( int i{0}; i<pis.size(); ++i )
    {
      auto tt0 = kitty::cofactor0( tt, aig.pi_index(aig.get_node(pis[i]) ));
      auto tt1 = kitty::cofactor1( tt, aig.pi_index(aig.get_node(pis[i]) ));
      auto mk0 = kitty::cofactor0( mk, aig.pi_index(aig.get_node(pis[i]) ));
      auto mk1 = kitty::cofactor1( mk, aig.pi_index(aig.get_node(pis[i]) ));

      if( kitty::is_const0( tt0&mk0 ) ) // x & F1
      {
        aig_network::signal x = pis[i];
        pis.erase( pis.begin() + i );
        aig_network::signal f1 = synthesize_aig_rec( aig, pis, tt1, mk1 );
        return aig.create_and( x, f1 );
      }

      if( kitty::is_const0( tt1&mk1 ) ) // x' & F0
      {
        aig_network::signal x = pis[i];
        pis.erase( pis.begin() + i );
        aig_network::signal f0 = synthesize_aig_rec( aig, pis, tt0, mk0 );
        return aig.create_and( !x, f0 );
      }

      if( kitty::equal( tt0&mk0, mk0 ) ) // x' + F1
      {
        aig_network::signal x = pis[i];
        pis.erase( pis.begin() + i );
        aig_network::signal f1 = synthesize_aig_rec( aig, pis, tt1, mk1 );
        return aig.create_or( !x, f1 );
      }

      if( kitty::equal( tt1&mk1, mk1 ) ) // x + F0
      {
        aig_network::signal x = pis[i];
        pis.erase( pis.begin() + i );
        aig_network::signal f0 = synthesize_aig_rec( aig, pis, tt0, mk0 );
        return aig.create_or( x, f0 );
      }
      
      uint32_t n0 = kitty::count_ones( (~tt) & mk );
      uint32_t n1 = kitty::count_ones( tt & mk );
      impurity = n0*n1;
      if( impurity < best_impurity && ( n0 > 0 || n1 > 0) )
      {
        best_impurity = impurity;
        idx = i;
      }

    }

    if( pis.size() <= 4u )
    {
      return match_aig( aig, pis, tt, mk );
    }

    aig_network::signal x = pis[idx];
    pis.erase( pis.begin() + idx );
    aig_network::signal f1 = synthesize_aig_rec( aig, pis, kitty::cofactor1( tt, aig.pi_index(aig.get_node(pis[idx]) ) ), kitty::cofactor1( mk, aig.pi_index(aig.get_node(pis[idx]) ) ) );
    aig_network::signal f0 = synthesize_aig_rec( aig, pis, kitty::cofactor0( tt, aig.pi_index(aig.get_node(pis[idx]) ) ), kitty::cofactor0( mk, aig.pi_index(aig.get_node(pis[idx]) ) ) );

    return aig.create_ite( x, f1, f0 );
  }

  std::tuple<kitty::static_truth_table<4u>, kitty::static_truth_table<4u>> extract_4functionality( kitty::dynamic_truth_table const& tt, kitty::dynamic_truth_table const& mk )
  {

    std::vector<kitty::dynamic_truth_table> xs;
    std::vector<kitty::static_truth_table<4u>> x4;
    for ( uint32_t i{ 0 }; i < 4u; ++i )
    {
      xs.emplace_back( tt.num_vars() );
      x4.emplace_back( );
      kitty::create_nth_var( xs[i], i );
      kitty::create_nth_var( x4[i], i );
    }

    kitty::static_truth_table<4u> func_s;
    kitty::static_truth_table<4u> care_s;
    auto temp = tt.construct();
    auto temp_s = func_s.construct();

    for ( uint32_t m{ 0u }; m < 16u; ++m )
    {
      temp = temp | ~temp;
      temp_s = temp_s | ~temp_s;

      for ( uint32_t l{ 0u }; l < 4; ++l )
      {
        if ( ( m >> l ) & 0x1 == 0x1 )
        {
          temp &= xs[l];
          temp_s &= x4[l];
        }
        else
        {
          temp &= ~xs[l];
          temp_s &= ~x4[l];
        }
      }

      if ( kitty::count_ones( temp & mk ) > 0 ) // care value
      {
        care_s |= temp_s;

        if ( kitty::count_ones( temp & tt ) > 0 )
        {
          func_s |= temp_s;
        }
      }
    }
    auto rnd_tt = func_s.construct();
    kitty::create_random( rnd_tt, _seed++ );

    func_s |= ( rnd_tt & ~care_s );
    // kitty::print_binary( func_s );printf("\n");

    return std::make_tuple( func_s, care_s );
  }

  aig_network::signal match_aig( aig_network& aig, std::vector<signal<aig_network>> vars, kitty::dynamic_truth_table const& tt, kitty::dynamic_truth_table const& mk )
  {

    kitty::static_truth_table<4u> tt4;
    kitty::static_truth_table<4u> mk4;
    auto config4 = extract_4functionality( tt, mk );
    tt4 = std::get<0>(config4);
    mk4 = std::get<1>(config4);

    auto config = exact_npn_canonization( tt4 );

    auto func_npn = std::get<0>( config );
    auto neg = std::get<1>( config );
    auto perm = std::get<2>( config );

    auto dc_npn = ~apply_npn_transformation( mk4, neg & ~( 1 << 4u ), perm );

    auto const structures = _lib.get_supergates( func_npn, dc_npn, neg, perm );

    bool phase = ( neg >> 4 == 1 ) ? true : false;

    for( auto i{0}; i<vars.size(); ++i )
    {
      if( ( neg >> i ) & 0x1 == 0x1 )
        vars[i] = !vars[i];
    }
    std::array<signal<aig_network>, 4> leaves {aig.get_constant(false)};

    for( auto i{0}; i<4; ++i )
    {
      if( perm[i]<vars.size() )
        leaves[i] = vars[perm[i]];
    }
    auto & db = _lib.get_database();
    index_list_t mapped_index_list = index_list;

    auto res = create_aig( aig, db.get_node( structures->at(0).root ), leaves );

    bool is_output_negated = ( phase != db.is_complemented( structures->at(0).root ) );

    return is_output_negated ? !res : res;

  }

  template< class node_t >
  signal<aig_network> create_aig( aig_network& aig, node_t const& n, std::array<signal<aig_network>, 4u> const& leaves )
  {
    return create_aig_rec( aig, n, leaves );
  }

  template< class node_t >
  signal<aig_network> create_aig_rec( aig_network& aig, node_t const& n, std::array<signal<aig_network>, 4u> const& leaves )
  {
    auto& db = _lib.get_database();
 
    std::array<signal<aig_network>, 2u> node_data;

    int i{0};
    db.foreach_fanin( n, [&]( auto const& f, auto i ) 
    {
      node<aig_network> g = db.get_node( f );
      if( db.is_pi( g ) )
      {
        node_data[i] = db.is_complemented( f ) ? !leaves[f.index-1] : leaves[f.index-1];
        return;
      }
      if( db.is_and( g ) )
      {
        auto res = create_aig_rec( aig, g, leaves );
        node_data[i] = db.is_complemented( f ) ? !res : res;
        return;
      }
    } );

    if( db.is_and( n ) )
    {
      signal<aig_network> new_sig = aig.create_and( node_data[0], node_data[1] );
      return new_sig;
    }
    
    return aig.get_constant(false);
  }

#pragma endregion synthesize aig


  std::optional<uint32_t> _map_with_database( std::vector<uint32_t> const& supp, kitty::static_truth_table<4u> const& func, kitty::static_truth_table<4u> const& care, double max_inserts )
  {
    std::array<uint32_t,4> lits0={0,0,0,0};
    for ( int i{0}; i<supp.size(); ++i )
    {
      lits0[i] = (supp[i]+1) << 1u ;
      //printf("l0=%d\n", (supp[i]+1) << 1u );
    }

    auto dcset = ~care;

  //  printf("\n");
    //kitty::print_binary(func);
    //printf("\n");

    std::vector<uint32_t> dcs;
    for( int bit{0}; bit<16; ++bit )
    {
      if( kitty::get_bit( dcset, bit ) > 0 )
      {
        dcs.push_back(bit);
      }
    }


    uint64_t best_key;
    double best_area = max_inserts+1;
    std::vector<uint8_t> best_perm;

    for( uint32_t m{0}; m< (1<<dcs.size()); ++m )
    {
      auto tt = func;

      for( int i{0}; i<dcs.size(); ++i )
      {
        if( (m >> i)&0x1 == 0x1 )
        {
          kitty::flip_bit( tt, dcs[i] );
        }
      }
      //kitty::print_binary(tt);
      const auto support = kitty::min_base_inplace( tt );

     // printf("%d\n", support.size());
      if( support.size() != supp.size() )
        continue;
      /* p-canonize */
      auto config = exact_p_canonization( tt );

      auto func_p = std::get<0>( config );
      auto neg = std::get<1>( config );
      auto perm = std::get<2>( config );

      //kitty::print_binary(func_p);
      //printf("\n");

      uint64_t key = _pClassMap[ func_p._bits&0xFFFF ];
      if( _areas[key] <= best_area  )
      {
        best_key = key;
        best_area = _areas[key];
        best_perm = perm;
      }
    }

    //printf("%f >? %f\n", max_inserts, best_area);
    if( best_area <= max_inserts )
    {
      std::vector<uint32_t> lits = {0,0,0,0,0};
      for( auto i{0}; i<4; ++i )
      {
        //printf("p[%d]=%d\n", i, best_perm[i]);
        lits[i+1] = lits0[best_perm[i]];
        //printf("V[%d]=%d\n", i, lits0[best_perm[i]]);
      }

      //printf( "try match %f %f\n", best_area, max_inserts );
      auto entry = _idlists[ best_key ];
      int type = 0;
      int nFins = 0;
      uint32_t sc_id;
      std::vector<uint32_t> children;
      uint32_t lit;
      //printf("build index list from entry of size %d\n", entry.size());
      for( int i{0}; i<entry.size(); ++i )
      {
        if( type == 0 )
        {
          nFins = entry[i];
          //printf("nFins=%d\n", nFins);
          type=1;
        }
        else if( type == 1 )
        {
          children.push_back( lits[entry[i]] );//not accounting for 0
          //printf("%d ", lits[entry[i]] );
          if( children.size() == nFins )
          {
            //printf("\n");
            type = 2;
          }
        }
        else if( type == 2 )
        {
          type = 0;
          sc_id = entry[i];
          lit = index_list.add_function( children, _gates[sc_id].function, _gates[sc_id].area, _gates[sc_id].id );
          //kitty::print_binary( _gates[sc_id].function );
          //printf("\n");
          lits.push_back( lit );
          children.clear();
        }
      }
      //printf(" %f %f %f\n", index_list.get_area(), _areas[best_key], max_inserts );
      return lit;
    }



    return std::nullopt;
  }

 

  std::optional<uint32_t> _map_on_the_fly( std::vector<uint32_t> const& supp, kitty::dynamic_truth_table const& func, kitty::dynamic_truth_table const& care, double max_inserts )
  {

    std::vector<uint32_t> lits;
    aig_network aig;
    std::vector<signal<aig_network>> pis;
    for ( uint32_t x : supp )
    {
      lits.push_back( (x+1) << 1u );
      pis.push_back( aig.create_pi() );
    }

    auto sig_out = synthesize_aig_inplace( aig, pis, func, care );
    if( aig.is_constant( aig.get_node(sig_out) ) )
      return std::nullopt;

    scopt::emap2_params ps2;
    ps2.cut_enumeration_ps.minimize_truth_table = true;
    ps2.cut_enumeration_ps.cut_limit = 1;
    ps2.cut_enumeration_ps.cut_limit = 1;
    ps2.area_oriented_mapping = true;
    scopt::emap2_stats st2;

    scopt::scg_network scg = scopt::emap2_klut( aig, _tech_lib, ps2, &st2 );

    //printf("substitute if %f<=%f\n", scg.compute_area(), max_inserts );
    if( scg.compute_area() <= max_inserts )
    {      
      scg.foreach_pi( [&]( auto n, auto i ) { 
        scg.set_value( n, lits[i] );
      } );

      std::vector<uint32_t> children;
      uint32_t lit_out;
      scg.foreach_gate( [&]( auto n ) 
      { 
        children.clear();
        scg.foreach_fanin( n, [&]( auto const& fi ) {
          if( scg.is_complemented(fi) )
            children.push_back( scg.value( scg.get_node(fi) ) ^ 0x1 );
          else
            children.push_back( scg.value( scg.get_node(fi) ) );
        } );
        lit_out = index_list.add_function( children, scg.node_function( n ), scg.get_binding(n).area, scg.get_binding(n).id );

        scg.set_value( n, lit_out );
      } );
      return lit_out;
    }

    //  depth_view<scg_network> scg_d{ scg };

    //  if( lit_out && _decomposer.num_luts() <= max_inserts )
    //  {
    //    return _decomposer.to_index_list( index_list, lits );
    //  }

    return std::nullopt;
  }

  // index_list.add_function( lits, func );

#pragma endregion Synthesis

#pragma region support_selection

  template<class SCOREDIVS>
  static std::vector<uint32_t> find_greedy_from_unbalancing( const typename static_params::truth_table_storage_type* pTts, SCOREDIVS const& scored_divisors, std::vector<std::conditional_t<static_params::copy_tts, TT, typename static_params::node_type>> const& divs, spfd_covering_manager_t<truth_table_t, 1 << static_params::max_support_size>& uSPFD, uint32_t pivot, bool complement, bool use_pivot )
  {

    if ( pivot >= scored_divisors.size() )
      return std::vector<uint32_t>{};
    std::mt19937 scgrng;
    scgrng.seed( pivot );

    uint32_t cost, best_cost;
    std::vector<uint32_t> best_candidates;
    std::vector<uint32_t> supp;
    auto const& mask = ( *pTts )[divs[scored_divisors[pivot].div]];
    uSPFD.reset( mask, complement );

    if ( use_pivot )
    {
      supp.push_back( scored_divisors[pivot].div );
    }

    /* add recomputation of the support */
    int nAttempts = 0;
    while ( !uSPFD.is_covered() && nAttempts < static_params::max_support_size )
    {
      nAttempts++;
      best_cost = std::numeric_limits<uint32_t>::max();
      if ( uSPFD.is_saturated() )
        break;
      for ( uint32_t iCnd{ 0 }; iCnd < divs.size(); ++iCnd )
      {
        cost = uSPFD.evaluate( ( *pTts )[divs[iCnd]] );
        if ( cost < best_cost )
        {
          best_cost = cost;
          best_candidates = { iCnd };
        }
        else if ( cost == best_cost )
        {
          best_candidates.push_back( iCnd );
        }
      }
      if ( best_candidates.size() == 0 )
        break;

      std::uniform_int_distribution<> distrib( 0, best_candidates.size() - 1 );
      int idx = distrib( scgrng );
      supp.push_back( best_candidates[idx] );
      uSPFD.update( ( *pTts )[divs[best_candidates[idx]]] );
    }

    if ( uSPFD.is_covered() && supp.size() <= static_params::max_support_size )
    {

      uSPFD.reset();
      for ( auto x : supp )
        uSPFD.update( ( *pTts )[divs[x]] );

      if ( uSPFD.is_covered() )
      {
        std::sort( supp.begin(), supp.end() );
        return supp;
      }
    }
    return std::vector<uint32_t>{};
  }

  static std::vector<uint32_t> find_from_unbalancing( const typename static_params::truth_table_storage_type* pTts, std::vector<scored_div> const& scored_divisors, std::vector<std::conditional_t<static_params::copy_tts, TT, typename static_params::node_type>> const& divs, spfd_covering_manager_t<truth_table_t, 1 << static_params::max_support_size>& uSPFD, uint32_t pivot )
  {
    auto supp1p = find_greedy_from_unbalancing( pTts, scored_divisors, divs, uSPFD, pivot, false, true );
    if ( pivot < divs.size() && supp1p.size() > 0 )
    {
      return supp1p;
    }

    auto supp0p = find_greedy_from_unbalancing( pTts, scored_divisors, divs, uSPFD, pivot, true, true );
    if ( pivot < divs.size() && supp0p.size() > 0 )
    {
      return supp0p;
    }

    auto supp1f = find_greedy_from_unbalancing( pTts, scored_divisors, divs, uSPFD, pivot, false, false );
    if ( pivot < divs.size() && supp1f.size() > 0 )
    {
      return supp1f;
    }

    auto supp0f = find_greedy_from_unbalancing( pTts, scored_divisors, divs, uSPFD, pivot, true, false );
    if ( pivot < divs.size() && supp0f.size() > 0 )
    {
      return supp0f;
    }

    return std::vector<uint32_t>{};
  }

  static std::vector<uint32_t> find_from_funbalancing( const typename static_params::truth_table_storage_type* pTts, std::vector<fscored_div> const& fscored_divisors, std::vector<std::conditional_t<static_params::copy_tts, TT, typename static_params::node_type>> const& divs, spfd_covering_manager_t<truth_table_t, 1 << static_params::max_support_size>& uSPFD, uint32_t pivot )
  {
    auto supp1p = find_greedy_from_unbalancing( pTts, fscored_divisors, divs, uSPFD, pivot, false, true );
    if ( pivot < divs.size() && supp1p.size() > 0 )
    {
      return supp1p;
    }

    auto supp0p = find_greedy_from_unbalancing( pTts, fscored_divisors, divs, uSPFD, pivot, true, true );
    if ( pivot < divs.size() && supp0p.size() > 0 )
    {
      return supp0p;
    }

    auto supp1f = find_greedy_from_unbalancing( pTts, fscored_divisors, divs, uSPFD, pivot, false, false );
    if ( pivot < divs.size() && supp1f.size() > 0 )
    {
      return supp1f;
    }

    auto supp0f = find_greedy_from_unbalancing( pTts, fscored_divisors, divs, uSPFD, pivot, true, false );
    if ( pivot < divs.size() && supp0f.size() > 0 )
    {
      return supp0f;
    }

    return std::vector<uint32_t>{};
  }

  std::optional<std::vector<uint32_t>> find_support()
  {
    if ( static_params::support_selection == support_selection_t::GREEDY )
    {
      if( _idelays.size() > 0 )
      {
        auto supp = find_support_greedy( 0 );
        if ( supp )
        {
          return *supp;
        }
      }
      else
      {
        auto supp = find_support_greedy( 0 );
        if ( supp )
        {
          return *supp;
        }
      }
      return std::nullopt;
    }
    if ( static_params::support_selection == support_selection_t::NGREEDY )
    {
      if( _idelays.size() > 0 )
      {
        auto supp = find_support_ngreedy_with_delay( 0 );
        if ( supp )
        {
          return *supp;
        }
      }
      else
      {
        auto supp = find_support_ngreedy( 0 );
        if ( supp )
        {
          return *supp;
        }
      }
      return std::nullopt;
    }
    if ( static_params::support_selection == support_selection_t::PIVOT )
    {
      auto supp = find_support_greedy( 0 );
      if ( supp )
        return *supp;

      for ( uint32_t i{ 0 }; i < scored_divs.size() * static_params::fraction_of_10 / 10; ++i )
      {
        auto supp = find_from_unbalancing( i );
        if ( supp )
          return *supp;
      }

      return std::nullopt;
    }
  }

  /*! \brief find support greedy */
  std::optional<std::vector<uint32_t>> find_support_greedy( uint32_t start = 0, std::vector<uint32_t> supp0 = {} )
  {
    uint32_t cost, best_cost;
    std::vector<uint32_t> best_candidates;
    std::vector<uint32_t> supp;

    _uSPFD.reset();
    for ( auto x : supp0 )
    {
      _uSPFD.update( get_div( x ) );
      supp.push_back( x );
    }

    /* add recomputation of the support */
    while ( !_uSPFD.is_covered() && supp.size() < static_params::max_support_size )
    {
      best_cost = std::numeric_limits<uint32_t>::max();
      if ( _uSPFD.is_saturated() )
        break;
      for ( uint32_t iCnd{ start }; iCnd < divisors.size(); ++iCnd )
      {
        cost = _uSPFD.evaluate( get_div( iCnd ) );
        if ( cost < best_cost )
        {
          best_cost = cost;
          best_candidates = { iCnd };
        }
        else if ( cost == best_cost )
        {
          best_candidates.push_back( iCnd );
        }
      }
      if ( best_candidates.size() == 0 )
        break;

      std::uniform_int_distribution<> distrib( 0, best_candidates.size() - 1 );
      int idx = distrib( RIGRNG );
      supp.push_back( best_candidates[idx] );
      _uSPFD.update( get_div( best_candidates[idx] ) );
    }

    if ( _uSPFD.is_covered() && supp.size() <= static_params::max_support_size )
    {
      std::sort( supp.begin(), supp.end() );

      return supp;
    }
    return std::nullopt;
  }
 /*! \brief find support greedy */
 
  std::optional<std::vector<uint32_t>> find_support_ngreedy( uint32_t start = 0, std::vector<uint32_t> supp0 = {} )
  {
    uint32_t cost, best_cost0, best_cost1;
    std::vector<uint32_t> best_costs;
    std::vector<std::vector<uint32_t>> best_cands;

    int nxt=0;

    std::vector<uint32_t> supp;

    _uSPFD.reset();
    for ( auto x : supp0 )
    {
      _uSPFD.update( get_div( x ) );
      supp.push_back( x );
    }

    /* add recomputation of the support */
    while ( !_uSPFD.is_covered() && supp.size() < static_params::max_support_size )
    {
      best_costs.clear();
      best_cands.clear();
      for( int i{0}; i<static_params::nBest; ++i )
      {
        best_costs.push_back( std::numeric_limits<double>::max() );
        best_cands.push_back( {} );
      }

      if ( _uSPFD.is_saturated() )
        break;
      for ( uint32_t iCnd{ start }; iCnd < divisors.size(); ++iCnd )
      {
        cost = _uSPFD.evaluate( get_div( iCnd ) );
        int repl=-1;
        for( int j{0}; j<static_params::nBest; ++j )
        {
          if( best_costs[j]>=cost )
          {
            repl=j;
          }
        }
        if( repl >= 0 )
        {
          if( best_costs[repl]==cost )
          {
            best_cands[repl].push_back(iCnd);
          }
          else
          {
            for( int j{0}; j<repl; ++j )
            {
              best_cands[j]=best_cands[j+1];
              best_costs[j]=best_costs[j+1];
            }  
            best_cands[repl]={iCnd};
            best_costs[repl]=cost;
          }
        }       
      }

      std::vector<uint32_t> best_candidates;

      for( auto cands : best_cands )
      {
        for( auto cand : cands )
          best_candidates.push_back( cand );
      }

      std::uniform_int_distribution<> distrib( 0, best_candidates.size() - 1 );
      int idx = distrib( RIGRNG );
      supp.push_back( best_candidates[idx] );
      _uSPFD.update( get_div( best_candidates[idx] ) );

    }

    if ( _uSPFD.is_covered() && supp.size() <= static_params::max_support_size )
    {
      std::sort( supp.begin(), supp.end() );

      return supp;
    }
    return std::nullopt;
  }

 
  std::optional<std::vector<uint32_t>> find_support_ngreedy_with_delay( uint32_t start = 0, std::vector<uint32_t> supp0 = {} )
  {
    uint32_t cost, best_cost0, best_cost1;
    std::vector<uint32_t> best_costs;
    std::vector<std::vector<uint32_t>> best_cands;

    int nxt=0;

    std::vector<uint32_t> supp;

    _uSPFD.reset();
    for ( auto x : supp0 )
    {
      _uSPFD.update( get_div( x ) );
      supp.push_back( x );
    }

    /* add recomputation of the support */
    while ( !_uSPFD.is_covered() && supp.size() < static_params::max_support_size )
    {
      best_costs.clear();
      best_cands.clear();
      for( int i{0}; i<static_params::nBest; ++i )
      {
        best_costs.push_back( std::numeric_limits<double>::max() );
        best_cands.push_back( {} );
      }

      if ( _uSPFD.is_saturated() )
        break;
      for ( uint32_t iCnd{ start }; iCnd < divisors.size(); ++iCnd )
      {
        cost = _uSPFD.evaluate( get_div( iCnd ) );
        int repl=-1;
        for( int j{0}; j<static_params::nBest; ++j )
        {
          if( best_costs[j]>=cost )
          {
            repl=j;
          }
        }
        if( repl >= 0 )
        {
          if( best_costs[repl]==cost )
          {
            best_cands[repl].push_back(iCnd);
          }
          else
          {
            for( int j{0}; j<repl; ++j )
            {
              best_cands[j]=best_cands[j+1];
              best_costs[j]=best_costs[j+1];
            }  
            best_cands[repl]={iCnd};
            best_costs[repl]=cost;
          }
        }       
      }

      std::vector<uint32_t> best_candidates;
      std::vector<double> T;
      double Tmax=std::numeric_limits<double>::min();
      double Tmin=std::numeric_limits<double>::max();

      for( auto cands : best_cands )
      {
        for( auto cand : cands )
        {
          best_candidates.push_back( cand );
          if( _idelays[cand] > Tmax )
            Tmax=_idelays[cand];
          else if( _idelays[cand] < Tmin )
            Tmin=_idelays[cand];
          T.push_back( _idelays[cand] );
        }
      }

      double bestT=std::numeric_limits<double>::max();
      int idx=0;
      for( int i{0}; i<T.size(); ++i )
      {
        if( T[i]<bestT )
        {
          idx=i;
          bestT=T[i];
        }
      }
      supp.push_back( best_candidates[idx] );
      _uSPFD.update( get_div( best_candidates[idx] ) );

    }

    if ( _uSPFD.is_covered() && supp.size() <= static_params::max_support_size )
    {
      std::sort( supp.begin(), supp.end() );

      return supp;
    }
    return std::nullopt;
  }

  

  /*! \brief find support from unbalancing */
  std::optional<std::vector<uint32_t>> find_from_unbalancing( uint32_t pivot )
  {
    uint32_t div = scored_divs[pivot].div;
    auto tti = get_div( div );

    auto supp1p = find_greedy_from_unbalancing( pivot, false, true );
    if ( supp1p )
    {
      return *supp1p;
    }
    auto supp1f = find_greedy_from_unbalancing( pivot, false, false );
    if ( supp1f )
    {
      return *supp1f;
    }

    auto supp0p = find_greedy_from_unbalancing( pivot, true, true );
    if ( supp0p )
    {
      return *supp0p;
    }
    auto supp0f = find_greedy_from_unbalancing( pivot, true, false );
    if ( supp0f )
    {
      return *supp0f;
    }

    return std::nullopt;
  }

  std::optional<std::vector<uint32_t>> find_greedy_from_unbalancing( uint32_t pivot, bool complement, bool use_pivot )
  {
    uint32_t cost, best_cost;
    std::vector<uint32_t> best_candidates;
    std::vector<uint32_t> supp;
    auto const& mask = get_div( scored_divs[pivot].div );
    _uSPFD.reset( mask, complement );

    if ( use_pivot )
    {
      supp.push_back( scored_divs[pivot].div );
    }

    /* add recomputation of the support */
    int nAttempts = 0;
    while ( !_uSPFD.is_covered() && nAttempts < static_params::max_support_size )
    {
      nAttempts++;
      best_cost = std::numeric_limits<uint32_t>::max();
      if ( _uSPFD.is_saturated() )
        break;
      for ( uint32_t iCnd{ pivot + 1 }; iCnd < divisors.size(); ++iCnd )
      {
        cost = _uSPFD.evaluate( get_div( iCnd ) );
        if ( cost < best_cost )
        {
          best_cost = cost;
          best_candidates = { iCnd };
        }
        else if ( cost == best_cost )
        {
          best_candidates.push_back( iCnd );
        }
      }
      if ( best_candidates.size() == 0 )
        break;

      std::uniform_int_distribution<> distrib( 0, best_candidates.size() - 1 );
      int idx = distrib( RIGRNG );
      supp.push_back( best_candidates[idx] );
      _uSPFD.update( get_div( best_candidates[idx] ) );
    }

    if ( _uSPFD.is_covered() && supp.size() <= static_params::max_support_size )
    {

      _uSPFD.reset();
      for ( auto x : supp )
        _uSPFD.update( get_div( x ) );

      if ( _uSPFD.is_covered() )
      {
        std::sort( supp.begin(), supp.end() );
        return supp;
      }
    }
    return std::nullopt;
  }

#pragma endregion support_selection

  inline TT const& get_div( uint32_t idx ) const
  {
    if constexpr ( static_params::copy_tts )
    {
      return divisors[idx];
    }
    else
    {
      return ( *ptts )[divisors[idx]];
    }
  }

private:
  std::array<TT, 2> on_off_sets;
  std::array<uint32_t, 2> num_bits; /* number of bits in on-set and off-set */

  const typename static_params::truth_table_storage_type* ptts;
  std::vector<std::conditional_t<static_params::copy_tts, TT, typename static_params::node_type>> divisors;

  spfd_covering_manager_t<truth_table_t, 1 << static_params::max_support_size> _uSPFD;
  lut_resynthesis_t<static_params::max_fanin_size, static_params::max_support_size> _decomposer;

  index_list_t index_list;

  /* positive unate: not overlapping with off-set
     negative unate: not overlapping with on-set */
  std::vector<scored_div> scored_divs;
  std::vector<fscored_div> fscored_divs;

  stats& st;

  std::default_random_engine::result_type _seed = 1;

  xag_npn_resynthesis<aig_network, aig_network, xag_npn_db_kind::aig_complete> _resyn;

  exact_library_params _lps;
  exact_library<aig_network> _lib;

  tech_library_params _tps;
  std::vector<gate> _gates;
  tech_library<5, classification_type::np_configurations> _tech_lib;
  std::vector<std::vector<uint32_t>> _idlists;
  std::vector<double> _areas;
  std::unordered_map<uint64_t, uint32_t> _pClassMap;
  std::vector<double> _idelays;

  double _area_th;

}; /* xag_resyn_decompose */

}; /* namespace rils */

} /* namespace mockturtle */

