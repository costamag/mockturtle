/* mockturtle: C++ logic network library
 * Copylight (C) 2018-2022  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the lights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copylight notice and this permission notice shall be
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
  \file lig_resyn.hpp
  \brief Resynthesis by extraction of functional cuts

  \author Andrea Costamagna
*/

#pragma once

#include "../../utils/index_list.hpp"
#include "../../utils/node_map.hpp"
#include "../../utils/stopwatch.hpp"
#include "../../utils/spfd_utils.hpp"

#include <abcresub/abcresub.hpp>
#include <fmt/format.h>
#include <kitty/kitty.hpp>

#include <algorithm>
#include <optional>
#include <type_traits>
#include <vector>
#include <thread>
#include <mutex>

namespace mockturtle
{

namespace rils
{

bool VERBOSE{false};

enum support_selection_t
{
  GREEDY,
  PIVOT,
};

std::mt19937 RIGRNG(5);

struct lig_resyn_static_params
{
  using base_type = lig_resyn_static_params;

  /*! \brief Whether to copy truth tables. */
  static constexpr bool copy_tts{ false };

  /*! \brief Reserved capacity for divisor truth tables (number of divisors). */
  static constexpr uint32_t reserve{ 200u };

  /*! \brief Whether to preserve depth. */
  static constexpr bool preserve_depth{ false };

  /*! \brief Whether the divisors have uniform costs (size and depth, whenever relevant). */
  static constexpr bool uniform_div_cost{ true };

  static constexpr uint32_t max_support_size{6u};
  static constexpr uint32_t fraction_of_10{10};

  static constexpr int max_fanin_size = -1;
  static constexpr bool accept_worse{false};

  static constexpr support_selection_t support_selection{ GREEDY };

  using truth_table_storage_type = void;
  using node_type = void;
};

template<class TT>
struct lig_resyn_static_params_default : public lig_resyn_static_params
{
  using truth_table_storage_type = std::vector<TT>;
  using node_type = uint32_t;
};


template<class Ntk, support_selection_t SUP_SEL, uint32_t SUPP_SIZE, int K=-1, int NRELAX=0>
struct lig_resyn_static_params_for_sim_resub : public lig_resyn_static_params
{
  using truth_table_storage_type = incomplete_node_map<kitty::partial_truth_table, Ntk>;
  using node_type = typename Ntk::node;
  static constexpr support_selection_t support_selection = SUP_SEL;
  static constexpr uint32_t max_support_size = SUPP_SIZE;
  static constexpr int max_fanin_size = K;
  static constexpr bool accept_worse = NRELAX>0;
};

template<class Ntk, support_selection_t SUP_SEL, uint32_t NumVars, uint32_t SUPP_SIZE, int K=-1, int NRELAX=0>
struct lig_resyn_static_params_for_sim_resub_static : public lig_resyn_static_params
{
  using truth_table_storage_type = incomplete_node_map<kitty::static_truth_table<NumVars>, Ntk>;
  using node_type = typename Ntk::node;
  static constexpr support_selection_t support_selection = SUP_SEL;
  static constexpr uint32_t max_support_size = SUPP_SIZE;
  static constexpr int max_fanin_size = K;
  static constexpr bool accept_worse = NRELAX>0;
};

struct lig_resyn_stats
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

template<class TT, class static_params = lig_resyn_static_params_default<TT>, support_selection_t SUP_SEL=GREEDY>
class lig_resyn_decompose
{

public:
  using stats = lig_resyn_stats;
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

  explicit lig_resyn_decompose( stats& st ) noexcept
      : st( st )
  {
    static_assert( std::is_same_v<typename static_params::base_type, lig_resyn_static_params>, "Invalid static_params type" );
    static_assert( !( static_params::uniform_div_cost && static_params::preserve_depth ), "If depth is to be preserved, divisor depth cost must be provided (usually not uniform)" );
    divisors.reserve( static_params::reserve );
    RIGRNG.seed(5);

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
  std::optional<index_list_t> operator()( TT const& target, TT const& care, iterator_type begin, iterator_type end, typename static_params::truth_table_storage_type const& tts, uint32_t max_size )
  {
    static_assert( static_params::copy_tts || std::is_same_v<typename std::iterator_traits<iterator_type>::value_type, typename static_params::node_type>, "iterator_type does not dereference to static_params::node_type" );

    ptts = &tts;
    on_off_sets[0] = ~target & care;
    on_off_sets[1] = target & care;

    _uSPFD.init( target, care );

    divisors.resize( 1 ); /* clear previous data and reserve 1 dummy node for constant */
    scored_divs.clear();

    while ( begin != end )
    {
      if constexpr ( static_params::copy_tts )
      {
        divisors.emplace_back( ( *ptts )[*begin] );
        scored_divs.emplace_back( divisors.size()-1, _uSPFD.evaluate( get_div( divisors.size()-1 ) ) );
      }
      else
      {
        divisors.emplace_back( *begin );
        scored_divs.emplace_back( divisors.size()-1, _uSPFD.evaluate( get_div( divisors.size()-1 ) ) );
      }
      ++begin;
    }

    call_with_stopwatch( st.time_sort, [&]() {
      std::sort( scored_divs.begin(), scored_divs.end() ); 
    });

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
  std::optional<index_list_t> compute_function( uint32_t num_inserts )
  {
    index_list.clear();
    index_list.reset_area();
    index_list.add_inputs( divisors.size() - 1 );
    auto const lit = compute_function_rec( num_inserts );
    if ( lit )
    {
      index_list.add_output( *lit );
      return index_list;
    }
    return std::nullopt;
  }

  std::optional<uint32_t> compute_function_rec( uint32_t num_inserts )
  {

    /* try 0-resub and collect unate literals */
    auto const res0 = call_with_stopwatch( st.time_0resub, [&]() {
      return try_0resub( num_inserts );
    } );
    if ( res0 )
    {
      return *res0;
    }
    
    if ( num_inserts <= 0 )
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
  std::optional<uint32_t> try_0resub( uint32_t max_inserts = 0 )
  {
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

    for ( auto v = 1u; v < divisors.size(); ++v )
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
        return ( v << 1 );
      }
      if ( unateness[1] && unateness[2] )
      {
        return ( v << 1 ) + 1u;
      }
    }
    return std::nullopt;
  }

  /* See if we cna define a new function of the other divisors
   */
  std::optional<uint32_t> try_nresub( uint32_t max_inserts )
  {
    auto supp = find_support();

    if( supp )
    {
      auto const [func, care] = extract_functionality_from_signatures( *supp );
      return _1_node_synthesis( *supp, func, care, max_inserts );
    }
    /* resynthesis */
    return std::nullopt;
  }

  std::tuple<kitty::dynamic_truth_table, kitty::dynamic_truth_table> extract_functionality_from_signatures( std::vector<uint32_t> const& supp )
  {
    assert( supp.size() <= static_params::max_support_size );

    std::vector<kitty::dynamic_truth_table> xs;
    for( uint32_t i{0}; i<supp.size(); ++i )
    {
      xs.emplace_back(supp.size());
      kitty::create_nth_var( xs[i], i );
    }

    kitty::dynamic_truth_table func_s(supp.size());
    kitty::dynamic_truth_table care_s = func_s.construct();
    auto  temp = _uSPFD.care.construct();
    auto  temp_s = func_s.construct();

    for( uint32_t m{0u}; m < ( 1u << supp.size() ); ++m )
    {
      temp = temp | ~temp;
      temp_s = temp_s | ~temp_s;

      for( uint32_t l{0u}; l < supp.size(); ++l )
      {
        if( ( m >> l ) & 0x1 == 0x1 )
        {
          temp &= get_div(supp[l]);
          temp_s &= xs[l];
        }
        else
        {
          temp &= ~get_div(supp[l]);
          temp_s &= ~xs[l];
        }
      }

      if( kitty::count_ones( temp & _uSPFD.care ) > 0 ) // care value
      {
        care_s |= temp_s;

        if( kitty::count_ones( temp & _uSPFD.func[1] ) > 0 )
        {
          func_s |= temp_s;
        }
      }
    }
    auto rnd_tt = func_s.construct();
    kitty::create_random( rnd_tt, _seed++ );

    func_s |= ( rnd_tt & ~care_s );
    //kitty::print_binary( func_s );printf("\n");

    return std::make_tuple(func_s, care_s);
  }

#pragma region synthesis

std::optional<uint32_t> _1_node_synthesis( std::vector<uint32_t> const& supp, kitty::dynamic_truth_table const& func, kitty::dynamic_truth_table const& care, uint32_t max_inserts )
{
  std::vector<uint32_t> lits;
  for( uint32_t x : supp )
    lits.push_back( x << 1u );


  _decomposer.clear();
  auto lit_out = _decomposer.decompose( func, care, max_inserts );
  if( _decomposer.num_luts() <= max_inserts )
  {
    return _decomposer.to_index_list( index_list, lits );
  }

  return std::nullopt;
}

std::array<uint32_t, 4> compute_literals( std::vector<uint32_t> const& supp )
{
  std::array<uint32_t, 4> lits {0};
  for( int i{0}; i<supp.size(); ++i )
  {
    lits[i] = supp[i] << 1u;
  }
  return lits;
}


//index_list.add_function( lits, func );

#pragma endregion Synthesis

#pragma region support_selection

  template< class SCOREDIVS >
  static std::vector<uint32_t> find_greedy_from_unbalancing( const typename static_params::truth_table_storage_type* pTts, SCOREDIVS const& scored_divisors, std::vector<std::conditional_t<static_params::copy_tts, TT, typename static_params::node_type>> const& divs, spfd_covering_manager_t<truth_table_t, 1 << static_params::max_support_size> & uSPFD, uint32_t pivot, bool complement, bool use_pivot )
  {

    if( pivot >= scored_divisors.size() ) return std::vector<uint32_t>{};
    std::mt19937 ligrng;
    ligrng.seed( pivot );

    uint32_t cost, best_cost;
    std::vector<uint32_t> best_candidates;
    std::vector<uint32_t> supp;
    auto const& mask = ( *pTts )[divs[scored_divisors[pivot].div ]];
    uSPFD.reset( mask, complement );

    if( use_pivot )
    {
      supp.push_back( scored_divisors[pivot].div );
    }

    /* add recomputation of the support */
    int nAttempts=0;
    while( !uSPFD.is_covered() && nAttempts < static_params::max_support_size )
    {
      nAttempts++;
      best_cost = std::numeric_limits<uint32_t>::max();
      if( uSPFD.is_saturated() ) break;
      for( uint32_t iCnd{1}; iCnd<divs.size(); ++iCnd )
      {
        cost = uSPFD.evaluate( ( *pTts )[divs[iCnd]] );
        if( cost < best_cost  )
        {
          best_cost = cost;
          best_candidates = {iCnd};
        }
        else if( cost == best_cost )
        {
          best_candidates.push_back( iCnd );
        }
      }
      if( best_candidates.size() == 0 ) break;

      std::uniform_int_distribution<> distrib(0, best_candidates.size()-1);
      int idx = distrib(ligrng);
      supp.push_back( best_candidates[idx] );
      uSPFD.update( ( *pTts )[divs[best_candidates[idx]]] );
    }

    if( uSPFD.is_covered() && supp.size() <= static_params::max_support_size )
    {

      uSPFD.reset();
      for( auto x : supp )
        uSPFD.update(( *pTts )[divs[x]]);
      
      if( uSPFD.is_covered() )
      {
        std::sort( supp.begin(), supp.end() );
        return supp;
      }
    }
    return std::vector<uint32_t>{};
  }

  static std::vector<uint32_t> find_from_unbalancing( const typename static_params::truth_table_storage_type * pTts, std::vector<scored_div> const& scored_divisors, std::vector<std::conditional_t<static_params::copy_tts, TT, typename static_params::node_type>> const& divs, spfd_covering_manager_t<truth_table_t, 1 << static_params::max_support_size> & uSPFD, uint32_t pivot )
  {
    auto supp1p = find_greedy_from_unbalancing( pTts, scored_divisors, divs, uSPFD, pivot, false, true );
    if( pivot < divs.size() && supp1p.size() > 0 )
    {
      return supp1p;
    }

    auto supp0p = find_greedy_from_unbalancing( pTts, scored_divisors, divs, uSPFD, pivot, true, true );
    if( pivot < divs.size() && supp0p.size() > 0 )
    {
      return supp0p;
    }

    auto supp1f = find_greedy_from_unbalancing( pTts, scored_divisors, divs, uSPFD, pivot, false, false );
    if( pivot < divs.size() && supp1f.size() > 0 )
    {
      return supp1f;
    }

    auto supp0f = find_greedy_from_unbalancing( pTts, scored_divisors, divs, uSPFD, pivot, true, false );
    if( pivot < divs.size() && supp0f.size() > 0 )
    {
      return supp0f;
    }

    return std::vector<uint32_t>{};
  }

  static std::vector<uint32_t> find_from_funbalancing( const typename static_params::truth_table_storage_type * pTts, std::vector<fscored_div> const& fscored_divisors, std::vector<std::conditional_t<static_params::copy_tts, TT, typename static_params::node_type>> const& divs, spfd_covering_manager_t<truth_table_t, 1 << static_params::max_support_size> & uSPFD, uint32_t pivot )
  {
    auto supp1p = find_greedy_from_unbalancing( pTts, fscored_divisors, divs, uSPFD, pivot, false, true );
    if( pivot < divs.size() && supp1p.size() > 0 )
    {
      return supp1p;
    }

    auto supp0p = find_greedy_from_unbalancing( pTts, fscored_divisors, divs, uSPFD, pivot, true, true );
    if( pivot < divs.size() && supp0p.size() > 0 )
    {
      return supp0p;
    }

    auto supp1f = find_greedy_from_unbalancing( pTts, fscored_divisors, divs, uSPFD, pivot, false, false );
    if( pivot < divs.size() && supp1f.size() > 0 )
    {
      return supp1f;
    }

    auto supp0f = find_greedy_from_unbalancing( pTts, fscored_divisors, divs, uSPFD, pivot, true, false );
    if( pivot < divs.size() && supp0f.size() > 0 )
    {
      return supp0f;
    }

    return std::vector<uint32_t>{};
  }


  std::optional<std::vector<uint32_t>> find_support()
  {
    if( static_params::support_selection == support_selection_t::GREEDY )
    {
      auto supp = find_support_greedy(1);
      if( supp )
      {          
        return *supp;
      }
      return std::nullopt;
    }
    if( static_params::support_selection == support_selection_t::PIVOT )
    {
      auto supp = find_support_greedy();
      if( supp )
        return *supp;
      
      for( uint32_t i{0}; i<scored_divs.size()*static_params::fraction_of_10/10; ++i )
      {
        auto supp = find_from_unbalancing(i);
        if( supp )
          return *supp;
      }

      return std::nullopt;
    }
  }

  /*! \brief find support greedy */
  std::optional<std::vector<uint32_t>> find_support_greedy( uint32_t start=1, std::vector<uint32_t> supp0 = {} )
  {
    uint32_t cost, best_cost;
    std::vector<uint32_t> best_candidates;
    std::vector<uint32_t> supp;
    
    _uSPFD.reset();
    for( auto x : supp0 )
    {
      _uSPFD.update( get_div(x) );
      supp.push_back(x);
    }

    /* add recomputation of the support */
    while( !_uSPFD.is_covered() && supp.size() < static_params::max_support_size )
    {
      best_cost = std::numeric_limits<uint32_t>::max();
      if( _uSPFD.is_saturated() ) break;
      for( uint32_t iCnd{start}; iCnd<divisors.size(); ++iCnd )
      {
        cost = _uSPFD.evaluate( get_div( iCnd ) );
        if( cost < best_cost  )
        {
          best_cost = cost;
          best_candidates = {iCnd};
        }
        else if( cost == best_cost )
        {
          best_candidates.push_back( iCnd );
        }
      }
      if( best_candidates.size() == 0 ) break;

      std::uniform_int_distribution<> distrib(0, best_candidates.size()-1);
      int idx = distrib(RIGRNG);
      supp.push_back( best_candidates[idx] );
      _uSPFD.update( get_div( best_candidates[idx] ) );
    }

    if( _uSPFD.is_covered() && supp.size() <= static_params::max_support_size )
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
    if( supp1p )
    {
      return *supp1p;
    }
    auto supp1f = find_greedy_from_unbalancing( pivot, false, false );
    if( supp1f )
    {
      return *supp1f;
    }

    auto supp0p = find_greedy_from_unbalancing( pivot, true, true );
    if( supp0p )
    {
      return *supp0p;
    }
    auto supp0f = find_greedy_from_unbalancing( pivot, true, false );
    if( supp0f )
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

    if( use_pivot )
    {
      supp.push_back( scored_divs[pivot].div );
    }

    /* add recomputation of the support */
    int nAttempts=0;
    while( !_uSPFD.is_covered() && nAttempts < static_params::max_support_size )
    {
      nAttempts++;
      best_cost = std::numeric_limits<uint32_t>::max();
      if( _uSPFD.is_saturated() ) break;
      for( uint32_t iCnd{pivot+1}; iCnd<divisors.size(); ++iCnd )
      {
        cost = _uSPFD.evaluate( get_div( iCnd ) );
        if( cost < best_cost  )
        {
          best_cost = cost;
          best_candidates = {iCnd};
        }
        else if( cost == best_cost )
        {
          best_candidates.push_back( iCnd );
        }
      }
      if( best_candidates.size() == 0 ) break;

      std::uniform_int_distribution<> distrib(0, best_candidates.size()-1);
      int idx = distrib(RIGRNG);
      supp.push_back( best_candidates[idx] );
      _uSPFD.update( get_div( best_candidates[idx] ) );
    }

    if( _uSPFD.is_covered() && supp.size() <= static_params::max_support_size )
    {

      _uSPFD.reset();
      for( auto x : supp )
        _uSPFD.update(get_div(x));
      
      if( _uSPFD.is_covered() )
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

  std::default_random_engine::result_type _seed=1;

}; /* xag_resyn_decompose */

}; /* namespace rils */

} /* namespace mockturtle */


