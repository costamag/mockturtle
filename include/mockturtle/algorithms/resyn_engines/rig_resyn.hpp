/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
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
  \file xag_resyn.hpp
  \brief Resynthesis by recursive decomposition for AIGs or XAGs.
  (based on ABC's implementation in `giaResub.c` by Alan Mishchenko)

  \author Andrea Costamagna
*/

#pragma once

#include "../../utils/index_list.hpp"
#include "../../utils/node_map.hpp"
#include "../../utils/stopwatch.hpp"

#include <abcresub/abcresub.hpp>
#include <fmt/format.h>
#include <kitty/kitty.hpp>

#include <algorithm>
#include <optional>
#include <type_traits>
#include <vector>

namespace mockturtle
{

namespace rils
{

enum support_selection_t
{
  GREEDY,
  ALL,
  GREALL,
  GREEDY2,
  GREEDYN,
  GREEDY3,
  PIVOT,
  RANDOM
};

struct rig_resyn_static_params
{
  using base_type = rig_resyn_static_params;

  /*! \brief Maximum number of binate divisors to be considered. */
  static constexpr uint32_t max_binates{ 50u };

  /*! \brief Reserved capacity for divisor truth tables (number of divisors). */
  static constexpr uint32_t reserve{ 200u };

  /*! \brief Whether to consider single XOR gates (i.e., using XAGs instead of AIGs). */
  static constexpr bool use_xor{ true };

  /*! \brief Whether to copy truth tables. */
  static constexpr bool copy_tts{ false };

  /*! \brief Whether to preserve depth. */
  static constexpr bool preserve_depth{ false };

  /*! \brief Whether the divisors have uniform costs (size and depth, whenever relevant). */
  static constexpr bool uniform_div_cost{ true };

  /*! \brief Size cost of each AND gate. */
  static constexpr uint32_t size_cost_of_and{ 1u };

  /*! \brief Size cost of each XOR gate (only relevant when `use_xor = true`). */
  static constexpr uint32_t size_cost_of_xor{ 1u };

  /*! \brief Depth cost of each AND gate (only relevant when `preserve_depth = true`). */
  static constexpr uint32_t depth_cost_of_and{ 1u };

  /*! \brief Depth cost of each XOR gate (only relevant when `preserve_depth = true` and `use_xor = true`). */
  static constexpr uint32_t depth_cost_of_xor{ 1u };

  static constexpr uint32_t max_support_size{6u};
  static constexpr uint32_t fraction_of_10{10};

  static constexpr support_selection_t support_selection{ GREEDY };

  using truth_table_storage_type = void;
  using node_type = void;
};

template<class TT>
struct rig_resyn_static_params_default : public rig_resyn_static_params
{
  using truth_table_storage_type = std::vector<TT>;
  using node_type = uint32_t;
};


template<class Ntk, support_selection_t SUP_SEL, uint32_t K>
struct rig_resyn_static_params_for_sim_resub : public rig_resyn_static_params
{
  using truth_table_storage_type = incomplete_node_map<kitty::partial_truth_table, Ntk>;
  using node_type = typename Ntk::node;
  static constexpr support_selection_t support_selection = SUP_SEL;
  static constexpr uint32_t max_support_size = K;
};

struct rig_resyn_stats
{
  /*! \brief Time for finding 0-resub and collecting unate literals. */
  stopwatch<>::duration time_unate{ 0 };

  /*! \brief Time for finding 1-resub. */
  stopwatch<>::duration time_resub1{ 0 };

  /*! \brief Time for finding 2-resub. */
  stopwatch<>::duration time_resub2{ 0 };

  /*! \brief Time for finding 3-resub. */
  stopwatch<>::duration time_resub3{ 0 };

  /*! \brief Time for sorting unate literals and unate pairs. */
  stopwatch<>::duration time_sort{ 0 };

  /*! \brief Time for collecting unate pairs. */
  stopwatch<>::duration time_collect_pairs{ 0 };

  /*! \brief Time for dividing the target and recursive call. */
  stopwatch<>::duration time_divide{ 0 };

  void report() const
  {
    fmt::print( "[i]         <xag_resyn_decompose>\n" );
    fmt::print( "[i]             0-resub      : {:>5.2f} secs\n", to_seconds( time_unate ) );
    fmt::print( "[i]             1-resub      : {:>5.2f} secs\n", to_seconds( time_resub1 ) );
    fmt::print( "[i]             2-resub      : {:>5.2f} secs\n", to_seconds( time_resub2 ) );
    fmt::print( "[i]             3-resub      : {:>5.2f} secs\n", to_seconds( time_resub3 ) );
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
template<class TT, class static_params = rig_resyn_static_params_default<TT>, support_selection_t SUP_SEL=GREEDY>
class rig_resyn_decompose
{
private:
  std::mt19937 RIGRNG;//(5);


public:
  using stats = rig_resyn_stats;
  using index_list_t = large_rig_index_list;
  using truth_table_t = TT;

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

  struct fanin_pair
  {
    fanin_pair( uint32_t l1, uint32_t l2 )
        : lit1( l1 < l2 ? l1 : l2 ), lit2( l1 < l2 ? l2 : l1 )
    {}

    fanin_pair( uint32_t l1, uint32_t l2, bool is_xor )
        : lit1( l1 > l2 ? l1 : l2 ), lit2( l1 > l2 ? l2 : l1 )
    {
      (void)is_xor;
    }

    bool operator==( fanin_pair const& other ) const
    {
      return lit1 == other.lit1 && lit2 == other.lit2;
    }

    uint32_t lit1, lit2;
    uint32_t score{ 0 };
  };

  template<class LTT, uint32_t CAP>
  struct u_spfd_manager_t
  {
    u_spfd_manager_t(){}
    
    void init( LTT const& target, LTT const& careset )
    {
      care = careset;
      safe_care = careset;
      func[1] =  target & careset;
      func[0] = ~target & careset;
      reset();
    }

    void reset()
    {
      masks[0] = safe_care;
      nMasks = 1;
      nEdges = kitty::count_ones( func[1] ) * kitty::count_ones( func[0] );  
      killed[0] = nEdges > 0 ? false : true;
      nKills = nEdges > 0 ? 0 : 1;
    }

    void reset( LTT const& modified_care, bool complement )
    {
      masks[0] = complement ? safe_care & ~modified_care : safe_care & modified_care;
      nMasks = 1;
      nEdges = kitty::count_ones( func[1] & masks[0] ) * kitty::count_ones( func[0] & masks[0] );  
      killed[0] = nEdges > 0 ? false : true;
      nKills = nEdges > 0 ? 0 : 1;
    }

    bool update( LTT const& tt )
    {
      nEdges = 0;
      for( uint32_t iMask{0}; iMask < nMasks; ++iMask )
      {
        if( killed[iMask] )
        {
          killed[nMasks+iMask] = true;
          nKills++;
        }
        else
        {
          masks[nMasks+iMask] = masks[iMask] & tt;
          masks[iMask] &= ~tt;

          if( ( kitty::count_ones( masks[nMasks+iMask] & func[1] ) == 0 ) || ( kitty::count_ones( masks[nMasks+iMask] & func[0] ) == 0 ) )
          {
            killed[nMasks+iMask] = true;
            nKills++;
          }
          else
          {
            killed[nMasks+iMask] = false;
            nEdges += kitty::count_ones( func[1] & masks[nMasks+iMask] ) * kitty::count_ones( func[0] & masks[nMasks+iMask] );  
          }

          if( kitty::count_ones( masks[iMask] & func[1] ) == 0 || kitty::count_ones( masks[iMask] & func[0] ) == 0 )
          {
            killed[iMask] = true;
            nKills++;
          }
          else
          {
            killed[iMask] = false;
            nEdges += kitty::count_ones( func[1] & masks[iMask] ) * kitty::count_ones( func[0] & masks[iMask] );  
          }
        }
      }
      nMasks = nMasks * 2;
      return true;
    }

    uint32_t evaluate( LTT const& tt )
    {
      uint32_t res=0;
      for( auto iMask{0}; iMask<nMasks; ++iMask )
      {
        if( !killed[iMask] )
        {
          res+= kitty::count_ones( func[1] & masks[iMask] &  tt ) * kitty::count_ones( func[0] & masks[iMask] & tt );  
          res+= kitty::count_ones( func[1] & masks[iMask] & ~tt ) * kitty::count_ones( func[0] & masks[iMask] & ~tt );  
        }
      }
      return res;
    } 

    uint32_t evaluate( LTT const& tt1, LTT const& tt2 )
    {
      uint32_t res=0;
      for( auto iMask{0}; iMask<nMasks; ++iMask )
      {
        if( !killed[iMask] )
        {
          res+= kitty::count_ones( func[1] & masks[iMask] &  tt1 & tt2 ) * kitty::count_ones( func[0] & masks[iMask] & tt1 & tt2 );  
          res+= kitty::count_ones( func[1] & masks[iMask] & ~tt1 & tt2 ) * kitty::count_ones( func[0] & masks[iMask] & ~tt1 & tt2 );  
          res+= kitty::count_ones( func[1] & masks[iMask] & ~tt1 &~tt2 ) * kitty::count_ones( func[0] & masks[iMask] & ~tt1 & ~tt2 );  
          res+= kitty::count_ones( func[1] & masks[iMask] & tt1 & ~tt2 ) * kitty::count_ones( func[0] & masks[iMask] & tt1 &~tt2 );  
        }
      }
      return res;
    } 

    uint32_t evaluate( LTT const& tt1, LTT const& tt2, LTT const& tt3 )
    {
      uint32_t res=0;
      for( auto iMask{0}; iMask<nMasks; ++iMask )
      {
        if( !killed[iMask] )
        {
          res+= kitty::count_ones( func[1] & masks[iMask] &  tt1 & tt2 & tt3 ) * kitty::count_ones( func[0] & masks[iMask] & tt1 & tt2 & tt3 );  
          res+= kitty::count_ones( func[1] & masks[iMask] & ~tt1 & tt2 & tt3 ) * kitty::count_ones( func[0] & masks[iMask] & ~tt1 & tt2 & tt3 );  
          res+= kitty::count_ones( func[1] & masks[iMask] & ~tt1 &~tt2 & tt3) * kitty::count_ones( func[0] & masks[iMask] & ~tt1 & ~tt2 & tt3 );  
          res+= kitty::count_ones( func[1] & masks[iMask] & tt1 & ~tt2 & tt3) * kitty::count_ones( func[0] & masks[iMask] & tt1 &~tt2 & tt3);  
          res+= kitty::count_ones( func[1] & masks[iMask] &  tt1 & tt2 & ~tt3 ) * kitty::count_ones( func[0] & masks[iMask] & tt1 & tt2 & ~tt3 );  
          res+= kitty::count_ones( func[1] & masks[iMask] & ~tt1 & tt2 & ~tt3) * kitty::count_ones( func[0] & masks[iMask] & ~tt1 & tt2  & ~tt3);  
          res+= kitty::count_ones( func[1] & masks[iMask] & ~tt1 &~tt2 & ~tt3) * kitty::count_ones( func[0] & masks[iMask] & ~tt1 & ~tt2 & ~tt3);  
          res+= kitty::count_ones( func[1] & masks[iMask] & tt1 & ~tt2 & ~tt3) * kitty::count_ones( func[0] & masks[iMask] & tt1 &~tt2   & ~tt3);  
        }
      }
      return res;
    } 

    bool is_covered()
    {
      return nMasks <= nKills;
    }

    bool is_saturated()
    {
      return nMasks >= CAP;
    }

    std::array<LTT, CAP> masks;
    std::array<bool, CAP> killed;
    uint32_t nMasks;
    uint32_t nKills;
    uint32_t nEdges;
    LTT care;
    LTT safe_care;
    std::array<LTT, 2u> func;
  };

public:
  explicit rig_resyn_decompose( stats& st ) noexcept
      : st( st )
  {
    static_assert( std::is_same_v<typename static_params::base_type, rig_resyn_static_params>, "Invalid static_params type" );
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
  std::optional<index_list_t> operator()( TT const& target, TT const& care, iterator_type begin, iterator_type end, typename static_params::truth_table_storage_type const& tts, uint32_t max_size = std::numeric_limits<uint32_t>::max() )
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
    std::sort( scored_divs.begin(), scored_divs.end() );
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
    index_list.add_inputs( divisors.size() - 1 );
    auto const lit = compute_function_rec( num_inserts );
    if ( lit )
    {
      assert( index_list.num_gates() <= num_inserts );
      index_list.add_output( *lit );
      //std::cout << to_index_list_string(index_list) << std::endl;
      return index_list;
    }
    return std::nullopt;
  }

  std::optional<uint32_t> compute_function_rec( uint32_t num_inserts )
  {

    /* try 0-resub and collect unate literals */
    auto const res0 = call_with_stopwatch( st.time_unate, [&]() {
      return find_one_unate();
    } );
    if ( res0 )
    {
      return *res0;
    }
    if ( num_inserts == 0u )
    {
      return std::nullopt;
    }

    /* try 1-resub */
    auto const res1 = call_with_stopwatch( st.time_unate, [&]() {
      return try_1_resub();
    } );
    if ( res1 )
    {
      return *res1;
    }
    if ( num_inserts == 1u )
    {
      return std::nullopt;
    }

    return std::nullopt;
  }

  /* See if there is a constant or divisor covering all on-set bits or all off-set bits.
     1. Check constant-resub
     2. Collect unate literals
     3. Find 0-resub (both positive unate and negative unate) and collect binate (neither pos nor neg unate) divisors
   */
  std::optional<uint32_t> find_one_unate()
  {

    //return std::nullopt;

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
        return ( v << 1 ) + 1;
      }
    }
    return std::nullopt;
  }

  /* See if we cna define a new function of the other divisors
   */
  std::optional<uint32_t> try_1_resub()
  {
    /* support selection */
    auto supp = find_support();
    if( supp )
    {
      auto func = extract_functionality_from_signatures( *supp );

      std::vector<uint32_t> lits;
      for( uint32_t x : *supp )
        lits.push_back( x << 1u );

      return index_list.add_function( lits, func );

    }
    /* resynthesis */
    return std::nullopt;
  }

  kitty::dynamic_truth_table extract_functionality_from_signatures( std::vector<uint32_t> const& supp )
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

    return func_s;
  }

  std::optional<std::vector<uint32_t>> find_support()
  {
    if( static_params::support_selection == support_selection_t::GREEDY )
    {
      return find_support_greedy(1);
    }
    if( static_params::support_selection == support_selection_t::GREEDYN )
    {
      for( uint32_t i{0}; i<divisors.size(); ++i )
      {
        auto supp = find_support_greedy(i);
        if( supp )
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
    else if( static_params::support_selection == support_selection_t::ALL )
    {
      return find_support_all();
    }
    else if( static_params::support_selection == support_selection_t::GREALL )
    {
      auto supp_greedy = find_support_greedy(1);
      if( supp_greedy )
        return supp_greedy;
      else
        return find_support_all();
    }
    else if( static_params::support_selection == support_selection_t::GREEDY2 )
    {
      auto supp_greedy = find_support_greedy(1);
      if( supp_greedy )
        return supp_greedy;
      else
      {
        auto pair = find_pair_greedy();
        if( pair )
        {
          return find_support_greedy( 1, *pair );
        }
      }
    }
    else if( static_params::support_selection == support_selection_t::GREEDY3 )
    {
      auto supp_greedy = find_support_greedy(1);
      if( supp_greedy )
        return supp_greedy;

      auto pair = find_pair_greedy();
      if( pair )
      {
        auto supp = find_support_greedy( 1, *pair );
        if( supp )
          return *supp;
      }

      auto triplet = find_triplet_greedy();
      if( triplet )
      {
        auto supp = find_support_greedy( 1, *triplet );
        if( supp )
          return *supp;
      }
      return std::nullopt;
    }
  }

  std::optional<std::vector<uint32_t>> find_support_all()
  {
    auto supps3 = find_support3();
    if( supps3 )
    {
      std::uniform_int_distribution<> distrib(0, supps3->size()-1);
      int idx = distrib(RIGRNG);

      std::set<std::vector<uint32_t>>::iterator it = supps3->begin();
      std::advance(it, idx); 

      return *it;
    }
    return std::nullopt;
  }

  std::optional<std::set<std::vector<uint32_t>>> find_support3()
  {
    std::set<std::vector<uint32_t>> res;

    std::vector<uint32_t> supp;

    _uSPFD.reset();
    uint32_t nEdges = _uSPFD.nEdges;

    for( int i0{0}; i0<scored_divs.size()-1; ++i0 )
    {
      _uSPFD.reset();
      _uSPFD.update( get_div( scored_divs[i0].div ) );

      if( _uSPFD.is_covered() )
      {
        continue;
      }
      for( int i1{i0+1}; i1<scored_divs.size(); ++i1 )
      {
        if( (scored_divs[i1].score + scored_divs[i0].score) <= nEdges )
        {
          _uSPFD.reset();
          _uSPFD.update( get_div( scored_divs[i0].div ) );
          _uSPFD.update( get_div( scored_divs[i1].div ) );

          if( _uSPFD.is_covered() )
          {
            supp={scored_divs[i0].div, scored_divs[i1].div};
            std::sort( supp.begin(), supp.end() );
            res.insert( supp );
            return res;
            continue;
          }
        }

        for( int i2{i1+1}; i2<scored_divs.size(); ++i2 )
        {    
          if( (scored_divs[i2].score + scored_divs[i1].score + scored_divs[i0].score) > 2*nEdges ) break;

          _uSPFD.reset();
          _uSPFD.update( get_div( scored_divs[i0].div ) );
          _uSPFD.update( get_div( scored_divs[i1].div ) );
          _uSPFD.update( get_div( scored_divs[i2].div ) );

          if( _uSPFD.is_covered() )
          {
            supp={scored_divs[i0].div, scored_divs[i1].div, scored_divs[i2].div};
            std::sort( supp.begin(), supp.end() );
            res.insert( supp );
            return res;
            continue;
          }

        } 
      } 
    }
    if( res.size() > 0 )
    {
      return res;
    }
    else
    {
      return std::nullopt;
    }
  }

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

  std::optional<std::vector<uint32_t>> find_from_unbalancing( uint32_t pivot )
  {
    auto supp1p = find_greedy_from_unbalancing( pivot, false, true );
    if( supp1p )
    {
//      for( auto x : *supp1p )
//        std::cout << x << " ";
//      std::cout <<  " : pivot is " << pivot << " A" << std::endl;
      
      return *supp1p;
    }

    auto supp0p = find_greedy_from_unbalancing( pivot, true, true );
    if( supp0p )
    {
//      for( auto x : *supp0p )
//        std::cout << x << " ";
//      std::cout <<  " : pivot is " << pivot << " B" << std::endl;
      return *supp0p;
    }

    auto supp1f = find_greedy_from_unbalancing( pivot, false, false );
    if( supp1f )
    {
//      for( auto x : *supp1f )
//        std::cout << x << " ";
//      std::cout <<  " : pivot is " << pivot << " C" << std::endl;
      return *supp1f;
    }

    auto supp0f = find_greedy_from_unbalancing( pivot, true, false );
    if( supp0f )
    {
//      for( auto x : *supp0f )
//        std::cout << x << " ";
//      std::cout <<  " : pivot is " << pivot << " D" << std::endl;
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
    while( !_uSPFD.is_covered() && supp.size() < static_params::max_support_size )
    {
      best_cost = std::numeric_limits<uint32_t>::max();
      if( _uSPFD.is_saturated() ) break;
      for( uint32_t iCnd{1}; iCnd<divisors.size(); ++iCnd )
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

  std::optional<std::vector<uint32_t>> find_first_K( uint32_t K )
  {
    uint32_t cost, best_cost;
    std::vector<uint32_t> best_candidates;
    std::vector<uint32_t> supp;
    
    _uSPFD.reset();

    /* add recomputation of the support */
    while( supp.size() < static_params::max_support_size )
    {
      best_cost = std::numeric_limits<uint32_t>::max();
      if( _uSPFD.is_saturated() ) break;
      for( uint32_t iCnd{1}; iCnd<divisors.size(); ++iCnd )
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

    return supp;
  }

  std::optional<std::vector<uint32_t>> find_pair_greedy( std::vector<uint32_t> supp0={} )
  {
    uint32_t cost, best_cost;
    std::vector<std::array<uint32_t,2>> best_candidates;
    std::vector<uint32_t> supp=supp0;
    best_cost = std::numeric_limits<uint32_t>::max();

    _uSPFD.reset();
    for( auto x : supp0 )
    {
      _uSPFD.update( get_div( x ) );
    }

    for( uint32_t i{1}; i<divisors.size(); ++i )
    {
      for( uint32_t j{i+1}; j<divisors.size(); ++j )
      {
        cost = _uSPFD.evaluate( get_div(i), get_div(j) );
        if( cost < best_cost )
        {
          best_cost = cost;
          best_candidates = {{i,j}};
        }
        else if( cost == best_cost )
        {
          best_candidates.push_back({i,j});
        }
      }
    }

    if( best_candidates.size() == 0 )
    {
      return std::nullopt;
    }

    std::uniform_int_distribution<> distrib(0, best_candidates.size()-1);
    int idx = distrib(RIGRNG);
    supp.push_back( best_candidates[idx][0] );
    supp.push_back( best_candidates[idx][1] );
    
    return supp;
  }

  std::optional<std::vector<uint32_t>> find_pair_corr( std::vector<uint32_t> supp0={} )
  {
    uint32_t cost, best_cost;
    std::vector<std::array<uint32_t,2>> best_candidates;
    std::vector<uint32_t> supp=supp0;
    best_cost = std::numeric_limits<uint32_t>::max();

    _uSPFD.reset();
    for( auto x : supp0 )
    {
      _uSPFD.update( get_div( x ) );
    }

    for( uint32_t i{1}; i<divisors.size(); ++i )
    {
      for( uint32_t j{i+1}; j<divisors.size(); ++j )
      {
        uint32_t correlation = kitty::count_ones( get_div(i)^get_div(j) );
        correlation = std::max( correlation, on_off_sets[0].num_bits()-correlation );
        cost = _uSPFD.evaluate( get_div(i), get_div(j) )*correlation;
        if( cost < best_cost )
        {
          best_cost = cost;
          best_candidates = {{i,j}};
        }
        else if( cost == best_cost )
        {
          best_candidates.push_back({i,j});
        }
      }
    }

    if( best_candidates.size() == 0 )
    {
      return std::nullopt;
    }

    std::uniform_int_distribution<> distrib(0, best_candidates.size()-1);
    int idx = distrib(RIGRNG);
    supp.push_back( best_candidates[idx][0] );
    supp.push_back( best_candidates[idx][1] );
    
    return supp;
  }

  std::optional<std::vector<uint32_t>> find_triplet_greedy( std::vector<uint32_t> supp0={} )
  {
    uint32_t cost, best_cost;
    std::vector<std::array<uint32_t,3>> best_candidates;
    std::vector<uint32_t> supp=supp0;
    best_cost = std::numeric_limits<uint32_t>::max();

    _uSPFD.reset();
    for( auto x : supp0 )
    {
      _uSPFD.update( get_div( x ) );
    }

    for( uint32_t i{1}; i<divisors.size()-2; ++i )
    {
      for( uint32_t j{i+1}; j<divisors.size()-1; ++j )
      {
        for( uint32_t k{j+1}; k<divisors.size(); ++k )
        {
          cost = _uSPFD.evaluate( get_div(i), get_div(j), get_div(k) );
          if( cost < best_cost )
          {
            best_cost = cost;
            best_candidates = {{i,j,k}};
          }
          else if( cost == best_cost )
          {
            best_candidates.push_back({i,j,k});
          }
        }
      }
    }

    if( best_candidates.size() == 0 )
    {
      return std::nullopt;
    }

    std::uniform_int_distribution<> distrib(0, best_candidates.size()-1);
    int idx = distrib(RIGRNG);
    supp.push_back( best_candidates[idx][0] );
    supp.push_back( best_candidates[idx][1] );
    supp.push_back( best_candidates[idx][2] );
    
    return supp;
  }

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

  u_spfd_manager_t<truth_table_t, 1 << static_params::max_support_size> _uSPFD;


  index_list_t index_list;

  /* positive unate: not overlapping with off-set
     negative unate: not overlapping with on-set */
  std::vector<scored_div> scored_divs;

  stats& st;
}; /* xag_resyn_decompose */

}; /* namespace rils */

} /* namespace mockturtle */




