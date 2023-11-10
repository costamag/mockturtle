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
  \file aig_resyn.hpp
  \brief Resynthesis by recursive decomposition for AIGs or aigs.
  (based on ABC's implementation in `giaResub.c` by Alan Mishchenko)

  \author Siang-Yun Lee
*/

#pragma once

#include "../../../utils/index_list.hpp"
#include "../../../utils/node_map.hpp"
#include "../../../utils/stopwatch.hpp"

#include <abcresub/abcresub.hpp>
#include <fmt/format.h>
#include <kitty/kitty.hpp>

#include <algorithm>
#include <optional>
#include <type_traits>
#include <vector>

#include <random>
#include <limits>

namespace mockturtle
{

namespace spfd
{

std::mt19937 RNG(5);
int SEED{5};

struct aig_resyn_static_params
{
  using base_type = aig_resyn_static_params;

  /*! \brief Maximum number of binate divisors to be considered. */
  static constexpr uint32_t max_binates{ 50u };

  /*! \brief Reserved capacity for divisor truth tables (number of divisors). */
  static constexpr uint32_t reserve{ 200u };

  /*! \brief Whether to consider single XOR gates (i.e., using aigs instead of AIGs). */
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


  /*! \brief Maximum number of support variables. */
  static constexpr uint32_t max_support_size{ 4u };

  static constexpr uint32_t max_num_support_samplings{ 1u };

  using truth_table_storage_type = void;
  using node_type = void;
};

template<class TT>
struct aig_resyn_static_params_default : public aig_resyn_static_params
{
  using truth_table_storage_type = std::vector<TT>;
  using node_type = uint32_t;
  static constexpr bool use_xor = false;
};

template<class Ntk, uint32_t K, uint32_t S>
struct aig_resyn_static_params_for_sim_resub : public aig_resyn_static_params
{
  using truth_table_storage_type = incomplete_node_map<kitty::partial_truth_table, Ntk>;
  using node_type = typename Ntk::node;
  static constexpr bool use_xor = false;
  static constexpr uint32_t max_support_size = K;
  static constexpr uint32_t max_num_support_samplings = S;
};

struct aig_resyn_stats
{
  /*! \brief Time for finding 0-resub and collecting unate literals. */
  stopwatch<>::duration time_unate{ 0 };

  /*! \brief Number of 0-resub optimizations */
  uint32_t num_0resub{ 0 };

  /*! \brief Time for sorting the divisors. */
  stopwatch<>::duration time_sort{ 0 };

  /*! \brief Time for performing . */
  stopwatch<>::duration time_spfd{ 0 };

  void report() const
  {
    fmt::print( "[i]         <aig_resyn>\n" );
    fmt::print( "[i]             0-resub      : {:5d} {:>5.2f} secs\n", num_0resub, to_seconds( time_unate ) );
    fmt::print( "[i]             sort         : {:>5.2f} secs\n", to_seconds( time_sort ) );
  }
};

/*! \brief Logic resynthesis engine for AIGs or aigs.
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
      aig_resyn_stats st;
      aig_resyn<TT, node_map<TT, aig_network>, false, false, aig_network::node> resyn( st );
      auto result = resyn( target, care, divisors.begin(), divisors.end(), tts );
   \endverbatim
 */
template<class TT, class static_params = aig_resyn_static_params_default<TT>>
class aig_resyn
{
public:
  using stats = aig_resyn_stats;
  using index_list_t = large_xag_index_list;
  using truth_table_t = TT;
  using divisor_id_t = uint32_t;

private:
  struct scored_lit
  {
    scored_lit( uint32_t l, uint32_t s )
        : lit( l ), score( s )
    {}

    bool operator==( scored_lit const& other ) const
    {
      return lit == other.lit;
    }

    uint32_t lit;
    uint32_t score;
  };

  template<class LTT, uint32_t CAP>
  struct spfd_manager_t
  {
    spfd_manager_t(){}
    
    void init( LTT const& func, LTT const& careset )
    {
      care = careset;
      on_off_sets[0] = ~func & careset;
      on_off_sets[1] =  func & careset;
      reset();
    }

    void reset()
    {
      masks[0] = care;
      nMasks = 1;
      nEdges = kitty::count_ones( on_off_sets[1] ) * kitty::count_ones( on_off_sets[0] );  
      killed[0] = nEdges > 0 ? false : true;
      nKills = nEdges > 0 ? 0 : 1;
    }

    bool update( LTT const& tt )
    {
      if( is_saturated() )  return false;
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
          masks[nMasks+iMask] = masks[iMask] & ~tt;
          masks[iMask] &= tt;
          if( kitty::count_ones( masks[iMask] & on_off_sets[1] ) == 0 || kitty::count_ones( masks[iMask] & on_off_sets[0] ) == 0 )
          {
            killed[iMask] = true;
            nKills++;
          }
          else
          {
            nEdges += kitty::count_ones( on_off_sets[1] & masks[iMask] ) * kitty::count_ones( on_off_sets[0] & masks[iMask] );  
          }

          if( kitty::count_ones( masks[nMasks+iMask] & on_off_sets[1] ) == 0 || kitty::count_ones( masks[nMasks+iMask] & on_off_sets[0] ) == 0 )
          {
            killed[nMasks+iMask] = true;
            nKills++;
          }
          else
          {
            nEdges += kitty::count_ones( on_off_sets[1] & masks[nMasks+iMask] ) * kitty::count_ones( on_off_sets[0] & masks[nMasks+iMask] );  
          }
        }
      }
      nMasks+=nMasks;
      return true;
    }

    uint32_t evaluate( LTT tt )
    {
      uint32_t res=0;
      for( auto iMask{0}; iMask<nMasks; ++iMask )
      {
        if( !killed[iMask] )
        {
          res+= kitty::count_ones( on_off_sets[1] & masks[iMask] &  tt ) * kitty::count_ones( on_off_sets[0] & masks[iMask] & tt );  
          res+= kitty::count_ones( on_off_sets[1] & masks[iMask] & ~tt ) * kitty::count_ones( on_off_sets[0] & masks[iMask] & ~tt );  
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

    LTT masks[CAP];
    bool killed[CAP];
    uint32_t nMasks;
    uint32_t nKills;
    uint32_t nEdges;
    LTT care;
    LTT on_off_sets[2];
  };

  template<class LTT, uint32_t CAP>
  struct support_generator_t
  {
    support_generator_t(){};

    void init( LTT const& func, LTT const& careset )
    {
      analyzer.init( func, careset );
    }

    void reset()
    {
      past_supports.clear();
    }

    template<class DIV>
    std::optional<std::vector<uint32_t>> generate_support( DIV * pDiv, std::vector<uint32_t> const& candidates )
    {
      if( past_supports.size() == 0 )
        return generate_support_0( pDiv, candidates );
      else
        return generate_support_n( pDiv, candidates );
    }

    template<class DIV>
    std::optional<std::vector<uint32_t>> generate_support_0( DIV * pDiv, std::vector<uint32_t> const& candidates )
    {
      uint32_t cost, best_cost;
      std::vector<uint32_t> best_candidates;

      for( uint32_t i{0}; i<max_num_attempts; ++i )
      {
        analyzer.reset();
        best_candidates.clear();
        support.clear();
        best_cost = analyzer.nEdges;
        while( !analyzer.is_covered() )
        {
          if( analyzer.is_saturated() ) break;
          for( uint32_t iCnd{0}; iCnd<candidates.size(); ++iCnd )
          {
            cost = analyzer.evaluate( pDiv->get_div( candidates[iCnd] ) );
            if( cost < best_cost )
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
          int idx = distrib(RNG);
          support.push_back( candidates[best_candidates[idx]] );
          analyzer.update( pDiv->get_div( candidates[best_candidates[idx]] ) );
        }
        if( analyzer.is_covered() )
        {
          std::sort( support.begin(), support.end() );
          past_supports.insert( support );

          return support;
        }
      }
      return std::nullopt;
    }

    template<class DIV>
    std::optional<std::vector<uint32_t>> generate_support_n( DIV * pDiv, std::vector<uint32_t> const& candidates )
    {
      return std::nullopt;
    }

    uint32_t max_num_attempts{1u};
    std::vector<uint32_t> support;
    spfd_manager_t<LTT, 1u<<CAP> analyzer;
    std::set<std::vector<uint32_t>> past_supports;
  };

public:
  explicit aig_resyn( stats& st ) noexcept
      : st( st )
  {
    static_assert( std::is_same_v<typename static_params::base_type, aig_resyn_static_params>, "Invalid static_params type" );
    static_assert( !( static_params::uniform_div_cost && static_params::preserve_depth ), "If depth is to be preserved, divisor depth cost must be provided (usually not uniform)" );
    divisors.reserve( static_params::reserve );
    divisor_ids.reserve( static_params::reserve );
    lits.reserve( static_params::reserve );
  }

  /*! \brief Perform aig resynthesis.
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

    divisors.resize( 1 ); /* clear previous data and reserve 1 dummy node for constant */
    while ( begin != end )
    {
      if constexpr ( static_params::copy_tts )
      {
        divisors.emplace_back( ( *ptts )[*begin] );
      }
      else
      {
        divisors.emplace_back( *begin );
      }
      ++begin;
    }

    for( uint32_t i{1}; i < divisors.size(); ++i )
    {
      divisor_ids.emplace_back(i);
    }

    _support_generator.init( on_off_sets[1], care );

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

    auto const resS = call_with_stopwatch( st.time_spfd, [&]() {
      return find_resynthesis( num_inserts );
    } );
    if( resS )
    {
      return *resS;
    }

    return std::nullopt;
  }

  /* See if there is a Boolean cut and a candidate resynthesis function
     1. Perform iterative greedy support selection
     2. Extract the local functionality
     3. Resynthesize using Boolean matching with don't cares
   */
  std::optional<uint32_t> find_resynthesis( uint32_t max_num_gates )
  {
    _support_generator.reset();
    for( uint32_t i{0}; i<static_params::max_num_support_samplings; ++i )
    {
      RNG.seed(SEED++);
      auto const supp = _support_generator.generate_support( this, divisor_ids );
      if( supp )
      {
        for( auto x : *supp )
          printf("%d ", x);
        printf("\n");
      }
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
    num_bits[0] = kitty::count_ones( on_off_sets[0] ); /* off-set */
    num_bits[1] = kitty::count_ones( on_off_sets[1] ); /* on-set */
    num_edges = num_bits[0] + num_bits[1]; /* total */

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
  uint32_t num_edges;
  
  const typename static_params::truth_table_storage_type* ptts;
  std::vector<std::conditional_t<static_params::copy_tts, TT, typename static_params::node_type>> divisors;
 
  std::vector<divisor_id_t> divisor_ids;
  support_generator_t<truth_table_t, static_params::max_support_size> _support_generator;

  index_list_t index_list;
  std::vector<scored_lit> lits;

  stats& st;
}; /* aig_resyn */ 


} /* namespace spfd */

} /* namespace mockturtle */