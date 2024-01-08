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
  ALL
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


template<class Ntk>
struct rig_resyn_static_params_for_sim_resub : public rig_resyn_static_params
{
  using truth_table_storage_type = incomplete_node_map<kitty::partial_truth_table, Ntk>;
  using node_type = typename Ntk::node;
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
template<class TT, class static_params = rig_resyn_static_params_default<TT>>
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
      func[1] =  target & careset;
      func[0] = ~target & careset;
      reset();
    }

    void reset()
    {
      masks[0] = care;
      nMasks = 1;
      nEdges = kitty::count_ones( func[1] ) * kitty::count_ones( func[0] );  
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
      return find_support_greedy();
    }
    else if( static_params::support_selection == support_selection_t::ALL )
    {
      return find_support_all();
    }
  }

  std::optional<std::vector<uint32_t>> find_support_all()
  {
    auto supp3 = find_support3();
    if( supp3 )
      return *supp3;
    return std::nullopt;
  }

  std::optional<std::vector<uint32_t>> find_support3()
  {
    std::vector<uint32_t> supp;
    _uSPFD.reset();

    std::array<TT,2u> masks0;
    std::array<TT,4u> masks1;
    std::array<TT,8u> masks2;
    std::array<bool,8> is_killed;
    uint32_t nKills0;
    uint32_t nKills1;
    uint32_t nKills2;

    bool is_valid;
    for( int i0{0}; i0<scored_divs.size()-1; ++i0 )
    {
      supp.clear();
      nKills0 = 0;
      masks0[0] = _uSPFD.care & get_div( scored_divs[i0].div );
      masks0[1] = _uSPFD.care & ~get_div( scored_divs[i0].div );

      if( kitty::is_const0(masks0[0]) || kitty::equal( masks0[0] & _uSPFD.func[1], masks0[0] ) ) { is_killed[0] = true; nKills0++; } else { is_killed[0] = false; }
      if( kitty::is_const0(masks0[1]) || kitty::equal( masks0[1] & _uSPFD.func[1], masks0[1] ) ) { is_killed[1] = true; nKills0++; } else { is_killed[1] = false; }
      if( nKills0 == 2u )
      {
        continue;
      }
      for( int i1{i0+1}; i1<scored_divs.size(); ++i1 )
      {
        nKills1=0;
        for( int k1{0}; k1<2; k1++ )
        {
          masks1[k1] = masks0[k1] &  get_div( scored_divs[i1].div );
          masks1[k1+2] = masks0[k1] & ~get_div( scored_divs[i1].div );
          if( is_killed[k1] )
          { 
            is_killed[k1+2] = true; 
            nKills1+=2; 
          } 
          else 
          { 
            if( kitty::is_const0(masks1[k1]) || kitty::equal( masks1[k1] & _uSPFD.func[1], masks1[k1] ) ) { is_killed[k1] = true; nKills1++; } else { is_killed[k1] = false; }
            if( kitty::is_const0(masks1[k1+2]) || kitty::equal( masks1[k1+2] & _uSPFD.func[1], masks1[k1+2] ) ) { is_killed[k1+2] = true; nKills1++; } else { is_killed[k1+2] = false; }
          }
        }
        if( nKills1 == 4u )
        {
          supp={scored_divs[i0].div, scored_divs[i1].div};
          std::sort( supp.begin(), supp.end() );
          return supp;
        }

        for( int i2{i1+1}; i2<scored_divs.size(); ++i2 )
        {    
          nKills2=0;
          if( scored_divs[i2].score + scored_divs[i1].score + scored_divs[i0].score > _uSPFD.nEdges ) break;

          uint32_t count;

          for( int k2{0}; k2<4; k2++ )
          {
            masks2[k2] = masks1[k2] &  get_div( scored_divs[i2].div );
            masks2[k2+4] = masks1[k2] & ~get_div( scored_divs[i2].div );
            if( is_killed[k2] )
            { 
              is_killed[k2+4] = true; 
              nKills2+=2; 
            } 
            else 
            { 
              if( kitty::is_const0(masks2[k2]) || kitty::equal( masks2[k2] & _uSPFD.func[1], masks2[k2] ) ) { is_killed[k2] = true; nKills2++; } else { is_killed[k2] = false; }
              if( kitty::is_const0(masks2[k2+4]) || kitty::equal( masks2[k2+4] & _uSPFD.func[1], masks2[k2+4] ) ) { is_killed[k2+4] = true; nKills2++; } else { is_killed[k2+4] = false; }
            }
          }
          if( nKills2 == 8u )
          {
            supp={scored_divs[i0].div, scored_divs[i1].div, scored_divs[i2].div};
            std::sort( supp.begin(), supp.end() );
            return supp;
          }
        } 
      } 
    }
    return std::nullopt;
  }

  std::optional<std::vector<uint32_t>> find_support_greedy()
  {
    uint32_t cost, best_cost;
    std::vector<uint32_t> best_candidates;
    std::vector<uint32_t> supp;
    
    _uSPFD.reset();

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
      std::sort( supp.begin(), supp.end() );
      return supp;
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