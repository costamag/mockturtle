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
  \file xmg_resyn.hpp
  \brief Resynthesis by recursive decomposition for AIGs or xmgs.
  (based on ABC's implementation in `giaResub.c` by Alan Mishchenko)

  \author Siang-Yun Lee
*/

#pragma once

#include "../../../utils/index_list.hpp"
#include "../../../utils/node_map.hpp"
#include "../../../utils/stopwatch.hpp"
#include "../../../utils/tech_library.hpp"
#include "../../node_resynthesis/xmg_npn.hpp"

#include <abcresub/abcresub.hpp>
#include <fmt/format.h>
#include <kitty/kitty.hpp>
#include <kitty/constructors.hpp>

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

namespace xmg
{

bool VERBOSE{false};

template<class TT>
void print_tt_with_dcs( TT tt, TT mk )
{
  for( auto m = tt.num_bits()-1; m >= 0; --m )
  {
    if( kitty::get_bit( mk, m ) == 1 )
    {
      if( kitty::get_bit( tt, m ) == 1 )
      {
        printf("1");
      }
      else
      {
        printf("0");
      }
    }
    else
    {
      printf("*");
    }
  }
  printf("\n");
}

std::mt19937 RNG(5);

template<class TT> TT compute_buff( TT const& tt1, TT const& tt2, TT const& tt3 ){ return tt1; }
template<class TT> TT compute_m111( TT const& tt1, TT const& tt2, TT const& tt3 ){ return ( tt1 & tt2 ) | ( tt1 & tt3 ) | ( tt2 & tt3 ); }
template<class TT> TT compute_m110( TT const& tt1, TT const& tt2, TT const& tt3 ){ return ( tt1 & tt2 ) | ( tt1 &~tt3 ) | ( tt2 &~tt3 ); }
template<class TT> TT compute_m101( TT const& tt1, TT const& tt2, TT const& tt3 ){ return ( tt1 &~tt2 ) | ( tt1 & tt3 ) | (~tt2 & tt3 ); }
template<class TT> TT compute_m011( TT const& tt1, TT const& tt2, TT const& tt3 ){ return (~tt1 & tt2 ) | (~tt1 & tt3 ) | ( tt2 & tt3 ); }
template<class TT> TT compute_xor3( TT const& tt1, TT const& tt2, TT const& tt3 ){ return tt1 ^ tt2 ^ tt3; }

template<class LIST> uint32_t add_buff_to_list( LIST& list, uint32_t lit1, uint32_t lit2, uint32_t lit3 ){ return lit1; }
template<class LIST> uint32_t add_m111_to_list( LIST& list, uint32_t lit1, uint32_t lit2, uint32_t lit3 ){ return list.add_maj( lit1, lit2, lit3 ); }
template<class LIST> uint32_t add_m110_to_list( LIST& list, uint32_t lit1, uint32_t lit2, uint32_t lit3 ){ return list.add_maj( lit1, lit2, lit3 ^ 0x1 ); }
template<class LIST> uint32_t add_m101_to_list( LIST& list, uint32_t lit1, uint32_t lit2, uint32_t lit3 ){ return list.add_maj( lit1, lit2 ^ 0x1, lit3 ); }
template<class LIST> uint32_t add_m011_to_list( LIST& list, uint32_t lit1, uint32_t lit2, uint32_t lit3 ){ return list.add_maj( lit1 ^ 0x1, lit2, lit3 ); }
template<class LIST> uint32_t add_xor3_to_list( LIST& list, uint32_t lit1, uint32_t lit2, uint32_t lit3 ){ return list.add_xor3( lit1, lit2, lit3 ); }


struct xmg_resyn_static_params
{
  using base_type = xmg_resyn_static_params;

  /*! \brief Maximum number of binate divisors to be considered. */
  static constexpr uint32_t max_binates{ 50u };

  /*! \brief Reserved capacity for divisor truth tables (number of divisors). */
  static constexpr uint32_t reserve{ 200u };

  /*! \brief Whether to consider single XOR gates (i.e., using xmgs instead of AIGs). */
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
  static constexpr uint32_t max_num_support_samplings{ 20u };
  static constexpr uint32_t max_resyn_attempts{ 1 };
  static constexpr double beta_support{100};
  static constexpr double beta_synthesis{10000};

  static constexpr bool try_boolean_matching{false};
  static constexpr bool use_greedy_support{ false };
  static constexpr bool use_local_search{ true };

  using truth_table_storage_type = void;
  using node_type = void;
};

template<class TT>
struct xmg_resyn_static_params_default : public xmg_resyn_static_params
{
  using truth_table_storage_type = std::vector<TT>;
  using node_type = uint32_t;
  static constexpr bool use_xor = false;
};

template<class Ntk, uint32_t SUPP_SIZE, uint32_t N_SAMPL, uint32_t N_RESYN, bool IS_BMATCH, bool IS_GREEDY, bool IS_LSEARCH>
struct xmg_resyn_static_params_for_sim_resub : public xmg_resyn_static_params
{
  using truth_table_storage_type = incomplete_node_map<kitty::partial_truth_table, Ntk>;
  using node_type = typename Ntk::node;
  static constexpr bool use_xor = false;
  static constexpr uint32_t max_support_size = SUPP_SIZE;
  static constexpr uint32_t max_num_support_samplings = N_SAMPL;
  static constexpr uint32_t max_resyn_attempts = N_RESYN;
  static constexpr uint32_t try_boolean_matching = IS_BMATCH;
  static constexpr uint32_t use_greedy_support = IS_GREEDY;
  static constexpr uint32_t use_local_search = IS_LSEARCH;

};

struct xmg_resyn_stats
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
    fmt::print( "[i]         <xmg_resyn>\n" );
    fmt::print( "[i]             0-resub      : {:5d} {:>5.2f} secs\n", num_0resub, to_seconds( time_unate ) );
    fmt::print( "[i]             sort         : {:>5.2f} secs\n", to_seconds( time_sort ) );
  }
};

/*! \brief Logic resynthesis engine for AIGs or xmgs.
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
      const std::vector<xmg_network::node> divisors = ...;
      const node_map<TT, xmg_network> tts = ...;
      const TT target = ..., care = ...;
      xmg_resyn_stats st;
      xmg_resyn<TT, node_map<TT, xmg_network>, false, false, xmg_network::node> resyn( st );
      auto result = resyn( target, care, divisors.begin(), divisors.end(), tts );
   \endverbatim
 */
template<class TT, class static_params = xmg_resyn_static_params_default<TT>>
class xmg_resyn
{
public:
  using stats = xmg_resyn_stats;
  using index_list_t = xmg_index_list;
  using truth_table_t = TT;
  using truth_table4_t = kitty::static_truth_table<4u>;
  using truth_tableK_t = kitty::static_truth_table<static_params::max_support_size>;
  using divisor_id_t = uint32_t;

private:

  struct divisort_t;
  struct candidate_t;
  template<class LTT, uint32_t CAP> struct spfd_manager_t;
  struct functional_library_t;


  struct gate_t
  {
    gate_t(){};
    gate_t( uint32_t specs, truth_tableK_t (*pF)( truth_tableK_t const&, truth_tableK_t const&, truth_tableK_t const& ), uint32_t (*pG)( index_list_t&, uint32_t, uint32_t, uint32_t ) ) : specs(specs), pF(pF), pG(pG){}

    truth_tableK_t compute( truth_tableK_t const& tt1, truth_tableK_t const& tt2, truth_tableK_t const& tt3 ){ return pF( tt1, tt2, tt3 ); };
    truth_tableK_t compute( truth_tableK_t const& tt1 ){ return pF( tt1, tt1, tt1 ); };

    uint32_t add_to_list( index_list_t& list, uint32_t lit1, uint32_t lit2, uint32_t lit3 ){ return pG( list, lit1, lit2, lit3 ); }
    uint32_t add_to_list( index_list_t& list, uint32_t lit1 ){ return pG( list, lit1, lit1, lit1 ); }

    bool is_buffer()
    {
      return specs & 0x1 == 0x1;
    }

    uint32_t specs;
    truth_tableK_t (*pF)( truth_tableK_t const&, truth_tableK_t const&, truth_tableK_t const& );
    uint32_t (*pG)( index_list_t&, uint32_t, uint32_t, uint32_t );
  };

  struct functional_library_t
  {
    functional_library_t()
    {
      gates1[0] = gate_t{ 0x0, &compute_buff<truth_tableK_t>, &add_buff_to_list<index_list_t> };
      gates3[0] = gate_t{ 0x7, &compute_m111<truth_tableK_t>, &add_m111_to_list<index_list_t> };
      gates3[1] = gate_t{ 0x6, &compute_m110<truth_tableK_t>, &add_m110_to_list<index_list_t> };
      gates3[2] = gate_t{ 0x5, &compute_m101<truth_tableK_t>, &add_m101_to_list<index_list_t> };
      gates3[3] = gate_t{ 0x3, &compute_m011<truth_tableK_t>, &add_m011_to_list<index_list_t> };
      gates3[4] = gate_t{ 0x8, &compute_xor3<truth_tableK_t>, &add_xor3_to_list<index_list_t> };
    }

    std::array<gate_t, 1u> gates1;
    std::array<gate_t, 5u> gates3;
  };

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

  struct divisor_t
  {
    divisor_t( truth_tableK_t func, uint32_t lit ) : func(func), lit(lit) {}
    divisor_t( truth_tableK_t func ) : func( func ) {}
    divisor_t(){}

    truth_tableK_t func;
    uint32_t lit;
  };

  struct divisors_t
  {
    divisors_t(){}

    void emplace_back( truth_tableK_t func, uint32_t lit )
    {
      divs.emplace_back( func, lit );
    }

    uint32_t size()
    {
      return divs.size();
    }

    inline divisor_t const& operator[]( uint32_t idx ) const
    {
      return divs[idx];
    }

    inline truth_tableK_t const& get_div( uint32_t idx ) const
    {
      return divs[idx].func;
    } 

    void set_support( std::vector<uint32_t> const& supp, std::array<truth_tableK_t,static_params::max_support_size> const& funcs )
    {
      divs.clear();
      divs.emplace_back( funcs[0].construct(), 0 );

      for( auto i{0}; i<supp.size(); ++i )
      {
        divs.emplace_back( funcs[i], supp[i] << 1u );
      }
    }

    void set_target( truth_tableK_t const& func, truth_tableK_t const& care )
    {
      spfd.init( func, care );
    }

    void clear()
    {
      divs.clear();
      spfd.reset();
    }

    bool update( index_list_t & list, functional_library_t const& functional_library, uint32_t max_num_gates )
    {
      uint32_t num_buffers{0};
      std::vector<divisor_t> new_divs = {divs[0]};

      std::vector<candidate_t> candidates;
      uint32_t cand_id{0};
      for( uint32_t v1{0}; v1 < divs.size(); ++v1 )
      {
        for( auto gate : functional_library.gates1 )
        {
          if( v1 == 0 ) continue;
          candidates.emplace_back( cand_id++, gate, divs[v1] );
        }

        for( uint32_t v2 = v1 + 1; v2 < divs.size(); ++v2 )
        {
          for( uint32_t v3 = v2 + 1; v3 < divs.size(); ++v3 )
          {
            for( auto gate : functional_library.gates3 )
            {
              candidates.emplace_back( cand_id++, gate, divs[v1], divs[v2], divs[v3] );
            }
          }
        }
      }

      double cost;
      double min_cost = std::numeric_limits<double>::max();
      double max_cost = std::numeric_limits<double>::min();
      std::set<uint32_t> set_used;
      std::uniform_real_distribution<> U01(0,1);

      spfd.reset();

      while( !spfd.is_covered() && new_divs.size() < static_params::max_support_size+1 )
      {
        for( auto & cand : candidates )
        {
          cost = spfd.evaluate( cand.compute() );
          cand.cost = cost;
          if( cost < min_cost && set_used.find(cand.id) == set_used.end() ) min_cost = cost;
          if( cost > max_cost && set_used.find(cand.id) == set_used.end() ) max_cost = cost;
        }

        double Z{0};
        bool copy_previous;
        for( auto & cand : candidates )
        {
          copy_previous = set_used.find(cand.id) != set_used.end();
          copy_previous |= (cand.gate.is_buffer() && ( num_buffers >= divs.size()-2 ));
          Z = cand.update_cost( Z, min_cost, max_cost, copy_previous );
        }

        double rnd = U01(RNG);
        bool is_updated{false};

        for( auto & cand : candidates )
        {
          if( rnd*Z <= cand.cost )
          {
            set_used.insert( cand.id );
            if( cand.gate.is_buffer() ) num_buffers++;

            truth_tableK_t tt = cand.compute();
            new_divs.emplace_back( tt, cand.add_to_list( list ) );
            spfd.update( tt );
            is_updated = true;
            break;
          }
        }
        if( !is_updated )
        {
          return false;
        }
      }
      if( spfd.is_covered() )
      {
        divs = new_divs;
        return true;
      }

      return false;
    }

    std::vector<divisor_t> divs;
    spfd_manager_t<truth_tableK_t, 1<<static_params::max_support_size> spfd;

  };

  template<class LTT, uint32_t CAP>
  struct spfd_manager_t
  {
    spfd_manager_t(){}
    
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

    double evaluate( LTT const& tt )
    {
      double res=0;
      for( auto iMask{0}; iMask<nMasks; ++iMask )
      {
        if( !killed[iMask] )
        {
          res+= kitty::count_ones( func[1] & masks[iMask] &  tt ) * kitty::count_ones( func[0] & masks[iMask] & tt )/nEdges;  
          res+= kitty::count_ones( func[1] & masks[iMask] & ~tt ) * kitty::count_ones( func[0] & masks[iMask] & ~tt )/nEdges;  
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
    double nEdges;
    LTT care;
    std::array<LTT, 2u> func;
  };

  struct candidate_t
  {
    candidate_t(){}
    candidate_t( uint32_t id, gate_t gate, divisor_t const& div1 ) : id(id), gate(gate), div1(div1), div2(div1), div3(div1){}
    candidate_t( uint32_t id, gate_t gate, divisor_t const& div1, divisor_t const& div2, divisor_t const& div3 ) : id(id), gate(gate), div1(div1), div2(div2), div3(div3){}

    uint32_t add_to_list( index_list_t & list ){ return gate.add_to_list( list, div1.lit, div2.lit, div3.lit ); }

    truth_tableK_t compute(){ return gate.compute( div1.func, div2.func, div3.func ); }

    double update_cost( double const& cost_previous, double const& min_cost, double const& max_cost, bool copy_previous )
    {
      if( copy_previous )
      {
        cost = cost_previous;
      }
      else
      {
        cost = cost_previous + exp( -static_params::beta_synthesis*( cost - min_cost )/( max_cost - min_cost ) );
      }
      return cost;
    }

    uint32_t id;
    gate_t gate;
    double cost;
    divisor_t const& div1;
    divisor_t const& div2;
    divisor_t const& div3;
  };

public:
  explicit xmg_resyn( stats& st ) noexcept
      : st( st ), _database(_resyn, {})
  {
    static_assert( std::is_same_v<typename static_params::base_type, xmg_resyn_static_params>, "Invalid static_params type" );
    static_assert( !( static_params::uniform_div_cost && static_params::preserve_depth ), "If depth is to be preserved, divisor depth cost must be provided (usually not uniform)" );
    divisors.reserve( static_params::reserve );
    lits.reserve( static_params::reserve );
    _costs.reserve( static_params::reserve );
  }

  /*! \brief Perform xmg resynthesis.
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

    _gSPFD.init( target, care );

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

    for( auto div : divisors )
    {
      _costs.push_back(0);
    }

    for( uint32_t i{0}; i<4; ++i )
      kitty::create_nth_var( _xs4[i], i );

    for( uint32_t i{0}; i<static_params::max_support_size; ++i )
      kitty::create_nth_var( _xsK[i], i );

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
    _past_supports.clear();
  bool a=false;
    for( uint32_t i{0}; i<static_params::max_num_support_samplings; ++i )
    {
      RNG.seed(i++);
      const auto supp = find_support();

      if( supp )
      {
        if(VERBOSE)
        {
          for( auto s : *supp )
            printf("%d ", s);
          printf(" [%d]\n", max_num_gates );
        }

        a = true;

        if( static_params::try_boolean_matching )
        {
          if( supp->size() > 4u )
          {
            index_list_t index_list_copy;
            for( int iTry{0}; iTry<static_params::max_resyn_attempts; ++iTry )
            {
              index_list_copy=index_list;
              auto [funcK, careK] = extract_functionalityK_from_signatures( *supp );
              if( find_spfd_remapping( *supp, funcK, careK, max_num_gates ) )
              {
                auto [lits4, func4, care4] = extract_functionality4_from_Kdivs( funcK, careK );
                const auto res = find_boolean_matching( lits4, func4, care4, max_num_gates ); 
                if( res )
                {
                  return *res;
                }
              }
            }
          }
          else
          {
            auto [func4, care4] = extract_functionality4_from_signatures( *supp );
            std::array<uint32_t, 4> lits = compute_literals( *supp );
            const auto res =  find_boolean_matching( lits, func4, care4, max_num_gates ); 
            if( res )
            {
              return *res;
            }
          }
          return std::nullopt;
        }
        else       
        {
          if( supp->size() == 0u )  return std::nullopt;


          auto [funcK, careK] = extract_functionalityK_from_signatures( *supp );
          //print_tt_with_dcs( funcK, careK );

          const auto res = find_spfd_resynthesis( *supp, funcK, careK, max_num_gates );
          if( res )
          {
            return *res;
          }
        }
      }
    }
    return std::nullopt;
  }

  bool find_spfd_remapping( std::vector<uint32_t> supp, truth_tableK_t const& funcK, truth_tableK_t const& careK, uint32_t max_num_gates )
  {
    _divsK.clear();
    _divsK.set_target( funcK, careK );
    _divsK.set_support( supp, _xsK );

    while( _divsK.size() > 5 && index_list.num_gates() <= max_num_gates )
    {
      if( !_divsK.update( index_list, _functional_library, max_num_gates ) )  return false;
    }
    return _divsK.size() <= 5 ;
    
  }

  std::tuple<std::array<uint32_t, 4>, truth_table4_t, truth_table4_t> extract_functionality4_from_Kdivs( truth_tableK_t const& funcK, truth_tableK_t const& careK )
  {
    if( _divsK.size() > 5 ) printf("[w] divisors size exceeds the limit \n");
    if( _divsK[0].lit != 0 ) printf("[w] first divisor should be zero \n");

    std::array<uint32_t, 4> lits {0};
    for( int i{1}; i<_divsK.size(); ++i )
      lits[i-1] = _divsK[i].lit;
    
    truth_table4_t func4;
    truth_table4_t care4;
    truth_table4_t temp4;
    truth_tableK_t temp = _divsK[0].func.construct();

    for( uint32_t m{0u}; m < 16u; ++m )
    {
      if( m < ( 1 << (_divsK.size()-1) ) )
      {
        temp = temp | ~temp;
        temp4 = temp4 | ~temp4;

        for( uint32_t l{0u}; l < (_divsK.size()-1); ++l )
        {
          if( ( m >> l ) & 0x1 == 0x1 )
          {
            temp &= _divsK[l+1].func;
            temp4 &= _xs4[l];
          }
          else
          {
            temp &= ~_divsK[l+1].func;
            temp4 &= ~_xs4[l];
          }
        }

        if( kitty::count_ones( temp & careK ) > 0 ) // care value
        {
          care4 |= temp4;

          if( kitty::count_ones( temp & funcK ) > 0 )
          {
            func4 |= temp4;
          }
        }
      }
      else
        kitty::clear_bit( care4, m );
    }

    return std::make_tuple( lits, func4, care4 );
  }

  std::array<uint32_t, 4> compute_literals( std::vector<uint32_t> supp )
  {
    std::array<uint32_t, 4> lits {0};
    for( int i{0}; i<supp.size(); ++i )
    {
      lits[i] = supp[i] << 1u;
    }
    return lits;
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

  #pragma region support_sampling

  std::optional<std::vector<uint32_t>> find_support()
  {
    if( _past_supports.size() == 0 || !static_params::use_local_search ) 
    {
      auto supp = static_params::use_greedy_support ? find_support_greedy( {} ) : find_support_boltz( {} );
      if( supp )
      {
        return *supp;
      }
      else
      {
        return std::nullopt;
      }
    }
    else
    {
      std::vector<uint32_t> partial_support = _support;
      // randomly erase an element an element
      std::uniform_int_distribution<> distrib(0, partial_support.size()-1);
      int idx = distrib(RNG);
      uint32_t erased = partial_support[idx];
      partial_support.erase( partial_support.begin() + idx );
      // cover
      return static_params::use_greedy_support ? find_support_greedy( partial_support, erased ) : find_support_boltz( partial_support, erased );
    }
    return std::nullopt;
  }

  std::optional<std::vector<uint32_t>> find_support_greedy( std::vector<uint32_t>const& partial_support, uint32_t erased = std::numeric_limits<uint32_t>::max() )
  {
    double cost, best_cost;
    std::vector<uint32_t> best_candidates;
    std::vector<uint32_t> supp;
    
    _gSPFD.reset();
    for( auto div : partial_support )
    {
      if( _gSPFD.is_saturated() ) break;
      _gSPFD.update( get_div( div ) );
      supp.push_back( div );
    }
    /* add recomputation of the support */

    while( !_gSPFD.is_covered() )
    {
      best_cost = std::numeric_limits<double>::max();
      if( _gSPFD.is_saturated() ) break;
      for( uint32_t iCnd{1}; iCnd<divisors.size(); ++iCnd )
      {
        cost = _gSPFD.evaluate( get_div( iCnd ) );
        if( cost < best_cost && iCnd != erased )
        {
          best_cost = cost;
          best_candidates = {iCnd};
        }
        else if( cost == best_cost && iCnd != erased )
        {
          best_candidates.push_back( iCnd );
        }
      }
      if( best_candidates.size() == 0 ) break;

      std::uniform_int_distribution<> distrib(0, best_candidates.size()-1);
      int idx = distrib(RNG);
      supp.push_back( best_candidates[idx] );
      _gSPFD.update( get_div( best_candidates[idx] ) );
    }
    if( _gSPFD.is_covered() )
    {
      std::sort( supp.begin(), supp.end() );
      if( _past_supports.insert(supp).second )
      {
        _support = supp;
        return supp;
      }
    }
    return std::nullopt;
  }

  std::optional<std::vector<uint32_t>> find_support_boltz( std::vector<uint32_t> const& partial_support, uint32_t erased = std::numeric_limits<uint32_t>::max() )
  {
    _gSPFD.reset();
    double min_cost;
    double max_cost;
    std::vector<double> costs;
    std::vector<uint32_t> supp;
    for( auto div : partial_support )
    {
      if( _gSPFD.is_covered() ) break;
      supp.push_back( div );
      _gSPFD.update( get_div( div ) );
    }
    /* add recomputation of the support */

    int iter{0};
    while( !_gSPFD.is_covered() )
    {
      if( _gSPFD.is_saturated() ) break;
      costs.clear();
      costs.push_back(0);
      min_cost = std::numeric_limits<double>::max();
      max_cost = std::numeric_limits<double>::min();

      for( uint32_t iDiv{1}; iDiv<divisors.size(); ++iDiv )
      {
        costs.push_back(_gSPFD.evaluate( get_div( iDiv ) ));
        if( costs[iDiv] < min_cost ) min_cost = costs[iDiv];
        if( costs[iDiv] > max_cost ) max_cost = costs[iDiv];
      }

      for( uint32_t i{1}; i<costs.size(); ++i )
      {
        costs[i] = exp( -static_params::beta_support*(costs[i]-min_cost)/(max_cost-min_cost) );
      }
      
      for( auto div : supp )
        costs[div] = 0;
      if( erased != std::numeric_limits<uint32_t>::max() )
        costs[erased] = 0;

      for( uint32_t i{1}; i<costs.size(); ++i )
      {
        costs[i] += costs[i-1];
      }

      std::uniform_real_distribution<> distrib(0, 1);
      double rnd = distrib(RNG);

      bool found{false};
      for( uint32_t i{1}; i<costs.size(); ++i )
      {
        if( rnd*costs.back() <= costs[i] )
        {
          supp.push_back(i);
          _gSPFD.update( get_div(i) );
          found = true;
          break;
        }
      }

      if( !found ) return std::nullopt;
    }
    if( _gSPFD.is_covered() )
    {
      std::sort( supp.begin(), supp.end() );
      if( _past_supports.insert(supp).second )
      {
        _support = supp;
        return supp;
      }
    }

    return std::nullopt;
  }

  #pragma endregion support_sampling

  #pragma region function_extraction

  std::tuple<truth_table4_t,truth_table4_t> extract_functionality4_from_signatures( std::vector<uint32_t> const& supp )
  {
    truth_table4_t func4;
    truth_table4_t care4;
    truth_table4_t temp4;
    truth_table_t  temp = _gSPFD.care.construct();

    for( uint32_t m{0u}; m < 16u; ++m )
    {
      if( m < ( 1 << supp.size() ) )
      {
        temp = temp | ~temp;
        temp4 = temp4 | ~temp4;

        for( uint32_t l{0u}; l < supp.size(); ++l )
        {
          if( ( m >> l ) & 0x1 == 0x1 )
          {
            temp &= get_div(supp[l]);
            temp4 &= _xs4[l];
          }
          else
          {
            temp &= ~get_div(supp[l]);
            temp4 &= ~_xs4[l];
          }
        }

        if( kitty::count_ones( temp & _gSPFD.care ) > 0 ) // care value
        {
          care4 |= temp4;

          if( kitty::count_ones( temp & _gSPFD.func[1] ) > 0 )
          {
            func4 |= temp4;
          }
        }
      }
      else
        kitty::clear_bit( care4, m );
    }

    return { func4, care4 };
  }

  std::tuple<truth_tableK_t,truth_tableK_t> extract_functionalityK_from_signatures( std::vector<uint32_t> const& supp )
  {
    truth_tableK_t funcK;
    truth_tableK_t careK;
    truth_tableK_t tempK;
    truth_table_t  temp = _gSPFD.care.construct();

    for( uint32_t m{0u}; m < (1u << static_params::max_support_size); ++m )
    {
      if( m < ( 1 << supp.size() ) )
      {
        temp = temp | ~temp;
        tempK = tempK | ~tempK;

        for( uint32_t l{0u}; l < supp.size(); ++l )
        {
          if( ( m >> l ) & 0x1 == 0x1 )
          {
            temp &= get_div(supp[l]);
            tempK &= _xsK[l];
          }
          else
          {
            temp &= ~get_div(supp[l]);
            tempK &= ~_xsK[l];
          }
        }

        if( kitty::count_ones( temp & _gSPFD.care ) > 0 ) // care value
        {
          careK |= tempK;

          if( kitty::count_ones( temp & _gSPFD.func[1] ) > 0 )
          {
            funcK |= tempK;
          }
        }
      }
      else
        kitty::clear_bit( careK, m );
    }

    return { funcK, careK };
  }

  #pragma endregion function_extraction

  #pragma region boolean_matching_resynthesis

  std::optional<uint32_t> find_boolean_matching( std::array<uint32_t, 4> lits, truth_table4_t const& func4, truth_table4_t const& care4, uint32_t max_num_gates )
  {
    if(VERBOSE)
    {
      printf("TT(0):"); print_tt_with_dcs(func4, care4);
    }

    auto [func_npn, neg, perm] = exact_npn_canonization( func4 );
    if(VERBOSE)
    {
      printf("neg  = ");
      for( auto iBit=3; iBit>=0; iBit-- )
        printf("%d", (neg >> iBit) & 0x1);

      printf(" | perm  =");
      for( auto iBit=0; iBit<4; iBit++ )
        printf("%d ", perm[iBit]);
      printf("\n");
      
      for( auto iBit=0; iBit<4; iBit++ )
      {
        if( neg >> iBit & 0x1 == 0x1 )
          printf("%2d : ~X[%d] <= X[%d]  <<  X[%d] <= P[%d]\n", lits[iBit] ^ 0x1, iBit, iBit, perm[iBit], iBit );
        else
          printf("%2d :  X[%d] <= X[%d]  <<  X[%d] <= P[%d]\n", lits[iBit], iBit, iBit, perm[iBit], iBit );
      }
    }
    
    auto const care_npn = apply_npn_transformation( care4, neg & ~( 1 << 4u ), perm );
    if(VERBOSE)
    {
      printf("npn(TT)"); print_tt_with_dcs(func_npn, care_npn);
    }


    auto const structures = _database.get_supergates( func_npn, ~care_npn, neg, perm );
    if( structures == nullptr )
    {
      printf("[w] no structure ");
      return std::nullopt;
    }
    if(VERBOSE)
    {
      printf("neg* = ");
      for( auto iBit=3; iBit>=0; iBit-- )
        printf("%d", (neg >> iBit) & 0x1);

      printf(" | perm* =");
      for( auto iBit=0; iBit<4; iBit++ )
        printf("%d ", perm[iBit]);
      printf("\n");
    }
    bool phase = ( neg >> 4 == 1 ) ? true : false;

    for( auto i{0}; i<lits.size(); ++i )
    {
      if( ( neg >> i ) & 0x1 == 0x1 )
        lits[i] ^= 0x1;
    }

    std::array<uint32_t, 4> leaves {0};

    for( auto i{0}; i<4; ++i )
    {
      leaves[i] = lits[perm[i]];
    }


    auto & db = _database.get_database();
    db.incr_trav_id();

    std::unordered_map<uint64_t, uint32_t> existing_nodes;

    index_list_t index_list_copy = index_list;

    auto res = create_index_list( db.get_node( structures->at(0).root ), leaves, existing_nodes );

    //if( res )
    {
      if(VERBOSE)
      {
        printf(" || --> [%d <?= %d]\n", index_list.num_gates(), max_num_gates );
      }
      //printf("%d ", index_list.num_gates() );

      if( index_list.num_gates() <= max_num_gates )
      {
        return phase ? res ^ 0x1 : res;//create_index
      }
      else
        index_list = index_list_copy;
    }

    return std::nullopt;
  }

  uint32_t create_index_list( node<xmg_network> const& n, std::array<uint32_t, 4u> const& leaves, std::unordered_map<uint64_t, uint32_t> & existing_nodes )
  {
    auto new_lit = create_index_list_rec( n, leaves, existing_nodes );

    //if( new_lit )
    {
      return new_lit;
    }
    //return std::nullopt;
  }

  uint32_t create_index_list_rec( node<xmg_network> const& n, std::array<uint32_t, 4u> const& leaves, std::unordered_map<uint64_t, uint32_t> & existing_nodes )
  {
    auto& db = _database.get_database();

    std::array<uint32_t, 3> node_data;
    db.foreach_fanin( n, [&]( auto const& f, auto i ) 
    {
      node<xmg_network> g = db.get_node( f );
      if( db.is_pi( g ) )
      {
        node_data[i] = db.is_complemented( f ) ? leaves[f.index-1] ^ 0x1 : leaves[f.index-1];
        return;
      }
      else if( db.is_constant(g) )
      {
        node_data[i] = db.is_complemented( f ) ? 0x1 : 0x0;
        return;
      }
      else if( db.is_maj( g ) )
      {
        auto res = create_index_list_rec( g, leaves, existing_nodes );
        node_data[i] = db.is_complemented( f ) ? res ^ 0x1 : res;
      }
      else if( db.is_xor3( g ) )
      {
        auto res = create_index_list_rec( g, leaves, existing_nodes );
        node_data[i] = db.is_complemented( f ) ? res ^ 0x1 : res;
      }
    });

    if( db.is_maj(n) )
    {
      uint64_t key = get_key( node_data );
      uint32_t new_lit;
      if( auto search = existing_nodes.find( key ); search != existing_nodes.end() )
      {
        new_lit = search->second;
      }
      else
      {
        new_lit = index_list.add_maj( node_data[0], node_data[1], node_data[2] );
        existing_nodes[key] = new_lit;
      }
      return new_lit;
    }
    else if( db.is_xor3(n) )
    {
      uint64_t key = get_key_xor( node_data );
      uint32_t new_lit;
      if( auto search = existing_nodes.find( key ); search != existing_nodes.end() )
      {
        new_lit = search->second;
      }
      else
      {
        new_lit = index_list.add_xor3( node_data[0], node_data[1], node_data[2] );
        existing_nodes[key] = new_lit;
      }
      return new_lit;
    }
    else
    {
      printf("unknown recursion node\n");
      return 0;
    }
  }

  uint64_t get_key( std::array<uint32_t, 3> node_data )
  {
    std::vector<uint64_t> keys = { (uint64_t)node_data[0], (uint64_t)node_data[1], (uint64_t)node_data[2] };
    std::sort( keys.begin(), keys.end() );

    return keys[0] | ( keys[1] << 20u ) | ( keys[2] << 40u );
  }

  uint64_t get_key_xor( std::array<uint32_t, 3> node_data )
  {
    return get_key( node_data ) | ( 1u << 60u );
  }

  std::optional<uint32_t> create_index_list_rec_old( node<xmg_network> const& n, std::array<uint32_t, 4u> const& leaves, std::unordered_map<uint64_t, uint32_t> & existing_nodes )
  {
    auto& db = _database.get_database();

    std::array<uint32_t, 3u> node_data{0};

    db.foreach_fanin( n, [&]( auto const& f, auto i ) 
    {
      node<xmg_network> g = db.get_node( f );
      if( db.is_pi( g ) )
      {
        node_data[i] = db.is_complemented( f ) ? leaves[f.index-1] ^ 0x1 : leaves[f.index-1];
      }
      else if( db.is_maj(g) )
      {
        auto res = create_index_list_rec( g, leaves, existing_nodes );

        if( res )
        {
          node_data[i] = db.is_complemented( f ) ? (*res) ^ 0x1 : *res;
        }
        else
          return std::nullopt;
      }
      else if( db.is_constant( g ) )
      {
        node_data[i] = db.is_complemented( f ) ? 0x1 : 0x0;
      }
    } );

    if( db.is_maj( n ) )
    {
      uint32_t new_lit;
      uint64_t key0 = node_data[0];
      uint64_t key1 = node_data[1];
      uint64_t key2 = node_data[2];


      if( key0 < key1 && key1 < key2 ) // 0 1 2 
      {
        key0 = ( key0 << 40u ) | ( key1 << 20u ) | key2;
      }
      else if( key0 < key2 && key2 < key1 ) // 0 2 1
      {
        key0 = ( key0 << 40u ) | ( key2 << 20u ) | key1;
      }
      else if( key1 < key0 && key0 < key2 ) // 1 0 2
      {
        key0 = ( key1 << 40u ) | ( key0 << 20u ) | key2;
      }
      else if( key1 < key2 && key2 < key0 ) // 1 2 0
      {
        key0 = ( key1 << 40u ) | ( key2 << 20u ) | key0;
      }
      else if( key2 < key0 && key0 < key1 ) // 2 0 1
      {
        key0 = ( key2 << 40u ) | ( key0 << 20u ) | key1;
      }
      else if( key2 < key1 && key1 < key0 ) // 2 1 0
      {
        key0 = ( key2 << 40u ) | ( key1 << 20u ) | key0;
      }

      if( auto search = existing_nodes.find( key0 ); search != existing_nodes.end() )
      {
        new_lit = search->second;
      }
      else
      {
        new_lit = index_list.add_maj( node_data[0], node_data[1], node_data[2] );
        existing_nodes[key0] = new_lit;
      }
      return new_lit;
    }
    else
    {
      return std::nullopt;
    }
  }

  #pragma endregion boolean_matching_resynthesis

  #pragma region spfd_resynthesis
  std::optional<uint32_t> find_spfd_resynthesis( std::vector<uint32_t> const& supp, truth_tableK_t const& funcK, truth_tableK_t const& careK, uint32_t max_num_gates )
  {
    index_list_t index_list_copy = index_list;
    uint32_t max_num_gates_copy = max_num_gates;
    _divsK.set_target( funcK, careK );

    for( auto iTry{0}; iTry<static_params::max_resyn_attempts; ++iTry )
    {
      index_list = index_list_copy;
      max_num_gates = max_num_gates_copy;

      _divsK.set_support( supp, _xsK );
      while( _divsK.size() > 2 && index_list.num_gates() <= max_num_gates )
      {
        if( !_divsK.update( index_list, _functional_library, max_num_gates ) )  break;
      }
      if( _divsK.spfd.is_covered() && _divsK.size() == 2 )
      {
        if( kitty::equal( _divsK.get_div(1) & _divsK.spfd.care, _divsK.spfd.func[1] ) )
        {
          return _divsK[1].lit;
        }
        else if( kitty::equal( _divsK.get_div(1) & _divsK.spfd.care, _divsK.spfd.func[0] ) )
        {
          return _divsK[1].lit ^ 0x1;
        }
        else
        {
          printf( "[w]: one divisor not matching\n" );
        }
      }
    }
    return std::nullopt;
  }
  #pragma endregion spfd_resynthesis

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
  std::vector<double> _costs;
 
  index_list_t index_list;
  std::vector<scored_lit> lits;

  spfd_manager_t<truth_table_t, 1<<static_params::max_support_size> _gSPFD;
  std::array<truth_table4_t,4> _xs4;
  std::array<truth_tableK_t,static_params::max_support_size> _xsK;
  std::set<std::vector<uint32_t>> _past_supports;
  std::vector<uint32_t> _support;
  divisors_t _divsK;

  functional_library_t _functional_library;

  xmg_npn_resynthesis _resyn;
  exact_library<xmg_network, xmg_npn_resynthesis> _database;


  stats& st;
}; /* xmg_resyn */ 

} /* xmg */

} /* namespace spfd */

} /* namespace mockturtle */