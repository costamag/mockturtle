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
  \brief Resynthesis for AIGs or XAGs.
  The engine based on decomposition is based on ABC's implementation in `giaResub.c` by Alan Mishchenko

  \author Andrea Costamagna
*/

#pragma once

#include "../../utils/index_list.hpp"
#include "../../utils/node_map.hpp"
#include "../../utils/stopwatch.hpp"
#include "../../utils/tech_library.hpp"
#include "../node_resynthesis/xag_npn.hpp"

#include <abcresub/abcresub.hpp>
#include <fmt/format.h>
#include <kitty/kitty.hpp>

#include <algorithm>
#include <optional>
#include <type_traits>
#include <vector>

#include <thread>

namespace mockturtle
{

namespace bmatch
{

std::mt19937 RNG( 5 );

enum gate_t
{
  PA00,
  PA01,
  PA10,
  PA11,
  IA00,
  IA01,
  IA10,
  IA11,
  INV_,
  BUF_,
  EXOR,
  NONE
};


template<class DTT> DTT hpcompute_buf_( DTT const& a, DTT const& b ){ return  a; }
template<class DTT> DTT hpcompute_pa00( DTT const& a, DTT const& b ){ return ~a & ~b; }
template<class DTT> DTT hpcompute_pa01( DTT const& a, DTT const& b ){ return ~a &  b; }
template<class DTT> DTT hpcompute_pa10( DTT const& a, DTT const& b ){ return  a & ~b; }
template<class DTT> DTT hpcompute_pa11( DTT const& a, DTT const& b ){ return  a &  b; }
template<class DTT> DTT hpcompute_exor( DTT const& a, DTT const& b ){ return  a ^  b; }


template<class LIST> uint32_t add_buf__to_index_list_( LIST& list, uint32_t lit1, uint32_t lit2 ){ return  lit1; }
template<class LIST> uint32_t add_pa00_to_index_list_( LIST& list, uint32_t lit1, uint32_t lit2 ){ return list.add_and( lit1 ^ 0x1, lit2 ^ 0x1 ); }
template<class LIST> uint32_t add_pa01_to_index_list_( LIST& list, uint32_t lit1, uint32_t lit2 ){ return list.add_and( lit1 ^ 0x1, lit2 ); }
template<class LIST> uint32_t add_pa10_to_index_list_( LIST& list, uint32_t lit1, uint32_t lit2 ){ return list.add_and( lit1, lit2 ^ 0x1 ); }
template<class LIST> uint32_t add_pa11_to_index_list_( LIST& list, uint32_t lit1, uint32_t lit2 ){ return list.add_and( lit1, lit2 ); }
template<class LIST> uint32_t add_exor_to_index_list_( LIST& list, uint32_t lit1, uint32_t lit2 ){ return list.add_xor( lit1, lit2 ); }


template<class STT, uint32_t CAP>
struct spfd_manager_t
{
  spfd_manager_t(){}

  spfd_manager_t( STT careset, STT func ) : care( careset ), onset( func & careset ), offset( ~func & care )
  {
    reset();
  }

  void init( STT careset, STT func )
  {
    care = careset;
    onset = func & care;
    offset = ~func & care; 
    reset();
  }

  void reset()
  {
    masks[0]=care;
    killed[0]=false;
    nMasks = 1u;
    nEdges = kitty::count_ones( onset ) * kitty::count_ones( offset );
    nKilled = nEdges > 0 ? 0u : 1u;
  }

  void update( STT const& tt )
  {
    nEdges = 0;

    for( uint32_t iMask{0}; iMask < nMasks; ++iMask )
    {
      if( killed[iMask] )
      {
        killed[nMasks+iMask] = true;
        masks[nMasks+iMask] = masks[iMask];
        nKilled++;
      }
      else
      {
        masks[nMasks+iMask] = masks[iMask] & tt;
        killed[nMasks+iMask] = false;

        if( ( kitty::count_ones( onset & masks[nMasks+iMask] ) == 0 ) || ( kitty::count_ones( offset & masks[nMasks+iMask] ) == 0 ) )
        {
          killed[nMasks+iMask] = true;
          nKilled++;
        }
        else
        {
          nEdges += kitty::count_ones( masks[nMasks+iMask] & onset )*kitty::count_ones( masks[nMasks+iMask] & offset );
        }

        masks[iMask] = masks[iMask] & ~tt;
        if( ( kitty::count_ones( onset & masks[iMask] ) == 0 ) || ( kitty::count_ones( offset & masks[iMask] ) == 0 ) )
        {
          killed[iMask] = true;
          nKilled++;
        }
        else
        {
          nEdges += kitty::count_ones( masks[iMask] & onset )*kitty::count_ones( masks[iMask] & offset );
        }
      }
    }
    nMasks=nMasks*2;
  }

  double evaluate( STT const& tt )
  {
    double res = 0;
    for( auto m = 0; m < nMasks; ++m )
    {
      if( !killed[m] )
      {
        res += kitty::count_ones( masks[m] & tt & onset )*kitty::count_ones( masks[m] & tt & offset )/nEdges;
        res += kitty::count_ones( masks[m] & ~tt & onset )*kitty::count_ones( masks[m] & ~tt & offset )/nEdges;
      }
    }
    return res;
  }

  bool is_covered()
  {
    return nMasks <= nKilled;
  }

  /*! \brief original careset */
  STT care;
  STT onset;
  STT offset;

  std::array<STT, CAP> masks;
  std::array<bool, CAP> killed;
  uint32_t nMasks{1};
  uint32_t nKilled{0};
  double nEdges{1};
};

template<class TT>
struct divisor_t
{
  divisor_t( TT func, uint32_t lit ) : func(func), lit(lit) {}
  divisor_t( TT func ) : func(func) {}
  divisor_t(){}
  
  TT func;
  uint32_t lit;
};

template<class TT>
struct divisors_t
{
  divisors_t(){}

  void emplace_back( TT func, uint32_t lit )
  {
    divs.emplace_back(func, lit);
  }

  void clear()
  {
    divs.clear();
  }

  uint32_t size()
  {
    return divs.size();
  }

  inline divisor_t<TT> const& operator[]( uint32_t idx ) const
  {
      return divs[idx];
  }

  inline TT const& get_sign( uint32_t idx ) const
  {
      return divs[idx].func;
  }
 
  std::vector<divisor_t<TT>> divs;
};

struct xag_resyn_static_params
{
  using base_type = xag_resyn_static_params;

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

  /*! \brief Maximum support size */
  static constexpr uint32_t max_support_size{ 7u };
  static constexpr uint32_t max_num_spfds{ 10u };


  /* FOR BOOLEAN MATCHING RESUBSTITUTION */
  /*! \brief recursively decompose */
  static constexpr uint32_t max_support_attempts{ 10u };
  static constexpr uint32_t max_resynthesis_attempts{ 10u };
  static constexpr bool try_0resub{ true };
  static constexpr bool try_1resub{ false };
  static constexpr bool try_unateness_decomposition{ false };
  static constexpr bool use_boolean_matching{ false };
  static constexpr double beta_support{ 100 };
  static constexpr double beta_synthesis{ 100 };
  static constexpr bool use_greedy_support_selection{false};

  using truth_table_storage_type = void;
  using node_type = void;
};

template<class TT>
struct xag_resyn_static_params_default : public xag_resyn_static_params
{
  using truth_table_storage_type = std::vector<TT>;
  using node_type = uint32_t;
};

template<class TT>
struct aig_resyn_static_params_default : public xag_resyn_static_params_default<TT>
{
  static constexpr bool use_xor = false;
};

template<class Ntk, uint32_t K, uint32_t S, uint32_t I>
struct xag_resyn_static_params_for_sim_resub : public xag_resyn_static_params
{
  using truth_table_storage_type = incomplete_node_map<kitty::partial_truth_table, Ntk>;
  using node_type = typename Ntk::node;
  static constexpr uint32_t max_support_size = K;
  static constexpr uint32_t max_support_attempts = S;
  static constexpr uint32_t max_resynthesis_attempts = I;
  static constexpr uint32_t max_num_spfds = max_support_size + 2 ;
};

template<class Ntk, uint32_t K, uint32_t S, uint32_t I>
struct aig_resyn_static_params_for_sim_resub : public xag_resyn_static_params_for_sim_resub<Ntk, K, S, I>
{
  static constexpr bool use_xor = false;
};

template<class Ntk, uint32_t K, uint32_t S, uint32_t I>
struct bmatch_aig_resyn_static_params_for_sim_resub : public aig_resyn_static_params_for_sim_resub<Ntk, K, S, I>
{
  static constexpr bool use_boolean_matching = true;
};

template<class Ntk, uint32_t K, uint32_t S, uint32_t I>
struct bmatch_xag_resyn_static_params_for_sim_resub : public xag_resyn_static_params_for_sim_resub<Ntk, K, S, I>
{
  static constexpr bool use_boolean_matching = true;

};

#pragma region XAG_resyn

struct xag_resyn_stats
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
  stopwatch<>::duration time_boolean_matching{ 0 };
  stopwatch<>::duration time_spfd_synthesis{ 0 };


  void report() const
  {
    fmt::print( "[i]         <xag_resyn>\n" );
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
      xag_resyn<TT, node_map<TT, aig_network>, false, false, aig_network::node> resyn( st );
      auto result = resyn( target, care, divisors.begin(), divisors.end(), tts );
   \endverbatim
 */

template<class TT, class static_params = xag_resyn_static_params_default<TT>>
class xag_resyn
{
public:
  using stats = xag_resyn_stats;
  using index_list_t = large_xag_index_list;
  using truth_table_t = TT;
  using small_truth_table_t = kitty::static_truth_table<static_params::max_support_size>;
  using truth_table4_t = kitty::static_truth_table<4u>;

private:
  struct unate_lit
  {
    unate_lit( uint32_t l )
        : lit( l )
    {}

    bool operator==( unate_lit const& other ) const
    {
      return lit == other.lit;
    }

    uint32_t lit;
    uint32_t score{ 0 };
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

  struct xaig_gate_t
  {
    gate_t type;
    uint32_t nInputs;

    small_truth_table_t (*pF)( small_truth_table_t const&, small_truth_table_t const& );
    uint32_t (*pG)( index_list_t&, uint32_t, uint32_t );

    xaig_gate_t( gate_t type, uint32_t nInputs, small_truth_table_t (*pF)( small_truth_table_t const&, small_truth_table_t const& ), uint32_t (*pG)( index_list_t&, uint32_t, uint32_t ) ) : type(type), nInputs(nInputs), pF(pF), pG(pG){}
    xaig_gate_t(){}

    small_truth_table_t compute( small_truth_table_t const& a, small_truth_table_t const& b ){ return pF( a, b ); };
    small_truth_table_t compute( small_truth_table_t const& a ){ return pF( a, a ); };

    uint32_t add_to_list( index_list_t& list, uint32_t lit1, uint32_t lit2 ){ return pG( list, lit1, lit2 ); };
    uint32_t add_to_list( index_list_t& list, uint32_t lit1 ){ return pG( list, lit1, lit1 ); };
  };


  struct xaig_library_t
  {
    xaig_library_t(){
      gates1[0] = xaig_gate_t{ BUF_, 1, &hpcompute_buf_<small_truth_table_t>, &add_buf__to_index_list_<index_list_t> }; 
      gates2[0] = xaig_gate_t{ PA00, 2, &hpcompute_pa00<small_truth_table_t>, &add_pa00_to_index_list_<index_list_t> }; 
      gates2[1] = xaig_gate_t{ PA01, 2, &hpcompute_pa01<small_truth_table_t>, &add_pa01_to_index_list_<index_list_t> }; 
      gates2[2] = xaig_gate_t{ PA10, 2, &hpcompute_pa10<small_truth_table_t>, &add_pa10_to_index_list_<index_list_t> }; 
      gates2[3] = xaig_gate_t{ PA11, 2, &hpcompute_pa11<small_truth_table_t>, &add_pa11_to_index_list_<index_list_t> }; 
      gates2[4] = xaig_gate_t{ EXOR, 2, &hpcompute_exor<small_truth_table_t>, &add_exor_to_index_list_<index_list_t> };
    }

    std::array<xaig_gate_t, 1u> gates1;
    std::array<xaig_gate_t, 5u> gates2;
  };

  template<class LTT>
  struct xaig_candidate_t
  {
    xaig_gate_t gate;
    double cost;
    gate_t type;

    divisor_t<LTT> const& a;
    divisor_t<LTT> const& b;
    uint32_t id;

    xaig_candidate_t(){}
    xaig_candidate_t( uint32_t id, xaig_gate_t gate, double cost, divisor_t<LTT> const& a, divisor_t<LTT> const& b ) : id(id), gate(gate), type(gate.type), cost(cost), a(a), b(b){}
    xaig_candidate_t( uint32_t id, xaig_gate_t gate, double cost, divisor_t<LTT> const& a ) : id(id), gate(gate), type(gate.type), cost(cost), a(a), b(a){}

    uint32_t add_to_list( index_list_t& list, uint32_t lit1, uint32_t lit2 ){ return gate.add_to_list(list, lit1, lit2);};
    uint32_t add_to_list( index_list_t& list ){ return gate.add_to_list(list, a.lit, b.lit);};
    LTT compute( LTT const& tta, LTT const& ttb ){ return gate.compute(tta, ttb);};
    LTT compute(){ return gate.compute(a.func, b.func);};

    double update_cost( double const& costPrevious, double const& minCost, double const& maxCost, bool isNew )
    {
      if( isNew )
      {
        cost = costPrevious + exp(-static_params::beta_synthesis*(cost-minCost)/(maxCost-minCost));
      }
      else
      {
        cost = costPrevious;
      }
      return cost;
    }

  };

public:
  explicit xag_resyn( stats& st ) noexcept
      : st( st ), _database(_resyn, {})
  {
    static_assert( std::is_same_v<typename static_params::base_type, xag_resyn_static_params>, "Invalid static_params type" );
    static_assert( !( static_params::uniform_div_cost && static_params::preserve_depth ), "If depth is to be preserved, divisor depth cost must be provided (usually not uniform)" );
    divisors.reserve( static_params::reserve );
    
    for( uint32_t i{0}; i<static_params::max_support_size; ++i )
      kitty::create_nth_var(_xs[i], i );

    for( uint32_t i{0}; i<4u; ++i )
      kitty::create_nth_var(_xs4[i], i );

    //db = _database.get_database();
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

    _gSPFD.init( care, target );

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
    pos_unate_lits.clear();
    neg_unate_lits.clear();
    binate_divs.clear();
    pos_unate_pairs.clear();
    neg_unate_pairs.clear();

    /* try 0-resub and collect unate literals */
    if( static_params::try_0resub )
    {
      auto const res0 = call_with_stopwatch( st.time_unate, [&]() {
        return find_one_unate();
      } );
      if ( res0 )
      {
        return *res0;
      }
      if ( num_inserts <= 0u )
      {
        return std::nullopt;
      }
    }

    /* sort unate literals and try 1-resub */
    if( static_params::try_1resub )
    {
      call_with_stopwatch( st.time_sort, [&]() {
        sort_unate_lits( pos_unate_lits, 1 );
        sort_unate_lits( neg_unate_lits, 0 );
      } );
      auto const res1or = call_with_stopwatch( st.time_resub1, [&]() {
        return find_div_div( pos_unate_lits, 1 );
      } );
      if ( res1or )
      {
        return *res1or;
      }
      auto const res1and = call_with_stopwatch( st.time_resub1, [&]() {
        return find_div_div( neg_unate_lits, 0 );
      } );
      if ( res1and )
      {
        return *res1and;
      }
  
      if ( binate_divs.size() > static_params::max_binates )
      {
        binate_divs.resize( static_params::max_binates );
      }
  
      if constexpr ( static_params::use_xor )
      {
        /* collect XOR-type unate pairs and try 1-resub with XOR */
        auto const res1xor = find_xor();
        if ( res1xor )
        {
          return *res1xor;
        }
      }
      if ( num_inserts == 1u )
      {
        return std::nullopt;
      }
    }
    
    /* SPFD-based synthesis */
    auto const resi = call_with_stopwatch( st.time_boolean_matching, [&]() {
      return find_spfd_resynthesis( num_inserts );
    } );
    if ( resi )
    {
      return *resi;
    }

    /* collect AND-type unate pairs and sort (both types), then try 2- and 3-resub */
    if( static_params::try_unateness_decomposition )
    {
      call_with_stopwatch( st.time_collect_pairs, [&]() {
        collect_unate_pairs();
      } );
      call_with_stopwatch( st.time_sort, [&]() {
        sort_unate_pairs( pos_unate_pairs, 1 );
        sort_unate_pairs( neg_unate_pairs, 0 );
      } );
      auto const res2or = call_with_stopwatch( st.time_resub2, [&]() {
        return find_div_pair( pos_unate_lits, pos_unate_pairs, 1 );
      } );
      if ( res2or )
      {
        return *res2or;
      }
      auto const res2and = call_with_stopwatch( st.time_resub2, [&]() {
        return find_div_pair( neg_unate_lits, neg_unate_pairs, 0 );
      } );
      if ( res2and )
      {
        return *res2and;
      }
  
      if ( num_inserts >= 3u )
      {
        auto const res3or = call_with_stopwatch( st.time_resub3, [&]() {
          return find_pair_pair( pos_unate_pairs, 1 );
        } );
        if ( res3or )
        {
          return *res3or;
        }
        auto const res3and = call_with_stopwatch( st.time_resub3, [&]() {
          return find_pair_pair( neg_unate_pairs, 0 );
        } );
        if ( res3and )
        {
          return *res3and;
        }
      }
  
      /* choose something to divide and recursive call on the remainder */
      /* Note: dividing = AND the on-set (if using positive unate) or the off-set (if using negative unate)
                          with the *negation* of the divisor/pair (subtracting) */
      uint32_t on_off_div, on_off_pair;
      uint32_t score_div = 0, score_pair = 0;
  
      call_with_stopwatch( st.time_divide, [&]() {
        if ( pos_unate_lits.size() > 0 )
        {
          on_off_div = 1; /* use pos_lit */
          score_div = pos_unate_lits[0].score;
          if ( neg_unate_lits.size() > 0 && neg_unate_lits[0].score > pos_unate_lits[0].score )
          {
            on_off_div = 0; /* use neg_lit */
            score_div = neg_unate_lits[0].score;
          }
        }
        else if ( neg_unate_lits.size() > 0 )
        {
          on_off_div = 0; /* use neg_lit */
          score_div = neg_unate_lits[0].score;
        }
  
        if ( num_inserts > 3u )
        {
          if ( pos_unate_pairs.size() > 0 )
          {
            on_off_pair = 1; /* use pos_pair */
            score_pair = pos_unate_pairs[0].score;
            if ( neg_unate_pairs.size() > 0 && neg_unate_pairs[0].score > pos_unate_pairs[0].score )
            {
              on_off_pair = 0; /* use neg_pair */
              score_pair = neg_unate_pairs[0].score;
            }
          }
          else if ( neg_unate_pairs.size() > 0 )
          {
            on_off_pair = 0; /* use neg_pair */
            score_pair = neg_unate_pairs[0].score;
          }
        }
      } );
  
      if ( score_div > score_pair / 2 ) /* divide with a divisor */
      {
        /* if using pos_lit (on_off_div = 1), modify on-set and use an OR gate on top;
           if using neg_lit (on_off_div = 0), modify off-set and use an AND gate on top
         */ 
        uint32_t const lit = on_off_div ? pos_unate_lits[0].lit : neg_unate_lits[0].lit;
        call_with_stopwatch( st.time_divide, [&]() {
          on_off_sets[on_off_div] &= lit & 0x1 ? get_sign(lit >> 1) : ~get_sign( lit >> 1 );
        } );
  
        auto const res_remain_div = compute_function_rec( num_inserts - 1 );
        if ( res_remain_div )
        {
          auto const new_lit = index_list.add_and( ( lit ^ 0x1 ), *res_remain_div ^ on_off_div );
          return new_lit + on_off_div;
        }
      }
      else if ( score_pair > 0 ) /* divide with a pair */
      {
        fanin_pair const pair = on_off_pair ? pos_unate_pairs[0] : neg_unate_pairs[0];
        call_with_stopwatch( st.time_divide, [&]() {
          if constexpr ( static_params::use_xor )
          {
            if ( pair.lit1 > pair.lit2 ) /* XOR pair: ~(lit1 ^ lit2) = ~lit1 ^ lit2 */
            {
              on_off_sets[on_off_pair] &= ( pair.lit1 & 0x1 ? get_sign( pair.lit1 >> 1 ) : ~get_sign( pair.lit1 >> 1 ) ) ^ ( pair.lit2 & 0x1 ? ~get_sign( pair.lit2 >> 1 ) : get_sign( pair.lit2 >> 1 ) );
            }
            else /* AND pair: ~(lit1 & lit2) = ~lit1 | ~lit2 */
            {
              on_off_sets[on_off_pair] &= ( pair.lit1 & 0x1 ? get_sign( pair.lit1 >> 1 ) : ~get_sign( pair.lit1 >> 1 ) ) | ( pair.lit2 & 0x1 ? get_sign( pair.lit2 >> 1 ) : ~get_sign( pair.lit2 >> 1 ) );
            }
          }
          else /* AND pair: ~(lit1 & lit2) = ~lit1 | ~lit2 */
          {
            on_off_sets[on_off_pair] &= ( pair.lit1 & 0x1 ? get_sign( pair.lit1 >> 1 ) : ~get_sign( pair.lit1 >> 1 ) ) | ( pair.lit2 & 0x1 ? get_sign( pair.lit2 >> 1 ) : ~get_sign( pair.lit2 >> 1 ) );
          }
        } );
  
        auto const res_remain_pair = compute_function_rec( num_inserts - 2 );
        if ( res_remain_pair )
        {
          uint32_t new_lit1;
          if constexpr ( static_params::use_xor )
          {
            new_lit1 = ( pair.lit1 > pair.lit2 ) ? index_list.add_xor( pair.lit1, pair.lit2 ) : index_list.add_and( pair.lit1, pair.lit2 );
          }
          else
          {
            new_lit1 = index_list.add_and( pair.lit1, pair.lit2 );
          }
          auto const new_lit2 = index_list.add_and( new_lit1 ^ 0x1, *res_remain_pair ^ on_off_pair );
          return new_lit2 + on_off_pair;
        }
      }
    }

    return std::nullopt;
  }

  /* Perform SPFD-based resynthesis ( iterative )
     1. Sample a support
     2. Perform resynthesis using the support divisors
   */
   std::optional<uint32_t> find_spfd_resynthesis( uint32_t num_inserts )
   {
      std::set<std::vector<uint32_t>> explored_supports;      
      index_list_t index_list_copy = index_list;

      for( auto i{0u}; i < static_params::max_support_attempts; ++i )
      {
        RNG.seed(i);
        auto supp = static_params::use_greedy_support_selection ? find_support_greedy() : find_support();
                
        if( supp && ( explored_supports.find(*supp)==explored_supports.end() ) )
        { 
          if( static_params::use_boolean_matching )
          {
            auto const res = bmatch_resynthesis( *supp, num_inserts );
            if( res )
            {
              return *res;
            }
          }
          else
          {
            auto const res = spfd_resynthesis( *supp, num_inserts );
            if( res )
            {
              return *res;
            }
          }
          explored_supports.insert( *supp );
        }
        index_list = index_list_copy;
      }
      return std::nullopt;
   }

   std::optional<uint32_t> spfd_resynthesis( std::vector<uint32_t> const& supp, uint32_t& max_num_gates )
   {
      extract_local_functionality( this, _lSPFD, _gSPFD, supp );

      //printf(" %d ", supp.size());
      //kitty::print_binary(_lSPFD.onset);
      //printf(" f|c ");
      //kitty::print_binary(_lSPFD.care);
      //printf(": ");

      index_list_t index_list_copy = index_list;
      divisors_t<small_truth_table_t> divs;

      for( uint32_t iIter{0}; iIter < static_params::max_resynthesis_attempts; ++iIter )
      {
        index_list = index_list_copy;
        divs.clear();
        for( uint32_t i{0}; i<supp.size(); ++i )
          divs.emplace_back( _xs[i], supp[i] << 1u );

        while( divs.size() > 1 && index_list.num_gates() <= max_num_gates )
        {
          auto newDivs = update_divisors( _lSPFD, divs, max_num_gates );
          if( newDivs )
          {
            divs = *newDivs;
          }
          else  
            break;
        }
        //printf("%d<=%d\n", index_list.num_gates(), max_num_gates);

        if( divs.size() == 1 )
        {
          if( kitty::equal( divs[0].func & _lSPFD.care, _lSPFD.onset ) )
          {
            return divs[0].lit;
          }
          else if( kitty::equal( divs[0].func & _lSPFD.care, _lSPFD.offset ) )
          {
            return divs[0].lit ^ 0x1;
          }
          else
          {
            printf("[w]: One divisor not matching\n");
          }
        }
      }
      return std::nullopt;
   }

   std::optional<uint32_t> bmatch_resynthesis( std::vector<uint32_t> const& supp, uint32_t& max_num_gates )
   {

      divisors_t<small_truth_table_t> divs0;
      extract_local_functionality( this, _lSPFD, _gSPFD, supp );

      for( uint32_t i{0}; i<supp.size(); ++i )
        divs0.emplace_back( _xs[i], supp[i] << 1u );
      auto divs = divs0;

      index_list_t index_list_copy = index_list;

      divisors_t<truth_table4_t> div4;

      for( uint32_t iIter{0}; iIter < static_params::max_resynthesis_attempts; ++iIter )
      {
        index_list = index_list_copy;
        
        auto divs = divs0;

        while( divs.size() > 4 && index_list.num_gates() < max_num_gates )
        {
          auto newDivs = update_divisors( _lSPFD, divs, max_num_gates );
          if( newDivs )
          {
            divs = *newDivs;
          }
          else  
            break;
        }
        //printf(": %d<=%d ", index_list.num_gates(), max_num_gates);

          std::vector<uint32_t> sup4;
          for( uint32_t i{0}; i<divs.size(); ++i )
            sup4.push_back(i);
          //kitty::print_binary(_lSPFD.onset);
          //printf(" %d\n", divs.size());
          extract_local_functionality( &divs, _4SPFD, _lSPFD, sup4 );

          //kitty::print_binary(_4SPFD.onset);
          //printf(" ");
          //kitty::print_binary(_4SPFD.care);
          //printf("\n");

          //kitty::print_binary(_4SPFD.care);
          //printf("\n");

          div4.clear();
          for( uint32_t i{0}; i<divs.size(); ++i )
            div4.emplace_back( _xs4[i], divs[i].lit );

          //printf("a\n");

          auto res = boolean_match( div4, max_num_gates );
          
          //printf("b\n");

          //printf(": %d<=%d\n", index_list.num_gates(), max_num_gates);

          if( res )
          {
            return *res;
          }

      }
      return std::nullopt;
   }

std::optional<uint32_t> boolean_match( divisors_t<truth_table4_t> const& divs, uint32_t max_num_gates )
{
  //printf("supp=%d\n", divs.divs.size());
  std::array<uint32_t, 4u> leaves;
  std::array<uint8_t, 4u> permutation;

  //printf("\n");
  //kitty::print_binary(_4SPFD.onset);
  //printf(" f0|m0 ");
  //kitty::print_binary(_4SPFD.care);
  //printf("\n");

  auto specs_npn = exact_npn_canonization( _4SPFD.onset );
  auto tt_npn = std::get<0>( specs_npn );
  auto neg = std::get<1>( specs_npn );
  auto perm = std::get<2>( specs_npn );
  auto const care_npn = apply_npn_transformation( _4SPFD.care, neg & 0xF, perm );
  
  //printf("%d ", neg);
  //for( auto xx : perm )
  //  printf("%d ", xx);
  //printf("\n");
//
  //kitty::print_binary(tt_npn);
  //printf(" fR|mR ");
  //kitty::print_binary(care_npn);
  //printf(".\n");


  auto const structures = _database.get_supergates( tt_npn, ~care_npn, neg, perm );
  if ( structures == nullptr )
  {
    return std::nullopt;
  }
  uint32_t negation = 0;
  for ( auto j = 0u; j < 4u; ++j )
  {
    permutation[perm[j]] = j;
    negation |= ( ( neg >> perm[j] ) & 1 ) << j;
  }
  /* save output negation to apply */
  bool phase = ( neg >> 4u == 1 ) ? true : false;

  {
    auto j = 0u;
    for ( auto const leaf : divs.divs )
    {
      leaves[permutation[j++]] = leaf.lit;
    }

    while ( j < 4u )
      leaves[permutation[j++]] = 0u;
  }

  for ( auto j = 0u; j < 4u; ++j )
  {
    if ( ( negation >> j ) & 1 )
    {
      leaves[j] = leaves[j] ^ 0x1;
    }
  }

  std::unordered_map<uint64_t, uint32_t> existing_nodes; 

  auto& db = _database.get_database();
  
  db.incr_trav_id();

  auto res = create_index_list( db.get_node(structures->at(0).root), leaves, existing_nodes );

  //if( res )
  //  printf("%d %d < %d\n", res->second, index_list.num_gates(), max_num_gates);


  //std::vector<uint32_t> v;
  //v[3]=1;

  if(res && res->second <= max_num_gates)
  {
    if( phase ) 
      return res->first ^ 0x1;
    else
      return res->first;
  }
  return std::nullopt;
}

std::optional<std::pair<int32_t, uint32_t>> create_index_list( node<xag_network> const& n, std::array<uint32_t, 4u> const& leaves, std::unordered_map<uint64_t, uint32_t>& existing_nodes )
{
  auto& db = _database.get_database();
  db.incr_trav_id();

  return create_index_list_rec( n, leaves, existing_nodes );
}

std::optional<std::pair<int32_t, uint32_t>> create_index_list_rec( node<xag_network> const& n, std::array<uint32_t, 4u> const& leaves, std::unordered_map<uint64_t, uint32_t>& existing_nodes )
{
  auto& db = _database.get_database();
  //db.set_visited( n, db.trav_id() );
  if ( db.is_pi( n ) || db.is_constant( n ) )
    return std::nullopt;
  if ( db.visited( n ) == db.trav_id() )
    return std::nullopt;

  db.set_visited( n, db.trav_id() );

  int32_t area = 0;
  uint32_t level = 0;
  bool hashed = true;

  std::array<uint32_t, 2u> node_data;

  db.foreach_fanin( n, [&]( auto const& f, auto i ) {
    node<xag_network> g = db.get_node( f );
    if ( db.is_pi( g ) )
    {
      //printf("PI %d\n", leaves[f.index-1]>>1);
      node_data[i] = db.is_complemented( f ) ? leaves[f.index-1] ^ 0x1 : leaves[f.index-1];
    }
    else
    {
      auto res = create_index_list_rec( g, leaves, existing_nodes );
      if( res )
      {
        node_data[i] = db.is_complemented( f ) ? res->first ^ 0x1 : res->first;
        area += res->second;
      }
      else
        return std::nullopt;
    }
  } );

  uint64_t key0 = node_data[0];
  uint64_t key1 = node_data[1];

  uint32_t new_lit;
  if( db.is_and(n) )
  {
    //printf("and( %d.%d, %d.%d)\n", node_data[0]&1u, node_data[0]>>1, node_data[1]&1u, node_data[1]>>1 );

    if( key0<key1 )
      key0 = ( key0 & 0x00000000FFFFFFFF ) | ( key1 & 0x00000000FFFFFFFF << 32u );
    else
      key0 = ( key1 & 0x00000000FFFFFFFF ) | ( key0 & 0x00000000FFFFFFFF << 32u );

    if( auto search = existing_nodes.find(key0); search != existing_nodes.end() )
    {
      new_lit = search->second;
    }
    else
    {
      new_lit = index_list.add_and(node_data[0], node_data[1]);
      existing_nodes[key0] = new_lit;
      area+=1;
    }
  }
  else if( db.is_xor(n) )
  {
    //printf("xor( %d.%d, %d.%d)\n", node_data[0]&1u, node_data[0]>>1, node_data[1]&1u, node_data[1]>>1 );

    if( key0>key1 )
      key0 = ( key0 & 0x00000000FFFFFFFF ) | ( key1 & 0x00000000FFFFFFFF << 32u );
    else
      key0 = ( key1 & 0x00000000FFFFFFFF ) | ( key0 & 0x00000000FFFFFFFF << 32u );

    if( auto search = existing_nodes.find(key0); search != existing_nodes.end() )
    {
      new_lit = search->second;
    }
    else
    {
      new_lit = index_list.add_xor(node_data[0], node_data[1]);
      existing_nodes[key0] = new_lit;
      area+=1;
    }
  }

  return std::make_pair(new_lit, area );
}

std::optional<std::vector<uint32_t>> find_support_greedy()
{
  _gSPFD.reset();
  std::vector<uint32_t> supp;
  std::vector<uint32_t> candidates;
  double cost, minCost = std::numeric_limits<double>::max();

  while( !_gSPFD.is_covered() && ( supp.size() < static_params::max_support_size ) ) 
  {
    for( auto v = 0u; v < divisors.size(); ++v )
    {
      cost = _gSPFD.evaluate( get_sign(v) );
      if( cost < minCost ) 
      {
        minCost = cost;
        candidates = {v};
      }
      else if( cost == minCost )
      {
        candidates.push_back( v );
      }
    }

    std::uniform_int_distribution<> distrib(0, candidates.size()-1);
    double rnd = distrib(RNG);
    supp.push_back( candidates[rnd] );
    _gSPFD.update( get_sign(candidates[rnd]) );
  }

  if( _gSPFD.is_covered() )
  {
    std::sort(supp.begin(), supp.end());
    return supp;
  }
  return std::nullopt;
}

std::optional<std::vector<uint32_t>> find_support()
{
  _gSPFD.reset();
  std::vector<uint32_t> supp;

  double cost;
  double minCost = std::numeric_limits<double>::max();
  double maxCost = std::numeric_limits<double>::min();
  std::vector<double> costs;

  while( !_gSPFD.is_covered() && ( supp.size() < static_params::max_support_size ) ) 
  {
    costs.clear();
    for( auto v = 0u; v < divisors.size(); ++v )
    {
      cost = _gSPFD.evaluate( get_sign(v) );
      costs.push_back( cost );
      if( cost < minCost ) minCost = cost;
      if( cost > maxCost ) maxCost = cost;
    }

    for( uint32_t i{0}; i<costs.size(); ++i )
      costs[i] = exp(-static_params::beta_support*(costs[i]-minCost)/(maxCost-minCost));
    for( auto v : supp )
      costs[v]=0;
    for( uint32_t i{1}; i<costs.size(); ++i )
    {
      costs[i] += costs[i-1];
    }
    
    std::uniform_real_distribution<> distrib(0, 1);
    double rnd = distrib(RNG);

    for( uint32_t i{0}; i<costs.size(); ++i )
    {
      if( std::isnan( costs[i] ) )
      {
        printf("[w]: NAN\n");
        return std::nullopt;
      }
      if( rnd*costs.back() <= costs[i] )
      {
        supp.push_back(i);
        _gSPFD.update( get_sign(i) );
        break;
      }
    }
  }

  if( _gSPFD.is_covered() )
  {
    std::sort(supp.begin(), supp.end());
    return supp;
  }
  return std::nullopt;
}

template<class DIVS, class SPFDA, class SPFDB>
void extract_local_functionality( DIVS * pDivs, SPFDB & managerB, SPFDA const& managerA, std::vector<uint32_t> const& supp )
{
  auto func = managerB.onset.construct();
  auto care = managerB.care.construct();
  auto jolly = managerA.onset.construct();

  for( uint32_t m{0}; m < (1u << func.num_vars()); ++m )
  {
    if( m < ( 1u << supp.size() ) )
    {
      jolly = jolly | ~jolly;
      for( uint32_t v{0}; v < supp.size(); ++v )
      {
        if( ( m >> v ) & 1u == 1u )
          jolly &= pDivs->get_sign(supp[v]);
        else
          jolly &= ~pDivs->get_sign(supp[v]);
      }
      
      if( kitty::count_ones( jolly & managerA.care ) > 0 )
      {
        kitty::set_bit( care, m );
        if( kitty::count_ones( jolly & managerA.onset & managerA.care ) > 0 )
        {
          kitty::set_bit( func, m );
        }
        else
        {
          kitty::clear_bit( func, m );
        }
      }
      else
      {
        kitty::clear_bit( care, m );
      }
    }
    else
      break;
  }

  //kitty::print_binary(func);
  //printf("\n");
  //kitty::print_binary(care);
  //printf(",\n");

  auto mk0 = managerB.care.construct();
  auto mk1 = managerB.care.construct();
  auto tt0 = managerB.care.construct();
  auto tt1 = managerB.care.construct();
  auto var = managerB.care.construct();

  for( int i{supp.size()-1}; i>=0; --i )
  {
    mk0 = kitty::cofactor0(care, i);
    mk1 = kitty::cofactor1(care, i);
    tt0 = kitty::cofactor0(func, i);
    tt1 = kitty::cofactor1(func, i);

    if( kitty::equal(mk0&mk1&tt0, mk0&mk1&tt1) )
    {
      //printf("erase %d\n", i);
      care = mk0 | mk1;
      kitty::create_nth_var(var,i);
      care &= var;
      func = mk0&tt0 | mk1&tt1;
    }
  }
  //kitty::print_binary(func);
  //printf("\n");
  //kitty::print_binary(care);
  //printf(".\n");

  managerB.init( care, func );
}


template<class LTT, class SPFD>
std::optional<divisors_t<LTT>> update_divisors( SPFD& manager, divisors_t<LTT> & divs, uint32_t max_num_gates )
{
  manager.reset();
  uint32_t nBuffers{0};
  divisors_t<LTT> res;

  double cost;
  double minCost = std::numeric_limits<uint32_t>::max();
  double maxCost = std::numeric_limits<uint32_t>::min();

  std::vector<xaig_candidate_t<LTT>> candidates;

  std::set<uint32_t> USED;
  uint32_t idx{0};

  while( !manager.is_covered() && res.size() < static_params::max_num_spfds )
  {
    for( auto v1 = 0; v1 < divs.size(); ++v1 )
    {
      for( auto gate : _lib.gates1 )
      {
        if( nBuffers >= divs.size()-1 )  continue;

        cost = manager.evaluate( gate.compute( divs.get_sign(v1) ) );
        candidates.emplace_back( idx++, gate, cost, divs[v1] );
        if( cost < minCost ) minCost = cost;
        if( cost > maxCost ) maxCost = cost;
      }

      for( auto v2 = v1+1; v2 < divs.size(); ++v2 )
      {
        for( auto gate : _lib.gates2 )
        {
          cost = manager.evaluate( gate.compute( divs.get_sign(v1), divs.get_sign(v2) ) );
          candidates.emplace_back( idx++, gate, cost, divs[v1], divs[v2] );
          if( cost < minCost ) minCost = cost;
          if( cost > maxCost ) maxCost = cost;
        }
      }
    }

    double costPrevious{0};
    for( auto & cand : candidates )
    {
      costPrevious = cand.update_cost( costPrevious, minCost, maxCost, USED.find(cand.id) == USED.end() );
    }

    double sum = candidates[idx-1].cost;
    std::uniform_real_distribution<> distrib(0, 1);
    double rnd = distrib(RNG);
    bool isUpd{false};

    for( auto & cand : candidates )
    {
      if( rnd <= cand.cost/sum )
      {
        USED.insert(cand.id);
        if( cand.type == BUF_ )  nBuffers++;

        LTT tt = cand.compute();
        res.emplace_back( tt, cand.add_to_list( index_list ) );
        manager.update( tt );  
        isUpd=true;
        break;
      }
    }
    if( !isUpd || index_list.num_gates() > max_num_gates )//+1 )//NEW
      return std::nullopt;
  }

  if( manager.is_covered() )
    return res;
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
      if ( kitty::intersection_is_empty<TT, 1, 1>( get_sign( v ), on_off_sets[0] ) )
      {
        pos_unate_lits.emplace_back( v << 1 );
        unateness[0] = true;
      }
      else if ( kitty::intersection_is_empty<TT, 0, 1>( get_sign( v ), on_off_sets[0] ) )
      {
        pos_unate_lits.emplace_back( v << 1 | 0x1 );
        unateness[1] = true;
      }

      /* check intersection with on-set */
      if ( kitty::intersection_is_empty<TT, 1, 1>( get_sign( v ), on_off_sets[1] ) )
      {
        neg_unate_lits.emplace_back( v << 1 );
        unateness[2] = true;
      }
      else if ( kitty::intersection_is_empty<TT, 0, 1>( get_sign( v ), on_off_sets[1] ) )
      {
        neg_unate_lits.emplace_back( v << 1 | 0x1 );
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
      /* useless unate literal */
      if ( ( unateness[0] && unateness[2] ) || ( unateness[1] && unateness[3] ) )
      {
        pos_unate_lits.pop_back();
        neg_unate_lits.pop_back();
      }
      /* binate divisor */
      else if ( !unateness[0] && !unateness[1] && !unateness[2] && !unateness[3] )
      {
        binate_divs.emplace_back( v );
      }
    }
    return std::nullopt;
  }

  /* Sort the unate literals by the number of minterms in the intersection.
     - For `pos_unate_lits`, `on_off` = 1, sort by intersection with on-set;
     - For `neg_unate_lits`, `on_off` = 0, sort by intersection with off-set
   */
  void sort_unate_lits( std::vector<unate_lit>& unate_lits, uint32_t on_off )
  {
    for ( auto& l : unate_lits )
    {
      l.score = kitty::count_ones( ( l.lit & 0x1 ? ~get_sign( l.lit >> 1 ) : get_sign( l.lit >> 1 ) ) & on_off_sets[on_off] );
    }
    std::sort( unate_lits.begin(), unate_lits.end(), [&]( unate_lit const& l1, unate_lit const& l2 ) {
      return l1.score > l2.score; // descending order
    } );
  }

  void sort_unate_pairs( std::vector<fanin_pair>& unate_pairs, uint32_t on_off )
  {
    for ( auto& p : unate_pairs )
    {
      if constexpr ( static_params::use_xor )
      {
        p.score = ( p.lit1 > p.lit2 ) ? kitty::count_ones( ( ( p.lit1 & 0x1 ? ~get_sign( p.lit1 >> 1 ) : get_sign( p.lit1 >> 1 ) ) ^ ( p.lit2 & 0x1 ? ~get_sign( p.lit2 >> 1 ) : get_sign( p.lit2 >> 1 ) ) ) & on_off_sets[on_off] )
                                      : kitty::count_ones( ( p.lit1 & 0x1 ? ~get_sign( p.lit1 >> 1 ) : get_sign( p.lit1 >> 1 ) ) & ( p.lit2 & 0x1 ? ~get_sign( p.lit2 >> 1 ) : get_sign( p.lit2 >> 1 ) ) & on_off_sets[on_off] );
      }
      else
      {
        p.score = kitty::count_ones( ( p.lit1 & 0x1 ? ~get_sign( p.lit1 >> 1 ) : get_sign( p.lit1 >> 1 ) ) & ( p.lit2 & 0x1 ? ~get_sign( p.lit2 >> 1 ) : get_sign( p.lit2 >> 1 ) ) & on_off_sets[on_off] );
      }
    }
    std::sort( unate_pairs.begin(), unate_pairs.end(), [&]( fanin_pair const& p1, fanin_pair const& p2 ) {
      return p1.score > p2.score; // descending order
    } );
  }

  /* See if there are two unate divisors covering all on-set bits or all off-set bits.
     - For `pos_unate_lits`, `on_off` = 1, try covering all on-set bits by combining two with an OR gate;
     - For `neg_unate_lits`, `on_off` = 0, try covering all off-set bits by combining two with an AND gate
   */
  std::optional<uint32_t> find_div_div( std::vector<unate_lit>& unate_lits, uint32_t on_off )
  {
    for ( auto i = 0u; i < unate_lits.size(); ++i )
    {
      uint32_t const& lit1 = unate_lits[i].lit;
      if ( unate_lits[i].score * 2 < num_bits[on_off] )
      {
        break;
      }
      for ( auto j = i + 1; j < unate_lits.size(); ++j )
      {
        uint32_t const& lit2 = unate_lits[j].lit;
        if ( unate_lits[i].score + unate_lits[j].score < num_bits[on_off] )
        {
          break;
        }
        auto const ntt1 = lit1 & 0x1 ? get_sign( lit1 >> 1 ) : ~get_sign( lit1 >> 1 );
        auto const ntt2 = lit2 & 0x1 ? get_sign( lit2 >> 1 ) : ~get_sign( lit2 >> 1 );
        if ( kitty::intersection_is_empty( ntt1, ntt2, on_off_sets[on_off] ) )
        {
          auto const new_lit = index_list.add_and( ( lit1 ^ 0x1 ), ( lit2 ^ 0x1 ) );
          return new_lit + on_off;
        }
      }
    }
    return std::nullopt;
  }

  std::optional<uint32_t> find_div_pair( std::vector<unate_lit>& unate_lits, std::vector<fanin_pair>& unate_pairs, uint32_t on_off )
  {
    for ( auto i = 0u; i < unate_lits.size(); ++i )
    {
      uint32_t const& lit1 = unate_lits[i].lit;
      for ( auto j = 0u; j < unate_pairs.size(); ++j )
      {
        fanin_pair const& pair2 = unate_pairs[j];
        if ( unate_lits[i].score + pair2.score < num_bits[on_off] )
        {
          break;
        }
        auto const ntt1 = lit1 & 0x1 ? get_sign( lit1 >> 1 ) : ~get_sign( lit1 >> 1 );
        TT ntt2;
        if constexpr ( static_params::use_xor )
        {
          if ( pair2.lit1 > pair2.lit2 )
          {
            ntt2 = ( pair2.lit1 & 0x1 ? get_sign( pair2.lit1 >> 1 ) : ~get_sign( pair2.lit1 >> 1 ) ) ^ ( pair2.lit2 & 0x1 ? ~get_sign( pair2.lit2 >> 1 ) : get_sign( pair2.lit2 >> 1 ) );
          }
          else
          {
            ntt2 = ( pair2.lit1 & 0x1 ? get_sign( pair2.lit1 >> 1 ) : ~get_sign( pair2.lit1 >> 1 ) ) | ( pair2.lit2 & 0x1 ? get_sign( pair2.lit2 >> 1 ) : ~get_sign( pair2.lit2 >> 1 ) );
          }
        }
        else
        {
          ntt2 = ( pair2.lit1 & 0x1 ? get_sign( pair2.lit1 >> 1 ) : ~get_sign( pair2.lit1 >> 1 ) ) | ( pair2.lit2 & 0x1 ? get_sign( pair2.lit2 >> 1 ) : ~get_sign( pair2.lit2 >> 1 ) );
        }

        if ( kitty::intersection_is_empty( ntt1, ntt2, on_off_sets[on_off] ) )
        {
          uint32_t new_lit1;
          if constexpr ( static_params::use_xor )
          {
            if ( pair2.lit1 > pair2.lit2 )
            {
              new_lit1 = index_list.add_xor( pair2.lit1, pair2.lit2 );
            }
            else
            {
              new_lit1 = index_list.add_and( pair2.lit1, pair2.lit2 );
            }
          }
          else
          {
            new_lit1 = index_list.add_and( pair2.lit1, pair2.lit2 );
          }
          auto const new_lit2 = index_list.add_and( ( lit1 ^ 0x1 ), new_lit1 ^ 0x1 );
          return new_lit2 + on_off;
        }
      }
    }
    return std::nullopt;
  }

  std::optional<uint32_t> find_pair_pair( std::vector<fanin_pair>& unate_pairs, uint32_t on_off )
  {
    for ( auto i = 0u; i < unate_pairs.size(); ++i )
    {
      fanin_pair const& pair1 = unate_pairs[i];
      if ( pair1.score * 2 < num_bits[on_off] )
      {
        break;
      }
      for ( auto j = i + 1; j < unate_pairs.size(); ++j )
      {
        fanin_pair const& pair2 = unate_pairs[j];
        if ( pair1.score + pair2.score < num_bits[on_off] )
        {
          break;
        }
        TT ntt1, ntt2;
        if constexpr ( static_params::use_xor )
        {
          if ( pair1.lit1 > pair1.lit2 )
          {
            ntt1 = ( pair1.lit1 & 0x1 ? get_sign( pair1.lit1 >> 1 ) : ~get_sign( pair1.lit1 >> 1 ) ) ^ ( pair1.lit2 & 0x1 ? ~get_sign( pair1.lit2 >> 1 ) : get_sign( pair1.lit2 >> 1 ) );
          }
          else
          {
            ntt1 = ( pair1.lit1 & 0x1 ? get_sign( pair1.lit1 >> 1 ) : ~get_sign( pair1.lit1 >> 1 ) ) | ( pair1.lit2 & 0x1 ? get_sign( pair1.lit2 >> 1 ) : ~get_sign( pair1.lit2 >> 1 ) );
          }
          if ( pair2.lit1 > pair2.lit2 )
          {
            ntt2 = ( pair2.lit1 & 0x1 ? get_sign( pair2.lit1 >> 1 ) : ~get_sign( pair2.lit1 >> 1 ) ) ^ ( pair2.lit2 & 0x1 ? ~get_sign( pair2.lit2 >> 1 ) : get_sign( pair2.lit2 >> 1 ) );
          }
          else
          {
            ntt2 = ( pair2.lit1 & 0x1 ? get_sign( pair2.lit1 >> 1 ) : ~get_sign( pair2.lit1 >> 1 ) ) | ( pair2.lit2 & 0x1 ? get_sign( pair2.lit2 >> 1 ) : ~get_sign( pair2.lit2 >> 1 ) );
          }
        }
        else
        {
          ntt1 = ( pair1.lit1 & 0x1 ? get_sign( pair1.lit1 >> 1 ) : ~get_sign( pair1.lit1 >> 1 ) ) | ( pair1.lit2 & 0x1 ? get_sign( pair1.lit2 >> 1 ) : ~get_sign( pair1.lit2 >> 1 ) );
          ntt2 = ( pair2.lit1 & 0x1 ? get_sign( pair2.lit1 >> 1 ) : ~get_sign( pair2.lit1 >> 1 ) ) | ( pair2.lit2 & 0x1 ? get_sign( pair2.lit2 >> 1 ) : ~get_sign( pair2.lit2 >> 1 ) );
        }

        if ( kitty::intersection_is_empty( ntt1, ntt2, on_off_sets[on_off] ) )
        {
          uint32_t fanin_lit1, fanin_lit2;
          if constexpr ( static_params::use_xor )
          {
            if ( pair1.lit1 > pair1.lit2 )
            {
              fanin_lit1 = index_list.add_xor( pair1.lit1, pair1.lit2 );
            }
            else
            {
              fanin_lit1 = index_list.add_and( pair1.lit1, pair1.lit2 );
            }
            if ( pair2.lit1 > pair2.lit2 )
            {
              fanin_lit2 = index_list.add_xor( pair2.lit1, pair2.lit2 );
            }
            else
            {
              fanin_lit2 = index_list.add_and( pair2.lit1, pair2.lit2 );
            }
          }
          else
          {
            fanin_lit1 = index_list.add_and( pair1.lit1, pair1.lit2 );
            fanin_lit2 = index_list.add_and( pair2.lit1, pair2.lit2 );
          }
          uint32_t const output_lit = index_list.add_and( fanin_lit1 ^ 0x1, fanin_lit2 ^ 0x1 );
          return output_lit + on_off;
        }
      }
    }
    return std::nullopt;
  }

  std::optional<uint32_t> find_xor()
  {
    /* collect XOR-type pairs (d1 ^ d2) & off = 0 or ~(d1 ^ d2) & on = 0, selecting d1, d2 from binate_divs */
    for ( auto i = 0u; i < binate_divs.size(); ++i )
    {
      for ( auto j = i + 1; j < binate_divs.size(); ++j )
      {
        auto const tt_xor = get_sign( binate_divs[i] ) ^ get_sign( binate_divs[j] );
        bool unateness[4] = { false, false, false, false };
        /* check intersection with off-set; additionally check intersection with on-set is not empty (otherwise it's useless) */
        if ( kitty::intersection_is_empty<TT, 1, 1>( tt_xor, on_off_sets[0] ) && !kitty::intersection_is_empty<TT, 1, 1>( tt_xor, on_off_sets[1] ) )
        {
          pos_unate_pairs.emplace_back( binate_divs[i] << 1, binate_divs[j] << 1, true );
          unateness[0] = true;
        }
        if ( kitty::intersection_is_empty<TT, 0, 1>( tt_xor, on_off_sets[0] ) && !kitty::intersection_is_empty<TT, 0, 1>( tt_xor, on_off_sets[1] ) )
        {
          pos_unate_pairs.emplace_back( ( binate_divs[i] << 1 ) + 1, binate_divs[j] << 1, true );
          unateness[1] = true;
        }

        /* check intersection with on-set; additionally check intersection with off-set is not empty (otherwise it's useless) */
        if ( kitty::intersection_is_empty<TT, 1, 1>( tt_xor, on_off_sets[1] ) && !kitty::intersection_is_empty<TT, 1, 1>( tt_xor, on_off_sets[0] ) )
        {
          neg_unate_pairs.emplace_back( binate_divs[i] << 1, binate_divs[j] << 1, true );
          unateness[2] = true;
        }
        if ( kitty::intersection_is_empty<TT, 0, 1>( tt_xor, on_off_sets[1] ) && !kitty::intersection_is_empty<TT, 0, 1>( tt_xor, on_off_sets[0] ) )
        {
          neg_unate_pairs.emplace_back( ( binate_divs[i] << 1 ) + 1, binate_divs[j] << 1, true );
          unateness[3] = true;
        }

        if ( unateness[0] && unateness[2] )
        {
          return index_list.add_xor( ( binate_divs[i] << 1 ), ( binate_divs[j] << 1 ) );
        }
        if ( unateness[1] && unateness[3] )
        {
          return index_list.add_xor( ( binate_divs[i] << 1 ) + 1, ( binate_divs[j] << 1 ) );
        }
      }
    }

    return std::nullopt;
  }

  /* collect AND-type pairs (d1 & d2) & off = 0 or ~(d1 & d2) & on = 0, selecting d1, d2 from binate_divs */
  void collect_unate_pairs()
  {
    for ( auto i = 0u; i < binate_divs.size(); ++i )
    {
      for ( auto j = i + 1; j < binate_divs.size(); ++j )
      {
        collect_unate_pairs_detail<1, 1>( binate_divs[i], binate_divs[j] );
        collect_unate_pairs_detail<0, 1>( binate_divs[i], binate_divs[j] );
        collect_unate_pairs_detail<1, 0>( binate_divs[i], binate_divs[j] );
        collect_unate_pairs_detail<0, 0>( binate_divs[i], binate_divs[j] );
      }
    }
  }

  template<bool pol1, bool pol2>
  void collect_unate_pairs_detail( uint32_t div1, uint32_t div2 )
  {
    /* check intersection with off-set; additionally check intersection with on-set is not empty (otherwise it's useless) */
    if ( kitty::intersection_is_empty<TT, pol1, pol2>( get_sign( div1 ), get_sign( div2 ), on_off_sets[0] ) && !kitty::intersection_is_empty<TT, pol1, pol2>( get_sign( div1 ), get_sign( div2 ), on_off_sets[1] ) )
    {
      pos_unate_pairs.emplace_back( ( div1 << 1 ) + (uint32_t)( !pol1 ), ( div2 << 1 ) + (uint32_t)( !pol2 ) );
    }
    /* check intersection with on-set; additionally check intersection with off-set is not empty (otherwise it's useless) */
    else if ( kitty::intersection_is_empty<TT, pol1, pol2>( get_sign( div1 ), get_sign( div2 ), on_off_sets[1] ) && !kitty::intersection_is_empty<TT, pol1, pol2>( get_sign( div1 ), get_sign( div2 ), on_off_sets[0] ) )
    {
      neg_unate_pairs.emplace_back( ( div1 << 1 ) + (uint32_t)( !pol1 ), ( div2 << 1 ) + (uint32_t)( !pol2 ) );
    }
  }

public:

  inline TT const& operator[](uint32_t idx ) const
  {
      return get_sign( idx );
  }

  inline TT const& get_sign( uint32_t idx ) const
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

  uint32_t num_divisors()
  {
    return divisors.size();
  }

  inline TT const& get_onset() const
  {
    return _gSPFD.onset;
  }

  inline TT const& get_care() const
  {
    return _gSPFD.care;
  }

public:
  index_list_t index_list;

private:
  std::array<TT, 2> on_off_sets;
  std::array<uint32_t, 2> num_bits; /* number of bits in on-set and off-set */
  spfd_manager_t<truth_table_t, 1u << static_params::max_num_spfds> _gSPFD;
  spfd_manager_t<small_truth_table_t, 1u << static_params::max_num_spfds> _lSPFD;
  spfd_manager_t<truth_table4_t, 1u << static_params::max_num_spfds> _4SPFD;
  
  std::array<small_truth_table_t, static_params::max_support_size> _xs;
  std::array<truth_table4_t, 4u> _xs4;
  xaig_library_t _lib{};


  //xag_network db;

  xag_npn_resynthesis<xag_network, xag_network, xag_npn_db_kind::xag_complete> _resyn;
  exact_library<xag_network, xag_npn_resynthesis<xag_network, xag_network, xag_npn_db_kind::xag_complete>> _database;

  const typename static_params::truth_table_storage_type* ptts;
  std::vector<std::conditional_t<static_params::copy_tts, TT, typename static_params::node_type>> divisors;
  divisors_t<TT> _divisors;


  /* positive unate: not overlapping with off-set
     negative unate: not overlapping with on-set */
  std::vector<unate_lit> pos_unate_lits, neg_unate_lits;
  std::vector<uint32_t> binate_divs;
  std::vector<fanin_pair> pos_unate_pairs, neg_unate_pairs;

  stats& st;
}; /* xag_resyn */

#pragma endregion XAG_resyn


} /* namespace spfd */

} /* namespace mockturtle */






//|       c17 |    6 |      6 |      0.00 |      6 |   0.00 |      0.00 |         true |         true |
//|      c432 |  208 |    166 |      0.02 |    166 |   0.00 |      0.01 |         true |         true |
//|      c499 |  398 |    388 |      0.02 |    246 | -36.60 |      0.03 |         true |         true |
//|      c880 |  325 |    296 |      0.02 |    269 |  -9.12 |      0.02 |         true |         true |
//|     c1355 |  502 |    420 |      0.03 |    263 | -37.38 |      0.03 |         true |         true |
//|     c1908 |  341 |    280 |      0.04 |    181 | -35.36 |      0.02 |         true |         true |
//|     c2670 |  716 |    532 |      0.06 |    484 |  -9.02 |      0.04 |         true |         true |
//|     c3540 | 1024 |    787 |      0.33 |    750 |  -4.70 |      0.08 |         true |         true |
//|     c5315 | 1776 |   1277 |      0.15 |   1211 |  -5.17 |      0.10 |         true |         true |
//|     c6288 | 2337 |   1480 |      0.25 |   1426 |  -3.65 |      0.07 |         true |         true |
//|     c7552 | 1469 |   1291 |      0.20 |   1146 | -11.23 |      0.10 |         true |         true |
//  [ 0.00,    0.00,  -36.60,   -9.12,  -37.38,  -35.36,   -9.02,   -4.70,   -5.17,   -3.65,  -11.23 ]

//0000000000101100 m|c 0000000010101111.