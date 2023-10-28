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
  \file gen_resyn.hpp
  \brief Generic resynthesis.

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



struct gen_resyn_static_params
{
  using base_type = gen_resyn_static_params;

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
  static constexpr uint32_t max_support_size{ 5u };

  static constexpr bool parallel_support_selection{false};

  /*! \brief Use statistical support selection */
  static constexpr bool use_statistical_support{ false };

  static constexpr uint32_t max_resynthesis_attempts{ 1u };

  static constexpr uint32_t max_support_attempts{ 1u };

  static constexpr bool  is_xaig{ false };

  using truth_table_storage_type = void;
  using node_type = void;

};

template<class TT>
struct gen_resyn_static_params_default : public gen_resyn_static_params
{
  using truth_table_storage_type = std::vector<TT>;
  using node_type = uint32_t;
  using TT6_t = kitty::static_truth_table<max_support_size>;
};

template<class Ntk>
struct gen_resyn_static_params_for_sim_resub : public gen_resyn_static_params
{
  using truth_table_storage_type = incomplete_node_map<kitty::partial_truth_table, Ntk>;
  using node_type = typename Ntk::node;
  using TT6_t = kitty::static_truth_table<max_support_size>;

  static constexpr bool is_xaig = std::is_same_v<typename Ntk::base_type, xag_network>;

};

#pragma region SPFD_resyn

struct gen_resyn_stats
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
  stopwatch<>::duration time_bmatch{ 0 };


  void report() const
  {
    fmt::print( "[i]         <gen_resyn_decompose>\n" );
    fmt::print( "[i]             0-resub      : {:>5.2f} secs\n", to_seconds( time_unate ) );
    fmt::print( "[i]             1-resub      : {:>5.2f} secs\n", to_seconds( time_resub1 ) );
    fmt::print( "[i]             2-resub      : {:>5.2f} secs\n", to_seconds( time_resub2 ) );
    fmt::print( "[i]             3-resub      : {:>5.2f} secs\n", to_seconds( time_resub3 ) );
    fmt::print( "[i]             sort         : {:>5.2f} secs\n", to_seconds( time_sort ) );
    fmt::print( "[i]             collect pairs: {:>5.2f} secs\n", to_seconds( time_collect_pairs ) );
    fmt::print( "[i]             dividing     : {:>5.2f} secs\n", to_seconds( time_divide ) );
  }
};

template<class TT, class static_params = gen_resyn_static_params_default<TT>>
class gen_resyn_spfd
{
public:
  using stats = gen_resyn_stats;
  using index_list_t = large_xag_index_list;
  using truth_table_t = TT;

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

  struct divisor_s_t
  {
    divisor_s_t( typename static_params::TT6_t func, uint32_t lit ) : func(func), lit(lit) {}
    divisor_s_t( typename static_params::TT6_t func ) : func(func) {}
    typename static_params::TT6_t func;
    uint32_t lit;
  };

  enum best_t
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

public:
  explicit gen_resyn_spfd( stats& st ) noexcept
      : st( st )
  {
    static_assert( std::is_same_v<typename static_params::base_type, gen_resyn_static_params>, "Invalid static_params type" );
    static_assert( !( static_params::uniform_div_cost && static_params::preserve_depth ), "If depth is to be preserved, divisor depth cost must be provided (usually not uniform)" );
    divisors.reserve( static_params::reserve );
    for( uint32_t i{0}; i<static_params::max_support_size; ++i )
      kitty::create_nth_var( _s_xs[i], i );

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
    _care = care;

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
    for( auto i{0}; i<divisors.size(); ++i )
      _vNumEdges.push_back(0);

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

    /* try constant-resub */
    auto const resc = call_with_stopwatch( st.time_unate, [&]() {
      return find_cresub();
    } );
    if ( resc )
    {
      return *resc;
    }
    if ( num_inserts == 0u )
    {
      return std::nullopt;
    }

    /* try 0-resub */
    auto const res0 = call_with_stopwatch( st.time_unate, [&]() {
      return find_0resub();
    } );
    if ( res0 )
    {
      return *res0;
    }
    if ( num_inserts == 0u )
    {
      return std::nullopt;
    }

    /* try SPFD-resub */
    auto supp = find_support();
    if( supp && ( supp->size() <= static_params::max_support_size ) )
    {
      auto resS = synthesize( *supp, num_inserts );
      if( resS )
        return resS;
    }

    return std::nullopt;
  }

  /* See if the function is a constant.
   */
  std::optional<uint32_t> find_cresub()
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

    return std::nullopt;
  }

  /* See if there is a divisor covering all on-set bits or all off-set bits.
   */
  std::optional<uint32_t> find_0resub()
  {
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
      /* useless unate literal */

    }
    return std::nullopt;
  }

  std::optional<std::vector<uint32_t>> find_support( double beta = 1000 )
  { 
    std::mt19937 RNGLOC(5);

    std::vector<uint32_t> supp;
    supp.reserve(static_params::max_support_size+1);
    std::vector<uint32_t> vBestDivs;
    reset_masks();

    double numOnes, numEdge, divBest;
    double minEdge = std::numeric_limits<double>::max();
    double maxEdge = std::numeric_limits<double>::min();

    std::vector<double> costs;
    costs.reserve( divisors.size() );

    while( ( _nMasks > _nKilled ) && (supp.size() < static_params::max_support_size) ) 
    {
      double numEdgeTotal{0};
      for( auto m = 0; m < _nMasks; ++m )
      {
        if( !_killed[m] )
        {
          numOnes = kitty::count_ones( _masks[m] & on_off_sets[1] );
          numEdgeTotal += numOnes*(kitty::count_ones( _masks[m] )-numOnes);
        }
      }

      costs.clear();
      costs.push_back(0);
      for( auto v = 1u; v < divisors.size(); ++v )
      {
        numEdge = 0;
        for( auto m = 0; m < _nMasks; ++m )
        {
          if( !_killed[m] )
          {
            numOnes = kitty::count_ones( _masks[m] & get_div(v) & on_off_sets[1] );
            numEdge += numOnes*(kitty::count_ones( get_div(v) & _masks[m] )-numOnes)/numEdgeTotal;
            numOnes = kitty::count_ones( ~get_div(v) & on_off_sets[1] & _masks[m] );
            numEdge += numOnes*(kitty::count_ones( ~get_div(v) & _masks[m] )-numOnes)/numEdgeTotal;
          }
        }
        if( numEdge < minEdge ) minEdge = numEdge;
        if( numEdge > maxEdge ) maxEdge = numEdge;
        costs.push_back(numEdge);
      }

      costs[0]=0;
      for( uint32_t i{1}; i<costs.size(); ++i )
        costs[i] = exp(-beta*(costs[i]-minEdge)/(maxEdge-minEdge));
      for( auto v : supp )
        costs[v]=0;
      for( uint32_t i{1}; i<costs.size(); ++i )
        costs[i] += costs[i-1];

      double sum = costs[costs.size()-1];

      std::uniform_real_distribution<> distrib(0, 1);
      double rnd = distrib(RNGLOC);
      int v=-1;
      for( uint32_t i{1}; i<costs.size(); ++i )
      {
        if( rnd*sum <= costs[i] )
        {
          v = i;
          break;
        }
      }
      if(v<0)
        return std::nullopt;
      supp.push_back(v);
      update_masks(get_div(v));
    }

    if( _nMasks == _nKilled )
    {
      std::sort(supp.begin(), supp.end());
      return supp;
    }
    return std::nullopt;
  }


  std::optional<uint32_t> synthesize( std::vector<uint32_t> supp, uint32_t max_num_gates )
  {
    uint32_t delta, cnt{0};
    std::vector<divisor_s_t> divs;
    for( uint32_t v{0}; v<supp.size(); ++v )
      divs.emplace_back( _s_xs[v], supp[v] << 1 );
    
    uint32_t nSupp = supp.size();
    uint32_t nMnts = 1u << supp.size();
    
    TT jolly = on_off_sets[1];

    _s_care ^= _s_care;


    for( uint32_t m{0}; m < (1u<<static_params::max_support_size); ++m )
    {
      if( m < nMnts )
      {
        jolly = jolly | ~jolly;
        for( uint32_t v{0}; v < nSupp; ++v )
        {
          if( ( m >> v ) & 1u == 1u )
            jolly &= get_div(supp[v]);
          else
            jolly &= ~get_div(supp[v]);
        }
        
        if( kitty::count_ones( jolly ) > 0 )
        {
          kitty::set_bit( _s_care, m );
          jolly = jolly & _care & on_off_sets[1];
          if( kitty::count_ones(jolly) > 0 )
          {
            kitty::set_bit( _s_func, m );
          }
          else
          {
            kitty::clear_bit( _s_func, m );
          }
        }
        else
        {
          kitty::clear_bit( _s_care, m );
        }
      }
      else
        break;
    }
    reset_masks_s();
    
    index_list_t safe_copy_list = index_list;
    std::vector<divisor_s_t> safe_copy_divs = divs;
    uint32_t K = 0;

    while( ( max_num_gates > cnt ) && ( divs.size() > 1 ) && (K < static_params::max_resynthesis_attempts) )
    {
      auto UPD = update_divisors_s( divs, max_num_gates );
      if( UPD )
      {
        delta = UPD->first;
        divs = UPD->second;
        cnt+=delta;
      }
      if( !UPD || ((divs.size()>1)&&(cnt>=max_num_gates)) )
      {
          K++;
          index_list = safe_copy_list;
          divs = safe_copy_divs;
          cnt = 0;
      }
    }

    if( divs.size() == 1 )
    {
      if( kitty::equal( divs[0].func & _s_care, _s_func & _s_care ) )
      {
        return divs[0].lit;
      }
      else if( kitty::equal( ~divs[0].func & _s_care, _s_func & _s_care ) )
      {
        return divs[0].lit ^ 0x1; // MODIFY
      }
      else
      {
        printf("[w]\n");
        return std::nullopt;
      }
    }

    return std::nullopt;
  }

  std::optional<std::pair<uint32_t, std::vector<divisor_s_t>>> update_divisors_s( std::vector<divisor_s_t> const& divs, uint32_t max_num_gates )
  {
    std::vector<divisor_s_t> newDivs;
    uint32_t numGates{0};
    reset_masks_s();
    double numEdge{0};
    double numOnes;
    uint32_t KNT{0};

    uint32_t buffer_counter{0};

    std::set<uint32_t> USED;

    while( ( _s_nMasks > _s_nKilled ) && (newDivs.size() < _s_masks.size()) )
    {
      double numEdgeTotal{0};
      for( auto m = 0; m < _s_nMasks; ++m )
      {
        if( !_s_killed[m] )
        {
          numOnes = kitty::count_ones( _s_masks[m] & _s_func );
          numEdgeTotal += numOnes*(kitty::count_ones( _s_masks[m] )-numOnes);
        }
      }

      double minEdge = std::numeric_limits<uint32_t>::max();
      double maxEdge = std::numeric_limits<uint32_t>::min();
      std::vector<uint32_t> As = {0};
      std::vector<uint32_t> Bs = {0};
      std::vector<typename static_params::TT6_t> Tts = {_s_func};
      std::vector<best_t> Gates ={NONE};
      std::vector<double> Costs = {0};

      if( buffer_counter < divs.size()-1 )
      {
        for( auto v = 0; v<divs.size(); ++v )
        {
          numEdge = 0;
          for( auto m = 0; m < _s_nMasks; ++m )
          {
            if( !_s_killed[m] )
            {
              numOnes = kitty::count_ones( _s_masks[m] & divs[v].func & _s_func );
              numEdge += numOnes*(kitty::count_ones( divs[v].func & _s_masks[m] )-numOnes)/numEdgeTotal;
              numOnes = kitty::count_ones( ~divs[v].func & _s_func & _s_masks[m] );
              numEdge += numOnes*(kitty::count_ones( ~divs[v].func & _s_masks[m] )-numOnes)/numEdgeTotal;
            }
          }
          As.push_back(v);
          Bs.push_back(v);
          Tts.push_back(divs[v].func);
          Gates.push_back(BUF_);
          Costs.push_back(numEdge);
          if( numEdge < minEdge ) minEdge = numEdge;
          if( numEdge > maxEdge ) maxEdge = numEdge;
        }
      }

      typename static_params::TT6_t funcs[5];
      best_t gates[5] = { PA00, PA01, PA10, PA11, EXOR };

      uint32_t nFuncs = static_params::use_xor ? 5 : 4; 

      for( auto v1 = 0; v1<divs.size()-1; ++v1 )
      {
        for( auto v2 = v1+1; v2<divs.size(); ++v2 )
        {
          funcs[0] = ~divs[v1].func & ~divs[v2].func;
          funcs[1] = ~divs[v1].func & divs[v2].func; 
          funcs[2] = divs[v1].func & ~divs[v2].func;
          funcs[3] = divs[v1].func & divs[v2].func;
          if constexpr ( static_params::use_xor )
            funcs[4] = divs[v1].func ^ divs[v2].func;

          for( uint32_t iFn{0}; iFn < nFuncs; ++iFn )
          {
            numEdge = 0;
            for( auto m = 0; m < _s_nMasks; ++m )
            {
              if( !_s_killed[m] )
              {
                numOnes = kitty::count_ones( _s_masks[m] & funcs[iFn] & _s_func );
                numEdge += numOnes*(kitty::count_ones( funcs[iFn] & _s_masks[m] )-numOnes)/numEdgeTotal;
                numOnes = kitty::count_ones( ~funcs[iFn] & _s_func & _s_masks[m] );
                numEdge += numOnes*(kitty::count_ones( ~funcs[iFn] & _s_masks[m] )-numOnes)/numEdgeTotal;
              }
            }

            As.push_back(v1);
            Bs.push_back(v2);
            Tts.push_back(funcs[iFn]);
            Gates.push_back(gates[iFn]);
            Costs.push_back(numEdge);
            
            //printf("num %f\n", numEdge );

            if( numEdge < minEdge ) minEdge = numEdge;
            if( numEdge > maxEdge ) maxEdge = numEdge;
          }
        }
      }

      //printf("min %f max %f \n", minEdge, maxEdge);

      double beta{100};
      for( uint32_t i{1}; i<Costs.size(); ++i )
      {
        if( USED.find(i) == USED.end() )
        {
          Costs[i] = Costs[i-1] + exp(-beta*(Costs[i]-minEdge)/(maxEdge-minEdge));
        }
        else
          Costs[i] = Costs[i-1];
      }
      double sum = Costs[Costs.size()-1];
//
      std::uniform_real_distribution<> distrib(0, 1);
      double rnd = distrib(RNGSPFD);
      int v = -1;
      //printf("%f:\n", rnd );
      for( uint32_t i{1}; i<Costs.size(); ++i )
      {
        //printf("%f ", Costs[i]);
        if( rnd <= Costs[i]/sum )
        {
          v = i;
          break;
        }
      }

      if(v<0)
        return std::nullopt;
      uint32_t a = As[v];
      uint32_t b = Bs[v];
      best_t gate = Gates[v];
      typename static_params::TT6_t tt = Tts[v];
      USED.insert(v);
      
      switch (gate)
      {
      case PA00:
        newDivs.emplace_back( tt, index_list.add_and( divs[a].lit ^ 0x1, divs[b].lit ^ 0x1 ) );
        numGates+=1;
        break;
      case PA01:
        newDivs.emplace_back( tt, index_list.add_and( divs[a].lit ^ 0x1, divs[b].lit ) );
        numGates+=1;
        break;
      case PA10:
        newDivs.emplace_back( tt, index_list.add_and( divs[a].lit, divs[b].lit ^ 0x1) );
        numGates+=1;
        break;
      case PA11:
        newDivs.emplace_back( tt, index_list.add_and( divs[a].lit, divs[b].lit ) );
        numGates+=1;
        break;
      case EXOR:
        newDivs.emplace_back( tt, index_list.add_xor( divs[a].lit, divs[b].lit ) );
        numGates += 1 ;
        break;
      case BUF_:
        newDivs.emplace_back( tt, divs[a].lit );
        buffer_counter++;
        break;
      default:
        return std::nullopt;
        //printf("HERE\n");
        assert(0);
        break;
      }
      if( numGates > max_num_gates )//NEW
        return std::nullopt;//NEW

      update_masks_s( tt );
    }
    /* termination options */
    if( (numGates == 0) && (newDivs.size() == divs.size()))
    {
      printf("[w] no gates added\n");
      return std::nullopt;
    }
    return std::make_pair(numGates, newDivs);
  }

  inline void update_masks( TT const& tt )
  {
    for( uint32_t iMask{0}; iMask < _nMasks; ++iMask )
    {
      if( _killed[iMask] )
      {
        _killed[_nMasks+iMask] = true;
        _nKilled++;
      }
      else
      {
        _killed[_nMasks+iMask] = false;
        _masks[_nMasks+iMask] = _masks[iMask] & tt;
        _masks[iMask] = _masks[iMask] & ~tt;

        if( kitty::count_ones( on_off_sets[1] & _masks[_nMasks+iMask] ) == 0 )
        {
          _killed[_nMasks+iMask] = true;
          _nKilled++;
        }
        else if( kitty::equal( on_off_sets[1] & _masks[_nMasks+iMask], _masks[_nMasks+iMask] ) )
        {
          _killed[_nMasks+iMask] = true;
          _nKilled++;
        }

        if( kitty::count_ones( on_off_sets[1] & _masks[iMask] ) == 0 )
        {
          _killed[iMask] = true;
          _nKilled++;
        }
        else if( kitty::equal( on_off_sets[1] & _masks[iMask], _masks[iMask] ) )
        {
          _killed[iMask] = true;
          _nKilled++;
        }
      }
    }
    _nMasks=_nMasks*2;
  }

  void reset_masks()
  {
    _masks[0] = _care;
    _nMasks = 1u;
    _killed[0] = false;
    _nKilled = 0;
  }

  void update_masks_s( typename static_params::TT6_t tt )
  {
    for( uint32_t iMask{0}; iMask < _s_nMasks; ++iMask )
    {
      if( _s_killed[iMask] )
      {
        _s_killed[_s_nMasks+iMask] = true;
        _s_nKilled++;
      }
      else
      {
        _s_killed[_s_nMasks+iMask] = false;
        _s_masks[_s_nMasks+iMask] = _s_masks[iMask] & tt;
        _s_masks[iMask] = _s_masks[iMask] & ~tt;

        if( kitty::count_ones( _s_func & _s_masks[_s_nMasks+iMask] ) == 0 )
        {
          _s_killed[_s_nMasks+iMask] = true;
          _s_nKilled++;
        }
        else if( kitty::equal( _s_func & _s_masks[_s_nMasks+iMask], _s_masks[_s_nMasks+iMask] ) )
        {
          _s_killed[_s_nMasks+iMask] = true;
          _s_nKilled++;
        }

        if( kitty::count_ones( _s_func & _s_masks[iMask] ) == 0 )
        {
          _s_killed[iMask] = true;
          _s_nKilled++;
        }
        else if( kitty::equal( _s_func & _s_masks[iMask], _s_masks[iMask] ) )
        {
          _s_killed[iMask] = true;
          _s_nKilled++;
        }
      }
    }
    _s_nMasks=_s_nMasks*2;
  }

  void reset_masks_s()
  {
    _s_masks[0] = _s_care;
    _s_nMasks = 1u;
    _s_killed[0] = false;
    _s_nKilled = 0;
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
      if ( kitty::intersection_is_empty<TT, 1, 1>( get_div( v ), on_off_sets[0] ) )
      {
        pos_unate_lits.emplace_back( v << 1 );
        unateness[0] = true;
      }
      else if ( kitty::intersection_is_empty<TT, 0, 1>( get_div( v ), on_off_sets[0] ) )
      {
        pos_unate_lits.emplace_back( v << 1 | 0x1 );
        unateness[1] = true;
      }

      /* check intersection with on-set */
      if ( kitty::intersection_is_empty<TT, 1, 1>( get_div( v ), on_off_sets[1] ) )
      {
        neg_unate_lits.emplace_back( v << 1 );
        unateness[2] = true;
      }
      else if ( kitty::intersection_is_empty<TT, 0, 1>( get_div( v ), on_off_sets[1] ) )
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
      l.score = kitty::count_ones( ( l.lit & 0x1 ? ~get_div( l.lit >> 1 ) : get_div( l.lit >> 1 ) ) & on_off_sets[on_off] );
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
        p.score = ( p.lit1 > p.lit2 ) ? kitty::count_ones( ( ( p.lit1 & 0x1 ? ~get_div( p.lit1 >> 1 ) : get_div( p.lit1 >> 1 ) ) ^ ( p.lit2 & 0x1 ? ~get_div( p.lit2 >> 1 ) : get_div( p.lit2 >> 1 ) ) ) & on_off_sets[on_off] )
                                      : kitty::count_ones( ( p.lit1 & 0x1 ? ~get_div( p.lit1 >> 1 ) : get_div( p.lit1 >> 1 ) ) & ( p.lit2 & 0x1 ? ~get_div( p.lit2 >> 1 ) : get_div( p.lit2 >> 1 ) ) & on_off_sets[on_off] );
      }
      else
      {
        p.score = kitty::count_ones( ( p.lit1 & 0x1 ? ~get_div( p.lit1 >> 1 ) : get_div( p.lit1 >> 1 ) ) & ( p.lit2 & 0x1 ? ~get_div( p.lit2 >> 1 ) : get_div( p.lit2 >> 1 ) ) & on_off_sets[on_off] );
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
        auto const ntt1 = lit1 & 0x1 ? get_div( lit1 >> 1 ) : ~get_div( lit1 >> 1 );
        auto const ntt2 = lit2 & 0x1 ? get_div( lit2 >> 1 ) : ~get_div( lit2 >> 1 );
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
        auto const ntt1 = lit1 & 0x1 ? get_div( lit1 >> 1 ) : ~get_div( lit1 >> 1 );
        TT ntt2;
        if constexpr ( static_params::use_xor )
        {
          if ( pair2.lit1 > pair2.lit2 )
          {
            ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) ) ^ ( pair2.lit2 & 0x1 ? ~get_div( pair2.lit2 >> 1 ) : get_div( pair2.lit2 >> 1 ) );
          }
          else
          {
            ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) ) | ( pair2.lit2 & 0x1 ? get_div( pair2.lit2 >> 1 ) : ~get_div( pair2.lit2 >> 1 ) );
          }
        }
        else
        {
          ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) ) | ( pair2.lit2 & 0x1 ? get_div( pair2.lit2 >> 1 ) : ~get_div( pair2.lit2 >> 1 ) );
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
            ntt1 = ( pair1.lit1 & 0x1 ? get_div( pair1.lit1 >> 1 ) : ~get_div( pair1.lit1 >> 1 ) ) ^ ( pair1.lit2 & 0x1 ? ~get_div( pair1.lit2 >> 1 ) : get_div( pair1.lit2 >> 1 ) );
          }
          else
          {
            ntt1 = ( pair1.lit1 & 0x1 ? get_div( pair1.lit1 >> 1 ) : ~get_div( pair1.lit1 >> 1 ) ) | ( pair1.lit2 & 0x1 ? get_div( pair1.lit2 >> 1 ) : ~get_div( pair1.lit2 >> 1 ) );
          }
          if ( pair2.lit1 > pair2.lit2 )
          {
            ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) ) ^ ( pair2.lit2 & 0x1 ? ~get_div( pair2.lit2 >> 1 ) : get_div( pair2.lit2 >> 1 ) );
          }
          else
          {
            ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) ) | ( pair2.lit2 & 0x1 ? get_div( pair2.lit2 >> 1 ) : ~get_div( pair2.lit2 >> 1 ) );
          }
        }
        else
        {
          ntt1 = ( pair1.lit1 & 0x1 ? get_div( pair1.lit1 >> 1 ) : ~get_div( pair1.lit1 >> 1 ) ) | ( pair1.lit2 & 0x1 ? get_div( pair1.lit2 >> 1 ) : ~get_div( pair1.lit2 >> 1 ) );
          ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) ) | ( pair2.lit2 & 0x1 ? get_div( pair2.lit2 >> 1 ) : ~get_div( pair2.lit2 >> 1 ) );
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
        auto const tt_xor = get_div( binate_divs[i] ) ^ get_div( binate_divs[j] );
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
    if ( kitty::intersection_is_empty<TT, pol1, pol2>( get_div( div1 ), get_div( div2 ), on_off_sets[0] ) && !kitty::intersection_is_empty<TT, pol1, pol2>( get_div( div1 ), get_div( div2 ), on_off_sets[1] ) )
    {
      pos_unate_pairs.emplace_back( ( div1 << 1 ) + (uint32_t)( !pol1 ), ( div2 << 1 ) + (uint32_t)( !pol2 ) );
    }
    /* check intersection with on-set; additionally check intersection with off-set is not empty (otherwise it's useless) */
    else if ( kitty::intersection_is_empty<TT, pol1, pol2>( get_div( div1 ), get_div( div2 ), on_off_sets[1] ) && !kitty::intersection_is_empty<TT, pol1, pol2>( get_div( div1 ), get_div( div2 ), on_off_sets[0] ) )
    {
      neg_unate_pairs.emplace_back( ( div1 << 1 ) + (uint32_t)( !pol1 ), ( div2 << 1 ) + (uint32_t)( !pol2 ) );
    }
  }

private:
  std::array<TT, 2> on_off_sets;
  std::array<uint32_t, 2> num_bits; /* number of bits in on-set and off-set */

  TT _care;
  std::array<TT, 256> _masks;
  std::array<bool, 256> _killed;
  uint64_t _nMasks{1};
  uint32_t _nKilled{0};
  std::vector<uint32_t> _vNumEdges;


  std::array<typename static_params::TT6_t, static_params::max_support_size> _s_xs;
  typename static_params::TT6_t _s_care;
  typename static_params::TT6_t _s_func;
  std::array<typename static_params::TT6_t,256> _s_masks;
  std::array<bool, 256> _s_killed;
  uint64_t _s_nMasks{1};
  uint32_t _s_nKilled{0};


  const typename static_params::truth_table_storage_type* ptts;
  std::vector<std::conditional_t<static_params::copy_tts, TT, typename static_params::node_type>> divisors;

  index_list_t index_list;

  /* positive unate: not overlapping with off-set
     negative unate: not overlapping with on-set */
  std::vector<unate_lit> pos_unate_lits, neg_unate_lits;
  std::vector<uint32_t> binate_divs;
  std::vector<fanin_pair> pos_unate_pairs, neg_unate_pairs;

  stats& st;
}; /* gen_resyn_spfd */

#pragma endregion SPFD_resyn

} /* namespace mockturtle */


