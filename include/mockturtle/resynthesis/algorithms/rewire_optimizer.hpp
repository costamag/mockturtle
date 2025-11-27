/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2023  EPFL
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
  \file rewire_optimizer.hpp
  \brief Inplace rewrite

  \author Andrea Costamagna
*/

#pragma once

#include "../boolean/dependency_cut.hpp"
#include "../boolean/spfd_manager.hpp"
#include "../../algorithms/cut_enumeration.hpp"
#include "../../algorithms/cut_enumeration/rewrite_cut.hpp"

#include <fmt/format.h>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>


namespace mockturtle
{
  
  /*! \brief Parameters for Rewrite.
  *
  * The data structure `rewrite_params` holds configurable parameters with
  * default arguments for `rewrite`.
  */
 struct rewire_optimizer_params
 {
   rewire_optimizer_params()
   {
     /* 0 < Cut limit < 16 */
     cut_enumeration_ps.cut_limit = 8;
     cut_enumeration_ps.minimize_truth_table = true;
    }
    
  /*! \brief Cut enumeration parameters. */
  cut_enumeration_params cut_enumeration_ps{};

  /*! \brief If true, candidates are only accepted if they do not increase logic depth. */
  bool preserve_depth{ false };
  
  /*! \brief Allow rewrite with multiple structures */
  bool allow_multiple_structures{ true };
  
  /*! \brief Allow zero-gain substitutions */
  bool allow_zero_gain{ false };
  
  /*! \brief Use satisfiability don't cares for optimization. */
  bool use_dont_cares{ false };

  /*! \brief Window size for don't cares calculation. */
  uint32_t window_size{ 8u };
  
  /*! \brief Be verbose. */
  bool verbose{ false };
  
  /* 1 -> consider all the don't cares, 0 -> only one function */
  int effort{ 1 };
  bool greedy{false};
  /* maximum number of fanouts to be used as a new divisor */
  uint32_t fanout_limit{10};
};

/*! \brief Statistics for rewrite.
*
* The data structure `rewrite_stats` provides data collected by running
* `rewrite`.
*/
struct rewire_optimizer_stats
{
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{ 0 };
  
  /*! \brief Expected gain. */
  uint32_t estimated_gain{ 0 };
  
  /*! \brief Candidates */
  uint32_t candidates{ 0 };

  void report() const
  {
    std::cout << fmt::format( "[i] total time       = {:>5.2f} secs\n", to_seconds( time_total ) );
  }
};

template<class Ntk, class Analyzer_t, uint32_t num_vars, uint32_t max_leaves, uint32_t extra_vars = 0>
class rewire_optimizer
{
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;
  using cost = typename Analyzer_t::cost;
  static constexpr cost zero_cost = Analyzer_t::zero_cost;
  static constexpr uint32_t extra_space = 20;
  static constexpr uint32_t offset = 2;
  using signature_t = kitty::static_truth_table<max_leaves>;
  bool VERBOSEq = false;

public:
  rewire_optimizer( Ntk & ntk, Analyzer_t & analyzer, rewire_optimizer_params const& ps, rewire_optimizer_stats& st )
      : ntk( ntk ),
      analyzer( analyzer ),
      ps( ps ),
      st( st ),
      sims( max_leaves + extra_space ),
      node_to_sim( ntk ),
      list_simulator(),
      best_leaves(),
      best_list( num_vars ),
      best_reward( zero_cost )
  {
    init();
  }

  void init()
  {
    /* reserve and initialize the signatures for exhaustive window simulation */
    sims.resize( max_leaves + offset );
    sims[0] = sims[0]^sims[0];
    sims[1] = ~sims[0];
    for( int i = 0; i < max_leaves; ++i )
    {
      kitty::create_nth_var( sims[offset+i], i );
    }
  }

  void reset()
  {
    best_reward = zero_cost;
  }

  template<class WinManager_t>
  void analyze( node const& n, WinManager_t const& window_manager )
  {
if( VERBOSEq )
  window_manager.print();
    std::vector<node> leaves;
    best_reward = zero_cost;
    /* terminate if the window is not valid */
    if( !window_manager.is_valid() )
    {
      std::cout << "Window is not valid" << std::endl;
      return;
    }
    simulate_window( window_manager );
    std::vector<signal> children;
    std::vector<node> fanin;
    std::set<node> current;
    ntk.foreach_fanin( n, [&]( auto fi, auto ii ){
      children.push_back( fi );
      fanin.push_back( ntk.get_node( fi ) );
      current.insert( ntk.get_node( children.back() ) );
    } );

  double cost_curr = analyzer.evaluate( n, fanin, leaves );


    /* rewire the node itself */
    auto simn = get_sig( n );
    window_manager.foreach_divisor( [&]( auto d, auto i ){
if ( VERBOSEq )
  std::cout << "div " << d << std::endl;
if( analyzer.get_arrival( d ) <= analyzer.get_arrival( n ) )
{
  auto simd = get_sig( d );
  if ( kitty::equal( simn, simd ) )
  {
    /* Build the alternative list */
    mapped_list list( children.size() );
    list.create_po( 2 );
    leaves = { d };
    /* Evaluate current cost */
if ( VERBOSEq )
  std::cout << "cost curr " << cost_curr << std::endl;
    double cost_cand = analyzer.evaluate( list, leaves, n );
if ( VERBOSEq )
    std::cout << "cost_cand " << cost_cand << std::endl;
    double reward = cost_curr - cost_cand;
if ( VERBOSEq )
    std::cout << "reward " << reward << std::endl;
          if ( reward > best_reward )
          {
            best_reward = reward;
            best_list = list;
            best_leaves = leaves;
          }
        }
      }
    } );
    
    if ( ntk.fanin_size( n ) > 1 )
    {
      /* extract the SPF */
      auto tt = ntk.node_function( n );
  
      signature_t tmp;
      /* Set to avoid rewiring to current children */
      std::set<node> curr_children;

      ntk.foreach_fanin( n, [&]( auto fi, auto ii ){
        auto tt0 = kitty::cofactor0( tt, ii );
        auto tt1 = kitty::cofactor1( tt, ii );
        auto diff = tt0 ^ tt1;

        signature_t mask;
        int num_minterms = ( 1u << tt.num_vars() );
        for( int m = 0; m < num_minterms; ++m )
        {
          tmp |= ~tmp;
          if( kitty::get_bit( diff, m ) == 1 ) // Boolean difference one -> you must check
          {
            for( int im{0}; im < tt.num_vars(); ++im )
            {
              if( im == ii )
                continue;
              if( ( ( m >> im ) & 0x1 ) == 0x1 )
              {
                tmp &= get_sig( children[im].index );
              }
              else
              {
                tmp &= ~get_sig( children[im].index );
              }
            }
            mask |= tmp;
          }
        }
      
        // double cost_curr = analyzer.evaluate( n, fanin, new_leaves );
        window_manager.foreach_divisor( [&]( auto d, auto i ){
          if( ( current.find( d ) == current.end() )  )
          {
            auto simi = get_sig( fi.index );
            auto simd = get_sig( d );
            if ( d != fi.index )
            {
              if( kitty::equal( simi&mask, simd&mask ) )
              {
                mapped_list list( children.size() );
                std::vector<uint32_t> fanins;
                for( int i = 0; i < children.size(); ++i )
                {
                  fanins.push_back( i + 2 );
                }
                auto lit = list.create_node( fanins, ntk.get_binding_id( n ) );
                list.create_po( lit );
                std::vector<node> new_leaves;
                for( int i = 0; i < children.size(); ++i )
                {
                  new_leaves.push_back( fanin[i] );
                }
                new_leaves[ii] = d;
                for( int i = 0; i < children.size(); ++i )
                {
                  if ( VERBOSEq )
                    std::cout << new_leaves[i] << " ";
                }
                if ( VERBOSEq )
                  std::cout << "cost curr " << cost_curr << std::endl;
                double cost_cand = analyzer.evaluate( list, new_leaves, n );
                if ( VERBOSEq )
                {
                  list.print();
                  std::cout << "cost cand " << cost_cand << std::endl;
                }
                double reward = cost_curr - cost_cand;
                if ( reward > best_reward )
                {
                  best_reward = reward;
                  best_list = list;
                  best_leaves = new_leaves;
                }
              }
            }
          }
        } );
      } );
    }
  }

private:
  template<class WinManager_t>
  bool simulate_window( WinManager_t const& window_manager )
  {
    node_to_sim.resize();
    sims.resize( max_leaves + offset );
    window_manager.foreach_leaf( [&]( node const& l, uint32_t index ){
      node_to_sim[l] = offset + index; 
    } );

    window_manager.foreach_divisor( [&]( node const& d, uint32_t index ){
      if( index >= window_manager.leaves.size() )
      {
        std::vector<signature_t*> sim_ptrs;
        ntk.foreach_fanin( d, [&]( auto fi, auto ii ){
          node const ni = ntk.get_node( fi );
          uint32_t const index = node_to_sim[ni];
          sim_ptrs.push_back( &sims[index] );
        } );
        uint32_t const index = sims.size();
        node_to_sim[d] = index;
        sims.emplace_back();
        auto const alist = ntk.get_binding( d ).list;
        list_simulator( sims[index], alist, sim_ptrs );
      }
    } );

    window_manager.foreach_mffc( [&]( node const& m, uint32_t index ){
      if( !window_manager.is_leaf( m ) )
      {
        std::vector<signature_t*> sim_ptrs;
        ntk.foreach_fanin( m, [&]( auto fi, auto ii ){
          node const ni = ntk.get_node( fi );
          uint32_t const index = node_to_sim[ni];
          sim_ptrs.push_back( &sims[index] );
        } );
        uint32_t const index = sims.size();
        node_to_sim[m] = index;
        sims.emplace_back();
        auto const alist = ntk.get_binding( m ).list;
        list_simulator( sims[index], alist, sim_ptrs );
      }
    } );
    return true;
  }

  signature_t const& get_sig( node const& n )
  {
    uint32_t const index = node_to_sim[n];
    return sims[index];
  }

public:
  bool has_candidate()
  {
    return best_reward > zero_cost;
  }

  cost get_reward() const
  {
    return best_reward;
  }

  mapped_list get_list() const
  {
    return best_list;
  }

  std::vector<node> get_leaves() const
  {
    return best_leaves;
  }

  template<class WinManager_t>
  void print_simulation( WinManager_t const& window_manager )
  {
    std::cout << "leaves" << std::endl;
    window_manager.foreach_leaf( [&]( auto l, auto i ){
      std::cout << l << " ";
      kitty::print_binary( sims[node_to_sim[l]] );
      std::cout << std::endl;
    } );

    std::cout << "divs" << std::endl;
    window_manager.foreach_divisor( [&]( auto d, auto i ){
      std::cout << d << " ";
      kitty::print_binary( sims[node_to_sim[d]] );

      std::cout << " ";
      if( !ntk.is_pi( d ) )
      {
        auto g = ntk.get_binding( d );
        std::cout << g.expression << " ";
        ntk.foreach_fanin( d, [&]( auto fi, auto ii ){
          std::cout << ntk.get_node( fi ) << " ";
        } );
      }
      else
      {
        std::cout << " PI";
      }
      std::cout << std::endl;
    } );

    std::cout << "mffc" << std::endl;
    window_manager.foreach_mffc( [&]( auto m, auto i ){
      std::cout << m << " ";
      kitty::print_binary( sims[node_to_sim[m]] );

      std::cout << " ";
      if( !ntk.is_pi( m ) )
      {
        auto g = ntk.get_binding( m );
        std::cout << g.expression << " ";
        ntk.foreach_fanin( m, [&]( auto fi, auto ii ){
          std::cout << ntk.get_node( fi ) << " ";
        } );
      }
      else
      {
        std::cout << " PI";
      }
      std::cout << std::endl;
    } );
  }

private:
  Ntk & ntk;
  Analyzer_t & analyzer;
  rewire_optimizer_params const& ps;
  rewire_optimizer_stats & st;
  std::vector<signature_t> sims;
  incomplete_node_map<uint32_t, Ntk> node_to_sim;
  aig_list_simulator<signature_t> list_simulator;

  std::vector<node> best_leaves;
  mapped_list best_list;
  cost best_reward;

};

} /* namespace mockturtle */