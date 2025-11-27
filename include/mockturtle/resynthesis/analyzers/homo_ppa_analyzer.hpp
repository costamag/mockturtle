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
  \file homo_ppa_analyzer.hpp
  \brief Inplace rewrite

  \author Andrea Costamagna
*/

#pragma once

#include "../database/mapped_list.hpp"

#include "../windowing/window_manager.hpp"
#include "../../algorithms/simulation.hpp"
#include "../list_simulator.hpp"
#include "../timing/sensing_analyzer.hpp"
#include "../timing/arrival_analyzer.hpp"
#include "../power/node_simulator.hpp"
#include "../../algorithms/detail/resub_utils.hpp"
#include "../../algorithms/reconv_cut.hpp"

#include <fmt/format.h>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>

namespace mockturtle
{

class ppa
{
public:
  ppa( double power, double delay, double area )
  : power( power ), delay( delay ), area( area )
  {
  }

  ppa( double value )
  : power( value ), delay( value ), area( value )
  {
  }

  bool operator!=( ppa const& other )
  {
    bool cec_power = this->power != other.power;
    bool cec_delay = this->delay != other.delay;
    bool cec_area = this->area != other.area;
    return cec_area || cec_delay || cec_power;
  }

  double operator-( ppa const other ) const
  {
    double d1 = this->power ;//* this->delay * this->area;
    double d2 = other.power ;//* other.delay * other.area;
    return d1-d2;
  }


public:
  double power;
  double delay;
  double area;
};

template<class Ntk, uint32_t num_vars, uint32_t num_steps, uint32_t cube_size, uint32_t max_leaves>
class homo_ppa_analyzer
{
public:
  using cost = double;
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;
  static constexpr cost zero_cost = 0;
  static constexpr cost max_cost = std::numeric_limits<cost>::max();
  static constexpr cost min_cost = std::numeric_limits<cost>::min();

  using power_signature_t = kitty::static_truth_table<cube_size>;
  using node_simulator_t = node_simulator<Ntk, num_steps, cube_size>;
  using list_simulator_t = aig_list_simulator<power_signature_t>;

public:
  homo_ppa_analyzer( Ntk& ntk )
    : ntk( ntk ),
      lib( ntk.get_bindings() ),
      win_manager( ntk ),
      node_to_sim( ntk ),
      sensing( ntk ),
      arrival( ntk ),
      zero_delay( 0 ),
      dyn_switch( 0 )
  {
    size_t extra_space = 1000;
    nodes_sims.reserve( ntk.size() + extra_space );
    mlist_sims.reserve( 20 );
    init();
  }

  void init()
  {
    std::vector<double> input_sensing( ntk.num_pis(), -0.02 );
    std::vector<double> input_arrival( ntk.num_pis(), +0.00 );
    sensing.compute_sensing_times( input_sensing );
    arrival.compute_arrival_times( input_arrival );
    double const delay = arrival.worst_delay();
    power_signature_t const0;
    /* initialize the constants */
    std::array<power_signature_t, 2> const pair0 = {  const0,  const0 };
    std::array<power_signature_t, 2> const pair1 = { ~const0, ~const0 };
    node const n0 = ntk.get_node( ntk.get_constant( false ) );
    node const n1 = ntk.get_node( ntk.get_constant( true ) );
    double const start = 0;
    node_simulator_t sim0( n0, pair0, start, start, delay );
    node_simulator_t sim1( n1, pair1, start, start, delay );
    nodes_sims.push_back( sim0 );
    nodes_sims.push_back( sim1 );
    mlist_sims.push_back( sim0 );
    mlist_sims.push_back( sim1 );

    node_to_sim[ntk.get_node( ntk.get_constant( false ) )] = 0;
    node_to_sim[ntk.get_node( ntk.get_constant( true ) )] = 1;
    ntk.set_visited( 0, ntk.trav_id() );
    ntk.set_visited( 1, ntk.trav_id() );

    std::array<kitty::static_truth_table<cube_size>, 2> workload;
    std::default_random_engine::result_type seed1 = 1;
    std::default_random_engine::result_type seed2 = 42;

    for( int i = 0; i < max_leaves; ++i )
    {
      kitty::create_random( workload[0], seed1+2*i );
      kitty::create_random( workload[1], seed2+2*i+1 );
      node_simulator_t sim( 0, workload, -0.01, -0.02, 0 );
      nodes_sims.emplace_back( sim );
      seed1+=1;
      seed2+=2;
    }
  }

  void reset()
  {
    /* clear information */
    nodes_sims.resize( max_leaves + 2 );
    mlist_sims.resize( 2 );
    node_to_sim.resize();
    node_to_sim[ntk.get_node( ntk.get_constant( false ) )] = 0;
    node_to_sim[ntk.get_node( ntk.get_constant( true ) )] = 1;
  }

  void operator()( node const& n )
  {
    reset();
    dyn_switch = 0;
    zero_delay = 0;
    
    extract_side_windows( n );
    
    if( !win_manager.run( n ) )
    {
      std::cerr << "[e] window construction failed" << std::endl;
      return;
    }

    for( int i = 0; i < win_manager.leaves.size(); ++i )
    {
      node const l = win_manager.leaves[i];
      node_to_sim[l] = i + 2;
      nodes_sims[i+2].update_info( l, sensing.get_time( l ), arrival.get_time( l ) );
    }

    std::vector<node_simulator_t*> sim_ptrs;
    for( int i = 0; i < win_manager.divs.size() ; ++i )
    {
      if( i < win_manager.leaves.size() )
        continue; // this is a leaf node
      node const d = win_manager.divs[i];
      sim_ptrs.clear();
      ntk.foreach_fanin( d, [&]( auto f, auto j ){
        node const ni = ntk.get_node( f );
        sim_ptrs.push_back( &nodes_sims[node_to_sim[ni]] );
      } );
      uint32_t const index { nodes_sims.size() };
      node_to_sim[d] = index;
      nodes_sims.emplace_back( d, sensing.get_time( d ), arrival.get_time( d ) );
      auto const bind = ntk.get_binding( d );
      nodes_sims[index]( bind, list_simulator, sim_ptrs );
    }

    for( int i = 0; i < win_manager.mffc.size() ; ++i )
    {
      node const d = win_manager.mffc[i];
      sim_ptrs.clear();
      ntk.foreach_fanin( d, [&]( auto f, auto j ){
        node const ni = ntk.get_node( f );
        sim_ptrs.push_back( &nodes_sims[node_to_sim[ni]] );
      } );
      uint32_t const index { nodes_sims.size() };
      node_to_sim[d] = index;
      nodes_sims.emplace_back( d, sensing.get_time( d ), arrival.get_time( d ) );
      auto const bind = ntk.get_binding( d );
      nodes_sims[index]( bind, list_simulator, sim_ptrs );
    }
  }

  cost evaluate( mapped_list const& list, std::array<node, num_vars> const& leaves )
  {
    mlist_sims.reserve( list.num_gates() + list.num_pis() + 10 );
    mlist_sims.resize( 2 );

    for( int i = 0; i < num_vars; ++i )
    {
      node const l = leaves[i];
      if( ntk.is_null( l ) )
      {
        mlist_sims.push_back( nodes_sims[0] );
      }
      else
      {
        uint32_t const index = node_to_sim[l];
        mlist_sims.push_back( nodes_sims[index] );
      }
    }

    cost switching = zero_cost;
    cost zero_delay = zero_cost;
    list.foreach_gate( [&]( auto children, auto id ) {

      std::vector<node_simulator_t*> sim_ptrs;

      auto const g = lib[id];
      double sensing = std::numeric_limits<double>::max();
      double arrival = std::numeric_limits<double>::min();
      for( int i = 0; i < children.size(); ++i )
      {
        uint32_t const index = children[i];
        sim_ptrs.emplace_back( &mlist_sims[index] );
        double const rise = g.pins[i].rise_block_delay;
        double const fall = g.pins[i].fall_block_delay;
        auto psim = sim_ptrs[children[i]];
        double const child_sensing = mlist_sims[index].get_sensing();
        double const child_arrival = mlist_sims[index].get_arrival();
        sensing = std::min( sensing, child_sensing + std::min( rise, fall ) );
        arrival = std::max( arrival, child_arrival + std::max( rise, fall ) );
      }

      uint32_t const index = mlist_sims.size();
      mlist_sims.emplace_back( ntk.get_null(), sensing, arrival );
      mlist_sims[index]( g, list_simulator, sim_ptrs );
      auto contr = mlist_sims[index].get_dyn_switch();
      switching += contr;
      zero_delay +=  mlist_sims[index].get_zero_delay();
    } );

    return switching - zero_delay;
  }

  cost heuristic_cost( mapped_list const& list, std::array<node, num_vars> const& leaves )
  {
    std::vector<double> arrivals { 0, 0 };
    std::vector<double> sensings { 0, 0 };
    for( int i = 0; i < num_vars; ++i )
    {
      node const l = leaves[i];
      if( ntk.is_null( l ) )
      {
        mlist_sims.push_back( nodes_sims[0] );
        arrivals.push_back( 1.0*( 1u << 20u ) );
        sensings.push_back( 1.0*( 1u << 20u ) );
      }
      else
      {
        uint32_t const index = node_to_sim[l];
        mlist_sims.push_back( nodes_sims[index] );
        arrivals.push_back( arrival.get_time( l ) );
        sensings.push_back( sensing.get_time( l ) );
      }
    }

    list.foreach_gate( [&]( auto children, auto id ) {
      auto const g = lib[id];
      double ts = std::numeric_limits<double>::max();
      double ta = std::numeric_limits<double>::min();
      for( int i = 0; i < children.size(); ++i )
      {
        uint32_t const index = children[i];
        double const rise = g.pins[i].rise_block_delay;
        double const fall = g.pins[i].fall_block_delay;
        double const child_sensing = sensings[index];
        double const child_arrival = arrivals[index];
        ts = std::min( ts, child_sensing + std::min( rise, fall ) );
        ta = std::max( ta, child_arrival + std::max( rise, fall ) );
      }
      sensings.push_back( ts );
      arrivals.push_back( ta );
    } );
    uint32_t index = sensings.size() - 1;
    double const ts = sensings[index];
    double const ta = arrivals[index];

    cost non_sim = compute_non_simultaneity( ts, ta );

    return non_sim;
  }

  cost evaluate_rewiring( node const& n )
  {
    uint32_t index = node_to_sim[n];
    return nodes_sims[index].get_dyn_switch();
  }

  void extract_side_windows( node const& n )
  {
    /* how close is it to the desired time */
    cost non_simul = 0;
    side_windows.clear();

    ntk.foreach_fanout( n, [&]( auto fo, auto io ){
      node no = ntk.get_node( fo );
      ntk.foreach_fanin( no, [&]( auto fi, auto ii ){
        node ni = ntk.get_node( fi );
        double const arrival_i = arrival.get_time( ni );
        double const sensing_i = sensing.get_time( ni );
        side_windows.push_back( std::make_pair( sensing_i, arrival_i ) );
      } );
    } );

  }

  template<class Leaves_t>
  cost heuristic_cost( node const& n, Leaves_t const& leaves )
  {
    auto ts = sensing.get_time( n );
    auto ta = arrival.get_time( n );
    cost non_sim = compute_non_simultaneity( ts, ta );
    return non_sim;
  }

  cost compute_non_simultaneity( double tSn, double tAn )
  {
    double non_sim = 0;
    double cnt = 0;
    for( auto side : side_windows )
    {
      double tSi = side.first;
      double tAi = side.second;
      double contr = 0;
      cnt += 1;
      if( tAi >= tAn )
      {
        if( tAn >= tSi && tSn >= tSi )
        {
          non_sim += ( tAi - tAn ) + ( tSn - tSi );
        }
        else if( tAn <= tSi )
        {
          non_sim += ( tAn - tSn ) + ( tAi - tSi );
        }
        else
        {
          assert( ( ( tAn > tSi ) && ( tSi > tSn ) ) );
          non_sim += ( tAi - tAn ) + ( tSi - tSn );
        }
      }
      else
      {
        if( tAi >= tSn && tSi >= tSn )
        {
          non_sim += ( tAn - tAi ) + ( tSi - tSn );
        }
        else if( tAi <= tSn )
        {
          non_sim += ( tAi - tSi ) + ( tAn - tSn );
        }
        else
        {
          assert( ( ( tAi > tSn ) && ( tSn > tSi ) ) );
          non_sim += ( tAn - tAi ) + ( tSn - tSi );
        }
      }
    }

    return cnt <= 0 ? 0 : non_sim / cnt;
  }

  template<class Leaves_t>
  cost evaluate( node const& n, Leaves_t const& leaves )
  {
    cost cost_deref = measure_mffc_deref( n, leaves );
    cost cost_ref = measure_mffc_ref( n, leaves );
    if( cost_ref != cost_deref )
      std::cerr << "[e] different costs " << std::endl;

    return cost_deref;
  }

  template<class Leaves_t>
  cost measure_mffc_deref( node const& n, Leaves_t const& leaves )
  {
    /* reference cut leaves */
    for ( auto leaf : leaves )
    {
      if( !ntk.is_null( leaf ) )
        ntk.incr_fanout_size( leaf );
    }

    cost mffc_cost = recursive_deref( n );

    /* dereference leaves */
    for ( auto leaf : leaves )
    {
      ntk.decr_fanout_size( leaf );
    }

    return mffc_cost;
  }

  template<class Leaves_t>
  cost measure_mffc_ref( node const& n, Leaves_t const& leaves )
  {
    /* reference cut leaves */
    for ( auto leaf : leaves )
    {
      ntk.incr_fanout_size( leaf );
    }

    cost mffc_cost = recursive_ref( n );

    /* dereference leaves */
    for ( auto leaf : leaves )
    {
      ntk.decr_fanout_size( leaf );
    }

    return mffc_cost;
  }

  cost recursive_deref( node const& n )
  {
    /* terminate? */
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      return zero_cost;
    
    if( !node_to_sim.has( n ) )
    {
      std::cerr << "[e] mffc node not simulated" << std::endl;
    }

    /* recursively collect nodes */
    cost value{ nodes_sims[node_to_sim[n]].get_dyn_switch() - nodes_sims[node_to_sim[n]].get_zero_delay() };
    ntk.foreach_fanin( n, [&]( auto const& s ) {
      if ( ntk.decr_fanout_size( ntk.get_node( s ) ) == 0 )
      {
        value += recursive_deref( ntk.get_node( s ) );
      }
    } );
    return value;
  }

  cost recursive_ref( node const& n )
  {
    /* terminate? */
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      return zero_cost;

    /* recursively collect nodes */
    cost value{ nodes_sims[node_to_sim[n]].get_dyn_switch() - nodes_sims[node_to_sim[n]].get_zero_delay() };

    ntk.foreach_fanin( n, [&]( auto const& s ) {
      if ( ntk.incr_fanout_size( ntk.get_node( s ) ) == 0 )
      {
        value += recursive_ref( ntk.get_node( s ) );
      }
    } );
    return value;
  }

  double compare( cost const& current, cost const& candidate )
  {
    return current - candidate;
  }

  double get_time( node const& n )
  {
    return arrival.get_time( n );
  }

private:
  Ntk& ntk;
  bindings const lib;

  window_manager<Ntk> win_manager;
  std::array<node_simulator_t, max_leaves> input_sims;
  std::vector<node_simulator_t> nodes_sims;
  std::vector<node_simulator_t> mlist_sims;
  incomplete_node_map<uint32_t, Ntk> node_to_sim;
  sensing_analyzer<Ntk> sensing;
  arrival_analyzer<Ntk> arrival;
  std::vector<std::pair<double,double>> side_windows;

  double dyn_switch;
  double zero_delay;
  list_simulator_t list_simulator;

};

} /* namespace mockturtle */