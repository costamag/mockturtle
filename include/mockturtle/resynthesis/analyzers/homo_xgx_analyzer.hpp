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
  \file homo_xgx_analyzer.hpp
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
#include "../timing/switch_analyzer.hpp"
#include "../power/node_simulator.hpp"
#include "../power/analyzers/convo_analyzer.hpp"
#include "../../algorithms/detail/resub_utils.hpp"
#include "../../algorithms/reconv_cut.hpp"
#include "analyzers_common.hpp"

#include <fmt/format.h>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>

#include <set>


namespace mockturtle
{
  
  struct activity_interval
  {
    double sensing;
    double arrival;
  };

  template<class Ntk, uint32_t num_vars, uint32_t num_steps, uint32_t cube_size, uint32_t max_leaves>
  class homo_xgx_analyzer
{
  public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;
  using cost = double;
  static constexpr bool pass_manager = true;
  static constexpr cost zero_cost = 0;
  static constexpr cost min_cost = std::numeric_limits<cost>::min();
  static constexpr cost max_cost = std::numeric_limits<cost>::max() - 1000;
  bool VERBOSEg = false;

  using truth_table_t = kitty::static_truth_table<cube_size>;
  using node_simulator_t = node_simulator<Ntk, num_steps, cube_size>;
  using list_simulator_t = aig_list_simulator<truth_table_t>;

public:
  homo_xgx_analyzer( Ntk& ntk, analyzer_params ps = {} )
    : ntk( ntk ),
      lib( ntk.get_bindings() ),
      node_to_sim( ntk ),
      sensing( ntk ),
      arrival( ntk ),
      zero_delay( 0 ),
      dyn_switch( 0 ),
      ps( ps )
  {
    size_t extra_space = 1000;
    nodes_sims.reserve( ntk.size() + extra_space );
    nodes.reserve( ntk.size() + extra_space );
    mlist_sims.reserve( 20 );
    init();
  }

  homo_xgx_analyzer( Ntk& ntk, std::vector<double> arrival = {} )
    : ntk( ntk ),
      lib( ntk.get_bindings() ),
      node_to_sim( ntk ),
      sensing( ntk ),
      arrival( ntk ),
      zero_delay( 0 ),
      dyn_switch( 0 )
  {
    size_t extra_space = 1000;
    nodes_sims.reserve( ntk.size() + extra_space );
    nodes.reserve( ntk.size() + extra_space );
    mlist_sims.reserve( 20 );
    init();
  }

  void update()
  {
  }

  void init()
  {
    std::vector<double> input_sensing( ntk.num_pis(), -1 );
    std::vector<double> input_arrival( ntk.num_pis(), 0 );
    sensing.compute_sensing_times( input_sensing );
    arrival.compute_arrival_times( input_arrival );
    double const delay = arrival.worst_delay();
    truth_table_t const0;
    /* initialize the constants */
    std::array<truth_table_t, 2> const pair0 = {  const0,  const0 };
    std::array<truth_table_t, 2> const pair1 = { ~const0, ~const0 };
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
    std::default_random_engine::result_type seed1 = 12;
    std::default_random_engine::result_type seed2 = 7;
    for ( int i = 0; i < max_leaves; ++i )
    {
      kitty::create_random( workload[0], seed1*(i+1) );
      kitty::create_random( workload[1], seed2*(2*i+1) );
      node_simulator_t sim( 0, workload, -1, 0 );
      nodes_sims.emplace_back( sim );
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

  /*! \brief Run the analysis for node n
   * Takes as input the node and the window manager and run the analysis to simulate the 
   * window with the workload to evaluate the dynamic power.
  */
  void operator()( node const& n, window_manager<Ntk, max_leaves> const& win_manager )
  {
    reset();
    dyn_switch = 0;
    zero_delay = 0;
    //win_manager.print();
    for( int i = 0; i < win_manager.leaves.size(); ++i )
    {
      node const l = win_manager.leaves[i];

      node_to_sim[l] = i + 2;
      nodes_sims[i+2].set_sensing( sensing.get_time( l ) );
      nodes_sims[i+2].set_arrival( arrival.get_time( l ) );
    }

    window_nodes.clear();
    for ( auto l : win_manager.leaves )
    {
      window_nodes.insert(l);
    }
    for( auto d : win_manager.divs )
    {
      window_nodes.insert(d);
    }
    for( auto d : win_manager.mffc )
    {
      window_nodes.insert(d);
    }

    std::vector<node_simulator_t*> sim_ptrs;
    for ( int i = 0; i < win_manager.divs.size() ; ++i )
    {
      if ( i < win_manager.leaves.size() )
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

  bool has_node( node n )
  {
    return window_nodes.find( n ) != window_nodes.end();
  }

  cost evaluate( mapped_list const& list, std::array<double, num_vars> const& leaves )
  {
    std::cerr << "[e] Invalid cost function" << std::endl;
    return max_cost;
  }

  template<class Leaves_t>
  cost evaluate( mapped_list const& list, Leaves_t const& leaves, node const& n )
  {
    for ( auto const& l : leaves )
    {
      if ( !ntk.is_null( l ) )
      {
        if ( !has_node( l ) )
        {
          return max_cost;
        }
      }
    }

    auto act = evaluate_interval( list, leaves );
    if ( ps.preserve_depth )
    {
      if ( act.arrival > arrival.get_required( n ) )
      {
        return max_cost;
      }
    }

    if ( ps.preserve_interval )
    {
      auto req = evaluate_desired( n );
      if ( ( act.arrival > req.arrival ) )//|| ( act.sensing < req.sensing ) )
      {
        return max_cost;
      }
    }

    mlist_sims.resize( 2 );
    mlist_sims.reserve( list.num_gates() + list.num_pis() + 10 );

    std::set<uint32_t> null_set;

    for ( int i = 0; i < leaves.size(); ++i )
    {
      node const l = leaves[i];
      if( ntk.is_null( l ) )
      {
        mlist_sims.push_back( nodes_sims[0] );
        null_set.insert( i + 2 );
      }
      else
      {
        uint32_t const index = node_to_sim[l];
        mlist_sims.push_back( nodes_sims[index] );
        mlist_sims.back().set_arrival( arrival.get_time( l ) );
        mlist_sims.back().set_sensing( sensing.get_time( l ) );
      }
    }

    while ( mlist_sims.size() < ( num_vars + 2 ) )
    {
      null_set.insert( mlist_sims.size() );
      mlist_sims.push_back( nodes_sims[0] );
    }
    cost switching = 0;

    std::vector<uint32_t> fanout_sizes( list.size(), 0u );
    list.foreach_gate( [&]( auto children, auto id ) {
      for ( auto lit : children )
        fanout_sizes[lit] += 1;
    } );
    fanout_sizes[list.po_at(0)] = ntk.fanout_size( n );

    bool is_valid = true;
    list.foreach_gate( [&]( auto children, auto id ) {
      std::vector<node_simulator_t*> sim_ptrs( children.size() );

      auto const g = lib[id];
      double sensing = std::numeric_limits<double>::max();
      double arrival = std::numeric_limits<double>::min();
      for ( int i = 0; i < children.size(); ++i )
      {
        uint32_t const index = children[i];
        if ( null_set.find( index ) != null_set.end() )
        {
          if ( VERBOSEg )
            std::cout << index << " not valid " << std::endl;
          is_valid = false;
        }
        else
        {
          sim_ptrs[i] = &mlist_sims[index];
          double const rise = g.pins[i].rise_block_delay;
          double const fall = g.pins[i].fall_block_delay;
          double const child_sensing = mlist_sims[index].get_sensing();
          double const child_arrival = mlist_sims[index].get_arrival();
          sensing = std::min( sensing, child_sensing + std::min( rise, fall ) );
          arrival = std::max( arrival, child_arrival + std::max( rise, fall ) );
        }
      }

      uint32_t const index = mlist_sims.size();
      mlist_sims.emplace_back( ntk.get_null(), sensing, arrival );

      mlist_sims[index]( g, list_simulator, sim_ptrs );
      auto contr = mlist_sims[index].get_dyn_switch();
      switching += ps.use_fanouts ? contr*fanout_sizes[index] : contr;
      if ( VERBOSEg )
        std::cout << index << " contr " << contr << "x" << fanout_sizes[index] << "=" << contr*fanout_sizes[index] << std::endl;
    } );
    if ( !is_valid )
    {
      return max_cost;
    }

    for ( int i = 0; i < leaves.size(); ++i )
    {
      node const l = leaves[i];
      if ( !ntk.is_null( l ) )
      {
        if( ps.use_fanouts && !ntk.is_pi( l ) )
        {
          if ( VERBOSEg )
            std::cout << mlist_sims[i+2].get_dyn_switch() << "x" << fanout_sizes[2+i] << "=" << mlist_sims[i+2].get_dyn_switch()*fanout_sizes[2+i] << std::endl;
          switching += mlist_sims[i+2].get_dyn_switch()*fanout_sizes[2+i];
        }
      }
    }
    return switching;
  }

  double get_sensing( node n )
  {
    return sensing.get_time( n );
  }

  double get_arrival( node n )
  {
    return arrival.get_time( n );
  }

  template<class Cut_t>
  activity_interval evaluate_interval( mapped_list const& list, Cut_t const& leaves )
  {
    std::vector<double> arrival_times;
    std::vector<double> sensing_times;

    for ( int i = 0; i < 2; ++i )
    {
      arrival_times.push_back( 0 );
      sensing_times.push_back( 0 );
    }

    list.foreach_pi( [&]( auto v ) {
      if ( ( v - 2 ) < leaves.size() )
      {
        node n = leaves[ v - 2 ];
        if ( ntk.is_null( n ) )
        {
          arrival_times.push_back( max_cost );
          sensing_times.push_back( 0 );
        }
        else
        {
          arrival_times.push_back( arrival.get_time( n ) );
          sensing_times.push_back( sensing.get_time( n ) );
        }
      }
    } );
    while( arrival_times.size() < ( num_vars + 2 ) )
    {
      arrival_times.push_back( max_cost );
      sensing_times.push_back( 0 );
    }
    
    double node_arrival;
    double node_sensing;
    list.foreach_gate( [&]( auto children, auto id ){
      node_arrival = std::numeric_limits<double>::min();
      node_sensing = std::numeric_limits<double>::max();
      for( int i = 0; i < children.size(); ++i )
      {
        auto const g = lib[id];
        uint32_t const child = children[i];
        double const rise = g.pins[i].rise_block_delay;
        double const fall = g.pins[i].fall_block_delay;
        double const child_arrival = arrival_times[child];
        double const child_sensing = sensing_times[child];
        node_arrival = std::max( node_arrival, child_arrival + std::max( rise, fall ) );
        node_sensing = std::min( node_sensing, child_sensing + std::min( rise, fall ) );
      }
      arrival_times.push_back( node_arrival );
      sensing_times.push_back( node_sensing );
    } );
    activity_interval act;
    act.sensing = sensing_times[list.po_at(0)];
    act.arrival = arrival_times[list.po_at(0)];
    return act;
  }

  activity_interval evaluate_desired( node const& n )
  {
    activity_interval res;
    res.sensing = sensing.get_time( n );
    res.arrival = arrival.get_time( n );
    ntk.foreach_fanout( n, [&]( auto fo, auto io ){
      auto no = ntk.get_node( fo );
      auto const go = ntk.get_binding( no );
      int i = -1;
      ntk.foreach_fanin( no, [&]( auto fi, auto ii ){
        auto ni = ntk.get_node( fi );
        if ( ni == n )
        {
          i = ii;
        }
      } );
      assert( i >= 0 );
      double const rise = go.pins[i].rise_block_delay;
      double const fall = go.pins[i].fall_block_delay;
      res.arrival = std::min( res.arrival, arrival.get_time( no ) - std::max( rise, fall ) );
      res.sensing = std::max( res.sensing, sensing.get_time( no ) - std::min( rise, fall ) );
    } );
    if ( VERBOSEg )
      std::cout << res.sensing << "=tS " << res.arrival << "=tA " << std::endl;
    return res;
  }

  cost evaluate_rewiring( node const& n, node const& ncurr  )
  {
    if ( !has_node( n ) || !has_node( ncurr ) )
    {
      return max_cost;
    }

    if( ps.preserve_depth )
    {
      if ( arrival.get_arrival( n ) > arrival.get_required( ncurr ) )
      {
        return max_cost;
      }
    }

    if( ps.preserve_interval )
    {
      auto act = evaluate_desired( ncurr );
      if ( ( arrival.get_time( n ) > arrival.get_time( ncurr ) ) )//|| ( sensing.get_time( n ) < sensing.get_time( ncurr ) ) )
      {
        return max_cost;
      }
    }

    uint32_t const index = node_to_sim[n];
    double switching = 0;
    if ( ps.use_fanouts )
    {
      switching = mlist_sims[index].get_dyn_switch()*( ntk.fanout_size( n ) + ntk.fanout_size( ncurr ) );
    }
    else
    {
      switching = mlist_sims[index].get_dyn_switch();
    }

    return switching;
  }

  template<class Leaves_t>
  cost evaluate( node const& n, Leaves_t const& leaves )
  {
    for( auto l : leaves )
    {
      if( !ntk.is_null( l ) )
      {
        if ( !has_node( l ) )
        {
          if ( VERBOSEg )
            std::cout << "!has " << l << std::endl;
          return 0;
        }
      }
    }
    cost cost_deref = measure_mffc_deref( n, leaves );
    reference_mffc( n, leaves );
    return cost_deref;
  }

  template<class Leaves_t>
  cost measure_mffc_deref( node const& n, Leaves_t const& leaves )
  {
    if ( !has_node( n ) )
    {
      return 0;
    }
    /* reference cut leaves */
    for ( auto leaf : leaves )
    {
      if( ntk.is_null( leaf ) )
        continue;
      if ( !has_node( leaf ) )
      {
        if ( VERBOSEg )
          std::cout << "NOt has " << std::endl;
        return 0;
      }
      ntk.incr_fanout_size( leaf );
    }

    cost mffc_cost = 0;

    if ( ps.use_fanouts )
    {
      mffc_cost = nodes_sims[node_to_sim[n]].get_dyn_switch()*ntk.fanout_size( n );
    }
    else
    {
      mffc_cost = nodes_sims[node_to_sim[n]].get_dyn_switch();
    }
    mffc_cost += recursive_deref( n );

    /* dereference leaves */
    for ( auto leaf : leaves )
    {
      if( ntk.is_null( leaf ) )
        continue;
      ntk.decr_fanout_size( leaf );
    }

    return mffc_cost;
  }

  cost recursive_deref( node const& n )
  {
    /* terminate? */
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      return 0;

    /* recursively collect nodes */
    cost value = 0;

    ntk.foreach_fanin( n, [&]( auto const& fi ) {
      auto ni = ntk.get_node( fi );
      if ( ps.use_fanouts && !ntk.is_pi( ni ) )
      {
        value += nodes_sims[node_to_sim[ni]].get_dyn_switch();
      }
      if ( ntk.decr_fanout_size( ni ) == 0 )
      {
        value += recursive_deref( ni );
      }
    } );
    return value;
  }

  template<class Leaves_t>
  void reference_mffc( node const& n, Leaves_t const& leaves )
  {
    /* reference cut leaves */
    for ( auto leaf : leaves )
    {
      if( ntk.is_null( leaf ) )
        continue;
      ntk.incr_fanout_size( leaf );
    }

    recursive_ref( n );

    /* dereference leaves */
    for ( auto leaf : leaves )
    {
      if( ntk.is_null( leaf ) )
        continue;
      ntk.decr_fanout_size( leaf );
    }

    return;
  }

  void recursive_ref( node const& n )
  {
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      return;
    /* recursively collect nodes */
    ntk.foreach_fanin( n, [&]( auto const& s ) {
      if ( ntk.incr_fanout_size( ntk.get_node( s ) ) == 0 )
      {
        recursive_ref( ntk.get_node( s ) );
      }
    } );
  }

  double compare( cost const& current, cost const& candidate )
  {
    return current - candidate;
  }

  double get_time( node const& n )
  {
    return arrival.get_time( n );
  }

  template<typename Fn>
  void foreach_gate( Fn&& fn ) const
  {
    for ( int i = 0; i < nodes.size(); i++ )
    {
      fn( nodes[i] );
    }
  }

  bool is_interesting( node n )
  {
    return true;
  }

  void choose_gates()
  {

    if( ps.only_critical )
    {
      choose_hot_gates();
    }
    else
    {
      choose_simple();
    }
  }

  void choose_hot_gates()
  {

    convo_analyzer<Ntk, num_steps, cube_size, max_leaves> power( ntk );
    power( 1 );

    uint32_t const size = ntk.size();

    std::vector<std::pair<node, double>> nodes_with_hotness;

    ntk.foreach_gate( [&]( auto n, auto i ) {
      if( i > size )
      {
        return false;
      }
      if( ntk.is_pi( n ) || ntk.is_constant( n ) )
      {
        return true;
      }

      if ( ntk.is_dead( n ) )
      {
        return true;
      }
      
      if( ps.only_critical )
      {
        double cost = 0;
        cost += ps.use_fanouts ? power.get_dyn_switch( n )*ntk.fanout_size( n ) : power.get_dyn_switch( n );
        nodes_with_hotness.push_back( std::make_pair( n, cost ) );
      }
      return true;
    } );
    std::sort( nodes_with_hotness.begin(), nodes_with_hotness.end(), []( auto a, auto b ){
      return a.second > b.second;
    } );

    for( auto n : nodes_with_hotness )
      nodes.push_back( n.first );
  }

  void choose_simple()
  {

    uint32_t const size = ntk.size();
    ntk.foreach_gate( [&]( auto n, auto i ) {
      if( i > size )
      {
        return false;
      }
      if( ntk.is_pi( n ) || ntk.is_constant( n ) )
      {
        return true;
      }

      if ( ntk.is_dead( n ) )
      {
        return true;
      }
      
      nodes.push_back( n );
      return true;
    } );
  }

public:

private:
  Ntk & ntk;
  bindings const lib;

  std::array<node_simulator_t, max_leaves> input_sims;
  std::vector<node_simulator_t> nodes_sims;
  std::vector<node_simulator_t> mlist_sims;
  incomplete_node_map<uint32_t, Ntk> node_to_sim;
  sensing_analyzer<Ntk> sensing;
  arrival_analyzer<Ntk> arrival;
  std::vector<typename Ntk::node> nodes;
  std::set<typename Ntk::node> window_nodes;

  double dyn_switch;
  double zero_delay;
  list_simulator_t list_simulator;
  analyzer_params ps;

};

} /* namespace mockturtle */