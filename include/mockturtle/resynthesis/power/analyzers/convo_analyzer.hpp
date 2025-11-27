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
  \file convo_analyzer.hpp
  \brief Analyze power including dyn_switch.

  This engine can be used for power analysis of mapped network.
  For each node, the following information is stored:
  - The sensing time : the first time at which a transition can happen
  - The arrival time : the first time at which the output is stable
  - A vector of simulation patterns identifying quantized timesteps in this interval
  The power of a node is evaluated by computing a reconvergence driven cut and simulating it.
  The class has the following template parameters:
  \param N : number of timesteps used for simulating each signal
  \param I : number of inputs to be used for simulating the static truth table ( 2^I input pairs )

  \author Andrea Costamagna
*/

#pragma once

#include "../node_simulator.hpp"
#include "../../list_simulator.hpp"
#include "../../timing/arrival_analyzer.hpp"
#include "../../timing/sensing_analyzer.hpp"
#include "../../windowing/window_manager.hpp"
#include "../../../algorithms/detail/resub_utils.hpp"
#include "../../../algorithms/reconv_cut.hpp"
#include "../../../algorithms/simulation.hpp"

namespace mockturtle
{

template<class Ntk, uint32_t N = 16, uint32_t I = 12, uint32_t L = 14>
class convo_analyzer
{
public:
  using truth_table_t = kitty::static_truth_table<I>;
  using signal = typename Ntk::signal;
  using node = typename Ntk::node;
  using node_simulator_t = node_simulator<Ntk, N, I>;
  using list_simulator_t = aig_list_simulator<truth_table_t>;

  convo_analyzer( Ntk& ntk )
      : ntk( ntk ),
        win_manager( ntk ),
        node_to_sim( ntk ),
        node_to_zero_delay( ntk ),
        node_to_dyn_switch( ntk ),
        node_to_dyn_power( ntk ),
        node_to_count( ntk ),
        sensing( ntk ),
        arrival( ntk ),
        gen( 336 )
  {
    /* allocate enough memory for storing the simulation information */
    size_t extra_space = 1000;
    nodes_sims.reserve( ntk.size() + extra_space );
    init();

    /* define the events */
    add_event = ntk.events().register_add_event( [&]( const auto& n ) {
      node_to_sim.resize();
      node_to_dyn_switch.resize();
      node_to_dyn_power.resize();
      node_to_zero_delay.resize();
      node_to_count.resize();
      ntk.incr_trav_id();
      compute_power();
    } );

    deleted_event = ntk.events().register_delete_event( [&]( const auto& n ) {
    } );

    modified_event = ntk.events().register_modified_event( [&]( const auto& n, auto children ) {
      node_to_sim.resize();
      node_to_dyn_switch.resize();
      node_to_dyn_power.resize();
      node_to_zero_delay.resize();
      node_to_count.resize();
      ntk.incr_trav_id();
      compute_power();
    } );

  }

  ~convo_analyzer()
  {
    if ( add_event )
    {
      ntk.events().release_add_event( add_event );
    }

    if ( modified_event )
    {
      ntk.events().release_modified_event( modified_event );
    }

    if ( deleted_event )
    {
      ntk.events().release_delete_event( deleted_event );
    }
  }

  void init()
  {
    std::vector<double> input_sensing( ntk.num_pis(), +0.0 );
    std::vector<double> input_arrival( ntk.num_pis(), +0.0 );
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

    node_to_sim[ntk.get_node( ntk.get_constant( false ) )] = 0;
    node_to_sim[ntk.get_node( ntk.get_constant( true ) )] = 1;
    ntk.set_visited( 0, ntk.trav_id() );
    ntk.set_visited( 1, ntk.trav_id() );
  }

  void operator()( uint32_t num_workloads = 1 )
  {
    std::array<kitty::static_truth_table<I>, 2> workload;
    for( int w = 0; w < num_workloads; ++w )
    {
      reset();

      ntk.incr_trav_id();
      ntk.set_visited( ntk.get_node( ntk.get_constant( false ) ), ntk.trav_id() );
      ntk.set_visited( ntk.get_node( ntk.get_constant( true ) ), ntk.trav_id() );

      std::default_random_engine::result_type seed1 = 5;
      std::default_random_engine::result_type seed2 = 8;

      /* Initialize the PIs */
      int k = 1;
      for( int i = 0; i < L; ++i )
      {
        kitty::create_random( workload[0], seed1*( k ) );
        kitty::create_random( workload[1], seed2*( k+1 ) );
        node_simulator_t sim( 0, workload, 0, 0 );
        nodes_sims.emplace_back( sim );
        k += 2;
      }
      compute_power();
    }
  }

  void compute_power()
  {
    ntk.foreach_gate( [&]( auto n ){
      node_to_zero_delay[n] = 0;
      node_to_dyn_switch[n] = 0;
      node_to_dyn_power[n] = 0;
      node_to_count[n] = 0;
    } );

    ntk.foreach_gate( [&]( auto n ){
      if( !ntk.is_dead( n ) && win_manager.run( n ) )
      {
        node_to_sim.reset();
        nodes_sims.resize( 2 + L );
        node_to_sim[ntk.get_node( ntk.get_constant( false ) )] = 0;
        node_to_sim[ntk.get_node( ntk.get_constant( true ) )] = 1;
  
        std::set<node> visited;
        visited.insert( 0 );
        visited.insert( 1 );
        for( int i = 0; i < win_manager.leaves.size() ; ++i )
        {
          node const l = win_manager.leaves[i];
          node_to_sim[l] = i + 2;
          nodes_sims[i+2].update_info( l, sensing.get_time( l ), arrival.get_time( l ) );
          visited.insert( l );
        }
  
        std::vector<node_simulator_t*> sim_ptrs;
        for( int i = 0; i < win_manager.divs.size() ; ++i )
        {
          if( i < win_manager.leaves.size() )
          {
            continue; // this is a leaf node
          }
  
          node const d = win_manager.divs[i];
          sim_ptrs.clear();
          ntk.foreach_fanin( d, [&]( auto f, auto j ){
            node const ni = ntk.get_node( f );
            if( visited.find( ni ) == visited.end() )
              std::cerr << "not simulated" << std::endl;
            assert( visited.find( ni ) != visited.end() );
            sim_ptrs.push_back( &nodes_sims[node_to_sim[ni]] );
          } );
          uint32_t const index { nodes_sims.size() };
          node_to_sim[d] = index;
          nodes_sims.emplace_back( d, sensing.get_time( d ), arrival.get_time( d ) );
          auto const bind = ntk.get_binding( d );
          nodes_sims[index]( bind, list_simulator, sim_ptrs );   
          visited.insert( d ); 
          node_to_zero_delay[n] += nodes_sims[index].get_zero_delay();
          node_to_dyn_switch[n] += nodes_sims[index].get_dyn_switch();
          node_to_dyn_power[n] += node_to_dyn_switch[n]*ntk.fanout_size( n );
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
  
        uint32_t const index = node_to_sim[n];
      }
      else if ( !ntk.is_dead( n ) )
      {
        std::cout << "[w] window construction failed" << std::endl;
      }

      //return;
    } );
  }

  void reset()
  {
    /* clear information */
    nodes_sims.resize( 2 );
    node_to_sim.reset();
    node_to_sim[ntk.get_node( ntk.get_constant( false ) )] = 0;
    node_to_sim[ntk.get_node( ntk.get_constant( true ) )] = 1;
  }

  void compute_power_tfi( signal f )
  {
    node n = ntk.get_node( f );
    if ( ( ntk.visited( n ) == ntk.trav_id() ) || ( ntk.is_pi( n ) ) )
    {
      return;
    }
    ntk.set_visited( n, ntk.trav_id() );

    node_simulator_t sim( n, sensing.get_time( n ), arrival.get_time( n ) );

    std::vector<node_simulator_t*> sims_ptrs;
    ntk.foreach_fanin( n, [&]( auto fi, auto ii ) {
      compute_power_tfi( fi );
      node const ni = ntk.get_node( fi );
      uint32_t const index = node_to_sim[ni];
      sims_ptrs.push_back( &nodes_sims[index] );
    } );
    sim( ntk, list_simulator, sims_ptrs );
    
    if ( node_to_sim.has( n ) )
    {
      uint32_t const index = node_to_sim[n];
      node_simulator_t& old_sim = nodes_sims[index];
      old_sim = sim;
    }
    else
    {
      uint32_t const index = nodes_sims.size();
      node_to_sim[n] = index;
      nodes_sims.push_back( sim );
    }
  }

  double get_zero_delay( node n ) const
  {
    if( ntk.is_pi( n ) )
      return 0.5;
    if( ntk.is_constant( n ) )
      return 0;
    if( ntk.is_null( n ) )
      return 1u << 10;
    return node_to_zero_delay[n];
  }

  double get_dyn_switch( node n ) const
  {
    if( ntk.is_pi( n ) )
      return 0.5;
    if( ntk.is_constant( n ) )
      return 0;
    if( ntk.is_null( n ) )
      return 1u << 10;
    return node_to_dyn_switch[n];
  }

  double get_dyn_power( node n ) const
  {
    if( ntk.is_null( n ) )
      return 1u << 10;
    return get_dyn_switch( n ) * ntk.fanout_size( n );
  }

  double approximate_dyn_switch( node n )
  {
    auto sim = nodes_sims[node_to_sim[n]];
    return sim.approximate_dyn_switch();
  }

  uint32_t num_steps() const
  {
    return N;
  }

  auto get_probs( node n )
  {
    auto sim = nodes_sims[node_to_sim[n]];
    return sim.get_probs();
  }

  auto get_trans_probs( node n ) const
  {
    auto sim = nodes_sims[node_to_sim[n]];
    return sim.get_trans_probs();
  }

  auto get_arrival( node n ) const
  {
    auto sim = nodes_sims[node_to_sim[n]];
    return sim.get_arrival();
  }

  auto get_sensing( node n )
  {
    auto sim = nodes_sims[node_to_sim[n]];
    return sim.get_sensing();
  }

  void print( node n ) const
  {
    auto const& idx = node_to_sim[n];
    auto const& sim = nodes_sims[idx];
    std::cout << "node " << n << std::endl;
    sim.print();
  }

  void print()
  {
    for( auto const& sim : nodes_sims )
    {
      std::cout << "node " << sim.get_node() << std::endl;
      sim.print();
    }
  }

  void print_power()
  {
    ntk.foreach_gate( [&]( auto n ){
      std::string name = ntk.get_name( n );
      std::cout << n << " " << name << " ";
      ntk.foreach_fanin( n, [&]( auto fi, auto ii ){
        std::cout << fi.index << " ";
      } );
      std::cout << " FO: ";
      ntk.foreach_fanout( n, [&]( auto fo ){
        std::cout << fo << " ";
      } );
      std::cout << get_dyn_switch(n)*ntk.fanout_size( n ) << std::endl;
      return true;
    } );

    ntk.foreach_po( [&]( auto fo ){
      std::cout << fo.index << " ";
    } );
    std::cout << " <- POs" << std::endl;
  }

  node_simulator_t * get_simulator_ptr( node n )
  {
    auto const idx = node_to_sim[n];
    return &nodes_sims[idx];
  }

  auto get_trans_probs( node n )
  {
    auto const idx = node_to_sim[n];
    return nodes_sims[idx].get_trans_probs();
  }

//  double get_dyn_switch( node n )
//  {
//    auto const idx = node_to_sim[n];
//    return nodes_sims[idx].get_dyn_switch();
//  }

  auto get_trans_matrix( node n )
  {
    auto const idx = node_to_sim[n];
    return nodes_sims[idx].get_trans_matrix();
  }

  auto get_sum( node n ) const
  {
    auto const idx = node_to_sim[n];
    return nodes_sims[idx].get_sum();
  }

  uint32_t print_intervals( node n )
  {
    auto const idx = node_to_sim[n];
    return nodes_sims[idx].print_intervals();  
  }

  void print_distr( node n )
  {
    auto const idx = node_to_sim[n];
    nodes_sims[idx].print_distribution();  
  }

  void print_ntk()
  {
    auto binds = ntk.get_bindings();
    ntk.foreach_gate( [&]( auto n ){
      if( ntk.is_dead( n ) )
        return true;
      std::string name = ntk.get_name( n );
      std::cout << n << " ";
      std::cout << " GLITCH " << node_to_dyn_switch[n]- node_to_zero_delay[n] << " ";
      std::cout << " SWITCH " << node_to_dyn_switch[n] << " ";
      std::cout << name << " ";
      ntk.foreach_fanin( n, [&]( auto fi, auto ii ){
        std::cout << fi.index << " ";
      } );
      std::cout << " FO: ";
      ntk.foreach_fanout( n, [&]( auto fo ){
        std::cout << fo << " ";
      } );
      double ta = arrival.get_time( n );
      std::cout << " P0 " << node_to_zero_delay[n] << " PD " << node_to_dyn_switch[n] << " tA " << ta << std::endl;
      return true;
    } );

    ntk.foreach_po( [&]( auto fo ){
      std::cout << fo.index << " ";
    } );
    std::cout << " <- POs" << std::endl;
  }

  void print_worst_glitches( double thresh = 0.3 )
  {
    ntk.foreach_gate( [&]( auto n ){
      double glitch = node_to_dyn_switch[n]- node_to_zero_delay[n];
      if( glitch >= thresh )
      {
        ntk.incr_trav_id();
        printf(" NODE AT %5.2f glitches\n", glitch );
        print_ntk_rec( n );
        printf("\n");
      }
    } );
  }

  void print_ntk_rec( node const n )
  {
    if( ntk.visited( n ) == ntk.trav_id() )
      return;
    if( ntk.is_dead( n ) )
      return;
    if( ntk.is_pi( n ) )
    {
      return;
    }
    ntk.foreach_fanin( n, [&]( auto fi, auto ii ){
      node ni = ntk.get_node( fi );
      print_ntk_rec( ni );
    } );

    auto binds = ntk.get_bindings();

      std::string name = ntk.get_name( n );
      std::cout << n << " ";
      std::cout << " GLITCH " << node_to_dyn_switch[n]- node_to_zero_delay[n] << " ";
      std::cout << " SWITCH " << node_to_dyn_switch[n] << " ";
      std::cout << name << " ";
      ntk.set_visited( n, ntk.trav_id() );
      ntk.foreach_fanin( n, [&]( auto fi, auto ii ){
        std::cout << fi.index << " ";
      } );
      std::cout << " FO: ";
      ntk.foreach_fanout( n, [&]( auto fo ){
        std::cout << fo << " ";
      } );
      double ta = arrival.get_time( n );
      std::cout << " P0 " << node_to_zero_delay[n] << " PD " << node_to_dyn_switch[n] << " tA " << ta << std::endl;
      return;
  }

private:
  Ntk & ntk;

  window_manager<Ntk, L> win_manager;
  std::array<node_simulator_t, L> input_sims;
  std::vector<node_simulator_t> nodes_sims;
  incomplete_node_map<uint32_t, Ntk> node_to_sim;
  incomplete_node_map<double, Ntk> node_to_dyn_switch;
  incomplete_node_map<double, Ntk> node_to_zero_delay;
  incomplete_node_map<double, Ntk> node_to_dyn_power;
  incomplete_node_map<double, Ntk> node_to_count;
  sensing_analyzer<Ntk> sensing;
  arrival_analyzer<Ntk> arrival;
  std::mt19937 gen; // mersenne_twister_engine seeded with rd()
  list_simulator_t list_simulator;
  std::shared_ptr<typename network_events<Ntk>::add_event_type> add_event;
  std::shared_ptr<typename network_events<Ntk>::modified_event_type> modified_event;
  std::shared_ptr<typename network_events<Ntk>::delete_event_type> deleted_event;
};

} // namespace mockturtle