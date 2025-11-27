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
  \file power_analyzer.hpp
  \brief Analyze power including dyn_switch.

  This engine can be used for power analysis of mapped network.
  For each node, the following information is stored:
  - The sensing time : the first time at which a transition can happen
  - The arrival time : the first time at which the output is stable
  - A vector of simulation patterns identifying quantized timesteps in this interval
  The class has the following template parameters:
  \param N : number of timesteps used for simulating each signal
  \param I : number of inputs to be used for simulating the static truth table ( 2^I input pairs )

  \author Andrea Costamagna
*/

#pragma once

#include "../../../algorithms/simulation.hpp"
#include "../../list_simulator.hpp"
#include "../../timing/sensing_analyzer.hpp"
#include "../../timing/arrival_analyzer.hpp"
#include "../node_simulator.hpp"

namespace mockturtle
{

template<class Ntk, uint32_t N = 16, uint32_t I = 10>
class power_analyzer
{
public:
  using truth_table_t = kitty::static_truth_table<I>;
  using signal = typename Ntk::signal;
  using node = typename Ntk::node;
  using node_simulator_t = node_simulator<Ntk, N, I>;
  using list_simulator_t = aig_list_simulator<truth_table_t>;

  power_analyzer( Ntk& ntk )
      : ntk( ntk ),
        nodes_sims(),
        node_to_sim( ntk ),
        sensing( ntk ),
        arrival( ntk ),
        zero_delay( 0 ),
        dyn_switch( 0 ),
        list_simulator()
  {
    /* allocate enough memory for storing the simulation information */
    size_t extra_space = 1000;
    nodes_sims.reserve( ntk.size() + extra_space );
    topo_sort.reserve( ntk.size() + extra_space );
    _init();

    /* define the events */
    add_event = ntk.events().register_add_event( [&]( const auto& n ) {
      node_to_sim.resize();
      ntk.incr_trav_id();
      ntk.foreach_po( [&]( auto f ) {
        compute_power_tfi( f );
      } );
    } );

    modified_event = ntk.events().register_modified_event( [&]( const auto& n, auto children ) {
      node_to_sim.resize();
      ntk.incr_trav_id();
      ntk.foreach_po( [&]( auto f ) {
        compute_power_tfi( f );
      } );
    } );
  }

  void _init()
  {
  }

  void operator()( std::vector<node_simulator_t> const& input_sims )
  {
    nodes_sims.clear();
    node_to_sim.reset();
    dyn_switch = 0;
    zero_delay = 0;
    dyn_power = 0;
    ntk.incr_trav_id();

    assert( ntk.num_pis() == input_sims.size() );
    std::vector<double> input_sensing( ntk.num_pis() );
    std::vector<double> input_arrival( ntk.num_pis() );
    std::transform( input_sims.begin(), input_sims.end(), input_sensing.begin(), [&]( auto s ) {
      return s.get_sensing();
    } );
    std::transform( input_sims.begin(), input_sims.end(), input_arrival.begin(), [&]( auto s ) {
      return s.get_arrival();
    } );
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

    node_to_sim[ntk.get_node( ntk.get_constant( false ) )] = 0;
    node_to_sim[ntk.get_node( ntk.get_constant( true ) )] = 1;
    ntk.set_visited( 0, ntk.trav_id() );
    ntk.set_visited( 1, ntk.trav_id() );

    /* initialize the PIs */
    int index = 0;
    ntk.incr_trav_id();
    ntk.foreach_pi( [&]( auto n ) {
      node_to_sim[n] = nodes_sims.size();
      nodes_sims.emplace_back( input_sims[index++] );
      ntk.set_visited( n, ntk.trav_id() );
    } );
    /* compute the power for all the nodes */
    sort_topologically();
    propagate();
    compute_power();
    //print();
  }

  ~power_analyzer()
  {
    if ( add_event )
    {
      ntk.events().release_add_event( add_event );
    }

    if ( modified_event )
    {
      ntk.events().release_modified_event( modified_event );
    }
  }

  void sort_topologically( signal f )
  {
    node n = ntk.get_node( f );
    if ( ( ntk.visited(n) == ntk.trav_id() ) || ntk.is_pi( n ) )
      return;

    ntk.foreach_fanin( n, [&]( auto fi, auto ii ){
      sort_topologically( fi );
    } );
    topo_sort.push_back( n );
    ntk.set_visited( n, ntk.trav_id() );
  }

  void sort_topologically()
  {
    topo_sort.clear();
    ntk.incr_trav_id();
    ntk.foreach_po( [&]( auto fo ){
      sort_topologically( fo );
    } );
  }

  void propagate()
  {
    for ( auto n : topo_sort )
    {
      if ( ntk.is_dead( n ) )
        continue;
      if( !arrival.has( n ) )
      {
        arrival.compute_arrival_tfi( ntk.make_signal( n ) );
      }
      node_simulator_t sim( n, sensing.get_time( n ), arrival.get_time( n ) );

      std::vector<node_simulator_t*> sims_ptrs;
      ntk.foreach_fanin( n, [&]( auto fi, auto ii ) {
        node const ni = ntk.get_node( fi );
        uint32_t const index = node_to_sim[ni];
        sims_ptrs.push_back( &nodes_sims[index] );
      } );
      sim( ntk, list_simulator, sims_ptrs );
      dyn_switch += sim.get_dyn_switch();
      zero_delay += sim.get_zero_delay();

      if ( node_to_sim.has( n ) )
      {
        uint32_t const index = node_to_sim[n];
        node_simulator_t& old_sim = nodes_sims[index];
        dyn_switch -= old_sim.get_dyn_switch();
        zero_delay -= old_sim.get_zero_delay();
        old_sim = sim;
      }
      else
      {
        uint32_t const index = nodes_sims.size();
        node_to_sim[n] = index;
        nodes_sims.push_back( sim );
      }
    }
  }

  void compute_power_tfi( signal f )
  {
    node n = ntk.get_node( f );
    //std::cout << n << ": ";
    ntk.foreach_fanin( n, [&]( auto fi, auto ii ) {
      node const ni = ntk.get_node( fi );
      //std::cout << ni << " ";
    } );
    //std::cout << std::endl;

    if ( ( node_to_sim.has( n ) ) || ( ntk.visited( n ) == ntk.trav_id() ) || ( ntk.is_pi( n ) ) )
    {
      return;
    }
    if( !arrival.has( n ) )
    {
      arrival.compute_arrival_tfi( ntk.make_signal( n ) );
    }
    node_simulator_t sim( n, sensing.get_time( n ), arrival.get_time( n ) );
    ntk.set_visited( n, ntk.trav_id() );

    std::vector<node_simulator_t*> sims_ptrs;
    ntk.foreach_fanin( n, [&]( auto fi, auto ii ) {
      compute_power_tfi( fi );
      node const ni = ntk.get_node( fi );
      uint32_t const index = node_to_sim[ni];
      sims_ptrs.push_back( &nodes_sims[index] );
    } );
//if( sims_ptrs.size() == 1 && n == 17 )
//{
//  sims_ptrs[0]->print();
//  auto g  = ntk.get_binding( n );
//  kitty::print_binary( g.function );
//  std::cout << std::endl; 
//}
    sim( ntk, list_simulator, sims_ptrs );
//if( sims_ptrs.size() == 1 && n == 17 )
//{
//  ntk.foreach_fanin( n, [&]( auto fi, auto ii ) {
//    node const ni = ntk.get_node( fi );
//    std::cout << ni << ".fi " << std::endl;
//  } );
//  std::cout << "dyn switch " << sim.get_dyn_switch() << std::endl; 
//}
    dyn_switch += sim.get_dyn_switch();
    zero_delay += sim.get_zero_delay();
    
    if ( node_to_sim.has( n ) )
    {
      uint32_t const index = node_to_sim[n];
      node_simulator_t& old_sim = nodes_sims[index];
      dyn_switch -= old_sim.get_dyn_switch();
      zero_delay -= old_sim.get_zero_delay();
      old_sim = sim;
    }
    else
    {
      uint32_t const index = nodes_sims.size();
      node_to_sim[n] = index;
      nodes_sims.push_back( sim );
    }
  }

  void compute_power_tfo( node n )
  {
    if ( ntk.visited( n ) == ntk.trav_id() )
    {
      return;
    }

    node_simulator_t sim( n, sensing.get_time( n ), arrival.get_time( n ) );

    std::vector<node_simulator_t*> sims_ptrs;
    ntk.foreach_fanin( n, [&]( auto fi, auto ii ) {
      compute_power_tfi( fi );
      node const ni = ntk.get_node( fi );
      uint32_t const index = node_to_sim[ni];
      sims_ptrs.push_back( &nodes_sims[index] );
    } );
    sim( ntk, list_simulator, sims_ptrs );
    dyn_switch += sim.get_dyn_switch();
    zero_delay += sim.get_zero_delay();
    if ( node_to_sim.has( n ) )
    {
      uint32_t const index = node_to_sim[n];
      node_simulator_t& old_sim = nodes_sims[index];
      dyn_switch -= old_sim.get_dyn_switch();
      zero_delay -= old_sim.get_zero_delay();
      old_sim = sim;
    }
    else
    {
      uint32_t const index = nodes_sims.size();
      node_to_sim[n] = index;
      nodes_sims.push_back( sim );
    }
    ntk.set_visited( n, ntk.trav_id() );
    ntk.foreach_fanout( n, [&]( auto no ) {
      compute_power_tfo( no );
    } );
  }

  double get_zero_delay() const
  {
    return zero_delay;
  }

  double get_dyn_switch() const
  {
    return dyn_switch;
  }

  double get_zero_delay( node n ) const
  {
    auto sim = nodes_sims[node_to_sim[n]];
    return sim.get_zero_delay();
  }

  double get_dyn_switch( node n ) const
  {
    auto sim = nodes_sims[node_to_sim[n]];
    return sim.get_dyn_switch();
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

  double get_dyn_switch( node n )
  {
    auto const idx = node_to_sim[n];
    return nodes_sims[idx].get_dyn_switch();
  }

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


  void compute_power()
  {
    bool failed = false;
    dyn_power = 0;
    zero_delay = 0;
    dyn_switch = 0;

    for ( auto n : topo_sort ){
      if( failed )
      {
        continue;
      }
      if( ntk.is_pi( n ) || ntk.is_constant( n ) )
      {
        continue;
      }
      if( ntk.is_dead( n ) )
      {
        continue;
      }
      dyn_power += nodes_sims[node_to_sim[n]].get_dyn_switch() * fanout_size( n );
// std::cout << n << " " << dyn_power << std::endl;
      zero_delay += nodes_sims[node_to_sim[n]].get_zero_delay();
      dyn_switch += nodes_sims[node_to_sim[n]].get_dyn_switch();
    }
  }
public:
  double get_dyn_power() const
  {
    return dyn_power;
  }

  void print_power()
  {
    ntk.foreach_gate( [&]( auto n ){
      if ( ntk.is_dead( n ) || ntk.is_pi( n ) )
        return true;
      std::string name = ntk.get_name( n );
      std::cout << n << " " << name << " ";
      ntk.foreach_fanin( n, [&]( auto fi, auto ii ){
        std::cout << fi.index << " ";
      } );
      std::cout << " FO: {";
      ntk.foreach_fanout( n, [&]( auto fo ){
        std::cout << fo << " ";
      } );
      std::cout << "} DS = {" << get_dyn_switch(n) << 
                   "} ZD = {" << get_zero_delay(n) << 
                   "} PW = {" << get_dyn_switch(n)*fanout_size( n ) << "} tA=" << get_arrival( n ) << " tS=" << get_sensing( n ) << std::endl;
      return true;
    } );

    ntk.foreach_po( [&]( auto fo ){
      std::cout << fo.index << " ";
    } );
    std::cout << " <- POs" << std::endl;
  }

  uint32_t fanout_size( node const& n )
  {
    auto const fo_size = ntk.fanout_size( n ); 
    auto const& fanouts = ntk._storage->nodes[n].fanout;
    if ( fo_size != fanouts.size() )
    {
      assert( ntk.is_po( n ) );
    }
    return fo_size;
  }

private:
  Ntk & ntk;
  std::vector<node_simulator_t> nodes_sims;
  incomplete_node_map<uint32_t, Ntk> node_to_sim;
  sensing_analyzer<Ntk> sensing;
  arrival_analyzer<Ntk> arrival;
  std::vector<node> topo_sort;
  double dyn_switch;
  double zero_delay;
  double dyn_power;
  list_simulator_t list_simulator;

  /* events */
  std::shared_ptr<typename network_events<Ntk>::add_event_type> add_event;
  std::shared_ptr<typename network_events<Ntk>::modified_event_type> modified_event;
};

} // namespace mockturtle