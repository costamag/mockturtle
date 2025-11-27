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
  \file timing.hpp
  \brief Compute timing information of a network.

  This engine can be used for power analysis of mapped network.
  For each node, the following information is stored:
  - The sensing time : the first time at which a transition can happen
  - The cumwire time : the first time at which the output is stable
  - A vector of simulation patterns identifying quantized timesteps in this interval
  The class has the following template parameters:
  \param N : number of timesteps used for simulating each signal
  \param I : number of inputs to be used for simulating the static truth table ( 2^I input pairs )

  \author Andrea Costamagna
*/

#pragma once

#include "../../utils/node_map.hpp"
#include <set>

namespace mockturtle
{

template<class Ntk>
class cumwire_analyzer
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

  cumwire_analyzer( Ntk & ntk )
    : ntk( ntk ),
      cumwires( ntk ),
      required( ntk )
  {
    rev_topo_order.reserve( ntk.size() );
    compute_cumwire_times();
    compute_required_times();
    add_event = ntk.events().register_add_event( [&]( const auto& n ) {
      cumwires.resize();
      signal const f = ntk.make_signal( n );
//      compute_cumwire_tfi( f );
//      compute_cumwire_times();
    //  compute_cumwire_times();
    //  compute_required_times();

    } );

    modified_event = ntk.events().register_modified_event( [&]( const auto& n, auto children ) {
      cumwires.resize();
      ntk.incr_trav_id();
      //compute_cumwire_tfo( n );
      //if ( input_cumwires_.size() > 0 )
      //  compute_cumwire_times( input_cumwires_ );
      //else
      compute_cumwire_times();
      compute_required_times();
    } );
  }

  cumwire_analyzer( Ntk & ntk, std::vector<double> const& pattern )
    : ntk( std::make_shared<Ntk>( ntk ) ),
      cumwires( ntk ),
      required( ntk ),
      input_cumwires_( pattern )
  {
    rev_topo_order.reserve( ntk.size() );
    compute_cumwire_times( input_cumwires_ );
    compute_required_times();
    add_event = ntk.events().register_add_event( [&]( const auto& n ) {
      cumwires.resize();
      signal const f = ntk.make_signal( n );
      //compute_cumwire_times();
      //compute_required_times();
    } );

    modified_event = ntk.events().register_modified_event( [&]( const auto& n, auto children ) {
      cumwires.resize();
      ntk.incr_trav_id();
      //compute_cumwire_tfo( n );
      compute_cumwire_times();
      compute_required_times();
    } );
  }

  ~cumwire_analyzer()
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

#pragma region Arrival
void compute_cumwire_times( std::vector<double> const& input_cumwires )
{
  input_cumwires_ = input_cumwires;
  ntk.incr_trav_id();
  cumwires.reset();
  int index = 0;
  ntk.foreach_pi( [&]( auto n ){
    cumwires[n] = input_cumwires[index++];
  } );
  ntk.foreach_po( [&]( auto f ){
    compute_cumwire_tfi( f );
  } );
}

void compute_cumwire_times()
{
  cumwires.reset();
  std::vector<double> input_cumwires( ntk.num_pis(), 0 );
  if ( input_cumwires_.size() > 0 )
    input_cumwires = input_cumwires_;
  compute_cumwire_times( input_cumwires );
}

void topo_sort()
{
  rev_topo_order.clear();
  ntk.incr_trav_id();

  std::set<node> outputs;
  ntk.foreach_po( [&]( auto fo, auto io ){
    auto no = ntk.get_node( fo );
    outputs.insert( no );
  } );
  while( outputs.size() > 0 )
  {
    for( auto no : outputs )
    {
      bool dangling = true;
      ntk.foreach_fanout( no, [&]( auto fo, auto io ){
        if( outputs.find( fo ) != outputs.end() )
        {
          dangling = false;
        }
      } );
      if( dangling )
      {
        topo_sort_rec( ntk.make_signal( no ) );
        outputs.erase( no );
        break;
      }
    }
  }
  std::reverse( rev_topo_order.begin(), rev_topo_order.end() );
}

void topo_sort_rec( signal f )
{
  node n = ntk.get_node( f );
  if( ntk.visited( n ) == ntk.trav_id() )
    return;
  if( ntk.is_pi( n ) )
    return;
  
  ntk.set_visited( n, ntk.trav_id() );
  ntk.foreach_fanin( n, [&]( auto fi, auto ii ){
    topo_sort_rec( fi );
  } );
  rev_topo_order.push_back( n );
}

void compute_required_times()
{
  topo_sort();
  required.reset();

  ntk.foreach_node( [&]( auto no, auto io ){
    required[no] = worst_delay();
  } );

  for( auto n : rev_topo_order )
  {
    ntk.foreach_fanout( n, [&]( auto fo, auto io ){
      node no = ntk.get_node( fo );
//      std::cout << "doing " << no << std::endl;
      if( !ntk.is_dead( no ) )
      {
        auto go = ntk.get_binding( no );
        uint32_t index = -1;
        ntk.foreach_fanin( no, [&]( auto fi, auto ii ){
          if( ntk.get_node( fi ) == n )
          {
            index = ii;
          }
        } );
//        assert( index != -1 && "Node not found in the FIs of the FO node" );
        if( index != -1 )
        {
          double const rise_time = go.pins[index].rise_block_delay;
          double const fall_time = go.pins[index].fall_block_delay;
          double const delta = std::max( rise_time, fall_time );
  //        std::cout << delta << " delta" << std::endl;
  //        std::cout << required[no] << " required[" << no << "]" << std::endl;
  //        std::cout << no << " " << required[no] - delta << std::endl;
          required[n] = std::min( required[n], required[no] - delta );
        }
//        std::cout << "required[" << n << "]= " << required[n] << std::endl;
//        assert( ( required[no] >= delta ) );
      }
    } );
  }
}

void compute_cumwire_tfi( signal const& f )
{
  node n = ntk.get_node( f );

  if( n > ntk.size() )
  {
    return;
  }

  if( ntk.is_pi( n ) )
  {
    return;
  }
  if( ntk.visited( n ) == ntk.trav_id() )
  {
    return;
  }
  auto g = ntk.get_binding( n );
  cumwires[n] = 0;
  ntk.foreach_fanin( n, [&]( auto fi, auto ii ){
    node ni = ntk.get_node( fi );
    assert( ( ni != n ) );
    compute_cumwire_tfi( fi );
    double const rise_time = g.pins[ii].rise_block_delay;
    double const fall_time = g.pins[ii].fall_block_delay;
    double const max_delay = std::max( rise_time, fall_time );
    cumwires[n] = std::max( cumwires[n], cumwires[ni] + max_delay );    
  } );
  ntk.set_visited( n, ntk.trav_id() );
}

void compute_cumwire_tfo( node const& n )
{
  if( ntk.visited( n ) == ntk.trav_id() )
    return;
  ntk.set_visited( n, ntk.trav_id() );

  auto g = ntk.get_binding( n );
  cumwires[n] = 0;
  ntk.foreach_fanin( n, [&]( auto fi, auto ii ){
    node ni = ntk.get_node( fi );
    double const rise_time = g.pins[ii].rise_block_delay;
    double const fall_time = g.pins[ii].fall_block_delay;
    double const max_delay = std::max( rise_time, fall_time );
    cumwires[n] = std::max( cumwires[n], cumwires[ni] + max_delay );    
  } );

  ntk.foreach_fanout( n, [&]( auto no ){
    compute_cumwire_tfo( no );
  } );

}

double worst_delay()
{
  double worst_delay = 0;
  ntk.foreach_po( [&]( auto f ){
    worst_delay = std::max( worst_delay, cumwires[f] );
  } );
  return worst_delay;
}

double get_time( node n )
{
  if( ntk.is_null( n ) ) 
    return 1.0*(1u << 20u);

  if( !cumwires.has( n ) )
    compute_cumwire_times();//_tfi( ntk.make_signal( n ) );

  return cumwires[n];
}

double get_time( signal f )
{
  node const n = ntk.get_node( f );
  return get_time( n );
}

double get_required( signal f )
{
  node const n = ntk.get_node( f );
  return required[n];
}

double get_required( node n )
{
  if ( !required.has( n ) )
    compute_required_times();
  if ( !required.has( n ) )
    return 0;
  return required[n];
}

double get_cumwire( signal f )
{
  node const n = ntk.get_node( f );
  return cumwires[n];
}

bool has( node n ) const
{
  return cumwires.has( n )&&required.has(n);
}

bool is_critical( node n ) const
{
  return std::abs( required[n] - cumwires[n] ) <= 0.01;
}

#pragma endregion



private:
  Ntk & ntk;
  incomplete_node_map<double, Ntk> cumwires;
  incomplete_node_map<double, Ntk> required;
  std::vector<double> input_cumwires_;
  std::vector<node> rev_topo_order;
  /* events */
  std::shared_ptr<typename network_events<Ntk>::add_event_type> add_event;
  std::shared_ptr<typename network_events<Ntk>::modified_event_type> modified_event;
};

} // namespace mockturtle