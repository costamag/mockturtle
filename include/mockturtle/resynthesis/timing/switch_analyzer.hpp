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
  - The arrival time : the first time at which the output is stable
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
class switch_analyzer
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

  switch_analyzer( Ntk & ntk )
    : pntk( std::make_shared<Ntk>( ntk ) ),
      flow( ntk )
  {
    add_event = pntk->events().register_add_event( [&]( const auto& n ) {
      flow.resize();
      signal const f = pntk->make_signal( n );
      compute_flow_tfi( f );
    } );

    modified_event = pntk->events().register_modified_event( [&]( const auto& n, auto children ) {
      flow.resize();
      pntk->incr_trav_id();
      //compute_arrival_tfo( n );
      compute_flow();
    } );
  }

  ~switch_analyzer()
  {
    if ( add_event )
    {
      pntk->events().release_add_event( add_event );
    }

    if ( modified_event )
    {
      pntk->events().release_modified_event( modified_event );
    }
  }

#pragma region Arrival
void compute_flow( std::vector<double> const& input_switches )
{
  pntk->incr_trav_id();
  flow.reset();
  int index = 0;
  pntk->foreach_pi( [&]( auto n ){
    flow[n] = input_switches[index++];
  } );
  pntk->foreach_po( [&]( auto f ){
    compute_flow_tfi( f );
  } );
}

void compute_flow()
{
  flow.reset();
  std::vector<double> input_switches( pntk->num_pis(), 0 );
  compute_flow( input_switches );
}

void compute_flow_tfi( signal const& f )
{
  node n = pntk->get_node( f );

  if( n > pntk->size() )
  {
    return;
  }

  if( pntk->is_constant( n ) )
  {
    flow[n] = 0;
    return;
  }

  if( pntk->is_pi( n ) )
  {
    flow[n] = 0.5;
    return;
  }
  if( pntk->visited( n ) == pntk->trav_id() )
  {
    return;
  }
  flow[n] = pntk->get_power( n );

  pntk->foreach_fanin( n, [&]( auto fi, auto ii ){
    node ni = pntk->get_node( fi );
    compute_flow_tfi( fi );
    flow[n] += flow[ni];    
  } );
  pntk->set_visited( n, pntk->trav_id() );
}

double get_flow( node n )
{
  if( pntk->is_null( n ) ) 
    return 1.0*(1u << 20u);

  if( !flow.has( n ) )
    compute_flow();//_tfi( pntk->make_signal( n ) );

  return flow[n];
}

bool has( node n ) const
{
  return flow.has( n );
}

#pragma endregion

private:
  std::shared_ptr<Ntk> pntk;
  incomplete_node_map<double, Ntk> flow;
  /* events */
  std::shared_ptr<typename network_events<Ntk>::add_event_type> add_event;
  std::shared_ptr<typename network_events<Ntk>::modified_event_type> modified_event;
};

} // namespace mockturtle