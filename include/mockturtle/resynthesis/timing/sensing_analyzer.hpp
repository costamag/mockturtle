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
  - The sensing time : the first time at which the output is stable
  - A vector of simulation patterns identifying quantized timesteps in this interval
  The class has the following template parameters:
  \param N : number of timesteps used for simulating each signal
  \param I : number of inputs to be used for simulating the static truth table ( 2^I input pairs )

  \author Andrea Costamagna
*/

#pragma once

#include "../../utils/node_map.hpp"


namespace mockturtle
{

template<class Ntk>
class sensing_analyzer
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

  sensing_analyzer( Ntk & ntk )
    : pntk( std::make_shared<Ntk>( ntk ) ),
      sensings( ntk )
  {
    add_event = pntk->events().register_add_event( [&]( const auto& n ) {
      sensings.resize();
      signal const f = pntk->make_signal( n );
      compute_sensing_tfi( f );
    } );

    modified_event = pntk->events().register_modified_event( [&]( const auto& n, auto children ) {
      sensings.resize();
      pntk->incr_trav_id();
      compute_sensing_times();
    } );

  }

  ~sensing_analyzer()
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
void compute_sensing_times( std::vector<double> const& input_sensings )
{
  sensings.reset();
  int index = 0;
  pntk->foreach_pi( [&]( auto n ){
    sensings[n] = input_sensings[index++];
  } );
  pntk->foreach_po( [&]( auto f ){
    compute_sensing_tfi( f );
  } );
}

void compute_sensing_times()
{
  sensings.reset();
  std::vector<double> input_sensings( pntk->num_pis(), -1 );
  compute_sensing_times( input_sensings );
}

void compute_sensing_tfi( signal const& f )
{
  node n = pntk->get_node( f );
  if( sensings.has( n ) )
  {
    return;
  }
  auto g = pntk->get_binding( n );
  sensings[n] = std::numeric_limits<double>::max();;
  pntk->foreach_fanin( n, [&]( auto fi, auto ii ){
    node ni = pntk->get_node( fi );
    compute_sensing_tfi( fi );
    double const rise_time = g.pins[ii].rise_block_delay;
    double const fall_time = g.pins[ii].fall_block_delay;
    double const min_delay = std::min( rise_time, fall_time );
    sensings[n] = std::min( sensings[n], sensings[ni] + min_delay );    
  } );
}

void compute_sensing_tfo( node const& n )
{
  if( pntk->visited( n ) == pntk->trav_id() )
    return;
  pntk->set_visited( n, pntk->trav_id() );

  auto g = pntk->get_binding( n );
  sensings[n] = std::numeric_limits<double>::max();;
  pntk->foreach_fanin( n, [&]( auto fi, auto ii ){
    node ni = pntk->get_node( fi );
    double const rise_time = g.pins[ii].rise_block_delay;
    double const fall_time = g.pins[ii].fall_block_delay;
    double const min_delay = std::min( rise_time, fall_time );
    sensings[n] = std::min( sensings[n], sensings[ni] + min_delay );    
  } );

  pntk->foreach_fanout( n, [&]( auto no ){
    compute_sensing_tfo( no );
  } );

}

double get_time( node n )
{
  return sensings[n];
}

double get_time( signal f )
{
  return sensings[f.index];
}

#pragma endregion

private:
  std::shared_ptr<Ntk> pntk;
  unordered_node_map<double, Ntk> sensings;
  /* events */
  std::shared_ptr<typename network_events<Ntk>::add_event_type> add_event;
  std::shared_ptr<typename network_events<Ntk>::modified_event_type> modified_event;
};

} // namespace mockturtle