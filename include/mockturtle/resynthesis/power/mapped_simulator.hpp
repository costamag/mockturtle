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
  \file mapped_simulator.hpp
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

#include "../../algorithms/simulation.hpp"
#include "../list_simulator.hpp"
#include "../timing/sensing_analyzer.hpp"
#include "../timing/arrival_analyzer.hpp"
#include "node_simulator.hpp"

namespace mockturtle
{

template<class Ntk, uint32_t N = 16, uint32_t I = 10>
class mapped_simulator
{
public:
  using truth_table_t = kitty::static_truth_table<I>;
  using signal = typename Ntk::signal;
  using node = typename Ntk::node;
  using node_simulator_t = node_simulator<Ntk, N, I>;
  using list_simulator_t = aig_list_simulator<truth_table_t>;
  static constexpr node null_node = std::numeric_limits<uint32_t>::max();

  mapped_simulator( Ntk& ntk )
      : pntk( std::make_shared<Ntk>( ntk ) ),
        sims(),
        dyn_switch( 0 ),
        zero_delay( 0 ),
        list_simulator()
  {
    /* allocate enough memory for storing the simulation information */
    sims.reserve( 20 );
    truth_table_t const const0;
    /* initialize the constants */
    std::array<truth_table_t, 2> const pair0 = {  const0,  const0 };
    std::array<truth_table_t, 2> const pair1 = { ~const0, ~const0 };
    double const start = 0;
    double const delay = std::numeric_limits<double>::max() - 1;
    node_simulator_t sim0( null_node, pair0, start, start, delay );
    node_simulator_t sim1( null_node, pair1, start, start, delay );
  }

  double operator()( mapped_list const& mlist, std::vector<node_simulator_t*> const& input_sims )
  {
    sims.reserve( mlist.num_gates() + mlist.num_pis() + 2 );
    sims.resize( 2 );
    dyn_switch = 0;
    zero_delay = 0;

    int i = 2;    
    for( auto const psim : input_sims )
    {
      node_simulator_t const sim_new( *psim );
      sims.emplace_back( sim_new );
    }
    mlist.foreach_gate( [&]( auto children, auto id ){
      auto const g = pntk->get_binding_from_id( id );
      /* compute sensing time */
      double sensing = std::numeric_limits<double>::max();
      double arrival = std::numeric_limits<double>::min();
      for( int i = 0; i < children.size(); ++i )
      {
        double const rise = g.pins[i].rise_block_delay;
        double const fall = g.pins[i].fall_block_delay;
        sensing = std::min( sensing, sims[children[i]].get_sensing() + std::min( rise, fall ) );
        arrival = std::max( arrival, sims[children[i]].get_arrival() + std::max( rise, fall ) );
      }
      /* compute arrival time */
      node_simulator_t sim( null_node, sensing, arrival );

      std::vector<node_simulator_t*> sims_ptrs;
      for( int i = 0; i < children.size(); ++i )
      {
        auto const index = children[i];
        sims_ptrs.push_back( &sims[index] );
      }
      sim( g, list_simulator, sims_ptrs );

      dyn_switch += sim.get_dyn_switch();
      zero_delay += sim.get_zero_delay();
    } );
    return dyn_switch;
  }

  ~mapped_simulator()
  {
  }

private:
  std::share_ptr<Ntk> pntk;
  std::vector<node_simulator_t> sims;
  double dyn_switch;
  double zero_delay;
  list_simulator_t list_simulator;

};

} // namespace mockturtle