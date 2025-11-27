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
  \file trans_density_matrix.hpp
  \brief Stores and manipulate a transition density matrix.

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
#include "node_simulator.hpp"

namespace mockturtle
{

template<int N>
class transition_matrix
{

public:

  /*! \brief compute the transition matrix from the simulation vector */
  template<class Simulator_t>
  void operator()( Simulator_t & sim )
  {
    matrix[0][0] = 0;
    matrix[1][0] = 0;
    matrix[2][0] = 0;
    matrix[3][0] = 0;

    for( int step = 1; step < N; ++step )
    {
      auto const bef = *( sim.get_tt_ptr( step-1 ) );
      auto const now = *( sim.get_tt_ptr( step ) );
      double const num_pairs = sim.get_num_pairs();
      matrix[0][step] = kitty::count_ones( ( ~bef ) & ( ~now ) ) / num_pairs;
      matrix[1][step] = kitty::count_ones( ( ~bef ) & (  now ) ) / num_pairs;
      matrix[2][step] = kitty::count_ones( (  bef ) & ( ~now ) ) / num_pairs;
      matrix[3][step] = kitty::count_ones( (  bef ) & (  now ) ) / num_pairs;
    }
  }

  /*! \brief compute the transition matrix from the simulation vector */
//  template<class Binding_t, class Simulator_t>
//  void operator()( Binding_t const& binding, Simulator_t * node, std::vector<Simulator_t*> fanins )
//  {
//    matrix[0][0] = 0;
//    matrix[1][0] = 0;
//    matrix[2][0] = 0;
//    matrix[3][0] = 0;
//    std::vector<int> steps_bef_1( binding.num_vars );
//    std::vector<int> steps_bef_2( binding.num_vars );
//    /* simulate all the other timesteps */
//    for( int step = 1; step < N; ++step )
//    {
//      /* compute the time step */
//      double t = node->get_time( step );
//      /* iterate over all the inputs and store the simulation patterns */
//      for( int i = 0; i < binding.num_vars; ++i )
//      {
//        double const delta_i = g.pin_delays[i];
//        double const time_bef_1 = t - delta_i - 1;
//        double const time_bef_2 = t - delta_i;
//        /* get the input steps for the transition */
//        steps_bef_1[i] = fanins[i]->get_step( time_bef_1 );
//        steps_bef_2[i] = fanins[i]->get_step( time_bef_2 );
//      }
//      /* compute the transition entries */
//      auto const onset = binding.function;
//      uint32_t const num_minterms = 1u << binding.num_vars;
//      for( uint32_t m1 = 0; m1 < num_minterms; ++m1 )
//      {
//        int set1 = kitty::get_bit( onset, m1 );
//        for( uint32_t m2 = 0; m2 < num_minterms; ++m2 )
//        {
//          int set2 = kitty::get_bit( onset, m2 );
//          
//          if( ( set1 == 0 ) && ( set2 == 0 ) ) /* 0->0 */
//          {
//
//          }
//          else if( ( set1 == 0 ) && ( set2 == 1 ) ) /* 0->1 */
//          {
//
//          }
//          else if( ( set1 == 1 ) && ( set2 == 0 ) ) /* 1->0 */
//          {
//            
//          }
//          else if( ( set1 == 1 ) && ( set2 == 1 ) ) /* 1->1 */
//          {
//            
//          }
//        }
//      }
//
//    }
//  }

  double compute_switching() const
  {
    double switching = 0;
    for( int step = 1; step < N; ++step )
    {
      switching += matrix[1][step];
      switching += matrix[2][step];
    }
    return switching;
  }

  void print() const
  {
    std::vector<std::string> descr { "0->0 ", "0->1 ", "1->0 ", "1->1 " };
    for( int r = 0; r < 4; ++r )
    {
      std::cout << descr[r];
      for( int t = 0; t < N; ++t )
      {
        printf( "%3.2f ", matrix[r][t] );
      }
      std::cout << std::endl;
    }
    std::cout << "transition switching: " << compute_switching() << std::endl;
  }

private:
  std::array<std::array<double, N>, 4u> matrix;
};

} // namespace mockturtle
