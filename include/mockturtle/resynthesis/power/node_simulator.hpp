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
  - The sensing time : the first time at which a trans can happen
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

namespace mockturtle
{

template<class Ntk, int N = 16, int I = 10>
class node_simulator
{
  static constexpr double num_pairs = static_cast<double>( 1u << I );
public:
  using truth_table_t = kitty::static_truth_table<I>;
  using list_simulator_t = aig_list_simulator<truth_table_t>;
  using index_list_t = typename list_simulator_t::index_list_t;
  using node = typename Ntk::node;

public:

  node_simulator()
    :  n( 0 ),
      sensing( 0 ),
      arrival( 0 ),
      zero_delay( 0 ),
      dyn_switch( 0 ),
      gen( 336 )
  {}

  template<class Simulator_t>
  node_simulator( Simulator_t const& other )
    :  n( other.get_node() ), sensing( other.get_sensing() ), arrival( other.get_arrival() ), zero_delay( other.get_zero_delay() ), dyn_switch( other.get_dyn_switch() ), gen( other.gen )
  {
    sims = other.get_sims();
    probs = other.get_trans_probs();
  }

  explicit node_simulator( node const& n, double const& sensing, double const& arrival ) noexcept
    :  n( n ), sensing( sensing ), arrival( arrival ), zero_delay( 0 ), dyn_switch( 0 ), gen( 336 )
  {}

  explicit node_simulator( node const& n, std::array<truth_table_t, 2> const& pattern_pairs, double const& toggle, double const& sensing, double const& arrival ) noexcept
    :  n( n ), sensing( sensing ), arrival( arrival ), zero_delay( 0 ), dyn_switch( 0 ), gen( 336 )
  {
    sims[0] = pattern_pairs[0];
    sims[N-1] = pattern_pairs[1];
    int toggle_step = get_step( toggle );
    if ( ( toggle_step % ( N - 1 ) ) == 0 )
    {
      toggle_step = ( N - 1 ) / 2;
    }

    for( int step = 1; step < ( N - 1 ); ++step )
    {
      sims[step] = step < toggle_step ? pattern_pairs[0] : pattern_pairs[1];
    }
    probs[0] = 0;
    for( int step = 1; step < N; ++step )
    {
      auto const switches = kitty::count_ones( sims[ step ] ^ sims[ step - 1 ] );
      dyn_switch += switches;
      probs[step] = switches / num_pairs;
    }

    for( int step = 0; step < N; ++step )
    {
      auto const switches = kitty::count_ones( sims[ step ] );
      prob1[step] = 0.001 + switches / num_pairs;
    }

    dyn_switch /= num_pairs;
    zero_delay = kitty::count_ones( sims[ 0 ] ^ sims[ N - 1 ] ) / num_pairs;
  }

  /*! \brief Constructor without the toggle time 
   * This constuctor assumes that the toggle time could happen at any instant in the 
   * interval between the arrival time and the sensing time.
  */
  explicit node_simulator( node const& n, std::array<truth_table_t, 2> const& pattern_pairs, double const& sensing, double const& arrival ) noexcept
    :  n( n ), sensing( sensing ), arrival( arrival ), zero_delay( 0 ), dyn_switch( 0 ), gen( 336 )
  {
    std::uniform_int_distribution<> U( 1, N - 2 );
    for ( int b = 0; b < pattern_pairs[0].num_bits(); ++b )
    {
      int toggle_step = sensing == arrival ? ( ( N - 1 ) / 2 ) : U( gen );
      for( int step = 0; step < N; ++step )
      {
        if ( step < toggle_step )
        {
          if ( kitty::get_bit( pattern_pairs[0], b ) )
            kitty::set_bit( sims[step], b );
          else
            kitty::clear_bit( sims[step], b );
        }
        else
        {
          if ( kitty::get_bit( pattern_pairs[1], b ) )
            kitty::set_bit( sims[step], b );
          else
            kitty::clear_bit( sims[step], b );
        }
      }
    }
    for( int step = 1; step < N; ++step )
    {
      auto const switches = kitty::count_ones( sims[ step ] ^ sims[ step - 1 ] );
      dyn_switch += switches;
    }

    dyn_switch /= num_pairs;
    zero_delay = kitty::count_ones( sims[ 0 ] ^ sims[ N - 1 ] ) / num_pairs;
  }

  void update_info( node const& nd, double const& sensing_time, double const& arrival_time )
  {
    n = nd;
    sensing = sensing_time;
    arrival = arrival_time;
  }

  void operator()( Ntk const& ntk, list_simulator_t & simulator, std::vector<node_simulator*> const& sims_ptrs )
  {
    auto g = ntk.get_binding( n );
    simulate( g, simulator, sims_ptrs );
  }

  template<class Binding_t>
  void operator()( Binding_t const& g, list_simulator_t & simulator, std::vector<node_simulator*> const& sims_ptrs )
  {
    simulate( g, simulator, sims_ptrs );
  }

  template<class Binding_t>
  void operator()( Binding_t const& g, list_simulator_t & simulator, std::vector<node_simulator> const& sims_vals )
  {
    std::vector<node_simulator*> sims_ptrs;
    for ( int i = 0; i < sims_vals.size(); ++i )
    {
      sims_ptrs.push_back( &sims_vals[i] );
    }
    simulate( g, simulator, sims_ptrs );
  }

  template<class Binding_t>
  void simulate( Binding_t const& g, list_simulator_t & simulator, std::vector<node_simulator*> const& sims_ptrs )
  {
    index_list_t const& list = g.list;
    zero_delay = 0;
    dyn_switch = 0;
    std::vector<truth_table_t*> tts_ptrs( sims_ptrs.size() );

    /* initialize the simulation using the initial pattern */
    for( int i = 0; i < sims_ptrs.size(); ++i )
    {
      tts_ptrs[i] = sims_ptrs[i]->get_tt_ptr( 0 );
    }

    int toggle = ( arrival == sensing ) ? ( N - 1 ) / 2 : N - 1;
    for( int step = 0; step < toggle; ++step )
    {
      simulator( sims[step], list, tts_ptrs );
    }
    /* complete the simulation using the last pattern */
    for( int i = 0; i < sims_ptrs.size(); ++i )
    {
      tts_ptrs[i] = sims_ptrs[i]->get_tt_ptr( N - 1 );
    }
    for( int step = toggle; step < N; ++step )
    {
      simulator( sims[step], list, tts_ptrs );
    }

    /* simulate all the other timesteps */
    if ( toggle == ( N - 1 ) )
    {
      for( int step = 1; step < ( N - 1 ); ++step )
      {
        /* compute the time step */
        double t = get_time( step );
        /* iterate over all the inputs and store the simulation patterns */
        for( int i = 0; i < sims_ptrs.size(); ++i )
        {
          double delta_i = g.pin_delays[ i ];
          double time_bef = t - delta_i;
          int step_bef = sims_ptrs[i]->get_step( time_bef );
          tts_ptrs[i] = sims_ptrs[i]->get_tt_ptr( step_bef );
        }
        simulator( sims[step], list, tts_ptrs );
      }
    }

    probs[0] = 0;
    for( int step = 1; step < N; ++step )
    {
      auto const switches = kitty::count_ones( sims[ step ] ^ sims[ step - 1 ] );
      dyn_switch += switches;
      probs[step] = switches / num_pairs;
    }

    for( int step = 0; step < N; ++step )
    {
      auto const switches = kitty::count_ones( sims[ step ] );
      prob1[step] = 0.001 + switches / num_pairs;
    }

    dyn_switch /= num_pairs;
    zero_delay = kitty::count_ones( sims[ 0 ] ^ sims[ N - 1 ] ) / num_pairs;
  }

  /*! \brief get the simulation time given the simulation step in the activity window */
  double get_time( int step ) const
  {
    step = std::max( 0, std::min( N - 1, step ) );
    return sensing + step * ( arrival - sensing ) / ( N - 1 );
  }

  double get_num_pairs() const
  {
    return num_pairs;
  }

  /*! \brief get the simulation timestep of a node within the activity window */
  int get_step( double time ) const
  {
    auto const step_float = ( N - 1 ) * ( time - sensing ) / ( arrival - sensing );
    auto const step_round = std::ceil( step_float );
    int step = static_cast<int>( step_round );
    step = std::max( 0, std::min( N - 1, step ) );
    return step;
  }

  truth_table_t * get_tt_ptr( int step )
  {
    return &sims[step];
  }

  truth_table_t get_tt( int step ) const
  {
    return sims[step];
  }

  double get_prob( int step )
  {
    return probs[step];
  }

  double get_zero_delay() const
  {
    return zero_delay;
  }

  double get_dyn_switch() const
  {
    return dyn_switch;
  }

  void set_arrival( double arrival_val )
  {
    arrival = arrival_val;
  }

  void set_sensing( double sensing_val )
  {
    sensing = sensing_val;
  }

  double get_arrival() const
  {
    return arrival;
  }

  double get_sensing() const
  {
    return sensing;
  }

  double get_node() const
  {
    return n;
  }

  std::array<double, N> get_trans_probs() const
  {
    return probs;
  }

  std::array<double, N> get_probs() const
  {
    return prob1;
  }

  std::array<truth_table_t, N> get_sims() const
  {
    return sims;
  }

  std::array<double, N> get_distrs() const
  {
    return distr;
  }

  void print() const
  {
    int const num_patterns = 1u << I;
    for( int m = 0; m < num_patterns; ++m )
    {
      printf( "%4.2f ", sensing );
      for( int step = 0; step < N; ++step )
      {
        if( kitty::get_bit( sims[step], m ) )
          std::cout << "-";
        else
          std::cout << "_";
      }
      printf( " %4.2f\n", arrival );
    }
  }

  uint32_t print_intervals() const
  {
    uint32_t count = 0;
//    for( int step = 0; step < N; ++step )
//    {
//      bool const cec = ( step );
//      if( cec )
//      {
//        count += 1;
//        std::cout << "-";
//      }
//      else
//        std::cout << "_";
//    }
    return count;
  }

  void print_distribution() const
  {
    double sum=0;
    for( int step = 0; step < N; ++step )
    {
      printf("%5.2f ", distr[step] );
      sum += distr[step]*prob1[step];
    }
    // std::cout << " -> md " << sum << " dy " << dyn_switch << " zd " << zero_delay << std::endl;
    std::cout << std::endl;
  }

private:
  node n;
  std::array<truth_table_t, N> sims;
  std::array<double, N> probs;
  std::array<double, N> prob1;
  std::array<double, N> distr;
  double sensing;
  double arrival;
  std::mt19937 gen; // mersenne_twister_engine seeded with rd()
  double zero_delay;
  double dyn_switch;
};

} // namespace mockturtle
