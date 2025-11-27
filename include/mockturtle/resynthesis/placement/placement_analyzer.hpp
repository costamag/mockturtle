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
#include "../database/mapped_list.hpp"
#include "placement_storage.hpp"
#include "placement_utils.hpp"
#include "particle.hpp"
#include <math.h>
#include <set>
#include <queue>
#include <array>

namespace mockturtle
{

template<class Ntk, uint32_t Nx, uint32_t Ny>
struct placers_t
{
  placers_t( Ntk & ntk, placement_params & ps )
  {
    for ( int i = 0; i < ps.num_placers; ++i )
      data.emplace_back( placement_storage_t<Ntk, Nx, Ny>( ntk, ps, (1+i) * 5 ) );
  }
  std::vector<placement_storage_t<Ntk, Nx, Ny>> data;
};

template<class Ntk, int Nx, int Ny>
class placement_analyzer
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

  placement_analyzer( Ntk & ntk, placement_params & ops )
    : ntk( ntk ),
      placers( ntk, ops ),
      gen( 336u )
  {
    double Uinit = compute_cost() ;
    _init();
    double Ufire = compute_cost() ;
    add_event = ntk.events().register_add_event( [&]( const auto& n ) {
    } );

    modified_event = ntk.events().register_modified_event( [&]( const auto& n, auto old_children ) {
    } );

    delete_event = ntk.events().register_delete_event( [&]( const auto& n ) {
    } );
  }

  ~placement_analyzer()
  {
    if ( add_event )
    {
      ntk.events().release_add_event( add_event );
    }

    if ( modified_event )
    {
      ntk.events().release_modified_event( modified_event );
    }

    if ( delete_event )
    {
      ntk.events().release_delete_event( delete_event );
    }
  }

  /*! \brief Random initializer setting the positions at random places */
  void _init()
  {
  }

  void fire( placement_storage_t<Ntk, Nx, Ny> & placer )
  {
    double gamma0 = 0.01;
    double gamma = gamma0;
    double gamma_dec = 0.99;
    double dt_dec = 0.5;
    double dt_inc = 1.1;
    double dt = 0.01;
    double dt_max = 0.1;
    int Nh = 7;

    //if ( ps.anchor )
    placer.mark_anchored();
    //placer.thermal_fluctuation();
    placer.reset_moments();
//    placer.randomize_positions();


    double U_new = placer.compute_cost_old();
    double U_old = U_new + 100;
    int count = 0;
    while ( ( U_old - U_new ) > ps.eps )
    {
      U_old = U_new;
      U_old = placer.compute_cost_old();

      placer.reset_moments();
      placer.compute_forces();
      placer.integrate_moments( dt );

      if ( placer.is_descending() )
      {
        if ( count++ > Nh )
        {
          dt = std::min( dt_max, dt_inc * dt );
          gamma = gamma * gamma_dec;
        }
      }
      else
      {
        placer.reset_moments();
        gamma = gamma0;
        count = 0;
        dt = dt * dt_dec;
      }
      placer.rescale_moments( gamma );
      placer.integrate_positions( dt );
      placer.legalize_positions();
      U_new = placer.compute_cost_old();
      placer.check();

      if ( ( U_old - U_new ) <= ps.eps )
      {
        //placer.replace_extremes();
        U_new = placer.compute_cost_old();
        gamma0 = 0.01;
        gamma = gamma0;
        gamma_dec = 0.99;
        dt_dec = 0.5;
        dt_inc = 1.1;
        dt = 0.01;
        dt_max = 0.1;
        Nh = 7;
        placer.reset_moments();
      }
      //std::cout << U_old << "|->>" << U_new << std::endl;
    }
  }

  void fire()
  {
    for ( auto & placer : placers.data )
    {
      fire( placer );
      placer.compute_cumwire_times();
    }
  }

  double compute_cost()
  {
    double res = 0;
    for ( auto & placer : placers.data )
    {
      res += placer.compute_cost_old();
    }
    return res / ( (double)placers.data.size() );
  }

  void print_costs()
  {
    int i = 0;
    for ( auto & placer : placers.data )
    {
      std::cout << "placer #" << i++ << std::endl; 
      placer.print_cost();
    }
  }

  double compute_hpwl()
  {
    double res = 0;
    for ( auto & placer : placers.data )
    {
      res += placer.compute_hpwl();
    }
    return res / ( (double)placers.data.size() );
  }


  void print_params_placer()
  {
    placers.data[0].print_params();
  }

  void print_state()
  {
    ntk.foreach_node( [&]( node n ){
      if ( ntk.is_dead( n ) )
      {
        std::cout << n << " is dead" << std::endl;
        return true;
      }
      vector_t<double> q = placers.data[0].get_position( n );
      std::cout << n << " " << q[0] << " " << q[1] << std::endl;
      return true;
    } );

    //placers.data[i].save_density_matrix( file_name + "_matrix.txt" );
  }

  void save_result( std::string file_name, int i = 0 )
  {
    std::cout << file_name << std::endl;
    placers.data[i].compute_cost();
    std::ofstream p( file_name + ".txt" );
    p <<  " " << placers.data[i].get_num_bins_x() <<
          " " << placers.data[i].get_num_bins_y() <<
          " " << placers.data[i].get_slots_per_bin() << 
          " " << placers.data[i].get_unit_length() << std::endl;
    ntk.foreach_pi( [&]( node n ){
      p << n << " ";
    } );
    p << std::endl;
    ntk.foreach_po( [&]( signal f ){
      node n = ntk.get_node( f );
      p << n << " ";
    } );
    p << std::endl;

    ntk.foreach_node( [&]( node n ){
      if ( ntk.is_dead( n ) )
      {
        std::cout << n << " is dead" << std::endl;
        return true;
      }
      vector_t<double> q = placers.data[i].get_position( n );
      p << n << " " << ntk.get_range( n ) << " " << q[0] << " " << q[1] << " ";
      ntk.foreach_fanin( n, [&]( auto f, auto i ){
        node ni = ntk.get_node( f );
        p << ni << " ";
      } );
      p << std::endl;
      return true;
    } );
    p.close();

    //placers.data[i].save_density_matrix( file_name + "_matrix.txt" );
  }

  template<class Leaves_t, class LeavesNew_t>
  double evaluate( node const& x, Leaves_t const& leaves, std::set<node> const& erased, LeavesNew_t & new_leaves )
  {
  //  double avg = 0;
  //  for ( auto & placer : placers.data )
  //  {
  //    avg += placer.get_cum_wires( x );
  //  }
  //  return avg / ( ( double ) placers.data.size() ) ;

    auto lib = ntk.get_bindings();
    mapped_list list( ntk, lib, leaves, ntk.make_signal( x ) );
  //  list.print();
  bool one = false;
  for ( auto l : new_leaves )
    one |= ( l < ntk.size() );
    if ( !one )
      return std::numeric_limits<double>::min();
    auto cost = evaluate_l( list, leaves, x, erased, new_leaves );
    return cost;
  }

  template<class Leaves_t, class LeavesNew_t>
  double evaluate_l( mapped_list const& list, Leaves_t const& leaves, node const& n, std::set<node> const& erased, LeavesNew_t & new_leaves )
  {
   // std::cout << "START" << std::endl;
   // std::cout << "state: " << std::endl;
   // print_state();
   /// list.print();
  //  for ( auto l : new_leaves )
  //    std::cout << l << " ";
  //  std::cout << std::endl;
  //  std::cout << "Evaluate" << std::endl;
    double res_avg = 0;

    //auto & placer = placers.data[0];
    for ( auto & placer : placers.data )
    {
      std::vector<particle_t> parts;
      placer.save_state();
      placer.init_place_list_l( parts, list, leaves, n, erased, new_leaves );

      double num_times = 1;
      for ( int t = 0; t < num_times; ++t )
      {
        double gamma0 = 0.01;
        double gamma = gamma0;
        double gamma_dec = 0.99;
        double dt_dec = 0.5;
        double dt_inc = 1.1;
        double dt = 0.001;
        double dt_max = 0.1;
        int Nh = 7;

        placer.reset_moments( parts );

        double U_new = placer.compute_costs( parts, list, leaves, n ); //H
        placer.compute_forces( parts, list, leaves, n );

        // grads are up to date
        double U_old = U_new + 100;
        int count = 0;
        while ( ( U_old - U_new ) > ps.eps )
        {
          U_old = U_new;
          placer.integrate_moments( parts, dt );

          if ( placer.is_descending( parts ) )
          {
            if ( count++ > Nh )
            {
              dt = std::min( dt_max, dt_inc * dt );
              gamma = gamma * gamma_dec;
            }
          }
          else
          {
            placer.reset_moments( parts );
            gamma = gamma0;
            count = 0;
            dt = dt * dt_dec;
          }
          placer.rescale_moments( parts, gamma );
          placer.integrate_positions( parts, dt );
          placer.legalize_positions( parts );
          U_new = placer.compute_costs( parts, list, leaves, n );
          placer.compute_forces( parts, list, leaves, n );
        }
        double placer_cost = placer.compute_costs( parts, list, leaves, n )/num_times;
        res_avg += placer_cost;
      }

      placer.reset_place_list( parts, list, leaves, n, erased, new_leaves );
      placer.restore_state();
    }
    auto res_cost = res_avg/( (double)placers.data.size() );
 //   std::cout << "final cost " << res_cost << std::endl;

    return res_cost;
  }

  template<class Leaves_t, class LeavesNew_t>
  double evaluate( mapped_list const& list, Leaves_t const& leaves, node const& n, std::set<node> const& erased, LeavesNew_t & new_leaves )
  {
  //  std::cout << "START" << std::endl;
  //  std::cout << "state: " << std::endl;
  //  print_state();
  //  list.print();
  //  for ( auto l : new_leaves )
  //    std::cout << l << " ";
  //  std::cout << std::endl;
  //  std::cout << "Evaluate" << std::endl;
    double res_avg = 0;

    ///auto & placer = placers.data[0];
    
    for ( auto & placer : placers.data )
    {
      std::vector<particle_t> parts;
      placer.save_state();
      placer.init_place_list( parts, list, leaves, n, erased, new_leaves );

      double num_times = 1;
      for ( int t = 0; t < num_times; ++t )
      {
        double gamma0 = 0.01;
        double gamma = gamma0;
        double gamma_dec = 0.99;
        double dt_dec = 0.5;
        double dt_inc = 1.1;
        double dt = 0.001;
        double dt_max = 0.1;
        int Nh = 7;

        placer.reset_moments( parts );

        double U_new = placer.compute_costs( parts, list, leaves, n ); //H
        placer.compute_forces( parts, list, leaves, n );

        // grads are up to date
        double U_old = U_new + 100;
        int count = 0;
        while ( ( U_old - U_new ) > ps.eps )
        {
          U_old = U_new;
          placer.integrate_moments( parts, dt );

          if ( placer.is_descending( parts ) )
          {
            if ( count++ > Nh )
            {
              dt = std::min( dt_max, dt_inc * dt );
              gamma = gamma * gamma_dec;
            }
          }
          else
          {
            placer.reset_moments( parts );
            gamma = gamma0;
            count = 0;
            dt = dt * dt_dec;
          }
          placer.rescale_moments( parts, gamma );
          placer.integrate_positions( parts, dt );
          placer.legalize_positions( parts );
          U_new = placer.compute_costs( parts, list, leaves, n );
          placer.compute_forces( parts, list, leaves, n );
        }
        double placer_cost = placer.compute_costs( parts, list, leaves, n )/num_times;
        res_avg += placer_cost;
      }

      placer.reset_place_list( parts, list, leaves, n, erased, new_leaves );
      placer.restore_state();
    }
    auto res_cost = res_avg/( (double)placers.data.size() );
  //  std::cout << "final cost " << res_cost << std::endl;

    return res_cost;
  }

  /* TODO define */
  double evaluate( node const& nnew, node const& nold )
  {
    double res = std::numeric_limits<double>::max();
    return res;
  }

  #pragma endregion evaluate

  double get_cost( node n )
  {
    double res = 0;
    for ( int i = 0; i < placers.data.size(); ++i )
    {
      res += placers.data[i].compute_cost( n );
    }
    return res / ( (double)placers.data.size() );
  }

  void erase_from_grid( std::set<node> & removed )
  {
    for ( int i = 0; i < placers.data.size(); ++i )
    {
      for ( auto n : removed )
      {
        placers.data[i].erase_from_grid( n );
      }
    }
  }

  void save_nets_stats( std::string path )
  {
    std::ofstream p( path );
    ntk.foreach_gate( [&]( node n ){
      auto b = ntk.get_binding( n ); 
      if ( !ntk.is_dead( n ) )
      {
        ntk.foreach_fanin( n, [&]( auto fi, auto ii ){
          node ni = ntk.get_node( fi );
          for ( int i = 0; i < ps.num_placers; ++i )
          {
            p << pow(placers.data[i].get_distance( n, ni ), 2) << " ";
          }
          p << std::endl;
        } );
      }
    } );
    p.close();
  }

  void substitution( node n, signal fnew )
  {
    for ( int i = 0; i < placers.data.size(); ++i )
    {
      placers.data[i].substitution( n, fnew );
    }
  }

private:

  placement_stats st;
  placement_params ps;

  /*! \brief Pointer to the network to be placed */
  Ntk & ntk;

  placers_t<Ntk, Nx, Ny> placers;

  std::mt19937 gen; // mersenne_twister_engine seeded with rd()

  /* events */
  std::shared_ptr<typename network_events<Ntk>::add_event_type> add_event;
  std::shared_ptr<typename network_events<Ntk>::modified_event_type> modified_event;
  std::shared_ptr<typename network_events<Ntk>::delete_event_type> delete_event;
};


} // namespace mockturtle


#if 0
#pragma region evaluate




  double copy_cost( node const& n )
  {
    double res = 0;
    for ( auto & placer : placers.data )
    {
      res += ( placer.compute_cost_place(n)+placer.compute_cost_route( n ) );
    }
    res /= ( (double)placers.data.size() );
    return res;
  }

#pragma endregion evaluate

  void print( int i = 0 )
  {   
    auto & placer = placers.data[i]; 
    for ( int x = 0; x < Nx; ++x )
    {
      printf("%4d/%4d|", x, Nx-1 );
      for ( int y = Ny - 1; y >= 0; --y )
      {
        auto nodes = placer.copy_nodes( x, y );

        for ( auto n : nodes )
        {
          auto q = placer.copy_positions( n );
          std::cout << n << ":( " << q[0] << ", " << q[1] << ") ";
        }
      }
      printf("\n");      
    }
  }
#endif