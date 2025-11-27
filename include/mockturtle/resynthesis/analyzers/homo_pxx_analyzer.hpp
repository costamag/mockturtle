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
  \file homo_pxx_analyzer.hpp
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
#include "../power/node_simulator.hpp"
#include "../../algorithms/detail/resub_utils.hpp"
#include "../../algorithms/reconv_cut.hpp"
#include "../algorithms/resynthesize.hpp"

#include <fmt/format.h>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>

namespace mockturtle
{

struct pxx_resyn_params
{
  bool only_critical = true;
};

template<class Ntk, uint32_t num_vars, uint32_t num_steps, uint32_t cube_size, uint32_t max_leaves, uint32_t extra_vars = 0>
class homo_pxx_analyzer
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;
  using cost = double;
  static constexpr bool pass_manager = false;
  static constexpr cost zero_cost = 0;
  static constexpr cost max_cost = std::numeric_limits<cost>::max();
  static constexpr cost min_cost = std::numeric_limits<cost>::min();

  using truth_table_t = kitty::static_truth_table<cube_size>;
  using node_simulator_t = node_simulator<Ntk, num_steps, cube_size>;
  using list_simulator_t = aig_list_simulator<truth_table_t>;

public:
  homo_pxx_analyzer( Ntk& ntk, std::vector<double> pattern = {} )
    : pntk( std::make_shared<Ntk>( ntk ) ),
      lib( ntk.get_bindings() ),
      nodes( std::make_shared<std::vector<typename Ntk::node>>() ),
      arrival( ntk )
  {
    init( pattern );
  }

  void init( std::vector<double> pattern )
  {
    std::vector<double> input_arrival( pntk->num_pis(), +0.00 );
    if( pattern.size() != 0 )
      input_arrival = pattern;
    arrival.compute_arrival_times( input_arrival );
    double const delay = arrival.worst_delay();
    nodes->reserve( pntk->size() );
  }

  void operator()( node const& n )
  {
  }


  cost evaluate( mapped_list const& list, std::vector<node> const& leaves, node const& n )
  {
    (void)n;
    std::vector<double> times( list.size(), 0 );
    times[0] = 0;
    times[1] = 0;
    list.foreach_pi( [&]( auto v ) {
      node n = leaves[ v - 2 ];
      if( pntk->is_null( n ) )
        times[v] = 1 << 20;
      else
        times[v] = arrival.get_time( n );
    } );
    
    uint32_t index = 2 + list.num_pis();
    
    double node_arrival;
    list.foreach_gate( [&]( auto children, auto id ){
      node_arrival = std::numeric_limits<double>::min();
      for( int i = 0; i < children.size(); ++i )
      {
        auto const g = lib[id];
        uint32_t const child = children[i];
        double const rise = g.pins[i].rise_block_delay;
        double const fall = g.pins[i].fall_block_delay;
        double const child_arrival = times[child];
        node_arrival = std::max( node_arrival, child_arrival + std::max( rise, fall ) );
      }
      times[index++] = node_arrival;

    } );
    return times[list.po_at(0)];
  }

  cost evaluate( mapped_list const& list, std::array<double, num_vars> const& itimes )
  {
    std::vector<double> times( list.size(), 0 );
    times[0] = 0;
    times[1] = 0;
    list.foreach_pi( [&]( auto v ) {
      times[v] = itimes[v-2];
    } );
    
    uint32_t index = 2 + list.num_pis();
    
    double node_arrival;
    list.foreach_gate( [&]( auto children, auto id ){
      node_arrival = std::numeric_limits<double>::min();
      for( int i = 0; i < children.size(); ++i )
      {
        auto const g = lib[id];
        uint32_t const child = children[i];
        double const rise = g.pins[i].rise_block_delay;
        double const fall = g.pins[i].fall_block_delay;
        double const child_arrival = times[child];
        node_arrival = std::max( node_arrival, child_arrival + std::max( rise, fall ) );
      }
      times[index++] = node_arrival;

    } );
    return times[list.po_at(0)];
  }

  void update()
  {}

  template<class Leaves_t>
  cost evaluate( mapped_list const& list, Leaves_t const& leaves, node const& n )
  {
    (void)n;
    std::vector<double> times( list.size(), 0 );
    times[0] = 0;
    times[1] = 0;
    list.foreach_pi( [&]( auto v ) {
      node n = leaves[ v - 2 ];

      if( pntk->is_null( n ) )
        times[v] = 1 << 20;
      else
        times[v] = arrival.get_time( n );
    } );
    uint32_t index = 2 + list.num_pis();
    
    double node_arrival;
    list.foreach_gate( [&]( auto children, auto id ){
      node_arrival = std::numeric_limits<double>::min();
      for( int i = 0; i < children.size(); ++i )
      {
        auto const g = lib[id];
        uint32_t const child = children[i];
        double const rise = g.pins[i].rise_block_delay;
        double const fall = g.pins[i].fall_block_delay;
        double const child_arrival = times[child];
        node_arrival = std::max( node_arrival, child_arrival + std::max( rise, fall ) );
      }
      times[index++] = node_arrival;

    } );
    return times[list.po_at(0)];
  }

  template<class Leaves_t>
  cost evaluate( node const& n, Leaves_t const& leaves )
  {
    (void)leaves;
    return arrival.get_time( n );
  }

  cost evaluate_rewiring( node const& n, node const& ncurr )
  {
    (void)ncurr;
    return arrival.get_time( n );
  }

  bool is_interesting( node n )
  {
    if( !arrival.has( n ) )
      return false;
    return get_required( n ) <= (get_arrival( n )+0.1);
  }

  double compare( cost const& current, cost const& candidate )
  {
    return current - candidate;
  }

  double get_time( node const& n )
  {
    return arrival.get_time( n );
  }

  double get_required( node const& n )
  {
    return arrival.get_required( n );
  }

  double get_arrival( node const& n )
  {
    return arrival.get_time( n );
  }


  template<typename Fn>
  void foreach_gate( Fn&& fn ) const
  {
    for ( int i = 0; i < nodes->size(); i++ )
    {
      fn( (*nodes)[i] );
    }
  }

  void choose_gates()
  {
    if( ps.only_critical )
    {
      arrival.compute_required_times();
    }

    uint32_t const size = pntk->size();
    pntk->foreach_gate( [&]( auto n, auto i ) {
      if( i > size )
      {
        return false;
      }
      if( pntk->is_pi( n ) || pntk->is_constant( n ) )
      {
        return true;
      }

      if ( pntk->is_dead( n ) )
      {
        return true;
      }
      
      if( !ps.only_critical || ( arrival.is_critical( n ) ) )
      {
        nodes->push_back( n );
      }
      return true;
    } );
  }

  double worst_delay()
  {
    return arrival.worst_delay();
  }

private:
  std::shared_ptr<Ntk> pntk;
  bindings const lib;
  std::shared_ptr<std::vector<typename Ntk::node>> nodes;
  pxx_resyn_params ps;

  arrival_analyzer<Ntk> arrival;

};

} /* namespace mockturtle */