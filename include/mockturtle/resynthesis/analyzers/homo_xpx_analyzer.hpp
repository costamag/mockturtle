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
  \file homo_xpx_analyzer.hpp
  \brief Inplace rewrite

  \author Andrea Costamagna
*/

#pragma once

#include "../database/mapped_list.hpp"
#include "../timing/arrival_analyzer.hpp"
#include "../algorithms/resynthesize.hpp"

#include <fmt/format.h>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>

namespace mockturtle
{

struct xpx_resyn_params
{
  bool only_critical = false;
  bool preserve_depth = true;
};

template<class Ntk, uint32_t num_vars, uint32_t num_steps, uint32_t cube_size, uint32_t max_leaves>
class homo_xpx_analyzer
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
  homo_xpx_analyzer( Ntk& ntk )
      : pntk( std::make_shared<Ntk>(ntk) ),
      lib( ntk.get_bindings() ),
      nodes( std::make_shared<std::vector<typename Ntk::node>>() ),
      arrival( ntk )
  {
    init();
  }

  void init()
  {
    std::vector<double> input_arrival( pntk->num_pis(), +0.00 );
    arrival.compute_arrival_times( input_arrival );
    double const delay = arrival.worst_delay();
    nodes->reserve( pntk->size() );
  }

  void operator()( node const& n )
  {
  }

  cost evaluate( mapped_list const& list, std::array<node, num_vars> const& leaves, node const& n )
  {
    if( ps.preserve_depth )
    {
      auto time = evaluate_arrival( list, leaves );
      if( time > arrival.get_required( n ) )
        return std::numeric_limits<double>::max() - 1000;
    }

    return list.get_power( lib );
  }

  template<class Leaves_t>
  cost evaluate( node const& n, Leaves_t const& leaves )
  {
    cost cost_deref = measure_mffc_deref( n, leaves );
    cost cost_ref = measure_mffc_ref( n, leaves );
    if( std::abs( cost_ref - cost_deref ) > 0.01 )
    {
      std::cout << cost_ref << " != " << cost_deref << std::endl;
      std::cerr << "[e] different costs " << std::endl;
      assert( cost_ref == cost_deref );
    }
    return cost_deref;
  }


  cost evaluate_arrival( mapped_list const& list, std::array<node, num_vars> const& leaves )
  {
    std::vector<double> times( list.size(), 0 );
    times[0] = 0;
    times[1] = 0;
    list.foreach_pi( [&]( auto v ) {
      node n = leaves[ v - 2 ];
      if( pntk->is_null( n ) )
        times[v] = 1 << 20;
      else
        times[v] = arrival.get_time( n );
//std::cout << "t[" << v << ", " << leaves[v] << "]=" << times[v] << std::endl;
    } );
    
    uint32_t index = 2 + list.num_pis();
    
    double node_arrival;
    list.foreach_gate( [&]( auto children, auto id ){
      node_arrival = std::numeric_limits<double>::min();
      for( int i = 0; i < children.size(); ++i )
      {
        auto const g = lib[id];
        uint32_t const child = children[i];
//std::cout << " t_c" << i << " " << times[child] << std::endl;
        double const rise = g.pins[i].rise_block_delay;
        double const fall = g.pins[i].fall_block_delay;
        double const child_arrival = times[child];
        node_arrival = std::max( node_arrival, child_arrival + std::max( rise, fall ) );
      }
      times[index++] = node_arrival;
//std::cout << "t[" << index-1 << "]=" << times[index-1] << std::endl;

    } );
//std::cout << "t[" << list.po_at(0) << "]=" << times[list.po_at(0)] << std::endl;
    return times[list.po_at(0)];
  }



  cost evaluate_rewiring( node const& n, node const& ncurr  )
  {
    if( ps.preserve_depth )
    {
      if( arrival.get_arrival( n ) > arrival.get_required( ncurr ) )
      {
        return std::numeric_limits<double>::max() - 1000;
      }
    }
    return zero_cost;
  }


  template<class Leaves_t>
  cost measure_mffc_deref( node const& n, Leaves_t const& leaves )
  {
    /* reference cut leaves */
    for ( auto leaf : leaves )
    {
      if( !pntk->is_null( leaf ) )
        pntk->incr_fanout_size( leaf );
    }

    cost mffc_cost = recursive_deref( n );

    /* dereference leaves */
    for ( auto leaf : leaves )
    {
      pntk->decr_fanout_size( leaf );
    }

    return mffc_cost;
  }

  template<class Leaves_t>
  cost measure_mffc_ref( node const& n, Leaves_t const& leaves )
  {
    /* reference cut leaves */
    for ( auto leaf : leaves )
    {
      pntk->incr_fanout_size( leaf );
    }

    cost mffc_cost = recursive_ref( n );

    /* dereference leaves */
    for ( auto leaf : leaves )
    {
      pntk->decr_fanout_size( leaf );
    }

    return mffc_cost;
  }

  cost recursive_deref( node const& n )
  {
    /* terminate? */
    if ( pntk->is_constant( n ) || pntk->is_pi( n ) )
      return zero_cost;

    /* recursively collect nodes */
    cost value{ pntk->get_power( n ) };
    pntk->foreach_fanin( n, [&]( auto const& s ) {
      if ( pntk->decr_fanout_size( pntk->get_node( s ) ) == 0 )
      {
        value += recursive_deref( pntk->get_node( s ) );
      }
    } );
    return value;
  }

  cost recursive_ref( node const& n )
  {
    /* terminate? */
    if ( pntk->is_constant( n ) || pntk->is_pi( n ) )
      return zero_cost;

    /* recursively collect nodes */
    cost value{ pntk->get_power( n ) };

    pntk->foreach_fanin( n, [&]( auto const& s ) {
      if ( pntk->incr_fanout_size( pntk->get_node( s ) ) == 0 )
      {
        value += recursive_ref( pntk->get_node( s ) );
      }
    } );
    return value;
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


  double compare( cost const& current, cost const& candidate )
  {
    return current - candidate;
  }

  double get_time( node const& n )
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

private:
  std::shared_ptr<Ntk> pntk;
  std::shared_ptr<std::vector<typename Ntk::node>> nodes;
  bindings const lib;
  arrival_analyzer<Ntk> arrival;
  xpx_resyn_params ps;

};

} /* namespace mockturtle */