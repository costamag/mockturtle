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
  \file homo_pla_analyzer.hpp
  \brief Inplace rewrite

  \author Andrea Costamagna
*/

#pragma once

#include "../database/mapped_list.hpp"

#include "../networks/bindings.hpp"
#include "../windowing/window_manager.hpp"
#include "../../algorithms/simulation.hpp"
#include "../list_simulator.hpp"
#include "../timing/sensing_analyzer.hpp"
#include "../timing/arrival_analyzer.hpp"
#include "../timing/switch_analyzer.hpp"
#include "../power/node_simulator.hpp"
#include "../placement/placement_analyzer.hpp"
#include "../power/analyzers/convo_analyzer.hpp"
#include "../../algorithms/detail/resub_utils.hpp"
#include "../../algorithms/reconv_cut.hpp"
#include "analyzers_common.hpp"

#include <fmt/format.h>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>


namespace mockturtle
{
  
  template<class Ntk, uint32_t num_vars, uint32_t num_steps, uint32_t cube_size, uint32_t max_leaves, uint32_t Size1, uint32_t Size2, uint32_t extra_vars = 0>
  class homo_pla_analyzer
  {
  public:
    using node = typename Ntk::node;
    using signal = typename Ntk::signal;
    using cost = double;
    static constexpr bool pass_manager = true;
    static constexpr cost zero_cost = 0.01;
    static constexpr cost min_cost = std::numeric_limits<cost>::min();
    static constexpr cost max_cost = std::numeric_limits<cost>::max() - 1000;
    bool VERBOSEg = false;

  public:

    homo_pla_analyzer( Ntk& ntk, analyzer_params ps = {}  )
    : ntk( ntk ),
      lib( ntk.get_bindings() ),
      placement( ntk, ps.placement_ps ),
      arrival( ntk ),
      ps( ps )
    {
      init();
    }

    /* Performs initial placement and computes arrival times */
    void init()
    {
      double Unew = placement.compute_cost();
      double Uold = Unew + 10;
      while ( Uold > ( Unew + 1 ) )
      {
        Uold = Unew;
        placement.fire();
        Unew = placement.compute_cost();
      }
      arrival.compute_arrival_times();
    }

    void update( node n, signal fnew )
    {
      placement.substitution( n, fnew );
      double U_new = placement.compute_cost();
      double U_old = U_new + 10;
      while ( U_new < ( U_old - 1 ) )
      {
        U_old = U_new;
        placement.fire();
        U_new = placement.compute_cost();
      }
    }


    void reset()
    {
      mffc.clear();
    }

    template<class Leaves_t>
    void collect_mffc( node const& n, Leaves_t const& leaves )
    {
      mffc.clear();
      collect_mffc_deref( n, leaves );
      clean_mffc_ref( n, leaves );
    }

    void operator()( node const& n )
    {
      reset();
    }

    void operator()( node const& n, window_manager<Ntk, max_leaves> const& win_manager )
    {
      reset();
      window_nodes.clear();
      for ( auto l : win_manager.leaves )
      {
        window_nodes.insert(l);
      }
      for( auto d : win_manager.divs )
      {
        window_nodes.insert(d);
      }
      for( auto d : win_manager.mffc )
      {
        window_nodes.insert(d);
      }

    }

    /*! \brief Evaluate the cost increase due to inserting the list into the network */
    template<class Leaves_t>
    cost evaluate( mapped_list const& list, Leaves_t const& leaves, node const& n )
    {

      if( ps.preserve_depth )
      {
        auto time = evaluate_arrival( list, leaves );
        if( time > arrival.get_required( n ) )
        {
          return max_cost;
        }
      }

      collect_mffc( n, leaves );

      double area_mffc = 0;
      for( auto n : mffc )
      {
        auto b = ntk.get_binding( n );
        area_mffc += b.area;
      }
      double area_list = 0;
      list.foreach_gate( [&]( auto children, auto id ){
        area_list += lib[id].area;
      } );
      if( area_list > area_mffc )
      {
        return max_cost;
      }

//      if ( list.num_gates() > mffc.size() )
//        return max_cost;
      return placement.evaluate( list, leaves, n, mffc, leaves );
    }

    double get_arrival( node n )
    {
      return arrival.get_time( n );
    }

    template<class Cut_t>
    cost evaluate_arrival( mapped_list const& list, Cut_t const& leaves )
    {
      std::vector<double> times( list.size(), 0 );
      times[0] = 0;
      times[1] = 0;
      for ( int i = 0; i < leaves.size(); ++i )
      {
        node n = leaves[i];
        if( ntk.is_null( n ) )
          times[i+2] = 1 << 20;
        else
          times[i+2] = arrival.get_time( n );
      }
    
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

    bool has_node( node n )
    {
      return window_nodes.find( n ) != window_nodes.end();
    }

    /*! \brief Computes the cost obtained by removing the mffc nodes */
    template<class Leaves_t, class LeavesNew_t>
    cost evaluate( node const& x, Leaves_t const& leaves, LeavesNew_t const & new_leaves )
    {
      for( auto l : leaves )
      {
        if( !ntk.is_null( l ) )
        {
          if ( !has_node( l ) )
          {
            if ( VERBOSEg )
              std::cout << "!has " << l << std::endl;
            return 0;
          }
        }
      }
      collect_mffc( x, leaves );
      return placement.evaluate( x, leaves, mffc, new_leaves );
    }

    void remove_mffc_info()
    {
      placement.erase_from_grid( mffc );
    }

    template<class Times_t>
    cost evaluate( mockturtle::mapped_list& mlist, Times_t& times )
    {
      assert(0);
      return std::numeric_limits<cost>::max();
    } 

    cost evaluate_rewiring( node const& nnew, node const& nold )
    {
      if( ps.preserve_depth )
      {
        if( arrival.get_arrival( nnew ) > arrival.get_required( nold ) )
        {
          return std::numeric_limits<double>::max() - 1000;
        }
      }
      return placement.evaluate( nnew, nold );
    }

    template<class Leaves_t>
    void collect_mffc_deref( node const& n, Leaves_t const& leaves )
    {
      /* reference cut leaves */
      for ( auto leaf : leaves )
      {
        if( ntk.is_null( leaf ) )
          continue;
        ntk.incr_fanout_size( leaf );
      }

      recursive_deref( n );

      /* dereference leaves */
      for ( auto leaf : leaves )
      {
        if( ntk.is_null( leaf ) )
          continue;
        ntk.decr_fanout_size( leaf );
      }
    }

    template<class Leaves_t>
    void clean_mffc_ref( node const& n, Leaves_t const& leaves )
    {
      /* reference cut leaves */
      for ( auto leaf : leaves )
      {
        if( ntk.is_null( leaf ) )
          continue;
        ntk.incr_fanout_size( leaf );
      }

      recursive_ref( n );

      /* dereference leaves */
      for ( auto leaf : leaves )
      {
        if( ntk.is_null( leaf ) )
          continue;
        ntk.decr_fanout_size( leaf );
      }
    }

    void recursive_deref( node const& n )
    {
      /* terminate? */
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
        return;

      /* recursively collect nodes */
      node nc = n;
      mffc.insert( nc );

      ntk.foreach_fanin( n, [&]( auto const& s ) {
        if ( ntk.decr_fanout_size( ntk.get_node( s ) ) == 0 )
        {
          recursive_deref( ntk.get_node( s ) );
        }
      } );
    }

    void recursive_ref( node const& n )
    {
      /* terminate? */
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
        return;

      ntk.foreach_fanin( n, [&]( auto const& s ) {
        if ( ntk.incr_fanout_size( ntk.get_node( s ) ) == 0 )
        {
          recursive_ref( ntk.get_node( s ) );
        }
      } );
      return;
    }

    double compare( cost const& current, cost const& candidate )
    {
      return current - candidate;
    }

    double get_time( node const& n )
    {
      return arrival.get_time( n );
    }

    bool is_interesting( node const& n )
    {
      (void)n;
      return true;
    }

    double get_cost()
    {
      return placement.compute_cost();
    }

  /* Using the current placemet information */
  void choose_gates()
  {
    placement.compute_cost();

    uint32_t const size = ntk.size();
    nodes.clear();
    nodes.reserve( size );

    std::vector<std::pair<node, double>> nodes_with_hotness;

    ntk.foreach_gate( [&]( auto n, auto i ) {
      if( i > size )
      {
        return false;
      }
      if( ntk.is_pi( n ) || ntk.is_constant( n ) )
      {
        return true;
      }

      if ( ntk.is_dead( n ) )
      {
        return true;
      }
      
      nodes_with_hotness.push_back( std::make_pair( n, placement.get_cost( n ) ) );
      return true;
    } );
    std::sort( nodes_with_hotness.begin(), nodes_with_hotness.end(), []( auto a, auto b ){
      return a.second > b.second;
    } );

    for( auto n : nodes_with_hotness )
      nodes.push_back( n.first );
  }

  template<typename Fn>
  void foreach_gate( Fn&& fn ) const
  {
    for ( int i = 0; i < nodes.size(); i++ )
    {
      fn( nodes[i] );
    }
  }

  double compute_cost()
  {
    return placement.compute_cost();
  }

  private:
    Ntk & ntk;
    bindings lib;
    std::set<node> mffc;
    placement_analyzer<Ntk, Size1, Size2> placement;
    arrival_analyzer<Ntk> arrival;
    analyzer_params ps;
    std::vector<typename Ntk::node> nodes;
    std::set<typename Ntk::node> window_nodes;


    /* grid storing literals */
  };

} /* namespace mockturtle */