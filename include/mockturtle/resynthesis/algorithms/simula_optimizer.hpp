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
  \file simula_optimizer.hpp
  \brief Inplace rewrite

  \author Andrea Costamagna
*/

#pragma once

#include "../boolean/dependency_cut.hpp"
#include "../boolean/spfd_manager.hpp"
#include "../../algorithms/cut_enumeration.hpp"
#include "../../algorithms/cut_enumeration/rewrite_cut.hpp"

#include <fmt/format.h>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>

namespace mockturtle
{

/*! \brief Parameters for Rewrite.
 *
 * The data structure `rewrite_params` holds configurable parameters with
 * default arguments for `rewrite`.
 */
struct simula_optimizer_params
{
  simula_optimizer_params()
  {
    /* 0 < Cut limit < 16 */
    cut_enumeration_ps.cut_limit = 8;
    cut_enumeration_ps.minimize_truth_table = true;
  }

  /*! \brief Cut enumeration parameters. */
  cut_enumeration_params cut_enumeration_ps{};

  /*! \brief If true, candidates are only accepted if they do not increase logic depth. */
  bool preserve_depth{ false };

  /*! \brief Allow rewrite with multiple structures */
  bool allow_multiple_structures{ true };

  /*! \brief Allow zero-gain substitutions */
  bool allow_zero_gain{ false };

  /*! \brief Use satisfiability don't cares for optimization. */
  bool use_dont_cares{ false };

  /*! \brief Window size for don't cares calculation. */
  uint32_t window_size{ 8u };
  uint32_t fanout_limit{100};
  int decompose { 0 };

  /*! \brief Be verbose. */
  bool verbose{ false };

  uint64_t effort{ 1 };
  bool greedy{ true };
  bool time{ false };
};

/*! \brief Statistics for rewrite.
 *
 * The data structure `rewrite_stats` provides data collected by running
 * `rewrite`.
 */
struct simula_optimizer_stats
{
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{ 0 };

  /*! \brief Expected gain. */
  uint32_t estimated_gain{ 0 };

  /*! \brief Candidates */
  uint32_t candidates{ 0 };

  void report() const
  {
    std::cout << fmt::format( "[i] total time       = {:>5.2f} secs\n", to_seconds( time_total ) );
  }
};

template<class Ntk, class Database_t, class Analyzer_t, uint32_t num_vars, uint32_t cube_size, uint32_t extra_vars = 0>
class simula_optimizer
{
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;
  using cost = typename Analyzer_t::cost;
  static constexpr cost zero_cost = Analyzer_t::zero_cost;
  static constexpr uint32_t extra_space = 20;
  static constexpr uint32_t offset = 2;
  using cut_t = static_dependency_cut<Ntk, num_vars+extra_vars>;
  using network_cuts_t = std::vector<cut_t>;
  using truth_table_t = kitty::static_truth_table<num_vars+extra_vars>;
  using signature_t = kitty::static_truth_table<cube_size>;
  using word_t = kitty::static_truth_table<6u>;

public:
  simula_optimizer( Ntk & ntk, Database_t & database, Analyzer_t & analyzer, simula_optimizer_params const& ps, simula_optimizer_stats& st )
    : ntk( ntk ),
      database( database ),
      analyzer( analyzer ),
      ps( ps ),
      st( st ),
      cuts( extra_space ),
      spfd_manager(),
      best_reward( zero_cost ),
      node_to_sim( ntk ),
      node_to_word( ntk ),
      num_cex( 1 ),
      proj_functions(),
      best_leaves(),
      best_list( num_vars+extra_vars ),
      list_simulator_sign(),
      list_simulator_word()
  {
    add_event = ntk.events().register_add_event( [&]( const auto& n ) {
      signal f = ntk.make_signal( n );
      node_to_sim.resize();
      node_to_word.resize();
      simulate_rec( f, node_to_sim, list_simulator_sign );
    } );

    modified_event = ntk.events().register_modified_event( [&]( const auto& n, auto children ) {
      signal f = ntk.make_signal( n );
      node_to_sim.resize();
      node_to_word.resize();
      simulate_rec( f, node_to_sim, list_simulator_sign );
    } );
    init();
  }

  ~simula_optimizer()
  {
    ntk.events().release_add_event( add_event );
    ntk.events().release_modified_event( modified_event );
  }

  void init()
  {
    /* reserve and initialize the signatures for exhaustive window simulation */
    node_to_sim.resize();
    node_to_word.resize();
    node const n0 = ntk.get_node( ntk.get_constant( false ) );
    node const n1 = ntk.get_node( ntk.get_constant( true ) );
    signature_t sig0;
    word_t word0;
    node_to_sim[n0] = sig0;
    node_to_sim[n1] = ~sig0;
    node_to_word[n0] = word0;
    node_to_word[n1] = ~word0;
    
    std::default_random_engine::result_type seed = 42;

    ntk.foreach_pi( [&]( auto pi, auto ii ){
      kitty::create_random( node_to_sim[pi], seed++ );
      seed += ii;
      kitty::create_random( node_to_word[pi], seed++ );
      seed += ii;
    } );

    for( int i = 0; i < (num_vars+extra_vars); ++i )
    {
      kitty::create_nth_var( proj_functions[i], i );
    }

    simulate_network( node_to_sim, list_simulator_sign );
  }

  template<class Tt_t>
  void simulate_network( incomplete_node_map<Tt_t, Ntk> & sims, aig_list_simulator<Tt_t> & simulator )
  {
    ntk.incr_trav_id();
    ntk.foreach_po( [&]( auto fo, auto io ){
      simulate_rec( fo, sims, simulator );
    } );  
  }

  template<class Tt_t>
  void simulate_rec( signal f, incomplete_node_map<Tt_t, Ntk> & sims, aig_list_simulator<Tt_t> & simulator )
  {
    node n = ntk.get_node( f );

    if( sims.has( n ) )
      return;

    if( ntk.is_pi( n ) || ( ntk.visited( n ) == ntk.trav_id() ) )
      return;

    std::vector<Tt_t*> sim_ptrs;
    ntk.foreach_fanin( n, [&]( auto fi, auto ii ){
      node ni = ntk.get_node( fi );
      simulate_rec( fi, sims, simulator );

      sim_ptrs.push_back( &sims[ni] );
    } );
    
    auto const alist = ntk.get_binding( n ).list;
    simulator( sims[n], alist, sim_ptrs );
    ntk.set_visited( n, ntk.trav_id() );
  }

  template<class WinManager_t>
  void analyze( node const& n, WinManager_t const& window_manager )
  {
    std::array<node, num_vars+extra_vars> times;
    std::array<node, num_vars+extra_vars> leaves;
    std::array<node, num_vars+extra_vars> perm_leaves;
    std::array<uint8_t, num_vars+extra_vars> permutation;

    best_reward = zero_cost;
    /* terminate if the window is not valid */
    if( !window_manager.is_valid() )
    {
      return;
    }

    if( !compute_dependency_cuts( n, window_manager ) )
    {
      return;
    }

    for ( auto const& cut : cuts )
    {
      /* extract the cut leaves and functionality */
      auto const leaves = cut.get_leaves();

//for( auto l : leaves )
//  std::cout << l << " ";
//std::cout << std::endl;

      auto cost_curr = analyzer.evaluate( n, window_manager.leaves, leaves );
      auto const& func = cut.get_function();
      auto const& care = cut.get_careset();
//kitty::print_binary( func );
  if ( kitty::is_const0( func&care ) )
  {
    best_list.create_po( 0 );
    best_reward = 1;
    return;
  }
  else if ( kitty::equal( func&care, care ) )
  {
    best_list.create_po( 1 );
    best_reward = 1;
    return;
  }
//std::cout << std::endl;

      int cnt = 0;
      for( int j{0}; j < ( num_vars + extra_vars ); ++j )
      {
        if( !ntk.is_null( leaves[j] ) )
        {
          cnt++;
          auto tt_leaf = node_to_sim[leaves[j]];
          auto tt_root = node_to_sim[n];
          if( kitty::equal( tt_root, tt_leaf ) )
          {
            cost cost_cand = analyzer.evaluate_rewiring( leaves[j], n );
            cost reward = analyzer.compare( cost_curr, cost_cand );

            if( ( reward > best_reward ) )
            {
              mapped_list list( num_vars + extra_vars );
              list.create_po( list.pi_at( j ) );
              best_list = list;
              best_reward = reward;
              best_leaves = leaves;
            }
          }
        }
      }

      std::vector<uint64_t> dont_cares;
      for( int m = 0; m < func.num_bits(); ++m )
      {
        if( kitty::get_bit( care, m ) == 0 )
          dont_cares.push_back( m );
      }
      uint64_t num_funcs = 1u << dont_cares.size();
      num_funcs = std::min( num_funcs, ps.effort );

      for( uint64_t m = 0; m < num_funcs; ++m )
      {
       // std::cout << m << "/" << num_funcs << std::endl;
        truth_table_t onset = func & care;
        truth_table_t tmp = onset | ~onset;

        for( int dc = 0; dc < dont_cares.size(); ++dc )
        {
          if( (( m >> dc ) & 0x1) == 0x1 )
            kitty::set_bit( onset, dont_cares[dc] );
          else
            kitty::clear_bit( onset, dont_cares[dc] );
        }
// kitty::print_binary( onset ); std::cout << std::endl;
        auto ttnew = onset;
        int cnt = 0;

        if constexpr ( extra_vars == 0 )
        {
          /* extract the database entries */
          auto res = database.get_entry( onset );
          if( res )
          {
            auto perm = std::get<0>( *res );
            auto entry = std::get<1>( *res );

            for ( auto j = 0u; j < (num_vars+extra_vars); ++j )
            {
              permutation[perm[j]] = j;
            }

            for( int i = 0; i < (num_vars+extra_vars); ++i )
            {
              auto const l = leaves[i];
              perm_leaves[permutation[i]] = l;
              times[permutation[i]] = !ntk.is_null( l ) ? analyzer.get_time( l ) : std::numeric_limits<double>::max();
            }

            database.foreach_list( entry, [&]( auto list, auto index ){
              auto const symm = database.get_symm( entry );
              list.match_times( times, symm );

              cost cost_cand = analyzer.evaluate( list, perm_leaves, n );
              cost reward = analyzer.compare( cost_curr, cost_cand );
              if( ( reward > best_reward ) )
              {
                best_reward = reward;
                best_list = list;
                best_leaves = perm_leaves;
              }
            } );
            //std::cout << std::endl;
          }
        }
        else
        {
         /* In this method there are three levels of abstraction:
             - LUT LIST
             - MAPPED LIST
             - SMALL LISTS
           */
          /* contains the times at the  */
          std::vector<double> times_lut;
          std::vector<double> times_lut_det{ 0, 0 };
          std::vector<uint32_t> indices_lut;
          std::vector<uint32_t> indices_list;
          std::vector<uint32_t> lut_to_list;
          for ( int k = 0; k < num_vars + extra_vars + 2; ++k )
          {
            indices_list.push_back( k );
            lut_to_list.push_back( k );
          }

          int j = 0;
          for ( int i = 0; i < num_vars + extra_vars; ++i )
          {
            auto n = leaves[i];
            if ( !ntk.is_null( n ) )
            {
              double T = analyzer.get_time( n );
              times_lut.push_back( T );
              times_lut_det.push_back( T );
//std::cout << j+2 << " times_lut_det : node " << n << " : " << times_lut_det[j+2] << std::endl;
              indices_lut.push_back( j++ );
            }
          }
          while ( times_lut_det.size() < ( num_vars + extra_vars + 2 )  )
            times_lut_det.push_back( std::numeric_limits<double>::max() );

          mapped_list cand_list( num_vars + extra_vars );
          truth_table_list_t<num_vars, num_vars+extra_vars> tt_list;
// kitty::print_binary( func ); std::cout << std::endl;
          uint32_t lit = std::numeric_limits<uint32_t>::max();
          if ( ps.decompose == 0 )
            lit = tt_list.delay_oriented_shannon( func, times_lut, indices_lut );
          else if ( ps.decompose == 1 )
            lit = tt_list.delay_oriented_shannon3( func, times_lut, indices_lut );
            if ( tt_list.is_null( lit ) )
              continue;

            if ( lit == 1 )
            {
              mapped_list lout( num_vars + extra_vars );
              lout.create_po( 1 );
              best_reward = 10;
              return;
            }
            if ( lit == 0 )
            {
              mapped_list lout( num_vars + extra_vars );
              lout.create_po( 0 );
              best_reward = 10;
              return;
            }

          uint32_t out_lit;
          tt_list.foreach_gate( [&]( auto nd ){
            auto const& tt = nd.func;
            auto const& lit = nd.lit;
            auto const& children = nd.children;
            /* extract the database entry */
            auto res = database.get_entry( tt );
            if ( res )
            {
              auto perm = std::get<0>( *res );
              auto entry = std::get<1>( *res );

              std::array<uint32_t, num_vars> permutation{ std::numeric_limits<uint32_t>::max() };
              for ( int i = 0; i < perm.size(); ++i )
              {
                permutation[perm[i]] = i;
              }
              std::array<uint32_t, num_vars> perm_leaves;
              std::array<double, num_vars> times_cut;
              for ( int i = 0; i < num_vars; ++i )
              {
                /* literal in the truth table representation */
//  std::cout << "children[i] " << children[i] << std::endl;
                auto lit = tt_list.is_null( children[i] ) ? children[i] : lut_to_list[children[i]];
//  if ( !tt_list.is_null( children[i] ) )
//    std::cout << "lut to list " << lut_to_list[children[i]] << std::endl;
//  std::cout << "lit " << lit << std::endl;
                perm_leaves[permutation[i]] = lit;
                times_cut[permutation[i]] = tt_list.is_null( children[i] ) ? std::numeric_limits<double>::max() : times_lut_det[children[i]]; // remove the offset of the constant 0 and 1
//  if ( !tt_list.is_null( children[i] ) )
//    std::cout << "times_cut[permutation[i]] " << times_cut[permutation[i]] << std::endl;
              }

              double min_cost = std::numeric_limits<double>::max();
              mapped_list best_list_s( num_vars );
              database.foreach_list( entry, [&]( auto list, auto index ){
                auto const symm = database.get_symm( entry );
                list.match_times( times_cut, symm );
                cost cost_cand = analyzer.evaluate( list, times_cut ); // TODO define this function
                if ( cost_cand < min_cost )
                {
                  min_cost = cost_cand;
                  best_list_s = list;
                }
              } );

              times_lut_det.push_back( min_cost );
              /* map the small mapped list to the larger one */
              std::vector<uint32_t> small_to_large { 0, 1 };
              for ( int i = 0; i < num_vars; ++i )
                small_to_large.push_back( perm_leaves[i] );

              best_list_s.foreach_gate( [&]( auto fanin, auto id ){
                std::vector<uint32_t> new_children;
                for ( int i = 0; i < fanin.size(); ++i )
                  new_children.push_back( small_to_large[fanin[i]] );
                out_lit = cand_list.create_node( new_children, id );
                small_to_large.push_back( out_lit );
                indices_list.push_back( out_lit );
              } );
              if ( best_list_s.num_gates() == 0 )
              {
                out_lit = small_to_large[ best_list_s.po_at(0) ];
                indices_list.push_back( out_lit );
              }

//  std::cout << " out lit " << out_lit << std::endl;
              lut_to_list.push_back( out_lit );
            }
            else
            {
//  std::cout << "no match found" << std::endl; 
              best_reward = zero_cost;
              return;
            }
          } );
          cand_list.create_po( out_lit );
          best_list = cand_list;
          std::array<node, num_vars + extra_vars> large_leaves;
          best_leaves = leaves;
          cost cost_cand = analyzer.evaluate( best_list, best_leaves, n );
          cost reward = analyzer.compare( cost_curr, cost_cand );
          best_reward = reward;
        }
      }
    }
  }

  template<class WinManager_t>
  void print_simulation( WinManager_t const& window_manager )
  {
    std::cout << "leaves" << std::endl;
    window_manager.foreach_leaf( [&]( auto l, auto i ){
      std::cout << l << " ";
      kitty::print_binary( node_to_sim[l] );
      std::cout << std::endl;
    } );

    std::cout << "divs" << std::endl;
    window_manager.foreach_divisor( [&]( auto d, auto i ){
      std::cout << d << " ";
      kitty::print_binary( node_to_sim[d] );
      std::cout << std::endl;
    } );

    std::cout << "mffc" << std::endl;
    window_manager.foreach_mffc( [&]( auto m, auto i ){
      std::cout << m << " ";
      kitty::print_binary( node_to_sim[m] );
      std::cout << std::endl;
    } );
  }

  template<class WinManager_t>
  bool compute_dependency_cuts( node n, WinManager_t const& window_manager )
  {
    /* prepare the data */
    cuts.clear();
    /* call greedy support selection */
    /* call greedy support selection */
    if( ps.greedy == true )
      return greedy_support_selection( n, window_manager );
    else if( ps.time == true )
      return times_support_selection( n, window_manager );
    else
      return pivot_support_selection( n, window_manager );

  }

  template<class WinManager_t>
  bool pivot_support_selection( node const& n, WinManager_t const& window_manager )
  {
    bool found = false;
    window_manager.foreach_divisor( [&]( auto piv, auto i ){
      if( found )
        return true;

      /* initialize the spfd manager */
      signature_t const & sign = node_to_sim[n];
      spfd_manager.init( sign );
      
      signature_t const& sig = node_to_sim[piv];
      auto min_cost = spfd_manager.evaluate( sig );
      spfd_manager.update( sig );
      std::array<node, num_vars+extra_vars> leaves;
      leaves[0] = piv;
      int cnt = 1;

      /* prepare the data */
      node candidate = ntk.get_null();
      while( !spfd_manager.is_covered() && !spfd_manager.is_saturated() && ( cnt < ( num_vars + extra_vars ) ))
      {
        bool found_new = false;
        window_manager.foreach_divisor( [&]( auto d, auto i ){
          if( !ps.preserve_depth || ( analyzer.get_time( d ) < (analyzer.get_time( n )) ) )
          {
            signature_t const& sig = node_to_sim[d];
            auto const cost_div = spfd_manager.evaluate( sig );
            if( cost_div < min_cost )
            {
              found_new = true;
              min_cost = cost_div;
              candidate = d;
            }
          }
        } );
        if( !found_new )
        {
          return true;
        }
        signature_t const& sig = node_to_sim[candidate];
        spfd_manager.update( sig );
        leaves[cnt++] = candidate;
      }
      if( spfd_manager.is_covered() )
      {
        while( cnt < (num_vars+extra_vars) )
          leaves[cnt++] = ntk.get_null();
  
        auto [func, care] = compute_function( n, leaves );
        cuts.emplace_back( n, leaves, func, care );
        found = true;
      }
      return true;
    } );

    return found;
  }
  template<class WinManager_t>
  bool leaves_support_selection( node const& n, WinManager_t const& window_manager )
  {
    std::array<node, num_vars + extra_vars> leaves;
    bool found = false;
    if( window_manager.leaves.size() <= ( num_vars + extra_vars )  )
    {
      for ( int i = 0; i < num_vars + extra_vars; ++i )
      {
        if ( i < window_manager.leaves.size() )
          leaves[i] = window_manager.leaves[i];
        else
          leaves[i] = ntk.get_null();
      }
      auto [func, care] = compute_function( n, leaves );
      cuts.emplace_back( n, leaves, func, care );
      found = true;
    }
    return found;
  }


  template<class WinManager_t>
  bool times_support_selection( node const& n, WinManager_t const& window_manager )
  {
    bool found = false;
    double min_time = std::numeric_limits<double>::max();
    double max_time = analyzer.get_time( n );
    std::vector<node> sorted_divs;
    window_manager.foreach_divisor( [&]( auto d, auto i ){
      sorted_divs.push_back( d );
      auto const& d_time = analyzer.get_time( d );
      min_time = std::min( min_time, d_time );
      max_time = std::max( max_time, d_time );
    } );
    max_time = analyzer.get_time( n )-20;
    std::sort( sorted_divs.begin(), sorted_divs.end(), [&]( auto na, auto nb ){
      return analyzer.get_time( na ) < analyzer.get_time( nb );
    } );

    uint32_t steps = 10;
    double delta = ( max_time - min_time ) / ((double)steps);

    for ( uint32_t s = 0; s < steps; ++s )
    {
      double threshold = min_time + s * delta;

      for ( int ip = 0; ip < sorted_divs.size(); ++ip )
      {
        auto piv = sorted_divs[ip];
        if( found )
          continue;
        if( ntk.fanout_size( piv ) > ps.fanout_limit || ntk.fanout_size( piv ) == 0 )
          continue;

        if ( analyzer.get_time( piv ) >= threshold )
          break;

        /* initialize the spfd manager */
        signature_t const & sign = node_to_sim[n];

        spfd_manager.init( sign );

        signature_t const& sig = node_to_sim[piv];
        auto min_cost = spfd_manager.evaluate( sig );
        spfd_manager.update( sig );
        std::array<node, num_vars + extra_vars> leaves;
        leaves[0] = piv;
        int cnt = 1;

        /* prepare the data */
        node candidate = ntk.get_null();
        while( !spfd_manager.is_covered() && !spfd_manager.is_saturated() && ( cnt < ( num_vars+extra_vars) ) )
        {
          bool found_new = false;
          for ( int jp = ip+1; jp < sorted_divs.size(); ++jp )
          {
            auto d = sorted_divs[jp];
            if( ntk.fanout_size( d ) <= ps.fanout_limit && ntk.fanout_size( d ) > 0 )
            {
              if( !ps.preserve_depth || ( analyzer.get_time( d ) < threshold ) )
              {
                signature_t const & sign = node_to_sim[d];

                auto const cost_div = spfd_manager.evaluate( sig );
                if( cost_div < min_cost )
                {
                  found_new = true;
                  min_cost = cost_div;
                  candidate = d;
                }
              }
            }
          }
          if( !found_new )
          {
            break;
          }
          signature_t const& sig = node_to_sim[candidate];
          spfd_manager.update( sig );
          leaves[cnt++] = candidate;
        }
        if( spfd_manager.is_covered() )
        {
          while( cnt < num_vars + extra_vars )
            leaves[cnt++] = ntk.get_null();
    
          auto [func, care] = compute_function( n, leaves );
          cuts.emplace_back( n, leaves, func, care );
          found = true;
        }
        if ( found )
          break;
      }
      if ( found )
        break;
    }

    return found;
  }

  template<class WinManager_t>
  bool greedy_support_selection( node const& n, WinManager_t const& window_manager )
  {
    /* initialize the spfd manager */
    signature_t const & sign = node_to_sim[n];
    spfd_manager.init( sign );
    /* prepare the data */
    std::array<node, num_vars+extra_vars> leaves;
    auto min_cost = spfd_manager.get_num_edges();
    node candidate = ntk.get_null();
    int cnt = 0;
    while( !spfd_manager.is_covered() && !spfd_manager.is_saturated() )
    {
      bool found_new = false;
      window_manager.foreach_divisor( [&]( auto d, auto i ){
        if( !ps.preserve_depth || ( analyzer.get_time( d ) < analyzer.get_time( n ) ) )
        {
          signature_t const& sig = node_to_sim[d];
          auto const cost_div = spfd_manager.evaluate( sig );
          if( cost_div < min_cost )
          {
            found_new = true;
            min_cost = cost_div;
            candidate = d;
          }
        }
      } );
      if( !found_new )
      {
        return false;
      }
      signature_t const& sig = node_to_sim[candidate];
      spfd_manager.update( sig );
      leaves[cnt++] = candidate;
    }
    if( spfd_manager.is_covered() )
    {
      while( cnt < (num_vars+extra_vars) )
      {
        leaves[cnt++] = ntk.get_null();
      }
      auto [func, care] = compute_function( n, leaves );
      cuts.emplace_back( n, leaves, func, care );
    }
    return true;
  }



  std::tuple<truth_table_t, truth_table_t> compute_function( node const& n, std::array<node, num_vars+extra_vars> const& leaves )
  {
    signature_t const& root_signature = node_to_sim[n];
    uint32_t num_minterms = 1u << (num_vars+extra_vars);
    truth_table_t func;
    truth_table_t care;
    truth_table_t temp;
    func ^= func;
    care ^= care;
    signature_t mask;
    for( uint32_t m = 0; m < num_minterms; ++m )
    {
      mask |= ~mask;
      temp |= ~temp;
      for( uint32_t v = 0; v < (num_vars+extra_vars); ++v )
      {
        if( !ntk.is_null( leaves[v] ) )
        {
          signature_t const& sig = node_to_sim[leaves[v]];
          bool polarity = ( ( m >> v ) & 0x1 ) > 0;
          mask &= polarity ? sig : ~sig;
          temp &= polarity ? proj_functions[v] : ~proj_functions[v];
        }
      }
      if( kitty::count_ones( mask ) != 0 )
      {
        care |= temp;
      }
      if( kitty::count_ones( mask & root_signature ) > 0 )
      {
        func |= temp;
      }
    }
    return std::make_tuple( func, care );
  }

public:
  bool has_candidate()
  {
    return best_reward > zero_cost;
  }

  cost get_reward() const
  {
    return best_reward;
  }

  mapped_list get_list() const
  {
    return best_list;
  }

  std::array<node, num_vars+extra_vars> get_leaves() const
  {
    return best_leaves;
  }

  void add_pattern( std::vector<bool> cex )
  {
    uint32_t bit = num_cex - 1;
    ntk.foreach_pi( [&]( auto pi, auto ii ){
      if( cex[ii] )
        kitty::set_bit( node_to_word[pi], bit );
      else
        kitty::clear_bit( node_to_word[pi], bit );
    } );
    num_cex++;
    if( ( num_cex % 64u ) == 0 )
    {
      simulate_network( node_to_word, list_simulator_word );
    }
    uint32_t const num_blocks = 1u << ( cube_size - 6u );
    uint32_t block = ( num_cex / 64 ) % num_blocks;
  }

private:
  Ntk & ntk;
  Database_t & database;
  Analyzer_t & analyzer;
  simula_optimizer_params const& ps;
  simula_optimizer_stats & st;
  network_cuts_t cuts;
  spfd_manager_t<cube_size, num_vars+extra_vars> spfd_manager;

  incomplete_node_map<signature_t, Ntk> node_to_sim;
  incomplete_node_map<word_t, Ntk> node_to_word;

  aig_list_simulator<signature_t> list_simulator_sign;
  aig_list_simulator<word_t> list_simulator_word;

  uint32_t num_cex;

  std::array<truth_table_t, num_vars+extra_vars> proj_functions;
  std::array<node, num_vars+extra_vars> best_leaves;
  mapped_list best_list;
  cost best_reward;

  /* events */
  std::shared_ptr<typename network_events<Ntk>::add_event_type> add_event;
  std::shared_ptr<typename network_events<Ntk>::modified_event_type> modified_event;

};

} /* namespace mockturtle */