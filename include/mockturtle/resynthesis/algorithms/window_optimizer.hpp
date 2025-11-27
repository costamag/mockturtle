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
  \file window_optimizer.hpp
  \brief Inplace rewrite

  \author Andrea Costamagna
*/

#pragma once

#include "../boolean/dependency_cut.hpp"
#include "../networks/bindings.hpp"
#include "../mapped_list_simulator.hpp"
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
 struct window_optimizer_params
 {
   window_optimizer_params()
   {
     /* 0 < Cut limit < 16 */
     cut_enumeration_ps.cut_limit = 8;
    cut_enumeration_ps.minimize_truth_table = true;
  }

  /*! \brief Cut enumeration parameters. */
  cut_enumeration_params cut_enumeration_ps{};
  
  /*! \brief If true, candidates are only accepted if they do not increase logic depth. */
  bool preserve_depth{ true };
  
  /*! \brief Allow rewrite with multiple structures */
  bool allow_multiple_structures{ true };
  
  /*! \brief Allow zero-gain substitutions */
  bool allow_zero_gain{ false };
  
  int decompose { 0 };
  
  /*! \brief Use satisfiability don't cares for optimization. */
  bool use_dont_cares{ false };
  bool use_old {false};
  
  /*! \brief Window size for don't cares calculation. */
  uint32_t window_size{ 8u };
  
  /*! \brief Be verbose. */
  bool verbose{ false };
  
  /* 1 -> consider all the don't cares, 0 -> only one function */
  uint64_t effort{ 1 };
  bool greedy{ true };
  bool time{ false };
  /* maximum number of fanouts to be used as a new divisor */
  uint32_t fanout_limit{100};
};

/*! \brief Statistics for rewrite.
*
* The data structure `rewrite_stats` provides data collected by running
* `rewrite`.
*/
struct window_optimizer_stats
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

template<class Ntk, class Database_t, class Analyzer_t, uint32_t num_vars, uint32_t max_leaves, uint32_t extra_vars = 0>
class window_optimizer
{
  using node = typename Ntk::node;
  using cost = typename Analyzer_t::cost;
  static constexpr cost zero_cost = Analyzer_t::zero_cost;
  static constexpr uint32_t extra_space = 20;
  static constexpr uint32_t offset = 2;
  using cut_t = static_dependency_cut<Ntk, num_vars + extra_vars>;
  using network_cuts_t = std::vector<cut_t>;
  using truth_table_t = kitty::static_truth_table<num_vars + extra_vars>;
  using signature_t = kitty::static_truth_table<max_leaves>;
  using bindings_t = bindings;
  bool VERBOSEw = false;
  
public:
  window_optimizer( Ntk & ntk, Database_t & database, Analyzer_t & analyzer, window_optimizer_params const& ps, window_optimizer_stats& st )
    : ntk( ntk ),
      database( database ),
      analyzer( analyzer ),
      ps( ps ),
      st( st ),
      cuts( extra_space ),
      spfd_manager(),
      sims( max_leaves + extra_space ),
      cst(),
      node_to_sim( ntk ),
      list_simulator(),
      mlist_simulator( ntk.get_bindings() ),
      proj_functions(),
      best_leaves(),
      best_list( num_vars + extra_vars ),
      best_reward( zero_cost )
  {
    init();
  }

  void init()
  {
    /* reserve and initialize the signatures for exhaustive window simulation */
    sims.resize( max_leaves + offset );
    sims[0] = sims[0]^sims[0];
    sims[1] = ~sims[0];
    for( int i = 0; i < max_leaves; ++i )
    {
      kitty::create_nth_var( sims[offset+i], i );
    }

    for( int i = 0; i < num_vars + extra_vars; ++i )
    {
      kitty::create_nth_var( proj_functions[i], i );
    }

    for( int i = 0; i < num_vars; ++i )
    {
      kitty::create_nth_var( _xs[i], i );
    }
  }

  void reset()
  {
    best_reward = zero_cost;
  }

  template<class WinManager_t>
  void analyze( node const& n, WinManager_t const& window_manager )
  {
    std::array<node, num_vars + extra_vars> leaves;
    std::array<node, num_vars + extra_vars> perm_leaves;
    std::array<node, num_vars + extra_vars> times;
    std::array<node, num_vars + extra_vars> permutation;
    best_reward = zero_cost;
    /* terminate if the window is not valid */
    if( !window_manager.is_valid() )
    {
      return;
    }

    simulate_window( window_manager );

    if( !compute_dependency_cuts( n, window_manager ) )
    {
      return;
    }

    for ( auto const& cut : cuts )
    {
      /* extract the cut leaves and functionality */
      auto const leaves = cut.get_leaves();

      auto cost_curr = analyzer.evaluate( n, window_manager.leaves, leaves );
if( VERBOSEw )
  std::cout << "cost curr " << cost_curr << std::endl;

      auto const& func = cut.get_function();
      auto const& care = cut.get_careset();

if ( ps.decompose == 0 )
{
      std::vector<uint32_t> dont_cares;
      uint32_t supp_size = 0;
      for ( auto l : leaves )
      {
        supp_size += ( !ntk.is_null( l ) ) ? 1 : 0;
//std::cout << l << " ";
      }
//std::cout << std::endl;

//kitty::print_binary( care ); std::cout << std::endl;
      kitty::static_truth_table<num_vars+extra_vars> tmp;
      for( int m = 0; m < ( 1u << supp_size ); ++m )
      {
        tmp |= ~tmp;
        for ( int v = 0; v < leaves.size(); ++v )
        {
          if ( !ntk.is_null( leaves[v] ) )
          {
            tmp &= ( ( ( m >> v ) & 0x1 ) == 0x1 ) ? proj_functions[v] : ~proj_functions[v];
          }
        }

        if( kitty::is_const0( care&tmp ) )
          dont_cares.push_back( m );
      }
      uint64_t num_funcs = 1u << dont_cares.size();
      num_funcs = std::min( num_funcs, ps.effort );

      uint64_t num_ones = kitty::count_ones( func & care );
      uint64_t num_zeros = kitty::count_ones( (~func) & care );

      //std::cout << num_funcs << std::endl;
      for( uint64_t m = 0; m < num_funcs; ++m )
      {
        truth_table_t onset = func & care;

        for( int dc = 0; dc < dont_cares.size(); ++dc )
        {
          tmp |= ~tmp;
          for ( int v = 0; v < leaves.size(); ++v )
          {
            if ( !ntk.is_null( leaves[v] ) )
            {
              tmp &= ( ( ( dont_cares[dc] >> v ) & 0x1 ) == 0x1 ) ? proj_functions[v] : ~proj_functions[v];
            }
          }

          if( (( m >> dc ) & 0x1) == 0x1 )
            onset |= tmp;
          else
            onset &= ~tmp;
        }

        int cnt = 0;
        for( int j{0}; j < num_vars + extra_vars; ++j )
        {
          if( !ntk.is_null( leaves[j] ) )
          {
            cnt++;
            uint32_t index_leaf = node_to_sim[leaves[j]];
            uint32_t index_root = node_to_sim[n];
            auto tt_leaf = sims[index_leaf];
            auto tt_root = sims[index_root];
            if( kitty::equal( tt_root, tt_leaf ) )
            {
              cost cost_cand = analyzer.evaluate_rewiring( leaves[j], n );
              cost reward = analyzer.compare( cost_curr, cost_cand );
              if( ( reward > best_reward ) )
              {
                mapped_list list( num_vars + extra_vars );
                list.create_po( list.pi_at( j ) );
                best_reward = reward;
                best_list = list;
                best_leaves = leaves;
              }
            }
          }
        }

        if constexpr ( extra_vars > 0 )
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
              indices_lut.push_back( j++ );
              std::cout << n << " ";
            }
          }
          //std::cout << j << "/" << num_vars+extra_vars << std::endl;
          while ( times_lut_det.size() < ( num_vars + extra_vars + 2 )  )
            times_lut_det.push_back( std::numeric_limits<double>::max() );

          mapped_list cand_list( num_vars + extra_vars );
          truth_table_list_t<num_vars, num_vars+extra_vars> tt_list;
          uint32_t lit = std::numeric_limits<uint32_t>::max();
//  std::cout << "W" << std::endl;          
          lit = tt_list.delay_oriented_shannon( func, times_lut, indices_lut );

          if ( tt_list.is_null( lit ) )
            continue;


          uint32_t out_lit=lit;
          tt_list.foreach_gate( [&]( auto nd ){
            auto const& tt = nd.func;
//kitty::print_binary( tt ); std::cout << std::endl;
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
//std::cout << "children[i] " << children[i] << std::endl;
                auto lit = tt_list.is_null( children[i] ) ? children[i] : lut_to_list[children[i]];
//if ( !tt_list.is_null( children[i] ) )
//  std::cout << "lut to list " << lut_to_list[children[i]] << std::endl;
//std::cout << "lit " << lit << std::endl;
                perm_leaves[permutation[i]] = lit;
                times_cut[permutation[i]] = tt_list.is_null( children[i] ) ? std::numeric_limits<double>::max() : times_lut_det[children[i]]; // remove the offset of the constant 0 and 1
//if ( !tt_list.is_null( children[i] ) )
//  std::cout << "times_cut[permutation[i]] " << times_cut[permutation[i]] << std::endl;
              }

              double min_cost = std::numeric_limits<double>::max();
              mapped_list best_list_s( num_vars );
              database.foreach_list( entry, [&]( auto list, auto index ){
                auto const symm = database.get_symm( entry );
                list.match_times( times_cut, symm );
                //list.print();
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

//std::cout << " out lit " << out_lit << std::endl;
              lut_to_list.push_back( out_lit );
            }
            else
            {
//std::cout << "no match found" << std::endl; 
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
        else
        {
          /* extract the database entries */
          auto res = database.get_entry( onset );
          if( res )
          {
            auto perm = std::get<0>( *res );
            auto entry = std::get<1>( *res );

            std::array<uint8_t, num_vars + extra_vars> permutation;
            for ( auto j = 0u; j < num_vars + extra_vars; ++j )
            {
              permutation[perm[j]] = j;
            }

            for( int i = 0; i < num_vars + extra_vars; ++i )
            {
              auto const l = leaves[i];
              perm_leaves[permutation[i]] = l;
              times[permutation[i]] = !ntk.is_null( l ) ? analyzer.get_time( l ) : std::numeric_limits<double>::max();
            }

            database.foreach_list( entry, [&]( auto list, auto index ){
              
              auto const symm = database.get_symm( entry );
              list.match_times( times, symm );
  if( VERBOSEw )
    list.print();
              cost cost_cand = analyzer.evaluate( list, perm_leaves, n );
  if( VERBOSEw )
    std::cout << "(w) cost cand " << cost_cand << std::endl;
              cost reward = analyzer.compare( cost_curr, cost_cand );
  if( VERBOSEw )
    std::cout << "(w) reward    " << reward << std::endl;

              if( ( reward > best_reward ) )
              {
                best_reward = reward;
  if( VERBOSEw )
    std::cout << "Best reward " << best_reward << std::endl;
                best_list = list;
                best_leaves = perm_leaves;
              }
            } );
          }
  //        else
  //        {
  //          std::cout << "no entry implementing this function" << std::endl;
  //        }
        }
        //if( m >= ps.effort*num_funcs )
        //{
        //  m = num_funcs - 1;
        //}
      }
    }
else
{
    std::vector<double> times_lut;
    std::vector<double> times_lut_det{ 0, 0 };
    std::vector<uint32_t> indices_lut;
    std::vector<uint32_t> indices_list;
    std::vector<uint32_t> lut_to_list;
    std::vector<kitty::static_truth_table<max_leaves>> lut_sims;
    std::vector<kitty::static_truth_table<max_leaves>> isims;
    lut_sims.emplace_back(); // 0
    lut_sims.push_back( ~lut_sims[0] ) ; // 1
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
        indices_lut.push_back( j++ );

        uint32_t const index = node_to_sim[n];
        lut_sims.push_back( sims[index] );
        //std::cout << n << " ";
      }
    }
    //std::cout << j << "/" << num_vars+extra_vars << std::endl;
    
    while ( lut_sims.size() < ( num_vars + extra_vars + 2 ) )
      lut_sims.emplace_back();

    while ( times_lut_det.size() < ( num_vars + extra_vars + 2 )  )
      times_lut_det.push_back( std::numeric_limits<double>::max() );

    mapped_list cand_list( num_vars + extra_vars );
    truth_table_list_t<num_vars, num_vars+extra_vars> tt_list;
    uint32_t lit = std::numeric_limits<uint32_t>::max();
// std::cout << "W" << std::endl;           
    lit = tt_list.delay_oriented_shannon_dc( func, care, times_lut, indices_lut );//
// std::cout << "Wed " << lit << std::endl;          
    if ( tt_list.is_null( lit ) )
    continue;

    std::array<uint32_t, num_vars> perm_leaves;
    std::array<double, num_vars>   times_cut;
    uint32_t out_lit = lit;
    tt_list.foreach_gate( [&]( auto nd ){
      auto const& tt = nd.func;
      auto const& lit = nd.lit;
      auto const& children = nd.children;
      uint32_t nvars = 0;
      for ( int v = 0; v < num_vars; ++v )
        nvars += ( children[v] != std::numeric_limits<uint32_t>::max() ) ? 1 : 0;
      
      isims.clear();
      for ( int i = 0; i < num_vars; ++i )
      {
        if ( children[i] != std::numeric_limits<uint32_t>::max() )
          isims.push_back( lut_sims[children[i]] );
      }

      uint32_t num_minterms = 1u << nvars;
      kitty::static_truth_table<max_leaves> tmp_l;
      std::vector<uint32_t> dcs;
      for ( uint32_t m = 0; m < num_minterms; ++m )
      {
        tmp_l |= ~tmp_l;
        for ( int v = 0; v < nvars; ++v )
        {
          auto const c = children[v];
          tmp_l &= ( ( ( m >> v ) & 0x1 ) == 0x1 ) ? lut_sims[c] : ~lut_sims[c];
        }
        if ( kitty::count_ones( tmp_l ) == 0 )
        {
          dcs.push_back( m );
        }
      }
      uint32_t num_funcs = std::min( ps.effort, static_cast<uint64_t>(1u << dcs.size()) );
      kitty::static_truth_table<num_vars> onset;
      kitty::static_truth_table<num_vars> tmp;
      
      mapped_list best_list_s( num_vars );
      double min_cost = std::numeric_limits<double>::max();

      for ( uint64_t j = 0; j < num_funcs; ++j )
      {
        onset = tt;
        for ( int k = 0; k < dcs.size(); ++k )
        {
          tmp |= ~tmp;
          auto const& dc = dcs[k];
          for ( int v = 0; v < nvars; ++v )
          {
            tmp &= ( ( ( dc >> v ) & 0x1 ) == 0x1 ) ? _xs[v] : ~_xs[v];
          }
          if ( ( ( j >> k ) & 0x1 ) == 0x1 )
            onset ^= tmp;
        }
        auto res = database.get_entry( onset );
        if ( res )
        {
          auto perm = std::get<0>( *res );
          auto entry = std::get<1>( *res );
          std::array<uint32_t, num_vars> permutation{ std::numeric_limits<uint32_t>::max() };
          for ( int i = 0; i < perm.size(); ++i )
          {
            permutation[perm[i]] = i;
          }

          for ( int i = 0; i < num_vars; ++i )
          {
            auto lit = tt_list.is_null( children[i] ) ? children[i] : lut_to_list[children[i]];
            perm_leaves[permutation[i]] = lit;
            times_cut[permutation[i]] = tt_list.is_null( children[i] ) ? std::numeric_limits<double>::max() : times_lut_det[children[i]]; // remove the offset of the constant 0 and 1
          }

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
        }
        else
        {
          best_reward = zero_cost;
          return;
        }
      }
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

      lut_to_list.push_back( out_lit );
      lut_sims.emplace_back();
      mlist_simulator( lut_sims[lut_sims.size()-1], best_list_s, isims );
    } );
    cand_list.create_po( out_lit );
    best_list = cand_list;
    std::array<node, num_vars + extra_vars> large_leaves;
    best_leaves = leaves;
    cost cost_cand = analyzer.evaluate( best_list, best_leaves, n );
    cost reward = analyzer.compare( cost_curr, cost_cand );
    //std::cout << " " << cost_curr << "->" << cost_cand << std::endl;
    best_reward = reward;
  }
    
  }

  }

private:
  template<class WinManager_t>
  bool simulate_window( WinManager_t const& window_manager )
  {
    node_to_sim.resize();
    sims.resize( max_leaves + offset );
    window_manager.foreach_leaf( [&]( node const& l, uint32_t index ){
      node_to_sim[l] = offset + index; 
if ( VERBOSEw )
{
//  kitty::print_binary( sims[node_to_sim[l]] ); std::cout << " " << l << std::endl;
}
    } );

    window_manager.foreach_divisor( [&]( node const& d, uint32_t index ){
      if( index >= window_manager.leaves.size() )
      {
        std::vector<signature_t*> sim_ptrs;
        ntk.foreach_fanin( d, [&]( auto fi, auto ii ){
          node const ni = ntk.get_node( fi );
          uint32_t const index = node_to_sim[ni];
          sim_ptrs.push_back( &sims[index] );
        } );
        uint32_t const index = sims.size();
        node_to_sim[d] = index;
        sims.emplace_back();
        auto const alist = ntk.get_binding( d ).list;
        list_simulator( sims[index], alist, sim_ptrs );
if ( VERBOSEw )
{
//  kitty::print_binary( sims[index] ); std::cout << " " << d << std::endl;
}
      }
    } );

    window_manager.foreach_mffc( [&]( node const& m, uint32_t index ){
      if ( !window_manager.is_leaf( m ) )
      {
        std::vector<signature_t*> sim_ptrs;
        ntk.foreach_fanin( m, [&]( auto fi, auto ii ){
          node const ni = ntk.get_node( fi );
          uint32_t const index = node_to_sim[ni];
          sim_ptrs.push_back( &sims[index] );
        } );
        uint32_t const index = sims.size();
        node_to_sim[m] = index;
        sims.emplace_back();
        auto const alist = ntk.get_binding( m ).list;
        list_simulator( sims[index], alist, sim_ptrs );
      }
    } );
    return true;
  }

  template<class WinManager_t>
  bool compute_dependency_cuts( node n, WinManager_t const& window_manager )
  {
    /* prepare the data */
    cuts.clear();

    /* call greedy support selection */
    if( ps.greedy == true )
      return greedy_support_selection( n, window_manager );
    else if( ps.time == true )
      return times_support_selection( n, window_manager );
    else
      return pivot_support_selection( n, window_manager );
  }


  template<class WinManager_t>
  bool time_support_selection( node const& n, WinManager_t const& window_manager )
  {
    bool found = false;
    node piv;
    double min_delay = std::numeric_limits<double>::max();
    window_manager.foreach_divisor( [&]( auto d, auto i ){
      double const div_delay = analyzer.get_arrival( d );
      if( div_delay < min_delay )
      {
        min_delay = analyzer.get_arrival( d );
        piv = d;
      }
    });

      /* initialize the spfd manager */
      uint32_t index = node_to_sim[n];
      signature_t const & sign = sims[index];
      spfd_manager.init( sign );

      uint32_t const piv_index = node_to_sim[piv];
      signature_t const& sig = sims[piv_index];
      auto min_cost = spfd_manager.evaluate( sig );
      spfd_manager.update( sig );
      std::array<node, num_vars + extra_vars> leaves;
      leaves[0] = piv;
      int cnt = 1;

      /* prepare the data */
      node candidate = ntk.get_null();
      while( !spfd_manager.is_covered() && !spfd_manager.is_saturated() )
      {
        bool found_new = false;
        window_manager.foreach_divisor( [&]( auto d, auto i ){
          if( ntk.fanout_size( d ) <= ps.fanout_limit && ntk.fanout_size( d ) > 0 )
          {
            if( !ps.preserve_depth || ( analyzer.get_time( d ) < (analyzer.get_time( n )-100) ) )
            {
              uint32_t const index = node_to_sim[d];
              signature_t const& sig = sims[index];
              auto const cost_div = spfd_manager.evaluate( sig );
              if( cost_div < min_cost )
              {
                found_new = true;
                min_cost = cost_div;
                candidate = d;
              }
            }
          }
        } );
        uint32_t const index = node_to_sim[candidate];
        signature_t const& sig = sims[index];
        spfd_manager.update( sig );
        leaves[cnt++] = candidate;
      }
      if( spfd_manager.is_covered() )
      {
        for ( int i = 0; i < cnt; ++i )
          std::cout << analyzer.get_time( leaves[cnt] ) << " ";
        std::cout << std::endl;

        while( cnt < (num_vars + extra_vars) )
          leaves[cnt++] = ntk.get_null();
  
        auto [func, care] = compute_function( n, leaves );
        cuts.emplace_back( n, leaves, func, care );
        found = true;
      }

    return found;
  }


  template<class WinManager_t>
  bool greedy_support_selection( node const& n, WinManager_t const& window_manager )
  {
    /* initialize the spfd manager */
    uint32_t index = node_to_sim[n];
    signature_t const & sign = sims[index];
    spfd_manager.init( sign );
    /* prepare the data */
    std::array<node, num_vars + extra_vars> leaves;
    auto min_cost = spfd_manager.get_num_edges();
    node candidate = ntk.get_null();
    int cnt = 0;
    while( !spfd_manager.is_covered() && !spfd_manager.is_saturated() )
    {
      bool found_new = false;
      window_manager.foreach_divisor( [&]( auto d, auto i ){
        if( ntk.fanout_size( d ) <= ps.fanout_limit && ntk.fanout_size( d ) > 0 )
        {
          if( !ps.preserve_depth || ( analyzer.get_time( d ) < analyzer.get_time( n ) ) )
          {
            uint32_t const index = node_to_sim[d];
            signature_t const& sig = sims[index];
            auto const cost_div = spfd_manager.evaluate( sig );
            if( cost_div < min_cost )
            {
              found_new = true;
              min_cost = cost_div;
              candidate = d;
            }
          }
        }
      } );
      if( !found_new )
      {
        return false;
      }
      uint32_t const index = node_to_sim[candidate];
      signature_t const& sig = sims[index];
      spfd_manager.update( sig );
      leaves[cnt++] = candidate;
    }
    if( spfd_manager.is_covered() )
    {
//      std::vector<std::pair<node, double>> node_times;
//      window_manager.foreach_divisor( [&]( auto d, auto i ){
//        node_times.push_back( std::make_pair( d, analyzer.get_time( d ) ) );
//      } );
//      std::stable_sort( node_times.begin(), node_times.end(), [&](auto a, auto b)
//                                  {
//                                    return a.second < b.second;
//                                  });

      while( cnt < num_vars + extra_vars )
        leaves[cnt++] = ntk.get_null();

      auto [func, care] = compute_function( n, leaves );

//      int index = 0;
//      while( cnt < num_vars + extra_vars && index < node_times.size() )
//      {
//        if( std::find( leaves.begin(), leaves.end(), node_times[index].first ) == leaves.end() )
//        {
//          leaves[cnt++] = node_times[index].first;
//        }
//        index++;
//      }

      cuts.emplace_back( n, leaves, func, care );
    }

    return true;
  }

  template<class WinManager_t>
  bool pivot_support_selection( node const& n, WinManager_t const& window_manager )
  {
    bool found = false;
    window_manager.foreach_divisor( [&]( auto piv, auto i ){
      if( found )
        return true;
      if( ntk.fanout_size( piv ) > ps.fanout_limit || ntk.fanout_size( piv ) == 0 )
        return true;

      /* initialize the spfd manager */
      uint32_t index = node_to_sim[n];
      signature_t const & sign = sims[index];
      spfd_manager.init( sign );

      uint32_t const piv_index = node_to_sim[piv];
      signature_t const& sig = sims[piv_index];
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
        window_manager.foreach_divisor( [&]( auto d, auto i ){
          if( ntk.fanout_size( d ) <= ps.fanout_limit && ntk.fanout_size( d ) > 0 )
          {
            if( !ps.preserve_depth || ( analyzer.get_time( d ) < (analyzer.get_time( n ) ) ) )
            {
              uint32_t const index = node_to_sim[d];
              signature_t const& sig = sims[index];
              auto const cost_div = spfd_manager.evaluate( sig );
              if( cost_div < min_cost )
              {
                found_new = true;
                min_cost = cost_div;
                candidate = d;
              }
            }
          }
        } );
        if( !found_new )
        {
          return true;
        }
        uint32_t const index = node_to_sim[candidate];
        signature_t const& sig = sims[index];
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
      return true;
    } );

    return found;
  }
  
  template<class WinManager_t>
  bool times_support_selection( node const& n, WinManager_t const& window_manager )
  {
    std::array<node, num_vars + extra_vars> leaves;
    if ( window_manager.leaves.size() <= ( num_vars + extra_vars ) )
    {
      for ( int i = 0; i < window_manager.leaves.size(); ++i )
        leaves[i]=window_manager.leaves[i];
      for ( int i = window_manager.leaves.size(); i < ( num_vars+extra_vars ); ++i )
        leaves[i]=ntk.get_null();
      auto [func, care] = compute_function( n, leaves );
      cuts.emplace_back( n, leaves, func, care );
      return true;
    }
    bool found = false;
    double min_time = std::numeric_limits<double>::max();
    double max_time = analyzer.get_time( n );
    std::vector<node> sorted_divs;
    window_manager.foreach_divisor( [&]( auto d, auto i ){
      if ( ( analyzer.get_time( d ) < (analyzer.get_time( n ) - 20) ) )
      {
        sorted_divs.push_back( d );
        auto const& d_time = analyzer.get_time( d );
        min_time = std::min( min_time, d_time );
        max_time = std::max( max_time, d_time );
      }
    } );
    max_time = analyzer.get_time( n ) - 20;
    std::sort( sorted_divs.begin(), sorted_divs.end(), [&]( auto na, auto nb ){
      return analyzer.get_time( na ) < analyzer.get_time( nb );
    } );
    if ( sorted_divs.size() <= (num_vars + extra_vars) )
      return false;

    uint32_t steps = 10;
    double delta = ( max_time - min_time ) / ((double)steps);

    for ( uint32_t s = 0; s < steps; ++s )
    {
      double threshold = min_time + s * delta;
      //std::cout << min_time << " #" << s << "/" << sorted_divs.size() << " " << threshold << " " << max_time << std::endl;
      for ( int ip = 0; ip < sorted_divs.size(); ++ip )
      {
        auto piv = sorted_divs[ip];
        if( found )
          continue;
        if( (ntk.fanout_size( piv ) > ps.fanout_limit) || (ntk.fanout_size( piv ) == 0) )
          continue;

        if ( analyzer.get_time( piv ) >= threshold )
        {
          //std::cout << "Break at " << ip << std::endl;
          break;
        }

        /* initialize the spfd manager */
        uint32_t index = node_to_sim[n];
        signature_t const & sign = sims[index];
        spfd_manager.init( sign );

        uint32_t const piv_index = node_to_sim[piv];
        signature_t const& sig = sims[piv_index];
        auto min_cost = spfd_manager.evaluate( sig );
        spfd_manager.update( sig );
        leaves[0] = piv;
        int cnt = 1;

        /* prepare the data */
        node candidate = ntk.get_null();
        while( !spfd_manager.is_covered() && !spfd_manager.is_saturated() && ( cnt < ( num_vars+extra_vars) ) )
        {
          bool found_new = false;
          for ( int jp = 0; jp < sorted_divs.size(); ++jp )
          {
            auto d = sorted_divs[jp];
            if( ntk.fanout_size( d ) <= ps.fanout_limit && ntk.fanout_size( d ) > 0 )
            {
              if( !ps.preserve_depth || ( analyzer.get_time( d ) < threshold ) )
              {
                uint32_t const index = node_to_sim[d];
                signature_t const& sig = sims[index];
                auto const cost_div = spfd_manager.evaluate( sig );
                if( cost_div < min_cost )
                {
                  found_new = true;
                  min_cost = cost_div;
                  candidate = d;
                }
              }
             // if ( ( analyzer.get_time( d ) >= threshold ) )
             //   std::cout << "stop at " << jp << std::endl;
            }
          }
          if( !found_new )
          {
            break;
          }
          uint32_t const index = node_to_sim[candidate];
          signature_t const& sig = sims[index];
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



  std::tuple<truth_table_t, truth_table_t> compute_function( node const& n, std::array<node, num_vars + extra_vars> const& leaves )
  {
    uint32_t const index = node_to_sim[n];
    signature_t const& root_signature = sims[index];

    uint32_t num_minterms = 1u << num_vars + extra_vars;
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
      for( uint32_t v = 0; v < num_vars + extra_vars; ++v )
      {
        if( !ntk.is_null( leaves[v] ) )
        {
          uint32_t index = node_to_sim[leaves[v]];
          signature_t const& sig = sims[index];
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

  std::array<node, num_vars + extra_vars> get_leaves() const
  {
    return best_leaves;
  }

  template<class WinManager_t>
  void print_simulation( WinManager_t const& window_manager )
  {
    std::cout << "leaves" << std::endl;
    window_manager.foreach_leaf( [&]( auto l, auto i ){
      std::cout << l << " ";
      kitty::print_binary( sims[node_to_sim[l]] );
      std::cout << std::endl;
    } );

    std::cout << "divs" << std::endl;
    window_manager.foreach_divisor( [&]( auto d, auto i ){
      std::cout << d << " ";
      kitty::print_binary( sims[node_to_sim[d]] );

      std::cout << " ";
      if( !ntk.is_pi( d ) )
      {
        auto g = ntk.get_binding( d );
        std::cout << g.expression << " ";
        ntk.foreach_fanin( d, [&]( auto fi, auto ii ){
          std::cout << ntk.get_node( fi ) << " ";
        } );
      }
      else
      {
        std::cout << " PI";
      }
      std::cout << std::endl;
    } );

    std::cout << "mffc" << std::endl;
    window_manager.foreach_mffc( [&]( auto m, auto i ){
      std::cout << m << " ";
      kitty::print_binary( sims[node_to_sim[m]] );

      std::cout << " ";
      if( !ntk.is_pi( m ) )
      {
        auto g = ntk.get_binding( m );
        std::cout << g.expression << " ";
        ntk.foreach_fanin( m, [&]( auto fi, auto ii ){
          std::cout << ntk.get_node( fi ) << " ";
        } );
      }
      else
      {
        std::cout << " PI";
      }
      std::cout << std::endl;
    } );
  }

private:
  Ntk & ntk;
  Database_t & database;
  Analyzer_t & analyzer;
  window_optimizer_params const& ps;
  window_optimizer_stats & st;
  network_cuts_t cuts;
  spfd_manager_t<max_leaves, num_vars + extra_vars> spfd_manager;
  std::vector<signature_t> sims;
  cut_enumeration_stats cst;
  incomplete_node_map<uint32_t, Ntk> node_to_sim;
  aig_list_simulator<signature_t> list_simulator;
  mapped_list_simulator<bindings_t, signature_t> mlist_simulator;

  std::array<truth_table_t, num_vars + extra_vars> proj_functions;
  std::array<kitty::static_truth_table<num_vars>, num_vars> _xs;
  std::array<node, num_vars + extra_vars> best_leaves;
  mapped_list best_list;
  cost best_reward;

};

} /* namespace mockturtle */