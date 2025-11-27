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
  \file struct_optimizer.hpp
  \brief Inplace rewrite

  \author Andrea Costamagna
*/

#pragma once

#include "../../algorithms/cut_enumeration.hpp"
#include "../../algorithms/cut_enumeration/rewrite_cut.hpp"
#include "../networks/bindings.hpp"
#include "../mapped_list_simulator.hpp"
#include "../database/mapped_list.hpp"
#include "../boolean/spfd_manager.hpp"
#include "../boolean/truth_table_list.hpp"

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
struct struct_optimizer_params
{
  struct_optimizer_params()
  {
    /* 0 < Cut limit < 16 */
    cut_enumeration_ps.cut_limit = 8;//12
    /* otherwise creates problem in the presence of reconvergence (A^B)^(B^C) becomes {A,C} */
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

  /*! \brief Be verbose. */
  bool verbose{ false };
  int decompose { 0 };

  /*! \brief Use satisfiability don't cares for optimization. */
  bool dynamic_database{ false };
  bool collect_mode{ false };
  bool use_old{ false };
  uint64_t effort{ 128 };
};

/*! \brief Statistics for rewrite.
 *
 * The data structure `rewrite_stats` provides data collected by running
 * `rewrite`.
 */
struct struct_optimizer_stats
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

template<class Ntk, class Database_t, class Analyzer_t, uint32_t num_vars, uint32_t extra_vars = 0>
class struct_optimizer
{
  bool VERBOSEr = false;
  using node = typename Ntk::node;
  using cost = typename Analyzer_t::cost;
  static constexpr cost zero_cost = Analyzer_t::zero_cost;
  static constexpr cost max_cost = Analyzer_t::max_cost;
  using network_cuts_t = dynamic_network_cuts<Ntk, num_vars+extra_vars, true, cut_enumeration_rewrite_cut>;
  using cut_manager_t = detail::dynamic_cut_enumeration_impl<Ntk, num_vars+extra_vars, true, cut_enumeration_rewrite_cut>;
  using cut_t = typename network_cuts_t::cut_t;
  using node_data = typename Ntk::storage::element_type::node_type;
  using truth_table_t = kitty::static_truth_table<num_vars+extra_vars>;
  using bindings_t = bindings;
public:
  struct_optimizer( Ntk & ntk, Database_t & database, Analyzer_t & analyzer, struct_optimizer_params const& ps, struct_optimizer_stats& st )
    : ntk( ntk ),
      database( database ),
      analyzer( analyzer ),
      ps( ps ),
      st( st ),
      cuts( ntk.size() + ( ntk.size() >> 1 ) ),
      cst(),
      mlist_simulator( ntk.get_bindings() ),
      spfd_manager(),
      cut_manager( ntk, ps.cut_enumeration_ps, cst, cuts ),
      best_leaves(),
      best_list( num_vars + extra_vars ),
      best_reward( zero_cost )
    {
      init();
      for ( int i = 0; i < num_vars+extra_vars; ++i )
        kitty::create_nth_var( _xl[i], i );
      for ( int i = 0; i < num_vars; ++i )
        kitty::create_nth_var( _xs[i], i );
    }

  void init()
  {
    cut_manager.init_cuts2();
  }

  void analyze( node const& n )
  {
    /* structural optimization */
    best_reward = zero_cost;

    if ( ntk.fanout_size( n ) == 0u )
       return;

    cut_manager.clear_cuts( n );
    cut_manager.compute_cuts( n );
    std::array<node, num_vars+extra_vars> leaves;

    for ( auto & cut : cuts.cuts( ntk.node_to_index( n ) ) )
    {
      cost worst_reward = zero_cost;
      /* skip trivial cut */
      if ( ( cut->size() == 1 ) && *cut->begin() == ntk.node_to_index( n ) )
      {
        continue;
      }
      /* extract the cut leaves and functionality */
      fill_leaves( leaves, *cut );

//for ( auto l : leaves )
//  std::cout << l << " ";
//std::cout << std::endl;

if( VERBOSEr )
  print_subntk( n, leaves );
      auto const cost_curr = analyzer.evaluate( n, leaves, leaves );
//std::cout << "cost curr struct "<<cost_curr << " " << std::endl;

      auto const& func = cuts.truth_table( *cut );

//kitty::print_binary( func ); std::cout << std::endl;
if ( kitty::is_const0( func ) || kitty::is_const0(~func) )
  continue;

      if ( ps.decompose == 0 )
      {
        if constexpr ( extra_vars == 0 )
        {
          kitty::static_truth_table<num_vars> tt;
          kitty::extend_to_inplace( tt, func );

          /* extract the database entries */
          auto res = database.get_entry( tt );
          if ( res )
          {
            auto perm = std::get<0>( *res );
            auto entry = std::get<1>( *res );
  
            assert( perm.size() == num_vars );
  
            std::array<uint32_t, num_vars> permutation{ std::numeric_limits<uint32_t>::max() };
            for ( int i = 0; i < perm.size(); ++i )
            {
              permutation[perm[i]] = i;
            }
  
            std::array<node, num_vars> perm_leaves;
            std::array<double, num_vars> times;
            // assuming each variable has a permutation
            for ( int i = 0; i < num_vars; ++i )
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

              if ( ( reward > best_reward ) )
              {
    if( VERBOSEr )
      std::cout << "New best reward " << reward << " " << cost_curr << " " << cost_cand << std::endl;
                best_reward = reward;
                best_list = list;
                best_leaves = perm_leaves;
              }
              else if( reward < worst_reward )
              {
                worst_reward = reward;
              }
            } );
  
          if ( cut->size() == 0 )//|| ( cut->size() == 1 && *cut->begin() != ntk.node_to_index( n ) ) )
            break;
          //  break;
        }
        if ( ps.collect_mode || ps.dynamic_database )
        {
          database.add_entry( ntk, leaves, ntk.make_signal( n ) );
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
        for ( int i = 0; i < num_vars + extra_vars + 2; ++i )
        {
          indices_list.push_back( i );
          lut_to_list.push_back( i );
        }

        int j = 0;
        for ( auto leaf : *cut )
        {
          auto n = ntk.index_to_node( leaf );
          times_lut.push_back( analyzer.get_time( n ) );
          times_lut_det.push_back( analyzer.get_time( n ) );
//std::cout << j+2 << " times_lut_det : node " << n << " : " << times_lut_det[j+2] << std::endl;
          indices_lut.push_back( j++ );
        }
        while ( times_lut_det.size() < ( num_vars + extra_vars + 2 )  )
          times_lut_det.push_back( std::numeric_limits<double>::max() );

        mapped_list cand_list( num_vars + extra_vars );
        truth_table_list_t<num_vars, num_vars+extra_vars> tt_list;
//kitty::print_binary( func ); 
        uint32_t lit = std::numeric_limits<uint32_t>::max();
          lit = tt_list.delay_oriented_shannon( func, times_lut, indices_lut );
//std::cout << "S" << std::endl;
//tt_list.print();
          if ( tt_list.is_null( lit ) )
            continue;
        uint32_t out_lit;
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
//list.print();
//for( auto t : times_cut )
//  std::cout << t << " ";
//std::cout << std::endl;
              list.match_times( times_cut, symm );
list.print();
              cost cost_cand = analyzer.evaluate( list, times_cut ); // TODO define this function
              if ( cost_cand < min_cost )
              {
                min_cost = cost_cand;
                best_list_s = list;
              }
            } );
//std::cout << "Ss" << std::endl;
//best_list_s.print();

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
//best_list.print();
        std::array<node, num_vars + extra_vars> large_leaves;
        fill_leaves( large_leaves, *cut );
        best_leaves = large_leaves;
        cost cost_cand = analyzer.evaluate( best_list, large_leaves, n );
std::cout << cost_cand << std::endl;
        cost reward = analyzer.compare( cost_curr, cost_cand );
        best_reward = reward;
std::cout << "reward " << reward << std::endl;
      }
}
else
{
{
    std::vector<double> times_lut;
    std::vector<double> times_lut_det{ 0, 0 };
    std::vector<uint32_t> indices_lut;
    std::vector<uint32_t> indices_list;
    std::vector<uint32_t> lut_to_list;
    std::vector<kitty::static_truth_table<num_vars+extra_vars>> lut_sims;
    std::vector<kitty::static_truth_table<num_vars+extra_vars>> isims;
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
//std::cout << "put " << n << std::endl;
        lut_sims.push_back( _xl[i] );
      }
    }
    while ( lut_sims.size() < ( num_vars + extra_vars + 2 ) )
      lut_sims.emplace_back();

    while ( times_lut_det.size() < ( num_vars + extra_vars + 2 )  )
      times_lut_det.push_back( std::numeric_limits<double>::max() );

    mapped_list cand_list( num_vars + extra_vars );
    truth_table_list_t<num_vars, num_vars+extra_vars> tt_list;
    uint32_t lit = std::numeric_limits<uint32_t>::max();
    //std::cout << " decompose" << std::endl;
    lit = tt_list.delay_oriented_shannon_dc( func, func|~func, times_lut, indices_lut );//
    //std::cout << lit << " done decompose" << std::endl;
    //tt_list.print();
    if ( tt_list.is_null( lit ) )
      continue;

    std::array<uint32_t, num_vars> perm_leaves;
    std::array<double, num_vars>   times_cut;
    uint32_t out_lit=lit;
    tt_list.foreach_gate( [&]( auto nd ){
      auto const& tt = nd.func;
      auto const& lit = nd.lit;
      auto const& children = nd.children;
      uint32_t nvars = 0;
      for ( int v = 0; v < num_vars; ++v )
      {
// std::cout << v << std::endl;
// std::cout << children[v] << std::endl;
        nvars += ( children[v] != std::numeric_limits<uint32_t>::max() ) ? 1 : 0;
      }
      
      isims.clear();
      for ( int i = 0; i < num_vars; ++i )
      {
        if ( children[i] != std::numeric_limits<uint32_t>::max() )
          isims.push_back( lut_sims[children[i]] );
      }

      uint32_t num_minterms = 1u << nvars;
      kitty::static_truth_table<num_vars+extra_vars> tmp_l;
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
    best_reward = reward;
  }
}
    }

    if( has_candidate() )
    {
//      for( auto l : best_leaves )
//        std::cout << l << " ";
//      std::cout << std::endl;
    }
  }

  void clear_cuts_fanout_rec( node const& n )
  {
    ntk.foreach_fanout( n, [&]( auto const& g ) {
      auto const index = ntk.node_to_index( g );
      if ( cuts.cuts( index ).size() > 0 )
      {
        cut_manager.clear_cuts( g );
        clear_cuts_fanout_rec( g );
      }
    } );
  }

  void fill_leaves( std::array<node, num_vars+extra_vars>& leaves, cut_t const& cut )
  {
    int i=0;
    for ( auto const leaf : cut )
    {
      leaves[i++] = ntk.index_to_node( leaf );
    }

    while( i < num_vars+extra_vars )
    {
      leaves[i++] = ntk.get_null();
    }
  }

  void reset()
  {
  }

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

void print_subntk( node n, std::array<node, num_vars+extra_vars>& leaves )
{
  ntk.incr_trav_id();
  for( auto l : leaves )
  {
    std::cout << l << "_l ";
    if ( !ntk.is_null( l ) )
      ntk.set_visited( l, ntk.trav_id() );
  }
  std::cout << std::endl;
  print_rec( n );
}

void print_rec( node n )
{
  if( ntk.visited( n ) == ntk.trav_id() )
    return;
  ntk.set_visited( n, ntk.trav_id() );
  ntk.foreach_fanin( n, [&]( auto fi, auto ii ){
    auto ni = ntk.get_node( fi );
    print_rec( ni );
  } );

  ntk.foreach_fanin( n, [&]( auto fi, auto ii ){
    std::cout << ntk.get_node( fi ) << " ";
  } );
  if( !ntk.is_pi( n ) )
    std::cout << " : " << ntk.get_name( n ) << " " << n << std::endl;

}

private:
  Ntk & ntk;
  Database_t & database;
  Analyzer_t & analyzer;
  struct_optimizer_params const& ps;
  struct_optimizer_stats & st;
  network_cuts_t cuts;
  cut_manager_t cut_manager;
  cut_enumeration_stats cst;
  spfd_manager_t<num_vars + extra_vars, num_vars + extra_vars> spfd_manager;
  mapped_list_simulator<bindings_t, truth_table_t> mlist_simulator;
  std::vector<truth_table_t> sims;
  std::array<truth_table_t, num_vars+extra_vars> _xl;
  std::array<kitty::static_truth_table<num_vars>, num_vars> _xs;

  std::array<node, num_vars+extra_vars> best_leaves;
  mapped_list best_list;
  cost best_reward;

};

} /* namespace mockturtle */