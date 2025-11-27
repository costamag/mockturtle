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
struct resub_optimizer_params
{
  resub_optimizer_params()
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

  /*! \brief Be verbose. */
  bool verbose{ false };

  /* 1 -> consider all the don't cares, 0 -> only one function */
  int effort{ 1 };
  bool greedy{false};
};

/*! \brief Statistics for rewrite.
 *
 * The data structure `rewrite_stats` provides data collected by running
 * `rewrite`.
 */
struct resub_optimizer_stats
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

template<class Ntk, class Database_t, class Analyzer_t, uint32_t num_vars, uint32_t max_leaves>
class resub_optimizer
{
public:
  using node = typename Ntk::node;
  using cost = typename Analyzer_t::cost;
  static constexpr cost zero_cost = Analyzer_t::zero_cost;
  static constexpr uint32_t extra_space = 20;
  static constexpr uint32_t offset = 2;
  using cut_t = static_dependency_cut<Ntk, num_vars>;
  using network_cuts_t = std::vector<cut_t>;
  using truth_table_t = kitty::static_truth_table<num_vars>;
  using signature_t = kitty::static_truth_table<max_leaves>;

  resub_optimizer( std::shared_ptr<Ntk> pntk, std::shared_ptr<Database_t> pdatabase, std::shared_ptr<Analyzer_t> panalyzer, resub_optimizer_params const& ps, resub_optimizer_stats& st )
    : pntk( pntk ),
      pDatabase( pdatabase ),
      pAnalyzer( panalyzer ),
      ps( ps ),
      st( st ),
      cuts( extra_space ),
      spfd_manager(),
      sims( max_leaves + extra_space ),
      cst(),
      node_to_sim( *pntk ),
      list_simulator(),
      proj_functions(),
      best_leaves(),
      best_list( num_vars ),
      best_reward( zero_cost )
  {
    init();
//    add_event = pntk->events().register_add_event( [&]( const auto& n ) {
//      node_to_sim.resize();
//    } );
//
//    modified_event = pntk->events().register_modified_event( [&]( const auto& n, auto children ) {
//      node_to_sim.resize();
//    } );
  }

  ~resub_optimizer()
  {
//    if ( add_event )
//    {
//      pntk->events().release_add_event( add_event );
//    }
//
//    if ( modified_event )
//    {
//      pntk->events().release_modified_event( modified_event );
//    }
  }


  void init()
  {
    /* reserve and initialize the signatures for exhaustive resub simulation */
    sims.resize( max_leaves + offset );
    sims[0] = sims[0]^sims[0];
    sims[1] = ~sims[0];
    for( int i = 0; i < max_leaves; ++i )
    {
      kitty::create_nth_var( sims[offset+i], i );
    }

    for( int i = 0; i < num_vars; ++i )
    {
      kitty::create_nth_var( proj_functions[i], i );
    }
  }

  void reset()
  {
    best_reward = zero_cost;
  }

  template<class WinManager_t>
  void analyze( node const& n, WinManager_t const& window_manager )
  {
    std::array<node, num_vars> leaves;
    std::array<node, num_vars> perm_leaves;
    std::array<node, num_vars> times;
    std::array<node, num_vars> permutation;
    best_reward = zero_cost;
    /* terminate if the window is not valid */
    if( !window_manager.is_valid() )
    {
      return;
    }
    if( window_manager.leaves.size() < 2 )
    {
      return;
    }

//  std::cout << "simulate" << std::endl;
    simulate_window( window_manager );
//window_manager.print();
//print_simulation( window_manager );

    if( !compute_dependency_cuts( n, window_manager ) )
    {
// std::cout << "no dependency cut" << std::endl;
      return;
    }

    for ( auto const& cut : cuts )
    {
      /* extract the cut leaves and functionality */
      auto const leaves = cut.get_leaves();
//for( int i = 0; i < num_vars; ++i )
//  std::cout << leaves[i] << " ";
//std::cout << std::endl;

      auto cost_curr = pAnalyzer->evaluate( n, window_manager.leaves );
//std::cout << "cost curr " << cost_curr << std::endl;

      auto const& func = cut.get_function();
      auto const& care = cut.get_careset();

      std::vector<uint32_t> dont_cares;
      for( int m = 0; m < func.num_bits(); ++m )
      {
        if( kitty::get_bit( care, m ) == 0 )
          dont_cares.push_back( m );
      }
      uint64_t const num_funcs = 1u << dont_cares.size();

      uint64_t num_ones = kitty::count_ones( func & care );
      uint64_t num_zeros = kitty::count_ones( (~func) & care );

      for( uint64_t m = 0; m < num_funcs; ++m )
      {
        truth_table_t onset = func & care;

        for( int dc = 0; dc < dont_cares.size(); ++dc )
        {
          if( (( m >> dc ) & 0x1) == 0x1 )
            kitty::set_bit( onset, dont_cares[dc] );
          else
            kitty::clear_bit( onset, dont_cares[dc] );
        }

        int cnt = 0;
        for( int j{0}; j < num_vars; ++j )
        {
          if( !pntk->is_null( leaves[j] ) )
          {
            cnt++;
            uint32_t index_leaf = node_to_sim[leaves[j]];
            uint32_t index_root = node_to_sim[n];
            auto tt_leaf = sims[index_leaf];
            auto tt_root = sims[index_root];
            if( kitty::equal( tt_root, tt_leaf ) )
            {
              cost cost_cand = pAnalyzer->evaluate_rewiring( leaves[j], n );
              cost reward = pAnalyzer->compare( cost_curr, cost_cand );
              if( ( reward > best_reward ) )
              {
                mapped_list list( num_vars );
                list.create_po( list.pi_at( j ) );
                best_reward = reward;
                best_list = list;
                best_leaves = leaves;
              }
            }
          }
        }

        /* extract the database entries */
        auto res = pDatabase->get_entry( onset );
        if( res )
        {
          auto perm = std::get<0>( *res );
          auto entry = std::get<1>( *res );

          std::array<uint8_t, num_vars> permutation;
          for ( auto j = 0u; j < num_vars; ++j )
          {
            permutation[perm[j]] = j;
          }

          for( int i = 0; i < num_vars; ++i )
          {
            auto const l = leaves[i];
            perm_leaves[permutation[i]] = l;
            times[permutation[i]] = !pntk->is_null( l ) ? pAnalyzer->get_time( l ) : std::numeric_limits<double>::max();
          }

          pDatabase->foreach_list( entry, [&]( auto list, auto index ){
            
            auto const symm = pDatabase->get_symm( entry );
            list.match_times( times, symm );

//list.print();
            cost cost_cand = pAnalyzer->evaluate( list, perm_leaves, n );
//std::cout << "cost cand " << cost_cand << std::endl;
            cost reward = pAnalyzer->compare( cost_curr, cost_cand );
//std::cout << "reward    " << reward << std::endl;

            if( ( reward > best_reward ) )
            {
              best_reward = reward;
              best_list = list;
              best_leaves = perm_leaves;
            }
          } );
        }
//        else
//        {
//          std::cout << "no entry implementing this function" << std::endl;
//        }
        if( ps.effort == 0 )
        {
          m = num_funcs - 1;
        }
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
    } );

    window_manager.foreach_divisor( [&]( node const& d, uint32_t index ){
      if( !window_manager.is_leaf( d ) )
      {
        std::vector<signature_t> sim_ptrs;
        pntk->foreach_fanin( d, [&]( auto fi, auto ii ){
          node const ni = pntk->get_node( fi );
          uint32_t const index = node_to_sim[ni];
          sim_ptrs.push_back( sims[index] );
        } );
        uint32_t const index = sims.size();
        node_to_sim[d] = index;
        sims.emplace_back();
        auto const alist = pntk->get_binding( d ).list;
        list_simulator( sims[index], alist, sim_ptrs );
      }
    } );

    window_manager.foreach_mffc( [&]( node const& m, uint32_t index ){
      if( !window_manager.is_leaf( m ) )
      {
        std::vector<signature_t> sim_ptrs;
        pntk->foreach_fanin( m, [&]( auto fi, auto ii ){
          node const ni = pntk->get_node( fi );
          uint32_t const index = node_to_sim[ni];
          sim_ptrs.push_back( sims[index] );
        } );
        uint32_t const index = sims.size();
        node_to_sim[m] = index;
        sims.emplace_back();
        auto const alist = pntk->get_binding( m ).list;
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
    else
      return pivot_support_selection( n, window_manager );
  }

  template<class WinManager_t>
  bool greedy_support_selection( node const& n, WinManager_t const& window_manager )
  {
    /* initialize the spfd manager */
    uint32_t index = node_to_sim[n];
    signature_t const & sign = sims[index];
    spfd_manager.init( sign );
    /* prepare the data */
    std::array<node, num_vars> leaves;
    auto min_cost = spfd_manager.get_num_edges();
    node candidate = pntk->get_null();
    int cnt = 0;
    while( !spfd_manager.is_covered() && !spfd_manager.is_saturated() )
    {
      bool found_new = false;
      window_manager.foreach_divisor( [&]( auto d, auto i ){
        if( !ps.preserve_depth || ( pAnalyzer->get_time( d ) < pAnalyzer->get_time( n ) ) )
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
//        node_times.push_back( std::make_pair( d, pAnalyzer->get_time( d ) ) );
//      } );
//      std::stable_sort( node_times.begin(), node_times.end(), [&](auto a, auto b)
//                                  {
//                                    return a.second < b.second;
//                                  });

      while( cnt < num_vars )
        leaves[cnt++] = pntk->get_null();

      auto [func, care] = compute_function( n, leaves );

//      int index = 0;
//      while( cnt < num_vars && index < node_times.size() )
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

      /* initialize the spfd manager */
      uint32_t index = node_to_sim[n];
      signature_t const & sign = sims[index];
      spfd_manager.init( sign );

      uint32_t const piv_index = node_to_sim[piv];
      signature_t const& sig = sims[piv_index];
      auto min_cost = spfd_manager.evaluate( sig );
      spfd_manager.update( sig );
      std::array<node, num_vars> leaves;
      leaves[0] = piv;
      int cnt = 1;

      /* prepare the data */
      node candidate = pntk->get_null();
      while( !spfd_manager.is_covered() && !spfd_manager.is_saturated() )
      {
        bool found_new = false;
        window_manager.foreach_divisor( [&]( auto d, auto i ){
          if( !ps.preserve_depth || ( pAnalyzer->get_time( d ) < pAnalyzer->get_time( n ) ) )
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
        while( cnt < num_vars )
          leaves[cnt++] = pntk->get_null();
  
        auto [func, care] = compute_function( n, leaves );
        cuts.emplace_back( n, leaves, func, care );
        found = true;
      }
      return true;
    } );

    return found;
  }

  std::tuple<truth_table_t, truth_table_t> compute_function( node const& n, std::array<node, num_vars> const& leaves )
  {
    uint32_t const index = node_to_sim[n];
    signature_t const& root_signature = sims[index];

    uint32_t num_minterms = 1u << num_vars;
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
      for( uint32_t v = 0; v < num_vars; ++v )
      {
        if( !pntk->is_null( leaves[v] ) )
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

  std::array<node, num_vars> get_leaves() const
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
      if( !pntk->is_pi( d ) )
      {
        auto g = pntk->get_binding( d );
        std::cout << g.expression << " ";
        pntk->foreach_fanin( d, [&]( auto fi, auto ii ){
          std::cout << pntk->get_node( fi ) << " ";
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
      if( !pntk->is_pi( m ) )
      {
        auto g = pntk->get_binding( m );
        std::cout << g.expression << " ";
        pntk->foreach_fanin( m, [&]( auto fi, auto ii ){
          std::cout << pntk->get_node( fi ) << " ";
        } );
      }
      else
      {
        std::cout << " PI";
      }
      std::cout << std::endl;
    } );
  }

public:
  signature_t get_sig( node n )
  {
    return sims[node_to_sim[n]];
  }

private:
  std::shared_ptr<Ntk> pntk;
  std::shared_ptr<Database_t> pDatabase;
  std::shared_ptr<Analyzer_t> pAnalyzer;
  resub_optimizer_params const& ps;
  resub_optimizer_stats & st;
  network_cuts_t cuts;
  spfd_manager_t<max_leaves, num_vars> spfd_manager;
  std::vector<signature_t> sims;
  cut_enumeration_stats cst;
  incomplete_node_map<uint32_t, Ntk> node_to_sim;
  aig_list_simulator<signature_t> list_simulator;
public:
  std::array<truth_table_t, num_vars> proj_functions;
  std::array<node, num_vars> best_leaves;
  mapped_list best_list;
  cost best_reward;
//  std::shared_ptr<typename network_events<Ntk>::add_event_type> add_event;
//  std::shared_ptr<typename network_events<Ntk>::modified_event_type> modified_event;
};

} /* namespace mockturtle */