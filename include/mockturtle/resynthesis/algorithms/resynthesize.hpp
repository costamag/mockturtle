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
  \file resynthesize.hpp
  \brief Inplace rewrite

  \author Andrea Costamagna
*/

#pragma once

#include "../windowing/window_manager.hpp"
#include "struct_optimizer.hpp"
#include "window_optimizer.hpp"
#include "simula_optimizer.hpp"
#include "rewire_optimizer.hpp"
#include "../analyzers/analyzers_common.hpp"
//#include "../analyzers/homo_pxx_analyzer.hpp"
//#include "../analyzers/homo_xpx_analyzer.hpp"
#include "../analyzers/homo_xgx_analyzer.hpp"
#include "../placement/placement_analyzer.hpp"
//#include "../analyzers/homo_xxa_analyzer.hpp"
#include "../analyzers/homo_pla_analyzer.hpp"
//#include "../analyzers/homo_ppa_analyzer.hpp"
#include "../database/homo_database.hpp"
#include "../../algorithms/circuit_validator.hpp"
#include "../../views/depth_view.hpp"

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
struct resynthesis_params
{
  resynthesis_params()
  {
    /* 0 < Cut limit < 16 */
    cut_enumeration_ps.cut_limit = 8;
    cut_enumeration_ps.minimize_truth_table = true;
  }

  analyzer_params analyzer_ps;

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

  bool collect_mode { false };
  int decompose {0};

  uint32_t max_clauses{ 1000 };
  uint32_t max_nodes{ 10000 };
  int32_t odc_levels{ 0 };
  uint32_t conflict_limit{ 10000 };
  uint32_t random_seed{5};
  uint32_t fanout_limit{100};
  double min_reward {10};//0.3 for iscas
  uint32_t max_level_simula{20};

  /*! \brief If true try cut-rewriting with structural cuts */
  bool try_struct{ false };

  /*! \brief If true try window-based rewriting with non-structural cuts */
  bool try_window{ false };

  /*! \brief If true try simulation-guided rewriting with non-structural cuts */
  bool try_simula{ false };

  /*! \brief If true try fanin rewiring */
  bool try_rewire{ false };

  bool dynamic_database{ false };
  std::vector<double> pattern;
};

/*! \brief Statistics for rewrite.
 *
 * The data structure `rewrite_stats` provides data collected by running
 * `rewrite`.
 */
struct resynthesis_stats
{
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{ 0 };

  /*! \brief Expected gain. */
  uint32_t estimated_gain{ 0 };

  /*! \brief Candidates */
  uint32_t candidates{ 0 };

  uint32_t num_struct{ 0 };
  uint32_t num_window{ 0 };
  uint32_t num_simula{ 0 };
  uint32_t num_rewire{ 0 };

  void report() const
  {
    std::cout << fmt::format( "[i] total time       = {:>5.2f} secs\n", to_seconds( time_total ) );
    std::cout << fmt::format( "    num struct       = {:5d}\n", num_struct );
    std::cout << fmt::format( "    num window       = {:5d}\n", num_window );
    std::cout << fmt::format( "    num simula       = {:5d}\n", num_simula );
    std::cout << fmt::format( "    num rewire       = {:5d}\n", num_rewire );
  }
};

namespace detail
{

template<class Ntk, class Database_t, class Analyzer_t, uint32_t num_vars = 4, uint32_t cube_size = 10, uint32_t max_leaves = 10, uint32_t extra_vars = 0>
class resynthesize_impl
{
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;
  using network_cuts_t = dynamic_network_cuts<Ntk, num_vars+extra_vars, true, cut_enumeration_rewrite_cut>;
  using node_data = typename Ntk::storage::element_type::node_type;
  using truth_table_t = kitty::static_truth_table<num_vars+extra_vars>;
  using validator_t = circuit_validator<Ntk, bill::solvers::bsat2, false, true, false>;//last is false

  using struct_optimizer_t = struct_optimizer<Ntk, Database_t, Analyzer_t, num_vars, extra_vars>;
  using window_optimizer_t = window_optimizer<Ntk, Database_t, Analyzer_t, num_vars, max_leaves, extra_vars>;
  using simula_optimizer_t = simula_optimizer<Ntk, Database_t, Analyzer_t, num_vars, cube_size, extra_vars>;
  using rewire_optimizer_t = rewire_optimizer<Ntk, Analyzer_t, num_vars, max_leaves>;

public:
  resynthesize_impl( Ntk & ntk, Database_t & database, resynthesis_params ps, resynthesis_stats& st )
      : ntk( ntk ),
        database( database ),
        analyzer( ntk, ps.analyzer_ps ),
        win_manager( ntk ),
        ps( ps ),
        st( st )
  {
  }

  ~resynthesize_impl()
  {
  }

  void run()
  {
    /* specifications for structural exploration */
    struct_optimizer_params struct_ps;
    struct_ps.cut_enumeration_ps.cut_size = num_vars+extra_vars;
    struct_ps.dynamic_database = ps.dynamic_database;
    struct_ps.collect_mode = ps.collect_mode;
    struct_ps.decompose = ps.decompose;
    struct_ps.preserve_depth = ps.preserve_depth;
    struct_optimizer_stats struct_st;
    struct_optimizer_t struct_opto( ntk, database, analyzer, struct_ps, struct_st );

    /* specifications for window-based exploration */
    window_optimizer_params window_ps;
    window_ps.preserve_depth = ps.preserve_depth;
    window_ps.decompose = ps.decompose;
    window_optimizer_stats window_st;

    window_optimizer_t window_opto( ntk, database, analyzer, window_ps, window_st );

    /* specifications for window-based exploration */
    simula_optimizer_params simula_ps;
    simula_ps.preserve_depth = ps.preserve_depth;
    simula_ps.decompose = ps.decompose;
    simula_optimizer_stats simula_st;

    simula_optimizer_t simula_opto( ntk, database, analyzer, simula_ps, simula_st );
    
    /* specifications for rewire-based exploration */
    rewire_optimizer_params rewire_ps;
    rewire_ps.preserve_depth = ps.preserve_depth;
    rewire_optimizer_stats rewire_st;

    rewire_optimizer_t rewire_opto( ntk, analyzer, rewire_ps, rewire_st );

    analyzer.choose_gates();
    validator_t validator( ntk, { ps.max_clauses, ps.odc_levels, ps.conflict_limit, ps.random_seed } );


int cnt = 0;
int cnt_max = ps.max_nodes;
//std::cout << "cnt_max  = ";
//std::cin >> cnt_max;

    analyzer.foreach_gate( [&]( auto n ) {
if ( cnt++ >= cnt_max )
  return true;
      /* Skip nodes which cannot result in optimization */
      if ( skip_node_analysis( n ) )
      {
        return true;
      }
      /* Build and run analysis on window, if required by any heuristic */
      if ( !window_analysis( n ) )
      {
        return true;
      }

      if( ps.collect_mode )
      {
        struct_opto.analyze( n );
        return true;
      }

      if( ( ps.try_window || ps.try_simula || ps.try_rewire ) && !win_manager.is_valid() )
      {
        window_opto.reset();
        rewire_opto.reset();
      }
      else 
      {
        if( ps.try_window )
        {
          window_opto.analyze( n, win_manager );
        }
        if( ps.try_rewire )
        {
          rewire_opto.analyze( n, win_manager );
        }
      }

      if( ps.try_rewire && win_manager.is_valid() )
      {
        rewire_opto.analyze( n, win_manager );
        if( rewire_opto.has_candidate() )
        {
          auto const list = rewire_opto.get_list();
          auto const leaves = rewire_opto.get_leaves();
          double const reward = rewire_opto.get_reward();
          if( reward > ps.min_reward )
          {
            signal fnew = list.insert( ntk, leaves );
            //double Ubef = analyzer.compute_cost();
            ntk.substitute_node( n, fnew );
            analyzer.remove_mffc_info();
            //double Uaft = analyzer.compute_cost();
            //cnt++;
            update( n, fnew );
            //double Uaft2 = analyzer.compute_cost();
            //std::cout << Ubef << "->" << Uaft << "->" << Uaft2 << std::endl;
            st.num_rewire += 1;
            return true;
          }
        }
      }

      if( ps.try_struct )
      {
        struct_opto.analyze( n );
      }

      if( window_opto.has_candidate() || struct_opto.has_candidate() )
      {
        bool do_struct = ps.try_struct && struct_opto.has_candidate();
        bool do_window = ps.try_window && window_opto.has_candidate();
        std::vector<std::pair<int, double>> cands;
        cands.push_back( std::make_pair( 0, do_struct ? struct_opto.get_reward() : -1 ) );
        cands.push_back( std::make_pair( 1, do_window ? window_opto.get_reward() : -1 ) );
        std::stable_sort( cands.begin(), cands.end(), [&]( auto a, auto b ){
          return a.second > b.second;
        } );

        signal fnew;

        if ( cands[0].second > ps.min_reward )
        {
          //std::cout << cands[0].second << std::endl;
          if ( cands[0].first == 0 )
          {
            auto const list = struct_opto.get_list();
            auto const leaves = struct_opto.get_leaves();
            fnew = list.insert( ntk, leaves );
            st.num_struct += 1;
          }
          else if ( cands[0].first == 1 )
          {
            auto const list = window_opto.get_list();
            auto const leaves = window_opto.get_leaves();
            fnew = list.insert( ntk, leaves );
            st.num_window += 1;
          }
          ntk.substitute_node( n, fnew );
          update( n, fnew );
          if( cands[0].first == 0 )
          {
            struct_opto.clear_cuts_fanout_rec( ntk.get_node( fnew ) );
          }
  
          return true;
        }
      }

      if( ps.try_simula && ( ntk.level( n ) < ps.max_level_simula ) )
      {
        simula_opto.analyze( n, win_manager );
        if( simula_opto.has_candidate() && ( simula_opto.get_reward() > ps.min_reward ) )
        {
          auto const list = simula_opto.get_list();
          auto const leaves = simula_opto.get_leaves();
          auto const new_index = ntk.size();
          signal fnew = list.insert( ntk, leaves );
          node nnew = ntk.get_node( fnew );
          {
            auto valid = validator.validate( ntk.make_signal( n ), fnew );
            if( valid )
            {
              if( *valid && n != fnew.index )
              {
                ntk.substitute_node( n, fnew );
                st.num_simula += 1;
              }
              else
              {
                simula_opto.add_pattern( validator.cex );
              }
              update( n, fnew );
              if( ps.try_simula )
                validator.update();
            }
          }
          return true;
        }
      }

      return true;
    } );
  }

  /*! \brief Checks if the node should be analyzed for optimization or skipped */
  bool skip_node_analysis( node const& n )
  {
    if( ntk.fanout_size( n ) > ps.fanout_limit )
    {
      return true;
    }
    if( ntk.fanout_size( n ) <= 0 )
    {
      return true;
    }
    if( ntk.is_pi( n ) || ntk.is_constant( n ) )
    {
      return true;
    }
    if ( ntk.is_dead( n ) )
    {
      return true;
    }
    if ( !analyzer.is_interesting( n ) )
    {
      return true;
    }
    return false;
  }

  /*! \brief Perform window analysis if required by any heuristic. Return false if failure */
  bool window_analysis( node const& n )
  {
    win_manager.reset();
    win_manager.run( n );

    if constexpr( Analyzer_t::pass_manager )
    {
      if( !win_manager.is_valid() )
      {
        return false;
      }
      analyzer( n, win_manager );
    }
    else
      analyzer( n );
    return true;
  }

  void update( node n, signal fnew )
  {
    analyzer.update( n, fnew );
  }

private:
  Ntk & ntk;
  Database_t & database;
  Analyzer_t analyzer;
  window_manager<Ntk, max_leaves> win_manager;
  //validator_t validator;
  resynthesis_params ps;
  resynthesis_stats & st;

};

} /* namespace detail */

/*! \brief Boolean rewrite.
 *
 * This algorithm rewrites enumerated cuts using new network structures from a database.
 * The algorithm performs changes in-place and keeps the substituted structures dangling
 * in the network.
 *
 * **Required network functions:**
 * - `get_node`
 * - `size`
 * - `make_signal`
 * - `foreach_gate`
 * - `substitute_node`
 * - `clear_visited`
 * - `clear_values`
 * - `fanout_size`
 * - `set_value`
 * - `foreach_node`
 *
 * \param ntk Input network (will be changed in-place)
 * \param library Exact library containing pre-computed structures
 * \param ps Rewrite params
 * \param pst Rewrite statistics
 * \param cost_fn Node cost function (a functor with signature `uint32_t(Ntk const&, node<Ntk> const&)`)
 */
#if 0
template<class Ntk, class Database_t, class Analyzer_t>
void resynthesize( Ntk& ntk, Database_t & database )
{
  resynthesis_params ps;
  using dNtk = depth_view<Ntk>;
  dNtk dntk{ ntk };
  detail::resynthesize_impl<dNtk, Database_t, Analyzer_t> p( dntk, database, ps );
  p.run();
}

template<class Ntk, class Database_t, uint32_t num_vars>
void power_resynthesize( Ntk& ntk, Database_t & database, resynthesis_params ps = {} )
{
  using dNtk = depth_view<Ntk>;
  dNtk dntk{ ntk };
  using Analyzer_t = homo_xpx_analyzer<dNtk, num_vars>;
  detail::resynthesize_impl<dNtk, Database_t, Analyzer_t, num_vars> p( dntk, database, ps );
  p.run();
}
#endif
//template<class Ntk, class Database_t, uint32_t num_vars, uint32_t num_steps, uint32_t cube_size, uint32_t max_leaves, uint32_t extra_vars = 0>
//void delay_resynthesize( Ntk& ntk, Database_t & database, resynthesis_params ps = {} )
//{
//  using dNtk = depth_view<Ntk>;
//  dNtk dntk{ ntk };
//
//  using Analyzer_t = homo_pxx_analyzer<dNtk, num_vars, num_steps, cube_size, max_leaves, extra_vars>;
//  resynthesis_stats st;
//  std::cout << "initialize" << std::endl;
//  detail::resynthesize_impl<dNtk, Database_t, Analyzer_t, num_vars, cube_size, max_leaves, extra_vars> p( dntk, database, ps, st );
//  std::cout << "run" << std::endl;
//  p.run();
//  std::cout << "report" << std::endl;
//  st.report();
//}
//
//template<class Ntk, class Database_t, uint32_t num_vars, uint32_t num_steps, uint32_t cube_size, uint32_t max_leaves, uint32_t extra_vars=0>
//void area_resynthesize( Ntk& ntk, Database_t & database, resynthesis_params ps = {} )
//{
//  using dNtk = depth_view<Ntk>;
//  dNtk dntk{ ntk };
//  using Analyzer_t = homo_xxa_analyzer<dNtk, num_vars, num_steps, cube_size, max_leaves>;
//  resynthesis_stats st;
//  detail::resynthesize_impl<dNtk, Database_t, Analyzer_t, num_vars, cube_size, max_leaves, extra_vars> p( dntk, database, ps, st );
//  p.run();
//  st.report();
//}
//
//template<class Ntk, class Database_t, uint32_t num_vars, uint32_t num_steps, uint32_t cube_size, uint32_t max_leaves>
//void power_resynthesize( Ntk& ntk, Database_t & database, resynthesis_params ps = {} )
//{
//  using dNtk = depth_view<Ntk>;
//  dNtk dntk{ ntk };
//  using Analyzer_t = homo_xpx_analyzer<dNtk, num_vars, num_steps, cube_size, max_leaves>;
//  resynthesis_stats st;
//  detail::resynthesize_impl<dNtk, Database_t, Analyzer_t, num_vars, cube_size, max_leaves> p( dntk, database, ps, st );
//  p.run();
//  st.report();
//}

template<class Ntk, class Database_t, uint32_t num_vars, uint32_t num_steps, uint32_t cube_size, uint32_t max_leaves>
void glitch_resynthesize( Ntk& ntk, Database_t & database, resynthesis_params ps = {} )
{
  using dNtk = depth_view<Ntk>;
  dNtk dntk{ ntk };
  using Analyzer_t = homo_xgx_analyzer<dNtk, num_vars, num_steps, cube_size, max_leaves>;
  resynthesis_stats st;
  detail::resynthesize_impl<dNtk, Database_t, Analyzer_t, num_vars, cube_size, max_leaves> p( dntk, database, ps, st );
  p.run();
  st.report();
}

template<class Ntk, class Database_t, uint32_t num_vars, uint32_t num_steps, uint32_t cube_size, uint32_t max_leaves, uint32_t Size1, uint32_t Size2, uint32_t extra_vars=0>
void place_resynthesize( Ntk& ntk, Database_t & database, resynthesis_params ps = {} )
{
  std::cout << "Start" << std::endl;
  using dNtk = depth_view<Ntk>;
  dNtk dntk{ ntk };
  using Analyzer_t = homo_pla_analyzer<dNtk, num_vars, num_steps, cube_size, max_leaves, Size1, Size2, extra_vars>;
  resynthesis_stats st;
  detail::resynthesize_impl<dNtk, Database_t, Analyzer_t, num_vars, cube_size, max_leaves, extra_vars> p( dntk, database, ps, st );
  p.run();
  st.report();
  std::cout << "Done" << std::endl;
}

#if 0
template<class Ntk, class Database_t, uint32_t num_vars>
void area_resynthesize( Ntk& ntk, Database_t & database, resynthesis_params ps = {} )
{
  using dNtk = depth_view<Ntk>;
  dNtk dntk{ ntk };
  using Analyzer_t = homo_xxa_analyzer<dNtk, num_vars>;
  detail::resynthesize_impl<dNtk, Database_t, Analyzer_t, num_vars> p( dntk, database, ps );
  p.run();
}

template<class Ntk, class Database_t, uint32_t num_vars, uint32_t num_steps, uint32_t cube_size, uint32_t max_leaves>
void ppa_resynthesize( Ntk& ntk, Database_t & database, resynthesis_params ps = {} )
{
  using dNtk = depth_view<Ntk>;
  dNtk dntk{ ntk };
  using Analyzer_t = homo_ppa_analyzer<dNtk, num_vars, num_steps, cube_size, max_leaves>;
  detail::resynthesize_impl<dNtk, Database_t, Analyzer_t, num_vars, cube_size, max_leaves> p( dntk, database, ps );
  p.run();
}

template<class Ntk, class Database_t, uint32_t num_vars, uint32_t num_steps, uint32_t cube_size, uint32_t max_leaves>
void ppa_resynthesize( Ntk& ntk, Database_t & database, std::vector<typename Ntk::node> const& nodes, resynthesis_params ps = {} )
{
  using dNtk = depth_view<Ntk>;
  dNtk dntk{ ntk };
  using Analyzer_t = homo_ppa_analyzer<dNtk, num_vars, num_steps, cube_size, max_leaves>;
  detail::resynthesize_impl<dNtk, Database_t, Analyzer_t, num_vars, cube_size, max_leaves> p( dntk, database, ps );
  p.run( nodes );
}
#endif


} /* namespace mockturtle */