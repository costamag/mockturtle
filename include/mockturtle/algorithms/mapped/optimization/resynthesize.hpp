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

#include "../dependencies/rewire_dependencies.hpp"
#include "../evaluators/area_resyn_evaluator.hpp"
#include "../evaluators/evaluators_utils.hpp"
#include "../windowing/window_manager.hpp"
#include "../windowing/window_simulator.hpp"
// #include "../dependencies/struct_dependencies.hpp"
// #include "../dependencies/window_dependencies.hpp"
// #include "../dependencies/simula_dependencies.hpp"
// #include "../../../algorithms/circuit_validator.hpp"
#include "../database/mapped_database.hpp"

#include <fmt/format.h>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>
#include <optional>

namespace mockturtle
{

/*! \brief Statistics for rewrite.
 *
 * The data structure `rewrite_stats` provides data collected by running
 * `rewrite`.
 */
struct resynthesis_stats
{

  window_manager_stats window_st;
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

struct default_resynthesis_params
{
  evaluator_params evaluator_ps;

  /*! \brief If true, candidates are only accepted if they do not increase logic depth. */

  struct window_manager_params : default_window_manager_params
  {
    static constexpr uint32_t max_num_leaves = 6u;
    bool preserve_depth = false;
    int32_t odc_levels = 0u;
    uint32_t skip_fanout_limit_for_divisors = 100u;
    uint32_t max_num_divisors{ 128 };
  };
  window_manager_params window_manager_ps;

  static constexpr bool do_strashing = true;
  /*! \brief Use satisfiability don't cares for optimization. */
  static constexpr bool use_dont_cares = false;

  /*! \brief If true try fanin rewiring */
  static constexpr bool try_rewire = false;

  /*! \brief If true try cut-rewriting with structural cuts */
  static constexpr bool try_struct = false;

  /*! \brief If true try window-based rewriting with non-structural cuts */
  static constexpr bool try_window = false;

  /*! \brief If true try simulation-guided rewriting with non-structural cuts */
  static constexpr bool try_simula = false;

  /*! \brief Activates lazy man's synthesis when set to true */
  static constexpr bool dynamic_database = false;

  /*! \brief Maximum number of leaves of the window */
  static constexpr uint32_t max_num_leaves = 6u;
  /*! \brief Cube size for the signatures in simulation-guided resubstitution */
  static constexpr uint32_t num_vars_sign = 10u;
  /*! \brief Maximum number of leaves in the dependency cuts */
  static constexpr uint32_t max_cuts_size = 6u;
  /*! \brief Maximum cube size exactly represented with SPFDs */
  static constexpr uint32_t max_cube_spfd = 12u;

  /*! \brief Maximum fanout size for a node to be optimized*/
  static constexpr uint32_t fanout_limit = 12u;
};

namespace detail
{

template<class Ntk, typename Database, typename Evaluator, typename Params = default_resynthesis_params>
class resynthesize_impl
{
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;
  using cut_t = dependency_cut_t<Ntk, Params::max_cuts_size>;

  // using validator_t = circuit_validator<Ntk, bill::solvers::bsat2, false, true, false>; // last is false

  struct rewire_params : default_rewire_params
  {
    static constexpr uint32_t max_cuts_size = Params::max_cuts_size;
  };

  //  struct window_params : default_window_params
  //  {
  //    static constexpr uint32_t num_vars_sign = Params::num_vars_sign;
  //    static constexpr uint32_t max_cuts_size = Params::max_cuts_size;
  //    static constexpr uint32_t max_cube_spfd = Params::max_cube_spfd;
  //  };
  //
  //  struct struct_params : default_struct_params
  //  {
  //    static constexpr uint32_t num_vars_sign = Params::num_vars_sign;
  //    static constexpr uint32_t max_cuts_size = Params::max_cuts_size;
  //  };

  using rewire_dependencies_t = rewire_dependencies<Ntk, rewire_params>;
  // using struct_dependencies_t = struct_dependencies<Ntk, custom_struct_params>;
  // using window_dependencies_t = window_dependencies<Ntk, custom_window_params>;
  // using simula_dependencies_t = simula_dependencies<Ntk, custom_simula_params>;

public:
  resynthesize_impl( Ntk& ntk, Database& database, Params ps, resynthesis_stats& st )
      : ntk_( ntk ),
        database_( database ),
        evaluator_( ntk, ps.evaluator_ps ),
        win_manager_( ntk, ps.window_manager_ps, st.window_st ),
        win_simulator_( ntk ),
        ps_( ps ),
        st_( st )
  {
  }

  void run()
  {
    /* specifications for rewire-based exploration */
    rewire_dependencies_t rewire( ntk_ );

    /* specifications for structural exploration */
    //    auto [struct_ps, struct_st] = init_struct_specs();
    //    struct_dependencies_t struct_dependencies( ntk, evaluator, struct_ps, struct_st );
    //
    //    /* specifications for window-based exploration */
    //    auto [window_ps, window_st] = init_window_specs();
    //    window_dependencies_t window_dependencies( ntk, evaluator, window_ps, window_st );
    //
    //    /* specifications for window-based exploration */
    //    auto [simula_ps, simula_st] = init_simula_specs();
    //    simula_dependencies_t simula_dependencies( ntk, evaluator, simula_ps, simula_st );
    //    validator_t validator( ntk, { ps_.max_clauses, ps_.odc_levels, ps_.conflict_limit, ps_.random_seed } );

    evaluator_.foreach_gate( [&]( auto n ) {
      /* Skip nodes which cannot result in optimization */
      if ( skip_node( n ) )
        return true;

      /* Build and run analysis on window */
      window_analysis( n );
      if ( ps_.try_rewire && win_manager_.is_valid() )
      {
        auto const& win_leaves = win_manager_.get_leaves();
        rewire.run( win_manager_, win_simulator_ );
        std::optional<cut_t> best_cut;
        double best_reward = 0;

        rewire.foreach_cut( [&]( auto& cut, auto i ) {
          auto const& cut_leaves = cut.leaves;
          auto const reward = evaluator_.evaluate_rewiring( n, cut_leaves, win_leaves );
          if ( reward > best_reward )
          {
            best_reward = reward;
            best_cut = std::make_optional( cut );
          }
        } );
        if ( best_cut )
        {
          auto const ids = ntk_.get_binding_ids( n );
          auto const fnew = ntk_.template create_node<Params::do_strashing>( ( *best_cut ).leaves, ids );
          auto const nnew = ntk_.get_node( fnew );
          std::vector<signal> fs;
          ntk_.foreach_output( nnew, [&]( auto f ) {
            fs.push_back( f );
          } );
          ntk_.substitute_node( n, fs );
          return true;
        }
      }
      //
      //      if ( ps.try_struct )
      //      {
      //        struct_dependencies.analyze( n );
      //        identify_dependencies( struct_dependencies, cands );
      //      }
      //
      //      if ( ps_.try_window && win_manager_.is_valid() )
      //      {
      //        window_dependencies.analyze( n, win_manager_ );
      //        identify_dependencies( window_dependencies, cands );
      //      }
      //
      //      if ( cands.size() > 0 )
      //      {
      //        auto cand = evaluator_.choose( cands );
      //        if ( cand )
      //        {
      //          auto const n = (*cand).cut->root;
      //          auto const& fnew = insert( ntk_, (*cand).cut.leaves, (*cand).cut.list );
      //          ntk_.substitute_node( n, fnew );
      //        }
      //      }
      //      cands.clear();
      //
      //      if ( ps_.try_simula && ( ntk_.level( n ) < ps_.max_level_simula ) )
      //      {
      //        simula_dependencies.analyze( n, win_manager_ );
      //        identify_dependencies( simula_dependencies, cands );
      //        if ( cands.size() > 0 )
      //        {
      //          while( cands.size() > 0 )
      //          {
      //            auto cand = evaluator_.choose( cands );
      //            if ( cand )
      //            {
      //              auto const n = (*cand).cut->root;
      //              auto const& fnew = insert( ntk_, (*cand).cut.leaves, (*cand).cut.list );
      //              auto [valid, new_tt, abort] = validate( fnew, n );
      //              if ( validator.equal( fnew, n ) )
      //              {
      //                ntk_.substitute_node( n, fnew );
      //                return true;
      //              }
      //              else
      //              {
      //                ntk_.take_out_node( fnew );
      //                if ( abort )
      //                  delete_cand
      //                else
      //                  update_tt
      //              }
      //            }
      //          }
      //        }
      //      }

      return true;
    } );
  }

private:
  /*! \brief Checks if the node should be analyzed for optimization or skipped */
  bool skip_node( node const& n )
  {
    if ( ntk_.fanout_size( n ) > ps_.fanout_limit )
      return true;
    if ( ntk_.fanout_size( n ) <= 0 )
      return true;
    if ( ntk_.is_pi( n ) || ntk_.is_constant( n ) )
      return true;
    if ( ntk_.is_dead( n ) )
      return true;
    return false;
  }

  /*! \brief Perform window analysis if required by any heuristic. Return false if failure */
  void window_analysis( node const& n )
  {
    win_manager_.run( n );
    win_simulator_.run( win_manager_ );

    if constexpr ( Evaluator::pass_window )
    {
      if ( !win_manager_.is_valid() )
      {
        return;
      }
      if constexpr ( Evaluator::node_depend )
        evaluator_( n, win_manager_ );
      else
        evaluator_( win_manager_ );
    }
    else
    {
      if constexpr ( Evaluator::node_depend )
        evaluator_( n, win_manager_ );
    }
  }

  std::vector<double> get_times( std::vector<signal> const& leaves )
  {
    assert( Evaluator::has_arrival && "[e] The evaluator does not have the arrival tracker" );
    std::vector<double> times( leaves.size() );
    std::transform( leaves.begin(), leaves.end(), times.begin(),
                    [&]( auto const& f ) { return evaluator_.get_arrival( f ); } );
  }

private:
  Ntk& ntk_;
  Database& database_;
  Evaluator evaluator_;
  window_manager<Ntk, typename Params::window_manager_params> win_manager_;
  window_simulator<Ntk, Params::max_num_leaves> win_simulator_;
  Params ps_;
  resynthesis_stats& st_;
};

} /* namespace detail */

template<class Ntk, class Database, typename Params = default_resynthesis_params>
void area_resynthesize( Ntk& ntk, Database& database, Params ps = {}, resynthesis_stats* pst = nullptr )
{
  using Evaluator = area_resyn_evaluator<Ntk>;
  resynthesis_stats st;
  detail::resynthesize_impl<Ntk, Database, Evaluator, Params> p( ntk, database, ps, st );
  p.run();
  if ( pst != nullptr )
    *pst = st;
}

#if 0
template<class Ntk, class Database, uint32_t MaxNumVars, uint32_t num_steps, uint32_t CubeSize, uint32_t MaxNumLeaves>
void glitch_resynthesize( Ntk& ntk, Database & database, resynthesis_params ps = {} )
{
  using dNtk = depth_view<Ntk>;
  dNtk dntk{ ntk };
  using Evaluator = homo_xgx_evaluator<dNtk, MaxNumVars, num_steps, CubeSize, MaxNumLeaves>;
  resynthesis_stats st;
  detail::resynthesize_impl<dNtk, Database, Evaluator, MaxNumVars, CubeSize, MaxNumLeaves> p( dntk, database, ps, st );
  p.run();
  st.report();
}

template<class Ntk, class Database, uint32_t MaxNumVars, uint32_t num_steps, uint32_t CubeSize, uint32_t MaxNumLeaves>
void ppa_resynthesize( Ntk& ntk, Database & database, resynthesis_params ps = {} )
{
  using dNtk = depth_view<Ntk>;
  dNtk dntk{ ntk };
  using Evaluator = homo_ppa_evaluator<dNtk, MaxNumVars, num_steps, CubeSize, MaxNumLeaves>;
  detail::resynthesize_impl<dNtk, Database, Evaluator, MaxNumVars, CubeSize, MaxNumLeaves> p( dntk, database, ps );
  p.run();
}

template<class Ntk, class Database, uint32_t MaxNumVars, uint32_t num_steps, uint32_t CubeSize, uint32_t MaxNumLeaves>
void ppa_resynthesize( Ntk& ntk, Database & database, std::vector<typename Ntk::node> const& nodes, resynthesis_params ps = {} )
{
  using dNtk = depth_view<Ntk>;
  dNtk dntk{ ntk };
  using Evaluator = homo_ppa_evaluator<dNtk, MaxNumVars, num_steps, CubeSize, MaxNumLeaves>;
  detail::resynthesize_impl<dNtk, Database, Evaluator, MaxNumVars, CubeSize, MaxNumLeaves> p( dntk, database, ps );
  p.run( nodes );
}
#endif

} /* namespace mockturtle */