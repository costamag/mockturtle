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
#include "resub_optimizer.hpp"
#include "../analyzers/homo_pxx_analyzer.hpp"
#include "../analyzers/homo_xpx_analyzer.hpp"
#include "../analyzers/homo_xgx_analyzer.hpp"
#include "../analyzers/homo_xxa_analyzer.hpp"
//#include "../analyzers/homo_ppa_analyzer.hpp"
#include "../database/homo_database.hpp"
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
struct resub_params
{
  resub_params()
  {
    /* 0 < Cut limit < 16 */
    cut_enumeration_ps.cut_limit = 12;
    cut_enumeration_ps.minimize_truth_table = true;
  }

  /*! \brief Cut enumeration parameters. */
  cut_enumeration_params cut_enumeration_ps{};

  /*! \brief If true, candidates are only accepted if they do not increase logic depth. */
  bool preserve_depth{ false };

  /*! \brief Allow rewrite with multiple structures */
  bool allow_multiple_structures{ true };

  /*! \brief Allow zero-gain resubs */
  bool allow_zero_gain{ false };

  /*! \brief Use satisfiability don't cares for optimization. */
  bool use_dont_cares{ false };

  /*! \brief Window size for don't cares calculation. */
  uint32_t window_size{ 8u };

  /*! \brief Be verbose. */
  bool verbose{ false };

  bool collect_mode { false };

  uint32_t max_clauses{ 1000 };
  uint32_t max_nodes{ 1000 };
  int32_t odc_levels{ 0 };
  uint32_t conflict_limit{ 1000 };
  uint32_t random_seed{5};

  bool try_struct{ false };
  bool try_window{ false };
  bool try_simula{ false };

  bool dynamic_database{ false };
//  bool only_critical{true};

};

/*! \brief Statistics for rewrite.
 *
 * The data structure `rewrite_stats` provides data collected by running
 * `rewrite`.
 */
struct resub_stats
{
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{ 0 };

  /*! \brief Expected gain. */
  uint32_t estimated_gain{ 0 };

  /*! \brief Candidates */
  uint32_t candidates{ 0 };

  uint32_t num_resub{ 0 };

  void report() const
  {
    std::cout << fmt::format( "[i] total time       = {:>5.2f} secs\n", to_seconds( time_total ) );
    std::cout << fmt::format( "    num resub        = {:5d}\n", num_resub );
  }
};

namespace detail
{

template<class Ntk, class Database_t, class Analyzer_t, uint32_t num_vars = 4, uint32_t cube_size = 10, uint32_t max_leaves = 10>
class resub_impl
{
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;
  using network_cuts_t = dynamic_network_cuts<Ntk, num_vars, true, cut_enumeration_rewrite_cut>;
  using cut_manager_t = detail::dynamic_cut_enumeration_impl<Ntk, num_vars, true, cut_enumeration_rewrite_cut>;
  using cut_t = typename network_cuts_t::cut_t;
  using node_data = typename Ntk::storage::element_type::node_type;
  using truth_table_t = kitty::static_truth_table<num_vars>;

  using struct_optimizer_t = struct_optimizer<Ntk, Database_t, Analyzer_t, num_vars>;
  using window_optimizer_t = window_optimizer<Ntk, Database_t, Analyzer_t, num_vars, max_leaves>;
  using resub_optimizer_t = resub_optimizer<Ntk, Database_t, Analyzer_t, num_vars, max_leaves>;
  using simula_optimizer_t = simula_optimizer<Ntk, Database_t, Analyzer_t, num_vars, cube_size>;

public:
  resub_impl( Ntk& ntk, Database_t & database, resub_params ps, resub_stats& st )
      : pntk( std::make_shared<Ntk>( ntk ) ),
        pDatabase( std::make_shared<Database_t>( database ) ),
        pAnalyzer( std::make_shared<Analyzer_t>( ntk ) ),
        pWin_manager( std::make_shared<window_manager<Ntk, max_leaves>>( ntk ) ),
        ps( ps ),
        st( st )
  {
  }

  void run()
  {
    /* specifications for window-based exploration */
    resub_optimizer_params resub_ps;
    resub_ps.preserve_depth = ps.preserve_depth;
    resub_optimizer_stats resub_st;
    resub_optimizer_t resub_opto( pntk, pDatabase, pAnalyzer, resub_ps, resub_st );

    convo_analyzer<Ntk, 10, cube_size, max_leaves> power( *pntk );
    power( 1 );
    pAnalyzer->choose_gates();

    int cnt = 0;
    bool found = false;
    pAnalyzer->foreach_gate( [&]( auto n ) {

      if( pntk->is_pi( n ) || pntk->is_constant( n ) )
      {
        return true;
      }

      if ( pntk->is_dead( n ) )
      {
        return true;
      }

      if( cnt++ > ps.max_nodes )
        return true;

      if constexpr( Analyzer_t::pass_manager )
      {
        if( !pWin_manager->run( n ) )
        {
          return true;
        }

        (*pAnalyzer)( n, *pWin_manager );
      }
      else
        (*pAnalyzer)( n );


      if( !pWin_manager->run( n ) )
      {
        std::cerr << "[e] window construction failed" << std::endl;
        resub_opto.reset();
      }
      else
      {

        resub_opto.analyze( n, *pWin_manager );

        typename resub_optimizer_t::signature_t tmp;
        std::vector<signal> children;
        pntk->foreach_fanin( n, [&]( auto fi, auto ii ){
          children.push_back( fi );
        } );

        auto tt = pntk->node_function( n );

        /* extract the SPF */
        std::vector<std::tuple<signal, int, double>> candidates;

        pntk->foreach_fanin( n, [&]( auto fi, auto ii ){
          auto tt0 = kitty::cofactor0( tt, ii );
          auto tt1 = kitty::cofactor1( tt, ii );
          auto diff = tt0 ^ tt1; // when 1 it matters

          typename resub_optimizer_t::signature_t mask;
          int num_minterms = ( 1u << tt.num_vars() );
          for( int m = 0; m < num_minterms; ++m )
          {
            tmp |= ~tmp;
            if( kitty::get_bit( diff, m ) == 1 ) // Boolean difference one -> you must check
            {
              for( int im{0}; im < tt.num_vars(); ++im )
              {
                if( im == ii )
                  continue;
                if( ( ( m >> im ) & 0x1 ) == 0x1 )
                {
                  tmp &= resub_opto.get_sig( children[im].index );
                }
                else
                {
                  tmp &= ~resub_opto.get_sig( children[im].index );
                }
              }
//              kitty::print_binary( tmp ); std::cout << " tmp " << m << std::endl;
              mask |= tmp;
            }
          }

          auto curr = power.get_dyn_switch( fi.index )*pntk->fanout_size( fi.index );
          pWin_manager->foreach_divisor( [&]( auto d, auto i ){
            if( pAnalyzer->get_arrival( d ) <= pAnalyzer->get_arrival( fi.index  ) )
            {
              auto cand = power.get_dyn_switch( d )*pntk->fanout_size( d );
              auto simi = resub_opto.get_sig( fi.index );
              auto simd = resub_opto.get_sig( d );
              if( (cand < curr) && ( d != fi.index ))
              {
                if( kitty::equal( simi&mask, simd&mask ) )
                {
                  candidates.push_back( std::make_tuple<signal, int, double>( pntk->make_signal( d ), ii, curr - cand ) );
                }
              }
            }
          } );
        } );

        if( candidates.size() > 0 )
        {
          double max_reward = 0;
          std::tuple<signal, int, double> chosen;
          for( auto cand : candidates )
          {
            //std::cout << std::get<0>( cand ) << " " << std::get<1>( cand ) << " " << std::get<2>( cand ) << " " << std::endl;
            if( std::get<2>( cand ) > max_reward )
            {
              max_reward = std::get<2>( cand );
              chosen = cand;
            }
          }

          children[std::get<1>( chosen )] = std::get<0>( chosen );
          auto g = pntk->get_binding( n );
          auto fnew = pntk->create_node( children, g.name )[0];
          pntk->substitute_node( n, fnew );
          st.num_resub++;
          update();

          //power.reset();
          //power( 1 );
          cnt++;
        }
      }
      return true;
    } );
  }

  void update()
  {
  }

private:
  std::shared_ptr<Ntk> pntk;
  std::shared_ptr<Database_t> pDatabase;
  std::shared_ptr<Analyzer_t> pAnalyzer;
  std::shared_ptr<window_manager<Ntk, max_leaves>> pWin_manager;
  resub_params ps;
  resub_stats & st;
 
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

template<class Ntk, class Database_t, uint32_t num_vars, uint32_t num_steps, uint32_t cube_size, uint32_t max_leaves>
void glitch_resub( Ntk& ntk, Database_t & database, resub_params ps = {} )
{
  using dNtk = depth_view<Ntk>;
  dNtk dntk{ ntk };
  using Analyzer_t = homo_xgx_analyzer<dNtk, num_vars, num_steps, cube_size, max_leaves>;
  resub_stats st;
  detail::resub_impl<dNtk, Database_t, Analyzer_t, num_vars, cube_size, max_leaves> p( dntk, database, ps, st );
  p.run();
  st.report();
}


} /* namespace mockturtle */