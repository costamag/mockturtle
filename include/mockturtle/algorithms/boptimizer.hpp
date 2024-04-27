/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
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
  \file boptimizer.hpp
  \brief Boolean optimizer

  \author Andrea Costamagna
*/

#pragma once

#include "../traits.hpp"
#include "../utils/progress_bar.hpp"
#include "../utils/stopwatch.hpp"
#include "../views/depth_view.hpp"
#include "../views/fanout_view.hpp"
#include "circuit_validator.hpp"
#include "pattern_generation.hpp"
#include "simulation.hpp"
#include "detail/resub_utils.hpp"
#include "cut_enumeration.hpp"
#include "cut_enumeration/rewrite_cut.hpp"

#include "dont_cares.hpp"
#include "reconv_cut.hpp"
#include "resyn_engines/lig_resyn.hpp"
#include "resyn_engines/scg_resyn.hpp"

#include <vector>

namespace mockturtle
{

std::mt19937 rng_opt(2);

//namespace rils
//{

/*! \brief Parameters for rboolean optimization.
 *
 * The data structure `boptimizer_params` holds configurable parameters with
 * default arguments for `boptimizer`.
 */
struct boptimizer_params
{
  boptimizer_params()
  {
    /* 0 < Cut limit < 16 */
    cut_enumeration_ps.cut_limit = 8;
    cut_enumeration_ps.minimize_truth_table = true;
  }
  /*! \brief Maximum number of PIs of reconvergence-driven cuts. */
  uint32_t max_pis{ 8 };

  /*! \brief Maximum number of divisors to consider. */
  uint32_t max_divisors{ 150 };

  /*! \brief Maximum number of nodes added by resubstitution. */
  double max_inserts{ 100 };

  /*! \brief Maximum fanout of a node to be considered as root. */
  uint32_t skip_fanout_limit_for_roots{ 1000 };

  /*! \brief Maximum fanout of a node to be considered as divisor. */
  uint32_t skip_fanout_limit_for_divisors{ 100 };

  /*! \brief Show progress. */
  bool progress{ false };

  /*! \brief Be verbose. */
  bool verbose{ false };

  bool verify_with_sim{false};

  bool timing_aware{false};

  bool use_wings{true};

  /****** window-based resub engine ******/

  /*! \brief Use don't cares for optimization. Only used by window-based resub engine. */
  bool use_dont_cares{ false };

  bool add_random_divs{false};

  /*! \brief Window size for don't cares calculation. Only used by window-based resub engine. */
  uint32_t window_size{ 12u };

  bool use_delay_constraints{false};
  bool high_effort_delay{false};

  /*! \brief Whether to prevent from increasing depth. Currently only used by window-based resub engine. */
  bool preserve_depth{ false };

  /****** simulation-based resub engine ******/

  /*! \brief Whether to use pre-generated patterns stored in a file.
   * If not, by default, 1024 random pattern + 1x stuck-at patterns will be generated. Only used by simulation-based resub engine.
   */
  std::optional<std::string> pattern_filename{};

  /*! \brief Whether to save the appended patterns (with CEXs) into file. Only used by simulation-based resub engine. */
  std::optional<std::string> save_patterns{};

  /*! \brief Maximum number of clauses of the SAT solver. Only used by simulation-based resub engine. */
  uint32_t max_clauses{ 1000 };

  /*! \brief Conflict limit for the SAT solver. Only used by simulation-based resub engine. */
  uint32_t conflict_limit{ 1000 };

  /*! \brief Random seed for the SAT solver (influences the randomness of counter-examples). Only used by simulation-based resub engine. */
  uint32_t random_seed{ 1 };

  /*! \brief Whether to utilize ODC, and how many levels. 0 = no. -1 = Consider TFO until PO. Only used by simulation-based resub engine. */
  int32_t odc_levels{ 0 };

  /*! \brief Maximum number of trials to call the resub functor. Only used by simulation-based resub engine. */
  uint32_t max_trials{ 100 };

  /* k-resub engine specific */
  /*! \brief Maximum number of divisors to consider in k-resub engine. Only used by `abc_resub_functor` with simulation-based resub engine. */
  uint32_t max_divisors_k{ 50 };

    cut_enumeration_params cut_enumeration_ps{};

};

/*! \brief Statistics for resubstitution.
 *
 * The data structure `boptimizer_stats` provides data collected by running
 * `resubstitution`.
 */
struct boptimizer_stats
{
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{ 0 };

  /*! \brief Accumulated runtime of the divisor collector. */
  stopwatch<>::duration time_divs{ 0 };

  /*! \brief Accumulated runtime of the divisor collector. */
  stopwatch<>::duration time_explore{ 0 };

  /*! \brief Accumulated runtime of the resub engine. */
  stopwatch<>::duration time_resub{ 0 };

  /*! \brief Accumulated runtime of the callback function. */
  stopwatch<>::duration time_callback{ 0 };

  /*! \brief Total number of divisors. */
  uint64_t num_total_divisors{ 0 };

  /*! \brief Total number of gain. */
  int estimated_gain{ 0 };

  /*! \brief Initial network size (before resubstitution). */
  uint64_t initial_size{ 0 };

  void report() const
  {
    // clang-format off
    fmt::print( "[i] <Top level>\n" );
    fmt::print( "[i]     ========  Stats  ========\n" );
    fmt::print( "[i]     #divisors = {:8d}\n", num_total_divisors );
    fmt::print( "[i]     est. gain = {:8d} ({:>5.2f}%)\n", estimated_gain, ( 100.0 * estimated_gain ) / initial_size );
    fmt::print( "[i]     ======== Runtime ========\n" );
    fmt::print( "[i]     total         : {:>5.2f} secs\n", to_seconds( time_total ) );
    fmt::print( "[i]       DivCollector: {:>5.2f} secs\n", to_seconds( time_divs ) );
    fmt::print( "[i]       Exploration : {:>5.2f} secs\n", to_seconds( time_explore ) );
    fmt::print( "[i]       ResubEngine : {:>5.2f} secs\n", to_seconds( time_resub ) );
    fmt::print( "[i]       callback    : {:>5.2f} secs\n", to_seconds( time_callback ) );
    fmt::print( "[i]     =========================\n\n" );
    // clang-format on
  }
};

namespace detail
{

template<typename Ntk>
bool substitute_fn( Ntk& ntk, typename Ntk::node const& n, typename Ntk::signal const& g )
{
  ntk.substitute_node( n, g );
  return true;
}

template<typename Ntk>
bool report_fn( Ntk& ntk, typename Ntk::node const& n, typename Ntk::signal const& g )
{
  fmt::print( "[i] Substitute node {} with signal {}{}\n", n, ntk.is_complemented( g ) ? "!" : "", ntk.get_node( g ) );
  return false;
}

#pragma region divisor_collection
struct collector_stats
{
  /*! \brief Total number of leaves. */
  uint64_t num_total_leaves{ 0 };

  /*! \brief Accumulated runtime for cut computation. */
  stopwatch<>::duration time_cuts{ 0 };

  /*! \brief Accumulated runtime for mffc computation. */
  stopwatch<>::duration time_mffc{ 0 };

  /*! \brief Accumulated runtime for divisor computation. */
  stopwatch<>::duration time_divs{ 0 };

  stopwatch<>::duration time_rand{ 0 };

  void report() const
  {
    // clang-format off
    fmt::print( "[i] <DivCollector: rils_divisor_collector>\n" );
    fmt::print( "[i]     #leaves = {:6d}\n", num_total_leaves );
    fmt::print( "[i]     ======== Runtime ========\n" );
    fmt::print( "[i]     reconv. cut : {:>5.2f} secs\n", to_seconds( time_cuts ) );
    fmt::print( "[i]     MFFC        : {:>5.2f} secs\n", to_seconds( time_mffc ) );
    fmt::print( "[i]     divs collect: {:>5.2f} secs\n", to_seconds( time_divs ) );
    fmt::print( "[i]     divs collect: {:>5.2f} secs\n", to_seconds( time_rand ) );
    fmt::print( "[i]     =========================\n\n" );
    // clang-format on
  }
};

/*! \brief Prepare the three public data members `leaves`, `divs` and `mffc`
 * to be ready for usage.
 *
 * `leaves`: sufficient support for all divisors
 * `divs`: divisor nodes that can be used for resubstitution
 * `mffc`: MFFC nodes which are needed to do simulation from
 * `leaves`, through `divs` and `mffc` until the root node,
 * but should be excluded from resubstitution.
 * The last element of `mffc` is always the root node.
 *
 * `divs` and `mffc` are in topological order.
 *
 * \param MffcMgr Manager class to compute the potential gain if a
 * resubstitution exists (number of MFFC nodes when the cost function is circuit size).
 * \param MffcRes Typename of the return value of `MffcMgr`.
 * \param cut_comp Manager class to compute reconvergence-driven cuts.
 */

template<class Ntk, class MffcMgr = rils_node_mffc_inside<Ntk>, typename MffcRes = double, typename cut_comp = detail::reconvergence_driven_cut_impl<Ntk>>
class rils_divisor_collector
{
public:
  using stats = collector_stats;
  using mffc_result_t = MffcRes;
  using node = typename Ntk::node;

  using cut_comp_parameters_type = typename cut_comp::parameters_type;
  using cut_comp_statistics_type = typename cut_comp::statistics_type;

public:
  explicit rils_divisor_collector( Ntk const& ntk, boptimizer_params const& ps, stats& st )
      : ntk( ntk ), ps( ps ), st( st ), cuts( ntk, cut_comp_parameters_type{ ps.max_pis }, cuts_st )
  {
  }

  bool run( node const& n, mffc_result_t& potential_gain )
  {
    /* skip nodes with many fanouts */
    if ( ntk.fanout_size( n ) > ps.skip_fanout_limit_for_roots )
    {
      return false;
    }

    /* compute a reconvergence-driven cut */
    leaves = call_with_stopwatch( st.time_cuts, [&]() {
      return cuts.run( { n } ).first;
    } );
    st.num_total_leaves += leaves.size();

    /* collect the MFFC */
    MffcMgr mffc_mgr( ntk );
    potential_gain = call_with_stopwatch( st.time_mffc, [&]() {
      return mffc_mgr.run( n, leaves, mffc );
    } );

    /* collect the divisor nodes in the cut */
    bool div_comp_success = call_with_stopwatch( st.time_divs, [&]() {
      return collect_divisors( n );
    } );

    if ( !div_comp_success )
    {
      return false;
    }

    return true;
  }

private:
  void collect_divisors_rec( node const& n )
  {
    /* skip visited nodes */
    if ( ntk.visited( n ) == ntk.trav_id() )
    {
      return;
    }
    ntk.set_visited( n, ntk.trav_id() );

    ntk.foreach_fanin( n, [&]( const auto& f ) {
      collect_divisors_rec( ntk.get_node( f ) );
    } );

    /* collect the internal nodes */
    if ( ntk.value( n ) == 0 && n != 0 ) /* ntk.fanout_size( n ) */
    {
      divs.emplace_back( n );
    }
  }

  bool collect_divisors( node const& root )
  {
    auto max_depth = std::numeric_limits<uint32_t>::max();
    if ( ps.preserve_depth )
    {
      max_depth = ntk.level( root ) - 1;
    }
    /* add the leaves of the cuts to the divisors */
    divs.clear();
    desp.clear();
    /*add the 0 divisor for constant resub*/
    ntk.incr_trav_id();
    for ( const auto& l : leaves )
    {
      divs.emplace_back( l );
      ntk.set_visited( l, ntk.trav_id() );
    }

    /* mark nodes in the MFFC */
    for ( const auto& t : mffc )
    {
      ntk.set_value( t, 1 );
    }

    /* collect the cone (without MFFC) */
    collect_divisors_rec( root );

    /* unmark the current MFFC */
    for ( const auto& t : mffc )
    {
      ntk.set_value( t, 0 );
    }

    /* check if the number of divisors is not exceeded */
    if ( divs.size() + mffc.size() - leaves.size() > ps.max_divisors - ps.max_pis )
    {
      return false;
    }
    uint32_t limit = ps.max_divisors - ps.max_pis - mffc.size() + leaves.size();

    if( ps.use_wings )
    {
      /* explore the fanouts, which are not in the MFFC */
      bool quit = false;
      for ( auto i = 0u; i < divs.size(); ++i )
      {
        auto const d = divs.at( i );

        if ( ntk.fanout_size( d ) > ps.skip_fanout_limit_for_divisors )
        {
          continue;
        }
        if ( divs.size() >= limit )
        {
          break;
        }

        /* if the fanout has all fanins in the set, add it */
        ntk.foreach_fanout( d, [&]( node const& p ) {
          if ( ntk.visited( p ) == ntk.trav_id() || ntk.level( p ) > max_depth )
          {
            return true; /* next fanout */
          }

          bool all_fanins_visited = true;
          ntk.foreach_fanin( p, [&]( const auto& g ) {
            if ( ntk.visited( ntk.get_node( g ) ) != ntk.trav_id() )
            {
              all_fanins_visited = false;
              return false; /* terminate fanin-loop */
            }
            return true; /* next fanin */
          } );

          if ( !all_fanins_visited )
            return true; /* next fanout */

          bool has_root_as_child = false;
          ntk.foreach_fanin( p, [&]( const auto& g ) {
            if ( ntk.get_node( g ) == root )
            {
              has_root_as_child = true;
              return false; /* terminate fanin-loop */
            }
            return true; /* next fanin */
          } );

          if ( has_root_as_child )
          {
            return true; /* next fanout */
          }

          divs.emplace_back( p );
          ntk.set_visited( p, ntk.trav_id() );

          /* quit computing divisors if there are too many of them */
          if ( divs.size() >= limit )
          {
            quit = true;
            return false; /* terminate fanout-loop */
          }

          return true; /* next fanout */
        } );

        if ( quit )
        {
          break;
        }
      }
    }

    /* note: different from the previous version, now we do not add MFFC nodes into divs */
    assert( root == mffc.at( mffc.size() - 1u ) );
    /* note: this assertion makes sure window_simulator does not go out of bounds */
    assert( divs.size() + mffc.size() - leaves.size() <= ps.max_divisors - ps.max_pis );

    for( auto nd : mffc )
    {
      if( nd == root ) continue;
      bool is_extr=true;
      ntk.foreach_fanin( nd, [&]( const auto& g ) {
        auto ng = ntk.get_node(g);
        if ( std::find( leaves.begin(), leaves.end(), ng ) == leaves.end() )
        {
          is_extr = false;
        }
      } );
      if( is_extr )
      {
        ntk.foreach_fanin( root, [&]( const auto& g ) {
        auto ng = ntk.get_node(g);
        if( ng == nd ) 
          is_extr=false;
      } );
      }
      if( is_extr )
        desp.push_back( nd );
    }

    return true;
  }

private:
  Ntk const& ntk;
  boptimizer_params ps;
  stats& st;

  cut_comp cuts;
  cut_comp_statistics_type cuts_st;

public:
  std::vector<node> leaves;
  std::vector<node> divs;
  std::vector<node> mffc;
  std::vector<node> desp;
};

#pragma endregion divisor_collection

#pragma region window_boptimizer

template<typename ResubFnSt>
struct window_boptimizer_stats
{
  /*! \brief Number of successful resubstitutions. */
  uint32_t num_resub{ 0 };

  /*! \brief Time for simulation. */
  stopwatch<>::duration time_sim{ 0 };

  /*! \brief Time for don't-care computation. */
  stopwatch<>::duration time_dont_care{ 0 };

  /*! \brief Time of the resub functor. */
  stopwatch<>::duration time_compute_function{ 0 };

  /*! \brief Time for pattern generation. */
  stopwatch<>::duration time_patgen{ 0 };

  /*! \brief Time for saving patterns. */
  stopwatch<>::duration time_patsave{ 0 };

  /*! \brief Time for simulation. */
  stopwatch<>::duration time_lsim{ 0 };

  /*! \brief Time for SAT solving. */
  stopwatch<>::duration time_sat{ 0 };
  stopwatch<>::duration time_sat_restart{ 0 };

  /*! \brief Time for computing ODCs. */
  stopwatch<>::duration time_odc{ 0 };

  /*! \brief Time for finding dependency function. */
  stopwatch<>::duration time_resyn{ 0 };

  /*! \brief Time for translating from index lists to network signals. */
  stopwatch<>::duration time_interface{ 0 };

  /*! \brief Number of patterns used. */
  uint32_t num_pats{ 0 };

  /*! \brief Number of counter-examples. */
  uint32_t num_cex{ 0 };

  /*! \brief Number of successful resubstitutions. */
  //uint32_t num_resub{ 0 };

  /*! \brief Number of SAT solver timeout. */
  uint32_t num_timeout{ 0 };

  /*! \brief Number of calls to the resynthesis engine. */
  uint32_t num_resyn{ 0 };

  void report() const
  {
    fmt::print( "[i] <ResubEngine: simulation_based_resub_engine>\n" );
    fmt::print( "[i]     #resub = {:6d}\n", num_resub );
    fmt::print( "[i]     ========  Stats  ========\n" );
    fmt::print( "[i]     #pat        = {:6d}\n", num_pats );
    fmt::print( "[i]     #resyn call = {:6d}\n", num_resyn );
    fmt::print( "[i]     #valid      = {:6d}\n", num_resub );
    fmt::print( "[i]     #CEX        = {:6d}\n", num_cex );
    fmt::print( "[i]     #timeout    = {:6d}\n", num_timeout );
    fmt::print( "[i]     ======== Runtime ========\n" );
    fmt::print( "[i]     generate pattern: {:>5.2f} secs [excluded]\n", to_seconds( time_patgen ) );
    fmt::print( "[i]     save pattern    : {:>5.2f} secs [excluded]\n", to_seconds( time_patsave ) );
    fmt::print( "[i]     g-simulation    : {:>5.2f} secs\n", to_seconds( time_sim ) );
    fmt::print( "[i]     l-simulation    : {:>5.2f} secs\n", to_seconds( time_sim ) );
    fmt::print( "[i]     don't care      : {:>5.2f} secs\n", to_seconds( time_dont_care ) );
    fmt::print( "[i]     functor         : {:>5.2f} secs\n", to_seconds( time_compute_function ) );
    fmt::print( "[i]     SAT solve       : {:>5.2f} secs\n", to_seconds( time_sat ) );
    fmt::print( "[i]     SAT restart     : {:>5.2f} secs\n", to_seconds( time_sat_restart ) );
    fmt::print( "[i]     compute ODCs    : {:>5.2f} secs\n", to_seconds( time_odc ) );
    fmt::print( "[i]     interfacing     : {:>5.2f} secs\n", to_seconds( time_interface ) );
    fmt::print( "[i]     compute function: {:>5.2f} secs\n", to_seconds( time_resyn ) );
    fmt::print( "[i]     ======== Details ========\n" );
    resyn_st.report();
    fmt::print( "[i]     =========================\n\n" );
  }

  ResubFnSt resyn_st;

};

template<class Ntk, typename validator_t, class ResynEngine, uint32_t SizeSupp, uint32_t nPisLoc=4u, uint32_t nPisGlb=4u, typename MffcRes = double>
class window_boptimizer
{
public:
  static constexpr bool require_leaves_and_mffc = true;
  using TTsig = kitty::static_truth_table<nPisGlb>;
  using TTcut = kitty::static_truth_table<nPisLoc>;
  using TTtmp = kitty::static_truth_table<6u>;
  using stats = window_boptimizer_stats<typename ResynEngine::stats>;
  using mffc_result_t = MffcRes;


  static constexpr uint32_t num_vars = SizeSupp;
  static constexpr uint32_t max_window_size = 8u;
  using network_cuts_t = dynamic_network_cuts<Ntk, num_vars, true, cut_enumeration_rewrite_cut>;
  using cut_manager_t = detail::dynamic_cut_enumeration_impl<Ntk, num_vars, true, cut_enumeration_rewrite_cut>;
  using cut_t = typename network_cuts_t::cut_t;
  using node_data = typename Ntk::storage::element_type::node_type;


  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

  explicit window_boptimizer( Ntk& ntk, boptimizer_params const& ps, stats& st )
      : ntk( ntk ), ps( ps ), st( st ), tts( ntk ), tt6( ntk ), _arr_times(ntk), _req_times(ntk), _lSim( ntk, ps.max_divisors, nPisLoc ), engine( ntk._library, st.resyn_st ), validator( ntk, { ps.max_clauses, ps.odc_levels, ps.conflict_limit, ps.random_seed } )//,
        //_cuts( ntk.size() + ( ntk.size() >> 1 ) ),
        //_cut_manager( ntk, ps.cut_enumeration_ps, _cst, _cuts )
  {
    add_event = ntk.events().register_add_event( [&]( const auto& n ) {
      tts.resize();
      tt6.resize();
      _arr_times.resize();
      _req_times.resize();
      call_with_stopwatch( st.time_sim, [&]() {
        simulate_node_static<Ntk, nPisGlb>( ntk, n, tts, _gSim );
        simulate_node_static<Ntk, 6>( ntk, n, tt6, _6Sim );
      } );
    } );
    //v1 asap_timing_information();
    kitty::static_truth_table<nPisGlb> tts0;
    kitty::static_truth_table<nPisLoc> tt60;
    tts[0] = tts0;
    tt6[0] = tt60;
    if( ps.use_delay_constraints )
    {
      timing_information();
    }

  }

  explicit window_boptimizer( Ntk& ntk, boptimizer_params const& ps, stats& st, std::vector<gate> const& gates )
      : ntk( ntk ), ps( ps ), st( st ), tts( ntk ), tt6( ntk ), _lSim( ntk, ps.max_divisors, ps.max_pis ), engine( st.resyn_st, gates ), validator( ntk, { ps.max_clauses, ps.odc_levels, ps.conflict_limit, ps.random_seed } )
      //, _cuts( ntk.size() + ( ntk.size() >> 1 ) ),
      //  _cut_manager( ntk, ps.cut_enumeration_ps, _cst, _cuts )
  {
    add_event = ntk.events().register_add_event( [&]( const auto& n ) {
      tts.resize();
      tt6.resize();
      _arr_times.resize();
      _req_times.resize();
      call_with_stopwatch( st.time_sim, [&]() {
        simulate_node_static<Ntk, nPisGlb>( ntk, n, tts, _gSim );
        simulate_node_static<Ntk, 6>( ntk, n, tt6, _6Sim );
      } );
    } );

    kitty::static_truth_table<nPisGlb> tts0;
    kitty::static_truth_table<nPisLoc> tt60;
    tts[0] = tts0;
    tt6[0] = tt60;
    if( ps.use_delay_constraints )
    {
      timing_information();
    }

  }

  ~window_boptimizer()
  {
//    if ( ps.save_patterns )
//    {
//      call_with_stopwatch( st.time_patsave, [&]() {
//        write_patterns( _gSim, *ps.save_patterns );
//      } );
//    }

    if ( add_event )
    {
      ntk.events().release_add_event( add_event );
    }
  }

  void init() 
  {
    /* prepare simulation patterns */
    call_with_stopwatch( st.time_patgen, [&]() {
      _gSim = static_simulator<nPisGlb>( ntk.num_pis() );
      _6Sim = static_simulator<6>( ntk.num_pis() );
    } );
    st.num_pats = _gSim.num_bits();
    assert( _gSim.num_bits() > 0 );
    assert( _6Sim.num_bits() > 0 );

    /* first simulation: the whole circuit; from 0 bits. */
    call_with_stopwatch( st.time_sim, [&]() {
      simulate_nodes_static<Ntk, nPisGlb>( ntk, tts, _gSim, true );
      simulate_nodes_static<Ntk, 6>( ntk, tt6, _6Sim, true );
    } );


    /* initialize cuts for constant nodes and PIs */
    //_cut_manager.init_cuts();

  }

  void init_topo_order()
  {
    topo_order.reserve( ntk.size() );

    topo_view<Ntk>( ntk ).foreach_node( [this]( auto n ) {
      topo_order.push_back( n );
    } );
  }

  double prop_arr_rec( node n )
  {
    if( _arr_times.has(n) )
    {
      return _arr_times[n];
    }

    if ( ntk.has_binding( n ) )
    {
      auto const& g = ntk.get_binding( n );
      double gate_delay = 0;
      ntk.foreach_fanin( n, [&]( auto const& f, auto i ) {
        double arr_fanin = prop_arr_rec( ntk.get_node(f) );
        gate_delay = std::max( gate_delay, (double)( arr_fanin + std::max( g.pins[i].rise_block_delay, g.pins[i].fall_block_delay ) ) );
      } );
      _arr_times[n] = gate_delay;
    }
    else
    {
      double gate_delay = 1;
      ntk.foreach_fanin( n, [&]( auto const& f, auto i ) {
        double arr_fanin = prop_arr_rec( ntk.get_node(f) );
        gate_delay = std::max( gate_delay, (double)( arr_fanin + 1u ) );
      } );
      _arr_times[n] = gate_delay;
    }
    return _arr_times[n];
  }

  double propagate_arrival_times()
  {
    _arr_times.reset();
    ntk.foreach_pi( [&]( auto const& n, auto i ) {
      _arr_times[n]=0.0;
    } );
    _arr_times[0]=0.0;

    double max_delay=0;
    ntk.foreach_po( [&]( auto const& no, auto i ) {
      auto out_del = prop_arr_rec(ntk.get_node(no));
      if( out_del > max_delay )
      {
        max_delay = out_del;
      }
    } );
    return max_delay;
  }

//  double propagate_arrival_times()
//  {
//    ntk.foreach_pi( [&]( auto const& n, auto i ) {
//      _arr_times[n]=0.0;
//    } );
//    _arr_times[0]=0.0;
//
//    double max_delay;
//    for ( auto it = topo_order.begin(); it != topo_order.end(); ++it )
//    {
//      if ( ntk.is_constant( *it ) || ntk.is_pi( *it ) )
//      {
//        _arr_times[*it]=0.0;
//        continue;
//      }
//      if ( ntk.has_binding( *it ) )
//      {
//        auto const& g = ntk.get_binding( *it );
//        double gate_delay = 0;
//        ntk.foreach_fanin( *it, [&]( auto const& f, auto i ) {
//          gate_delay = std::max( gate_delay, (double)( _arr_times[f] + std::max( g.pins[i].rise_block_delay, g.pins[i].fall_block_delay ) ) );
//        } );
//        _arr_times[*it] = gate_delay;
//        if( gate_delay > max_delay )
//        {
//          max_delay = gate_delay;
//        }
//      }
//      else
//      {
//        double gate_delay = 1;
//        ntk.foreach_fanin( *it, [&]( auto const& f, auto i ) {
//          gate_delay = std::max( gate_delay, (double)( _arr_times[f] + 1u ) );
//        } );
//        _arr_times[*it] = gate_delay;
//
//        if( gate_delay > max_delay )
//        {
//          max_delay = gate_delay;
//        }
//      }
//    }
//    return max_delay;
//  }

  void propagate_required_times( double worst_delay )
  {
    init_topo_order();
    for ( auto it = topo_order.begin(); it != topo_order.end(); ++it )
    {
      _req_times[*it]=worst_delay+1;
    }

    ntk.foreach_pi( [&]( auto const& n, auto i ) {
      _req_times[n]=0.0;
    } );

    ntk.foreach_po( [&]( auto const& n, auto i ) {
      _req_times[n]=worst_delay;
    } );


    for ( auto it = topo_order.rbegin(); it != topo_order.rend(); ++it )
    {
      if( ntk.is_pi(*it) )
      {
        continue;
      }
      auto const& fos = ntk.fanout( *it );
      if( ntk.has_binding( *it ) )
      {
        for( auto fo : fos )
        {
          auto const& g = ntk.get_binding( fo );
          uint32_t idx;
          ntk.foreach_fanin( fo, [&]( auto const& f, auto i ) {
            if( *it == ntk.get_node(f) )
            {
              idx=i;
              return;
            }
          } );
          _req_times[*it] = std::min( _req_times[*it], (double)( _req_times[fo] - std::max( g.pins[idx].rise_block_delay, g.pins[idx].fall_block_delay ) ) );
        }
      }
      else
      {
        for( auto fo : fos )
        {
          uint32_t idx;
          ntk.foreach_fanin( fo, [&]( auto const& f, auto i ) {
            if( *it == ntk.get_node(f) )
            {
              idx=i;
              return;
            }
          } );
          _req_times[*it] = std::min( _req_times[*it], (double)( _req_times[fo] - 1 ) );
        }
      }
    }

  }

  void timing_information()
  {
    _W_REQ_NODES.clear();
    _arr_times.reset();
    _req_times.reset();

    //compute_required_time();
    double worst_delay = propagate_arrival_times();
    propagate_required_times( worst_delay );

  }

  void update( node n, node nn ) 
  {
    //ntk.clear_marks();
    //timing_information();
    //init_topo_order();
    //propagate_arrival_times( );
              //reset_arrival_from( n );
              //_arr_times.reset();
              //compute_required_time();
    //if( _DO_ARR )
    //  timing_information();
    //else
    //if( ps.use_delay_constraints )
    //]{
    //  if( _DO_ARR )
    //timing_information();
    //  else
    //}
    if( ps.use_delay_constraints )
    {
      //timing_information();
      propagate_arrival_times();
      
    }

    if constexpr ( validator_t::use_odc_ || has_EXODC_interface_v<Ntk> )
    {
      call_with_stopwatch( st.time_sat_restart, [&]() {
        validator.update();
      } );
      tts.reset();
      call_with_stopwatch( st.time_sim, [&]() {
        simulate_nodes_static<Ntk, nPisGlb>( ntk, tts, _gSim, true );
      } );
    }
  }

  template<class LIST, class LIB>
  double compute_worst_delay( LIST list, std::vector<double> divs_delays, LIB const& lib )
  {
    //std::cout << to_index_list_string(list) << std::endl;
    double delay;
    if constexpr( std::is_same_v<typename Ntk::base_type, scopt::scg_network> )
    {
      list.foreach_gate( [&]( std::vector<uint32_t> children, uint32_t func_lit ) {
        auto g = lib[list.ids[func_lit]];
        int i{0};
        double delay=0;
        for( uint32_t child : children )
        {
          delay = std::max( delay, ( double)( divs_delays[(child >> 1u)] + std::max( g.pins[i].rise_block_delay, g.pins[i].fall_block_delay ) ) );
          i++;
          //printf("%d del %f->%f\n", func_lit, divs_delays[(child >> 1u)], delay );
        }
        divs_delays.push_back( delay );
      } );
    }
    else
    {
      list.foreach_gate( [&]( std::vector<uint32_t> children, uint32_t func_lit ) {
        int i{0};
        double delay=0;
        for( uint32_t child : children )
        {
          delay = std::max( delay, ( double)( divs_delays[(child >> 1u)] + 1 ));
          i++;
          //printf("%d del %f->%f\n", func_lit, divs_delays[(child >> 1u)], delay );
        }
        divs_delays.push_back( delay );
      } );
    }
//    if( list.num_gates() == 0 )
//    {
//      //std::cout << to_index_list_string(list) << std::endl;
//      printf( "return %f\n",divs_delays[list.values.back()>>1] );
//      std::vector<int> verr;
//      std::cout << verr[-1];
//    }

    return divs_delays[list.values.back()>>1];
  }

  void recursively_mark( node n )
  {
    if( ntk.is_pi(n) || ntk.is_constant(n) || ntk.is_marked(n) )
    {
      return;
    }

    ntk.foreach_fanin( n, [&]( const auto& f ) {
      recursively_mark( ntk.get_node(f) );
    } );

    ntk.set_mark( n );
  }

  std::optional<signal> run( node const& n, std::vector<node> const& leaves, std::vector<node> const& divs,  std::vector<node> const& desps, std::vector<node> const& mffc, mffc_result_t potential_gain, double& last_gain )
  {

    if( ps.use_delay_constraints && ntk.is_marked(n) )
    {
      timing_information();
      ntk.clear_marked();
    }

    /* make valid the simulation at each divisor node */
    check_tts( n );
    for ( auto const& d : divs )
    {
      check_tts( d );
    }

    /* compute the observability don't cares */
    kitty::static_truth_table<nPisGlb> const care = _gSim.compute_constant( true );//~observability_dont_cares( ntk, n, _gSim, tts, 10 );

    std::vector<double> divs_delays, divs_delays2;
    divs_delays.push_back(0);
    divs_delays2.push_back(0);
    _arr_times[0]=0;
    int i=0;
    //printf("%d|%d|%f\n", i++, 0, _arr_times[0] );
    for( auto div : divs )
    {
      divs_delays.push_back( _arr_times[div] );
      divs_delays2.push_back( _arr_times[div] );
     // printf("%d %f\n", div, _arr_times[div] );
    }
    for( auto div : desps )
    {
      divs_delays2.push_back( _arr_times[div] );
    }
    //printf("\n");

    //for( int i{0}; i<divs_delays.size(); ++i )
    //{
    //  printf("%d %d %f\n", i, divs[i], divs_delays[i] );
    //}
    //printf("\n");

/////////////////////////////////////////////////////////////



  //  _cut_manager.clear_cuts( n );
  //  _cut_manager.compute_cuts( n );
//
  //  uint32_t cut_index = 0;
  //  for ( auto& cut : _cuts.cuts( ntk.node_to_index( n ) ) )
  //  {
  //    /* skip trivial cut */
  //    if ( ( cut->size() == 1 && *cut->begin() == ntk.node_to_index( n ) ) )
  //    {
  //      ++cut_index;
  //      continue;
  //    }
  //    printf("we have a structural cut\n");
  //    /* Boolean matching */
  //    kitty::print_binary( _cuts.truth_table( *cut ) ); //printf("\n");
  //    /* measure the MFFC contained in the cut */
  //    double mffc_size = measure_mffc_deref( n, cut );
  //    double gain = mffc_size;// - nodes_added;
  //    printf("->%f\n", gain );
  //    const auto res_struct = call_with_stopwatch( st.time_resyn, [&]() {
  //      ++st.num_resyn;
  //      return engine( _cuts. _cuts.truth_table( *cut ), std::min( struct_gain, ps.max_inserts ) );
  //    } );
  //        measure_mffc_ref( n, cut );
  //  }

    for ( auto j = 0u; j < ps.max_trials; ++j )
    {
      /* do resynthesis */
      const auto res = call_with_stopwatch( st.time_resyn, [&]() {
        ++st.num_resyn;
        return engine( tts[n], care, std::begin( divs ), std::end( divs ), tts, std::min( potential_gain, ps.max_inserts ), j );
      } );
      if( res )
      {
        auto const& id_list = *res;
        assert( id_list.num_pos() == 1u );
        last_gain = potential_gain - id_list.get_area();

        double delay_candidate = ps.use_delay_constraints ? compute_worst_delay( id_list, divs_delays, ntk.get_library() ) : 0;
        if( !ps.use_delay_constraints || delay_candidate < _req_times[n] )
        {
          auto valid = call_with_stopwatch( st.time_sat, [&]() {
            return validator.validate( n, divs, id_list );
          } );
          if ( valid )
          {
            if ( *valid )
            {
              _stats_gen1[id_list.num_gates()]++;
              _stats_genT[id_list.num_gates()]++;
              ++st.num_resub;

              signal out_sig;
              std::vector<signal> divs_sig( divs.size() );
              std::vector<node> upd_req_nodes;
              call_with_stopwatch( st.time_interface, [&]() {
                std::transform( divs.begin(), divs.end(), divs_sig.begin(), [&]( const node n ) {
                  return ntk.make_signal( n );
                } );

                  insert( ntk, divs_sig.begin(), divs_sig.end(), id_list, [&]( signal const& s ) {
                  out_sig = s;
                  _NNEW = ntk.get_node(out_sig);
                } );
              } );
              
              if( ps.use_delay_constraints )
              {
                recursively_mark( ntk.get_node(out_sig) );
              }
              _DELAY_NEW = delay_candidate;
              return out_sig;
            }
            else
            {
              _stats_genT[id_list.num_gates()]++;
              found_cex();
              continue;
            }
          }
        }
        else
        {
          continue;
        }

      }
      else /* functor can not find any potential resubstitution */
      {
        return std::nullopt;
      }
    }

    bool try_desp = false;
    if( try_desp )
    {
      auto divs2 = divs;
      uint32_t nZero = divs.size();
      for( auto desp : desps )
      {
        divs2.push_back( desp );
      }

      for ( auto j = 0u; j < ps.max_trials; ++j )
      {
        /* do resynthesis */
        const auto res = call_with_stopwatch( st.time_resyn, [&]() {
          ++st.num_resyn;
          return engine( tts[n], care, std::begin( divs2 ), std::end( divs2 ), nZero, tts, std::min( potential_gain, ps.max_inserts ) );
        } );
        if( res )
        {
          auto const& id_list = *res;
          assert( id_list.num_pos() == 1u );
          last_gain = potential_gain - id_list.get_area();

          double delay_candidate = ps.use_delay_constraints ? compute_worst_delay( id_list, divs_delays2, ntk.get_library() ) : 0;
          if( !ps.use_delay_constraints || delay_candidate < _req_times[n] )
          {
            auto valid = call_with_stopwatch( st.time_sat, [&]() {
              return validator.validate( n, divs2, id_list );
            } );
            if ( valid )
            {
              if ( *valid )
              {
                _stats_gen1[id_list.num_gates()]++;
                _stats_genT[id_list.num_gates()]++;
                ++st.num_resub;

                signal out_sig;
                std::vector<signal> divs_sig( divs2.size() );
                std::vector<node> upd_req_nodes;
                call_with_stopwatch( st.time_interface, [&]() {
                  std::transform( divs2.begin(), divs2.end(), divs_sig.begin(), [&]( const node n ) {
                    return ntk.make_signal( n );
                  } );

                    insert( ntk, divs_sig.begin(), divs_sig.end(), id_list, [&]( signal const& s ) {
                    out_sig = s;
                    _NNEW = ntk.get_node(out_sig);
                  } );
                } );
                
                if( ps.use_delay_constraints )
                {
                  recursively_mark( ntk.get_node(out_sig) );
                }
                _DELAY_NEW = delay_candidate;
                return out_sig;
              }
              else
              {
                _stats_genT[id_list.num_gates()]++;
                found_cex();
                continue;
              }
            }
          }
          else
          {
            continue;
          }

        }
        else /* functor can not find any potential resubstitution */
        {
          return std::nullopt;
        }
      }
    }

    return std::nullopt;
  }


  TTcut simulate_subnet( typename Ntk::signal sig, incomplete_node_map<TTcut, Ntk>& loc_map )
  {
    auto nd = ntk.get_node(sig);
    if( loc_map.has( nd ) )
    {
      return ntk.is_complemented( sig ) ? ~loc_map[nd] : loc_map[nd];
    }

    std::vector<TTcut> tti;
    ntk.foreach_fanin( nd, [&]( const auto& f ) {
      tti.push_back( simulate_subnet( f, loc_map ) );
    } );

    loc_map[nd] = ntk.compute( nd, tti );
    return loc_map[nd];
  }

  void found_cex()
  {
    sig_pointer = (sig_pointer+1)%(1<<nPisGlb);
    ++st.num_cex;

    _6Sim.add_pattern( validator.cex );
    if ( sig_pointer % 64 == 0 )
    {
      tt6.reset();
      call_with_stopwatch( st.time_sim, [&]() {
        simulate_nodes_static<Ntk>( ntk, tt6, _6Sim, true );
      } );
      ntk.foreach_pi( [&]( auto const& n, auto i ) {
        *(tts[n].begin() + _block) = *(tt6[n].begin());
      } );

      ntk.foreach_gate( [&]( auto const& n, auto i ) {
        *(tts[n].begin() + _block) = *(tt6[n].begin());
      } );

      _block = nPisGlb == 6u ? 0u : ( _block + 1u ) % ( ( 1u << ( nPisGlb - 6u ) ) - 1u ) ;
    }
  }

private:
  void simulate( std::vector<node> const& leaves, std::vector<node> const& divs, std::vector<node> const& mffc )
  {
    _lSim.resize();
    for ( auto i = 0u; i < divs.size() + mffc.size(); ++i )
    {
      const auto d = i < divs.size() ? divs.at( i ) : mffc.at( i - divs.size() );

      /* skip constant 0 */
      if ( d == 0 )
        continue;

      /* assign leaves to variables */
      if ( i < leaves.size() )
      {
        _lSim.assign( d, i + 1 );
        continue;
      }

      /* compute truth tables of inner nodes */
      _lSim.assign( d, i - uint32_t( leaves.size() ) + ps.max_pis + 1 );
      std::vector<TTcut> tts;
      ntk.foreach_fanin( d, [&]( const auto& s ) {
        tts.emplace_back( _lSim.get_tt( ntk.make_signal( ntk.get_node( s ) ) ) ); /* ignore sign */
      } );

      auto const tt = ntk.compute( d, tts.begin(), tts.end() );
      _lSim.set_tt( i - uint32_t( leaves.size() ) + ps.max_pis + 1, tt );
    }

    /* normalize truth tables */
    _lSim.normalize( divs );
    _lSim.normalize( mffc );
  }

  void check_tts( node const& n )
  {
    if ( tts[n].num_bits() != _gSim.num_bits() )
    {
      call_with_stopwatch( st.time_sim, [&]() {
        simulate_node_static<Ntk, nPisGlb>( ntk, n, tts, _gSim );
      } );
    }
  }

  double measure_mffc_ref( node const& n, cut_t const* cut )
  {
    /* reference cut leaves */
    for ( auto leaf : *cut )
    {
      ntk.incr_fanout_size( ntk.index_to_node( leaf ) );
    }

    double mffc_size = static_cast<double>( recursive_ref( n ) );

    /* dereference leaves */
    for ( auto leaf : *cut )
    {
      ntk.decr_fanout_size( ntk.index_to_node( leaf ) );
    }

    return mffc_size;
  }

  double measure_mffc_deref( node const& n, cut_t const* cut )
  {
    /* reference cut leaves */
    for ( auto leaf : *cut )
    {
      ntk.incr_fanout_size( ntk.index_to_node( leaf ) );
    }

    double mffc_size = static_cast<double>( recursive_deref( n ) );

    /* dereference leaves */
    for ( auto leaf : *cut )
    {
      ntk.decr_fanout_size( ntk.index_to_node( leaf ) );
    }

    return mffc_size;
  }

  double recursive_deref( node const& n )
  {
    /* terminate? */
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      return 0;

    /* recursively collect nodes */
    double value{ ntk.get_area( n ) };
    ntk.foreach_fanin( n, [&]( auto const& s ) {
      if ( ntk.decr_fanout_size( ntk.get_node( s ) ) == 0 )
      {
        value += recursive_deref( ntk.get_node( s ) );
      }
    } );
    return value;
  }

  double recursive_ref( node const& n )
  {
    /* terminate? */
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      return 0;

    /* recursively collect nodes */
    double value{ ntk.get_area( n ) };
    ntk.foreach_fanin( n, [&]( auto const& s ) {
      if ( ntk.incr_fanout_size( ntk.get_node( s ) ) == 0 )
      {
        value += recursive_ref( ntk.get_node( s ) );
      }
    } );
    return value;
  }

public:
  TTsig get_signature( node n )
  {
    return tts[n];
  }

public:
  Ntk& ntk;
  boptimizer_params const& ps;
  stats& st;
  uint32_t _block{0};
  incomplete_node_map<TTsig, Ntk> tts;
  incomplete_node_map<TTtmp, Ntk> tt6;
  incomplete_node_map<double, Ntk> _arr_times;
  incomplete_node_map<double, Ntk> _req_times;
  std::set<signal> _W_REQ_NODES;
  std::set<signal> _W_ARR_NODES;

  std::vector<node> topo_order;
  double _DELAY_NEW{0};
  node _NNEW;
  bool _DO_ARR{false};
  double dT{0};
  window_simulator<Ntk, TTcut> _lSim;
  static_simulator<nPisGlb> _gSim;
  static_simulator<6u> _6Sim;
  uint32_t sig_pointer{0};
  std::default_random_engine::result_type _seed=1;
  validator_t validator;
  ResynEngine engine;

  /* events */
  std::shared_ptr<typename network_events<Ntk>::add_event_type> add_event;

  std::array<double,10> _stats_gen1{0};
  std::array<double,10> _stats_genT{0};


    /* initialize the cut manager */
    //cut_enumeration_stats _cst;
    //network_cuts_t _cuts;
    //cut_manager_t _cut_manager;

}; /* simulation_based_resub_engine_for_lig */


#pragma endregion window_boptimizer

#pragma region boptimizer_impl
template<class Ntk, class ResubEngine>
class boptimizer_impl
{
public:
  using DivCollector = detail::rils_divisor_collector<Ntk>;
  using engine_st_t = typename ResubEngine::stats;
  using collector_st_t = typename DivCollector::stats;
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;
  using resub_callback_t = std::function<bool( Ntk&, node const&, signal const& )>;
  using mffc_result_t = typename ResubEngine::mffc_result_t;

  /*! \brief Constructor of the top-level boptimizer framework.
   *
   * \param ntk The network to be optimized.
   * \param ps boptimizer parameters.
   * \param st Top-level boptimizer statistics.
   * \param engine_st Statistics of the boptimizer engine.
   * \param collector_st Statistics of the divisor collector.
   * \param callback Callback function when a boptimizer is found.
   */
  explicit boptimizer_impl( Ntk& ntk, boptimizer_params const& ps, boptimizer_stats& st, engine_st_t& engine_st, collector_st_t& collector_st )
      : ntk( ntk ), ps( ps ), st( st ), engine_st( engine_st ), collector_st( collector_st )
  {
    static_assert( std::is_same_v<typename ResubEngine::mffc_result_t, typename DivCollector::mffc_result_t>, "MFFC result type of the engine and the collector are different" );

    st.initial_size = ntk.num_gates();

    register_events();
  }

  ~boptimizer_impl()
  {
    ntk.events().release_add_event( add_event );
    ntk.events().release_modified_event( modified_event );
    ntk.events().release_delete_event( delete_event );
  }

  void run( resub_callback_t const& callback = substitute_fn<Ntk> )
  {
    stopwatch t( st.time_total );

    /* start the managers */
    DivCollector collector( ntk, ps, collector_st );
    ResubEngine resub_engine( ntk, ps, engine_st );
    call_with_stopwatch( st.time_resub, [&]() {
      resub_engine.init( );
    } );

    progress_bar pbar{ ntk.size(), "resub |{0}| node = {1:>4}   cand = {2:>4}   est. gain = {3:>5}", ps.progress };

    auto const size = ntk.num_gates();
    ntk.foreach_gate( [&]( auto const& n, auto i ) {
      if ( i >= size )
      {
        return false; /* terminate */
      }
      if ( (ntk.fanin_size(n) == 1) && ntk.po_index(n) != -1 )
      {
        return true; /* terminate */
      }

      if ( ntk.is_constant(n) )
      {
        return true; /* terminate */
      }

      if ( (ntk.fanin_size(n) == 1) && (ntk.is_pi( ntk.get_children( n, 0 ) ) || ntk.is_constant( ntk.get_children( n, 0 ) )) )
      {
        return true; /* terminate */
      }

      //kitty::print_binary( resub_engine.get_signature( n )); printf("@\n");

      pbar( i, i, candidates, st.estimated_gain );

      /* compute cut, collect divisors, compute MFFC */
      mffc_result_t potential_gain;
      const auto collector_success = call_with_stopwatch( st.time_divs, [&]() {
        return collector.run( n, potential_gain );
      } );
      if ( !collector_success )
      {
        return true; /* next */
      }

      /* update statistics */
      last_gain = 0;
      st.num_total_divisors += collector.divs.size();

      /* try to find a boptimizer with the divisors */
      auto g = call_with_stopwatch( st.time_resub, [&]() {
          return resub_engine.run( n, collector.leaves, collector.divs, collector.desp, collector.mffc, potential_gain, last_gain );
      } );
      if ( !g )
      {
        return true; /* next */
      }

      /* update progress bar */
      candidates++;
      st.estimated_gain += last_gain;

      /* update network */
      //double dbef = ntk.compute_worst_delay();

      bool updated = call_with_stopwatch( st.time_callback, [&]() {
        return callback( ntk, n, *g );
      } );

      if ( updated )
      {
        resub_engine.update( n, ntk.get_node(*g) );
      }
      //double daft = ntk.compute_worst_delay();
    //  if( daft > dbef )
    //  {
    //    std::vector<int> ver;
    //    std::cout << daft <<">"<< dbef << std::endl;
    //    std::cout << ver[-1] << std::endl;
    //  }

      return true; /* next */
    } );

//    for( int i{0}; i<resub_engine._stats_genT.size(); ++i )
//    {
//      if( resub_engine._stats_genT[i] > 0 )
//      {
//        printf("%d %f (%f)\n", i, resub_engine._stats_gen1[i]/resub_engine._stats_genT[i], resub_engine._stats_genT[i]);
//      }
//    }
  }

private:
  void register_events()
  {
    auto const update_level_of_new_node = [&]( const auto& n ) {
      ntk.resize_levels();
      update_node_level( n );
    };

    auto const update_level_of_existing_node = [&]( node const& n, const auto& old_children ) {
      (void)old_children;
      ntk.resize_levels();
      update_node_level( n );
    };

    auto const update_level_of_deleted_node = [&]( const auto& n ) {
      ntk.set_level( n, -1 );
    };

    add_event = ntk.events().register_add_event( update_level_of_new_node );
    modified_event = ntk.events().register_modified_event( update_level_of_existing_node );
    delete_event = ntk.events().register_delete_event( update_level_of_deleted_node );
  }

  /* maybe should move to depth_view */
  void update_node_level( node const& n, bool top_most = true )
  {
    uint32_t curr_level = ntk.level( n );

    uint32_t max_level = 0;
    ntk.foreach_fanin( n, [&]( const auto& f ) {
      auto const p = ntk.get_node( f );
      auto const fanin_level = ntk.level( p );
      if ( fanin_level > max_level )
      {
        max_level = fanin_level;
      }
    } );
    ++max_level;

    if ( curr_level != max_level )
    {
      ntk.set_level( n, max_level );

      /* update only one more level */
      if ( top_most )
      {
        ntk.foreach_fanout( n, [&]( const auto& p ) {
          update_node_level( p, false );
        } );
      }
    }
  }

private:
  Ntk& ntk;

  boptimizer_params const& ps;
  boptimizer_stats& st;
  engine_st_t& engine_st;
  collector_st_t& collector_st;

  /* temporary statistics for progress bar */
  uint32_t candidates{ 0 };
  double last_gain{ 0 };

  /* events */
  std::shared_ptr<typename network_events<Ntk>::add_event_type> add_event;
  std::shared_ptr<typename network_events<Ntk>::modified_event_type> modified_event;
  std::shared_ptr<typename network_events<Ntk>::delete_event_type> delete_event;

  std::vector<gate> _gates;
};


#pragma endregion boptimizer_impl

} /* namespace detail */

/*! \brief Window-based Boolean optimizer. */
template<support_selection_t SuppSel_t = support_selection_t::GRE, uint32_t SizeSupp=6u, uint32_t SizeFanin=SizeSupp>
void boptimize_klut( lig_network& ntk, boptimizer_params const& ps = {}, boptimizer_stats* pst = nullptr )
{
  using Ntk = lig_network;
  static constexpr uint32_t nPisLoc = 16u;
  static constexpr uint32_t nPisGlb = 10u;

  using bopt_view_t = fanout_view<depth_view<Ntk>>;
  depth_view<Ntk> depth_view{ ntk };
  bopt_view_t bopt_view{ depth_view };

  using resyn_params_t = rils::lig_resyn_static_params_for_sim_resub_static<bopt_view_t, SuppSel_t, nPisGlb, SizeSupp, SizeFanin>;
  using signature_t = kitty::static_truth_table<nPisGlb>;
  
  exact_library_params eps;
  eps.np_classification = false;
  eps.compute_dc_classes = true;

  using resyn_engine_t = rils::lig_resyn_decompose<bopt_view_t, signature_t, resyn_params_t, SuppSel_t>;

  if ( ps.max_pis <= nPisLoc )
  {
    /* only non odc optimized so far */

    using validator_t = circuit_validator<bopt_view_t, bill::solvers::bsat2, false, true, false>;//last is false
    using WindowEngine_t = typename detail::window_boptimizer<bopt_view_t, validator_t, resyn_engine_t, SizeSupp, nPisLoc, nPisGlb>;
    using bopt_impl_t = typename detail::boptimizer_impl<bopt_view_t, WindowEngine_t>;

    boptimizer_stats st;
    typename bopt_impl_t::engine_st_t engine_st;
    typename bopt_impl_t::collector_st_t collector_st;

    bopt_impl_t p( bopt_view, ps, st, engine_st, collector_st );
    p.run();
    st.time_resub -= engine_st.time_patgen;
    st.time_total -= engine_st.time_patgen + engine_st.time_patsave;

    if ( ps.verbose )
    {
      st.report();
      collector_st.report();
      engine_st.report();
    }

    if ( pst )
    {
      *pst = st;
    }
  }
  else
  {
    printf("ERROR\n");
  }

}

/*! \brief Window-based Boolean optimizer. */
template<support_selection_t SuppSel_t = support_selection_t::GRE, uint32_t SizeSupp=6u, uint32_t SizeFanin=SizeSupp>
void boptimize_sc( scopt::scg_network& ntk, boptimizer_params const& ps = {}, boptimizer_stats* pst = nullptr )
{
  using Ntk = scopt::scg_network;
  static constexpr uint32_t nPisLoc = 16u;
  static constexpr uint32_t nPisGlb = 11u;

  using bopt_view_t = fanout_view<depth_view<Ntk>>;
  depth_view<Ntk> depth_view{ ntk };
  bopt_view_t bopt_view{ depth_view };

  using resyn_params_t = scopt::scg_resyn_static_params_for_sim_resub_static<bopt_view_t, SuppSel_t, nPisGlb, SizeSupp, SizeFanin>;
  using signature_t = kitty::static_truth_table<nPisGlb>;
  
  exact_library_params eps;
  eps.np_classification = false;
  eps.compute_dc_classes = true;

  using resyn_engine_t = scopt::scg_resyn_decompose<bopt_view_t, signature_t, resyn_params_t, SuppSel_t>;

  if ( ps.max_pis <= nPisLoc )
  {
    /* only non odc optimized so far */

    using validator_t = circuit_validator<bopt_view_t, bill::solvers::bsat2, false, true, false>;//last is false
    using WindowEngine_t = typename detail::window_boptimizer<bopt_view_t, validator_t, resyn_engine_t, SizeSupp, nPisLoc, nPisGlb>;
    using bopt_impl_t = typename detail::boptimizer_impl<bopt_view_t, WindowEngine_t>;

    boptimizer_stats st;
    typename bopt_impl_t::engine_st_t engine_st;
    typename bopt_impl_t::collector_st_t collector_st;

    bopt_impl_t p( bopt_view ,ps, st, engine_st, collector_st );
    p.run();
    st.time_resub -= engine_st.time_patgen;
    st.time_total -= engine_st.time_patgen + engine_st.time_patsave;

    if ( ps.verbose )
    {
      st.report();
      collector_st.report();
      engine_st.report();
    }

    if ( pst )
    {
      *pst = st;
    }
  }
  else
  {
    printf("ERROR\n");
  }

}

//} /* namespace rils */

} /* namespace mockturtle */


