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
  \file analyzer.hpp
  \brief Generic analyzer framework

  \author Eleonora Testa
  \author Heinz Riener
  \author Mathias Soeken
  \author Shubham Rai
  \author Siang-Yun (Sonia) Lee
*/

#pragma once

#include "../traits.hpp"
#include "../utils/progress_bar.hpp"
#include "../utils/stopwatch.hpp"
#include "../views/depth_view.hpp"
#include "../views/fanout_view.hpp"

#include "detail/resub_utils.hpp"
#include "dont_cares.hpp"
#include "reconv_cut.hpp"

#include <vector>

namespace mockturtle
{

/*! \brief Parameters for analyzer.
 *
 * The data structure `analyzer_params` holds configurable parameters with
 * default arguments for `analyzer`.
 */
struct analyzer_params
{
  /*! \brief Maximum number of PIs of reconvergence-driven cuts. */
  uint32_t max_pis{ 8 };

  /*! \brief Maximum number of divisors to consider. */
  uint32_t max_divisors{ 150 };

  /*! \brief Maximum number of nodes added by analyzer. */
  uint32_t max_inserts{ 2 };

  /*! \brief Maximum fanout of a node to be considered as root. */
  uint32_t skip_fanout_limit_for_roots{ 1000 };

  /*! \brief Maximum fanout of a node to be considered as divisor. */
  uint32_t skip_fanout_limit_for_divisors{ 100 };

  /*! \brief Show progress. */
  bool progress{ false };

  /*! \brief Be verbose. */
  bool verbose{ false };

  /*! \brief IG. */
  bool useInfo{ false };

  /****** window-based resub engine ******/

  /*! \brief Use don't cares for optimization. Only used by window-based resub engine. */
  bool use_dont_cares{ false };

  /*! \brief Window size for don't cares calculation. Only used by window-based resub engine. */
  uint32_t window_size{ 12u };

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
};

/*! \brief Statistics for analyzer.
 *
 * The data structure `analyzer_stats` provides data collected by running
 * `analyzer`.
 */
struct analyzer_stats
{
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{ 0 };

  uint32_t nXXLMFFC{0};

  /*! \brief Accumulated runtime of the divisor collector. */
  stopwatch<>::duration time_divs{ 0 };

  /*! \brief Accumulated runtime of the resub engine. */
  stopwatch<>::duration time_resub{ 0 };

  /*! \brief Accumulated runtime of the callback function. */
  stopwatch<>::duration time_callback{ 0 };

  /*! \brief Total number of divisors. */
  uint64_t num_total_divisors{ 0 };

  /*! \brief Total number of gain. */
  uint64_t estimated_gain{ 0 };

  /*! \brief Initial network size (before analyzer). */
  uint64_t initial_size{ 0 };

  std::vector<uint32_t> hist;

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
    fmt::print( "[i]       ResubEngine : {:>5.2f} secs\n", to_seconds( time_resub ) );
    fmt::print( "[i]       callback    : {:>5.2f} secs\n", to_seconds( time_callback ) );
    fmt::print( "[i]     =========================\n\n" );
    // clang-format on
  }
};

namespace detail
{

struct analyzer_collector_stats
{
  /*! \brief Total number of leaves. */
  uint64_t num_total_leaves{ 0 };

  /*! \brief Accumulated runtime for cut computation. */
  stopwatch<>::duration time_cuts{ 0 };

  /*! \brief Accumulated runtime for mffc computation. */
  stopwatch<>::duration time_mffc{ 0 };

  /*! \brief Accumulated runtime for divisor computation. */
  stopwatch<>::duration time_divs{ 0 };

  void report() const
  {
    // clang-format off
    fmt::print( "[i] <DivCollector: analyzer_divisor_collector>\n" );
    fmt::print( "[i]     #leaves = {:6d}\n", num_total_leaves );
    fmt::print( "[i]     ======== Runtime ========\n" );
    fmt::print( "[i]     reconv. cut : {:>5.2f} secs\n", to_seconds( time_cuts ) );
    fmt::print( "[i]     MFFC        : {:>5.2f} secs\n", to_seconds( time_mffc ) );
    fmt::print( "[i]     divs collect: {:>5.2f} secs\n", to_seconds( time_divs ) );
    fmt::print( "[i]     =========================\n\n" );
    // clang-format on
  }
};

/*! \brief Prepare the three public data members `leaves`, `divs` and `mffc`
 * to be ready for usage.
 *
 * `leaves`: sufficient support for all divisors
 * `divs`: divisor nodes that can be used for analyzer
 * `mffc`: MFFC nodes which are needed to do simulation from
 * `leaves`, through `divs` and `mffc` until the root node,
 * but should be excluded from analyzer.
 * The last element of `mffc` is always the root node.
 *
 * `divs` and `mffc` are in topological order.
 *
 * \param MffcMgr Manager class to compute the potential gain if a
 * analyzer exists (number of MFFC nodes when the cost function is circuit size).
 * \param MffcRes Typename of the return value of `MffcMgr`.
 * \param cut_comp Manager class to compute reconvergence-driven cuts.
 */
template<class Ntk, class MffcMgr = node_mffc_inside<Ntk>, typename MffcRes = uint32_t, typename cut_comp = detail::reconvergence_driven_cut_impl<Ntk>>
class analyzer_divisor_collector
{
public:
  using stats = analyzer_collector_stats;
  using mffc_result_t = MffcRes;
  using node = typename Ntk::node;

  using cut_comp_parameters_type = typename cut_comp::parameters_type;
  using cut_comp_statistics_type = typename cut_comp::statistics_type;

public:
  explicit analyzer_divisor_collector( Ntk const& ntk, analyzer_params const& ps, stats& st )
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

    /* note: different from the previous version, now we do not add MFFC nodes into divs */
    assert( root == mffc.at( mffc.size() - 1u ) );
    /* note: this assertion makes sure window_simulator does not go out of bounds */
    assert( divs.size() + mffc.size() - leaves.size() <= ps.max_divisors - ps.max_pis );

    return true;
  }

private:
  Ntk const& ntk;
  analyzer_params ps;
  stats& st;

  cut_comp cuts;
  cut_comp_statistics_type cuts_st;

public:
  std::vector<node> leaves;
  std::vector<node> divs;
  std::vector<node> mffc;
};

/*! \brief The top-level analyzer framework.
 *
 * \param ResubEngine The engine that computes the resubtitution for a given root
 * node and divisors. One can choose from `window_based_resub_engine` which
 * does complete simulation within small windows, or `simulation_based_resub_engine`
 * which does partial simulation on the whole circuit.
 *
 * \param DivCollector Collects divisors near a given root node, and compute
 * the potential gain (MFFC size or its variants).
 * Currently only `analyzer_divisor_collector` is implemented, but
 * a frontier-based approach may be integrated in the future.
 * When using `window_based_resub_engine`, the `DivCollector` should prepare
 * three public data members: `leaves`, `divs`, and `mffc` (see documentation
 * of `analyzer_divisor_collector` for details). When using `simulation_based_resub_engine`,
 * only `divs` is needed.
 */
template<class Ntk, class DivCollector = analyzer_divisor_collector<Ntk>>
class analyzer_impl
{
public:
  using collector_st_t = typename DivCollector::stats;
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;
  using resub_callback_t = std::function<bool( Ntk&, node const&, signal const& )>;
  using mffc_result_t = uint32_t;

  /*! \brief Constructor of the top-level analyzer framework.
   *
   * \param ntk The network to be optimized.
   * \param ps Resubstitution parameters.
   * \param st Top-level analyzer statistics.
   * \param engine_st Statistics of the analyzer engine.
   * \param collector_st Statistics of the divisor collector.
   * \param callback Callback function when a analyzer is found.
   */
  explicit analyzer_impl( Ntk& ntk, analyzer_params const& ps, analyzer_stats& st )
      : ntk( ntk ), ps( ps ), st( st )
  {
    st.initial_size = ntk.num_gates();

  }

  ~analyzer_impl()
  {
    ntk.events().release_add_event( add_event );
    ntk.events().release_modified_event( modified_event );
    ntk.events().release_delete_event( delete_event );
  }

  void run()
  {
    std::vector<uint32_t> hist;
    
    stopwatch t( st.time_total );
    collector_st_t collector_st;
    /* start the managers */
    DivCollector collector( ntk, ps, collector_st );

    auto const size = ntk.num_gates();
    ntk.foreach_gate( [&]( auto const& n, auto i ) {
      if ( i >= size )
      {
        return false; /* terminate */
      }


      /* compute cut, collect divisors, compute MFFC */
      mffc_result_t potential_gain;
      const auto collector_success = call_with_stopwatch( st.time_divs, [&]() {
        return collector.run( n, potential_gain );
      } );
      if ( !collector_success )
      {
        return true; /* next */
      }
      else
      {
        for( uint32_t i{hist.size()}; i<potential_gain+1; ++i )
          hist.push_back(0);
        hist[potential_gain]+=1;
      }

      return true; /* next */
    } );
    
    uint32_t S4{0};
    for( uint32_t i{0}; i<hist.size(); ++i )
    {
      //printf("[%d:%d]", i, hist[i]);
      if( i>4 )
        S4+=hist[i];
    }
    st.nXXLMFFC = S4;
  }

private:
  Ntk& ntk;

  analyzer_params const& ps;
  analyzer_stats& st;

  /* temporary statistics for progress bar */
  uint32_t candidates{ 0 };
  uint32_t last_gain{ 0 };

  /* events */
  std::shared_ptr<typename network_events<Ntk>::add_event_type> add_event;
  std::shared_ptr<typename network_events<Ntk>::modified_event_type> modified_event;
  std::shared_ptr<typename network_events<Ntk>::delete_event_type> delete_event;
};

} /* namespace detail */

/*! \brief Window-based Boolean analyzer with default resub functor (only div0). */
template<class Ntk>
void default_analyzer( Ntk& ntk, analyzer_params const& ps = {}, analyzer_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_clear_values_v<Ntk>, "Ntk does not implement the clear_values method" );
  static_assert( has_fanout_size_v<Ntk>, "Ntk does not implement the fanout_size method" );
  static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin method" );
  static_assert( has_foreach_gate_v<Ntk>, "Ntk does not implement the foreach_gate method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_get_constant_v<Ntk>, "Ntk does not implement the get_constant method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement the is_complemented method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_make_signal_v<Ntk>, "Ntk does not implement the make_signal method" );
  static_assert( has_set_value_v<Ntk>, "Ntk does not implement the set_value method" );
  static_assert( has_set_visited_v<Ntk>, "Ntk does not implement the set_visited method" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_substitute_node_v<Ntk>, "Ntk does not implement the substitute_node method" );
  static_assert( has_value_v<Ntk>, "Ntk does not implement the value method" );
  static_assert( has_visited_v<Ntk>, "Ntk does not implement the visited method" );

  analyzer_stats st;
 using resub_view_t = fanout_view<depth_view<Ntk>>;
  depth_view depth_view{ntk};
  fanout_view fanout_view{depth_view};
    detail::analyzer_impl<resub_view_t> p( fanout_view, ps, st );
    p.run();
  *pst = st;

}

} /* namespace mockturtle */
