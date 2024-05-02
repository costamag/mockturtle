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
  \file rewrub.hpp
  \brief Inplace rewrub

  \author Andrea Costamagna
*/

#pragma once

#include "../traits.hpp"
#include "../utils/cost_functions.hpp"
#include "../utils/index_list.hpp"
#include "../utils/spfd_utils.hpp"
#include "../utils/node_map.hpp"
#include "detail/resub_utils.hpp"
#include "../utils/stopwatch.hpp"
#include "../views/color_view.hpp"
#include "../views/depth_view.hpp"
#include "../views/fanout_view.hpp"
#include "../views/window_view.hpp"
#include "cleanup.hpp"
#include "cut_enumeration.hpp"
#include "cut_enumeration/rewrite_cut.hpp"
#include "reconv_cut.hpp"
#include "simulation.hpp"

#include <fmt/format.h>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>

namespace mockturtle
{

std::mt19937 RNGRWS( 5 );


struct pLibrary_t
{
  pLibrary_t( std::string const& library )
  {
    kitty::static_truth_table<4u> tt;
    int i{0};
    std::string line;
    std::ifstream fTts ( library + ".tts" );
    if (fTts.is_open())
    {
      while ( std::getline (fTts,line) )
      {
        kitty::create_from_binary_string( tt, line );
        _pClassMap[tt._bits]=i++;
      }
      fTts.close();
    }
    else
    {
      printf("not found\n");
    }

    std::ifstream fAreas (library + ".area");
    if (fAreas.is_open())
    {
      while ( std::getline (fAreas,line) )
      {
        _areas.push_back( std::stof( line ) );
      }
      fAreas.close();
    }
    else
    {
      printf("not found\n");
    }


    std::ifstream fLists ( library + ".list" );
    if (fLists.is_open())
    {
      while ( std::getline (fLists,line) )
      {
        std::vector<uint32_t> list;
        uint32_t number;
        std::istringstream line_stream(line);
        while (line_stream >> number)
        {
            list.push_back(number);
        }
        _idlists.push_back( list );
      }
      fLists.close();
    }
    else
    {
      printf("not found\n");
    }
  }

  template<class TT>
  std::optional<uint32_t> get_key( TT const& tt )
  {
    uint64_t repr = tt._bits&0xFFFF;
    if( _pClassMap.find(repr) != _pClassMap.end() )
      return _pClassMap[tt._bits];
    return std::nullopt;
  }

  template<class TT>
  std::optional<float> get_area( TT const& tt )
  {
    auto key = get_key( tt );
    if( key ) 
      return _areas[*key];
    return std::nullopt;
  }

  /* objects */
  std::vector<std::vector<uint32_t>> _idlists;
  std::vector<float> _areas;
  std::unordered_map<uint64_t, uint32_t> _pClassMap;
};

/*! \brief Parameters for rewrub.
 *
 * The data structure `rewrub_sc_params` holds configurable parameters with
 * default arguments for `rewrub`.
 */
struct rewrub_sc_params
{
  rewrub_sc_params()
  {
    /* 0 < Cut limit < 16 */
    cut_enumeration_ps.cut_limit = 8;
    cut_enumeration_ps.minimize_truth_table = true;
  }

  /*! \brief Cut enumeration parameters. */
  cut_enumeration_params cut_enumeration_ps{};

  /*! \brief If true, candidates are only accepted if they do not increase logic depth. */
  bool preserve_depth{ false };

  /*! \brief Allow rewrub with multiple structures */
  bool allow_multiple_structures{ true };

  /*! \brief Allow zero-gain substitutions */
  bool allow_zero_gain{ false };

  /*! \brief Use satisfiability don't cares for optimization. */
  bool use_dont_cares{ false };

  /*! \brief Window size for don't cares calculation. */
  uint32_t window_size{ 16u };

  /*! \brief Maximum number of PIs of reconvergence-driven cuts. */
  uint32_t max_pis{ 16 };

  /*! \brief Maximum number of divisors to consider. */
  uint32_t max_divisors{ 150 };

  /*! \brief Maximum number of nodes added by resubstitution. */
  uint32_t max_inserts{ 2 };

  /*! \brief Maximum fanout of a node to be considered as root. */
  uint32_t skip_fanout_limit_for_roots{ 1000 };

  /*! \brief Maximum fanout of a node to be considered as divisor. */
  uint32_t skip_fanout_limit_for_divisors{ 100 };

  /*! \brief Number of sampling of the functional cuts. */
  uint32_t num_samplings{1};

  /*! \brief Be verbose. */
  bool verbose{ false };
};

/*! \brief Statistics for rewrub.
 *
 * The data structure `rewrub_sc_stats` provides data collected by running
 * `rewrub`.
 */
struct rewrub_sc_stats
{
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{ 0 };

  /*! \brief Expected gain. */
  uint32_t estimated_gain{ 0 };

  /*! \brief Candidates */
  uint32_t candidates{ 0 };

  void report() const
  {
    //std::cout << fmt::format( "[i] total time       = {:>5.2f} secs\n", to_seconds( time_total ) );
  }
};

namespace detail
{

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

template<typename Ntk>
class node_mffc_inside2
{
public:
  using node = typename Ntk::node;

public:
  explicit node_mffc_inside2( Ntk const& ntk )
      : ntk( ntk )
  {
    static_assert( has_incr_fanout_size_v<Ntk>, "Ntk does not implement the incr_fanout_size method" );
    static_assert( has_decr_fanout_size_v<Ntk>, "Ntk does not implement the decr_fanout_size method" );
    static_assert( has_fanout_size_v<Ntk>, "Ntk does not implement the fanout_size method" );
    static_assert( has_incr_trav_id_v<Ntk>, "Ntk does not implement the incr_trav_id method" );
    static_assert( has_trav_id_v<Ntk>, "Ntk does not implement the trav_id method" );
    static_assert( has_set_visited_v<Ntk>, "Ntk does not implement the set_visited method" );
    static_assert( has_visited_v<Ntk>, "Ntk does not implement the visited method" );
    static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin method" );
    static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
    static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  }

  template<typename Fn>
  float call_on_mffc_and_count( node const& n, std::vector<node> const& leaves, Fn&& fn )
  {
    /* increment the fanout counters for the leaves */
    ntk.incr_trav_id();
    //printf("TRAVID %d\n", ntk.trav_id());
    for ( const auto& l : leaves )
    {
      ntk.incr_fanout_size( l );
      ntk.set_visited( l, ntk.trav_id() ); 
    }

    //printf("DEREF: ");
    /* dereference the node */
    auto count1 = node_deref_rec( n );
    //printf("\n");

    /* call `fn` on MFFC nodes */
    node_mffc_cone_rec( n, true, fn );

    /* reference it back */
    auto count2 = node_ref_rec( n );
    (void)count2;

    float eps = 0.1;
    assert( abs(count1-count2) <= eps );

    for ( const auto& l : leaves )
      ntk.decr_fanout_size( l );

    return count1;
  }

  float run( node const& n, std::vector<node> const& leaves, std::vector<node>& inside )
  {
    //printf("COLLECTING..\n");
    inside.clear();
    return call_on_mffc_and_count( n, leaves, [&]( node const& m ) { inside.push_back( m ); } );
  }

private:
  /* ! \brief Dereference the node's MFFC */
  float node_deref_rec( node const& n )
  {
    //printf("REC @ %d ", n );
    if ( ntk.is_pi( n ) )
      return 0;

    float counter = ntk.get_area( n );
    ntk.foreach_fanin( n, [&]( const auto& f ) {
      auto const& p = ntk.get_node( f );
      //printf( "fanin %d FO[0]=%d", p, ntk.fanout_size( p ) );
      ntk.decr_fanout_size( p );
      //printf( "=>FO[1]=%d\n", ntk.fanout_size( p ) );

      if ( ntk.fanout_size( p ) == 0 )
      {
        counter += node_deref_rec( p );
      }
    } );

    return counter;
  }

  /* ! \brief Reference the node's MFFC */
  float node_ref_rec( node const& n )
  {
    if ( ntk.is_pi( n ) )
      return 0;

    float counter = ntk.get_area( n );
    ntk.foreach_fanin( n, [&]( const auto& f ) {
      auto const& p = ntk.get_node( f );

      auto v = ntk.fanout_size( p );
      ntk.incr_fanout_size( p );
      if ( v == 0 )
      {
        counter += node_ref_rec( p );
      }
    } );

    return counter;
  }

  template<typename Fn>
  void node_mffc_cone_rec( node const& n, bool top_most, Fn&& fn )
  {
    /* skip visited nodes */
    if ( ntk.visited( n ) == ntk.trav_id() )
    {
      return;
    }
    ntk.set_visited( n, ntk.trav_id() );

    if ( !top_most && ( ntk.is_pi( n ) || ntk.fanout_size( n ) > 0 ) )
    {
      return;
    }

    /* recurse on children */
    ntk.foreach_fanin( n, [&]( const auto& f ) {
      node_mffc_cone_rec( ntk.get_node( f ), false, fn );
    } );

    /* collect the internal nodes */
    fn( n );
  }

private:
  Ntk const& ntk;
}; /* rils_node_mffc_inside */

template<class Ntk, class MffcMgr = node_mffc_inside2<Ntk>, typename MffcRes = float, typename cut_comp = detail::reconvergence_driven_cut_impl<Ntk>>
class divisor_collector2_t
{
public:
  using stats = collector_stats;
  using mffc_result_t = MffcRes;
  using node = typename Ntk::node;

  using cut_comp_parameters_type = typename cut_comp::parameters_type;
  using cut_comp_statistics_type = typename cut_comp::statistics_type;

public:
  explicit divisor_collector2_t( Ntk const& ntk, rewrub_sc_params const& ps, stats& st )
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

  //  std::cout << "L: " << std::endl;
  //  for( auto x : leaves )
  //  {
  //    printf("l) %d ", x );
  //    std::cout << " F" << ntk.get_binding(x).id << "(";
  //    ntk.foreach_fanin( x, [&]( auto const& f ) {
  //      std::cout << " " << ntk.get_node(f);
  //    } );
  //
  //    printf(")\n");
  //  }

    /* collect the MFFC */
    MffcMgr mffc_mgr( ntk );
    potential_gain = call_with_stopwatch( st.time_mffc, [&]() {
      return mffc_mgr.run( n, leaves, mffc );
    } );

//    std::cout << "M: " << std::endl;
//    for( auto x : mffc )
//    {
//      printf("l) %d ", x );
//      std::cout << " F" << ntk.get_binding(x).id << "(";
//      ntk.foreach_fanin( x, [&]( auto const& f ) {
//        std::cout << " " << ntk.get_node(f);
//      } );
//
//      printf(")\n");
//    }

    /* collect the divisor nodes in the cut */
    bool div_comp_success = collect_divisors( n );

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
    if ( ntk.value( n ) != ntk.trav_id() ) /*  */
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
      ntk.set_value( t, ntk.trav_id() );
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

    //if( ps.use_wings )
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
        desp.emplace_back( nd );
    }

    return true;
  }

private:
  Ntk const& ntk;
  rewrub_sc_params ps;
  stats& st;

  cut_comp cuts;
  cut_comp_statistics_type cuts_st;

public:
  std::vector<node> leaves;
  std::vector<node> divs;
  std::vector<node> mffc;
  std::vector<node> desp;
};

template<class Ntk, uint32_t W>
class rewrub_sc_impl
{
  static constexpr uint32_t num_vars = 4u;
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

  using DivCollector = detail::divisor_collector2_t<Ntk>;
  using collector_st_t = typename DivCollector::stats;
  using mffc_result_t = float;

  using network_cuts_t = dynamic_network_cuts<Ntk, num_vars, true, cut_enumeration_rewrite_cut>;
  using cut_manager_t = detail::dynamic_cut_enumeration_impl<Ntk, num_vars, true, cut_enumeration_rewrite_cut>;
  using cut_t = typename network_cuts_t::cut_t;
  using node_data = typename Ntk::storage::element_type::node_type;

public:
  rewrub_sc_impl( Ntk& ntk, pLibrary_t& database, rewrub_sc_params const& ps, rewrub_sc_stats& st )// Library&& library, NodeCostFn const& cost_fn
      : ntk( ntk ), _database(database), ps( ps ), st( st ), _required( ntk, std::numeric_limits<float>::max() ), _arrival( ntk, 0 ), _ttW(ntk) //library( library ), cost_fn( cost_fn ),
  {
    // initialize reference simulation patterns
    for( int i{0}; i<W; ++i )
    {
      kitty::create_nth_var( _xsW[i], i );
    }
    for( int i{0}; i<4; ++i )
    {
      kitty::create_nth_var( _xs4[i], i );
    }
    //_tts[0] = _xsW[0].construct();

    /* timing information */
    _max_delay = ps.preserve_depth ? ntk.compute_worst_delay() : std::numeric_limits<float>::max();
  }

  ~rewrub_sc_impl()
  {
  }

  void run()
  {
    stopwatch t( st.time_total );

    if ( ps.preserve_depth )
    {
      compute_arrival();
      compute_required( _max_delay );
    }

    perform_rewrubbing();

    st.estimated_gain = _estimated_gain;
    st.candidates = _candidates;
  }

  #pragma region timing 

  float compute_arrival_rec( node n )
  {
    if( ntk.visited(n) > 0 )
    {
      return _arrival[n];
    }


    if ( ntk.has_binding( n ) )
    {
      auto const& g = ntk.get_binding( n );
      float gate_delay = 0;
      ntk.foreach_fanin( n, [&]( auto const& f, auto i ) {
        float arr_fanin = compute_arrival_rec( ntk.get_node(f) );
        gate_delay = std::max( gate_delay, (float)( arr_fanin + std::max( g.pins[i].rise_block_delay, g.pins[i].fall_block_delay ) ) );
      } );
      _arrival[n] = gate_delay;
    }
    else
    {
      float gate_delay = 1;
      ntk.foreach_fanin( n, [&]( auto const& f, auto i ) {
        float arr_fanin = compute_arrival_rec( ntk.get_node(f) );
        gate_delay = std::max( gate_delay, (float)( arr_fanin + 1u ) );
      } );
      _arrival[n] = gate_delay;
    }
    ntk.set_visited( n, 1u );
    return _arrival[n];
  }

  float compute_arrival()
  {
    ntk.clear_visited();
    _arrival.reset();
    ntk.foreach_pi( [&]( auto const& n, auto i ) {
      _arrival[n]=0.0;
      ntk.set_visited( n, 1u );
    } );
    _arrival[0]=0.0;
    ntk.set_visited( 0, 1u );


    float max_delay=0;
    ntk.foreach_po( [&]( auto const& no, auto i ) {
      auto out_del = compute_arrival_rec( ntk.get_node(no) );
      if( out_del > max_delay )
      {
        max_delay = out_del;
      }
    } );
    
    ntk.clear_visited();

    return max_delay;
  }

  float compute_required_rec( node n, float max_delay )
  {
    if( ntk.visited(n) > 0  )
    {
      //std::cout << n << " " << _required[n] << std::endl;
      return _required[n];
    }

    if ( ntk.has_binding( n ) )
    {
      float gate_required = max_delay;

      ntk.foreach_fanout( n, [&]( auto const& f, auto i ) {
        float req_fanout = compute_required_rec( ntk.get_node(f), max_delay );
        auto nfo = ntk.get_node(f);
        auto const& g = ntk.get_binding( nfo );

        uint32_t ig;
        ntk.foreach_fanin( nfo, [&]( auto const& fo, auto io ) {
          if( ntk.get_node(fo) == n )
          {
            ig = io;
            return;
          }
        });

        //std::cout << "[" << n << "]@fo" << ntk.get_node(f) << " " << req_fanout << std::endl;
        gate_required = std::min( gate_required, (float)( req_fanout - std::max( g.pins[ig].rise_block_delay, g.pins[ig].fall_block_delay ) ) );
      } );
      _required[n] = gate_required;
    }
    else
    {
      float gate_required = max_delay;
      ntk.foreach_fanout( n, [&]( auto const& f, auto i ) {
        float req_fanout = compute_required_rec( ntk.get_node(f), max_delay );
        gate_required = std::min( gate_required, (float)( req_fanout - 1u ) );
      } );
      _required[n] = gate_required;
    }
    ntk.set_visited( n, 1u );

    //std::cout << "r " << n << " " << _required[n] << std::endl;
    return _required[n];
  }

  void compute_required( float max_delay )
  {

    //ntk.print();
    //std::cout << "max del " <<max_delay << std::endl;

    ntk.clear_visited();
    _required.reset();
    ntk.foreach_po( [&]( auto const& fo, auto i ) {
      auto n = ntk.get_node(fo);
      _required[n]=max_delay;
      ntk.set_visited( n, 1u );
    } );


    ntk.foreach_pi( [&]( auto const& ni, auto i ) {
      auto req = compute_required_rec( ntk.get_node(ni), max_delay );
      //std::cout << "PI" << ni << " " << req << std::endl;
    } );
    
    ntk.clear_visited();

  }

  void print_slack()
  {
    ntk.foreach_gate( [&]( auto const& n, auto i ) { 
      printf("%4d %f\n", n, _required[n]-_arrival[n]);
    });
  }

  #pragma endregion timing

  #pragma region opto
  struct opto_candidate_t
  {
    uint32_t id;
    std::array<signal, num_vars> leaves;
    std::array<uint8_t, num_vars> permutation;
    float reward;
  };

  float measure_mffc_ref( node const& n, std::array<signal,num_vars> const& cut )
  {
    /* reference cut leaves */
    for ( auto leaf : cut )
    {
      ntk.incr_fanout_size( ntk.get_node(leaf) );
    }

    float mffc_size = static_cast<float>( recursive_ref( n ) );

    /* dereference leaves */
    for ( auto leaf : cut )
    {
      ntk.decr_fanout_size( ntk.get_node(leaf) );
    }

    return mffc_size;
  }

  float measure_mffc_deref( node const& n, std::array<signal,num_vars> const& cut )
  {
    /* reference cut leaves */
    for ( auto leaf : cut )
    {
      ntk.incr_fanout_size( ntk.get_node(leaf) );
    }

    float mffc_size = static_cast<float>( recursive_deref( n ) );

    /* dereference leaves */
    for ( auto leaf : cut )
    {
      ntk.decr_fanout_size( ntk.get_node(leaf) );
    }

    return mffc_size;
  }

  float recursive_deref( node const& n )
  {
    /* terminate? */
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      return 0;

    /* recursively collect nodes */
    float value{ ntk.get_area( n ) };
    ntk.foreach_fanin( n, [&]( auto const& s ) {
      if ( ntk.decr_fanout_size( ntk.get_node( s ) ) == 0 )
      {
        value += recursive_deref( ntk.get_node( s ) );
      }
    } );
    return value;
  }

  float recursive_ref( node const& n )
  {
    /* terminate? */
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      return 0;

    /* recursively collect nodes */
    float value{ ntk.get_area( n ) };
    ntk.foreach_fanin( n, [&]( auto const& s ) {
      if ( ntk.incr_fanout_size( ntk.get_node( s ) ) == 0 )
      {
        value += recursive_ref( ntk.get_node( s ) );
      }
    } );
    return value;
  }

  float area_contained_mffc( node n, std::array<signal, num_vars> const& leaves )
  {
    /* measure the MFFC contained in the cut */
    float mffc_size = measure_mffc_deref( n, leaves );
    /* restore contained MFFC */
    measure_mffc_ref( n, leaves );

    return mffc_size;
  }  

  std::optional<opto_candidate_t> find_structural_rewriting( cut_manager_t & cut_manager, network_cuts_t & cuts, node const& n )
  {
    cut_manager.clear_cuts( n );
    cut_manager.compute_cuts( n );

    std::vector<opto_candidate_t> cands;
    float best_reward = std::numeric_limits<float>::min();

    uint32_t cut_index = 0;
    for ( auto& cut : cuts.cuts( ntk.node_to_index( n ) ) )
    {
      opto_candidate_t cand;
      /* skip trivial cut */
      if ( ( cut->size() == 1 && *cut->begin() == ntk.node_to_index( n ) ) )
      {
        ++cut_index;
        continue;
      }

      /* Boolean matching */
      auto config = kitty::exact_p_canonization( cuts.truth_table( *cut ) );
      auto repr = std::get<0>( config );
      auto nega = std::get<1>( config );
      auto perm = std::get<2>( config );
      std::array<uint8_t, num_vars> permutation;

      assert( nega == 0u );
      for ( auto j = 0u; j < num_vars; ++j )
      {
        permutation[perm[j]] = j;
      }

      /* save output negation to apply */

      {
        auto j = 0u;
        for ( auto const leaf : *cut )
        {
          cand.leaves[permutation[j++]] = ntk.make_signal( ntk.index_to_node( leaf ) );
        }

        while ( j < num_vars )
          cand.leaves[permutation[j++]] = ntk.get_constant( false );
      }

      /* resynthesis cost */
      auto cost = _database.get_area( repr );
      if( !cost )
      {
        continue;
      }
      else
      {
        cand.id = *_database.get_key( repr );
      }

      /* resynthesis reward */
      float area_mffc = area_contained_mffc( n, cand.leaves );
      //std::cout << "S(amffc " << area_mffc << std::endl;

      cand.reward = area_mffc - *cost;
      if( cand.reward > 0 && *cost > 0 )
      {
      //  std::cout << "cost:" << *cost << "removed:" << area_mffc << "=> reward " << cand.reward << std::endl;
        if( cand.reward > best_reward )
        {
          best_reward = cand.reward;
          cands = {cand};
        }
        else if( cand.reward == best_reward )
        {
          cands.push_back( cand );
        }
      }
    }
    /* sample from the solutions */
    if( cands.size() > 0 )
    {
      std::uniform_int_distribution<> distrib( 0, cands.size() - 1 );
      int idx = distrib( RNGRWS );
      return cands[idx];
    }

    return std::nullopt;
  }

  void simulate_window( std::vector<node> const& leaves, std::vector<node> const& divs, std::vector<node> const& mffc, node n )
  {
    _ttW.reset();
    _ttW[0]=_xsW[0].construct();
    int i{0};
    for( auto l : leaves )
    {
      _ttW[l] = _xsW[i++];
      //printf("l %3d:",l); kitty::print_binary(_ttW[l]); printf("\n");

    }

    std::vector<kitty::static_truth_table<W>> children; 
    for( auto d : divs )
    {
      if( !ntk.is_constant(d) && std::find( leaves.begin(), leaves.end(), d ) == leaves.end() )
      {
        children.clear();
        ntk.foreach_fanin( d, [&]( auto const& f ) {
          children.push_back( _ttW[ntk.get_node(f)] );
        } );
        _ttW[d] = ntk.compute( d, children.begin(), children.end() );
        //printf("d %3d:",d); kitty::print_binary(_ttW[d]); 
//        ntk.foreach_fanin( d, [&]( auto const& f ) {
//          std::cout << " " << ntk.get_node(f);
//        } );
//        std::cout << " id" << ntk.get_binding(d).id;
//        printf("\n");
      }
    }

    for( auto d : mffc )
    {
      if( !ntk.is_constant(d) && std::find( leaves.begin(), leaves.end(), d ) == leaves.end() )
      {
        children.clear();
        ntk.foreach_fanin( d, [&]( auto const& f ) {
          children.push_back( _ttW[ntk.get_node(f)] );
        } );

        _ttW[d] = ntk.compute( d, children.begin(), children.end() );

//        printf("m %3d:",d); kitty::print_binary(_ttW[d]); 
//        ntk.foreach_fanin( d, [&]( auto const& f ) {
//          std::cout << " " << ntk.get_node(f);
//        } );
//        std::cout << " id" << ntk.get_binding(d).id;
//
//        printf("\n");

      }
    }
    children.clear();
    ntk.foreach_fanin( n, [&]( auto const& f ) {
        children.push_back( _ttW[ntk.get_node(f)] );
    } );

    _ttW[n] = ntk.compute( n, children.begin(), children.end() );
//    printf("n %3d:",n); kitty::print_binary(_ttW[n]); printf("\n");


  }

  template<class SIM>
  std::optional<std::vector<node>> find_support_greedy( std::vector<node> const& divs, SIM const& tts )
  {
    uint32_t cost, best_cost;
    std::vector<node> best_candidates;
    std::vector<node> supp;
    _spfd.reset();


    /* add recomputation of the support */
    while ( !_spfd.is_covered() && supp.size() < 4u )
    {
      best_cost = std::numeric_limits<uint32_t>::max();
      if ( _spfd.is_saturated() )
        return std::nullopt;
      for ( uint32_t v{0}; v < divs.size(); ++v )
      {
        cost = _spfd.evaluate( tts[divs[v]] );
        if ( cost < best_cost )
        {
          best_cost = cost;
          best_candidates = { divs[v] };
        }
        else if ( cost == best_cost )
        {
          best_candidates.push_back( divs[v] );
        }
      }
      if ( best_candidates.size() == 0 )
        return std::nullopt;
      std::uniform_int_distribution<> distrib( 0, best_candidates.size() - 1 );
      int idx = distrib( RNGRWS );
      supp.push_back( best_candidates[idx] );
      _spfd.update( tts[best_candidates[idx]] );
    }

    if ( _spfd.is_covered() && supp.size() <= 4u )
    {
      std::sort( supp.begin(), supp.end() );
      return supp;
    }
    return std::nullopt;
  }

  template<class SIM>
  std::optional<std::vector<node>> find_support( std::vector<node> const& divs, node n, SIM const& tts )
  {
    _spfd.init( tts[n] );
    auto supp = find_support_greedy( divs, tts );
    return supp;
  }

  template<class SIM>
  std::tuple<kitty::static_truth_table<4u>, kitty::static_truth_table<4u>> extract_functionality( std::vector<node> const & cut, node n, SIM const& tts )
  {
    kitty::static_truth_table<4u> tt;
    kitty::static_truth_table<4u> mk;

    auto tmp = tts[n].construct();
    
    for( uint32_t m{0}; m < ( 1u << cut.size() ); ++m )
    {
      auto tmp4 = _xs4[0]|~_xs4[0];
      tmp = tmp | ~tmp;
      for( int i{0}; i<cut.size(); ++i )
      {
        if( ( ( m >> i ) & 0x1 ) == 0x1 )
        {
          tmp &= tts[cut[i]];
          tmp4 &= _xs4[i];
        }
        else
        {
          tmp &= ~tts[cut[i]];
          tmp4 &= ~_xs4[i];
        }  
      }
      int n0 = kitty::count_ones( ~tts[n] & tmp );
      int n1 = kitty::count_ones(  tts[n] & tmp );
      if( n0 > 0 && n1 == 0 )
      {
        mk |= tmp4;
      }
      else if( n1 > 0 && n0 == 0 )
      {
        tt |= tmp4;
        mk |= tmp4;
      }
      else if( n1 != 0 || n0 != 0 )
      {
        printf("not a valid support\n");
        assert( false && "Not valid support\n" );
      }
    }
    return std::make_tuple( tt, mk );
  }

  std::optional<opto_candidate_t> find_functional_rewriting( node const& n )
  {
    float mffc_area;
    collector_st_t collector_st;

    DivCollector collector( ntk, ps, collector_st );
    const auto collector_success = collector.run( n, mffc_area );
    if ( !collector_success )
    {
      return std::nullopt; /* next */
    }

    //std::cout << "F(amffc " << area_mffc << std::endl;
    std::vector<signal> leaves_sig(collector.leaves.size());
    std::transform( collector.leaves.begin(), collector.leaves.end(), leaves_sig.begin(), [&]( const node n ) {
      auto f = ntk.make_signal( n );
      return f;
    } );

//    std::cout << "RC:";
//    for( auto x : collector.leaves )
//    {
//      std::cout << x << " ";
//    }
//    std::cout << std::endl;

    //collect_mffc( mffc, n, leaves_sig );

    simulate_window( collector.leaves, collector.divs, collector.mffc, n );

    // find functional cut
    auto supp = find_support( collector.divs, n, _ttW );
    if( supp )
    {
//      std::cout << "SUPP:";
//      for( auto x : *supp )
//      {
//        std::cout << x << " ";
//      }
//      std::cout << std::endl;

      auto [func, care] = extract_functionality( *supp, n, _ttW );
      auto dontcare = ~care;

//      kitty::print_binary( func );
//      printf("\n");
//      kitty::print_binary( care );
//      printf("\n");

      std::vector<uint32_t> dcs;
      for( int bit{0}; bit<16; ++bit )
      {
        if( kitty::get_bit( dontcare, bit ) > 0 )
        {
          dcs.push_back(bit);
        }
      }

      uint64_t best_key;
      float best_area = mffc_area;
      std::vector<uint8_t> best_perm;


      for( uint32_t m{0}; m< (1<<dcs.size()); ++m )
      {
        auto tt = func;

        for( int i{0}; i<dcs.size(); ++i )
        {
          if( (m >> i)&0x1 == 0x1 )
          {
            kitty::flip_bit( tt, dcs[i] );
          }
        }

        /* p-canonize */
        auto config = kitty::exact_p_canonization( tt );

        auto repr = std::get<0>( config );
        auto neg = std::get<1>( config );
        auto perm = std::get<2>( config );

        //kitty::print_binary(func_p);
        //printf("\n");

        auto key = _database.get_key( repr );
//        printf("%f <? %f\n", _database._areas[*key], mffc_area );
        if( _database._areas[*key] < best_area )
        {
          best_key = *key;
          best_area = _database._areas[*key];
          best_perm = perm;
        }
      }

      //printf("%f >? %f\n", max_inserts, best_area);
      if( best_area < mffc_area )
      {
//        printf("]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]\n");
//        printf("]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]\n");
//        printf("]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]\n");
        opto_candidate_t cand;

        std::array<uint8_t, num_vars> permutation;

        for ( auto j = 0u; j < num_vars; ++j )
        {
          permutation[best_perm[j]] = j;
        }

        /* save output negation to apply */

        {
          auto j = 0u;
          for ( auto const leaf : *supp )
          {
            cand.leaves[permutation[j++]] = ntk.make_signal( leaf );
          }

          while ( j < num_vars )
            cand.leaves[permutation[j++]] = ntk.get_constant( false );
        }

        /* resynthesis cost */
        cand.id = best_key;

        /* resynthesis reward */
        float area_mffc = area_contained_mffc( n, cand.leaves );
        //std::cout << "S(amffc " << area_mffc << std::endl;

        //cand.reward = area_mffc - *cost;
        return cand;
      }

      }

    return std::nullopt;
  }


  large_lig_index_list resynthesize_index_list( opto_candidate_t const& cand )
  {
    large_lig_index_list index_list(4u);

    std::vector<uint32_t> lits = {0,2,4,6,8};

    auto entry = _database._idlists[ cand.id ];
    int type = 0;
    int nFins = 0;
    uint32_t sc_id;
    std::vector<uint32_t> children;
    uint32_t lit;
    
    for( int i{0}; i<entry.size(); ++i )
    {
      if( type == 0 )
      {
        nFins = entry[i];
        type=1;
      }
      else if( type == 1 )
      {
        children.push_back( lits[entry[i]] );//not accounting for 0
        if( children.size() == nFins )
        {
          type = 2;
        }
      }
      else if( type == 2 )
      {
        type = 0;
        sc_id = entry[i];
        lit = index_list.add_function( children, ntk._library[sc_id].function, ntk._library[sc_id].area, ntk._library[sc_id].id );
//        std::cout << lit << std::endl;
        lits.push_back( lit );
        children.clear();
      }
    }
    index_list.add_output( lit );
//    std::cout << "id list area " << index_list.get_area() << std::endl;
    return index_list;
  }

  signal resynthesize_sub_network( large_lig_index_list const& index_list, std::array<uint8_t,4> const& perm, std::array<signal, 4> const& leaves )
  {
    std::vector<signal> divs_sig;
    for( int i{0}; i<4u; ++i ) 
    {
      divs_sig.push_back( leaves[i] );
    }
//    std::transform( leaves.begin(), leaves.end(), divs_sig.begin(), [&]( const node n ) {
//      auto f = ntk.make_signal( n );
//      return f;
//    } );

    signal res;
    insert( ntk, divs_sig.begin(), divs_sig.end(), index_list, [&]( signal const& s ) {
      res = s;
    } );
    return res;
  }

  void perform_rewrubbing()
  {  

    /* structural cuts */
    cut_enumeration_stats cst;
    network_cuts_t cuts( ntk.size() + ( ntk.size() >> 1 ) );
    cut_manager_t cut_manager( ntk, ps.cut_enumeration_ps, cst, cuts );
    cut_manager.init_cuts();

    /* window cut */
    reconvergence_driven_cut_parameters rcuts_ps;
    rcuts_ps.max_leaves = W;
    reconvergence_driven_cut_statistics rcuts_st;

    std::array<signal, num_vars> leaves;
    std::array<signal, num_vars> best_leaves;
    std::array<uint8_t, num_vars> permutation;
    //signal best_signal;


    /*
    auto& db = library.get_database();*/

    const auto size = ntk.size(); 
    ntk.foreach_gate( [&]( auto const& n, auto i ) {

      /* exit condition */
      if (  i >= size ) //ntk.fanout_size( n ) == 0u ||
        return false;

      if ( ntk.is_constant(n) )
      {
        return true; /* terminate */
      }

      if ( (ntk.fanin_size(n) == 1) && (ntk.is_pi( ntk.get_children( n, 0 ) ) || ntk.is_constant( ntk.get_children( n, 0 ) )) )
      {
        return true; /* terminate */
      }

      /* verify if there is the need to update the required times */
      if( ps.preserve_depth && ntk.is_marked( n ) )
      {
        compute_required( _max_delay );
        ntk.clear_marked();
      }

      /* find structural optimization opportunities */
      auto str_opto = find_structural_rewriting( cut_manager, cuts, n );

      /* best resub is structural */
      {
        if( str_opto )
        {
          auto index_list = resynthesize_index_list( *str_opto );
//          std::cout << "idlist string = " << to_index_list_string( index_list ) << std::endl;
//          std::cout << "idlist weird  = " ;
//          for( auto x : _database._idlists[(*str_opto).id] )
//            std::cout << x << " ";
//          std::cout << std::endl;

          auto fnew = resynthesize_sub_network( index_list, (*str_opto).permutation, (*str_opto).leaves );
//          std::cout << fnew <<std::endl;
          ntk.substitute_node( n, fnew );

          if( ps.preserve_depth )
          {
            compute_arrival();
          }
          return true;
        }
      }


      auto fun_opto = find_functional_rewriting( n );

      /* best resub is structural */
      {
        if( fun_opto )
        {
          auto index_list = resynthesize_index_list( *fun_opto );
          if( index_list.num_gates() > 0 )
          {
          //std::cout << "F: "<< to_index_list_string( index_list ) << std::endl;
          //for( auto x : _database._idlists[(*fun_opto).id] )
          //  std::cout << x << " ";
          //std::cout << std::endl;

          auto fnew = resynthesize_sub_network( index_list, (*fun_opto).permutation, (*fun_opto).leaves );
          //std::cout << "new signal " << fnew <<std::endl;
          //std::cout << "substitute " << fnew << " " << n <<std::endl;
          ntk.substitute_node( n, fnew );
          std::cout << "substituted " << fnew  <<std::endl;

          if( ps.preserve_depth )
          {
            compute_arrival();
          }
          return true;
          }
        }
      }
      return true;

    });
  }

  #pragma endregion opto

private:
  Ntk& ntk;
  rewrub_sc_params const& ps;
  rewrub_sc_stats& st;

  node_map<float, Ntk> _required;
  node_map<float, Ntk> _arrival;
  float _max_delay;

  pLibrary_t& _database;

  std::array<kitty::static_truth_table<W>,W> _xsW;
  std::array<kitty::static_truth_table<4u>,4u> _xs4;
  unordered_node_map<kitty::static_truth_table<W>, Ntk> _ttW;
  spfd_covering_manager_t<kitty::static_truth_table<W>, 16u> _spfd;
  std::vector<node> _leaves;
  std::vector<node> _divs;
  std::vector<node> _mffc;

  uint32_t _candidates{ 0 };
  uint32_t _estimated_gain{ 0 };


};

} /* namespace detail */

void rewrub_sc( scopt::scg_network& ntk, pLibrary_t& database, rewrub_sc_params const& ps = {}, rewrub_sc_stats* pst = nullptr )
{

  rewrub_sc_stats st;

  using opto_view_t = fanout_view<depth_view<scopt::scg_network>>;
  depth_view<scopt::scg_network> depth_view{ ntk };
  opto_view_t opto_view{ depth_view };

  static constexpr uint32_t W = 16; 

  detail::rewrub_sc_impl<opto_view_t, W> p( opto_view, database, ps, st );
  p.run();
  

  if ( ps.verbose )
  {
    st.report();
  }

  if ( pst )
  {
    *pst = st;
  }

  ntk = cleanup_dangling( ntk );
}


} /* namespace mockturtle */

 