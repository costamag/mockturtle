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
  \file rewrub_mapped.hpp
  \brief Rewrub for mapped networks

  \author Andrea Costamagna
*/

#pragma once

#include "../networks/mapped.hpp"
#include "../utils/stopwatch.hpp"
#include "../utils/spfd_utils.hpp"
#include "../views/depth_view.hpp"
#include "../views/fanout_view.hpp"
#include "cut_enumeration.hpp"
#include "cut_enumeration/rewrite_cut.hpp"
#include "circuit_validator.hpp"


#include <fmt/format.h>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>

namespace mockturtle
{

std::mt19937 RNGRWS( 5 );
bool VERBOSE{false};

/*! \brief Parameters for rewrub.
 *
 * The data structure `rewrub_mapped_params` holds configurable parameters with
 * default arguments for `rewrub`.
 */
struct rewrub_mapped_params
{
  rewrub_mapped_params()
  {
    /* 0 < Cut limit < 16 */
    cut_enumeration_ps.cut_limit = 8;
    cut_enumeration_ps.minimize_truth_table = true;
  }

  /*! \brief Cut enumeration parameters. */
  cut_enumeration_params cut_enumeration_ps{};


  /*! \brief Allow rewrub with multiple structures */
  bool allow_multiple_structures{ true };

  /*! \brief Allow zero-gain substitutions */
  bool allow_zero_gain{ false };

  /*! \brief Use satisfiability don't cares for optimization. */
  bool use_dont_cares{ false };

  /*! \brief Window size for don't cares calculation. */

  /*! \brief Maximum number of PIs of reconvergence-driven cuts. */

  /*! \brief Maximum number of divisors to consider. */
  uint32_t max_divisors{ 256 };

  /*! \brief Maximum number of nodes added by resubstitution. */
  uint32_t max_inserts{ 2 };

  double required_time{ std::numeric_limits<double>::max() };

  /*! \brief Maximum fanout of a node to be considered as root. */
  uint32_t skip_fanout_limit_for_roots{ 1000 };

  /*! \brief Maximum fanout of a node to be considered as divisor. */
  uint32_t skip_fanout_limit_for_divisors{ 100 };

  /*! \brief Number of sampling of the functional cuts. */
  uint32_t num_samplings{1};

  /*! \brief Be verbose. */
  bool verbose{ false };

  double eps_str{0.001};
  double eps_fun{0.001};
  double eps_time{0.001};

  bool try_struct{true};
  bool try_window{true};
  bool try_simula{true};
  bool delay_awareness{true};

  uint32_t max_clauses{ 1000 };
  int32_t odc_levels{ 0 };
  uint32_t conflict_limit{ 1000 };
  uint32_t random_seed{5};

};

/*! \brief Statistics for rewrub.
 *
 * The data structure `rewrub_mapped_stats` provides data collected by running
 * `rewrub`.
 */
struct rewrub_mapped_stats
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

namespace detail
{

template<class Ntk, class Database, uint32_t W, uint32_t S>
class rewrub_mapped_impl
{
  static constexpr uint32_t num_vars = 4u;
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

//  using DivCollector = detail::divisor_collector2_t<Ntk, W>;
//  using collector_st_t = typename DivCollector::stats;
//  using mffc_result_t = double;

  using signature_t = kitty::static_truth_table<S>;
  using word_t = kitty::static_truth_table<6u>;
//  using validator_t = circuit_validator<Ntk, bill::solvers::bsat2, false, true, false>;//last is false


  using network_cuts_t = dynamic_network_cuts<Ntk, num_vars, true, cut_enumeration_rewrite_cut>;
  using cut_manager_t = detail::dynamic_cut_enumeration_impl<Ntk, num_vars, true, cut_enumeration_rewrite_cut>;
  using cut_t = typename network_cuts_t::cut_t;
  using node_data = typename Ntk::storage::element_type::node_type;

public:
  rewrub_mapped_impl( Ntk& ntk, Database& database, rewrub_mapped_params const& ps, rewrub_mapped_stats& st )
      : ntk( ntk ), _database(database), ps( ps ), st( st ), _required( ntk ), _arrival( ntk ), _ttW(ntk), _ttG(ntk), _ttC(ntk)//, _validator( ntk, { ps.max_clauses, ps.odc_levels, ps.conflict_limit, ps.random_seed } )
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


    _is_delay_constrained = ( ps.required_time != std::numeric_limits<double>::max() );

    /* initialize the simulators */
//    _gSim = static_simulator<S>( ntk.num_pis() );
//    _cSim = static_simulator<6>( ntk.num_pis() );
//    simulate_nodes_static<Ntk, S>( ntk, _ttG, _gSim, true );
//    simulate_nodes_static<Ntk, 6>( ntk, _ttC, _cSim, true );

    add_event = ntk.events().register_add_event( [&]( const auto& n ) {
//      _ttG.resize();
//      _ttC.resize();
//      simulate_node_static<Ntk, S>( ntk, n, _ttG, _gSim );
//      simulate_node_static<Ntk, 6>( ntk, n, _ttC, _cSim );
    } );
  }

  ~rewrub_mapped_impl()
  {
    if ( add_event )
    {
      ntk.events().release_add_event( add_event );
    }
  }

  void run()
  {
    stopwatch t( st.time_total );

    perform_rewrub_mapped();

    if( ps.verbose )
    {
      printf("struct %f\n", aStr );
      printf("window %f\n", aWin );
      printf("simula %f\n", aSim );
    }

    printf("  struct %f\n", nStr/NTOT );
    printf("  window %f\n", nWin/NTOT );
    printf("  simula %f\n", nSim/NTOT );
    printf("<#DCS >  %f\n", aDCS );
    printf("#CASES   %f\n", NDCS );

    st.estimated_gain = _estimated_gain;
    st.candidates = _candidates;
  }

  struct opto_candidate_t
  {
    mapped_index_list list;
    std::vector<signal> leaves;
    double reward;
  };

  void perform_rewrub_mapped()
  {  
    /* structural cuts */
    cut_enumeration_stats cst;
    network_cuts_t cuts( ntk.size() + ( ntk.size() >> 1 ) );
    cut_manager_t cut_manager( ntk, ps.cut_enumeration_ps, cst, cuts );
    cut_manager.init_cuts();

    /* window cut */
//    reconvergence_driven_cut_parameters rcuts_ps;
//    rcuts_ps.max_leaves = W;
//    reconvergence_driven_cut_statistics rcuts_st;
//
//    std::array<signal, num_vars> leaves;
//    std::array<signal, num_vars> best_leaves;
//    std::array<uint8_t, num_vars> permutation;
//
//
//    tech_library_params tps;
//    tech_library<5, classification_type::np_configurations> tech_lib( ntk._library, tps );
//    auto [buf_area, buf_delay, buf_id] = tech_lib.get_buffer_info();
//    auto [inv_area, inv_delay, inv_id] = tech_lib.get_inverter_info();
//
//    _BUF_ID = buf_id;
//    
//
//    if ( _is_delay_constrained )
//    {
//      compute_arrival();
//      compute_required( ps.required_time );
//    }
//
    const auto size = ntk.size(); 

    ntk.foreach_gate( [&]( auto const& n, auto i ) {
      
      /* check for trivial cases */
      auto gate = ntk.get_binding( n );
      if( ntk.is_buffer( n ) || ntk.is_constant( n ) || ntk.is_dead( n ) )
      {
        return true;
      } 

      if (  i >= size )
      {
        return false;
      }

      std::vector<opto_candidate_t> candidates;

//      /* verify if there is the need to update the required times */
//      if( _is_delay_constrained && ntk.is_marked( n ) )
//      {
//        compute_arrival();
//        compute_required( ps.required_time );
//        ntk.clear_marked();
//      }
//
//      /* find structural optimization opportunities */
      find_structural_rewriting( candidates, cut_manager, cuts, n );
//      auto win_opto = find_functional_rewriting_exhaustive( n );

      double best_reward = 0;
      int best_index = 0;
      if( candidates.size() > 0 )
      {
        for( int i{0}; i<candidates.size(); ++i )
        {
          if( candidates[i].reward > best_reward )
          {
            printf( "%f\n", candidates[i].reward );
            best_reward = candidates[i].reward;
            best_index = i;
          }
        }
        signal fnew;

        insert( ntk, candidates[best_index].leaves.begin(), candidates[best_index].leaves.end(), candidates[best_index].list, [&]( auto f ){
          fnew = f;
        } );

        ntk.substitute_node( n, fnew );
        return true;
      }


//
//      int choice{-1};
//
//      if( ps.try_struct && str_opto )
//      {
//        bool win_struct = true; 
//        win_struct |= ps.try_window && win_opto && ( (*str_opto).reward > (*win_opto).reward );
//        //bool win_simula = ps.try_simula && ( (*str_opto).reward > (*sim_opto).reward );
//        if( win_struct )//&& win_simula )
//        {
//          choice=0;
//
//        }
//      }
//      
//      if( choice == -1 && ps.try_window && win_opto )
//      {
//        bool win_window = true;
//        win_window |= ps.try_struct && str_opto && ( (*win_opto).reward >= (*str_opto).reward );
//        //bool win_simula = ps.try_simula && ( (*win_opto).reward > (*sim_opto).reward );
//        if( win_window )//&& win_simula )
//          choice=1;
//      }
//
//      /* best resub is structural */
//      if( choice == 0 )
//      {
//        {
//          auto index_list = resynthesize_index_list( *str_opto );
//          auto fnew = resynthesize_sub_network( index_list, (*str_opto).permutation, (*str_opto).leaves );
//          if( index_list.num_gates() == 1 && index_list.ids[0]==buf_id )
//          {
//            fnew = ntk.get_children_signal( ntk.get_node(fnew), 0 );
//          }
//          
//          if( !_is_delay_constrained || is_timing_acceptable( (*str_opto).leaves, fnew, n ) )
//          {
//            nStr+=1;
//            NTOT+=1;
//
//            ntk.substitute_node( n, fnew );
//            if( _is_delay_constrained )
//            {
//              compute_arrival();
//              compute_required( ps.required_time );
//              //ntk.recursively_mark(n);
//            }
//            aStr += (*str_opto).reward;
//          }
//        }
//      }
//      else if( choice == 1 )
//      {
//        {
//          auto index_list = resynthesize_index_list( *win_opto );
//          if( index_list.num_gates() > 0 )
//          {
//            if( VERBOSE )
//            {
//              std::cout << "F: "<< to_index_list_string( index_list ) << std::endl;
//              for( auto x : _database._idlists[(*win_opto).id] )
//                std::cout << x << " ";
//              std::cout << std::endl;
//            }
//
//            auto fnew = resynthesize_sub_network( index_list, (*win_opto).permutation, (*win_opto).leaves );
//
//            if( !_is_delay_constrained || is_timing_acceptable( (*win_opto).leaves, fnew, n ) )
//            {
//              ntk.substitute_node( n, fnew );
//              nWin+=1;
//              NTOT+=1;
//
//              if( _is_delay_constrained )
//              {
//                compute_arrival();
//                compute_required( ps.required_time );
//                //ntk.recursively_mark(n);
//              }
//              aWin += (*win_opto).reward;
//            }
//          }
//        }
//      }
//      else 
//      {
//        auto sim_opto = find_functional_rewriting_signatures( n );
//        if( sim_opto )
//        {
//          auto index_list = resynthesize_index_list( *sim_opto );
//          if( index_list.num_gates() > 0 )
//          {
//            if( VERBOSE )
//            {
//              std::cout << "S: "<< to_index_list_string( index_list ) << std::endl;
//              for( auto x : _database._idlists[(*sim_opto).id] )
//                std::cout << x << " ";
//              std::cout << std::endl;
//            }
//
//            /* check equivalence */
//            std::vector<node> divs;
//            for( int i{0}; i<4; ++i )
//            {
//              divs.push_back( ( *sim_opto ).leaves[i] );
//            }
//
//            auto fnew = resynthesize_sub_network( index_list, (*sim_opto).permutation, (*sim_opto).leaves );
//            
//            if( !_is_delay_constrained || is_timing_acceptable( (*sim_opto).leaves, fnew, n ) )
//            {
//              auto valid = _validator.validate( ntk.make_signal(n), fnew );
//              if( valid )
//              {
//                if( *valid )
//                { 
//                  nSim+=1;
//                  NTOT+=1;
//                  ntk.substitute_node( n, fnew );
//                  aSim += (*sim_opto).reward;
//
//                  if( _is_delay_constrained )
//                  {
//                    compute_arrival();
//                    compute_required( ps.required_time );
//                    //ntk.recursively_mark(n);
//
//                  }
//                }
//                else
//                {
//                  found_cex();
//                }
//              }
//            }
//          }
//        }
//        //printf("simula has a candidate\n");
//      }
//      
      return true;

    });
  }

  void find_structural_rewriting( std::vector<opto_candidate_t>& candidates, cut_manager_t & cut_manager, network_cuts_t & cuts, node const& n )
  {
    cut_manager.clear_cuts( n );
    cut_manager.compute_cuts( n );

    uint32_t cut_index = 0;
    double best_reward = -1;

    for ( auto& cut : cuts.cuts( ntk.node_to_index( n ) ) )
    {
      /* skip trivial cut */
      if ( ( cut->size() == 1 && *cut->begin() == ntk.node_to_index( n ) ) )
      {
        ++cut_index;
        continue;
      }

      /* Boolean matching */
      kitty::dynamic_truth_table tt( 4 );
      kitty::create_from_binary_string( tt, kitty::to_binary( cuts.truth_table( *cut ) ) );

      std::vector<node> nodes;

      for ( auto const leaf : *cut )
      {
        nodes.push_back( ntk.index_to_node( leaf ) );
      }

      while ( nodes.size() < tt.num_vars() )
      {
        nodes.push_back( ntk.get_node( ntk.get_constant( false ) ) );
      }

      auto match = _database.get_matches( tt, nodes );

      if( match )
      {
        auto ordered_nodes = (*match).first;
        _database.foreach_entry( (*match).second, [&]( auto entry ){      
          double cut_area = area_contained( n, ordered_nodes );
          if( entry.area < cut_area )
          {
            opto_candidate_t cand;
            cand.leaves = ordered_nodes;
            cand.list = entry.list;
            cand.reward = cut_area - entry.area;
            candidates.push_back( cand );
          }
        } );
      }
    }
  }

  double area_contained_rec( node n )
  {
    if( ntk.visited( n ) == ntk.trav_id() )
    {
      return 0;
    }

    double area{ 0.0 };
    ntk.foreach_fanin( n, [&]( auto f ){
      auto ni = ntk.get_node( f );
      area += area_contained_rec( ni );
    } );
    ntk.set_visited( n, ntk.trav_id() );
    area += ntk.get_binding( n ).area;
    return area;
  }

  double area_contained( node n, std::vector<node> const& leaves )
  {
    ntk.incr_trav_id();
    for( auto leaf : leaves )
    {
      ntk.set_visited( leaf, ntk.trav_id() );
    }

    /* measure the MFFC contained in the cut */
    return area_contained_rec( n );
  } 

//
//  #pragma region timing 
//
//  double compute_arrival_rec( node n )
//  {
//    if( _arrival.has(n) && ( ntk.visited( n ) == 1u ) )
//    {
//      return _arrival[n];
//    }
//
//    auto const& g = ntk.get_binding( n );
//    double arrival = 0;
//
//    ntk.foreach_fanin( n, [&]( auto const& f, auto i ) {
//      double arr_fanin = compute_arrival_rec( ntk.get_node(f) ); // normalized by default
//      double pin_delay = std::max( g.pins[i].rise_block_delay, g.pins[i].fall_block_delay );
//      pin_delay = std::ceil(pin_delay * 100.0) / 100.0;
//      arrival = std::max( arrival, (double)( arr_fanin + pin_delay ) );
//    } );
//    _arrival[n] = std::ceil( arrival * 100.0) / 100.0;
//    ntk.set_visited( n, 1u );
//
//    return _arrival[n];
//  }
//
//  double compute_arrival()
//  {
//    ntk.clear_visited();
//    _arrival.reset();
//
//    ntk.foreach_pi( [&]( auto const& n, auto i ) {
//      _arrival[n]=0.0;
//      ntk.set_visited( n, 1u );
//    } );
//    _arrival[0]=0.0;
//    ntk.set_visited( 0, 1u );
//
//    double max_delay=0;
//    ntk.foreach_po( [&]( auto const& fo, auto i ) {
//      node no = ntk.get_node(fo);
//      double out_del = compute_arrival_rec( no );
//      if( out_del > max_delay )
//      {
//        max_delay = out_del;
//      }
//    } );
//
//    ntk.foreach_gate( [&]( auto const& n, auto i ) { 
//      compute_arrival_rec( n );
//    });
//    
//    ntk.clear_visited();
//    return max_delay;
//  }
//
//  double compute_required_rec( node n, double max_delay )
//  {
//    if( _required.has(n) && ( ntk.visited(n) == 1u ) )
//    {
//      return _required[n];
//    }
//
//    double gate_required = max_delay;
//
//    ntk.foreach_fanout( n, [&]( auto const& f, auto i ) {
//      auto nfo = ntk.get_node(f);
//      double req_fanout = compute_required_rec( nfo, max_delay );
//      auto const& g = ntk.get_binding( nfo );
//
//      uint32_t ig;
//      ntk.foreach_fanin( nfo, [&]( auto const& fi, auto ii ) {
//        if( ntk.get_node(fi) == n )
//        {
//          ig = ii;
//          return;
//        }
//      });
//
//      gate_required = std::min( gate_required, (double)( req_fanout - std::max( g.pins[ig].rise_block_delay, g.pins[ig].fall_block_delay ) ) );
//    } );
//    _required[n] = std::ceil(gate_required * 100.0) / 100.0;
// 
//    ntk.set_visited( n, 1u );
//
//    return _required[n];
//  }
//
//  void compute_required( double max_delay )
//  {
//    ntk.clear_visited();
//    _required.reset();
//
//    ntk.foreach_po( [&]( auto const& fo, auto i ) 
//    {
//      node no = ntk.get_node(fo);
//      if( ntk.fanout_size(no) == 1 )
//      {
//        _required[no]= std::ceil(max_delay * 100.0) / 100.0;
//        ntk.set_visited( no, 1u );
//      }
//    });
//
//    ntk.foreach_pi( [&]( auto const& ni, auto i ) {
//      auto req = compute_required_rec( ntk.get_node(ni), max_delay );
//    } );
//    
//    ntk.clear_visited();
//
//  }
//
//  void print_slack()
//  {
//    ntk.foreach_gate( [&]( auto const& n, auto i ) { 
//      if( ntk.po_index(n) != -1 )
//        printf("po %4d a=%f r=%f s=%f\n", n, _arrival[n], _required[n], _required[n]-_arrival[n]);
//      else if( ntk.is_pi(n) )
//        printf("pi %4d a=%f r=%f s=%f\n", n, _arrival[n], _required[n], _required[n]-_arrival[n]);
//      else
//        printf("nd %4d a=%f r=%f s=%f\n", n, _arrival[n], _required[n], _required[n]-_arrival[n]);
//
//
//    });
//  }
//
//  #pragma endregion timing
//
//  #pragma region opto
//  struct opto_candidate_t
//  {
//    uint32_t id;
//    std::array<signal, num_vars> leaves;
//    std::array<uint8_t, num_vars> permutation;
//    double reward;
//  };
//
//  double measure_mffc_ref( node const& n, std::array<signal,num_vars> const& cut )
//  {
//    /* reference cut leaves */
//    for ( auto leaf : cut )
//    {
//      ntk.incr_fanout_size( ntk.get_node(leaf) );
//    }
//
//    double mffc_size = static_cast<double>( recursive_ref( n ) );
//
//    /* dereference leaves */
//    for ( auto leaf : cut )
//    {
//      ntk.decr_fanout_size( ntk.get_node(leaf) );
//    }
//
//    return std::ceil(mffc_size*100.0)/100.0;
//  }
//
//  double measure_mffc_deref( node const& n, std::array<signal,num_vars> const& cut )
//  {
//    /* reference cut leaves */
//    for ( auto leaf : cut )
//    {
//      ntk.incr_fanout_size( ntk.get_node(leaf) );
//    }
//
//    double mffc_size = static_cast<double>( recursive_deref( n ) );
//
//    /* dereference leaves */
//    for ( auto leaf : cut )
//    {
//      ntk.decr_fanout_size( ntk.get_node(leaf) );
//    }
//
//    return std::ceil(mffc_size*100.0)/100.0;
//  }
//
//  double recursive_deref( node const& n )
//  {
//    /* terminate? */
//    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
//      return 0;
//
//    /* recursively collect nodes */
//    double value{ ntk.get_area( n ) };
//    ntk.foreach_fanin( n, [&]( auto const& s ) {
//      if ( ntk.decr_fanout_size( ntk.get_node( s ) ) == 0 )
//      {
//        value += recursive_deref( ntk.get_node( s ) );
//      }
//    } );
//    return std::ceil(value*100.0)/100.0;
//  }
//
//  double recursive_ref( node const& n )
//  {
//    /* terminate? */
//    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
//      return 0;
//
//    /* recursively collect nodes */
//    double value{ ntk.get_area( n ) };
//    ntk.foreach_fanin( n, [&]( auto const& s ) {
//      if ( ntk.incr_fanout_size( ntk.get_node( s ) ) == 0 )
//      {
//        value += recursive_ref( ntk.get_node( s ) );
//      }
//    } );
//    return std::ceil(value*100.0)/100.0;
//  }
//
//  double area_contained_mffc( node n, std::array<signal, num_vars> const& leaves )
//  {
//    /* measure the MFFC contained in the cut */
//    double mffc_size = measure_mffc_deref( n, leaves );
//    /* restore contained MFFC */
//    measure_mffc_ref( n, leaves );
//
//    return std::ceil(mffc_size*100.0)/100.0;
//  }  
//
//  void print_cone( node n, std::array<signal, num_vars> const&  leaves, std::set<node> & nodes )
//  {
//    if( std::find( leaves.begin(), leaves.end(), ntk.make_signal(n) ) == leaves.end() )
//    {
//      if( nodes.find( n ) == nodes.end() )
//      {
//        std::vector<node> children;
//        ntk.foreach_fanin( n, [&]( auto const& f ) {
//          children.push_back( ntk.get_node(f) );
//          print_cone( ntk.get_node(f), leaves, nodes );
//        } );
//        nodes.insert(n);
//        printf("n %d - id %d :", n, ntk.get_binding(n).id );
//        for( auto ch : children )
//        {
//          printf(" %d", ch );
//        }
//        printf("\n");
//      }
//
//    }
//    else
//    {
//      if( nodes.find( n ) == nodes.end() )
//      {
//        std::vector<node> children;
//        ntk.foreach_fanin( n, [&]( auto const& f ) {
//          children.push_back( ntk.get_node(f) );
//        } );
//        nodes.insert(n);
//        printf("n %d - id %d :", n, ntk.get_binding(n).id );
//        for( auto ch : children )
//        {
//          printf(" %d", ch );
//        }
//        printf("\n");
//      }
//    }
//    return;
//  }
//

//
//  void simulate_window( std::vector<node> const& leaves, std::vector<node> const& divs, std::vector<node> const& mffc, node n )
//  {
//    _ttW.reset();
//    _ttW[0]=_xsW[0].construct();
//    int i{0};
//    for( auto l : leaves )
//    {
//      _ttW[l] = _xsW[i++];
//      if( VERBOSE )
//      {
//        printf("[l %3d]",l); kitty::print_binary(_ttW[l]); 
//      }
//
//    }
//    if(VERBOSE)
//        printf("\n");
//
//    std::vector<kitty::static_truth_table<W>> children; 
//    for( auto d : divs )
//    {
//      if( !ntk.is_constant(d) && std::find( leaves.begin(), leaves.end(), d ) == leaves.end() )
//      {
//        children.clear();
//        ntk.foreach_fanin( d, [&]( auto const& f ) {
//          children.push_back( _ttW[ntk.get_node(f)] );
//        } );
//        _ttW[d] = ntk.compute( d, children.begin(), children.end() );
//        if( VERBOSE )
//        {
//          printf("d %3d:",d); kitty::print_binary(_ttW[d]); 
//          ntk.foreach_fanin( d, [&]( auto const& f ) {
//            std::cout << " " << ntk.get_node(f);
//          } );
//          std::cout << " id" << ntk.get_binding(d).id;
//          printf("\n");
//        }
//      }
//    }
//
//    for( auto d : mffc )
//    {
//      if( !ntk.is_constant(d) && std::find( leaves.begin(), leaves.end(), d ) == leaves.end() )
//      {
//        children.clear();
//        ntk.foreach_fanin( d, [&]( auto const& f ) {
//          children.push_back( _ttW[ntk.get_node(f)] );
//        } );
//
//        _ttW[d] = ntk.compute( d, children.begin(), children.end() );
//
//          if( VERBOSE )
//          {
//            printf("m %3d:",d); kitty::print_binary(_ttW[d]); 
//            ntk.foreach_fanin( d, [&]( auto const& f ) {
//              std::cout << " " << ntk.get_node(f);
//            } );
//            std::cout << " id" << ntk.get_binding(d).id;
//
//            printf("\n");
//          }
//
//      }
//    }
//    children.clear();
//    ntk.foreach_fanin( n, [&]( auto const& f ) {
//        children.push_back( _ttW[ntk.get_node(f)] );
//    } );
//
//    _ttW[n] = ntk.compute( n, children.begin(), children.end() );
//    if( VERBOSE )
//    {
//      printf("n %3d:",n); kitty::print_binary(_ttW[n]); printf("\n");
//    }
//
//  }
//
//  void print_window( std::vector<node> const& leaves, std::vector<node> const& divs, std::vector<node> const& mffc, node n )
//  {
//
//    int i{0};
//    for( auto l : leaves )
//    {
//      printf("[l %3d]",l); //kitty::print_binary(_ttW[l]); 
//    }
//    printf("\n");
//
//    std::vector<kitty::static_truth_table<W>> children; 
//    for( auto d : divs )
//    {
//      if( !ntk.is_constant(d) && std::find( leaves.begin(), leaves.end(), d ) == leaves.end() )
//      {
//          printf("d %3d id %d:",d, ntk.get_binding(d).id); 
//          //kitty::print_binary(_ttW[d]); 
//          ntk.foreach_fanin( d, [&]( auto const& f ) {
//            std::cout << " " << ntk.get_node(f);
//          } );
//          std::cout << " id" << ntk.get_binding(d).id;
//          printf("\n");
//      }
//    }
//
//    for( auto d : mffc )
//    {
//      if( !ntk.is_constant(d) && std::find( leaves.begin(), leaves.end(), d ) == leaves.end() )
//      {
//
//            printf("m %3d id %d:",d, ntk.get_binding(d).id ); 
//            //kitty::print_binary(_ttW[d]); 
//            ntk.foreach_fanin( d, [&]( auto const& f ) {
//              std::cout << " " << ntk.get_node(f);
//            } );
//            std::cout << " id" << ntk.get_binding(d).id;
//
//            printf("\n");
//
//      }
//    }
//
//
//      printf("n %3d id %d:",n, ntk.get_binding(n).id ); kitty::print_binary(_ttW[n]); printf("\n");
//
//  }
//
//  template<class SIM, class SPFD>
//  std::optional<std::vector<node>> find_support_greedy( std::vector<node> const& divs, SIM const& tts, SPFD & spfd )
//  {
//    uint32_t cost, best_cost;
//    std::vector<node> best_candidates;
//    std::vector<node> supp;
//    spfd.reset();
//
//    /* add recomputation of the support */
//    while ( !spfd.is_covered() && supp.size() < 4u )
//    {
//      best_cost = std::numeric_limits<uint32_t>::max();
//      if ( spfd.is_saturated() )
//        return std::nullopt;
//      for ( uint32_t v{0}; v < divs.size(); ++v )
//      {
//        cost = spfd.evaluate( tts[divs[v]] );
//
//        if ( cost < best_cost )
//        {
//          best_cost = cost;
//          best_candidates = { divs[v] };
//        }
//        else if ( cost == best_cost )
//        {
//          best_candidates.push_back( divs[v] );
//        }
//      }
//      if ( best_candidates.size() == 0 )
//        return std::nullopt;
//
//      std::uniform_int_distribution<> distrib( 0, best_candidates.size() - 1 );
//      int idx = distrib( RNGRWS );
//      supp.push_back( best_candidates[idx] );
//      spfd.update( tts[best_candidates[idx]] );
//    }
//
//    if ( spfd.is_covered() && supp.size() <= 4u )
//    {
//      std::sort( supp.begin(), supp.end() );
//      return supp;
//    }
//    return std::nullopt;
//  }
//
//  template<class SIM, class SPFD>
//  std::optional<std::vector<node>> find_support_greedy_delay( std::vector<node> const& divs, SIM const& tts, SPFD & spfd )
//  {
//    uint32_t cost, best_cost;
//    std::vector<node> best_candidates;
//    std::vector<node> supp;
//    spfd.reset();
//
//    /* add recomputation of the support */
//    while ( !spfd.is_covered() && supp.size() < 4u )
//    {
//      double best_delay = std::numeric_limits<double>::max();
//
//      best_cost = std::numeric_limits<uint32_t>::max();
//      if ( spfd.is_saturated() )
//        return std::nullopt;
//      for ( uint32_t v{0}; v < divs.size(); ++v )
//      {
//        cost = spfd.evaluate( tts[divs[v]] );
//
//        if ( cost < best_cost || ( cost == best_cost && _arrival[divs[v]]<best_delay ) )
//        {
//          best_cost = cost;
//          best_delay = _arrival[divs[v]];
//          best_candidates = { divs[v] };
//        }
//        else if ( cost == best_cost && _arrival[divs[v]] == best_delay )
//        {
//          best_candidates.push_back( divs[v] );
//        }
//      }
//      if ( best_candidates.size() == 0 )
//        return std::nullopt;
//
//      std::uniform_int_distribution<> distrib( 0, best_candidates.size() - 1 );
//      int idx = distrib( RNGRWS );
//      supp.push_back( best_candidates[idx] );
//      spfd.update( tts[best_candidates[idx]] );
//    }
//
//    if ( spfd.is_covered() && supp.size() <= 4u )
//    {
//      std::sort( supp.begin(), supp.end() );
//      return supp;
//    }
//    return std::nullopt;
//  }
//
//  template<class SIM, class SPFD>
//  std::optional<std::vector<node>> find_support( std::vector<node> const& divs, node n, SIM const& tts, SPFD & spfd )
//  {
//    spfd.init( tts[n] );
//    auto supp = ps.delay_awareness ? find_support_greedy_delay( divs, tts, spfd ) : find_support_greedy( divs, tts, spfd );
//    return supp;
//  }
//
//  template<class SIM>
//  std::tuple<kitty::static_truth_table<4u>, kitty::static_truth_table<4u>> extract_functionality( std::vector<node> const & cut, node n, SIM const& tts )
//  {
//    kitty::static_truth_table<4u> tt;
//    kitty::static_truth_table<4u> mk;
//
//    auto tmp = tts[n].construct();
//    
//    for( uint32_t m{0}; m < ( 1u << cut.size() ); ++m )
//    {
//      auto tmp4 = _xs4[0]|~_xs4[0];
//      tmp = tmp | ~tmp;
//      for( int i{0}; i<cut.size(); ++i )
//      {
//        if( ( ( m >> i ) & 0x1 ) == 0x1 )
//        {
//          tmp &= tts[cut[i]];
//          tmp4 &= _xs4[i];
//        }
//        else
//        {
//          tmp &= ~tts[cut[i]];
//          tmp4 &= ~_xs4[i];
//        }  
//      }
//      int n0 = kitty::count_ones( ~tts[n] & tmp );
//      int n1 = kitty::count_ones(  tts[n] & tmp );
//      if( n0 > 0 && n1 == 0 )
//      {
//        mk |= tmp4;
//      }
//      else if( n1 > 0 && n0 == 0 )
//      {
//        tt |= tmp4;
//        mk |= tmp4;
//      }
//      else if( n1 != 0 || n0 != 0 )
//      {
//        printf("not a valid support\n");
//        assert( false && "Not valid support\n" );
//      }
//    }
//    return std::make_tuple( tt, mk );
//  }
//
//  std::optional<opto_candidate_t> find_functional_rewriting_exhaustive( node const& n )
//  {
//
//    if( !ps.try_window )
//    {
//      return std::nullopt;
//    }
//    
//    double mffc_area;
//    collector_st_t collector_st;
//
//    DivCollector collector( ntk, _arrival, _required, ps, collector_st );
//    const auto collector_success = collector.run( n, mffc_area );
//    if ( !collector_success )
//    {
//      return std::nullopt; /* next */
//    }
//    else if( ps.required_time != std::numeric_limits<double>::max() )
//    {
//
//    }
//
//    std::vector<signal> leaves_sig(collector.leaves.size());
//    std::transform( collector.leaves.begin(), collector.leaves.end(), leaves_sig.begin(), [&]( const node n ) {
//      auto f = ntk.make_signal( n );
//      return f;
//    } );
//
//    simulate_window( collector.leaves, collector.divs, collector.mffc, n );
//
//    // find functional cut
//    auto supp = find_support( collector.divs, n, _ttW, _wSpfd );
//    if( supp )
//    {
//      if( VERBOSE )
//      {
//        std::cout << "SUPP|w:";
//        for( auto x : *supp )
//        {
//          std::cout << x << " ";
//        }
//        std::cout << std::endl;
//      }
//
//      auto [func, care] = extract_functionality( *supp, n, _ttW );
//      auto dontcare = ~care;
//
//      std::vector<uint32_t> dcs;
//      for( int bit{0}; bit<16; ++bit )
//      {
//        if( kitty::get_bit( dontcare, bit ) > 0 )
//        {
//          dcs.push_back(bit);
//        }
//      }
//      aDCS = aDCS*(NDCS-1)/NDCS+((double)dcs.size())/NDCS;
//      NDCS+=1;
//      
//
//      uint64_t best_key;
//      double best_area = mffc_area;
//      std::vector<uint8_t> best_perm;
//
//
//      for( uint32_t m{0}; m< (1<<dcs.size()); ++m )
//      {
//        auto tt = func;
//
//        for( int i{0}; i<dcs.size(); ++i )
//        {
//          if( (m >> i)&0x1 == 0x1 )
//          {
//            kitty::flip_bit( tt, dcs[i] );
//          }
//        }
//
//        /* p-canonize */
//        auto config = kitty::exact_p_canonization( tt );
//
//        auto repr = std::get<0>( config );
//        auto neg = std::get<1>( config );
//        auto perm = std::get<2>( config );
//
//        auto key = _database.get_key( repr );
//
//        if( key )
//        {
//          double area = _database._areas[*key];
//          if( *key == _BUF_ID )
//            area = 0;
//          if( area < best_area )
//          {
//            best_key = *key;
//            best_area = _database._areas[*key];
//            best_perm = perm;
//          }
//        }
//      }
//
//      if( mffc_area - best_area > ps.eps_fun )
//      {
//
//        print_window( collector.leaves, collector.divs, collector.mffc, n );
//
//        std::cout << "rwq"<<  (mffc_area - best_area) << "{}" << best_key << "^SUPP|w:";
//        for( auto x : *supp )
//        {
//          std::cout << x << " ";
//        }
//        std::cout << std::endl;
//        std::cout << std::endl;
//
//
//
//        opto_candidate_t cand;
//
//        std::array<uint8_t, num_vars> permutation;
//
//        for ( auto j = 0u; j < num_vars; ++j )
//        {
//          permutation[best_perm[j]] = j;
//        }
//
//        /* save output negation to apply */
//
//        {
//          auto j = 0u;
//          for ( auto const leaf : *supp )
//          {
//            cand.leaves[permutation[j++]] = ntk.make_signal( leaf );
//          }
//
//          while ( j < num_vars )
//            cand.leaves[permutation[j++]] = ntk.get_constant( false );
//        }
//
//        /* resynthesis cost */
//        cand.id = best_key;
//
//        /* resynthesis reward */
//        double area_mffc = area_contained_mffc( n, cand.leaves );
//
//        cand.reward = area_mffc - best_area;
//        return cand;
//      }
//
//      }
//
//    return std::nullopt;
//  }
//
//  void check_tts( node const& n )
//  {
//    if( !_ttG.has(n) )
//    {
//      _ttG.resize();
//      _ttC.resize();
//      simulate_node_static<Ntk, S>( ntk, n, _ttG, _gSim );
//      simulate_node_static<Ntk, 6>( ntk, n, _ttC, _cSim );
//    }
//    else if ( _ttG[n].num_bits() != _gSim.num_bits() )
//    {
//      simulate_node_static<Ntk, S>( ntk, n, _ttG, _gSim );
//      simulate_node_static<Ntk, 6>( ntk, n, _ttC, _cSim );
//    }
//  }
//
//  std::optional<opto_candidate_t> find_functional_rewriting_signatures( node const& n )
//  {
//
//    if( !ps.try_simula )
//    {
//      return std::nullopt;
//    }
//
//    double mffc_area;
//    collector_st_t collector_st;
//
//    DivCollector collector( ntk, _arrival, _required, ps, collector_st );
//    const auto collector_success = collector.run( n, mffc_area );
//    if ( !collector_success )
//    {
//      return std::nullopt; /* next */
//    }
//
//    std::vector<signal> leaves_sig(collector.leaves.size());
//    std::transform( collector.leaves.begin(), collector.leaves.end(), leaves_sig.begin(), [&]( const node n ) {
//      auto f = ntk.make_signal( n );
//      return f;
//    } );
//
//    // verify that all the signatures are valid
//    check_tts(n);
//    for( auto d : collector.divs )
//    {
//      check_tts(d);
//    }
//
//    // find functional cut
//    auto supp = find_support( collector.divs, n, _ttG, _gSpfd );
//    if( supp )
//    {
//      if( VERBOSE )
//      {
//        std::cout << "SUPP|s:";
//        for( auto x : *supp )
//        {
//          std::cout << x << " ";
//        }
//        std::cout << std::endl;
//      }
//      std::nullopt;
//
//      auto [func, care] = extract_functionality( *supp, n, _ttG );
//      auto dontcare = ~care;
//
//      std::vector<uint32_t> dcs;
//      for( int bit{0}; bit<16; ++bit )
//      {
//        if( kitty::get_bit( dontcare, bit ) > 0 )
//        {
//          dcs.push_back(bit);
//        }
//      }
//      aDCS = aDCS*(NDCS-1)/NDCS+((double)dcs.size())/NDCS;
//      NDCS+=1;
//
//      uint64_t best_key;
//      double best_area = mffc_area;
//      std::vector<uint8_t> best_perm;
//
//
//      for( uint32_t m{0}; m< (1<<dcs.size()); ++m )
//      {
//        auto tt = func;
//
//        for( int i{0}; i<dcs.size(); ++i )
//        {
//          if( (m >> i)&0x1 == 0x1 )
//          {
//            kitty::flip_bit( tt, dcs[i] );
//          }
//        }
//
//        /* p-canonize */
//        auto config = kitty::exact_p_canonization( tt );
//
//        auto repr = std::get<0>( config );
//        auto neg = std::get<1>( config );
//        auto perm = std::get<2>( config );
//
//        auto key = _database.get_key( repr );
//
//        if( key )
//        {
//          double area = _database._areas[*key];
//          if( *key == _BUF_ID )
//            area = 0;
//          if( area < best_area )
//          {
//            best_key = *key;
//            best_area = _database._areas[*key];
//            best_perm = perm;
//          }
//        }
//      }
//
//      if( mffc_area - best_area > ps.eps_fun )
//      {
//        opto_candidate_t cand;
//
//        std::array<uint8_t, num_vars> permutation;
//
//        for ( auto j = 0u; j < num_vars; ++j )
//        {
//          permutation[best_perm[j]] = j;
//        }
//
//        /* save output negation to apply */
//
//        {
//          auto j = 0u;
//          for ( auto const leaf : *supp )
//          {
//            cand.leaves[permutation[j++]] = ntk.make_signal( leaf );
//          }
//
//          while ( j < num_vars )
//            cand.leaves[permutation[j++]] = ntk.get_constant( false );
//        }
//
//        /* resynthesis cost */
//        cand.id = best_key;
//
//        /* resynthesis reward */
//        double area_mffc = area_contained_mffc( n, cand.leaves );
//
//        cand.reward = area_mffc - best_area;
//        return cand;
//      }
//
//      }
//
//    return std::nullopt;
//  }
//
//  large_lig_index_list resynthesize_index_list( opto_candidate_t const& cand )
//  {
//    large_lig_index_list index_list(4u);
//
//    std::vector<uint32_t> lits = {0,2,4,6,8};
//
//    auto entry = _database._idlists[ cand.id ];
//    int type = 0;
//    int nFins = 0;
//    uint32_t sc_id;
//    std::vector<uint32_t> children;
//    uint32_t lit;
//    
//    for( int i{0}; i<entry.size(); ++i )
//    {
//      if( type == 0 )
//      {
//        nFins = entry[i];
//        type=1;
//      }
//      else if( type == 1 )
//      {
//        children.push_back( lits[entry[i]] );//not accounting for 0
//        if( children.size() == nFins )
//        {
//          type = 2;
//        }
//      }
//      else if( type == 2 )
//      {
//        type = 0;
//        sc_id = entry[i];
//        lit = index_list.add_function( children, ntk._library[sc_id].function, ntk._library[sc_id].area, ntk._library[sc_id].id );
//        lits.push_back( lit );
//        children.clear();
//      }
//    }
//    index_list.add_output( lit );
//    return index_list;
//  }
//
//  signal resynthesize_sub_network( large_lig_index_list const& index_list, std::array<uint8_t,4> const& perm, std::array<signal, 4> const& leaves )
//  {
//    std::vector<signal> divs_sig;
//    for( int i{0}; i<4u; ++i ) 
//    {
//      divs_sig.push_back( leaves[i] );
//    }
//
//    signal res;
//    insert( ntk, divs_sig.begin(), divs_sig.end(), index_list, [&]( signal const& s ) {
//      res = s;
//    } );
//    return res;
//  }
//
//  bool is_timing_acceptable( std::array<signal, 4u> const& leaves, signal fnew, node nold )
//  {
//
//    node nnew = ntk.get_node( fnew );
//
//    /* necessary setup for evaluating arrival time */
//    ntk.clear_visited();
//    for( auto x : leaves )
//    {
//      ntk.set_visited( ntk.get_node(x), 1u );
//    }
//    double new_arrival = compute_arrival_rec( nnew );
//    double new_required = _required.has(nnew) ? _required[nnew] : std::numeric_limits<double>::max();
//
//    ntk.clear_visited();
//
//    if( ( new_arrival < (_required[nold]-ps.eps_time) ) && ( new_arrival < (new_required-ps.eps_time) ) && (new_arrival < (ps.required_time-ps.eps_time) ) )// && (new_arrival < ps.required_time ) && new_required  _required[nold] )
//    {
//      return true;
//    }
//
//    return false;
//  }
//
//  void found_cex()
//  {
//    sig_pointer = (sig_pointer+1)%(1<<S);
//
//    _cSim.add_pattern( _validator.cex );
//    if ( sig_pointer % 64 == 0 )
//    {
//      _ttC.reset();
//      simulate_nodes_static<Ntk>( ntk, _ttC, _cSim, true );
//
//      ntk.foreach_pi( [&]( auto const& n, auto i ) {
//        *(_ttG[n].begin() + _block) = *(_ttC[n].begin());
//      } );
//
//      ntk.foreach_gate( [&]( auto const& n, auto i ) {
//        *(_ttG[n].begin() + _block) = *(_ttC[n].begin());
//      } );
//
//      _block = S == 6u ? 0u : ( _block + 1u ) % ( ( 1u << ( S - 6u ) ) - 1u ) ;
//    }
//  }
//

//
//  #pragma endregion opto
//
private:
  Ntk& ntk;
  rewrub_mapped_params const& ps;
  rewrub_mapped_stats& st;
  Database _database;

  unordered_node_map<double, Ntk> _required;
  unordered_node_map<double, Ntk> _arrival;
  double _BUF_AREA{0};
  int _BUF_ID{0};

  std::array<kitty::static_truth_table<W>,W> _xsW;
  std::array<kitty::static_truth_table<4u>,4u> _xs4;
  unordered_node_map<kitty::static_truth_table<W>, Ntk> _ttW;

  incomplete_node_map<signature_t, Ntk> _ttG;
  incomplete_node_map<word_t, Ntk> _ttC;
  //validator_t _validator;

  //static_simulator<S> _gSim;
  //static_simulator<6> _cSim;
  uint32_t _block{0};
  double nStr{0};
  double nSim{0};
  double nWin{0};
  double aStr{0};
  double aSim{0};
  double aWin{0};
  double aDCS{0};
  double NTOT{0};
  double NDCS{1};

  uint32_t sig_pointer{0};


  spfd_covering_manager_t<kitty::static_truth_table<W>, 16u> _wSpfd;
  spfd_covering_manager_t<kitty::static_truth_table<S>, 16u> _gSpfd;
  std::vector<node> _leaves;
  std::vector<node> _divs;
  std::vector<node> _mffc;

  uint32_t _candidates{ 0 };
  uint32_t _estimated_gain{ 0 };

  bool _is_delay_constrained{false};

  std::shared_ptr<typename network_events<Ntk>::add_event_type> add_event;


};


} /* namespace detail */

template<class Database, uint32_t W=10u, uint32_t S=10u>
void rewrub_mapped( mapped_network& ntk, Database& database, rewrub_mapped_params const& ps = {}, rewrub_mapped_stats* pst = nullptr )
{
  rewrub_mapped_stats st;

  using opto_view_t = depth_view<mapped_network>;
  depth_view<mapped_network> depth_view{ ntk };

  detail::rewrub_mapped_impl<opto_view_t, Database, W, S> p( depth_view, database, ps, st );
  p.run();
  
  if ( ps.verbose )
  {
    st.report();
  }

  if ( pst )
  {
    *pst = st;
  }

  ntk = cleanup_mapped( ntk );
}

} /* namespace mockturtle */

 