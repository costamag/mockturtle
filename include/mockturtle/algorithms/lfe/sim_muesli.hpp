/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2021  EPFL
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
  \file sim_muesli.hpp
  \brief sim_muesli algorithm

  \author Andrea Costamagna
*/

#pragma once

#include <cstdint>
#include <vector>

#include "../../traits.hpp"
#include "simulation_view.hpp"
#include "sim_patterns.hpp"
#include "sim_decomposition_fastS.hpp"


#include <fmt/format.h>
#include <kitty/constructors.hpp>
#include <kitty/decomposition.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/statistics.hpp>
#include <kitty/operators.hpp>
#include <kitty/print.hpp>

#include "create_candidates.hpp"
#include "chatterjee_method.hpp"


namespace mockturtle
{
/*! \brief Parameters for sim_muesli algorithm */
struct sim_muesli_params
{
  size_t init_sup{2};
  size_t max_sup{3};
  size_t max_act{3};
  double eps_th{1};
  bool verbose{false};
  bool try_accuracy_recovery{false};
  bool use_min_mi{false};
  bool one_by_one{true};
  bool re_initialize{true};
};

namespace detail
{

  template<typename Ntk>
  class sim_muesli_impl
  {
    public:
      using TT = kitty::partial_truth_table;

    public:
      sim_muesli_impl( simulation_view<Ntk> & ntk, TT target, sim_muesli_params & ps )
      : ntk(ntk),
        ps(ps),
        target(target)
      {
        active_list = std::make_pair( std::vector{ &ntk.sim_patterns[0].pat }, std::vector{ &ntk.sim_patterns[0] } );
      }
    
      void fill_active_list( size_t const& act )
      {
        assert( act > 0 );
        active_list.first.clear();
        active_list.second.clear();
        std::pair<std::vector<TT*>, std::vector<sim_pattern<Ntk> * >> active_list_tmp;

        Iactive.clear();
        ntk.clear_flag();


        for( size_t i = 0; i < std::min( act, (size_t)ntk.sim_patterns.size()); ++i )
        {

          active_list.first.emplace_back( &(ntk.sim_patterns[0].pat) );
          active_list.second.emplace_back( &ntk.sim_patterns[0]);
          active_list_tmp.first.emplace_back( &(ntk.sim_patterns[0].pat) );
          active_list_tmp.second.emplace_back( &ntk.sim_patterns[0]);

          Iactive.emplace_back( (double)0 );

          double Imax{.0};
          double Inew{.0};

          for( size_t j = 0; j < ntk.sim_patterns.size(); ++j )
          { 
            if( !ntk.sim_patterns[j].flag )
            {
              active_list_tmp.first[i] = &(ntk.sim_patterns[j].pat);
              active_list_tmp.second[i] = &(ntk.sim_patterns[j]);

              if( i != 0 )
                Inew = kitty::mutual_information( active_list_tmp.first, &target );
              else
              {
                if( ntk.sim_patterns[j].weight < 0 )
                {
                  Inew = kitty::mutual_information( active_list_tmp.first, &target );
                  ntk.sim_patterns[j].weight = Inew;
                }
                else
                 Inew = ntk.sim_patterns[j].weight;
              }

              if( Inew >= Imax )
              {
                active_list.second[i] = active_list_tmp.second[i];
                active_list.first[i] = active_list_tmp.first[i];
                Iactive[i] = Inew;
                Imax = Inew;
              }
            }
          }

          active_list_tmp.first[i] = active_list.first[i];
          active_list_tmp.second[i] = active_list.second[i];
          
          ntk.sim_patterns[ntk.nodes_to_patterns[ ntk.get_node( (*active_list.second[i]).sig ) ]].flag = true;
        }
        if( ps.verbose )
        {
          for( auto x : ntk.sim_patterns )
          {       
            std::cout << x.sig << ":" << x.weight << " ";
          }
          std::cout << "\nact " << act << std::endl;
          std::cout << "mi(A;y) =mi([ ";
          for( auto x : active_list.second )
            std::cout << (*x).sig << " " ;
          std::cout << "])=" << Iactive[act-1] << std::endl;
        }
      }

      bool not_done( )
      {
        fill_active_list(1);
        double eps_nd = Iactive[0]/kitty::entropy( std::vector{&target} );
        if( ps.verbose )
          std::cout << "E " << eps_nd << std::endl;
        return ( eps_nd < ps.eps_th );
      }

      struct best_function_res
      {
        std::vector<signal<Ntk>> children = {};
        double mi;
        std::string tt;
        kitty::partial_truth_table pat;
        kitty::dynamic_truth_table dtt;
      };

      best_function_res best_function_2( size_t const& act )
      {
        std::vector< TT * > patterns_support = { active_list.first[ act  - 1 ], &ntk.sim_patterns[0].pat };
        std::vector< TT * > active_list_tt = active_list.first;
        create_candidates_result<TT> candidates;
        best_function_res res;
        res.children.push_back( (*active_list.second[ act - 1 ]).sig );
        res.children.push_back( (*active_list.second[ 0 ]).sig );

        double Imax = 0.0;
        bool none_considered = true;


        for( size_t i = 0 ; i < ntk.sim_patterns.size(); ++i )
        {
          if( active_list.first[ act - 1 ] == &ntk.sim_patterns[i].pat )
            continue;
          
          patterns_support[1] = &ntk.sim_patterns[i].pat;
          candidates = create_candidates_method( patterns_support, &target );
          for( size_t j = 0; j < candidates.tt_v.size(); ++j )
          {
            active_list_tt[ act - 1 ] = &candidates.pat_v[j];
            double Inew = kitty::mutual_information( active_list_tt, &target );

            if( ( Inew > Imax ) || ( none_considered && ( Inew >= Imax ) ) )
            {
              none_considered = false;
              Imax = Inew;
              res.mi = Inew;
              res.pat = candidates.pat_v[j];
              res.tt = candidates.tt_v[j];
              res.dtt = candidates.dtt_v[j];
              res.children[1] = ntk.sim_patterns[i].sig;
            }
          }
        }
        return res;
      }

      best_function_res best_function( size_t const& act, size_t const& sup )
      {
        assert( sup == 2 );
        best_function_res res;

        switch( sup )
        {
          case 2 :
            res = best_function_2( act );
            break;
          default:
            std::cerr << "[e] method for support of size " << sup << " is not implemented" << std::endl;
        }
        return res;
      }

      void add_node( best_function_res const& best_fn )
      {
        kitty::dynamic_truth_table tt(2u);
        kitty::create_from_binary_string( tt, best_fn.tt );
        auto fnew = ntk.create_node( best_fn.children, tt );
        if( ps.verbose )
        {
          std::cout <<" select: " << fnew << "= " << best_fn.children[1] << " " << best_fn.children[0] << " " << best_fn.tt << std::endl;  
          std::cout << "* "; kitty::print_binary(best_fn.pat);
          std::cout << std::endl;
        }

      }

      bool is_trivial( std::string tt )
      {
        if( tt == "0000" || tt == "1111" || tt == "1100" || tt == "0011" || tt == "1010" || tt == "0101" )
          return true;
        else
          return false;
      }

      bool improve_mi( size_t const& act, size_t const& sup )
      {
        fill_active_list( act );
        best_function_res best_fn = best_function( act, sup );
        double Iold = Iactive[ act - 1 ];
        double Inew = best_fn.mi;
        if( is_trivial( best_fn.tt ) || ( best_fn.dtt.num_bits() == 1 ) )
        {
          if( ps.verbose )
          {
            std::cout << "Fails to find f(" << (*active_list.second[act - 1]).sig << ",?) with mi([ ";
            for( size_t i = 0; i < active_list.second.size()-1; ++i )
              std::cout << (*active_list.second[i]).sig << " ";   
            std::cout << "f ]) > " << Iactive[act-1] << std::endl;
          }
          return false;
        }
        
        if( Inew > Iold )
        {
          add_node( best_fn );
          return true;
        }
        else
        {
          if( ps.verbose )
            std::cout << "Fails to find f(" << (*active_list.second[act - 1]).sig << ",?) with mi([f]) > " << Iactive[act-1] << std::endl;
          return false;
        }
      }

      signal<Ntk> run()
      {
        size_t act = 0;
        bool success{ false };

        fill_active_list( 1 );
        size_t sup = 2;
        while( not_done() && ( sup < ps.max_sup ) )
        {
          act = 0;
          do
          {
            act++;
            success = improve_mi( act, sup );
          } while ( success == false && ( act < ps.max_act ) );
          if( success == true )
          {
            sup = 2;
            while( success == true )
              success = improve_mi( act, sup );
          }
          else
            sup++;
        }
        
        fill_active_list( 1 );
        return ntk.sim_patterns[ ntk.nodes_to_patterns[ ntk.get_node( (*active_list.second[0]).sig ) ] ].sig;
      }
    double accuracy( kitty::partial_truth_table const& A, kitty::partial_truth_table const& B )
    {
      return 100*(double)kitty::count_ones( ~(A^B) )/A.num_bits();
    }
    
    private:
      simulation_view<Ntk>& ntk;
      sim_muesli_params ps;
      kitty::partial_truth_table target;
    public:
      std::pair<std::vector<kitty::partial_truth_table*>, std::vector<sim_pattern<Ntk> * >> active_list;
      std::vector<double> Iactive;
}; /* class sim_muesli_impl */
    
} /* namespace detail */



#pragma region sim_muesli
/*! \brief sim_muesli algorithm assembles a network bottom up.
 *
 * This method iteratively creates and adds new informative nodes. 
 * It takes an empty network, the simulation patterns of its input nodes and the target functon. 
 * Then, the rest of the network is created using the strategy described in the paper
 * "Learning complex boolean functions: Algorithms and applications." by Oliveira, A., & Sangiovanni-Vincentelli, A. (1993).
 * \param ntk the network wrapped in a simulation view.
 * \param examples the input patterns to be used for assembling the network
 * \param target target simulation pattern
 * \param ps parameters of the sim_muesli method
 *
 */
template<typename Ntk>
signal<Ntk> sim_muesli( simulation_view<Ntk>& ntk, std::vector<kitty::partial_truth_table> & examples, kitty::partial_truth_table const & target, sim_muesli_params& ps )
{
  if( ps.re_initialize )
    ntk.initialize_network( examples );

  if( ps.verbose )
  {
    std::cout << "  ";
    for( size_t i = 0; i < ntk.sim_patterns[0].pat.num_bits(); ++i )
      std::cout << "-";
    std::cout << std::endl;
    for( auto x : ntk.sim_patterns )
    {
      std::cout << x.sig << " "; kitty::print_binary( x.pat ); std::cout << std::endl;
    }
    std::cout << "  ";
    for( size_t i = 0; i < ntk.sim_patterns[0].pat.num_bits(); ++i )
      std::cout << "-";
    std::cout << std::endl;
    std::cout << "y "; kitty::print_binary( target ); std::cout << std::endl;
    std::cout << "  ";
    for( size_t i = 0; i < ntk.sim_patterns[0].pat.num_bits(); ++i )
      std::cout << "-";
    std::cout << std::endl;
  }

  detail::sim_muesli_impl impl( ntk, target, ps );
  auto osignal = impl.run();
  if(ps.verbose)
    std::cout << "num gates = " << ntk.num_gates() << " accuracy = " << impl.accuracy( target, ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat ) << "%" << std::endl;

  if( ps.try_accuracy_recovery && ( ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat != target ) )
  {

    sim_decomposition_fastS_params decps;
    decps.verbose = ps.verbose;
    osignal = sim_decomposition_fastS( ntk, examples, target, decps );
    if( ps.verbose )
    {
      std::cout << "num gates = " << ntk.num_gates() << " accuracy = " << impl.accuracy( target, ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat ) << "%" << std::endl;
    }
  }
  
  return osignal;
}
#pragma endregion sim_muesli
} // namespace mockturtle