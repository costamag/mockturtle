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
  \file decomposition.hpp
  \brief decomposition algorithm

  \author Andrea Costamagna
*/

#pragma once

#include <cstdint>
#include <vector>

#include "../../traits.hpp"
#include "../../views/depth_view.hpp"
#include "../simulation.hpp"
#include "simulation_view.hpp"
#include "chatterjee_method.hpp"
#include "sim_decomposition_fast_checks.hpp"
#include "sim_operations.hpp"
#include "sim_utils.hpp"
#include <stdlib.h> 

#include <fmt/format.h>
#include <kitty/constructors.hpp>
#include <kitty/decomposition.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/statistics.hpp>
#include <kitty/operators.hpp>
#include <kitty/print.hpp>

#include <mockturtle/algorithms/lfe/create_candidates.hpp>

namespace mockturtle
{
/*! \brief Parameters for sim_decomposition_xor algorithm */
struct sim_decomposition_xor_params
{
  bool verbose{true};
  uint32_t max_sup{2};
  bool is_informed{true};
  bool is_size_aware{false};
  bool try_top_decomposition{true};
  bool try_bottom_decomposition{false};
  bool use_correlation{false};
  bool branch_on_all{true};
  bool try_xor{false};
  bool is_relaxed{false};
  bool is_dc{false};
  int nImpurity{0};
};

namespace detail
{

  template<typename Ntk>
  class sim_decomposition_xor_impl
  {
    public:
      using TT = kitty::partial_truth_table;

    public:
      sim_decomposition_xor_impl( simulation_view<Ntk> & ntk, TT target, sim_decomposition_xor_params & ps )
      : ntk(ntk),
        ps(ps),
        target(target),
        Y(target),
        n_bits(target.num_bits())
      {
        if( ps.branch_on_all )
        {
          X = ntk.sim_patterns;
          X.erase( X.begin() );
          X.erase( X.begin() );
        }
        else
        {
          X = ntk.input_patterns;
        }
        assert( X[0].pat.num_bits() == target.num_bits() );
        for( uint32_t i = 0; i<ps.max_sup; ++i )
        {
          kitty::dynamic_truth_table x(ps.max_sup);
          kitty::create_nth_var( x, i ); 
          InSims.push_back( x );
          kitty::print_binary( InSims[i] );
          std::cout << std::endl;
        }
      }

      std::vector<sim_pattern<Ntk>*> GenerateSupport()
      {
        std::vector<sim_pattern<Ntk>*> res;
        std::vector<TT*> Xptr;
        double Inew, Imax ;
        std::vector<uint32_t> vSelectedVars;

        uint32_t idx = rand()%X.size();
        Xptr.push_back( &X[ idx ].pat );
        res.push_back( &X[idx] );
        vSelectedVars.push_back(2*X.size());

        uint32_t iVarBest;
        for( uint32_t iVar=1; iVar < ps.max_sup; ++iVar )
        {
          Xptr.push_back( &X[0].pat );
          res.push_back( &X[0] );
          vSelectedVars.push_back(2*X.size());
          Imax = -std::numeric_limits<double>::max();
          for( uint32_t i = 0; i < X.size(); ++i )
          {
            if( std::find(vSelectedVars.begin(), vSelectedVars.end(), i) != vSelectedVars.end() )
              continue;
            
            Xptr[iVar] = &X[i].pat;
            Inew = kitty::mutual_information( Xptr, &target );
            if( Inew >= Imax )
            {
              Imax = Inew;
              vSelectedVars[iVar] = i;
              res[iVar] = &X[i];
            }
          }
        }

        for( uint32_t i = 1; i < res.size(); ++i )
        {
          sim_pattern<Ntk> * key = res[i];
          int j = i-1;
          bool UPD = true;
          while( j >= 0 && UPD )
          { 
            if( key->sig < res[j]->sig )
            {
              res[j+1] = res[j];
              j-=1;
            }
            else
              UPD = false;
          } 
          res[j+1]=key;
        }

        return res;
      }

      std::vector<sim_pattern<Ntk>*> RandomGenerateSupport()
      {
        std::vector<sim_pattern<Ntk>*> res;
        std::vector<TT*> Xptr;
        double Inew, Imax ;
        std::vector<uint32_t> vSelectedVars;


        uint32_t idx = rand()%X.size();
        Xptr.push_back( &X[ idx ].pat );
        res.push_back( &X[idx] );
        vSelectedVars.push_back(idx);
        uint32_t iVar = 0;

        while( iVar < ps.max_sup )
        {
          idx = rand()%X.size();
          if( std::find(vSelectedVars.begin(), vSelectedVars.end(), idx ) != vSelectedVars.end() )
              continue;
          Xptr.push_back( &X[ idx ].pat );
          res.push_back( &X[idx] );
          vSelectedVars.push_back(idx);
          
          iVar++;
        }

        for( uint32_t i = 1; i < res.size(); ++i )
        {
          sim_pattern<Ntk> * key = res[i];
          uint32_t j = i-1;
          bool UPD = true;
          while( j >= 0 || UPD )
          { 
            if( key->sig < res[j]->sig )
              res[j+1]=res[j];
            else
              UPD = false;
            j--;
          } 
          res[j+1]=key;
        }
          
        return res;
      }

      create_candidates_result<TT> GenerateApproximation( std::vector<sim_pattern<Ntk>*> Support )
      {
        sim_pattern<Ntk> res;

        
        
        std::vector<signal<Ntk>> children;
        for( auto s : Support )
          children.push_back( s->sig );
        
        std::vector<TT*> sim_pats_ptr;
        for( uint32_t i = 0; i < Support.size(); ++i )
          sim_pats_ptr.push_back( &(Support)[i]->pat);

        create_candidates_result<TT> cand_res = create_candidates_method( sim_pats_ptr, &target );
        return cand_res;
      }      

  template< typename TT >
  struct DTSyn_result
  {
    TT pat;
    kitty::dynamic_truth_table dtt;
  };


      double information( TT & on_xi, TT & off_xi, TT & on_f, TT & off_f )
      {
        double n0, n1, n00, n01, n10, n11;
        n0 = (double)kitty::count_ones( off_xi );
        n0 = ( n0 == 0 ? 0 : n0*log2(n0)); 
        n1 = (double)kitty::count_ones( on_xi );
        n1 = ( n1 == 0 ? 0 : n1*log2(n1)); 

        n00 = (double)kitty::count_ones( off_xi & off_f );
        n00 = ( n00 == 0 ? 0 : n00*log2(n00)); 
        n01 = (double)kitty::count_ones( off_xi & on_f );
        n01 = ( n01 == 0 ? 0 : n01*log2(n01)); 
        n10 = (double)kitty::count_ones( on_xi & off_f );
        n10 = ( n10 == 0 ? 0 : n10*log2(n10)); 
        n11 = (double)kitty::count_ones( on_xi & on_f );
        n11 = ( n11 == 0 ? 0 : n11*log2(n11)); 

        return ( n00 + n01 + n10 + n11 - n0 - n1 );
      }

DTSyn_result<TT> simple_decompose( TT func, TT mask, std::vector<uint32_t> iSupport, std::vector<sim_pattern<Ntk>*> Support )
{
    double Inew, Imax = -std::numeric_limits<double>::max();
    TT on_f, off_f, on_x, off_x;
    DTSyn_result<TT> FR;

    on_f  = func & mask ;
    off_f = ~func & mask ;

    if( kitty::count_ones(on_f) == 0 )
    {
      FR.pat = mask& ~mask;
      FR.dtt = InSims[0] & ~InSims[0];
      return FR;
    }
  
    if( kitty::count_ones(off_f) == 0 )
    {
      FR.pat = mask | ~mask;
      FR.dtt = InSims[0] | ~InSims[0];
      return FR;
    }

    if( iSupport.size() == 1 )
    {
      if( kitty::count_ones( ( func ^ Support[iSupport[0]]->pat ) & mask ) > kitty::count_ones( ( func ^ ~Support[iSupport[0]]->pat ) & mask ) )
      {
        FR.pat = ~Support[iSupport[0]]->pat;
        FR.dtt = ~InSims[iSupport[0]];
        return FR;
      }
      else
      {
        FR.pat = Support[iSupport[0]]->pat;
        FR.dtt = InSims[iSupport[0]];
        return FR;     
      }
    }
    
    uint32_t idx, Id;
    for( uint32_t i = 0; i<iSupport.size(); ++i )
    {
      on_x  =  Support[iSupport[i]]->pat & mask ;
      off_x = ~Support[iSupport[i]]->pat & mask ;

      Inew = information( on_x, off_x, on_f, off_f );
      if( Inew > Imax )
      {
        Imax = Inew;
        idx = iSupport[i];
        Id = i;
      }
    }

    std::vector<uint32_t> iRedSupport = iSupport;
    iRedSupport.erase( iRedSupport.begin()+Id );

    DTSyn_result<TT> F0 = simple_decompose( func, mask &  Support[idx]->pat, iRedSupport,Support );

    DTSyn_result<TT> F1 = simple_decompose( func, mask & ~Support[idx]->pat, iRedSupport, Support );
  
    FR.pat = ( Support[idx]->pat & F1.pat ) | ( ~Support[idx]->pat & F0.pat );
    FR.dtt = ( InSims[idx] & F1.dtt ) | ( ~InSims[idx] & F0.dtt );
    return FR;
}

      DTSyn_result<TT> GenerateApproximationDT( std::vector<sim_pattern<Ntk>*> Support )
      {
        sim_pattern<Ntk> res;
        
        std::vector<signal<Ntk>> children;
        for( auto s : Support )
          children.push_back( s->sig );
        
        std::vector<TT*> sim_pats_ptr;
        std::vector<uint32_t> iSupport;
        for( uint32_t i = 0; i < Support.size(); ++i )
        {
          iSupport.push_back(i);
          sim_pats_ptr.push_back( &(Support)[i]->pat);
        }

        TT mask = ~target.construct();
        DTSyn_result<TT> cand_res = simple_decompose( target, mask, iSupport, Support );

        return cand_res;
      } 

      signal<Ntk> xdec_step()
      {
        
        bool isUpd = false;
        uint32_t n_ones = kitty::count_ones( target );
        std::cout << "#ones: " << n_ones << std::endl;
        if( n_ones <= ps.nImpurity )
          return ntk.get_constant( false ); 
        if( n_ones >= target.num_bits()-ps.nImpurity )
          return ntk.get_constant( true );

        while( !isUpd )
        {
          std::vector<sim_pattern<Ntk>*> support = GenerateSupport();
          DTSyn_result<TT> P = GenerateApproximationDT( support );

          std::vector<signal<Ntk>> children;
          for( auto s : support )
            children.push_back( s->sig );

          uint32_t idx;
          uint32_t newError = 0;
          /*for( uint32_t i = 0; i < P.sC_v.size(); ++i )
          {
            newError = std::min( kitty::count_ones( P.pat_v[i] ^ target ), kitty::count_ones( ~ P.pat_v[i] ^ target ) );
            if( newError < error-5 )
            {
              error = newError;
              idx = i;
              isUpd = true;
            }
          }   */  
          newError = std::min( kitty::count_ones( P.pat ^ target ), kitty::count_ones( ~ P.pat ^ target ) );
          std::cout << newError << std::endl;
          if( newError < error-5 )
          {
            error = newError;
            isUpd = true;
          }   
          if( isUpd )
          {
            signal<Ntk> fc = ntk.create_node( children, P.dtt );
            //if( ps.verbose )
            {
              std::cout << fc << " = ";
              for( uint32_t i = 0; i < children.size(); ++i )
                std::cout << children[i] << " ";
              kitty::print_binary( P.dtt );
              std::cout << std::endl;
            }
            
            if( kitty::count_ones(P.pat ^ target) > target.num_bits()/2 )
            {  
              target = ~P.pat ^ target;
              return ntk.create_xor( !fc, xdec_step() );
            }
            else
            {  
              target = P.pat ^ target;
              return ntk.create_xor( fc, xdec_step() );
            }
          }
          if( cnter++ > 1000 ) return ntk.get_constant( false ); 
        }
      }

      signal<Ntk> run()
      {
        std::cout << "run" << std::endl;
        cnter=0;
        return xdec_step();
      }
    
    private:
      simulation_view<Ntk>& ntk;
      sim_decomposition_xor_params & ps;
      kitty::partial_truth_table target;
      sim_pattern<Ntk> Y;
      uint32_t n_bits;
      bool branch_on_all;
    public:
      std::vector<double> Iactive;
      std::vector<kitty::dynamic_truth_table> InSims;
      std::vector<sim_pattern<Ntk>> X;
      uint32_t cStag=0;
      std::vector<sim_pattern<Ntk>*> old_support;
      uint32_t error;
      uint32_t cnter;

}; /* class sim_decomposition_xor_impl */
    
} /* namespace detail */




#pragma region sim_decomposition_xor
/*! \brief sim_decomposition_xor algorithm assembles a network bottom up.
 *
 * This method iteratively creates and adds new informative nodes. 
 * It takes an empty network, the simulation patterns of its input nodes and the target functon. 
 * Then, the rest of the network is created using the strategy described in the paper
 * "Learning complex boolean functions: Algorithms and applications." by Oliveira, A., & Sangiovanni-Vincentelli, A. (1993).
 * \param ntk the network wrapped in a simulation view.
 * \param examples the input patterns to be used for assembling the network
 * \param target target simulation pattern
 * \param ps parameters of the sim_decomposition_xor method
 *
 */
template<typename Ntk>
signal<Ntk> sim_decomposition_xor( simulation_view<Ntk>& ntk, std::vector<kitty::partial_truth_table> & examples, kitty::partial_truth_table const & target, sim_decomposition_xor_params & ps, bool re_initialize = true )
{
  if( re_initialize )
    ntk.initialize_network( examples );

  if( ps.verbose )
  {
    std::cout << "  ";
    for( uint32_t i = 0; i < ntk.sim_patterns[0].pat.num_bits(); ++i )
      std::cout << "-";
    std::cout << std::endl;
    for( auto x : ntk.sim_patterns )
    {
      std::cout << x.sig << " "; kitty::print_binary( x.pat ); std::cout << std::endl;
    }
    std::cout << "  ";
    for( uint32_t i = 0; i < ntk.sim_patterns[0].pat.num_bits(); ++i )
      std::cout << "-";
    std::cout << std::endl;
    std::cout << "y "; kitty::print_binary( target ); std::cout << std::endl;
    std::cout << "  ";
    for( uint32_t i = 0; i < ntk.sim_patterns[0].pat.num_bits(); ++i )
      std::cout << "-";
    std::cout << std::endl;
  }
  detail::sim_decomposition_xor_impl impl( ntk, target, ps );
  auto osignal = impl.run();
  return osignal;
}

#pragma endregion sim_decomposition_xor


} // namespace mockturtle