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
/*! \brief Parameters for sim_decomposition_fastS algorithm */
struct sim_decomposition_fastS_params
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
  class sim_decomposition_fastS_impl
  {
    public:
      using TT = kitty::partial_truth_table;

    public:
      sim_decomposition_fastS_impl( simulation_view<Ntk> & ntk, TT target, sim_decomposition_fastS_params & ps )
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
      }


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

      signal<Ntk> synthesize_leaf( std::vector<uint32_t> & support, std::vector<sim_pattern<Ntk>> & X, TT & amask, TT & on_f )
      {
        std::vector<TT> sim_pats;
        std::vector<TT*> sim_pats_ptr;

        for( uint32_t i = 0; i < support.size(); ++i )
        {
          TT sim_pat ;
          for( uint32_t j{0}; j < amask.num_bits(); ++j )
          {
            if( kitty::get_bit( amask, j ) == 1 )
              sim_pat.add_bit( kitty::get_bit( X[support[i]].pat, j ) );
          }
          sim_pats.push_back(sim_pat);
        
        }
        
        for( uint32_t i = 0; i < sim_pats.size(); ++i )
          sim_pats_ptr.push_back( &sim_pats[i] );

        TT target;
        for( uint32_t j{0}; j < amask.num_bits(); ++j )
        {
          if( kitty::get_bit( amask, j ) == 1 )
            target.add_bit( kitty::get_bit( on_f, j ) );
        }

        TT * target_ptr = &target;

        std::vector<signal<Ntk>> children;
        for( auto s : support )
          children.push_back( X[s].sig );

        chj_result chj_res = chatterjee_method( sim_pats_ptr, target_ptr, 123 );

        signal<Ntk> fc = ntk.create_node( children, chj_res.dtt );
        if( ps.verbose )
        {
          std::cout << fc << " = ";
          for( uint32_t i = 0; i < children.size(); ++i )
            std::cout << children[i] << " ";
          kitty::print_binary(chj_res.dtt);
          std::cout << std::endl;
        }

        return fc;
      }

      bool try_bottom_decomposition( std::vector<uint32_t>& support, TT & amask, TT & on_f, TT & off_f, double Imax )
      {
        std::vector<double> vect_I;
        std::vector<uint32_t> sorted_indeces;
        double Inew;
        bool is_success = false;
        TT on_xi;
        TT off_xi;
        double Llim = ps.is_relaxed ? (support.size() < 256 ? 0.95 : 1.00) : 1.00;//;
        double Rlim = ps.is_relaxed ? (support.size() < 256 ? 1.05 : 1.00) : 1.00;//1.00;//support.size() < 256 ? 1.05 : 1.00;
        double RTIO = ps.is_relaxed ? (support.size() < 256 ? 0.01 : 0.00) : 0.00;//0.00;//support.size() < 256 ? 0.01 : 0.0;

        chj_result chj_new_node;
        /*{
          std::string tt;
          kitty::partial_truth_table pat;
          bool both_not0_and_eq {false};
          bool both_not0 {false};
          kitty::dynamic_truth_table dtt;
        };*/

        for( uint32_t i = 0; i < support.size(); ++i )
        {
          on_xi = amask & X[support[i]].pat;
          off_xi = amask & ~X[support[i]].pat;
          Inew = information( on_xi, off_xi, on_f, off_f );
          if( (vect_I.size() == 0 ) || ( Inew < vect_I[vect_I.size()-1] ) )
          {
            vect_I.push_back( Inew );
            sorted_indeces.push_back( i );
          }
          else
          {
            for( uint32_t j = 0; j < vect_I.size(); ++j )
            {
              if( Inew >= vect_I[j] )
              {
                sorted_indeces.insert( sorted_indeces.begin()+j, i );
                vect_I.insert( vect_I.begin()+j, Inew );
                break;
              }
            }
          }
        }

        //uint32_t fanin_size;
        //uint32_t min_fanin_size = std::numeric_limits<uint32_t>::max();
        std::vector<TT*> support_pat_pointers;
        support_pat_pointers.push_back( &X[0].pat );
        support_pat_pointers.push_back( &X[0].pat );
        std::vector<uint32_t> new_node_support_indeces = {0, 0};

        for( uint32_t i = 0; i < sorted_indeces.size()-1; ++i )
        {
          support_pat_pointers[0] = &X[support[sorted_indeces[i]]].pat;
          for( uint32_t j = i+1; j < sorted_indeces.size(); ++j )
          {
            support_pat_pointers[1] = &X[support[sorted_indeces[j]]].pat;
            //if( vect_I[j] < vect_I[i] )

            double rtio = vect_I[j] == 0 ? 0 : abs( (vect_I[i] - vect_I[j] )/vect_I[i]);
            //if( vect_I[j] >= vect_I[i] )
            //  std::cout << rtio << std::endl;
            if( ( (vect_I[i] > vect_I[j]) && ( rtio > RTIO ) )  )
              break;
            else
            {
              //\std::cout << i << " " << vect_I[i] << " " << j << " " << vect_I[j] << std::endl; 

              create_candidates_result<TT> candidates = create_candidates_method( support_pat_pointers, &on_f );
              for( uint32_t k = 0; k < candidates.dtt_v.size(); ++k )
              {
                on_xi = amask & candidates.pat_v[k];
                off_xi = amask & ~candidates.pat_v[k];
                Inew = information( on_xi, off_xi, on_f, off_f );
////////////////////////////////////////////////////////////////


                if( Inew > Imax )
                {
                  TT xl = amask & X[support[sorted_indeces[i]]].pat;
                  TT xr = amask & X[support[sorted_indeces[j]]].pat;
                  TT ym = amask & on_f;
                  TT xn = on_xi;

                  std::vector<TT*> Xptr;
                  //A Xptr = std::vector{&xl, &xn};
                  //A double Iin = kitty::mutual_information( Xptr, &ym );
                  //A Xptr = std::vector{&xr, &xn};
                  //A double Ijn = kitty::mutual_information( Xptr, &ym );
                  //A if( Iin == Ijn )
                  {
                  Xptr = std::vector{&xl, &xr};
                  double Iij = kitty::mutual_information( Xptr, &ym );
                  //A if( Iij == Iin )
                  {
                  Xptr = std::vector{&xn};
                  double In = kitty::mutual_information( Xptr, &ym );
                  Xptr = std::vector{&xl,&xn};
                  if( ( In >= Llim*Iij && In <= Rlim*Iij  ) )
                  {
                  //  Xptr = std::vector{&xl, &xr, &xn};
                  //  double Iijn = kitty::mutual_information( Xptr, &ym );
                 //if( Iijn == In && Iin == In && Ijn == In && Iij == In )
                 //if(  Iijn == Iin )//&& Iin == In && Ijn == In && Iij == In )
///////////////////////////////////////////////////////
                {
                  std::vector<signal<Ntk>> children;
                  children.push_back( X[support[sorted_indeces[i]]].sig );
                  children.push_back( X[support[sorted_indeces[j]]].sig );
                  auto cand = std::make_pair( children, candidates.tt_v[k] );
                  if( ntk.available_nodes.find(cand) == ntk.available_nodes.end() )
                  {
                    Imax = Inew;
                    new_node_support_indeces[0] = sorted_indeces[i];
                    new_node_support_indeces[1] = sorted_indeces[j];
                    chj_new_node.tt = candidates.tt_v[k];
                    chj_new_node.dtt = candidates.dtt_v[k];
                    chj_new_node.pat = candidates.pat_v[k];
                    //min_fanin_size = fanin_size;
                    ntk.available_nodes.insert(cand);
                    is_success = true;
                  }
                }
                }
                  }
                  }
                }
///////////////////////////////////////////////////////
              }
            }
          }
        }
        if( is_success )
        {
          std::vector<signal<Ntk>> children;
          children.push_back( X[support[new_node_support_indeces[0]]].sig );
          children.push_back( X[support[new_node_support_indeces[1]]].sig );
          signal<Ntk> fc = ntk.create_node( children, chj_new_node.dtt );
          support.push_back( X.size() );
          X.push_back( ntk.sim_patterns[ ntk.get_node_pattern( fc ) ] );
        }
        return is_success;
      }


      /*bool try_bottom_decomposition( std::vector<uint32_t>& support, TT & amask, TT & on_f, TT & off_f, double Imax )
      {
        std::vector<double> vect_I;
        std::vector<uint32_t> sorted_indeces;
        double Inew;
        bool is_success = false;
        TT on_xi;
        TT off_xi;
        for( uint32_t i = 0; i < support.size(); ++i )
        {
          on_xi = amask & X[support[i]].pat;
          off_xi = amask & ~X[support[i]].pat;
          Inew = information( on_xi, off_xi, on_f, off_f );
          if( (vect_I.size() == 0 ) || ( Inew < vect_I[vect_I.size()-1] ) )
          {
            vect_I.push_back( Inew );
            sorted_indeces.push_back( i );
          }
          else
          {
            for( uint32_t j = 0; j < vect_I.size(); ++j )
            {
              if( Inew >= vect_I[j] )
              {
                sorted_indeces.insert( sorted_indeces.begin()+j, i );
                vect_I.insert( vect_I.begin()+j, Inew );
                break;
              }
            }
          }
        }

        uint32_t fanin_size;
        uint32_t min_fanin_size = std::numeric_limits<uint32_t>::max();
        std::vector<TT*> support_pat_pointers;
        support_pat_pointers.push_back( &X[0].pat );
        support_pat_pointers.push_back( &X[0].pat );

        std::vector<uint32_t> new_node_support_indeces = {0, 0};

        for( uint32_t i = 0; i < sorted_indeces.size()-1; ++i )
        {
          support_pat_pointers[0] = &X[support[sorted_indeces[i]]].pat;
          for( uint32_t j = i+1; j < sorted_indeces.size(); ++j )
          {
            support_pat_pointers[1] = &X[support[sorted_indeces[j]]].pat;
            if( vect_I[i] != vect_I[j] )
              break;
            else
            {
              //chj_res = chatterjee_method( support_pat_pointers, &on_f );
              create_candidates_result<TT> candidates = create_candidates_method( support_pat_pointers, &on_f );
              for( auto chj_res : candidates )
              {
                on_xi = amask & chj_res.pat;
                off_xi = amask & chj_res.pat;

                Inew = information( on_xi, off_xi, on_f, off_f );
                if( ps.is_size_aware )
                {
                  fanin_size = ntk.nodes_to_size_fanin[ntk.get_node(X[support[new_node_support_indeces[0]]].sig)]+
                                 ntk.nodes_to_size_fanin[ntk.get_node(X[support[new_node_support_indeces[1]]].sig)]+1;
                }

                TT xl = amask & X[support[sorted_indeces[i]]].pat;
                TT xr = amask & X[support[sorted_indeces[j]]].pat;
                TT xn = on_xi;
                std::vector<TT*> Xptr;
                Xptr = std::vector{&xn};
                double In = kitty::mutual_information( Xptr, &on_f, &amask );
                Xptr = std::vector{&xl,&xn};
                double Iin = kitty::mutual_information( Xptr, &on_f, &amask );
                Xptr = std::vector{&xr, &xn};
                double Ijn = kitty::mutual_information( Xptr, &on_f, &amask );
                Xptr = std::vector{&xl, &xr, &xn};
                double Iijn = kitty::mutual_information( Xptr, &on_f, &amask );

                if( ( In == Iijn ) && ( In > Ii ) && ( Inew >= Imax ) && ( !ps.is_size_aware  || (fanin_size < min_fanin_size) ) )
                {
                  std::vector<signal<Ntk>> children;
                  children.push_back( X[support[sorted_indeces[i]]].sig );
                  children.push_back( X[support[sorted_indeces[j]]].sig );
                  auto cand = std::make_pair( children, chj_res.tt );
                  if( ntk.available_nodes.find(cand) == ntk.available_nodes.end() )
                  {
                    Imax = Inew;
                    new_node_support_indeces[0] = sorted_indeces[i];
                    new_node_support_indeces[1] = sorted_indeces[j];
                    chj_res_new_node = chj_res;
                    min_fanin_size = fanin_size;
                    ntk.available_nodes.insert(cand);
                    is_success = true;
                  }
                }
              }
            }
          }
        }

        if( is_success )
        {

          std::vector<signal<Ntk>> children;
          children.push_back( X[support[new_node_support_indeces[0]]].sig );
          children.push_back( X[support[new_node_support_indeces[1]]].sig );
          signal<Ntk> fc = ntk.create_node( children, chj_res_new_node.dtt );

          support.push_back( X.size() );
          X.push_back( ntk.sim_patterns[ ntk.get_node_pattern( fc ) ] );
          //support.erase( support.begin() + std::max<uint32_t>( new_node_support_indeces[0],new_node_support_indeces[1] ) );
          //support.erase( support.begin() + std::min<uint32_t>( new_node_support_indeces[0],new_node_support_indeces[1] ) );

        }
        return is_success;

      }*/

      void clear_fanin_size( signal<Ntk> & sig )
      {
        ntk.clear_network_fanin_size_from_node(ntk.get_node(sig));
        ntk.update_network_fanin_size();
      }
    
      signal<Ntk> idsd_step( std::vector<uint32_t> support, TT amask, TT xmask, bool branch_on_last = false )
      {

        uint32_t n_ones = kitty::count_ones(amask);

        if( n_ones == 0 )
        {
          return ntk.get_constant( false ); 
        }   
        if( support.size() == 0 )
        {
          return ntk.get_constant( false );
        }
        TT on_f(n_bits); 
        TT off_f(n_bits); 

        on_f = amask&( xmask^Y.pat );
        off_f = amask&(~(xmask^Y.pat ));

        if( kitty::count_ones( on_f ) <= ps.nImpurity ) // contraddiction
        {
          return ntk.get_constant( false );
        }
        else if( kitty::count_ones( on_f ) >= n_ones-ps.nImpurity ) // tautology
        {
          return ntk.get_constant( true );
        }

        uint32_t bidx = 0;
        signal<Ntk> bsig;
        double Inew = -std::numeric_limits<double>::max();
        double Imax = -std::numeric_limits<double>::max();
        uint32_t max_fanin_size = std::numeric_limits<uint32_t>::max();
        

        std::vector<uint32_t> to_be_deleted;
        std::vector<uint32_t> to_be_deleted_idx;
        std::vector<uint32_t> new_support;
        TT on_x(n_bits);
        TT off_x(n_bits);
        TT on_xi(n_bits); 
        TT off_xi(n_bits); 
        std::vector<uint32_t> maxI_indeces;
        
        if( ps.is_informed )
        {
          if( branch_on_last )
          {
            bidx = uint32_t(support.size()-1);
            on_x = amask & X[support[bidx]].pat;
            off_x = amask & ~X[support[bidx]].pat;
            bsig = X[support[bidx]].sig;

            Imax = information( on_x, off_x, on_f, off_f );

            if( on_x == on_f )
              return X[support[bidx]].sig;
            if( on_x == off_f )
            {
              signal<klut_network>fo = ntk.create_not( X[support[bidx]].sig );
              if( ps.verbose )
              {
                std::cout << fo << "=" << X[support[bidx]].sig  << "'" << std::endl;
              }
              return fo;
            } 
          }
          else
          {
            for( uint32_t i = 0; i < support.size(); ++i )
            {
              on_xi = amask & X[support[i]].pat;
              off_xi = amask & ~X[support[i]].pat;

              if( on_xi == on_f )
                return X[support[i]].sig;
              if( on_xi == off_f )
              {
                signal<klut_network>fo = ntk.create_not( X[support[i]].sig );
                if( ps.verbose )
                  std::cout << fo << "=" << X[support[i]].sig << "'" << std::endl;
                return fo;
              } 

              if( ( on_xi == amask ) || ( off_xi == amask) ) 
              {
                to_be_deleted_idx.push_back( i );
                to_be_deleted.push_back( support[i] );
              }
              else
              { 
                Inew = information( on_xi, off_xi, on_f, off_f );
                signal<Ntk> sig = X[support[i]].sig;
                if(( Inew > Imax ) || ( (Inew == Imax) && ( !ps.is_size_aware || (ntk.nodes_to_size_fanin[ntk.get_node(sig)] <= max_fanin_size )) ) )
                {
                  Imax = Inew;
                  bidx = i;
                  bsig = sig;
                  max_fanin_size = ntk.nodes_to_size_fanin[ntk.get_node(sig)];
                }
              }
            }
          }
        }
        else
        {
          for( uint32_t i = 0; i < support.size(); ++i )
          {
            on_xi = amask & X[support[i]].pat;
            off_xi = amask & ~X[support[i]].pat;

            if( on_xi == on_f )
            {
              return X[support[i]].sig;
            }
            if( on_xi == off_f )
            {
              signal<klut_network>fo = ntk.create_not( X[support[i]].sig );
              if( ps.verbose )
                std::cout << fo << "=" << X[support[i]].sig << "'" << std::endl;
              return fo;
            } 

            if( ( on_xi == amask ) || ( off_xi == amask) ) 
            {
              to_be_deleted_idx.push_back( i );
              to_be_deleted.push_back(support[i]);
            }
          }
          bidx = 0;
        }


        if( to_be_deleted.size() > 0 )
        {
          for( auto x : support )
            if( std::find(to_be_deleted.begin(), to_be_deleted.end(), x ) == to_be_deleted.end() )
              new_support.push_back(x);
          //std::cout << "H" << BFcnt++ << std::endl;
          
          /*std::reverse(to_be_deleted.begin(), to_be_deleted.end());*/
          uint32_t count = 0;
          for( uint32_t i = 0; i < to_be_deleted.size() ; ++i )
          {
            //assert( support.size() > to_be_deleted[i] );
            //support.erase( support.begin()+to_be_deleted[i] ) ;
            if( to_be_deleted_idx[i] < bidx )
              count++;
          }
          bidx = uint32_t(bidx-count);
          //bidx -= count;

          //std::cout << "bidx " << bidx << std::endl;
          support=new_support;
        }

        if( support.size() == 0 )
          return ntk.get_constant( false );

        if( support.size() <= ps.max_sup )
          return synthesize_leaf( support, X, amask, on_f );

        on_x = amask & X[support[bidx]].pat;
        off_x = amask & ~X[support[bidx]].pat;
        bsig = X[support[bidx]].sig;

        TT amask1(n_bits);
        TT amask0(n_bits);

        amask1 = on_x;
        amask0 = off_x;

        TT xmask1(n_bits);
        TT xmask0(n_bits);

        xmask1 = on_x & xmask;
        xmask0 = off_x & xmask;
        
        std::vector<uint32_t> reduced_support = support;

        reduced_support.erase( reduced_support.begin() + bidx );
        //std::cout << "A4" << std::endl;

        std::vector<uint32_t> pis_support;
        if( ps.try_xor )
        {
          for( uint32_t k = 0; k < reduced_support.size(); ++k )
          {
            if( ntk.is_pi( ntk.get_node(X[reduced_support[k]].sig) ) )
              pis_support.push_back( reduced_support[k] );
          }
        }
        else
          pis_support = reduced_support;
        
        if( ps.is_informed && ps.try_top_decomposition )
        {
          sim_top_decomposition_fast res = is_top_decomposable_fast( X, pis_support, on_f, amask1, amask0, ps.try_xor, ps.is_dc );
        
          if ( res != sim_top_decomposition_fast::none )
          {
            switch ( res )
            {
              default:
                assert( false );
              case sim_top_decomposition_fast::and_:
              {
                signal<klut_network> F1 = idsd_step( reduced_support, amask1, xmask1 );
                signal<klut_network> Fnew = ntk.create_and( bsig, F1 );

                if( ps.is_size_aware ) clear_fanin_size( Fnew );
                if( ps.verbose ) std::cout << Fnew << "=" << bsig << " AND " << F1 << std::endl;
                return Fnew;
              }
              case sim_top_decomposition_fast::or_:
              {
                signal<klut_network> F0 = idsd_step( reduced_support, amask0, xmask0 );
                signal<klut_network> Fnew = ntk.create_or( bsig, F0 );
                if( ps.is_size_aware ) clear_fanin_size( Fnew );
                if( ps.verbose ) std::cout << Fnew << "=" << bsig << " OR " << F0 << std::endl;
                return Fnew;
              }
              case sim_top_decomposition_fast::lt_:
              {
                signal<klut_network> F0 = idsd_step( reduced_support, amask0, xmask0 );
                signal<klut_network> Fnew = ntk.create_lt( bsig, F0 );
                if( ps.is_size_aware ) clear_fanin_size( Fnew );
                if( ps.verbose ) std::cout << Fnew << "=" << bsig << "' AND " << F0 << std::endl;
                return Fnew;

              }
              case sim_top_decomposition_fast::le_:
              {  
                signal<klut_network> F1 = idsd_step( reduced_support, amask1, xmask1 );
                signal<klut_network> Fnew = ntk.create_le( bsig, F1 );
                if( ps.is_size_aware ) clear_fanin_size( Fnew );
                if( ps.verbose ) std::cout << Fnew << "=" << bsig << "' OR " << F1 << std::endl;
                return Fnew;
              }
              case sim_top_decomposition_fast::xor_:
              {
                xmask = xmask ^ on_x;
                signal<klut_network> Fxor = idsd_step( reduced_support, amask, xmask );
                signal<klut_network> Fnew = ntk.create_xor( bsig, Fxor );
                if( ps.is_size_aware ) clear_fanin_size( Fnew );
                if( ps.verbose ) std::cout << Fnew << "=" << bsig << " XOR " << Fxor << std::endl;
                return Fnew;
              }
            }
          }
        }

        if( !branch_on_last && ps.try_bottom_decomposition ) //support.size() < 256 && 
        {
          if( ps.is_informed )
          {
            if ( try_bottom_decomposition( support, amask, on_f, off_f, Imax ) )
            {
              //num_nobranch += 1;
              return idsd_step( support, amask, xmask, true );
            }
          }
          else
          { 
            std::cout << "don't care-based not yet implemented" << std::endl;
          }
        }
        if( ps.is_size_aware ) 
        {
          clear_fanin_size( bsig );
        }
      //std::cout << "A10-1" << std::endl;
        //std::cout << kitty::count_ones( amask0 ) << std::endl;
  //      std::cout << kitty::count_ones( ~amask0 ) << std::endl;
        //std::cout << reduced_support.size() << std::endl;
        //for( auto x : reduced_support )
        //  std::cout << x << " ";
        //std::cout << std::endl;
        
        klut_network::signal F0 = idsd_step( reduced_support, amask0, xmask0 );//reduced_support, amask0, xmask0 );
        //std::cout << "A10+1" << std::endl;

        klut_network::signal f0 = ntk.create_and( ntk.create_not( bsig ), F0 );
    //    std::cout << "A11" << std::endl;

        signal<klut_network> F1 = idsd_step( reduced_support, amask1, xmask1 );//( reduced_support, amask1, xmask1 );
    //    std::cout << "A12" << std::endl;

        signal<klut_network> f1 = ntk.create_and( bsig, F1 );
    //    std::cout << "A13" << std::endl;


        signal<klut_network> Fnew = ntk.create_or( f1, f0 );

        if( ps.verbose )
          std::cout << Fnew << "= ite(" << bsig << "," << F1 << "," << F0 << ")" << std::endl;

      //  std::cout << "b0" << std::endl;


        return Fnew;
      }

      signal<Ntk> run()
      {
        std::vector<uint32_t> support;

        for( uint32_t i = 0; i < X.size(); ++i )
          support.push_back( i );
        
        if( X.size() > 256 )
          size_thresh = uint32_t( 0.5*X.size() );
      

        TT xmask( Y.pat.num_bits() );
        TT amask = ~xmask;
    //    std::cout << "a1" << std::endl;
        return idsd_step( support, amask, xmask );
    //    std::cout << "a0" << std::endl;

      }
    
    private:
      simulation_view<Ntk>& ntk;
      sim_decomposition_fastS_params & ps;
      kitty::partial_truth_table target;
      sim_pattern<Ntk> Y;
      uint32_t n_bits;
      uint32_t size_thresh;
      bool branch_on_all;
    public:
      std::vector<double> Iactive;
      std::vector<sim_pattern<Ntk>> X;
}; /* class sim_decomposition_fastS_impl */
    
} /* namespace detail */




#pragma region sim_decomposition_fastS
/*! \brief sim_decomposition_fastS algorithm assembles a network bottom up.
 *
 * This method iteratively creates and adds new informative nodes. 
 * It takes an empty network, the simulation patterns of its input nodes and the target functon. 
 * Then, the rest of the network is created using the strategy described in the paper
 * "Learning complex boolean functions: Algorithms and applications." by Oliveira, A., & Sangiovanni-Vincentelli, A. (1993).
 * \param ntk the network wrapped in a simulation view.
 * \param examples the input patterns to be used for assembling the network
 * \param target target simulation pattern
 * \param ps parameters of the sim_decomposition_fastS method
 *
 */
template<typename Ntk>
signal<Ntk> sim_decomposition_fastS( simulation_view<Ntk>& ntk, std::vector<kitty::partial_truth_table> & examples, kitty::partial_truth_table const & target, sim_decomposition_fastS_params & ps, bool re_initialize = true )
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
  detail::sim_decomposition_fastS_impl impl( ntk, target, ps );
  auto osignal = impl.run();
  return osignal;
}

#pragma endregion sim_decomposition_fastS

#pragma region forest_fastS
/*! \brief sim_decomposition_fastS algorithm assembles a network bottom up.
 *
 * This method iteratively creates and adds new informative nodes. 
 * It takes an empty network, the simulation patterns of its input nodes and the target functon. 
 * Then, the rest of the network is created using the strategy described in the paper
 * "Learning complex boolean functions: Algorithms and applications." by Oliveira, A., & Sangiovanni-Vincentelli, A. (1993).
 * \param ntk the network wrapped in a simulation view.
 * \param examples the input patterns to be used for assembling the network
 * \param target target simulation pattern
 * \param ps parameters of the sim_decomposition_fastS method
 *
 */
/*
template<typename Ntk>
signal<Ntk> forest_fastS( simulation_view<Ntk>& ntk, std::vector<kitty::partial_truth_table> & examples, kitty::partial_truth_table const & target, sim_decomposition_fastS_params & ps, uint32_t ndecs = 3 )
{
  std::vector<signal<Ntk>> signal_voters;
  uint32_t num_elements = targets.num_bits()/ndecs;
  for( uint32_t i = 0; i < ndecs; ++i )
  {
    std::vector<kitty::partial_truth_table> examples_p(examples.size());
    kitty::partial_truth_table target_p;
    for( uint32_t j = i*num_elements ; j < std::min( (i+1)*num_elements, targets.num_bits() ); ++j  )
    {
      for( uint32_t c = 0; c < examples.size(); ++c )
        examples_p[c].add_bit( kitty::get_bit( examples[c], j ) );
      target_p.add_bit( kitty::get_bit( target, j ) );
    }
    ntk.initialize_network( examples_p );
    detail::sim_decomposition_fastS_impl impl( ntk, target_p, ps );

    signal<Ntk> osignal = impl.run();
    signal_voters.push_back( osignal ); 
  }

  
  return ntk.create_maj( osignal[0], osignal[1], osignal[2] );
}
*/
#pragma endregion forest_fastS

} // namespace mockturtle