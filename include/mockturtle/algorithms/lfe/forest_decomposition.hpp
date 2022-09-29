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

namespace mockturtle
{
/*! \brief Parameters for forest_decomposition algorithm */
struct forest_decomposition_params
{
  bool verbose{false};
  uint32_t max_sup{2};
  bool is_informed{true};
  bool is_size_aware{false};
  bool try_top_decomposition{true};
  bool try_bottom_decomposition{false};
  bool use_correlation{false};
  bool branch_on_all{true};
  bool try_xor{false};
  uint32_t num_trees{3};
};

namespace detail
{

  template<typename Ntk>
  class forest_decomposition_impl
  {
    public:
      using TT = kitty::partial_truth_table;

    public:
      forest_decomposition_impl( simulation_view<Ntk> & ntk, TT target, forest_decomposition_params & ps )
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
        chj_result chj_res;
        chj_result chj_res_new_node;
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
              chj_res = chatterjee_method( support_pat_pointers, &on_f );
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
              std::vector<TT*> Xptr = std::vector{&xn};
              double In = kitty::mutual_information( Xptr, &on_f );
              Xptr = std::vector{&xl};
              double Ii = kitty::mutual_information( Xptr, &on_f );
              Xptr = std::vector{&xl, &xr};
              double Iij = kitty::mutual_information( Xptr, &on_f );
              Xptr = std::vector{&xl, &xr, &xn};
              double Iijn = kitty::mutual_information( Xptr, &on_f );

              if( ( In == Iij ) && ( In = Iijn ) && ( In > Ii ) && ( Inew >= Imax ) && ( !ps.is_size_aware  || (fanin_size < min_fanin_size) ) )
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

      }

      void clear_fanin_size( signal<Ntk> & sig )
      {
        ntk.clear_network_fanin_size_from_node(ntk.get_node(sig));
        ntk.update_network_fanin_size();
      }
    
      signal<Ntk> idsd_step( std::vector<uint32_t> support, TT amask, TT xmask, bool branch_on_last = false )
      {
        uint32_t n_ones = kitty::count_ones(amask);

        if( n_ones == 0 )
          return ntk.get_constant( false );    
        if( support.size() == 0 )
          return ntk.get_constant( false );

        TT on_f(n_bits); 
        TT off_f(n_bits); 

        on_f = amask&( xmask^Y.pat );
        off_f = amask&(~(xmask^Y.pat ));

        if( kitty::count_ones( on_f ) == 0 ) // contraddiction
          return ntk.get_constant( false );
        else if( kitty::count_ones( on_f ) == n_ones ) // tautology
          return ntk.get_constant( true );

        uint32_t bidx = 0;
        signal<Ntk> bsig;
        double Inew = -std::numeric_limits<double>::max();
        double Imax = -std::numeric_limits<double>::max();
        uint32_t max_fanin_size = std::numeric_limits<uint32_t>::max();

        std::vector<uint32_t> to_be_deleted;
        TT on_x(n_bits);
        TT off_x(n_bits);
        TT on_xi(n_bits); 
        TT off_xi(n_bits); 
        std::vector<uint32_t> maxI_indeces;
        
        if( ps.is_informed )
        {
          if( branch_on_last )
          {
            bidx = support.size()-1;
            on_x = amask & X[support[bidx]].pat;
            off_x = amask & ~X[support[bidx]].pat;
            bsig = X[support[bidx]].sig;

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
                to_be_deleted.push_back(i);
              else
              { 
                Inew = information( on_xi, off_xi, on_f, off_f );
                signal<Ntk> sig = X[support[i]].sig;
                if(( Inew > Imax ) || ( (Inew == Imax) && ( !ps.is_size_aware || (ntk.nodes_to_size_fanin[ntk.get_node(sig)] < max_fanin_size )) ) )
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
              return X[support[i]].sig;
            if( on_xi == off_f )
            {
              signal<klut_network>fo = ntk.create_not( X[support[i]].sig );
              if( ps.verbose )
                std::cout << fo << "=" << X[support[i]].sig << "'" << std::endl;
              return fo;
            } 

            if( ( on_xi == amask ) || ( off_xi == amask) ) 
              to_be_deleted.push_back(i);
          }
          bidx = 0;
        }

          
        if( to_be_deleted.size() > 0 )
        {
          std::reverse(to_be_deleted.begin(), to_be_deleted.end());
          uint32_t count = 0;
          for( uint32_t i = 0; i < to_be_deleted.size() ; ++i )
          {
            support.erase( support.begin()+to_be_deleted[i] ) ;
            if( to_be_deleted[i] < bidx )
              count++;
          }
          bidx -= count;
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
          sim_top_decomposition_fast res = is_top_decomposable_fast( X, pis_support, on_f, amask1, amask0, ps.try_xor );
        
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

        if( ps.try_bottom_decomposition )
        {
          if( ps.is_informed )
          {
            if ( try_bottom_decomposition( support, amask, on_f, off_f, Imax ) )
              return idsd_step( support, amask, xmask, true );
            
          }
          else
          {
            std::cout << "don't care-based not yet implemented" << std::endl;
          }
        }
        if( ps.is_size_aware ) 
          clear_fanin_size( bsig );

        signal<klut_network> F0 = idsd_step( reduced_support, amask0, xmask0 );
        signal<klut_network> f0 = ntk.create_and( ntk.create_not( bsig ), F0 );

        signal<klut_network> F1 = idsd_step( reduced_support, amask1, xmask1 );
        signal<klut_network> f1 = ntk.create_and( bsig, F1 );

        signal<klut_network> Fnew = ntk.create_or( f1, f0 );

        if( ps.verbose )
          std::cout << Fnew << "= ite(" << bsig << "," << F1 << "," << F0 << ")" << std::endl;

        return Fnew;
      }

      signal<Ntk> run()
      {
        std::vector<uint32_t> support;

        for( uint32_t i = 0; i < X.size(); ++i )
          support.push_back( i );
    

        TT xmask( Y.pat.num_bits() );
        TT amask = ~xmask;

        signal<Ntk> fout;

        if( ps.num_trees == 3 )
        {
          uint32_t edge1 = X[0].pat.num_bits()/3;
          uint32_t edge2 = X[0].pat.num_bits()*2/3;
          
          TT amask1 = amask;
          TT amask2 = amask;
          TT amask3 = amask;
          for( uint32_t i = edge2; i < X[0].pat.num_bits(); ++i )
            kitty::clear_bit( amask1, i );
          for( uint32_t i = 0; i < edge1; ++i )
            kitty::clear_bit( amask2, i );
          for( uint32_t i = edge1; i < edge2; ++i )
            kitty::clear_bit( amask3, i );
          
          signal<Ntk> f1 = idsd_step( support, amask1, xmask );
          signal<Ntk> f2 = idsd_step( support, amask2, xmask );
          signal<Ntk> f3 = idsd_step( support, amask3, xmask );
          fout = ntk.create_maj( f1, f2, f3 );
        }
        else if( ps.num_trees == 5 )
        {
          uint32_t edge1 = X[0].pat.num_bits()/5;
          uint32_t edge2 = X[0].pat.num_bits()*2/5;
          uint32_t edge3 = X[0].pat.num_bits()*3/5;
          uint32_t edge4 = X[0].pat.num_bits()*4/5;
          
          TT amask1 = amask;
          TT amask2 = amask;
          TT amask3 = amask;
          TT amask4 = amask;
          TT amask5 = amask;

          for( uint32_t i = edge4; i < X[0].pat.num_bits(); ++i )
            kitty::clear_bit( amask1, i );

          for( uint32_t i = 0; i < edge1; ++i )
            kitty::clear_bit( amask2, i );

          for( uint32_t i = edge1; i < edge2 ; ++i )
            kitty::clear_bit( amask3, i );

          for( uint32_t i = edge2; i < edge3; ++i )
            kitty::clear_bit( amask4, i );

          for( uint32_t i = edge3; i < edge4 ; ++i )
            kitty::clear_bit( amask5, i );

          std::vector<signal<Ntk>> children;

          children.push_back( idsd_step( support, amask1, xmask ) );
          children.push_back( idsd_step( support, amask2, xmask ) );
          children.push_back( idsd_step( support, amask3, xmask ) );
          children.push_back( idsd_step( support, amask4, xmask ) );
          children.push_back( idsd_step( support, amask5, xmask ) );

          kitty::dynamic_truth_table maj5(5u);
          std::string stt = "11111110111010001110100010000000";
          kitty::create_from_binary_string( maj5, stt );
          fout = ntk.create_node( children, maj5 );
        }
        else
          std::cerr << "not implemented yet" << std::endl;

        return fout;

      }
    
    private:
      simulation_view<Ntk>& ntk;
      forest_decomposition_params & ps;
      kitty::partial_truth_table target;
      sim_pattern<Ntk> Y;
      uint32_t n_bits;
      bool branch_on_all;
    public:
      std::vector<double> Iactive;
      std::vector<sim_pattern<Ntk>> X;
}; /* class forest_decomposition_impl */
    
} /* namespace detail */




#pragma region forest_decomposition
/*! \brief forest_decomposition algorithm assembles a network bottom up.
 *
 * This method iteratively creates and adds new informative nodes. 
 * It takes an empty network, the simulation patterns of its input nodes and the target functon. 
 * Then, the rest of the network is created using the strategy described in the paper
 * "Learning complex boolean functions: Algorithms and applications." by Oliveira, A., & Sangiovanni-Vincentelli, A. (1993).
 * \param ntk the network wrapped in a simulation view.
 * \param examples the input patterns to be used for assembling the network
 * \param target target simulation pattern
 * \param ps parameters of the forest_decomposition method
 *
 */
template<typename Ntk>
signal<Ntk> forest_decomposition( simulation_view<Ntk>& ntk, std::vector<kitty::partial_truth_table> & examples, kitty::partial_truth_table const & target, forest_decomposition_params & ps, bool re_initialize = true )
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

  detail::forest_decomposition_impl impl( ntk, target, ps );
  auto osignal = impl.run();
  return osignal;
}

#pragma endregion forest_decomposition


} // namespace mockturtle