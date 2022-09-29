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
/*! \brief Parameters for sim_decomposition_fast algorithm */
struct sim_decomposition_fast_params
{
  bool verbose{false};
  uint32_t max_sup{2};
  bool is_informed{true};
  bool try_top_decomposition{true};
  bool try_bottom_decomposition{false};
  bool use_correlation{false};
  bool try_xor{false};
};

namespace detail
{

  template<typename Ntk>
  class sim_decomposition_fast_impl
  {
    public:
      using TT = kitty::partial_truth_table;

    public:
      sim_decomposition_fast_impl( simulation_view<Ntk> & ntk, TT target, sim_decomposition_fast_params & ps )
      : ntk(ntk),
        ps(ps),
        target(target),
        Y(target),
        n_bits(target.num_bits())
      {
        X = ntk.sim_patterns;
        X.erase( X.begin() );
        X.erase( X.begin() );
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
          //uint32_t k = 0;
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
    
      signal<Ntk> idsd_step( std::vector<uint32_t> support, TT amask, TT xmask )
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
        std::vector<uint32_t> to_be_deleted;
        TT on_x(n_bits);
        TT off_x(n_bits);
        TT on_xi(n_bits); 
        TT off_xi(n_bits); 

          for( uint32_t i = 0; i < support.size(); ++i )
          {
            on_xi = amask & X[support[i]].pat;
            off_xi = amask & ~X[support[i]].pat;

            if( on_xi == on_f )
              return X[support[i]].sig;
            if( on_xi == off_f )
              return ntk.create_not( X[support[i]].sig );

            if( ( on_xi == amask ) || ( off_xi == amask) ) 
              to_be_deleted.push_back(i);
            else
            { 
              Inew = information( on_xi, off_xi, on_f, off_f );
             
              if( Inew > Imax )
              {
                bidx = i;
                on_x = on_xi;
                off_x = off_xi;
                Imax = Inew;
                bsig = X[support[i]].sig;
              }
            }
          }
          
          if( to_be_deleted.size() > 0 )
          {
            std::reverse(to_be_deleted.begin(), to_be_deleted.end());
            for( uint32_t i = 0; i < to_be_deleted.size() ; ++i )
            {
              support.erase( support.begin()+to_be_deleted[i] ) ;
              uint32_t count = 0;
              if( to_be_deleted[i] < bidx )
                count++;
              bidx -= count;
            }
          }

        if( support.size() == 0 )
          return ntk.get_constant( false );

        if( support.size() <= ps.max_sup )
          return synthesize_leaf( support, X, amask, on_f );

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

        if( ps.try_top_decomposition )
        {
          sim_top_decomposition_fast res = is_top_decomposable_fast( X, pis_support, on_f, amask1, amask0 );
        
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

                if( ps.verbose )
                  std::cout << Fnew << "=" << bsig << " AND " << F1 << std::endl;
                return Fnew;
              }
              case sim_top_decomposition_fast::or_:
              {
                signal<klut_network> F0 = idsd_step( reduced_support, amask0, xmask0 );
                signal<klut_network> Fnew = ntk.create_or( bsig, F0 );
                if( ps.verbose )
                  std::cout << Fnew << "=" << bsig << " OR " << F0 << std::endl;
                return Fnew;
              }
              case sim_top_decomposition_fast::lt_:
              {
                signal<klut_network> F0 = idsd_step( reduced_support, amask0, xmask0 );
                signal<klut_network> Fnew = ntk.create_lt( bsig, F0 );
                if( ps.verbose )
                  std::cout << Fnew << "=" << bsig << "' AND " << F0 << std::endl;
                return Fnew;

              }
              case sim_top_decomposition_fast::le_:
              {  
                signal<klut_network> F1 = idsd_step( reduced_support, amask1, xmask1 );
                signal<klut_network> Fnew = ntk.create_le( bsig, F1 );
                if( ps.verbose )
                  std::cout << Fnew << "=" << bsig << "' OR " << F1 << std::endl;
                return Fnew;
              }
              case sim_top_decomposition_fast::xor_:
              {
                xmask = xmask ^ on_x;
                signal<klut_network> Fxor = idsd_step( reduced_support, amask, xmask );
                signal<klut_network> Fnew = ntk.create_xor( bsig, Fxor );
                if( ps.verbose )
                  std::cout << Fnew << "=" << bsig << " XOR " << Fxor << std::endl;
                return Fnew;
              }
            }
          }
        }
        
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
        for( uint32_t i = 0; i < ntk.sim_patterns.size()-2; ++i )
          support.push_back( i );
        
        TT xmask( Y.pat.num_bits() );
        TT amask = ~xmask;
        return idsd_step( support, amask, xmask );
      }
    
    private:
      simulation_view<Ntk>& ntk;
      sim_decomposition_fast_params & ps;
      kitty::partial_truth_table target;
      sim_pattern<Ntk> Y;
      uint32_t n_bits;
    public:
      std::vector<double> Iactive;
      std::vector<sim_pattern<Ntk>> X;
}; /* class sim_decomposition_fast_impl */
    
} /* namespace detail */




#pragma region sim_decomposition_fast
/*! \brief sim_decomposition_fast algorithm assembles a network bottom up.
 *
 * This method iteratively creates and adds new informative nodes. 
 * It takes an empty network, the simulation patterns of its input nodes and the target functon. 
 * Then, the rest of the network is created using the strategy described in the paper
 * "Learning complex boolean functions: Algorithms and applications." by Oliveira, A., & Sangiovanni-Vincentelli, A. (1993).
 * \param ntk the network wrapped in a simulation view.
 * \param examples the input patterns to be used for assembling the network
 * \param target target simulation pattern
 * \param ps parameters of the sim_decomposition_fast method
 *
 */
template<typename Ntk>
signal<Ntk> sim_decomposition_fast( simulation_view<Ntk>& ntk, std::vector<kitty::partial_truth_table> & examples, kitty::partial_truth_table const & target, sim_decomposition_fast_params & ps, bool re_initialize = true )
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
  //depth_view ntk_dw{ntk};
  detail::sim_decomposition_fast_impl impl( ntk, target, ps );
  auto osignal = impl.run();
  return osignal;
}

#pragma endregion sim_decomposition_fast

} // namespace mockturtle