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
/*! \brief Parameters for sim_decomposition_fastS algorithm */
struct dc_decomposition_fastS_params
{
  bool verbose{false};
  uint32_t max_sup{2};
  bool is_size_aware{false};
  bool try_top_decomposition{true};
  bool try_bottom_decomposition{true};
  bool use_correlation{false};
  bool branch_on_all{true};
  bool try_xor{false};
};

namespace detail
{

  template<typename Ntk>
  class dc_decomposition_fastS_impl
  {
    public:
      using TT = kitty::partial_truth_table;

    public:
      dc_decomposition_fastS_impl( simulation_view<Ntk> & ntk, TT target, dc_decomposition_fastS_params & ps )
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
        //uint32_t max_fanin_size = std::numeric_limits<uint32_t>::max();

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
        }

        if( support.size() == 0 )
          return ntk.get_constant( false );

        if( support.size() <= ps.max_sup )
          return synthesize_leaf( support, X, amask, on_f );
        
        std::vector<uint32_t> reduced_support;
        TT amask1(n_bits);
        TT amask0(n_bits);
        TT xmask1(n_bits);
        TT xmask0(n_bits);

        for( uint32_t bidx = 0; bidx < support.size(); ++bidx )
        {
          reduced_support = support;
          on_x = amask & X[support[bidx]].pat;
          off_x = amask & ~X[support[bidx]].pat;
          bsig = X[support[bidx]].sig;

          amask1 = on_x;
          amask0 = off_x;

          xmask1 = on_x & xmask;
          xmask0 = off_x & xmask;
        

          reduced_support.erase( reduced_support.begin() + bidx );
        
          if( ps.try_top_decomposition )
          {
            sim_top_decomposition_fast res = is_dc_top_decomposable_fast( X, reduced_support, on_f, amask1, amask0, ps.try_xor );
        
            if ( res != sim_top_decomposition_fast::none )
            {
              ntk.clear_network_fanin_size_from_node(ntk.get_node(bsig));
              ntk.update_network_fanin_size();
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
        }

        bidx = 0;
        reduced_support = support;
        reduced_support.erase( reduced_support.begin() + bidx );

        on_x = amask & X[support[bidx]].pat;
        off_x = amask & ~X[support[bidx]].pat;
        bsig = X[support[bidx]].sig;

        amask1 = on_x;
        amask0 = off_x;
        xmask1 = on_x & xmask;
        xmask0 = off_x & xmask;
        
        if( ps.is_size_aware )
        {
          ntk.clear_network_fanin_size_from_node(ntk.get_node(bsig));
          ntk.update_network_fanin_size();        
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
        if( branch_on_all )
        {
          for( uint32_t i = 0; i < ntk.sim_patterns.size()-2; ++i )
            support.push_back( i );
        }
        else
        {
          for( uint32_t i = 0; i < ntk.input_patterns.size(); ++i )
            support.push_back( i );
        }

        TT xmask( Y.pat.num_bits() );
        TT amask = ~xmask;
        return idsd_step( support, amask, xmask );
      }
    
    private:
      simulation_view<Ntk>& ntk;
      dc_decomposition_fastS_params & ps;
      kitty::partial_truth_table target;
      sim_pattern<Ntk> Y;
      uint32_t n_bits;
      bool branch_on_all;
    public:
      std::vector<double> Iactive;
      std::vector<sim_pattern<Ntk>> X;
}; /* class sim_decomposition_fastS_impl */
    
} /* namespace detail */




#pragma region dc_decomposition_fastS
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
signal<Ntk> dc_decomposition_fastS( simulation_view<Ntk>& ntk, std::vector<kitty::partial_truth_table> & examples, kitty::partial_truth_table const & target, dc_decomposition_fastS_params & ps, bool re_initialize = true )
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
  detail::dc_decomposition_fastS_impl impl( ntk, target, ps );
  auto osignal = impl.run();
  return osignal;
}

#pragma endregion dc_decomposition_fastS

#pragma region forest_fastS
/*! \brief dc_decomposition_fastS algorithm assembles a network bottom up.
 *
 * This method iteratively creates and adds new informative nodes. 
 * It takes an empty network, the simulation patterns of its input nodes and the target functon. 
 * Then, the rest of the network is created using the strategy described in the paper
 * "Learning complex boolean functions: Algorithms and applications." by Oliveira, A., & Sangiovanni-Vincentelli, A. (1993).
 * \param ntk the network wrapped in a simulation view.
 * \param examples the input patterns to be used for assembling the network
 * \param target target simulation pattern
 * \param ps parameters of the dc_decomposition_fastS method
 *
 */
/*
template<typename Ntk>
signal<Ntk> forest_fastS( simulation_view<Ntk>& ntk, std::vector<kitty::partial_truth_table> & examples, kitty::partial_truth_table const & target, dc_decomposition_fastS_params & ps, uint32_t ndecs = 3 )
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
    detail::dc_decomposition_fastS_impl impl( ntk, target_p, ps );

    signal<Ntk> osignal = impl.run();
    signal_voters.push_back( osignal ); 
  }

  
  return ntk.create_maj( osignal[0], osignal[1], osignal[2] );
}
*/
#pragma endregion forest_fastS

} // namespace mockturtle