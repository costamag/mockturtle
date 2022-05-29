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
#include "sim_decomposition_checks.hpp"
#include "sim_operations.hpp"
#include "sim_utils.hpp"


#include <fmt/format.h>
#include <kitty/constructors.hpp>
#include <kitty/decomposition.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/statistics.hpp>
#include <kitty/operators.hpp>
#include <kitty/print.hpp>

namespace mockturtle
{
/*! \brief Parameters for sim_decomposition algorithm */
struct sim_decomposition_params
{
  bool verbose{false};
  size_t max_sup{2};
  bool is_informed{true};
  bool try_top_decomposition{true};
  bool try_bottom_decomposition{false};
};

namespace detail
{

  template<typename Ntk>
  class sim_decomposition_impl
  {
    public:
      using TT = typename kitty::partial_truth_table;

    public:
      sim_decomposition_impl( simulation_view<Ntk> & ntk, TT target, sim_decomposition_params & ps )
      : ntk(ntk),
        target(target),
        ps(ps),
        F(target)
      {
      }
    
      signal<Ntk> idsd_step( std::vector<sim_pattern<Ntk>>& X, sim_pattern<Ntk> & Y )
      {
        
        if( X.size() == 0 )
          return ntk.get_constant( false );    
        if( X[0].pat.num_bits() == 0 )
          return ntk.get_constant( false );

        assert( ( X[0].pat.num_bits() == Y.pat.num_bits() ) );

        if( kitty::count_ones(Y.pat ) == 0 ) // contraddiction
          return ntk.get_constant( false );
        else if( kitty::count_ones(Y.pat) == Y.pat.num_bits() ) // tautology
          return ntk.get_constant( true );


        uint32_t bidx = 0;
        sim_pattern<Ntk> bpat;
        uint32_t idx_min = 0;
        double Inew = 0;
        double Imax = 0;
        std::vector<double> mi_vect;
        std::vector<uint32_t> idx_vect;

        std::vector<size_t> to_be_deleted;

        if( ps.is_informed )
        {
          for( size_t i = 0; i < X.size(); ++i )
          {
            if( ( kitty::count_ones(X[i].pat) == X[i].pat.num_bits() ) || ( kitty::count_ones(X[i].pat) == 0 ) )
            {
              to_be_deleted.push_back(i);
            }
            else
            {
              Inew = kitty::mutual_information( std::vector{&(X[i].pat)}, &(Y.pat) );
              mi_vect.push_back( Inew );
              idx_vect.push_back( i-to_be_deleted.size() );
              if( Inew >= Imax )
              {
                bidx = i;
                bpat = X[i];
                Imax = Inew;
              }
            }
          }
          if( to_be_deleted.size() > 0 )
          {
            std::reverse(to_be_deleted.begin(), to_be_deleted.end());
            
            size_t count = 0;
            for( size_t i = 0; i < to_be_deleted.size() ; ++i )
            {
              if( to_be_deleted[i] < bidx )
                count++;
              X.erase( X.begin()+to_be_deleted[i] ) ;
            }
            bidx -= count;
          }
        }

        if( X.size() <= ps.max_sup )
        {
          std::vector<TT*> ipatterns;
          TT * opatterns = & Y.pat;
          std::vector<signal<Ntk>> children;
          for( size_t i = 0; i < X.size(); ++i )
          {
            children.push_back( X[i].sig );
            ipatterns.push_back( &(X[i].pat) );
          }

          chj_result chj_res = chatterjee_method( ipatterns, &(Y.pat) );
          signal<Ntk> fc = ntk.create_node( children, chj_res.dtt );
          return fc;
        }

        std::pair< std::vector<sim_pattern<Ntk>>, sim_pattern<Ntk>> XY0 = compute_cofactor0( X, Y, bidx );
        std::pair< std::vector<sim_pattern<Ntk>>, sim_pattern<Ntk>> XY1 = compute_cofactor1( X, Y, bidx ); 

        if( ps.try_top_decomposition )
        {
          sim_top_decomposition res = is_top_decomposable( XY0, XY1 );
          if ( res != sim_top_decomposition::none )
          {
            switch ( res )
            {
              default:
                assert( false );
              case sim_top_decomposition::and_:
              {
                signal<klut_network> F1 = idsd_step( XY1.first, XY1.second );
                return ntk.create_and( bpat.sig, F1 );
              }
              case sim_top_decomposition::or_:
              {
                signal<klut_network> F0 = idsd_step( XY0.first, XY0.second );
                return ntk.create_or( bpat.sig, F0 );
              }
              case sim_top_decomposition::lt_:
              {
                signal<klut_network> F0 = idsd_step( XY0.first, XY0.second );
                return ntk.create_lt( bpat.sig, F0 );
              }
              case sim_top_decomposition::le_:
              {  
                signal<klut_network> F1 = idsd_step( XY1.first, XY1.second );
                return ntk.create_le( bpat.sig, F1 );
              }
              case sim_top_decomposition::xor_:
              {
                remove_column_and_invert( X, Y, bidx ); 
                return ntk.create_xor( bpat.sig, idsd_step( X, Y ) );
              }
            }
          }
        }
        if( ps.try_bottom_decomposition )
        {
          bottom_res<Ntk> bres = is_bottom_decomposable( X, Y, Imax, mi_vect, idx_vect );
          if( bres.found )
          {
            std::cout << "+" ; 
            ntk.create_node( bres.children, bres.chj.dtt );
            return idsd_step( X, Y );
          }
        }
        
        signal<klut_network> F0 = idsd_step( XY0.first, XY0.second );
        signal<klut_network> f0 = ntk.create_and( ntk.create_not( bpat.sig ), F0 );

        signal<klut_network> F1 = idsd_step( XY1.first, XY1.second );
        signal<klut_network> f1 = ntk.create_and( bpat.sig, F1 );

        return ntk.create_or( f1, f0 );
      }

      signal<Ntk> run()
      {
        std::vector<sim_pattern<Ntk>> X = ntk.sim_patterns;

        X.erase( X.begin() );
        X.erase( X.begin() );

        sim_pattern<Ntk> Y = F;

        return idsd_step( X, Y );
      }
    
    private:
      simulation_view<Ntk>& ntk;
      sim_decomposition_params & ps;
      kitty::partial_truth_table target;
    public:
      sim_pattern<Ntk> F; 
      std::vector<double> Iactive;
}; /* class sim_decomposition_impl */
    
} /* namespace detail */




#pragma region sim_decomposition
/*! \brief sim_decomposition algorithm assembles a network bottom up.
 *
 * This method iteratively creates and adds new informative nodes. 
 * It takes an empty network, the simulation patterns of its input nodes and the target functon. 
 * Then, the rest of the network is created using the strategy described in the paper
 * "Learning complex boolean functions: Algorithms and applications." by Oliveira, A., & Sangiovanni-Vincentelli, A. (1993).
 * \param ntk the network wrapped in a simulation view.
 * \param examples the input patterns to be used for assembling the network
 * \param target target simulation pattern
 * \param ps parameters of the sim_decomposition method
 *
 */
template<typename Ntk>
signal<Ntk> sim_decomposition( simulation_view<Ntk>& ntk, std::vector<kitty::partial_truth_table> & examples, kitty::partial_truth_table const & target, sim_decomposition_params & ps )
{
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
  //depth_view ntk_dw{ntk};
  detail::sim_decomposition_impl impl( ntk, target, ps );
  auto osignal = impl.run();
  return osignal;
}

template<typename Ntk>
std::vector<signal<Ntk>> sim_decomposition( simulation_view<Ntk>& ntk, std::vector<kitty::partial_truth_table> & examples, std::vector<kitty::partial_truth_table> const & targets, sim_decomposition_params & ps )
{
  std::vector<signal<Ntk>> osignals;
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
    for( auto y : targets )
    {
      std::cout << "  "; kitty::print_binary( y ); std::cout << std::endl;
    }
    std::cout << "  ";
    for( size_t i = 0; i < ntk.sim_patterns[0].pat.num_bits(); ++i )
      std::cout << "-";
    std::cout << std::endl;
  }

  for( uint32_t i = 0; i < targets.size(); ++i )
  {
    detail::sim_decomposition_impl impl( ntk, targets[i], ps );
    osignals.push_back(impl.run());
  }

  return osignals;
}

#pragma endregion sim_decomposition

} // namespace mockturtle