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
/*! \brief Parameters for sim_decomposition algorithm */
struct fgsim_decomposition_params
{
  bool verbose{false};
  size_t max_sup{2};
  bool is_informed{true};
  bool try_top_decomposition{true};
};

namespace detail
{

  template<typename Ntk>
  class fgsim_decomposition_impl
  {
    public:
      using TT = typename kitty::partial_truth_table;

    public:
      fgsim_decomposition_impl( simulation_view<Ntk> & ntk, TT target, TT global, fgsim_decomposition_params & ps )
      : ntk(ntk),
        target(target),
        global(global),
        ps(ps),
        F(target)
      {
      }
    
      signal<Ntk> idsd_step( std::vector<sim_pattern<Ntk>>& X, sim_pattern<Ntk> & Y, sim_pattern<Ntk> & G )
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
        double I1 = 0;
        double I2 = 0;
        double Imax = 0;
        std::vector<double> mi_vect;
        int64_t Nbits = int64_t(X[0].pat.num_bits());

        std::vector<uint32_t> idx_vect;

        std::vector<size_t> to_be_deleted;
          //kitty::print_binary( Y.pat ); std::cout << std::endl;
          for( size_t i = 0; i < X.size(); ++i )
          {
            if( ( kitty::count_ones(X[i].pat) == Nbits ) || ( kitty::count_ones(X[i].pat) == 0 ) )
            {
              to_be_deleted.push_back(i);
            }
            else
            {
              //Inew = std::abs(kitty::count_ones( X[i].pat&G.pat )-(double)kitty::count_ones( X[i].pat )*kitty::count_ones( G.pat )/Nbits); 
              Inew = std::abs(Nbits-int64_t(2*kitty::count_ones( X[i].pat^G.pat) )); 
              //Inew = std::max( I1, I2 );
              mi_vect.push_back( Inew );
            }

            idx_vect.push_back( i-to_be_deleted.size() );
            if( Inew >= Imax )
            {
              bidx = i;
              bpat = X[i];
              Imax = Inew;
            }
          }
          //std::cout << "-------------------- END --------------------" << std::endl;
          
          if( to_be_deleted.size() > 0 )
          {
            std::reverse(to_be_deleted.begin(), to_be_deleted.end());
            
            size_t count = 0;
            for( size_t i = 0; i < to_be_deleted.size() ; ++i )
            {
              if( to_be_deleted[i] <= bidx )
                count++;
              X.erase( X.begin()+to_be_deleted[i] ) ;
            }
            bidx -= count;
          }
          
        if( X.size() == 0 )
          return ntk.get_constant( false );  

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

        //std::cout << X.size() << " " << bidx << std::endl;
        std::pair< std::vector<sim_pattern<Ntk>>, sim_pattern<Ntk>> XY0 = compute_cofactor0( X, Y, bidx );
        std::pair< std::vector<sim_pattern<Ntk>>, sim_pattern<Ntk>> XY1 = compute_cofactor1( X, Y, bidx );
        std::pair< std::vector<sim_pattern<Ntk>>, sim_pattern<Ntk>> G0 = compute_cofactor0( X, G, bidx );
        std::pair< std::vector<sim_pattern<Ntk>>, sim_pattern<Ntk>> G1 = compute_cofactor1( X, G, bidx ); 


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
                signal<klut_network> F1 = idsd_step( XY1.first, XY1.second, G1.second );
                return ntk.create_and( bpat.sig, F1 );
              }
              case sim_top_decomposition::or_:
              {
                signal<klut_network> F0 = idsd_step( XY0.first, XY0.second, G0.second );
                return ntk.create_or( bpat.sig, F0 );
              }
              case sim_top_decomposition::lt_:
              {
                signal<klut_network> F0 = idsd_step( XY0.first, XY0.second, G0.second );
                return ntk.create_lt( bpat.sig, F0 );
              }
              case sim_top_decomposition::le_:
              {  
                signal<klut_network> F1 = idsd_step( XY1.first, XY1.second, G1.second );
                return ntk.create_le( bpat.sig, F1 );
              }
              case sim_top_decomposition::xor_:
              {
                auto Xf = X;
                remove_column_and_invert( X, Y, bidx );
                remove_column_and_invert( Xf, G, bidx );

                return ntk.create_xor( bpat.sig, idsd_step( X, Y, G ) );
              }
            }
          }
        }
        
        signal<klut_network> F0 = idsd_step( XY0.first, XY0.second, G0.second );
        signal<klut_network> f0 = ntk.create_and( ntk.create_not( bpat.sig ), F0 );

        signal<klut_network> F1 = idsd_step( XY1.first, XY1.second, G1.second );
        signal<klut_network> f1 = ntk.create_and( bpat.sig, F1 );

        return ntk.create_or( f1, f0 );
      }

      signal<Ntk> run()
      {
        std::vector<sim_pattern<Ntk>> X = ntk.sim_patterns;

        X.erase( X.begin() );
        X.erase( X.begin() );

        sim_pattern<Ntk> Y = F;
        sim_pattern<Ntk> G = global;

        return idsd_step( X, Y, G );
      }
    
    private:
      simulation_view<Ntk>& ntk;
      fgsim_decomposition_params & ps;
      kitty::partial_truth_table target;
    public:
      sim_pattern<Ntk> F; 
      sim_pattern<Ntk> global; 
      std::vector<double> Iactive;
}; /* class sim_decomposition_impl */
    
} /* namespace detail */




#pragma region fgsim_decomposition
/*! \brief fgsim_decomposition algorithm assembles a network bottom up.
 *
 * This method iteratively creates and adds new informative nodes. 
 * It takes an empty network, the simulation patterns of its input nodes and the target functon. 
 * Then, the rest of the network is created using the strategy described in the paper
 * "Learning complex boolean functions: Algorithms and applications." by Oliveira, A., & Sangiovanni-Vincentelli, A. (1993).
 * \param ntk the network wrapped in a simulation view.
 * \param examples the input patterns to be used for assembling the network
 * \param target target simulation pattern
 * \param ps parameters of the fgsim_decomposition method
 *
 */
template<typename Ntk>
signal<Ntk> fgsim_decomposition( simulation_view<Ntk>& ntk, std::vector<kitty::partial_truth_table> & examples, kitty::partial_truth_table const & target, fgsim_decomposition_params & ps )
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
  detail::fgsim_decomposition_impl impl( ntk, target, target, ps );
  auto osignal = impl.run();
  return osignal;
}

template<typename Ntk>
std::vector<signal<Ntk>> fgsim_decomposition( simulation_view<Ntk>& ntk, std::vector<kitty::partial_truth_table> & examples, std::vector<kitty::partial_truth_table> const & targets, fgsim_decomposition_params & ps )
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
  kitty::partial_truth_table global;
  kitty::partial_truth_table tmp;

  for( uint32_t i = 0; i < targets.size(); ++i )
  {
    global = targets[i];
    double Best = 0;
    double New = 0;
    for( uint32_t j = i+1; j < targets.size(); ++j )
    {
      New = std::abs( int64_t(targets[i].num_bits())-int64_t(2*kitty::count_ones( targets[i]^targets[j]) )); 
      if( New >= Best )
      {
        Best = New;
        tmp = targets[j];
      }
    }
      
    global^=tmp;
    
    detail::fgsim_decomposition_impl impl( ntk, targets[i],global, ps );
    osignals.push_back(impl.run());
  }

  return osignals;
}

#pragma endregion fgsim_decomposition

} // namespace mockturtle