/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
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
  \file boolean_chain.hpp
  \brief data structure for storing boolean function representations

  \author Andrea Costamagna
*/
#pragma once
#include "../../../networks/xag.hpp"
#include <kitty/partial_truth_table.hpp>
#include <kitty/bit_operations.hpp>
#include <kitty/operators.hpp>
#include <kitty/constructors.hpp>
#include <kitty/print.hpp>
#include <stdio.h>
#include <stack>

namespace mockturtle
{

namespace ccgame
{

using TT = kitty::partial_truth_table;

struct cusco_rem_ps
{
  /* method */
  int nIters;
  cusco_rem_ps( int nIters ) : nIters( nIters ) {}
};

template<class Ntk>
class cusco_rem
{
public:
    /* problem definition */
    std::vector<TT> X; // input simulations
    std::vector<TT> Y; // output simulations

public:
  /* creation and destruction */
  cusco_rem( std::vector<TT> const&, std::vector<TT> const& );
  ~cusco_rem();
  /* solve */
  Ntk solve_random( cusco_rem_ps const& );
};

/* creation and destruction */
template<class Ntk>
cusco_rem<Ntk>::cusco_rem( std::vector<TT> const& X, std::vector<TT> const& Y ): X(X), Y(Y) {}
template<class Ntk>
cusco_rem<Ntk>::~cusco_rem(){}

template<class Ntk>
Ntk cusco_rem<Ntk>::solve_random( cusco_rem_ps const& ps )
{
  std::random_device rd;  // a seed source for the random number engine
  std::mt19937 gen(5); // mersenne_twister_engine seeded with rd()
  Ntk ntk;
  Ntk ntk_best;
  int nVars = ceil(log2(X[0].num_bits()));
  DTT func( nVars );
  DTT mask( nVars );
  kitty::create_from_binary_string( func, kitty::to_binary( Y[0] ) );
  kitty::create_from_binary_string( mask, kitty::to_binary( Y[0] | ~Y[0] ) );
  int nBest = 10000u;
  std::vector<DTT> xs;
  for( int i{0}; i < nVars; ++i )
  {
    xs.emplace_back( nVars );
    kitty::create_nth_var( xs[i], i );
  }
      analyzer_t analyzer;

  for( int i{0}; i < ps.nIters; ++i )
  {
    net_t net( X, Y );
    net.print();
    net.cuts[net.nCuts-1].set_func( func );
    net.cuts[net.nCuts-1].set_mask( mask );
    int idBound = 2;
    int bestRwd = -1;
    symmetry_t bestSym;
    //while( net.nHunging > 0 )
    {
      std::vector<symmetry_t> candidates = net.symmetry_analysis( xs, idBound );
      for( int i{0}; i<candidates.size(); ++i )
      {
        if( candidates[i].rwd > bestRwd )
        {
          bestSym = candidates[i];
          bestRwd = candidates[i].rwd;
        }
      }
      analyzer.print_symmetries( {bestSym} );
      //cut_t closed_c = net.check_closure( candidates );
      //net.add_cut( closed_c );

      if( net.nHunging == 0 )
        break;

      //tab_t table( candidates, net.outCut );
      //table.greedy_set_covering( );

      //std::uniform_int_distribution<> distrib(0, table.subsets.size()-1);
      //std::vector<int> SelIds = table.subsets[distrib(gen)];
      
      //cut_t new_c;
      //for( int i{0}; i < SelIds.size(); ++i )
      //  new_c.add_node( candidates.nodes[SelIds[i]] );

      //net.complete_cut( new_c );
       
    }
    ntk = net.convert<Ntk>();
    if( ntk.num_gates() < nBest )
    {
      nBest = ntk.num_gates();
      ntk_best = ntk;
    }
    printf("|gates| = %d\n", ntk.num_gates());
  }
  printf("END\n", ntk.num_gates());
  printf("|gates*|= %d\n", ntk_best.num_gates());

  return ntk_best;
}


} // namespace ccgame

} // namespace mockturtle