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
#include "../utils/ccg_net.hpp"
#include <kitty/partial_truth_table.hpp>
#include <kitty/bit_operations.hpp>
#include <kitty/operators.hpp>
#include <kitty/print.hpp>
#include <stdio.h>
#include <stack>

namespace mockturtle
{

namespace ccgame
{

using TT = kitty::partial_truth_table;

struct cusco_cov_ps
{
  /* method */
  int nIters;
  cusco_cov_ps( int nIters ) : nIters( nIters ) {}
};

template<class Ntk>
class cusco_cov
{
public:
    /* problem definition */
    std::vector<TT> X; // input simulations
    std::vector<TT> Y; // output simulations

public:
  /* creation and destruction */
  cusco_cov( std::vector<TT> const&, std::vector<TT> const& );
  ~cusco_cov();
  /* solve */
  Ntk solve_random( cusco_cov_ps const& );
};

/* creation and destruction */
template<class Ntk>
cusco_cov<Ntk>::cusco_cov( std::vector<TT> const& X, std::vector<TT> const& Y ): X(X), Y(Y) {}
template<class Ntk>
cusco_cov<Ntk>::~cusco_cov(){}

template<class Ntk>
Ntk cusco_cov<Ntk>::solve_random( cusco_cov_ps const& ps )
{
  std::random_device rd;  // a seed source for the random number engine
  std::mt19937 gen(5); // mersenne_twister_engine seeded with rd()
  Ntk ntk;
  int nBest = 10000u;
  for( int i{0}; i < ps.nIters; ++i )
  {
    net_t net( X, Y );
//net.print(); // @
    while( net.nHunging > 0 )
    {
      cut_t candidates = net.list_candidate_divs();
//candidates.print(); // @
      bool ClosedSome = net.check_closure( candidates );
      if( net.nHunging == 0 )
        break;

      tab_t table( candidates, net.output_c );
//table.print();
      table.greedy_set_covering( );

      std::uniform_int_distribution<> distrib(0, table.subsets.size()-1);
      std::vector<int> SelIds = table.subsets[distrib(gen)];
      
      cut_t new_c;
      for( int i{0}; i < SelIds.size(); ++i )
        new_c.add_divisor( candidates.divisors[SelIds[i]] );

      if( ClosedSome )
        net.complete_cut( new_c );
      else
        net.add_cut( new_c ); 
//net.print(); // @
    }
//net.print(); // @
  ntk = net.convert<Ntk>();
  printf("|gates| = %d\n", ntk.num_gates());
  }
  return ntk;
}


} // namespace ccgame

} // namespace mockturtle