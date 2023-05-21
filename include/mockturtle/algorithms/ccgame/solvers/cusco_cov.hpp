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
#include <mockturtle/algorithms/cleanup.hpp>
#include <stdio.h>
#include <stack>

namespace mockturtle
{

namespace ccgame
{

using TT = kitty::dynamic_truth_table;


template<class Ntk>
struct report_cov_t
{
  int nIt0;
  int nMin;
  int nMax;
  Ntk ntk;
};

struct cusco_cov_ps
{
  /* method */
  /*! \brief number of iterations */
  int nIters;
  /*! \brief capacity: #candidates considered */
  int nCap;
  /*! \brief capacity: #candidates considered */
  bool dcMax;
  cusco_cov_ps( int nIters, int nCap ) : nIters( nIters ), nCap(nCap) { dcMax=false; }
  cusco_cov_ps( int nIters, int nCap, bool dc_maximize ) : nIters( nIters ), nCap(nCap), dcMax(dc_maximize) {}
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
  report_cov_t<Ntk> solve_random( cusco_cov_ps const& );
};

/* creation and destruction */
template<class Ntk>
cusco_cov<Ntk>::cusco_cov( std::vector<TT> const& X, std::vector<TT> const& Y ): X(X), Y(Y) {}
template<class Ntk>
cusco_cov<Ntk>::~cusco_cov(){}

template<class Ntk>
report_cov_t<Ntk> cusco_cov<Ntk>::solve_random( cusco_cov_ps const& ps )
{
  report_cov_t<Ntk> rep;
  rep.nMax = -1;
  rep.nMin = 10000;
  std::random_device rd;  // a seed source for the random number engine
  std::mt19937 gen(5); // mersenne_twister_engine seeded with rd()
  for( int iIt{0}; iIt < ps.nIters; ++iIt )
  {
    Ntk ntk;
    net_t net( X, Y );
    while( net.nHunging > 0 )
    {
      cut_t candidates = net.enumerate_divs();
      cut_t closed_c = net.check_closure( candidates );

      net.add_cut( closed_c );

      if( net.nHunging == 0 )
        break;

      tab_t table( candidates, net.outCut );
      table.greedy_set_covering( ps.nCap );
      //if( ps.dcMax )  
      //  table.select_dc_maximizers( );

      std::uniform_int_distribution<> distrib(0, table.subsets.size()-1);
      int rnum = distrib(gen);
      std::vector<int> SelIds = table.subsets[rnum];
      
      cut_t new_c;
      for( int i{0}; i < SelIds.size(); ++i )
        new_c.add_node( candidates.nodes[SelIds[i]] );

      net.complete_cut( new_c );
       
    }
    ntk = net.convert<Ntk>();
    ntk = cleanup_dangling( ntk );
    if( iIt == 0 )  rep.nIt0 = ntk.num_gates();
    if( ntk.num_gates() < rep.nMin )
    {
      rep.ntk = ntk;
      rep.nMin = ntk.num_gates();
    }
    else if ( ntk.num_gates() > rep.nMax )  rep.nMax = ntk.num_gates();
  }
  //printf("END\n", ntk_best.num_gates());
  //printf("|gates*|= %d\n", ntk_best.num_gates());

  return rep;
}


} // namespace ccgame

} // namespace mockturtle