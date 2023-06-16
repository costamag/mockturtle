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
#include "cusco_cov.hpp"
#include "cusco_rem.hpp"
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

using TT = kitty::dynamic_truth_table;

/*! \brief methods implemented
 */
enum solver_t
{
  _SYM_1SH,
  _SYM_RND,
  _SYM_ENT,
  _COV_RND,
  _COV_DCM,
  _COV_MCTS,
  _COV_GEN
};

template<class Ntk>
struct report_t
{
  int nIt0;
  int nMin;
  int nMax;
  Ntk ntk;
  double time;
};

struct cusco_ps
{
  /* method */
  /*! \brief solver type */
  solver_t type;
  /*! \brief number of iterations */
  int nIters;
  /*! \brief capacity [only for covering -1 to let unbounded] */
  int nCap;

  cusco_ps( solver_t type, int nIters ) : type( type ), nIters( nIters ) 
  {
    nCap=-1;
  }
  cusco_ps( solver_t type, int nIters, int nCap ) : type( type ), nIters( nIters ), nCap(nCap) {}
};

template<class Ntk>
class cusco
{
public:
    /* problem definition */
    std::vector<TT> X; // input simulations
    std::vector<TT> Y; // output simulations

public:
  /* creation and destruction */
  cusco( std::vector<TT> const&, std::vector<TT> const& );
  ~cusco();
  /* solve */
  report_t<Ntk> solve( cusco_ps const& );
};

/* creation and destruction */
template<class Ntk>
cusco<Ntk>::cusco( std::vector<TT> const& X, std::vector<TT> const& Y ): X(X), Y(Y) {}
template<class Ntk>
cusco<Ntk>::~cusco(){}

template<class Ntk>
report_t<Ntk> cusco<Ntk>::solve( cusco_ps const& ps )
{
  std::clock_t start;
  double duration;
  start = std::clock();
  report_t<Ntk> rp;

  switch ( ps.type )
  {
    case _SYM_1SH : 
    {
      assert( Y.size() == 1 );
      cusco_rem<Ntk> solver0( X, Y );
      cusco_rem_ps ps0( 1 );
      report_rem_t<Ntk> rp0 = solver0.solve_1shot( ps0 );
      rp.nIt0 = rp0.nIt0;
      rp.nMin = rp0.nIt0;
      rp.nMax = rp0.nIt0;
      rp.ntk  = rp0.ntk;
      break;
    }
    case _SYM_RND : 
    {
      assert( Y.size() == 1 );
      cusco_rem<Ntk> solver1( X, Y );
      cusco_rem_ps ps1( ps.nIters );
      report_rem_t<Ntk> rp1 = solver1.solve_random( ps1 );
      rp.nIt0 = rp1.nIt0;
      rp.nMin = rp1.nMin;
      rp.nMax = rp1.nMax;
      rp.ntk  = rp1.ntk;
      break;
    }
    case _SYM_ENT: 
    {
      assert( Y.size() == 1 );
      cusco_rem<Ntk> solver1( X, Y );
      cusco_rem_ps ps1( ps.nIters );
      report_rem_t<Ntk> rp1 = solver1.solve_entropic( ps1 );
      rp.nIt0 = rp1.nIt0;
      rp.nMin = rp1.nMin;
      rp.nMax = rp1.nMax;
      rp.ntk  = rp1.ntk;
      break;
    }
    case _COV_RND :
    {
      cusco_cov<Ntk> solver2( X, Y );
      cusco_cov_ps ps2( ps.nIters, ps.nCap );
      report_cov_t<Ntk> rp2 = solver2.solve_random( ps2 );
      rp.nIt0 = rp2.nIt0;
      rp.nMin = rp2.nMin;
      rp.nMax = rp2.nMax;
      rp.ntk  = rp2.ntk;
      break;
    }
    case _COV_DCM :
    {
      cusco_cov<Ntk> solver3( X, Y );
      cusco_cov_ps ps3( ps.nIters, ps.nCap, true );
      report_cov_t<Ntk> rp3 = solver3.solve_random( ps3 );
      rp.nIt0 = rp3.nIt0;
      rp.nMin = rp3.nMin;
      rp.nMax = rp3.nMax;
      rp.ntk  = rp3.ntk;
      break;
    }
    case _COV_MCTS :
    {
      cusco_cov<Ntk> solver4( X, Y );
      cusco_cov_ps ps4( ps.nIters, ps.nCap, false );
      report_cov_t<Ntk> rp4 = solver4.solve_mcts( ps4 );
      rp.nIt0 = rp4.nIt0;
      rp.nMin = rp4.nMin;
      rp.nMax = rp4.nMax;
      rp.ntk  = rp4.ntk;
      printf("%d %d %d %d\n", rp.nIt0, rp.nMin, rp.nMax, rp.ntk );
      break;
    }
    case _COV_GEN :
    {
      cusco_cov<Ntk> solver5( X, Y );
      cusco_cov_ps ps5( ps.nIters, ps.nCap, false );
      report_cov_t<Ntk> rp5 = solver5.solve_genetic( ps5 );
      rp.nIt0 = rp5.nIt0;
      rp.nMin = rp5.nMin;
      rp.nMax = rp5.nMax;
      rp.ntk  = rp5.ntk;
      break;
    }
  }
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  rp.time = duration;
  //printf("SUMMARY:\n");
  //printf( "ngates = %d  time = %.2f\n", ntk.num_gates(), duration );

  return rp;
}


} // namespace ccgame

} // namespace mockturtle