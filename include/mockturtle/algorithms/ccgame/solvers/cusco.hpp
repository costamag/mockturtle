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
  _SYM_1DE,
  _SYM_RND,
  _SYM_RDE,
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
  signal<Ntk> osig;
  double levels{0u};
  double time;
  double area;
  bool Esl{false};
  void set_ntk( Ntk ntk_new ){ ntk=ntk_new; };

  void print()
  {
    printf( "nIt0=%d nMin=%d nMax=%d ntk.size()=%d time=%f\n", nIt0, nMin, nMax, ntk.size() );
  }
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
  /*! \brief input arrival patterns */
  std::vector<double> T;

  library__t lib;

  cusco_ps( solver_t type, int nIters ) : type( type ), nIters( nIters ) 
  {
    nCap=-1;
  }
  cusco_ps( solver_t type, int nIters, int nCap ) : type( type ), nIters( nIters ), nCap(nCap) {}
  cusco_ps( solver_t type, int nIters, library__t Lib ) : type( type ), nIters( nIters ), lib(Lib) {}
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
  report_t<Ntk> solve( cusco_ps const&, std::vector<signal<Ntk>>, Ntk * );
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
    case _SYM_1DE : 
    {
      assert( Y.size() == 1 );
      cusco_rem<Ntk> solver0( X, Y );
      cusco_rem_ps ps0_de( 1, ps.lib );
      ps0_de.T = ps.T;
      assert( ps.T.size() == X.size() );
      report_rem_t<Ntk> rp0_de = solver0.solve_1delay( ps0_de );
      rp.nIt0 = rp0_de.nIt0;
      rp.nMin = rp0_de.nIt0;
      rp.nMax = rp0_de.nIt0;
      rp.ntk  = rp0_de.ntk;
      rp.Esl = rp0_de.E_solution;
      rp.levels = rp0_de.levels;
      rp.area = rp0_de.area;
      break;
    }
    case _SYM_RDE : 
    {
      assert( Y.size() == 1 );
      cusco_rem<Ntk> solverR( X, Y );
      cusco_rem_ps psR_de( ps.nIters, ps.lib );
      psR_de.T = ps.T;
      assert( ps.T.size() == X.size() );
      report_rem_t<Ntk> rpR_de = solverR.solve_Rdelay( psR_de );
      rp.nIt0 = rpR_de.nIt0;
      rp.nMin = rpR_de.nIt0;
      rp.nMax = rpR_de.nIt0;
      rp.ntk  = rpR_de.ntk;
      rp.Esl = rpR_de.E_solution;
      rp.levels = rpR_de.levels;
      rp.area = rpR_de.area;

      break;
    }
  }
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  rp.time = duration;

  return rp;
}


template<class Ntk>
report_t<Ntk> cusco<Ntk>::solve( cusco_ps const& ps, std::vector<signal<Ntk>> inSigs, Ntk * pNtk )
{
  std::clock_t start;
  double duration;
  start = std::clock();
  report_t<Ntk> rp;

  switch ( ps.type )
  {
    case _SYM_1DE : 
    {
      assert( Y.size() == 1 );
      cusco_rem<Ntk> solver0( X, Y );
      cusco_rem_ps ps0_de( 1, ps.lib );
      ps0_de.T = ps.T;
      assert( ps.T.size() == X.size() );
      report_rem_t<Ntk> rp0_de = solver0.solve_1delay( ps0_de, pNtk, inSigs );
      rp.nIt0 = rp0_de.nIt0;
      rp.nMin = rp0_de.nIt0;
      rp.nMax = rp0_de.nIt0;
      rp.ntk  = rp0_de.ntk;
      rp.Esl = rp0_de.E_solution;
      rp.levels = rp0_de.levels;
      rp.area = rp0_de.area;
      rp.osig = rp0_de.osig;

      break;
    }
    case _SYM_RDE : 
    {
      assert( Y.size() == 1 );
      cusco_rem<Ntk> solver1( X, Y );
      cusco_rem_ps ps1_de( ps.nIters, ps.lib );
      ps1_de.T = ps.T;
      assert( ps.T.size() == X.size() );
      report_rem_t<Ntk> rp1_de = solver1.solve_Rdelay( ps1_de, pNtk, inSigs );
      rp.nIt0 = rp1_de.nIt0;
      rp.nMin = rp1_de.nIt0;
      rp.nMax = rp1_de.nIt0;
      rp.ntk  = rp1_de.ntk;
      rp.Esl  = rp1_de.E_solution;
      rp.levels = rp1_de.levels;
      rp.area = rp1_de.area;
      rp.osig = rp1_de.osig;
      break;
    }
  }
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  rp.time = duration;

  return rp;
}

} // namespace ccgame

} // namespace mockturtle