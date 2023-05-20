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

using TT = kitty::partial_truth_table;

/*! \brief methods implemented
 */
enum solver_t
{
  _SYM_RND,
  _COV_RND
};

struct cusco_ps
{
  /* method */
  solver_t type;
  int nIters;

  cusco_ps( solver_t type, int nIters ) : type( type ), nIters( nIters ) {}
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
  Ntk solve( cusco_ps const& );
};

/* creation and destruction */
template<class Ntk>
cusco<Ntk>::cusco( std::vector<TT> const& X, std::vector<TT> const& Y ): X(X), Y(Y) {}
template<class Ntk>
cusco<Ntk>::~cusco(){}

template<class Ntk>
Ntk cusco<Ntk>::solve( cusco_ps const& ps )
{
  std::clock_t start;
  double duration;
  start = std::clock();
  Ntk ntk;
  switch ( ps.type )
  {
    case _SYM_RND : 
    {
      assert( Y.size() == 1 );
      cusco_rem<Ntk> solver0( X, Y );
      cusco_rem_ps ps0( ps.nIters );
      ntk = solver0.solve_random( ps0 );
      break;
    }
    case _COV_RND :
    {
      cusco_cov<Ntk> solver1( X, Y );
      cusco_cov_ps ps1( ps.nIters );
      ntk = solver1.solve_random( ps1 );
      break;
    }
  }
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

  printf("SUMMARY:\n");
  printf( "ngates = %d  time = %.2f\n", ntk.num_gates(), duration );

  return ntk;
}


} // namespace ccgame

} // namespace mockturtle