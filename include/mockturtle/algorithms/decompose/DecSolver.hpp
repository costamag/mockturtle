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
  \file DecSims.hpp
  \brief data structure for synthesizing an network

  \author Andrea Costamagna
*/

#pragma once

#include <stdio.h>
#include <kitty/print.hpp>
#include "DecNet.hpp"
#include "DecAnalyzer.hpp"

namespace mockturtle
{

template<class TT, class Ntk>
class DecSolver
{
private:
    std::vector<TT>  vTruths;
    std::vector<TT>  vMasks;
    /* solver view */
    std::vector<signal_t> X; // ( remapped ) inputs signals |X| = n
    std::vector<int>      V;
    std::vector<signal_t> Y; // targets signals             |Y| = m

public:
    DecSolver( const std::vector<TT>&, const std::vector<TT>& );
    ~DecSolver();
    /* solve */
    Ntk solve();
    /* visualize */
    void PrintSpecs();

};

template<class TT, class Ntk>
DecSolver<TT, Ntk>::DecSolver( const std::vector<TT>& vTruths, const std::vector<TT>& vMasks ) :
vTruths(vTruths),
vMasks( vMasks )
{
}

template<class TT, class Ntk>
DecSolver<TT, Ntk>::~DecSolver()
{
}

#pragma region solve

template<class TT, class Ntk>
Ntk DecSolver<TT, Ntk>::solve()
{
  Ntk ntk;
  Ntk ntkBest;
  int nBest = 100000;
  for( int it{0}; it<1; ++it )
  {
    /* init */
    DecNet<TT, Ntk>  net;
    net.init( vTruths, vMasks );
    X = net.getPIs();
    for( int i{0}; i < X.size(); ++i )  V.push_back( i );
    Y = net.getTargets();

    /* solve */
    //while( net.numTargets() > 0 )
    {
      DecAnalyzer<TT, Ntk> checker( &net, &X, &Y, &V );
      checker.check1();
      checker.check2();
      std::vector<action_t<TT>> T1 = checker.get_topdec();
      std::vector<action_t<TT>> RM = checker.get_remove();
      checker.print_actions( T1 );
      checker.print_actions( RM );
      std::vector<action_t<TT>> RP = checker.get_remap(); 
      checker.print_actions( RP );

    }
    /* convert result */
    //DecChsToGraph<TT, Ntk> conv( net );  
    //ntk = conv.convert();
    //if( ntk.num_gates() < nBest )
    //    ntkBest = ntk;
  }
  return ntk;
}

#pragma endregion solve

#pragma region visualize
template<class TT, class Ntk>
void DecSolver<TT, Ntk>::PrintSpecs()
{
    printf("TRUTHS:\n");
    for( int i{0}; i<vTruths.size(); ++i ){ printf("%d ", i ); kitty::print_binary( vTruths[i] ); printf("\n");}
    printf("MASKS:\n");
    for( int i{0}; i<vMasks.size(); ++i ){ printf("%d ", i ); kitty::print_binary( vMasks[i] ); printf("\n");}
}
#pragma endregion visualize


} // namespace mockturtle