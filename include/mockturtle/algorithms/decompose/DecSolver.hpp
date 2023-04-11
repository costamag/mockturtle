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
#include "DecChsToGraph.hpp"

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

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
    Ntk man_sym_solve();
    void remap( DecNet<TT, Ntk> *, action_t<TT> );
    void close( DecNet<TT, Ntk> *, std::vector<action_t<TT>> );
    /* visualize */
    void PrintSpecs();
    void show_state( DecNet<TT, Ntk> *, std::vector<signal_t> * );

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
Ntk DecSolver<TT, Ntk>::man_sym_solve()
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
    while( Y.size() > 0 )
    {
      show_state( &net, &Y );

      DecAnalyzer<TT, Ntk> checker( &net, &X, &Y, &V );

      checker.check0();
      std::vector<action_t<TT>> CS = checker.get_closure();
      if( CS.size() > 0 )
        close( &net, CS );
      printf("|CS|=%d\n", CS.size());
      if( Y.size() > 0 )
      {
        checker.check2();
        std::vector<action_t<TT>> RP = checker.get_remap();
        checker.print_actions( RP );
        int iEnd = RP.size()-1;
        printf( "Choose the move[%2d-%3d]:", 0, iEnd  );
        int MV;
        std::cin >> MV;
        remap( &net, RP[MV] );
      }
    }
    /* convert result */
    DecChsToGraph<TT, Ntk> conv( net );  
    ntk = conv.convert();
    //if( ntk.num_gates() < nBest )
    //    ntkBest = ntk;
  }
  return ntk;
}

#pragma endregion solve

#pragma region closure
template<class TT, class Ntk>
void DecSolver<TT, Ntk>::close( DecNet<TT, Ntk> * pNet, std::vector<action_t<TT>> actions )
{
  printf("i\n");
  std::vector<int> rm_targs;
  for( int i{0}; i < Y.size(); ++i )
    rm_targs.push_back(0);

  for( auto act : actions )
  {
    if( rm_targs[act.sigs[0]] == 0 )
    {
      if( act.type == DecAct_t::BUF_ )
        pNet->close_target( Y[act.sigs[0]], X[act.sigs[1]], 0 );
      else if( act.type == DecAct_t::INV_ )
        pNet->close_target( Y[act.sigs[0]], X[act.sigs[1]], 1 );
      printf("a\n");
        rm_targs[act.sigs[0]] = 1;
    }
  }
  for( int i{Y.size()-1}; i >= 0; --i )
  {
    if( rm_targs[i] == 1 )
      Y.erase( Y.begin() + i );  
  }
  printf("e\n");
  
}
#pragma endregion closure

#pragma region remap
template<class TT, class Ntk>
 void DecSolver<TT, Ntk>::remap( DecNet<TT, Ntk> * pNet, action_t<TT> act )
  {
    pNet->change_sim_info( Y[act.sigs[0]], act.func, act.mask );

    int i = act.sigs[1];  
    int j = act.sigs[2];
    signal_t rj = X[j];
    signal_t ri = X[i];
    switch( act.type )
    {
      case DecAct_t::NES_:
        if( act.id_ord == 0 )     
        { 
          rj = pNet->create_or ( X[i], X[j] );  
          ri = pNet->create_and( X[i], X[j] ); 
        }
        else if( act.id_ord == 1 ){ rj = pNet->create_and( X[i], X[j] );  ri = pNet->create_or ( X[i], X[j] ); }
        else  std::cerr << "id_ord not valid_ord" << std::endl;
        
        X[j] = rj;
        X[i] = ri;
        break; 
  
      case DecAct_t::ES_:
        if( act.id_ord == 0 )     
        { 
          rj = pNet->create_le( X[i], X[j] ); 
          ri = pNet->create_le( X[j], X[i] ); 
        }
        else if( act.id_ord == 1 )
        { 
          rj = pNet->create_lt( X[i], X[j] ); 
          ri = pNet->create_lt( X[j], X[i] ); 
        }
        else  
          std::cerr << "id_ord not valid_ord" << std::endl;
        
        X[j] = rj;
        X[i] = ri;
        break; 

      case DecAct_t::SVS_:
        if( act.id_sym == 0 ) // SVS0X_ { SVS Xj } Xi'
        {
          if( act.id_ord == 0 )
            rj = pNet->create_le( X[i], X[j] );  // ( Xi & Xj' )'
          else
            rj = pNet->create_and( X[i], X[j] ); // ( Xi & Xj )
          ri = X[i];
        }
        else if( act.id_sym == 1 ) // SVS1X  { SVS Xj } Xi
        {
          if( act.id_ord == 0 )
            rj = pNet->create_or( X[i], X[j] ); 
          else
            rj = pNet->create_lt( X[i], X[j] );
          ri = X[i];
        }
        else if( act.id_sym == 2 ) // SVSX0_  { SVS Xi } Xj'
        {
          rj = X[j];
          if( act.id_ord == 0 )
            ri = pNet->create_le( X[j], X[i] );
          else
            ri = pNet->create_and( X[i], X[j] );
        }
        else if( act.id_sym == 3 ) //SVSX1_  { SVS Xi } Xj
        {
          rj = X[j];
          if( act.id_ord == 0 )
            ri = pNet->create_or( X[i], X[j] );
          else
            ri = pNet->create_lt( X[j], X[i] );
        }
        else
          std::cerr << "wrong symmetry identifier for SVS" ;
        
        X[j] = rj;
        X[i] = ri;
        break;

      case DecAct_t::MS_:
        if( act.id_ord == 0 )
        {
          X[i] = pNet->create_xnor( X[i], X[j] ) ;
          X.erase( X.begin() + j );
          V.erase( V.begin() + j );
        }
        else if( act.id_ord == 1 )
        {
          X[j] = pNet->create_xor( X[i], X[j] );
          X.erase( X.begin() + i );
          V.erase( V.begin() + i );
        }        
        else if( act.id_ord == 2 )
        {
          X[j] = pNet->create_xnor( X[i], X[j] );
          X.erase( X.begin() + i );
          V.erase( V.begin() + i );
        }        
        else if( act.id_ord == 3 )
        {
          X[i] = pNet->create_xor( X[i], X[j] );
          X.erase( X.begin() + j );
          V.erase( V.begin() + j );
        }  
        break; 
      case DecAct_t::CSVS_:
      {
        if( act.id_sym == 0 )
        {
          if( act.id_ord == 0 )
          {
            X[i] = pNet->create_and( X[i], X[j] );
            X.erase( X.begin() + j );
            V.erase( V.begin() + j );
          }
          else
          {
            X[j] = pNet->create_and( X[i], X[j] );
            X.erase( X.begin() + i );
            V.erase( V.begin() + i );
          }
        }
        else if( act.id_sym == 1 )
        {
          if( act.id_ord == 0 )
          {
            X[i] = pNet->create_le( X[j], X[i] );
            X.erase( X.begin() + j );
            V.erase( V.begin() + j );
          }
          else
          {
            X[j] = pNet->create_lt( X[i], X[j] );
            X.erase( X.begin() + i );
            V.erase( V.begin() + i );
          }
        }
        else if( act.id_sym == 2 )
        {
          if( act.id_ord == 0 )
          {
            X[j] = pNet->create_le( X[i], X[j] );
            X.erase( X.begin() + i );
            V.erase( V.begin() + i );
          }
          else
          {
            X[i] = pNet->create_lt( X[j], X[i] );
            X.erase( X.begin() + j );
            V.erase( V.begin() + j );
          }
        }
        else if( act.id_sym == 3 )
        {
          if( act.id_ord == 0 )
          {
            X[j] =  pNet->create_or( X[i], X[j] );
            X.erase( X.begin() + i );
            V.erase( V.begin() + i );
          }
          else
          {
            X[i] =  pNet->create_or( X[i], X[j] );
            X.erase( X.begin() + j );
            V.erase( V.begin() + j );
          }
        }
        else
          std::cerr << "wrong symmetry identifier for CSVS" ;
      }
      break; 
    }
  }
#pragma endregion remap

#pragma region visualize
template<class TT, class Ntk>
void DecSolver<TT, Ntk>::PrintSpecs()
{
    printf("TRUTHS:\n");
    for( int i{0}; i<vTruths.size(); ++i ){ printf("%d ", i ); kitty::print_binary( vTruths[i] ); printf("\n");}
    printf("MASKS:\n");
    for( int i{0}; i<vMasks.size(); ++i ){ printf("%d ", i ); kitty::print_binary( vMasks[i] ); printf("\n");}
}

template<class TT, class Ntk>
void DecSolver<TT, Ntk>::show_state( DecNet<TT, Ntk> * pNet, std::vector<signal_t> * pY )
{
  for( int i{0}; i<(*pY).size(); ++i )
  {
    printf( ANSI_COLOR_YELLOW " TARGET #%d" ANSI_COLOR_RESET "", i );
    printf( ANSI_COLOR_YELLOW );
    TT F = *pNet->getFuncP( (*pY)[i] );
    TT M = *pNet->getMaskP( (*pY)[i] );
    printf("\n");
    kitty::print_binary( F );
    printf("\n");
    auto km_tt = kitty::karnaugh_map( F );
    km_tt.print( M );
    printf( ANSI_COLOR_RESET );
  }
}
#pragma endregion visualize


} // namespace mockturtle