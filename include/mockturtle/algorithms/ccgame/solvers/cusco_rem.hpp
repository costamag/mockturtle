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

using TT = kitty::dynamic_truth_table;


template<class Ntk>
struct report_rem_t
{
  int nIt0;
  int nMin;
  int nMax;
  Ntk ntk;
};

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
  report_rem_t<Ntk> solve_random( cusco_rem_ps const& );
  report_rem_t<Ntk> solve_entropic( cusco_rem_ps const& );
  report_rem_t<Ntk> solve_1shot( cusco_rem_ps const& );
};

/* creation and destruction */
template<class Ntk>
cusco_rem<Ntk>::cusco_rem( std::vector<TT> const& X, std::vector<TT> const& Y ): X(X), Y(Y) {}
template<class Ntk>
cusco_rem<Ntk>::~cusco_rem(){}

template<class Ntk>
int sym_cost( symmetry_t sym )
{
  int costAND = 1u;
  int costINV = 0u;
  int costXOR = std::is_same<Ntk, xag_network>::value ? 1 : 3;
  switch ( sym.type )
  {
      case 0x33: return 2*costAND + 4*costINV; break; // nand( l', r )   nand( l , r')
      case 0xCC: return 2*costAND + 2*costINV; break; //  and( l , r')    and( l', r )
      case 0x66: return 2*costAND + 3*costINV; break; //   or( l , r )    and( l , r )
      case 0x99: return 2*costAND + 3*costINV; break; //  and( l , r )     or( l , r )
      case 0x44: return costAND;               break; // l                and( l , r )
      case 0x11: return costAND   + 2*costINV; break; // l               nand( l , r')
      case 0x77: return costAND   + 3*costINV; break; //   or( l , r )   r            
      case 0xDD: return costAND   +   costINV; break; //  and( l , r')   r            
      case 0x88: return costAND              ; break; //  and( l , r )   r            
      case 0x22: return costAND   + 2*costINV; break; // nand( l', r )   r            
      case 0xBB: return costAND   + 3*costINV; break; // l                 or( l , r )
      case 0xEE: return costAND   +   costINV; break; // l                and( l', r )
      case 0x36: return costXOR   +   costINV; break; // ]               xnor( l , r )
      case 0x6C: return costXOR              ; break; //  xor( l , r )   ]            
      case 0x9C: return costXOR              ; break; // ]                xor( l , r )
      case 0x39: return costXOR   +   costINV; break; // xnor( l , r )   ]            
      case 0x19: return costAND              ; break; //  and( l , r )   ]            
      case 0x26: return costAND              ; break; // ]                and( l , r )
      case 0x37: return costAND   + 2*costINV; break; // ]               nand( l , r')
      case 0x4C: return costAND   +   costINV; break; //  and( l , r')   ]            
      case 0x8C: return costAND   +   costINV; break; // ]                and( l', r )
      case 0x3B: return costAND   + 2*costINV; break; // nand( l', r )   ]            
      case 0x6E: return costAND   + 3*costINV; break; //   or( l , r )   ]            
      case 0x9D: return costAND   + 3*costINV; break; //   or( l , r ) 
    break;
  }
}

template<class Ntk>
report_rem_t<Ntk> cusco_rem<Ntk>::solve_1shot( cusco_rem_ps const& ps )
{
  report_rem_t<Ntk> rep;
  std::random_device rd;  // a seed source for the random number engine
  std::mt19937 gen(5); // mersenne_twister_engine seeded with rd()
  int nVars = ceil(log2(X[0].num_bits()));
  TT func( nVars );
  TT mask( nVars );
  kitty::create_from_binary_string( func, kitty::to_binary( Y[0] ) );
  kitty::create_from_binary_string( mask, kitty::to_binary( Y[0] | ~Y[0] ) );
  std::vector<TT> xs;
  for( int i{0}; i < nVars; ++i )
  {
    xs.emplace_back( nVars );
    kitty::create_nth_var( xs[i], i );
  }
  analyzer_t analyzer;

  net_t net( X, Y );
  net.cuts[net.nCuts-1].set_func( func );
  net.cuts[net.nCuts-1].set_mask( mask );
  int idBound = 1;
  int bestRwd = -1;
  symmetry_t bestSym;
  while( net.nHunging > 0 )
  {
    //std::vector<symmetry_t> candidates = net.symmetry_analysis( &xs );
    std::vector<symmetry_t> candidates = net.symmetry_analysis( &xs, idBound );
    int bestCost = 1000;
    if( candidates.size() == 0 ) break;
    for( int iCand{0}; iCand<candidates.size(); ++iCand )
    {
      int cost = sym_cost<Ntk>( candidates[iCand] );
      if( candidates[iCand].rwd > bestRwd || ( candidates[iCand].rwd == bestRwd && cost < bestCost ) )
      {
        bestSym = candidates[iCand];
        bestRwd = candidates[iCand].rwd;
        bestCost = cost;
      }
    }
    if( ( bestSym.idL == idBound ) || ( bestSym.idR == idBound  ) )
      idBound += 2;
    net.add_cut( &bestSym );
    net.check_sym_closure();
  }
  if( net.nHunging > 0 )
  {
    rep.nMax = -1;
    rep.nMin = -1;
  }
  else
  {
    rep.ntk  = net.convert<Ntk>();
    rep.nIt0 = rep.ntk.num_gates();
    rep.nMax = rep.nIt0;
    rep.nMin = rep.nIt0;
  }
  return rep;
}


template<class Ntk>
report_rem_t<Ntk> cusco_rem<Ntk>::solve_random( cusco_rem_ps const& ps )
{
  report_rem_t<Ntk> rep;
  rep.nMax = 0;
  rep.nMin = 0x00FFFFFF;

  std::random_device rd;  // a seed source for the random number engine
  std::mt19937 gen(5); // mersenne_twister_engine seeded with rd()
  Ntk ntk;
  int nVars = ceil(log2(X[0].num_bits()));
  TT func( nVars );
  TT mask( nVars );
  kitty::create_from_binary_string( func, kitty::to_binary( Y[0] ) );
  kitty::create_from_binary_string( mask, kitty::to_binary( Y[0] | ~Y[0] ) );
  int nBest = 10000u;
  std::vector<TT> xs;
  for( int i{0}; i < nVars; ++i )
  {
    xs.emplace_back( nVars );
    kitty::create_nth_var( xs[i], i );
  }
      analyzer_t analyzer;

  for( int iIt{0}; iIt < ps.nIters; ++iIt )
  {
    bool notFound = false;
    net_t net( X, Y );
    net.cuts[net.nCuts-1].set_func( func );
    net.cuts[net.nCuts-1].set_mask( mask );
    int idBound = 1;
    int bestRwd = -1;
    while( net.nHunging > 0 )
    {
      std::vector<symmetry_t> candidates = net.symmetry_analysis( &xs, idBound );
      int bestCost = 1000;
      std::vector<symmetry_t> bestCands;
      if( candidates.size() == 0 )
      {
        notFound = true;
        break;
      }
      for( int i{0}; i<candidates.size(); ++i )
      {
        int cost = sym_cost<Ntk>( candidates[i] );
        if( candidates[i].rwd > bestRwd || ( candidates[i].rwd == bestRwd && cost < bestCost ) )
          {
            bestCands = { candidates[i] };
            bestRwd = candidates[i].rwd;
            bestCost = cost;
          }
          else if( iIt > 0 && ( candidates[i].rwd == bestRwd ) && ( cost == bestCost ) )
          {
            bestCands.push_back( candidates[i] );
          }
      }
      std::uniform_int_distribution<> distrib(0, bestCands.size()-1);
      int iSym = distrib(gen);

      symmetry_t bestSym = bestCands[iSym];
      net.add_cut( &bestSym );
      cut_t lastCut = net.get_last_cut();

      if( ( bestSym.idL >= idBound ) || ( bestSym.idR >= idBound  ) )
      {
        idBound += 2;
      }
      net.check_sym_closure();
    }
    if( !notFound )
    {
      ntk = net.convert<Ntk>();
      if( iIt == 0 )  rep.nIt0 = ntk.num_gates();;
      if( ntk.num_gates() < rep.nMin )
      {
        rep.ntk = ntk;
        rep.nMin = ntk.num_gates();
      }
      if ( ntk.num_gates() > rep.nMax )  rep.nMax = ntk.num_gates();
    }
  }
  return rep;
}

template<class Ntk>
report_rem_t<Ntk> cusco_rem<Ntk>::solve_entropic( cusco_rem_ps const& ps )
{
  report_rem_t<Ntk> rep;
  rep.nMax = 0;
  rep.nMin = 0x00FFFFFF;

  std::random_device rd;  // a seed source for the random number engine
  std::mt19937 gen(5); // mersenne_twister_engine seeded with rd()
  Ntk ntk;
  int nVars = ceil(log2(X[0].num_bits()));
  TT func( nVars );
  TT mask( nVars );
  kitty::create_from_binary_string( func, kitty::to_binary( Y[0] ) );
  kitty::create_from_binary_string( mask, kitty::to_binary( Y[0] | ~Y[0] ) );
  int nBest = 10000u;
  std::vector<TT> xs;
  for( int i{0}; i < nVars; ++i )
  {
    xs.emplace_back( nVars );
    kitty::create_nth_var( xs[i], i );
  }
      analyzer_t analyzer;

  for( int iIt{0}; iIt < ps.nIters; ++iIt )
  {
    bool notFound = false;
    net_t net( X, Y );
    net.cuts[net.nCuts-1].set_func( func );
    net.cuts[net.nCuts-1].set_mask( mask );
    int oldRwd = 0;
    while( net.nHunging > 0 )
    {      
      std::vector<symmetry_t> candidates = net.symmetry_analysis( &xs );
      int bestCost = 1000;
      std::vector<symmetry_t> bestCands;
      if( candidates.size() == 0 )
      {
        notFound = true;
        break;
      }
      for( int i{0}; i<candidates.size(); ++i )
      {
        if( candidates[i].rwd >= oldRwd )
          bestCands.push_back( candidates[i] );
      }
      std::uniform_int_distribution<> distrib(0, bestCands.size()-1);
      int iSym = distrib(gen);
      symmetry_t bestSym = bestCands[iSym];
      oldRwd = bestSym.rwd;
      net.add_cut( &bestSym );
      cut_t lastCut = net.get_last_cut();

      net.check_sym_closure();
    }
    if( !notFound )
    {
      ntk = net.convert<Ntk>();
      if( iIt == 0 )  rep.nIt0 = ntk.num_gates();;
      if( ntk.num_gates() < rep.nMin )
      {
        rep.ntk = ntk;
        rep.nMin = ntk.num_gates();
      }
      if ( ntk.num_gates() > rep.nMax )  rep.nMax = ntk.num_gates();
    }
  }
  return rep;
}

} // namespace ccgame

} // namespace mockturtle
