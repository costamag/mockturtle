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
#include "../utils/ccg_rng.hpp"
#include "../utils/mct_utils.hpp"
#include "../utils/ccg_analyzer.hpp"
#include "../utils/ccg_net.hpp"
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
  double area;
  double levels;
  bool E_solution{false};
  std::vector<signal<Ntk>> S;
  Ntk * pNtk;
  signal<Ntk> osig;
};

struct cusco_rem_ps
{
  /* method */
  int nIters;
  std::vector<double> T;
  library__t lib;


  cusco_rem_ps( int nIters ) : nIters( nIters ) {}
  cusco_rem_ps( int nIters, library__t Lib ) : nIters( nIters ), lib(Lib) {}
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
  report_rem_t<Ntk> solve_1delay( cusco_rem_ps const& );
  report_rem_t<Ntk> solve_1delay( cusco_rem_ps const&, Ntk *, std::vector<signal<Ntk>> );
  report_rem_t<Ntk> solve_Rdelay( cusco_rem_ps const& );
  report_rem_t<Ntk> solve_Rdelay( cusco_rem_ps const&, Ntk *, std::vector<signal<Ntk>> );
};

/* creation and destruction */
template<class Ntk>
cusco_rem<Ntk>::cusco_rem( std::vector<TT> const& X, std::vector<TT> const& Y ): X(X), Y(Y) {}
template<class Ntk>
cusco_rem<Ntk>::~cusco_rem(){}

template<class Ntk>
int sym_cost( symmetry_t sym )
{
  double costAND = 1u;
  double costINV = 0u;
  double costXOR = std::is_same<Ntk, xag_network>::value ? 1 : 3;
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
report_rem_t<Ntk> cusco_rem<Ntk>::solve_1delay( cusco_rem_ps const& ps, Ntk * pNtk, std::vector<signal<Ntk>> S )
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
  std::vector<double> T;
  for( int i{0}; i < nVars; ++i )
  {
    xs.emplace_back( nVars );
    kitty::create_nth_var( xs[i], i );
    T.push_back( ps.T[i] );
  }
  analyzer_t analyzer;

  net_t net( X, T, Y, ps.lib );

  net.cuts[net.nCuts-1].set_func( func );
  net.cuts[net.nCuts-1].set_mask( mask );
  int idBound = 1;
  int bestRwd = -1;
  double bestLevel=0;
  symmetry_t bestSym;
  while( net.nHunging > 0 )
  {
    std::vector<symmetry_t> candidates = net.symmetry_analysis( &xs, idBound );
    // process candidate to give levels
    double bestCost = 1000;
    int iCHOSEN;
    if( candidates.size() == 0 ) break;
    for( int iCand{0}; iCand<candidates.size(); ++iCand )
    {
      int size_ = sym_cost<Ntk>( candidates[iCand] );
      double cost = net.predelay_cost(&candidates[iCand]);
      if( candidates[iCand].rwd > bestRwd || ( candidates[iCand].rwd == bestRwd && cost < bestCost ) )
      {
        bestSym = candidates[iCand];
        bestRwd = candidates[iCand].rwd;
        bestCost = cost;
        bestLevel = cost;
        iCHOSEN=iCand;
        rep.levels = cost;
      }
    }
    if( ( bestSym.idL == idBound ) || ( bestSym.idR == idBound  ) )
      idBound += 2;
    net.add_cut( &bestSym );

    cut_t lcut = net.get_last_cut();

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
    if( S.size() > 0 )
      rep.osig = net.create_in_ntk<Ntk>( pNtk, S );
    rep.nIt0 = rep.ntk.num_gates();
    rep.nMax = rep.nIt0;
    rep.nMin = rep.nIt0;
    rep.E_solution = true;
    node_t nd_out = net.outCut.nodes[0];
    rep.levels = nd_out.level;
  }
  return rep;
}

template<class Ntk>
report_rem_t<Ntk> cusco_rem<Ntk>::solve_1delay( cusco_rem_ps const& ps )
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
  std::vector<double> T;
  for( int i{0}; i < nVars; ++i )
  {
    xs.emplace_back( nVars );
    kitty::create_nth_var( xs[i], i );
    T.push_back( ps.T[i] );
  }
  analyzer_t analyzer;

  net_t net( X, T, Y, ps.lib );

  auto rt = net.cuts[net.nCuts-1].nodes;
  //for( auto nd : rt ) 
  //  printf("%f \n", nd.level );
  //printf("\n");

  net.cuts[net.nCuts-1].set_func( func );
  net.cuts[net.nCuts-1].set_mask( mask );
  int idBound = 1;
  int bestRwd = -1;
  double bestLevel=0;
  symmetry_t bestSym;
  while( net.nHunging > 0 )
  {
    std::vector<symmetry_t> candidates = net.symmetry_analysis( &xs, idBound );
    // process candidate to give levels
    double bestCost = 1000;
    int iCHOSEN;
    if( candidates.size() == 0 ) break;
    bool NONEFOUND{true};
    for( int iCand{0}; iCand<candidates.size(); ++iCand )
    {
      int size_ = sym_cost<Ntk>( candidates[iCand] );
      double cost = net.predelay_cost(&candidates[iCand]);
      //printf("%f ", cost);
      if( candidates[iCand].rwd > bestRwd || ( candidates[iCand].rwd == bestRwd && cost < bestCost ) )
      {
        bestSym = candidates[iCand];
        bestRwd = candidates[iCand].rwd;
        bestCost = cost;
        bestLevel = cost;
        iCHOSEN=iCand;
        rep.levels = cost;
        NONEFOUND = false;
      }
    }
    //printf(":: cost = %f\n", bestCost);
    if( NONEFOUND ) break;
    if( ( bestSym.idL == idBound ) || ( bestSym.idR == idBound  ) )
      idBound += 2;
    net.add_cut( &bestSym );

    cut_t lcut = net.get_last_cut();

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
    rep.nMin = rep.ntk.num_gates();
    rep.area = net.compute_area<Ntk>();
    rep.nIt0 = rep.nMin;
    rep.nMax = rep.nMin;
    rep.E_solution = true;
    node_t nd_out = net.outCut.nodes[0];
    rep.levels = nd_out.level;
    printf("d=%f a=%f\n", rep.levels, rep.area);
  }
  //net.print();
  return rep;
}


template<class Ntk>
report_rem_t<Ntk> cusco_rem<Ntk>::solve_Rdelay( cusco_rem_ps const& ps )
{
  std::random_device rd;  // a seed source for the random number engine
  std::mt19937 gen(5); // mersenne_twister_engine seeded with rd()
  int nVars = ceil(log2(X[0].num_bits()));
  TT func( nVars );
  TT mask( nVars );
  kitty::create_from_binary_string( func, kitty::to_binary( Y[0] ) );
  kitty::create_from_binary_string( mask, kitty::to_binary( Y[0] | ~Y[0] ) );
  std::vector<TT> xs;
  std::vector<double> T;
  for( int i{0}; i < nVars; ++i )
  {
    xs.emplace_back( nVars );
    kitty::create_nth_var( xs[i], i );
    T.push_back( ps.T[i] );
  }
  analyzer_t analyzer;

  report_rem_t<Ntk> REP;
  double BEST_SIZE = std::numeric_limits<double>::max();
  double BEST_DEPTH = std::numeric_limits<double>::max();
  for( int iIter{0}; iIter < ps.nIters; ++iIter )
  {
    report_rem_t<Ntk> rep;
    net_t net( X, T, Y, ps.lib );

    net.cuts[net.nCuts-1].set_func( func );
    net.cuts[net.nCuts-1].set_mask( mask );
    int idBound = 1;
    int bestRwd = -1;
    double bestLevel=0;
    symmetry_t bestSym;
    int old_reward = 0;
    while( net.nHunging > 0 )
    {
      std::vector<symmetry_t> candidates = net.symmetry_analysis( &xs, idBound );
      // process candidate to give levels
      double bestCost = 1000;
      if( candidates.size() == 0 ) break;
      std::vector<int> selected;
      for( int iCand{0}; iCand<candidates.size(); ++iCand )
      {
        int size_ = sym_cost<Ntk>( candidates[iCand] );
        double cost = net.predelay_cost(&candidates[iCand]);
        if( candidates[iCand].rwd > bestRwd || ( candidates[iCand].rwd == bestRwd && cost <= bestCost ) )
        {
          if( (candidates[iCand].rwd > bestRwd) || ( candidates[iCand].rwd == bestRwd && cost < bestCost ) )
          {
            // reset
            selected={iCand};
            bestRwd = candidates[iCand].rwd;
            bestCost = cost;
            bestLevel = cost;
          }
          else if( candidates[iCand].rwd == bestRwd && cost == bestCost )
          {
            selected.push_back( iCand );
            bestRwd = candidates[iCand].rwd;
            bestCost = cost;
            bestLevel = cost;
          }
          else
            assert(0);
        }
      }

      std::uniform_int_distribution<> distrib(0, selected.size()-1);
      uint32_t IDX = distrib(ccg_gen);
      uint32_t iCHOSEN = selected[IDX];
      rep.levels = net.predelay_cost(&candidates[iCHOSEN]);;
      bestSym = candidates[iCHOSEN];
      bestRwd = candidates[iCHOSEN].rwd;
      bestCost = rep.levels;
      bestLevel = rep.levels;

      if( ( bestSym.idL == idBound ) || ( bestSym.idR == idBound  ) )
        idBound += 2;
      net.add_cut( &bestSym );

      cut_t lcut = net.get_last_cut();

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
      rep.ntk = cleanup_dangling( rep.ntk );
      rep.area = net.compute_area<Ntk>();
      rep.nIt0 = rep.ntk.num_gates();
      rep.nMax = rep.nIt0;
      rep.nMin = rep.nIt0;
      rep.E_solution = true;
      node_t nd_out = net.outCut.nodes[0];
      rep.levels = nd_out.level;
      if( rep.levels < BEST_DEPTH || ( rep.levels == BEST_DEPTH && rep.area < BEST_SIZE  ))
      {
        BEST_DEPTH     = rep.levels;
        BEST_SIZE      = rep.area;
        REP.area       = rep.area;
        REP.ntk        = rep.ntk         ;
        REP.nIt0       = rep.nIt0        ;
        REP.nMax       = rep.nMax        ;
        REP.nMin       = rep.nMin        ;
        REP.E_solution = true ;
        REP.levels     = rep.levels     ;
      }
    }
  }
  if( REP.E_solution )
    printf("d=%f a=%f\n", REP.levels, REP.area);

  return REP;
}

template<class Ntk>
report_rem_t<Ntk> cusco_rem<Ntk>::solve_Rdelay( cusco_rem_ps const& ps, Ntk * pNtk, std::vector<signal<Ntk>> S )
{
  std::random_device rd;  // a seed source for the random number engine
  std::mt19937 gen(5); // mersenne_twister_engine seeded with rd()
  int nVars = ceil(log2(X[0].num_bits()));
  TT func( nVars );
  TT mask( nVars );
  kitty::create_from_binary_string( func, kitty::to_binary( Y[0] ) );
  kitty::create_from_binary_string( mask, kitty::to_binary( Y[0] | ~Y[0] ) );
  std::vector<TT> xs;
  std::vector<double> T;
  for( int i{0}; i < nVars; ++i )
  {
    xs.emplace_back( nVars );
    kitty::create_nth_var( xs[i], i );
    T.push_back( ps.T[i] );
  }
  analyzer_t analyzer;

  report_rem_t<Ntk> REP;
  net_t NET( X, T, Y, ps.lib );
  bool FOUND_ONE = false;
  uint32_t BEST_SIZE = std::numeric_limits<uint32_t>::max();
  uint32_t BEST_DEPTH = std::numeric_limits<uint32_t>::max();
  for( int iIter{0}; iIter < ps.nIters; ++iIter )
  {
    report_rem_t<Ntk> rep;
    net_t net( X, T, Y, ps.lib );

    net.cuts[net.nCuts-1].set_func( func );
    net.cuts[net.nCuts-1].set_mask( mask );
    int idBound = 1;
    int bestRwd = -1;
    double bestLevel=0;
    symmetry_t bestSym;
    int old_reward = 0;

    while( net.nHunging > 0 )
    {
      std::vector<symmetry_t> candidates = net.symmetry_analysis( &xs, idBound );
      // process candidate to give levels
      double bestCost = 1000;
      if( candidates.size() == 0 ) break;
      std::vector<int> selected;
      for( int iCand{0}; iCand<candidates.size(); ++iCand )
      {
        int size_ = sym_cost<Ntk>( candidates[iCand] );
        double cost = net.predelay_cost(&candidates[iCand]);

        if( candidates[iCand].rwd > bestRwd || ( candidates[iCand].rwd == bestRwd && cost <= bestCost ) )
        {
          if( (candidates[iCand].rwd > bestRwd) || ( candidates[iCand].rwd == bestRwd && cost < bestCost ) )
          {
            // reset
            selected={iCand};
            bestRwd = candidates[iCand].rwd;
            bestCost = cost;
            bestLevel = cost;
          }
          else if( candidates[iCand].rwd == bestRwd && cost == bestCost )
          {
            selected.push_back( iCand );
            bestRwd = candidates[iCand].rwd;
            bestCost = cost;
            bestLevel = cost;
          }
          else
            assert(0);
        }
      }
      if( selected.size() <= 0 )
        break;
      std::uniform_int_distribution<> distrib(0, selected.size()-1);
      uint32_t IDX = distrib(ccg_gen);
      uint32_t iCHOSEN = selected[IDX];
      rep.levels = net.predelay_cost(&candidates[iCHOSEN]);
      bestSym = candidates[iCHOSEN];
      bestRwd = candidates[iCHOSEN].rwd;
      bestCost = rep.levels;
      bestLevel = rep.levels;

      if( ( bestSym.idL == idBound ) || ( bestSym.idR == idBound  ) )
        idBound += 2;
      net.add_cut( &bestSym );
      cut_t lcut = net.get_last_cut();
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
      rep.ntk = cleanup_dangling( rep.ntk );
      rep.nIt0 = rep.ntk.num_gates();
      rep.nMax = rep.nIt0;
      rep.nMin = rep.nIt0;
      rep.E_solution = true;
      rep.levels = bestLevel;
      if( rep.levels < BEST_DEPTH || ( rep.levels == BEST_DEPTH && rep.nMin < BEST_SIZE  ))
      {
        FOUND_ONE = true;
        BEST_DEPTH     = rep.levels;
        BEST_SIZE      = rep.nMin;
        REP.ntk        = rep.ntk ;
        REP.nIt0       = rep.nIt0;
        REP.nMax       = rep.nMax;
        REP.nMin       = rep.nMin;
        REP.E_solution = true ;
        REP.levels     = rep.levels;
        NET = net;
      }
    }
  }
  if( S.size() > 0 )
    REP.osig = NET.create_in_ntk<Ntk>( pNtk, S );

  return REP;
}

} // namespace ccgame

} // namespace mockturtle
