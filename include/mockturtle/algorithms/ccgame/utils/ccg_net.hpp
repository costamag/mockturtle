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
  \file ccg_cut.hpp
  \brief data structure for storing the cuts for the ccgame

  \author Andrea Costamagna
*/
#pragma once

#include "ccg_cut.hpp"
#include "ccg_tab.hpp"
#include "ccg_analyzer.hpp"
#include <stdio.h>
#include <stack>

namespace mockturtle
{

namespace ccgame
{

using TT = kitty::dynamic_truth_table;

class net_t
{
private:
  analyzer_t analyzer;

public:
  std::vector<cut_t> cuts;
  cut_t outCut;
  int nHunging;
  uint32_t nCuts{0u};
  uint32_t nNodes{0u};
  uint32_t cost_xor{1u};

  net_t();
  net_t( std::vector<TT> const&, std::vector<TT> const& );
  net_t( std::vector<TT> const&, std::vector<uint32_t> const&, std::vector<TT> const& );
  ~net_t();

  /* manipulate */
  cut_t enumerate_divs();
  std::vector<symmetry_t> symmetry_analysis( std::vector<TT> *, int );
  std::vector<symmetry_t> symmetry_analysis( std::vector<TT> * );
  void add_cut( cut_t );
  void add_cut( symmetry_t * );
  uint32_t predelay_cost( symmetry_t * );
  void add_node_symL( cut_t *, symmetry_t * );
  void add_node_symR( cut_t *, symmetry_t * );
  void complete_cut( cut_t );
  bool check_closure( cut_t *, node_t *, node_t, node_t );
  cut_t check_closure( cut_t );
  cut_t check_closure( );
  bool check_sym_closure( );
  cut_t get_last_cut();

  template<class Ntk> Ntk convert();
  template<class Ntk> signal<Ntk> create_in_ntk( Ntk * , std::vector<signal<Ntk>> );

  /* visualize */
  void print();

};

#pragma region constructors
net_t::net_t(){};

net_t::net_t( std::vector<TT> const& X, std::vector<TT> const& Y )
{
  uint32_t INFTY = 0xFFFFFFFF;

  /* inputs */
  cut_t cut;
  cut.set_id(0x00000000);
  
  for( uint32_t i{0u}; i < X.size(); ++i )
  {
    cut.add_node( X[i], gate_t::PIS, 0, i, i );
    nNodes++;
  }
  cuts.push_back(cut);
  nCuts++;

  outCut.set_id( 0xFFFF );
  /* output */
  for( int i{0u}; i < Y.size(); ++i )
    outCut.add_node( Y[i], gate_t::POS, INFTY, INFTY );
  nHunging = Y.size();
}

net_t::net_t( std::vector<TT> const& X, std::vector<uint32_t> const& T, std::vector<TT> const& Y )
{
  uint32_t INFTY = 0xFFFFFFFF;

  /* inputs */
  cut_t cut;
  cut.set_id(0x00000000);
  
  for( uint32_t i{0u}; i < X.size(); ++i )
  {
    cut.add_node( X[i], gate_t::PIS, T[i], i, i );
    nNodes++;
  }
  cuts.push_back(cut);
  nCuts++;

  outCut.set_id( 0xFFFF );
  /* output */
  for( int i{0u}; i < Y.size(); ++i )
    outCut.add_node( Y[i], gate_t::POS, INFTY, INFTY );
  nHunging = Y.size();
}

net_t::~net_t(){}
#pragma endregion

#pragma region inquiring
/*! \brief get the last cut */
cut_t net_t::get_last_cut(){  return cuts[cuts.size()-1]; }
/*! \brief list essential candidate nodes */
cut_t net_t::enumerate_divs(){  return analyzer.enumerate_divs( get_last_cut() );}
std::vector<symmetry_t> net_t::symmetry_analysis( std::vector<TT> * pXs, int idBound )
{   
  cut_t cut = get_last_cut();
  std::vector<int> ancestor_to_node;
  std::vector<int> node_to_ancestor;
  for( int iVar{0}; iVar < cut.tt.num_vars(); ++iVar )
    ancestor_to_node.push_back(-1);
  int idNext = 0;
  bool notYet = true;
  for( int iNd{0}; iNd < cut.size(); ++iNd )
  {
    if( cut.nodes[iNd].is_remapped() )
    {
      if( cut.nodes[iNd].remapped_pi() <= idBound )
        node_to_ancestor.push_back(cut.nodes[iNd].remapped_pi());
      if( cut.nodes[iNd].remapped_pi() > idBound && notYet )
      {
        idNext = cut.nodes[iNd].remapped_pi();
        notYet = false;
      }
      ancestor_to_node[ cut.nodes[iNd].remapped_pi() ] = iNd;
    }
  }
  std::vector<symmetry_t> sym1, sym2;
  sym1 = analyzer.find_symmetries( pXs, &cut.tt, &cut.mk, &node_to_ancestor );
  if( !notYet )
  {
    std::vector<int> vNext { idBound, idNext };
    sym2 = analyzer.find_symmetries( pXs, &cut.tt, &cut.mk, &vNext );
    for( auto x : sym2 )
      sym1.push_back( x );
  }
  /* now in sym1 there should be all the symmetries */
  return sym1;
}

std::vector<symmetry_t> net_t::symmetry_analysis( std::vector<TT> * pXs )
{   
  cut_t cut = get_last_cut();
  std::vector<int> ancestor_to_node;
  std::vector<int> node_to_ancestor;
  for( int iVar{0}; iVar < cut.tt.num_vars(); ++iVar )
    ancestor_to_node.push_back(-1);
  for( int iNd{0}; iNd < cut.size(); ++iNd )
  {
    if( cut.nodes[iNd].is_remapped() )
    {
        node_to_ancestor.push_back(cut.nodes[iNd].remapped_pi());
        ancestor_to_node[ cut.nodes[iNd].remapped_pi() ] = iNd;
    }
  }
  std::vector<symmetry_t> sym;
  sym = analyzer.find_symmetries( pXs, &cut.tt, &cut.mk, &node_to_ancestor );
  /* now in sym1 there should be all the symmetries */
  return sym;
}
#pragma endregion inquiring

#pragma region manipulate

/*! \brief Add a cut to the network after adjusting its identifier. */
void net_t::add_cut( cut_t cut )
{
  cut_t newCut;
  newCut.set_id( nCuts++ );
  for( uint32_t i{0}; i < cut.size(); ++i )
    newCut.add_node( cut.nodes[i] );
  nNodes += newCut.size();
  cuts.push_back( newCut );
}

void net_t::add_node_symL( cut_t * pCut, symmetry_t * pSym )
{
  cut_t oldCut = get_last_cut();
  node_t xL, xR;
  for( int i{0}; i < oldCut.nodes.size(); ++i )
  {
    if( oldCut.nodes[i].is_remapped() )
    {
      if( oldCut.nodes[i].remapped_pi() == pSym->idL )
        xL = oldCut.nodes[i];
      else if( oldCut.nodes[i].remapped_pi() == pSym->idR )
        xR = oldCut.nodes[i];
    }
  }
  switch ( pSym->type )
  {
    case 0x33: pCut->add_node( ~( ~xL.tt &  xR.tt ), OI01, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;  //nand( l', r )
    case 0xCC: pCut->add_node(  (  xL.tt & ~xR.tt ), AI10, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;  // and( l , r')
    case 0x66: pCut->add_node( ~( ~xL.tt & ~xR.tt ), OI00, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;  //or( l , r )
    case 0x99: pCut->add_node(  (  xL.tt &  xR.tt ), AI11, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;  //and( l , r )
    case 0x44: pCut->add_node(                xL.tt, PRJL, xL.level, xL.id, xL.id ); break;  // l            
    case 0x11: pCut->add_node(                xL.tt, PRJL, xL.level, xL.id, xL.id ); break;  // l            
    case 0x77: pCut->add_node( ~( ~xL.tt & ~xR.tt ), OI00, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;  //   or( l , r )
    case 0xDD: pCut->add_node(  (  xL.tt & ~xR.tt ), AI10, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;  //  and( l , r')
    case 0x88: pCut->add_node(  (  xL.tt &  xR.tt ), AI11, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;  //  and( l , r )
    case 0x22: pCut->add_node( ~( ~xL.tt &  xR.tt ), OI01, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;  // nand( l', r )
    case 0xBB: pCut->add_node(                xL.tt, PRJL, xL.level, xL.id, xL.id ); break;  // l            
    case 0xEE: pCut->add_node(                xL.tt, PRJL, xL.level, xL.id, xL.id ); break;  // l            
    case 0x36: break;  // ]            
    case 0x6C: pCut->add_node(  (  xL.tt ^  xR.tt ), EXOR, std::max(xL.level, xR.level)+cost_xor, xL.id, xR.id ); break;  //  xor( l , r )
    case 0x9C: break;  // ]            
    case 0x39: pCut->add_node( ~(  xL.tt ^  xR.tt ), XNOR, std::max(xL.level, xR.level)+cost_xor, xL.id, xR.id ); break;  // xnor( l , r )
    case 0x19: pCut->add_node(  (  xL.tt &  xR.tt ), AI11, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;  //  and( l , r )
    case 0x26: break;  // ]            
    case 0x37: break;  // ]            
    case 0x4C: pCut->add_node(  (  xL.tt & ~xR.tt ), AI10, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;  //  and( l , r')
    case 0x8C: break;  // ]            
    case 0x3B: pCut->add_node(  ( ~xL.tt &  xR.tt ), AI01, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;  // nand( l', r )
    case 0x6E: pCut->add_node( ~( ~xL.tt & ~xR.tt ), OI00, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;  //   or( l , r )
    case 0x9D: break;  // ]            
  }
}

void net_t::add_node_symR( cut_t * pCut, symmetry_t * pSym )
{
  cut_t oldCut = get_last_cut();
  node_t xL, xR;
  for( int i{0}; i < oldCut.nodes.size(); ++i )
  {
    if( oldCut.nodes[i].is_remapped() )
    {
      if( oldCut.nodes[i].remapped_pi() == pSym->idL )
        xL = oldCut.nodes[i];
      else if( oldCut.nodes[i].remapped_pi() == pSym->idR )
        xR = oldCut.nodes[i];
    }
  }
  switch ( pSym->type )
  {
    case 0x33: pCut->add_node( ~(  xL.tt & ~xR.tt ), OI10, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;// nand( l , r')
    case 0xCC: pCut->add_node(  ( ~xL.tt &  xR.tt ), AI01, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;//  and( l', r )
    case 0x66: pCut->add_node(  (  xL.tt &  xR.tt ), AI11, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;//  and( l , r )
    case 0x99: pCut->add_node( ~( ~xL.tt & ~xR.tt ), OI00, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;//   or( l , r )
    case 0x44: pCut->add_node(  (  xL.tt &  xR.tt ), AI11, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;//  and( l , r )
    case 0x11: pCut->add_node( ~(  xL.tt & ~xR.tt ), OI10, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;// nand( l , r')
    case 0x77: pCut->add_node(                xR.tt, PRJR, xR.level, xR.id, xR.id ); break;// r            
    case 0xDD: pCut->add_node(                xR.tt, PRJR, xR.level, xR.id, xR.id ); break;// r            
    case 0x88: pCut->add_node(                xR.tt, PRJR, xR.level, xR.id, xR.id ); break;// r            
    case 0x22: pCut->add_node(                xR.tt, PRJR, xR.level, xR.id, xR.id ); break;// r            
    case 0xBB: pCut->add_node( ~( ~xL.tt & ~xR.tt ), OI00, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;//   or( l , r )
    case 0xEE: pCut->add_node(  ( ~xL.tt &  xR.tt ), AI01, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;//  and( l', r )
    case 0x36: pCut->add_node( ~(  xL.tt ^  xR.tt ), XNOR, std::max(xL.level, xR.level)+cost_xor, xL.id, xR.id ); break;// xnor( l , r )
    case 0x6C: break;// ]            
    case 0x9C: pCut->add_node(  (  xL.tt ^  xR.tt ), EXOR, std::max(xL.level, xR.level)+cost_xor, xL.id, xR.id ); break;//  xor( l , r )
    case 0x39: break;// ]            
    case 0x19: break;// ]            
    case 0x26: pCut->add_node(  (  xL.tt &  xR.tt ), AI11, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;//  and( l , r )
    case 0x37: pCut->add_node( ~(  xL.tt & ~xR.tt ), OI10, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;// nand( l , r')
    case 0x4C: break;// ]            
    case 0x8C: pCut->add_node(  ( ~xL.tt &  xR.tt ), AI01, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;//  and( l', r )
    case 0x3B: break;// ]            
    case 0x6E: break;// ]            
    case 0x9D: pCut->add_node( ~( ~xL.tt & ~xR.tt ), OI00, std::max(xL.level, xR.level)+1, xL.id, xR.id ); break;//   or( l , r )
  }
}

/*! \brief Add a cut to the network after adjusting its identifier. */
void net_t::add_cut( symmetry_t * pSym )
{
  cut_t oldCut = cuts[nCuts-1];
  cut_t newCut;
  newCut.set_id( nCuts++ );

  for( uint32_t i{0}; i < oldCut.size(); ++i )
  {
    int nNodes0 = newCut.nNodes;
    bool isRem = oldCut.nodes[i].is_remapped();
    int idPi   = oldCut.nodes[i].remapped_pi();
    if( isRem && ( idPi == pSym->idL || idPi == pSym->idR ) )
    {
      if( idPi == pSym->idL ) 
      {
        add_node_symL( &newCut, pSym );
      }
      else
        add_node_symR( &newCut, pSym );
      if( newCut.nNodes > nNodes0 ) newCut.nodes[nNodes0].idPi = idPi;
    }
    else
    {
      newCut.add_node( oldCut.nodes[i].tt, gate_t::PRJL, oldCut.nodes[i].level, oldCut.nodes[i].id, oldCut.nodes[i].id );
      if( oldCut.nodes[i].is_remapped() ) newCut.nodes[nNodes0].idPi = idPi;
    }
  }
  nNodes += newCut.size();
  newCut.tt = pSym->tt;
  newCut.mk = pSym->mk;

  cuts.push_back( newCut );
}

/*! \brief Add a cut to the network after adjusting its identifier. */
uint32_t net_t::predelay_cost( symmetry_t * pSym )
{
  cut_t oldCut = cuts[nCuts-1];
  cut_t newCut;
  uint32_t level=0;

  for( uint32_t i{0}; i < oldCut.size(); ++i )
  {
    int nNodes0 = newCut.nNodes;
    bool isRem = oldCut.nodes[i].is_remapped();
    int idPi   = oldCut.nodes[i].remapped_pi();
    if( isRem && ( idPi == pSym->idL || idPi == pSym->idR ) )
    {
      if( idPi == pSym->idL ) 
      {
        add_node_symL( &newCut, pSym );
      }
      else
        add_node_symR( &newCut, pSym );
      if( newCut.nNodes > nNodes0 ) newCut.nodes[nNodes0].idPi = idPi;
    }
    else
    {
      newCut.add_node( oldCut.nodes[i].tt, gate_t::PRJL, oldCut.nodes[i].level, oldCut.nodes[i].id, oldCut.nodes[i].id );
      if( oldCut.nodes[i].is_remapped() ) newCut.nodes[nNodes0].idPi = idPi;
    }
  }
  for( auto nd : newCut.nodes )
  {
    if( nd.level > level )
      level = nd.level;
  }

  return level;
}

/*! \brief Add nodes to the last cut. */
void net_t::complete_cut( cut_t cut )
{
  cut_t * pCut = &cuts[cuts.size()-1];
  for( int i{0}; i < cut.size(); ++i )
    pCut->add_node( cut.nodes[i] );
  nNodes += cut.size();
}

/*! \brief Check if there is a node in the input cut synthesizing an output*/
cut_t net_t::check_closure( cut_t candidates )
{
  cut_t newCut;
  newCut.set_id( nCuts );
  std::vector<int> res;
  int iDiv;
  bool found;
  for( int iOut{0}; iOut < outCut.nodes.size(); ++iOut )
  {
    node_t * pOut = &outCut.nodes[iOut];
    node_t * pDiv;
    if( pOut->gate == gate_t::POS )
    {
      found = false;
      for( int i{0}; i < newCut.nodes.size(); ++i )
      {
        if( found ) break;
        pDiv = &newCut.nodes[i];
        if( kitty::equal( pOut->tt, pDiv->tt ) )
        {
          pOut->gate = gate_t::PRJL;
          pOut->idL  = pDiv->id;
          pOut->idR  = pDiv->id;
          nHunging--;
          found = true;
        }
        else if( kitty::equal( ~pOut->tt, pDiv->tt ) )
        {
          pOut->gate = gate_t::CMPL;
          pOut->idL  = pDiv->id;
          pOut->idR  = pDiv->id;
          nHunging--;
          found = true;
        }
      }
      for( int i{0}; i < candidates.size(); ++i )
      {
        if( found ) break;
        pDiv = &candidates.nodes[i];
        if( kitty::equal( pOut->tt, pDiv->tt ) )
        {
          node_t node = newCut.add_node( *pDiv );
          pOut->gate = gate_t::PRJL;
          pOut->idL  = node.id;
          pOut->idR  = node.id;
          nHunging--;
          found = true;
        }
        else if( kitty::equal( ~pOut->tt, pDiv->tt ) )
        {
          node_t node = newCut.add_node( *pDiv );
          pOut->gate = gate_t::CMPL;
          pOut->idL  = node.id;
          pOut->idR  = node.id;
          nHunging--;
          found = true;
        }
      }
    }
  }
  return newCut;
}

bool net_t::check_closure( cut_t * pCut, node_t * pOut, node_t xL, node_t xR )
{
  if( kitty::equal( xL.tt & xR.tt, pOut->tt ) ) 
  {
    node_t node = pCut->add_node( xL.tt & xR.tt, AI11, std::max(xL.level, xR.level)+1, xL.id, xR.id );
    pOut->gate = gate_t::PRJL;
    pOut->idL  = node.id;
    pOut->idR  = node.id;
    nHunging--;
    return true;
  }
  else if( kitty::equal( xL.tt & xR.tt, ~pOut->tt ) ) 
  {
    node_t node = pCut->add_node( xL.tt & xR.tt, OI11, std::max(xL.level, xR.level)+1, xL.id, xR.id );
    pOut->gate = gate_t::PRJL;
    pOut->idL  = node.id;
    pOut->idR  = node.id;
    nHunging--;
    return true;
  }
  else if( kitty::equal( xL.tt & ~xR.tt, pOut->tt ) ) 
  {
    node_t node = pCut->add_node( xL.tt & ~xR.tt, AI10, std::max(xL.level, xR.level)+1, xL.id, xR.id );
    pOut->gate = gate_t::PRJL;
    pOut->idL  = node.id;
    pOut->idR  = node.id;
    nHunging--;
    return true;
  }
  else if( kitty::equal( xL.tt & ~xR.tt, ~pOut->tt ) ) 
  {
    node_t node = pCut->add_node( ~(xL.tt & ~xR.tt), OI10, std::max(xL.level, xR.level)+1, xL.id, xR.id );
    pOut->gate = gate_t::PRJL;
    pOut->idL  = node.id;
    pOut->idR  = node.id;
    nHunging--;
    return true;
  }
  else if( kitty::equal( ~xL.tt & xR.tt, pOut->tt ) ) 
  {
    node_t node = pCut->add_node( ~xL.tt & xR.tt, AI01, std::max(xL.level, xR.level)+1, xL.id, xR.id );
    pOut->gate = gate_t::PRJL;
    pOut->idL  = node.id;
    pOut->idR  = node.id;
    nHunging--;
    return true;
  }
  else if( kitty::equal( ~(~xL.tt & xR.tt), pOut->tt ) ) 
  {
    node_t node = pCut->add_node( ~(~xL.tt & xR.tt), OI01, std::max(xL.level, xR.level)+1, xL.id, xR.id );
    pOut->gate = gate_t::PRJL;
    pOut->idL  = node.id;
    pOut->idR  = node.id;
    nHunging--;
    return true;
  }
  else if( kitty::equal( ~xL.tt & ~xR.tt, pOut->tt ) ) 
  {
    node_t node = pCut->add_node( ~xL.tt & ~xR.tt, AI00, std::max(xL.level, xR.level)+1, xL.id, xR.id );
    pOut->gate = gate_t::PRJL;
    pOut->idL  = node.id;
    pOut->idR  = node.id;
    nHunging--;
    return true;
  }
  else if( kitty::equal( ~(~xL.tt & xR.tt), pOut->tt ) ) 
  {
    node_t node = pCut->add_node( ~(~xL.tt & ~xR.tt), OI00, std::max(xL.level, xR.level)+1, xL.id, xR.id );
    pOut->gate = gate_t::PRJL;
    pOut->idL  = node.id;
    pOut->idR  = node.id;
    nHunging--;
    return true;
  }
  else if( kitty::equal( xL.tt ^ xR.tt, pOut->tt ) ) 
  {
    node_t node = pCut->add_node( xL.tt ^ xR.tt, EXOR, std::max(xL.level, xR.level)+cost_xor, xL.id, xR.id );
    pOut->gate = gate_t::PRJL;
    pOut->idL  = node.id;
    pOut->idR  = node.id;
    nHunging--;
    return true;
  }
  else if( kitty::equal( ~( xL.tt ^ xR.tt), pOut->tt ) ) 
  {
    node_t node = pCut->add_node( ~( xL.tt ^  xR.tt), XNOR, std::max(xL.level, xR.level)+cost_xor, xL.id, xR.id );
    pOut->gate = gate_t::PRJL;
    pOut->idL  = node.id;
    pOut->idR  = node.id;
    nHunging--;
    return true;
  }
}

/*! \brief Check if there is a node in the input cut synthesizing an output*/
cut_t net_t::check_closure()
{
  cut_t newCut;
  newCut.set_id( nCuts );
  std::vector<int> res;
  int iDiv;
  bool found;
  cut_t lastCut = get_last_cut();
  for( int iOut{0}; iOut < outCut.nodes.size(); ++iOut )
  {
    node_t * pOut = &outCut.nodes[iOut];
    node_t * pDiv;
    found = false;
    if( pOut->gate == gate_t::POS )
    {    
      for( int iNdLst{0}; iNdLst<lastCut.size(); ++iNdLst )
      {
        if( found ) break;
        for( int iNdLst2{iNdLst+1}; iNdLst2<lastCut.size(); ++iNdLst2 )
        {
          if( found ) break;
          node_t xL = lastCut.nodes[iNdLst]; 
          node_t xR = lastCut.nodes[iNdLst2]; 
          found = check_closure( &newCut, pOut, xL, xR );
        }

        for( int iCut{0}; iCut < cuts.size()-1; ++iCut )
        {
          if( found ) break;
          for( int iNd{0}; iNd<cuts[iCut].size(); ++iCut )
          {
            if( found ) break;
            node_t xL = lastCut.nodes[iNd]; 
            node_t xR = lastCut.nodes[iNd]; 
            found = check_closure( &newCut, pOut, xL, xR );
          }
        }
      }
    }
  }
  return newCut;
}

/*! \brief Check if there is a node in the input cut synthesizing an output*/
bool net_t::check_sym_closure()
{
  node_t * pOut = &outCut.nodes[0];
  node_t * pDiv;
  cut_t lstCut = get_last_cut();
  for( int i{0}; i < lstCut.size(); ++i )
  {
    pDiv = &lstCut.nodes[i];
    if( kitty::equal( pOut->tt, pDiv->tt ) )
    {
      pOut->gate = gate_t::PRJL;
      pOut->idL  = pDiv->id;
      pOut->idR  = pDiv->id;
      nHunging--;
      return true;
    }
    else if( kitty::equal( ~pOut->tt, pDiv->tt ) )
    {
      pOut->gate = gate_t::CMPL;
      pOut->idL  = pDiv->id;
      pOut->idR  = pDiv->id;
      nHunging--;
      return true;
    }
  }
  return false;
}

#pragma endregion manipulate

#pragma region transform
template<class Ntk> 
Ntk net_t::convert()
{
  Ntk ntk;
  std::vector<std::vector<signal<Ntk>>> chain;
  chain.emplace_back();
  for( int iPi{0}; iPi < cuts[0].nodes.size(); ++iPi )
    chain[0].push_back( ntk.create_pi() );

  for( int iCut{1}; iCut < cuts.size(); ++iCut )
  {    
    chain.emplace_back();

    for( int iNode{0}; iNode < cuts[iCut].size(); ++iNode )
    {
      node_t * pNode = &cuts[iCut].nodes[iNode];
      assert( ( (iCut-1) == ( pNode->get_glb_idL() ) ) );
      assert( ( (iCut-1) == ( pNode->get_glb_idR() ) ) );
      signal<Ntk> xL = chain[pNode->get_glb_idL()][pNode->get_loc_idL()];
      signal<Ntk> xR = chain[pNode->get_glb_idR()][pNode->get_loc_idR()];

      switch (pNode->gate)
      {
      case gate_t::PIS :
        chain[iCut].push_back( ntk.create_pi() );
        break;
      case gate_t::AI00 :
        chain[iCut].push_back( ntk.create_and( ntk.create_not(xL), ntk.create_not(xR)) );
        break;
      case gate_t::AI01 :
        chain[iCut].push_back( ntk.create_and( ntk.create_not(xL), xR) );
        break;
      case gate_t::AI10 :
        chain[iCut].push_back( ntk.create_and( xL, ntk.create_not(xR)) );
        break;
      case gate_t::AI11 :
        chain[iCut].push_back( ntk.create_and( xL, xR ) );
        break;
      case gate_t::CMPL :
        chain[iCut].push_back( ntk.create_not(xL) );
        break;
      case gate_t::CMPR :
        chain[iCut].push_back( ntk.create_not(xR) );
        break;
      case gate_t::EXOR :
        chain[iCut].push_back( ntk.create_xor( xL, xR) );
        break;
      case gate_t::OI00 :
        chain[iCut].push_back( ntk.create_nand( ntk.create_not(xL), ntk.create_not(xR)) );
        break;
      case gate_t::OI01 :
        chain[iCut].push_back( ntk.create_nand( ntk.create_not(xL), xR) );
        break;
      case gate_t::OI10 :
        chain[iCut].push_back( ntk.create_nand( xL, ntk.create_not(xR)) );
        break;
      case gate_t::OI11 :
        chain[iCut].push_back( ntk.create_nand( xL, xR ) );
        break;
      case gate_t::PRJL :
        chain[iCut].push_back( ntk.create_buf(xL) );
        break;
      case gate_t::PRJR :
        chain[iCut].push_back( ntk.create_buf(xR) );
        break;
      case gate_t::XNOR :
        chain[iCut].push_back( ntk.create_not(ntk.create_xor( xL, xR ) ) );
        break;
      default:
        break;
      }
    }
  }

  for( int iOut{0}; iOut < outCut.nodes.size(); ++iOut )
  {
    node_t out = outCut.nodes[iOut];
    switch ( out.gate )
    {
    case gate_t::CMPL :
      ntk.create_po( ntk.create_not( chain[out.get_glb_idL()][out.get_loc_idL()] ) );
      break;
    case gate_t::CMPR :
      ntk.create_po( ntk.create_not( chain[out.get_glb_idR()][out.get_loc_idR()] ) );
      break;
    case gate_t::PRJL :
      ntk.create_po( chain[out.get_glb_idL()][out.get_loc_idL()] );
      break;
    case gate_t::PRJR :
      ntk.create_po( chain[out.get_glb_idR()][out.get_loc_idR()] );
      break;
    default:
      break;
    }
  }
  return ntk;
}


template<class Ntk> 
signal<Ntk> net_t::create_in_ntk( Ntk * pNtk, std::vector<signal<Ntk>> iSigs )
{
  std::vector<std::vector<signal<Ntk>>> chain;
  chain.emplace_back();
  assert( cuts[0].nodes.size() == iSigs.size() );
  for( int iPi{0}; iPi < cuts[0].nodes.size(); ++iPi )
    chain[0].push_back( iSigs[iPi] );

  for( int iCut{1}; iCut < cuts.size(); ++iCut )
  {    
    chain.emplace_back();

    for( int iNode{0}; iNode < cuts[iCut].size(); ++iNode )
    {
      node_t * pNode = &cuts[iCut].nodes[iNode];
      assert( ( (iCut-1) == ( pNode->get_glb_idL() ) ) );
      assert( ( (iCut-1) == ( pNode->get_glb_idR() ) ) );
      signal<Ntk> xL = chain[pNode->get_glb_idL()][pNode->get_loc_idL()];
      signal<Ntk> xR = chain[pNode->get_glb_idR()][pNode->get_loc_idR()];

      switch (pNode->gate)
      {
      case gate_t::PIS :
        chain[iCut].push_back( pNtk->create_pi() );
        break;
      case gate_t::AI00 :
        chain[iCut].push_back( pNtk->create_and( pNtk->create_not(xL), pNtk->create_not(xR)) );
        break;
      case gate_t::AI01 :
        chain[iCut].push_back( pNtk->create_and( pNtk->create_not(xL), xR) );
        break;
      case gate_t::AI10 :
        chain[iCut].push_back( pNtk->create_and( xL, pNtk->create_not(xR)) );
        break;
      case gate_t::AI11 :
        chain[iCut].push_back( pNtk->create_and( xL, xR ) );
        break;
      case gate_t::CMPL :
        chain[iCut].push_back( pNtk->create_not(xL) );
        break;
      case gate_t::CMPR :
        chain[iCut].push_back( pNtk->create_not(xR) );
        break;
      case gate_t::EXOR :
        chain[iCut].push_back( pNtk->create_xor( xL, xR) );
        break;
      case gate_t::OI00 :
        chain[iCut].push_back( pNtk->create_nand( pNtk->create_not(xL), pNtk->create_not(xR)) );
        break;
      case gate_t::OI01 :
        chain[iCut].push_back( pNtk->create_nand( pNtk->create_not(xL), xR) );
        break;
      case gate_t::OI10 :
        chain[iCut].push_back( pNtk->create_nand( xL, pNtk->create_not(xR)) );
        break;
      case gate_t::OI11 :
        chain[iCut].push_back( pNtk->create_nand( xL, xR ) );
        break;
      case gate_t::PRJL :
        chain[iCut].push_back( pNtk->create_buf(xL) );
        break;
      case gate_t::PRJR :
        chain[iCut].push_back( pNtk->create_buf(xR) );
        break;
      case gate_t::XNOR :
        chain[iCut].push_back( pNtk->create_not(pNtk->create_xor( xL, xR ) ) );
        break;
      default:
        break;
      }
    }
  }

  for( int iOut{0}; iOut < outCut.nodes.size(); ++iOut )
  {
    node_t out = outCut.nodes[iOut];
    switch ( out.gate )
    {
    case gate_t::CMPL :
      return pNtk->create_not( chain[out.get_glb_idL()][out.get_loc_idL()] );
      break;
    case gate_t::CMPR :
      return pNtk->create_not( chain[out.get_glb_idR()][out.get_loc_idR()] );
      break;
    case gate_t::PRJL :
      return chain[out.get_glb_idL()][out.get_loc_idL()];
      break;
    case gate_t::PRJR :
      return chain[out.get_glb_idR()][out.get_loc_idR()];
      break;
    default:
      break;
    }
  }
}

#pragma endregion transform

#pragma region visualize
void net_t::print()
{
  for( int i{0}; i < cuts.size(); ++i )
  {
    printf(" CUT %d\n", i );
    cuts[i].print();
  }
  printf("\n" );
  printf("OUTPUTS:\n" );
  outCut.print();

}
#pragma endregion visualize

} // namespace ccgame

} // namespace mockturtle