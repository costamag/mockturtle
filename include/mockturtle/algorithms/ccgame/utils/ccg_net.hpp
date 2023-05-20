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

using DTT = kitty::dynamic_truth_table;

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

  net_t( std::vector<TT> const&, std::vector<TT> const& );
  ~net_t();

  /* manipulate */
  cut_t enumerate_divs();
  std::vector<symmetry_t> symmetry_analysis( std::vector<DTT>, int );
  void add_cut( cut_t );
  void complete_cut( cut_t );
  cut_t check_closure( cut_t );
  cut_t get_last_cut();

  template<class Ntk> Ntk convert();

  /* visualize */
  void print();

};

#pragma region constructors
net_t::net_t( std::vector<TT> const& X, std::vector<TT> const& Y )
{
  uint32_t INFTY = 0xFFFFFFFF;

  /* inputs */
  cut_t cut;
  cut.set_id(0u);
  
  for( uint32_t i{0u}; i < X.size(); ++i )
  {
    cut.add_node( X[i], gate_t::PIS, INFTY, INFTY );
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
std::vector<symmetry_t> net_t::symmetry_analysis( std::vector<DTT> xs, int idBound )
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
  sym1 = analyzer.find_symmetries( xs, cut.tt, cut.mk, node_to_ancestor );
  if( !notYet )
  {
    sym2 = analyzer.find_symmetries( xs, cut.tt, cut.mk, { idBound, idNext } );
    for( auto x : sym2 )
      sym1.push_back( x );
  }
  /* now in sym1 there should be all the symmetries */
  analyzer.print_symmetries( sym1 );
  return sym1;
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