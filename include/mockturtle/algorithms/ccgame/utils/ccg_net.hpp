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
#include <stdio.h>
#include <stack>

namespace mockturtle
{

namespace ccgame
{

class net_t
{

public:
  std::vector<cut_t> cuts;
  cut_t output_c;
  int nNodes{0};
  int nHunging;

  net_t( std::vector<TT> const&, std::vector<TT> const& );
  ~net_t();

  /* manipulate */
  cut_t list_candidate_divs();
  void add_cut( cut_t );
  void complete_cut( cut_t );
  bool check_closure( cut_t );

  template<class Ntk> Ntk convert();

  /* visualize */
  void print();

};

#pragma region constructors
net_t::net_t( std::vector<TT> const& X, std::vector<TT> const& Y )
{
  /* inputs */
  cut_t cut;
  for( int i{0u}; i < X.size(); ++i )
  {
    cut.add_divisor( X[i], gate_t::PIS, i, i, i, 0u );
    nNodes++;
  }
  cuts.push_back(cut);
  /* output */
  for( int i{0u}; i < Y.size(); ++i )
    output_c.add_divisor( Y[i], gate_t::POS, 0xFFFFFFFF, 0xFFFFFFFF, i, 0u );
  nHunging = Y.size();
}

net_t::~net_t(){}
#pragma endregion

#pragma region manipulate
/*! \brief combines the gates in the last cut to propose all possible divisors*/
cut_t net_t::list_candidate_divs()
{
    cut_t divs;
    cut_t cut = cuts[cuts.size()-1]; 

    for( int i{0}; i < cut.divisors.size() ; ++i )
    {
        divisor_t d = cut.divisors[i];
        divs.add_divisor( d.tt, PRJL, d.inL, d.inL, 0xFFFFFFFF, 0x0 );
    }

    divisor_t dL;
    divisor_t dR;

    for( int iR{0}; iR < cut.divisors.size()-1 ; ++iR )
    {
        for( int iL{iR+1}; iL < cut.divisors.size(); ++iL ) // iL > iR
        {
            dL = cut.divisors[iL];
            dR = cut.divisors[iR];
            divs.add_divisor(  dL.tt &  dR.tt, AI11, dL.id, dR.id, 0xFFFFFFFF, 0x0 );
            divs.add_divisor(  dL.tt & ~dR.tt, AI10, dL.id, dR.id, 0xFFFFFFFF, 0x0 );
            divs.add_divisor( ~dL.tt &  dR.tt, AI01, dL.id, dR.id, 0xFFFFFFFF, 0x0 );
            divs.add_divisor( ~dL.tt & ~dR.tt, AI00, dL.id, dR.id, 0xFFFFFFFF, 0x0 );
            divs.add_divisor(  dL.tt ^  dR.tt, EXOR, dL.id, dR.id, 0xFFFFFFFF, 0x0 );
        }
    }
    return divs;
} 

void net_t::add_cut( cut_t cut )
{
  cuts.push_back( cut );
  for( int i{0}; i<cuts[ cuts.size() - 1 ].divisors.size(); ++i )
    cuts[ cuts.size() - 1 ].divisors[i].id = nNodes++;
}

void net_t::complete_cut( cut_t cut )
{
  cut_t * pCut = &cuts[cuts.size()-1];
  for( int i{0}; cut.size(); ++i )
  {
    pCut->divisors.push_back( cut.divisors[i] );
    pCut->divisors[pCut->size()-1].id = nNodes++;
  }
}

bool net_t::check_closure( cut_t candidates )
{
  cut_t new_c;
  std::vector<int> res;
  int iDiv;
  bool found;
  for( int iOut{0}; iOut < output_c.divisors.size(); ++iOut )
  {
    divisor_t * pOut = &output_c.divisors[iOut];
    divisor_t * pDiv;
    if( pOut->gate == gate_t::POS )
    {
      found = false;
      for( int i{0}; i < new_c.divisors.size(); ++i )
      {
        if( found ) break;
        pDiv = &new_c.divisors[i];
        if( kitty::equal( pOut->tt, pDiv->tt ) )
        {
          pOut->gate = gate_t::PRJL;
          pOut->inL  = pDiv->id;
          pOut->inR  = pDiv->id;
          nHunging--;
          found = true;
        }
        else if( kitty::equal( ~pOut->tt, pDiv->tt ) )
        {
          pOut->gate = gate_t::CMPL;
          pOut->inL  = pDiv->id;
          pOut->inR  = pDiv->id;
          nHunging--;
          found = true;
        }
      }
      for( int i{0}; i < candidates.size(); ++i )
      {
        if( found ) break;
        pDiv = &candidates.divisors[i];
        if( kitty::equal( pOut->tt, pDiv->tt ) )
        {
          new_c.add_divisor( *pDiv );
          pOut->gate = gate_t::PRJL;
          pOut->inL  = nNodes;
          pOut->inR  = nNodes;
          new_c.divisors[new_c.divisors.size()-1].id = nNodes++;
          nHunging--;
          found = true;
          res.push_back( i );
        }
        else if( kitty::equal( ~pOut->tt, pDiv->tt ) )
        {
          new_c.add_divisor( *pDiv );
          pOut->gate = gate_t::CMPL;
          pOut->inL  = nNodes;
          pOut->inR  = nNodes;
          new_c.divisors[new_c.divisors.size()-1].id = nNodes++;
          nHunging--;
          found = true;
          res.push_back( i );
        }
      }
    }
  }
  if( res.size() > 0 )
    cuts.push_back( new_c );
  return( res.size() > 0 );
}
#pragma endregion manipulate

#pragma region transform
template<class Ntk> 
Ntk net_t::convert()
{
  Ntk ntk;
  std::vector<signal<Ntk>> chain;
  for( int iPi{0}; iPi < cuts[0].divisors.size(); ++iPi )
    chain.push_back( ntk.create_pi() );
  for( int iCut{0}; iCut < cuts.size(); ++iCut )
  {    
    for( int iNode{0}; iNode < cuts[iCut].divisors.size(); ++iNode )
    {
      divisor_t * pNode = &cuts[iCut].divisors[iNode];
      signal<Ntk> xL = chain[pNode->inL];
      signal<Ntk> xR = chain[pNode->inR];

      switch (pNode->gate)
      {
      case gate_t::AI00 :
        chain.push_back( ntk.create_and( ntk.create_not(xL), ntk.create_not(xR)) );
        break;
      case gate_t::AI01 :
        chain.push_back( ntk.create_and( ntk.create_not(xL), xR) );
        break;
      case gate_t::AI10 :
        chain.push_back( ntk.create_and( xL, ntk.create_not(xR)) );
        break;
      case gate_t::AI11 :
        chain.push_back( ntk.create_and( xL, xR ) );
        break;
      case gate_t::CMPL :
        chain.push_back( ntk.create_not(xL) );
        break;
      case gate_t::CMPR :
        chain.push_back( ntk.create_not(xR) );
        break;
      case gate_t::EXOR :
        chain.push_back( ntk.create_xor( xL, xR) );
        break;
      case gate_t::OI00 :
        chain.push_back( ntk.create_nand( ntk.create_not(xL), ntk.create_not(xR)) );
        break;
      case gate_t::OI01 :
        chain.push_back( ntk.create_nand( ntk.create_not(xL), xR) );
        break;
      case gate_t::OI10 :
        chain.push_back( ntk.create_nand( xL, ntk.create_not(xR)) );
        break;
      case gate_t::OI11 :
        chain.push_back( ntk.create_nand( xL, xR ) );
        break;
      case gate_t::PRJL :
        chain.push_back( ntk.create_buf(xL) );
        break;
      case gate_t::PRJR :
        chain.push_back( ntk.create_buf(xR) );
        break;
      case gate_t::XNOR :
        chain.push_back( ntk.create_not(ntk.create_xor( xL, xR ) ) );
        break;
      default:
        break;
      }
    }
    chain.push_back( ntk.create_pi() );
  }

  for( int iOut{0}; iOut < output_c.divisors.size(); ++iOut )
  {
    divisor_t out = output_c.divisors[iOut];
    switch ( out.gate )
    {
    case gate_t::CMPL :
      ntk.create_po( ntk.create_not( chain[ out.inL ] ) );
      break;
    case gate_t::CMPR :
      ntk.create_po( ntk.create_not( chain[ out.inR ] ) );
      break;
    case gate_t::PRJL :
      ntk.create_po( chain[ out.inL ] );
      break;
    case gate_t::PRJR :
      ntk.create_po( chain[ out.inR ] );
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
  output_c.print();

}
#pragma endregion visualize

} // namespace ccgame

} // namespace mockturtle