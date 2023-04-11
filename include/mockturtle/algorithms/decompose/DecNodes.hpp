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
  \brief data structure for storing the Sims

  \author Andrea Costamagna
*/
#pragma once

#include <stdio.h>
#include <stack>
#include "../../networks/aig.hpp"

namespace mockturtle
{

enum DecFunc_t
{
  NONE_,
  PI_,
  PO_,
  HUNG_,
  NOT_,
  BUF_,
  AND_,
  NAND_,
  XOR_,
  XNOR_,
  OR_,
  NOR_,
  LT_,
  GE_,
  GT_,
  LE_
};

using node_t = uint32_t;
using sim_t  = uint32_t;

template<class Ntk>
class DecNodes
{
private:
  std::vector<std::vector<sim_t>> vInSims;    // fanins-simulation of each node
  std::vector<std::vector<node_t>>vFanIns;    // fanins-simulation of each node
  std::vector<sim_t>              vSims;       // simulation of the node
  std::vector<signal<Ntk>>        vSigs;      // vector of network signals
  std::vector<DecFunc_t>          vFnTypes;   // function type of each node
  std::vector<bool>               vUsed;      // function type of each node
  std::vector<int>                vSynt;      // function type of each node
  std::stack<node_t>              sFree;      // stack of free node
  int                             nNodes; 

public:
  DecNodes();
  ~DecNodes();
  /* modify */
  node_t addNode( std::vector<sim_t>, sim_t, DecFunc_t );
  node_t addHungNode( sim_t );
  void attachHunging( node_t, sim_t, node_t, DecFunc_t );

  void rmNode( node_t );
  void setSig( node_t, signal<Ntk> );
  /* read */
  bool isUsed( node_t );
  int isSynt( node_t );
  bool isPI( node_t );
  int  size();
  std::vector<sim_t>  * getInSimsP( node_t );
  std::vector<node_t> * getFanInsP( node_t );
  sim_t                getSim( node_t );
  DecFunc_t            getFunc( node_t );
  signal<Ntk>          getNtkSig( node_t );
};

#pragma region constructors
template<class Ntk>
DecNodes<Ntk>::DecNodes()
{
  nNodes = 0;
}

template<class Ntk>
DecNodes<Ntk>::~DecNodes()
{
}
#pragma endregion

#pragma region read
template<class Ntk> int DecNodes<Ntk>::size(){ return nNodes;  }
template<class Ntk> bool DecNodes<Ntk>::isUsed( node_t ref ){ return ( ref<vInSims.size() & vUsed[ref] ); };

template<class Ntk> int DecNodes<Ntk>::isSynt( node_t ref ){   return vSynt[ref]; };

template<class Ntk> bool DecNodes<Ntk>::isPI( node_t ref ){ return ( vFnTypes[ref] == DecFunc_t::PI_ ); };
template<class Ntk>
std::vector<sim_t> * DecNodes<Ntk>::getInSimsP( node_t ref )
{
  assert( ref < vInSims.size() ); 
  assert( vUsed[ref] ); 
  return &(vInSims[ref]); 
}

template<class Ntk>
std::vector<sim_t> * DecNodes<Ntk>::getFanInsP( node_t ref )
{
  assert( ref < vInSims.size() ); 
  assert( vUsed[ref] ); 
  return &(vFanIns[ref]); 
}

template<class Ntk>
sim_t DecNodes<Ntk>::getSim( node_t ref )
{ 
  assert( ref < vInSims.size() ); 
  assert( vUsed[ref] ); 
  return vSims[ref]; 
};
template<class Ntk>
DecFunc_t DecNodes<Ntk>::getFunc( node_t ref )
{ 
  assert( ref < vInSims.size() ); 
  assert( vUsed[ref] ); 
  return vFnTypes[ref]; 
}

template<class Ntk> signal<Ntk> DecNodes<Ntk>::getNtkSig( node_t node ){ return vSigs[node]; };

#pragma endregion

#pragma region modify
template<class Ntk>
node_t DecNodes<Ntk>::addNode( std::vector<node_t> Fanins, sim_t pSim, DecFunc_t pFunc )
{

  std::vector<sim_t> InSims;
  for( auto x : Fanins )
  {
    InSims.push_back( getSim(x) );
  }
  assert( vInSims.size() == vSims.size() );
  assert( vInSims.size() == vFnTypes.size() );
  assert( vInSims.size() == vUsed.size() );
  assert( vInSims.size() == vSynt.size() );
  int ref;
  if( sFree.size() == 0 )
  {
    assert( nNodes == vInSims.size() );
    ref = nNodes;
    vInSims.push_back( InSims );
    vFanIns.push_back( Fanins );
    vSims.push_back( pSim );
    vUsed.push_back( true );
    vSynt.push_back( 0 );
    vFnTypes.push_back( pFunc );
    vSigs.emplace_back();
  }
  else
  {
    ref = sFree.top();
    vInSims[ref] = InSims;
    vFanIns[ref] = Fanins;
    vSims[ref] = pSim;
    vUsed[ref] = true;
    vSynt[ref] = 0;
    vFnTypes[ref] = pFunc;
    sFree.pop();
  }
  nNodes++;
  return ref;
}

template<class Ntk>
void DecNodes<Ntk>::attachHunging( node_t inNode, sim_t inSim, node_t out, DecFunc_t FuncT )
{
  vFanIns[out].push_back(inNode);
  vInSims[out].push_back(inSim);
  vFnTypes[out] = FuncT;
}

template<class Ntk>
node_t DecNodes<Ntk>::addHungNode( sim_t pSim )
{
  assert( vInSims.size() == vSims.size() );
  assert( vInSims.size() == vFnTypes.size() );
  assert( vInSims.size() == vUsed.size() );
  assert( vInSims.size() == vSynt.size() );
  int ref;
  if( sFree.size() == 0 )
  {
    assert( nNodes == vInSims.size() );
    ref = nNodes;
    vInSims.emplace_back();
    vFanIns.emplace_back();
    vSims.push_back( pSim );
    vUsed.push_back( true );
    vSynt.push_back( 0 );
    vFnTypes.push_back( DecFunc_t::HUNG_ );
    vSigs.emplace_back();
  }
  else
  {
    ref = sFree.top();
    vInSims[ref] = {};
    vFanIns[ref] = {};
    vSims[ref] = pSim;
    vUsed[ref] = true;
    vSynt[ref] = 0;
    vFnTypes[ref] = DecFunc_t::HUNG_;
    sFree.pop();
  }
  nNodes++;
  return ref;
}

template<class Ntk>
void DecNodes<Ntk>::rmNode( node_t ref )
{
  assert(nNodes>0);
  assert(vInSims.size()>ref);
  assert(vFnTypes.size()>ref);
  assert(vSims.size()>ref);
  assert(vUsed.size()>ref);
  assert(vUsed[ref]);
  vSims[ref] = 0;
  vUsed[ref] = false;
  vSynt[ref] = 0;
  vFnTypes[ref] = DecFunc_t::NONE_;
  sFree.push(ref);
  nNodes--;
}

template<class Ntk>
void DecNodes<Ntk>::setSig( node_t node, signal<Ntk> ntk_sig )
{ 
  vSigs[node] = ntk_sig; 
  vSynt[node] = 1;
}
#pragma endregion modify

} // namespace mockturtle