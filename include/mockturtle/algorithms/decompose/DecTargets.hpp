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
  \file Dec_Nodes.hpp
  \brief data structure for the graph verteces

  \author Andrea Costamagna
*/

#include <stdio.h>
#include <stack>

namespace mockturtle
{

template<class TT>
class DecTargets
{
private:
  std::vector<TT>   vFuncs;
  std::vector<TT>   vMasks;
  std::vector<bool> vUsed;
  std::stack<int>   sFree;
  int nTargets;

public:
  DecTargets();
  ~DecTargets();
  int size();
  /* modify */
  int insert( const TT&, const TT& );
  void remove( int );
  /* read */
  TT * getFuncP( int );
  TT * getMaskP( int );

};

#pragma region constructors
template<class TT>
DecTargets<TT>::DecTargets()
{
  nTargets = 0;
}

template<class TT>
DecTargets<TT>::~DecTargets()
{
}
#pragma endregion

#pragma region read
template<class TT> int  DecTargets<TT>::size(){ return nTargets;  }
template<class TT> TT * DecTargets<TT>::getFuncP( int ref ){  return &vFuncs[ref];  }
template<class TT> TT * DecTargets<TT>::getMaskP( int ref ){  return &vMasks[ref];  }
#pragma endregion

#pragma region modify
template<class TT>
int DecTargets<TT>::insert( const TT& func, const TT& mask )
{
  assert( vFuncs.size() == vMasks.size() );
  assert( vUsed.size()  == vMasks.size() );
  int ref;
  if( sFree.size() == 0 )
  {
    assert( nTargets == vFuncs.size() );
    ref = nTargets;
    vFuncs.push_back( func );
    vMasks.push_back( mask );
    vUsed.push_back( true );
  }
  else
  {
    ref = sFree.top();
    vFuncs[ref] = func;
    vMasks[ref] = mask;
    vUsed[ref] = true;
    sFree.pop();
  }
  nTargets++;
  return ref;
}

template<class TT>
void DecTargets<TT>::remove( int ref )
{
  assert(nTargets>0);
  assert(vFuncs.size()>ref);
  assert(vMasks.size()>ref);
  assert(vUsed.size()>ref);
  assert(vUsed[ref]);
  vUsed[ref] = false;
  vFuncs[ref] &= ~vFuncs[ref];
  vMasks[ref] |= ~vMasks[ref];
  sFree.push(ref);
  nTargets--;
}
#pragma endregion modify





} // namespace mockturtle