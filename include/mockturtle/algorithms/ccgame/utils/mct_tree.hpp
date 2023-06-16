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
  \file mct_tree.hpp
  \brief generic monte carlo tree search engine

  \author Andrea Costamagna
*/
#pragma once

#include "ccg_net.hpp"
#include "ccg_mcnodes.hpp"
#include "ccg_rng.hpp"
#include "ccg_supportor.hpp"
#include <stdio.h>
#include <stack>
#include <iostream>

namespace mockturtle
{

namespace ccgame
{

using DTT = kitty::dynamic_truth_table;

template<class NODE, class METHOD>
class mct_tree_t
{
    public:
        std::vector<NODE> nodes;
        METHOD method;

        /* CONSTRUCT/DESCTRUCT */
        mct_tree_t( NODE, METHOD );
        mct_tree_t(){};
        ~mct_tree_t(){};

        /* GROW */
        NODE * select();
        NODE * expand( NODE * );
        NODE * simulate( NODE * );
        void backpropagate( NODE * );

        NODE * solve();

};

template<class NODE, class METHOD>
mct_tree_t<NODE, METHOD>::mct_tree_t( NODE root, METHOD method )
{
    nodes.push_back( root );
} 

template<class NODE, class METHOD>
NODE * mct_tree_t<NODE, METHOD>::select()
{
    return method.select( &nodes );
}

template<class NODE, class METHOD>
NODE * mct_tree_t<NODE, METHOD>::expand( NODE * pNdA )
{
    NODE NdB = method.expand( pNdA );
    nodes.push_back( NdB );
    return nodes.back();
}

template<class NODE, class METHOD>
NODE * mct_tree_t<class NODE, class METHOD>::simulate( NODE * pNdA )
{
    NODE NdB = method.simulate( &nodes, pNdA );
    nodes.push_back( NdB );
    return nodes.back();
}

template<class NODE, class METHOD>
void mct_tree_t<NODE, METHOD>::backpropagate( NODE * pNd )
{
    method.backpropagate( &nodes, pNd );
}

template<class NODE, class METHOD>
NODE * mct_tree_t<NODE, METHOD>::solve()
{
    NODE * pNdBest{nullptr};
    // your solution
    for( int it{0}; it<method.nIters; ++it )
    {
        NODE * pNdSel = select();
        NODE * pNdExp = expand( pNdSel );
        NODE * pNdEnd ;

        if( pNdExp->isLeaf )
            pNdEnd = pNdExp;
        else
            pNdEnd = simulate( pNdExp );

        backpropagate( pNdEnd );

        if( method.isBest( pNdEnd ) )
            pNdBest = pNdEnd;
    }
    assert( pNdBest != nullptr )
    return pNdBest;
}

} // namespace ccgame

} // namespace mockturtle