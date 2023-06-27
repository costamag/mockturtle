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

#include "ml_rng.hpp"
#include <stdio.h>
#include <stack>
#include <iostream>

namespace mockturtle
{

namespace mcts
{

using DTT = kitty::dynamic_truth_table;

struct mct_ps
{
    int nIters = 1;
    int nSims = 1;
    bool verbose = false;
};

template<class NODE, template<class> class METHOD>
class mct_tree_t
{
    public:
        std::vector<NODE> nodes;
        METHOD<NODE> method;
        mct_ps ps;

        /* CONSTRUCT/DESCTRUCT */
        mct_tree_t( NODE, METHOD<NODE>, mct_ps );
        mct_tree_t(){};
        ~mct_tree_t(){};

        /* GROW */
        int add_node( int, NODE );
        int select();
        int expand( int );
        int simulate( int );
        void backpropagate( int );

        int solve();
        void path_print( int );

        double evaluate( int );
};

template<class NODE, template<class> class METHOD>
mct_tree_t<NODE, METHOD>::mct_tree_t( NODE root, METHOD<NODE> method, mct_ps ps ) : method(method), ps(ps)
{
    nodes.push_back( root );
}

template<class NODE, template<class> class METHOD>
int mct_tree_t<NODE, METHOD>::select()
{
    return method.select( &nodes );
}

template<class NODE, template<class> class METHOD>
int mct_tree_t<NODE, METHOD>::add_node( int id, NODE child )
{
    if( child.is_null() )   return -1;

    child.id = nodes.size();
    nodes.push_back( child );
    nodes[id].add_child( child.id );
    nodes[child.id].idPar = id;
    return child.id;
}

template<class NODE, template<class> class METHOD>
int mct_tree_t<NODE, METHOD>::expand( int id )
{
    NODE nd = method.expand( &nodes[id] );
    return add_node( id, nd );
}

template<class NODE, template<class> class METHOD>
int mct_tree_t<NODE, METHOD>::simulate( int id )
{
    NODE nd;
    while( !nodes[id].is_leaf() && !nodes[id].is_null() )
    {
        nd = method.simulate( &nodes[id] );
        if( nd.is_null() ) return -1;
        id = add_node( id, nd );
    }
    return id;
}

template<class NODE, template<class> class METHOD>
void mct_tree_t<NODE, METHOD>::backpropagate( int id )
{
    method.backpropagate( &nodes, &nodes[id] );
}

template<class NODE, template<class> class METHOD>
int mct_tree_t<NODE, METHOD>::solve()
{
    int idBest = -1;
    double bestCost = 100000;
    for( int it{0}; it<ps.nIters; ++it )
    {
        bool FoundLeaf{false};
        int idEnd;
        if(ps.verbose)  printf("iter %d\n", it );
        int idSel = select();
        if( nodes[idSel].is_leaf() ) { FoundLeaf = true; idEnd = idSel;}
        if( nodes[idSel].is_null() ) { continue; }
        int idExp;
        if( !FoundLeaf )
        {
            idExp = expand( idSel );
            if( nodes[idExp].is_leaf() ) { FoundLeaf = true; idEnd = idExp; }
            if( nodes[idExp].is_null() ) { continue; }
        }
        if( FoundLeaf )
        {
            double cost = evaluate( idEnd );
            if(ps.verbose)  printf("cost %f\n", cost);
            if( cost >= 0 && cost < bestCost )
            {
                idBest = idEnd;
                bestCost = cost;
            }
        }
        else
        {
            for( int itSim{0}; itSim<ps.nSims; ++itSim )
            {
                idEnd = simulate( idExp );
                if( idEnd < 0 ) { continue; }
                //backpropagate( idEnd );
                double cost = evaluate( idEnd );
                if(ps.verbose)  printf("cost %f\n", cost);
                if( cost >= 0 && cost < bestCost )
                {
                    idBest = idEnd;
                    bestCost = cost;
                }
            }
        }
    }
    return idBest;
}

template<class NODE, template<class> class METHOD>
double mct_tree_t<NODE, METHOD>::evaluate( int id )
{
    for( auto trg : nodes[id].supportor.targets )
        if( !trg.isDone )   return -1;
    std::vector<NODE*> path;
    do
    {
        path.insert( path.begin(), &nodes[id] );
        id = nodes[id].idPar;
    } while ( id >= 0 );
    
    return method.evaluate( path );
}

template<class NODE, template<class> class METHOD>
void mct_tree_t<NODE, METHOD>::path_print( int id )
{
    if( id == -1 ) return;
    path_print( nodes[id].idPar );
    printf("=============================\n");
    printf("NODE %d\n", id);
    nodes[id].print();

} 

}// namespace mcts

} // namespace mockturtle