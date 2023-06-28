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
  \brief data structure for storing the cuts for the mcts

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

#pragma region parameters


#pragma endregion parameters

template<class NODE>
class mct_method_t
{
    public:
        mct_method_ps ps;

        /* CONSTRUCT/DESCTRUCT */
        mct_method_t( mct_method_ps ps ) : ps(ps) {};
        mct_method_t(){};
        ~mct_method_t(){};

        int select( std::vector<NODE> * );
        
        NODE expand( NODE * );
        NODE simulate( NODE * );
        void backpropagate( std::vector<NODE> *, int, double );

        double evaluate( std::vector<NODE *> );
};

#pragma region SELECT

template<class NODE>
int select_at_random( std::vector<NODE> * vNdPtrs  )
{
    std::uniform_int_distribution<> distrib(0, vNdPtrs->size()-1);
    return distrib(ml_gen);
}

template<class NODE>
int mct_method_t<NODE>::select( std::vector<NODE> * vNdPtrs )
{
    switch ( ps.sel_type )
    {
    case node_selection_t::NODE_RAND:
        return select_at_random( vNdPtrs );
        break;
    
    default:
        break;
    }
}
#pragma endregion SELECT

#pragma region EXPAND
template<class NODE>
NODE mct_method_t<NODE>::expand( NODE * pNd )
{
    return pNd->find_new();
}
#pragma endregion EXPAND

#pragma region SIMULATE
template<class NODE>
NODE mct_method_t<NODE>::simulate( NODE * pNd )
{
    return pNd->find_new();
}
#pragma endregion simulate

#pragma region BACKPROP
template<class NODE>
void mct_method_t<NODE>::backpropagate( std::vector<NODE> * vNdPtrs, int idEnd, double cost )
{
    NODE nd = (*vNdPtrs)[idEnd]; 
    (*vNdPtrs)[idEnd].add_cost( cost );

    while( !nd.isRoot )
    {
        nd = (*vNdPtrs)[nd.idPar];
        (*vNdPtrs)[idEnd].add_cost( cost );
    }
}
#pragma endregion BACKPROP

template<class NODE>
double mct_method_t<NODE>::evaluate( std::vector<NODE *> pNds )
{
    return pNds.back()->evaluate( pNds );
}

} // namespace mcts

} // namespace mockturtle