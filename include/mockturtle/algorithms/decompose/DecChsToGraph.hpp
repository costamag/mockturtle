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
  \file ChSimnet.hpp
  \brief data structure for combining simulations and nodes

  \author Andrea Costamagna
*/
#pragma once

#include <stdio.h>
#include <stack>
#include "DecNet.hpp"
#include "../../networks/aig.hpp"

namespace mockturtle
{

template<class TT, class Ntk>
class DecChsToGraph
{
private:
    DecNet<TT, Ntk> net;
    Ntk ntk;
    std::vector<signal<Ntk>> pis;
public:
    DecChsToGraph();
    DecChsToGraph( DecNet<TT, Ntk>& );
    ~DecChsToGraph();
    /* read */
    signal_t NodeToSig( node_t );
    /* actions */
    Ntk convert();
    signal<Ntk> reconvert( signal_t );
};

template<class TT, class Ntk>   DecChsToGraph<TT, Ntk>::DecChsToGraph()
{
}

template<class TT, class Ntk>   DecChsToGraph<TT, Ntk>::DecChsToGraph(DecNet<TT, Ntk>& net) : net(net)
{
}
template<class TT, class Ntk>   DecChsToGraph<TT, Ntk>::~DecChsToGraph(){}

#pragma region read
template<class TT, class Ntk>  signal_t DecChsToGraph<TT, Ntk>::NodeToSig( node_t node){ return net.NodeToSig( node ); }
#pragma endregion

#pragma region actions

template<class TT, class Ntk>
signal<Ntk> DecChsToGraph<TT, Ntk>::reconvert( signal_t sig )
{
    std::vector<signal<Ntk>> children_ntk;   
    if( net.isSynt( sig ) == 1 )
    {
      return net.getNtkSig( sig );
    }
    else
    {
        net.foreach_fanin( sig, [&]( node_t x ) 
        {   
            signal_t child = NodeToSig(x);
            children_ntk.push_back( reconvert( child ) );
        } );

        switch ( net.getFnType( sig ) )
        {
            case DecFunc_t::NOT_:
                net.setSig( sig.node, ntk.create_not( children_ntk[0] ) );
                break;
            case DecFunc_t::BUF_:
                net.setSig( sig.node, ntk.create_buf( children_ntk[0] ) );
                break;
            case DecFunc_t::AND_:
                net.setSig( sig.node, ntk.create_and( children_ntk[0], children_ntk[1] ));
                break;
            case DecFunc_t::NAND_:
                net.setSig( sig.node, ntk.create_not(ntk.create_and( children_ntk[0], children_ntk[1] )));
                break;
            case DecFunc_t::OR_:
                net.setSig( sig.node, ntk.create_or( children_ntk[0], children_ntk[1] ));
                break;
            case DecFunc_t::NOR_:
                net.setSig( sig.node, ntk.create_not(ntk.create_or( children_ntk[0], children_ntk[1] )));
                break;
            case DecFunc_t::XOR_:
                net.setSig( sig.node, ntk.create_xor( children_ntk[0], children_ntk[1]));
                break;
            case DecFunc_t::XNOR_:
                net.setSig( sig.node, ntk.create_not( ntk.create_xor( children_ntk[0], children_ntk[1])));
                break;
            case DecFunc_t::LT_:
                net.setSig( sig.node, ntk.create_and( ntk.create_not( children_ntk[0] ), children_ntk[1] ) );
                break;
            case DecFunc_t::GE_:
                net.setSig( sig.node, ntk.create_or( children_ntk[0], ntk.create_not( children_ntk[1] ) ) );
                break;
            case DecFunc_t::LE_:
                net.setSig( sig.node, ntk.create_or( ntk.create_not( children_ntk[0] ), children_ntk[1] ) );
                break;
            case DecFunc_t::GT_:
                net.setSig( sig.node, ntk.create_and( children_ntk[0], ntk.create_not( children_ntk[1] )));
                break;
            default:
                break;
        }
    }
    return net.getNtkSig( sig ); 

}

template<class TT, class Ntk>
Ntk DecChsToGraph<TT, Ntk>::convert( )
{
    net.foreach_pi( [&]( const auto& x, auto index ) 
    {    
        pis.push_back(ntk.create_pi()); 
        net.setSig( x.node, pis[index] );
    } );
    net.foreach_po( [&]( const auto& x, auto index ) 
    {   
        ntk.create_po( reconvert( x ) ); 
    } );
    return ntk;
}

#pragma endregion
}//namespace mockturtle