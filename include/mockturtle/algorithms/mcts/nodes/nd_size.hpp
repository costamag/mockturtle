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
  \file ccg_divisors.hpp
  \brief data structure for storing the divisorss for the mcts

  \author Andrea Costamagna
*/
#pragma once

#include "../ml_rng.hpp"
#include "../supportor.hpp"
#include "../mct_utils.hpp"
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <stdio.h>
#include <stack>
#include <iostream>

namespace mockturtle
{

namespace mcts
{

using DTT = kitty::dynamic_truth_table;

template<class NTK>
class nd_size_t
{
    public:
        support_generator_t supportor;
        std::vector<divisor_t> divisors;
        std::vector<target_t>  targets;
        std::vector<int> TargetsDoneHere;
        int id{0};
        int idPar{-1};
        std::vector<int> vKids;
        bool isNull;
        bool isRoot;
        bool isLeaf;
        node_ps ps;
        NTK ntk;

        /* CONSTRUCT/DESCTRUCT */
        /*! \brief generic node constructor */
        nd_size_t( std::vector<divisor_t>, std::vector<target_t>, node_ps );
        /*! \brief root node constructor*/
        nd_size_t( std::vector<DTT>, std::vector<double>, std::vector<DTT>, node_ps );
        /* default */
        nd_size_t(){ isNull = true; };
        ~nd_size_t(){};
        /* PROPERTIES */
        void print();
        bool is_root();
        bool is_leaf();
        bool is_null();
        /* GROW */
        nd_size_t find_new();
        nd_size_t null_node();
        void add_child( int );
        /* SYNTHESIZE */
        double evaluate( std::vector<nd_size_t<NTK> *> );
        bool check_closure();
};



template<class NTK>
nd_size_t<NTK>::nd_size_t( std::vector<divisor_t> X, std::vector<target_t> Y, node_ps ps ) : ps(ps)
{
    isNull = false;
    isRoot = false;
    for( int iTrg{0}; iTrg<Y.size(); ++iTrg )  
        targets.push_back(Y[iTrg]);
    for( int iDiv{0}; iDiv<X.size(); ++iDiv )
    {
        divisors.push_back( X[iDiv] );
        divisors[iDiv].id = iDiv;
    }
    isLeaf = check_closure();
    supportor = support_generator_t{ &divisors, &targets, ps };
}

template<class NTK>
nd_size_t<NTK>::nd_size_t( std::vector<DTT> X, std::vector<double> T, std::vector<DTT> Y, node_ps ps ):ps(ps)
{
    assert( X.size() == T.size() );
    isNull = false;
    isRoot = true;
    for( int iTrg{0}; iTrg<Y.size(); ++iTrg )  
        targets.emplace_back(iTrg, Y[iTrg]);
    
    for( int iDiv{0}; iDiv<X.size(); ++iDiv )
        divisors.emplace_back( iDiv, X[iDiv], 0.0, T[iDiv], gate_t::PIS );
    
    isLeaf = check_closure();
    supportor = support_generator_t{ &divisors, &targets, ps };
}

template<class NTK>
bool nd_size_t<NTK>::check_closure()
{
    bool isClosed{true};
    for( int iTrg{0}; iTrg<targets.size(); ++iTrg )
    {
        if( targets[iTrg].isDone ) continue;
        for( int iDiv{0}; iDiv<divisors.size(); ++iDiv )
        {
            targets[iTrg].isDone = true;
            targets[iTrg].div    = iDiv;
            if( kitty::equal( targets[iTrg].tt, divisors[iDiv].tt ) )
            {
                divisors[iDiv].isPo = true;
                targets[iTrg].type = gate_t::PRJL;
                targets[iTrg].isDone = true;
                break;
            }
            else if( kitty::equal( targets[iTrg].tt, ~divisors[iDiv].tt ) )
            {
                divisors[iDiv].isPo = true;
                targets[iTrg].type   = gate_t::CMPL;
                targets[iTrg].isDone = true;
                break;
            }
            else
                targets[iTrg].isDone = false;
        }
        isClosed &= targets[iTrg].isDone;
        if( targets[iTrg].isDone ) TargetsDoneHere.push_back(iTrg);
    }
    return isClosed;
}

#pragma region PROPERTIES
template<class NTK> bool nd_size_t<NTK>::is_null(){ return isNull; }
template<class NTK> bool nd_size_t<NTK>::is_root(){ return isRoot; }
template<class NTK> bool nd_size_t<NTK>::is_leaf(){ return isLeaf; }
#pragma endregion PROPERTIES

#pragma region GROW
template<class NTK>
nd_size_t<NTK> nd_size_t<NTK>::find_new()
{
    std::vector<int> supp = supportor.find_new( 1 );
    if( supp.size() == 0 )  {return null_node();}
    std::vector<divisor_t> divs;
    for( auto s : supp )    divs.push_back( supportor.divisors[s] );

    nd_size_t node( divs, supportor.targets, ps );
    return node;
}

template<class NTK>
void nd_size_t<NTK>::add_child( int idChild )
{
    vKids.push_back( idChild );
}
#pragma endregion GROW

template<class NTK>
void nd_size_t<NTK>::print()
{
    printf("=============================\n");
    for( int i{0}; i<divisors.size(); ++i )
        divisors[i].print();
}

template<class NTK>
nd_size_t<NTK> nd_size_t<NTK>::null_node( )
{
    nd_size_t<NTK> node0;
    node0.isLeaf=false;
    node0.isRoot=false;
    node0.isNull=true;
    return node0;
}

template<class NTK>
double nd_size_t<NTK>::evaluate( std::vector<nd_size_t<NTK>*> vPtrs )
{
    NTK net;
    std::vector<signal<NTK>> sigs_old;
    std::vector<signal<NTK>> sigs_new;
    std::vector<signal<NTK>> outSigs;
    /* deal with pis */
    nd_size_t<NTK> * pNd = vPtrs[0];
    assert( pNd->idPar == -1 );
    for( int iPi{0}; iPi<pNd->divisors.size(); ++iPi )
        sigs_old.push_back( net.create_pi() );
    
    if( pNd->TargetsDoneHere.size() > 0 )
    {
        for( auto iTrg : pNd->TargetsDoneHere )
        {
            int idDiv = pNd->supportor.targets[iTrg].div;
            if( pNd->targets[iTrg].type == gate_t::CMPL )
                outSigs.insert( outSigs.begin()+iTrg,net.create_not(sigs_old[idDiv]));
            else if( pNd->targets[iTrg].type == gate_t::PRJL )
                outSigs.insert( outSigs.begin()+iTrg, sigs_old[idDiv]);
            else
                assert(0);
        }
    }

    for( int iLev{1}; iLev < vPtrs.size(); ++iLev )
    {
        for( int iDiv{0}; iDiv < vPtrs[iLev]->divisors.size(); ++iDiv )
        {
            switch ( vPtrs[iLev]->divisors[iDiv].type )
            {
            case gate_t::AI00:
                sigs_new.push_back( net.create_and(!sigs_old[vPtrs[iLev]->divisors[iDiv].fanins[1]],!sigs_old[vPtrs[iLev]->divisors[iDiv].fanins[0]] ));
                break;
            case gate_t::AI01:
                sigs_new.push_back( net.create_and(!sigs_old[vPtrs[iLev]->divisors[iDiv].fanins[1]], sigs_old[vPtrs[iLev]->divisors[iDiv].fanins[0]] ));
                break;
            case gate_t::AI10:
                sigs_new.push_back( net.create_and( sigs_old[vPtrs[iLev]->divisors[iDiv].fanins[1]],!sigs_old[vPtrs[iLev]->divisors[iDiv].fanins[0]] ));
                break;
            case gate_t::AI11:
                sigs_new.push_back( net.create_and( sigs_old[vPtrs[iLev]->divisors[iDiv].fanins[1]], sigs_old[vPtrs[iLev]->divisors[iDiv].fanins[0]] ));
                break;
            case gate_t::PRJL:
                sigs_new.push_back( sigs_old[vPtrs[iLev]->divisors[iDiv].fanins[1]]);
                break;
            case gate_t::PRJR:
                sigs_new.push_back( sigs_old[vPtrs[iLev]->divisors[iDiv].fanins[0]]);
                break;
            case gate_t::CMPL:
                sigs_new.push_back( !sigs_old[vPtrs[iLev]->divisors[iDiv].fanins[1]]);
                break;
            case gate_t::CMPR:
                sigs_new.push_back( !sigs_old[vPtrs[iLev]->divisors[iDiv].fanins[0]]);
                break;
            default:
                break;
            }
        }
        sigs_old = sigs_new;
        sigs_new = {};
        if( vPtrs[iLev]->TargetsDoneHere.size() > 0 )
        {
            for( auto iTrg : vPtrs[iLev]->TargetsDoneHere )
            {
                int idDiv = vPtrs[iLev]->supportor.targets[iTrg].div;
                if( vPtrs[iLev]->targets[iTrg].type == gate_t::CMPL )
                    outSigs.insert( outSigs.begin()+iTrg,net.create_not(sigs_old[idDiv]));
                else if( vPtrs[iLev]->targets[iTrg].type == gate_t::PRJL )
                    outSigs.insert( outSigs.begin()+iTrg, sigs_old[idDiv]);
                else
                    assert(0);
            }
        }
    }
    /* synthesize the outouts */
    for( int iTrg{0}; iTrg<pNd->supportor.targets.size(); ++iTrg )
    {
        net.create_po(outSigs[iTrg]);
    }
    ntk = cleanup_dangling(net);

    return ntk.num_gates();
}

} // namespace mcts

} // namespace mockturtle