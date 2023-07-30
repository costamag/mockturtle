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
#include <limits>
#include <stack>
#include <iostream>

namespace mockturtle
{

namespace mcts
{

using DTT = kitty::dynamic_truth_table;

template<class NTK>
class nd_delay_t
{
    public:
        support_generator_t supportor;
        std::vector<divisor_t> divisors;
        std::vector<target_t>  targets;
        std::vector<int> TargetsDoneHere;
        std::vector<double> costs;
        double bestCost = std::numeric_limits<double>::max();
        int id{0};
        int idPar{-1};
        std::vector<int> vKids;
        bool isNull;
        bool isRoot;
        bool isLeaf;
        node_ps ps;
        NTK ntk;
        /* BACKPROP PARAMS */
        double ni = 0;
        double wi = 0;
        double Ni = 0; 
        NTK * pNtkOUT;

        /* CONSTRUCT/DESCTRUCT */
        /*! \brief generic node constructor */
        nd_delay_t( std::vector<divisor_t>, std::vector<target_t>, node_ps );
        /*! \brief root node constructor*/
        nd_delay_t( std::vector<DTT>, std::vector<double>, std::vector<DTT>, node_ps );
        nd_delay_t( std::vector<DTT>, std::vector<double>, std::vector<DTT>, node_ps, NTK * );
        /* default */
        nd_delay_t(){ isNull = true; };
        ~nd_delay_t(){};
        /* PROPERTIES */
        void print();
        bool is_root();
        bool is_leaf();
        bool is_null();
        /* GROW */
        nd_delay_t find_new();
        nd_delay_t null_node();
        void add_child( int );
        /* SYNTHESIZE */
        double evaluate( std::vector<nd_delay_t<NTK> *> );
        bool check_closure();
        void add_cost( double cost );
        void update_support_info( nd_delay_t<NTK>, double );
        signal<NTK> implant( std::vector<signal<NTK>>, std::vector<nd_delay_t<NTK>> );
};



template<class NTK>
nd_delay_t<NTK>::nd_delay_t( std::vector<divisor_t> X, std::vector<target_t> Y, node_ps eps )
{
    ps=eps;
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
nd_delay_t<NTK>::nd_delay_t( std::vector<DTT> X, std::vector<double> T, std::vector<DTT> Y, node_ps eps )
{
    ps = eps;
    assert( X.size() == T.size() );
    isNull = false;
    isRoot = true;
    for( int iTrg{0}; iTrg<Y.size(); ++iTrg )  
        targets.emplace_back( ps.use_inf_graph, iTrg, Y[iTrg]);
    
    for( int iDiv{0}; iDiv<X.size(); ++iDiv )
        divisors.emplace_back(  ps.use_inf_graph, iDiv, X[iDiv], 0.0, T[iDiv], gate_t::PIS );
    
    isLeaf = check_closure();
    supportor = support_generator_t{ &divisors, &targets, ps };
}

template<class NTK>
nd_delay_t<NTK>::nd_delay_t( std::vector<DTT> X, std::vector<double> T, std::vector<DTT> Y, node_ps eps, NTK * pNTKOUT )
{
    pNtkOUT = pNTKOUT;
    ps = eps;
    assert( X.size() == T.size() );
    isNull = false;
    isRoot = true;
    for( int iTrg{0}; iTrg<Y.size(); ++iTrg )  
        targets.emplace_back( ps.use_inf_graph, iTrg, Y[iTrg]);
    
    for( int iDiv{0}; iDiv<X.size(); ++iDiv )
        divisors.emplace_back(  ps.use_inf_graph, iDiv, X[iDiv], 0.0, T[iDiv], gate_t::PIS );
    
    isLeaf = check_closure();
    supportor = support_generator_t{ &divisors, &targets, ps };
}

template<class NTK>
bool nd_delay_t<NTK>::check_closure()
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
template<class NTK> bool nd_delay_t<NTK>::is_null(){ return isNull; }
template<class NTK> bool nd_delay_t<NTK>::is_root(){ return isRoot; }
template<class NTK> bool nd_delay_t<NTK>::is_leaf(){ return isLeaf; }
#pragma endregion PROPERTIES

#pragma region GROW
template<class NTK>
nd_delay_t<NTK> nd_delay_t<NTK>::find_new()
{
    std::vector<int> supp;
    
    switch (ps.sel_type)
    {
    case supp_selection_t::SUP_ENER:
        supp = supportor.find_new<supp_selection_t::SUP_ENER>( ps.nIters );
        break;
    case supp_selection_t::SUP_DECT:
        supp = supportor.find_new<supp_selection_t::SUP_DECT>( ps.nIters );
        break;
    default:
        break;
    }

        
    if( supp.size() == 0 )  {return null_node();}
    std::vector<divisor_t> divs;
    for( auto s : supp )    divs.push_back( supportor.divisors[s] );

    nd_delay_t node( divs, supportor.targets, ps );
    return node;
}

template<class NTK>
void nd_delay_t<NTK>::add_child( int idChild )
{
    vKids.push_back( idChild );
}
#pragma endregion GROW

template<class NTK>
void nd_delay_t<NTK>::print()
{
    printf("=============================\n");
    for( int i{0}; i<divisors.size(); ++i )
        divisors[i].print();
    printf( "costs: ");
    for( int i{0}; i<costs.size(); ++i )
        printf( "%f ", costs[i] );
    printf("\n");
}

template<class NTK>
nd_delay_t<NTK> nd_delay_t<NTK>::null_node( )
{
    nd_delay_t<NTK> node0;
    node0.isLeaf=false;
    node0.isRoot=false;
    node0.isNull=true;
    return node0;
}

template<class NTK>
signal<NTK> nd_delay_t<NTK>::implant( std::vector<signal<NTK>> S, std::vector<nd_delay_t<NTK>> path )
{
    NTK * pNet = path[0].pNtkOUT;
    std::vector<signal<NTK>> sigs_old = S;
    std::vector<signal<NTK>> sigs_new;
    std::vector<signal<NTK>> outSigs;
    /* deal with pis */
    nd_delay_t<NTK> * pNd = &path[0];
    assert( pNd->idPar == -1 );
    
    if( pNd->TargetsDoneHere.size() > 0 )
    {
        for( auto iTrg : pNd->TargetsDoneHere )
        {
            int idDiv = pNd->supportor.targets[iTrg].div;
            if( pNd->targets[iTrg].type == gate_t::CMPL )
                outSigs.insert( outSigs.begin()+iTrg,pNet->create_not(sigs_old[idDiv]));
            else if( pNd->targets[iTrg].type == gate_t::PRJL )
                outSigs.insert( outSigs.begin()+iTrg, sigs_old[idDiv]);
            else if( pNd->targets[iTrg].type == gate_t::CMPR )
                outSigs.insert( outSigs.begin()+iTrg,pNet->create_not(sigs_old[idDiv]));
            else if( pNd->targets[iTrg].type == gate_t::PRJR )
                outSigs.insert( outSigs.begin()+iTrg, sigs_old[idDiv]);
            else
                assert(0);
        }
    }

    for( int iLev{1}; iLev < path.size(); ++iLev )
    {
        for( int iDiv{0}; iDiv < path[iLev].divisors.size(); ++iDiv )
        {
            switch ( path[iLev].divisors[iDiv].type )
            {
            case gate_t::AI00:
                sigs_new.push_back( pNet->create_and(!sigs_old[path[iLev].divisors[iDiv].fanins[1]],!sigs_old[path[iLev].divisors[iDiv].fanins[0]] ));
                break;
            case gate_t::AI01:
                sigs_new.push_back( pNet->create_and(!sigs_old[path[iLev].divisors[iDiv].fanins[1]], sigs_old[path[iLev].divisors[iDiv].fanins[0]] ));
                break;
            case gate_t::AI10:
                sigs_new.push_back( pNet->create_and( sigs_old[path[iLev].divisors[iDiv].fanins[1]],!sigs_old[path[iLev].divisors[iDiv].fanins[0]] ));
                break;
            case gate_t::AI11:
                sigs_new.push_back( pNet->create_and( sigs_old[path[iLev].divisors[iDiv].fanins[1]], sigs_old[path[iLev].divisors[iDiv].fanins[0]] ));
                break;
            case gate_t::EXOR:
                sigs_new.push_back( pNet->create_xor( sigs_old[path[iLev].divisors[iDiv].fanins[1]], sigs_old[path[iLev].divisors[iDiv].fanins[0]] ));
                break;
            case gate_t::PRJL:
                sigs_new.push_back( sigs_old[path[iLev].divisors[iDiv].fanins[1]]);
                break;
            case gate_t::PRJR:
                sigs_new.push_back( sigs_old[path[iLev].divisors[iDiv].fanins[0]]);
                break;
            case gate_t::CMPL:
                sigs_new.push_back( !sigs_old[path[iLev].divisors[iDiv].fanins[1]]);
                break;
            case gate_t::CMPR:
                sigs_new.push_back( !sigs_old[path[iLev].divisors[iDiv].fanins[0]]);
                break;
            default:
                break;
            }
        }
        sigs_old = sigs_new;
        sigs_new = {};
        if( path[iLev].TargetsDoneHere.size() > 0 )
        {
            for( auto iTrg : path[iLev].TargetsDoneHere )
            {
                //printf("TYPE OUT %d\n", path[iLev].targets[iTrg].type );
                int idDiv = path[iLev].supportor.targets[iTrg].div;
                if( path[iLev].targets[iTrg].type == gate_t::CMPL )
                    outSigs.insert( outSigs.begin()+iTrg,pNet->create_not(sigs_old[idDiv]));
                else if( path[iLev].targets[iTrg].type == gate_t::PRJL )
                    outSigs.insert( outSigs.begin()+iTrg, sigs_old[idDiv]);
                else if( path[iLev].targets[iTrg].type == gate_t::CMPR )
                    outSigs.insert( outSigs.begin()+iTrg,pNet->create_not(sigs_old[idDiv]));
                else if( path[iLev].targets[iTrg].type == gate_t::PRJR )
                    outSigs.insert( outSigs.begin()+iTrg, sigs_old[idDiv]);
                else
                    assert(0);
            }
        }
    }
    /* synthesize the outouts */
    return outSigs[0];
}


template<class NTK>
double nd_delay_t<NTK>::evaluate( std::vector<nd_delay_t<NTK>*> vPtrs )
{
    NTK net;
    std::vector<signal<NTK>> sigs_old;
    std::vector<signal<NTK>> sigs_new;
    std::vector<signal<NTK>> outSigs;
    /* deal with pis */
    nd_delay_t<NTK> * pNd = vPtrs[0];
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
            else if( pNd->targets[iTrg].type == gate_t::CMPR )
                outSigs.insert( outSigs.begin()+iTrg,net.create_not(sigs_old[idDiv]));
            else if( pNd->targets[iTrg].type == gate_t::PRJR )
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
            case gate_t::EXOR:
                sigs_new.push_back( net.create_xor( sigs_old[vPtrs[iLev]->divisors[iDiv].fanins[1]], sigs_old[vPtrs[iLev]->divisors[iDiv].fanins[0]] ));
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
                //printf("TYPE OUT %d\n", vPtrs[iLev]->targets[iTrg].type );
                int idDiv = vPtrs[iLev]->supportor.targets[iTrg].div;
                if( vPtrs[iLev]->targets[iTrg].type == gate_t::CMPL )
                    outSigs.insert( outSigs.begin()+iTrg,net.create_not(sigs_old[idDiv]));
                else if( vPtrs[iLev]->targets[iTrg].type == gate_t::PRJL )
                    outSigs.insert( outSigs.begin()+iTrg, sigs_old[idDiv]);
                else if( vPtrs[iLev]->targets[iTrg].type == gate_t::CMPR )
                    outSigs.insert( outSigs.begin()+iTrg,net.create_not(sigs_old[idDiv]));
                else if( vPtrs[iLev]->targets[iTrg].type == gate_t::PRJR )
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

    return vPtrs.back()->divisors[vPtrs.back()->supportor.targets[0].div].delay;
}

template<class NTK>
void nd_delay_t<NTK>::add_cost( double cost )
{
    costs.push_back( cost );
    if( cost < bestCost )
        bestCost = cost;  
}

template<class NTK>
void nd_delay_t<NTK>::update_support_info( nd_delay_t<NTK> child, double cost )
{
    int idx;
    for( int i=0; i<vKids.size(); ++i  )
    {
        if( child.id == vKids[i] )
        {
            idx = i;
            break;
        }
    }
    supportor.history.update_cost( idx, cost );
}

} // namespace mcts

} // namespace mockturtle