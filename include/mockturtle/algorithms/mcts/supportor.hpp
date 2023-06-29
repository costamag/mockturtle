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
  \file ccg_supporter.hpp
  \brief data structure for generating new supports given a set of divisors

  \author Andrea Costamagna
*/
#pragma once

#include "ml_rng.hpp"
#include "mct_utils.hpp"

#include <kitty/print.hpp>
#include <stdio.h>
#include <stack>
#include <set>
#include <math.h>
#include <limits>

namespace mockturtle
{

namespace mcts
{

struct suppor_history_t
{
    std::set<std::vector<int>> set;
    std::vector<std::vector<int>> list;
    std::vector<double> costs;

    void insert( std::vector<int> gene )
    {
        set.insert( gene );
        list.push_back( gene );
        costs.push_back( std::numeric_limits<double>::max() );
    }

    std::set<std::vector<int>>::iterator find_in_set( std::vector<int> gene )
    {
        return set.find( gene );
    }

    std::set<std::vector<int>>::iterator end_of_set()
    {
        return set.end();
    }

    void update_cost( int idx, double cost )
    {
        if( cost < costs[idx] )
            costs[idx] = cost;
    }
};

class support_generator_t
{
    public:
        std::vector<divisor_t> divisors;
        std::vector<target_t> targets;
        std::vector<int> TargetsDoneHere;
        node_ps ps;
        int nIdentity;
        //std::set<std::vector<int>> history;
        suppor_history_t history;

        /* CONSTRUCT */
        support_generator_t( std::vector<divisor_t> *, std::vector<target_t> *, node_ps );
        support_generator_t( std::vector<divisor_t>, std::vector<target_t>, node_ps, int ); // will be removed
        support_generator_t(){};
        ~support_generator_t(){};

        /* GROW */
        template<supp_selection_t SEL_t>
        std::vector<int> find_new( int );

        void store_new( std::vector<int> );
        void print();
        void mark_closing_divisors();
        void next_layer( std::vector<divisor_t> *, std::vector<target_t> * );
        void add_cost( int, double );
};

support_generator_t::support_generator_t( std::vector<divisor_t> divisors, std::vector<target_t> targets, node_ps ps, int nIdentity ):
divisors(divisors), targets(targets), ps(ps), nIdentity(nIdentity)
{
    std::vector<int> identity;
    for( int i{0}; i < nIdentity; ++i ) identity.push_back( i );
    history.insert( identity );
    for( int i{nIdentity}; i < divisors.size(); ++i ) identity.push_back( i );

    std::vector<DTT> targets_tt;
    for( auto f : targets )
        targets_tt.push_back( f.graph );
};


void support_generator_t::next_layer( std::vector<divisor_t> * pDivs0, std::vector<target_t> * pTrgs0 )
{
    targets = *pTrgs0;

    std::vector<int> fanins;
    for( int i{0}; i < pDivs0->size(); ++i )
    {
        fanins = { (*pDivs0)[i].id, (*pDivs0)[i].id };
        divisors.emplace_back( divisors.size(), (*pDivs0)[i].tt, 0, (*pDivs0)[i].delay, gate_t::PRJL, fanins );
    }

    for( int iR{0}; iR < pDivs0->size(); ++iR )
    {
        for( int iL{iR+1}; iL < pDivs0->size(); ++iL )
        {
            fanins = { (*pDivs0)[iR].id, (*pDivs0)[iL].id };
            DTT ttL = (*pDivs0)[iL].tt;
            DTT ttR = (*pDivs0)[iR].tt;
            DTT tt;
            double maxDelay = std::max( (*pDivs0)[iL].delay, (*pDivs0)[iR].delay );
            tt = ~ttL & ~ttR;
            if( kitty::count_ones(tt)>0 && kitty::count_zeros(tt)>0 )
                divisors.emplace_back( divisors.size(), tt, 1, maxDelay+1, gate_t::AI00, fanins );
            for( int iTrg{0}; iTrg < targets.size(); ++iTrg )
            {
                if( targets[iTrg].isDone ) continue;
                if( kitty::equal( targets[iTrg].tt, divisors.back().tt ) || kitty::equal( targets[iTrg].tt, ~divisors.back().tt ) )
                {    divisors.back().isPo = true; divisors.back().id2 = iTrg; targets[iTrg].div = divisors.back().id;}
            }

            tt = ~ttL &  ttR;
            if( kitty::count_ones(tt)>0 && kitty::count_zeros(tt)>0 )
                divisors.emplace_back( divisors.size(), tt, 1, maxDelay+1, gate_t::AI01, fanins );
            for( int iTrg{0}; iTrg < targets.size(); ++iTrg )
            {
                if( targets[iTrg].isDone ) continue;
                if( kitty::equal( targets[iTrg].tt, divisors.back().tt ) || kitty::equal( targets[iTrg].tt, ~divisors.back().tt ) )
                {    divisors.back().isPo = true; divisors.back().id2 = iTrg; targets[iTrg].div = divisors.back().id; }
            }

            tt = ttL & ~ttR;
            if( kitty::count_ones(tt)>0 && kitty::count_zeros(tt)>0 )
                divisors.emplace_back( divisors.size(), tt, 1, maxDelay+1, gate_t::AI10, fanins );
            for( int iTrg{0}; iTrg < targets.size(); ++iTrg )
            {
                if( targets[iTrg].isDone ) continue;
                if( kitty::equal( targets[iTrg].tt, divisors.back().tt ) || kitty::equal( targets[iTrg].tt, ~divisors.back().tt ) )
                {    divisors.back().isPo = true; divisors.back().id2 = iTrg; targets[iTrg].div = divisors.back().id; }
            }

            tt =  ttL &  ttR;
            if( kitty::count_ones(tt)>0 && kitty::count_zeros(tt)>0 )
                divisors.emplace_back( divisors.size(), tt, 1, maxDelay+1, gate_t::AI11, fanins );
            for( int iTrg{0}; iTrg < targets.size(); ++iTrg )
            {
                if( targets[iTrg].isDone ) continue;
                if( kitty::equal( targets[iTrg].tt, divisors.back().tt ) || kitty::equal( targets[iTrg].tt, ~divisors.back().tt ) )
                {    divisors.back().isPo = true; divisors.back().id2 = iTrg; targets[iTrg].div = divisors.back().id; }
            }

            tt =  ttL ^  ttR;
            if( kitty::count_ones(tt)>0 && kitty::count_zeros(tt)>0 )
                divisors.emplace_back( divisors.size(), tt, 1, maxDelay+1, gate_t::EXOR, fanins );
            for( int iTrg{0}; iTrg < targets.size(); ++iTrg )
            {
                if( targets[iTrg].isDone ) continue;
                if( kitty::equal( targets[iTrg].tt, divisors.back().tt ) || kitty::equal( targets[iTrg].tt, ~divisors.back().tt ) )
                {    divisors.back().isPo = true; divisors.back().id2 = iTrg; targets[iTrg].div = divisors.back().id; }
            }

        }
    }  

}

support_generator_t::support_generator_t( std::vector<divisor_t> * pDivs0, std::vector<target_t> * pTrgs0, node_ps ps ): ps(ps)
{
    /* store trivial solution */
    std::vector<int> identity;
    for( int i{0}; i < pDivs0->size(); ++i ) identity.push_back( i );
    history.insert( identity );
    /*  */
    next_layer( pDivs0, pTrgs0 );
};

template<>
std::vector<int> support_generator_t::find_new<supp_selection_t::SUP_RAND>( int nIters )
{
    bool isEnd = true;
    for( int i{0}; i<targets.size(); ++i )
        isEnd &= targets[i].isDone;
    std::vector<int> support0;

    if( isEnd )
        return support0;
    for( int i{0}; i<divisors.size();++i )
    {
        if( divisors[i].isPo )
        {
            support0.push_back(i);
        }
    }
    std::vector<int> support;
    double BETA=0;

    for( int it{0}; it<nIters; ++it )
    {
        support = support0;
        std::vector<DTT> target_graphs;
        for( int i{0}; i<targets.size(); ++i )
        {
            //if( targets[i].div < 0 )
            target_graphs.push_back( targets[i].graph );
        }
        std::vector<int> divisors_id;    
        for( int i{0}; i<divisors.size(); ++i )
            divisors_id.push_back( i );

        for( int i{support.size()-1}; i>=0; --i )
        {
           target_graphs = cover_the_targets( &target_graphs, divisors[support[i]].graph ); 
           divisors_id.erase( divisors_id.begin() + i );
        }

        int iNdPar{0};
        int iNd;
        while( target_graphs.size() > 0 )
        {
            std::uniform_int_distribution<> distrib(0, divisors_id.size()-1);
            int iNew = distrib(ml_gen);
            // update the targets
            target_graphs = cover_the_targets( &target_graphs, divisors[divisors_id[iNew]].graph );
            support.push_back( divisors_id[iNew] );
            divisors_id.erase( divisors_id.begin() + iNew );

            for( int i{target_graphs.size()-1}; i>=0; --i )
            {
                if( kitty::count_ones( target_graphs[i] ) == 0 )
                    target_graphs.erase( target_graphs.begin() + i );
            }
        }
        std::sort( support.begin(), support.end() );
        if(support.size() > 1) support = erase_non_essential( &divisors, &targets, support );
        if( history.find_in_set( support ) == history.end_of_set() )
            return support;
        else
            support = {};
    }
    return support;   
}

template<>
std::vector<int> support_generator_t::find_new<supp_selection_t::SUP_ENER>( int nIters )
{
    bool isEnd = true;
    for( int i{0}; i<targets.size(); ++i )
        isEnd &= targets[i].isDone;
    std::vector<int> support0;

    if( isEnd )
        return support0;
    for( int i{0}; i<divisors.size();++i )
    {
        if( divisors[i].isPo )
        {
            support0.push_back(i);
        }
    }
    std::vector<int> support;
    double BETA;

    for( int it{0}; it<nIters; ++it )
    {
        BETA = nIters <= 1 ? ps.BETA0 : ps.BETA0 + it*(ps.BETAZ-ps.BETA0)/(nIters-1);
        support = support0;
        std::vector<DTT> target_graphs;
        for( int i{0}; i<targets.size(); ++i )
        {
            //if( targets[i].div < 0 )
            target_graphs.push_back( targets[i].graph );
        }
        std::vector<int> divisors_id;    
        for( int i{0}; i<divisors.size(); ++i )
            divisors_id.push_back( i );

        for( int i{support.size()-1}; i>=0; --i )
        {
           target_graphs = cover_the_targets( &target_graphs, divisors[support[i]].graph ); 
           divisors_id.erase( divisors_id.begin() + i );
        }

        int iNdPar{0};
        int iNd;
        while( target_graphs.size() > 0 )
        {
            std::vector<double> costs = compute_costs( ps, &divisors, &target_graphs, divisors_id );
            std::vector<double> CDF = compute_cdf( costs, BETA );
            int iNew = choose_divisor_from_cdf( CDF );
            // update the targets
            target_graphs = cover_the_targets( &target_graphs, divisors[divisors_id[iNew]].graph );
            support.push_back( divisors_id[iNew] );
            divisors_id.erase( divisors_id.begin() + iNew );

            for( int i{target_graphs.size()-1}; i>=0; --i )
            {
                if( kitty::count_ones( target_graphs[i] ) == 0 )
                    target_graphs.erase( target_graphs.begin() + i );
            }
        }
        std::sort( support.begin(), support.end() );
        if(support.size() > 1) support = erase_non_essential( &divisors, &targets, support );
        if( history.find_in_set( support ) == history.end_of_set() )
            return support;
        else
            support = {};
    }
    return support;   
}

template<>
std::vector<int> support_generator_t::find_new<supp_selection_t::SUP_GENE>( int nIters )
{
    bool isEnd = true;
    for( int i{0}; i<targets.size(); ++i )
        isEnd &= targets[i].isDone;

    std::vector<int> support0;

    if( isEnd )
        return support0;
    for( int i{0}; i<divisors.size();++i )
        if( divisors[i].isPo )
            support0.push_back(i);

    std::vector<int> support;

    for( int it{0}; it<nIters; ++it )
    {
        support = support0;
        std::vector<DTT> target_graphs;
        for( int i{0}; i<targets.size(); ++i )
            target_graphs.push_back( targets[i].graph );

        std::vector<int> divisors_id;    
        for( int i{0}; i<divisors.size(); ++i )
            divisors_id.push_back( i );

        for( int i{support.size()-1}; i>=0; --i )
        {
           target_graphs = cover_the_targets( &target_graphs, divisors[support[i]].graph ); 
           divisors_id.erase( divisors_id.begin() + i );
        }

        int iNdPar{0};
        int iNd;
        while( target_graphs.size() > 0 )
        {
            std::vector<double> costs = compute_costs( ps, &divisors, &target_graphs, divisors_id );
            std::vector<double> CDF = compute_cdf( costs, 0 );
            int iNew = choose_divisor_from_cdf( CDF );
            // update the targets
            target_graphs = cover_the_targets( &target_graphs, divisors[divisors_id[iNew]].graph );
            support.push_back( divisors_id[iNew] );
            divisors_id.erase( divisors_id.begin() + iNew );

            for( int i{target_graphs.size()-1}; i>=0; --i )
            {
                if( kitty::count_ones( target_graphs[i] ) == 0 )
                    target_graphs.erase( target_graphs.begin() + i );
            }
        }
        std::sort( support.begin(), support.end() );
        if(support.size() > 1) support = erase_non_essential( &divisors, &targets, support );
        if( history.find_in_set( support ) == history.end_of_set() )
            return support;
        else
            support = {};
    }
    return support;
    
}

void support_generator_t::mark_closing_divisors()
{
    for( auto trg : targets )
    {
        if( trg.isDone )
        int idDiv=-1;
        for( int iDiv{0}; iDiv < divisors.size(); ++iDiv )
        {
            if( divisors[iDiv].isPo != 0 ) continue;
            if( kitty::equal(divisors[iDiv].tt, trg.tt) )
            {
                divisors[iDiv].isPo = 1;
            }
            else if( kitty::equal(divisors[iDiv].tt,~trg.tt) )
                divisors[iDiv].isPo =-1;
            else
                divisors[iDiv].isPo = 0;
        }
    }
}

void support_generator_t::store_new( std::vector<int> support )
{
    std::sort( support.begin(), support.end() );
    history.insert( support );
}

void support_generator_t::print()
{
    printf("DIVISORS\n");
    for( auto div : divisors )
        div.print();
    printf("TARGETS\n");
    for( auto trg : targets )
        trg.print();
    
}

void support_generator_t::add_cost( int idSupp, double cost )
{
    if( history.costs[idSupp] > cost )  history.costs[idSupp] = cost;
}

} // namespace mcts

} // namespace mockturtle
