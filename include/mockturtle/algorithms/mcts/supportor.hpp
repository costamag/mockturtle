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


std::vector<double> bdd_compute_costs( std::vector<divisor_t> * pDivs, std::vector<DTT> * pFns, std::vector<DTT> * pMks, std::vector<int> ids )
{
    double edge_to_cover{0};
    for( int j{0}; j<pFns->size(); ++j )
    {
        double n = (double)kitty::count_ones( (*pMks)[j] );
        double n1 = (double)kitty::count_ones( (*pFns)[j]&(*pMks)[j] );
        double n0 = n-n1;
        edge_to_cover += n0*n1;
    }

    double min_cost = std::numeric_limits<double>::max();
    double max_cost = std::numeric_limits<double>::min();

    std::vector<double> costs;
    double nBits = (*pFns)[0].num_bits();
    for( int i{0}; i<ids.size(); ++i )
    {
        double cost{0};
        for( int j{0}; j<pFns->size(); ++j )
        {
            double nOn = kitty::count_ones( (*pMks)[j] );
            if( nOn > 0 )
            {
                DTT D = (*pDivs)[ids[i]].tt;
                DTT M0 = (*pMks)[j] & ~D;
                DTT M1 = (*pMks)[j] &  D;
                int nM0 = kitty::count_ones(M0);
                int c1M0 = kitty::count_ones(M0 & (*pFns)[j] );
                int c0M0 = nM0 - c1M0;

                int nM1 = kitty::count_ones(M1);
                int c1M1 = kitty::count_ones(M1 & (*pFns)[j] );
                int c0M1 = nM1 - c1M1;

                cost += c0M0*c1M0+c0M1*c1M1;
            }
        }
        costs.push_back(cost/edge_to_cover);

        min_cost = min_cost < costs[i] ? min_cost : costs[i];
        max_cost = max_cost > costs[i] ? max_cost : costs[i];
        //printf("cost %f\n", 1-max_reward/pFns->size());
    }

    for( int i{0}; i<costs.size(); ++i )
        costs[i] = (costs[i]-min_cost)/(max_cost-min_cost);

    return costs;
}

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
        divisors.emplace_back( ps.use_inf_graph, divisors.size(), (*pDivs0)[i].tt, 0, (*pDivs0)[i].delay, gate_t::PRJL, fanins );
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

            for( int i{0}; i<ps.lib.size(); ++i )
            {
                if( ps.lib[i].nInputs == 2 )
                {
                    tt = ps.lib[i].compute( {ttR, ttL} );
                    if( kitty::count_ones(tt)>0 && kitty::count_zeros(tt)>0 )
                        divisors.emplace_back( ps.use_inf_graph, divisors.size(), tt, 1, maxDelay+ps.lib[i].delay, ps.lib[i].type, fanins );
                    for( int iTrg{0}; iTrg < targets.size(); ++iTrg )
                    {
                        if( targets[iTrg].isDone ) continue;
                        if( kitty::equal( targets[iTrg].tt, divisors.back().tt ) || kitty::equal( targets[iTrg].tt, ~divisors.back().tt ) )
                        {    divisors.back().isPo = true; divisors.back().id2 = iTrg; targets[iTrg].div = divisors.back().id; }
                    }
                }
                else if( ps.lib[i].nInputs == 3 )
                {
                    for( int iZ{iL+1}; iZ < pDivs0->size(); ++iZ )
                    {
                        fanins = { (*pDivs0)[iR].id, (*pDivs0)[iL].id, (*pDivs0)[iZ].id };
                        DTT ttZ = (*pDivs0)[iZ].tt;
                        tt = ps.lib[i].compute( {ttR, ttL, ttZ} );
                        if( kitty::count_ones(tt)>0 && kitty::count_zeros(tt)>0 )
                            divisors.emplace_back( ps.use_inf_graph, divisors.size(), tt, 1, maxDelay+ps.lib[i].delay, ps.lib[i].type, fanins );
                        for( int iTrg{0}; iTrg < targets.size(); ++iTrg )
                        {
                            if( targets[iTrg].isDone ) continue;
                            if( kitty::equal( targets[iTrg].tt, divisors.back().tt ) || kitty::equal( targets[iTrg].tt, ~divisors.back().tt ) )
                            {    divisors.back().isPo = true; divisors.back().id2 = iTrg; targets[iTrg].div = divisors.back().id; }
                        }  
                    }  
                }
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
    for( int i{0}; i<divisors.size();++i )
    {
        if( divisors[i].isPo )
        {
            support0.push_back(i);
        }
    }
    if( isEnd )
        return support0;
        
    std::vector<int> support;
    double BETA;

    for( int it{0}; it<nIters; ++it )
    {
        BETA = nIters <= 1 ? ps.BETA0 : ps.BETA0 + it*(ps.BETAZ-ps.BETA0)/(nIters-1);
        printf("BETA %f\n", BETA);
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
            if( support.size() > ps.thresh )
            {
                support={};
                return support;
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
std::vector<int> support_generator_t::find_new<supp_selection_t::SUP_NORM>( int nIters )
{
    bool isEnd = true;
    for( int i{0}; i<targets.size(); ++i )
        isEnd &= targets[i].isDone;
    std::vector<int> support0;
    for( int i{0}; i<divisors.size();++i )
    {
        if( divisors[i].isPo )
        {
            support0.push_back(i);
        }
    }
    if( isEnd )
        return support0;
        
    std::vector<int> support;
    double BETA;

    for( int it{0}; it<nIters; ++it )
    {
        BETA = nIters <= 1 ? ps.BETA0 : ps.BETA0 + it*(ps.BETAZ-ps.BETA0)/(nIters-1);
        //printf("BETA %f\n", BETA);
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
            if( support.size() > ps.thresh )
            {
                support={};
                return support;
            }
        }
        std::sort( support.begin(), support.end() );
        if( ps.erase_not_essentials )
            if(support.size() > 1) support = erase_non_essential( &divisors, &targets, support );
        if( history.find_in_set( support ) == history.end_of_set() )
            return support;
        else
            support = {};
    }
    return support;   
}

template<>
std::vector<int> support_generator_t::find_new<supp_selection_t::SUP_BDD>( int nIters )
{
    bool isEnd = true;
    for( int i{0}; i<targets.size(); ++i )
        isEnd &= targets[i].isDone;
    std::vector<int> support0;
    for( int i{0}; i<divisors.size();++i )
    {
        if( divisors[i].isPo )
        {
            support0.push_back(i);
        }
    }
    if( isEnd )
        return support0;
        
    std::vector<int> support;
    double BETA;

    for( int it{0}; it<nIters; ++it )
    {
        BETA = nIters <= 1 ? ps.BETA0 : ps.BETA0 + it*(ps.BETAZ-ps.BETA0)/(nIters-1);
        //printf("BETA %f\n", BETA);
        support = support0;
        
        std::vector<DTT> target_grafs;
        std::vector<DTT> target_masks;

        for( int i{0}; i<targets.size(); ++i )
        {
            //if( targets[i].div < 0 )
            target_grafs.push_back( targets[i].tt );
            target_masks.push_back( ~(targets[i].tt.construct()) );
        }
        std::vector<int> divisors_id;    
        for( int i{0}; i<divisors.size(); ++i )
            divisors_id.push_back( i );

        for( int i{support.size()-1}; i>=0; --i )
        {
           cover_the_targets_bdd( &target_grafs, &target_masks, divisors[support[i]].tt ); 
           divisors_id.erase( divisors_id.begin() + i );
        }
        for( int i{target_grafs.size()-1}; i>=0; --i )
        {
            bool bErase{false};
            bErase |= kitty::count_ones(target_grafs[i]&target_masks[i])==0 ;
            bErase |= kitty::equal(target_grafs[i]&target_masks[i], target_masks[i] );
            if( bErase )
            {
                target_grafs.erase( target_grafs.begin() + i );
                target_masks.erase( target_masks.begin() + i );
            }
        }

        int iNdPar{0};
        int iNd;
        while( target_grafs.size() > 0 )
        {
            std::vector<double> costs = bdd_compute_costs( &divisors, &target_grafs, &target_masks, divisors_id );

            std::vector<double> CDF = compute_cdf( costs, BETA );
            
            int iNew = choose_divisor_from_cdf( CDF );

            //int knt{0};
            //for( int iD{0}; iD<divisors_id.size(); ++iD )
            //{
            //    printf("%d:", knt++);
            //    kitty::print_binary(divisors[divisors_id[iD]].tt);
            //    printf("\n");
            //}

            //printf("CHOOSE: %d\n", iNew);
            // update the targets
            cover_the_targets_bdd( &target_grafs, &target_masks, divisors[divisors_id[iNew]].tt );
            support.push_back( divisors_id[iNew] );
            divisors_id.erase( divisors_id.begin() + iNew );
            for( int i{target_grafs.size()-1}; i>=0; --i )
            {
                bool bErase{false};
                bErase |= (kitty::count_ones(target_grafs[i]&target_masks[i])==0) ;
                bErase |= kitty::equal(target_grafs[i]&target_masks[i], target_masks[i] );
                if( bErase )
                {
                    target_grafs.erase( target_grafs.begin() + i );
                    target_masks.erase( target_masks.begin() + i );
                }
            }
            //printf("|MK|=%d\n", target_masks.size() );
            //for( int i{0}; i<target_masks.size(); ++i )
            //{
            //    auto mk = target_masks[i];
            //    auto tt = target_grafs[i];
            //    printf("mk:");
            //    kitty::print_binary(mk);
            //    printf("\n");
            //    printf("tt:");
            //    kitty::print_binary(tt);
            //    printf("\n");
            //}

            if( support.size() > ps.thresh )
            {
                //printf("THR\n");
                support={};
                return support;
            }
        }
        std::sort( support.begin(), support.end() );
        //if( ps.erase_not_essentials )
        //    if(support.size() > 1) support = erase_non_essential( &divisors, &targets, support );
        if( history.find_in_set( support ) == history.end_of_set() )
            return support;
        else
            support = {};
    }
    return support;    
}

void update_targets( std::vector<DTT> * pFns, std::vector<DTT> * pMks, divisor_t div )
{
    int nFns = pFns->size();
    for( int iFn{0}; iFn < nFns; ++iFn )
    {
        DTT F0 = (*pFns)[iFn] & ~div.tt;
        DTT F1 = (*pFns)[iFn] &  div.tt;
        DTT M0 = (*pMks)[iFn] & ~div.tt; 
        DTT M1 = (*pMks)[iFn] &  div.tt; 
        pFns->push_back( F0 );
        pFns->push_back( F1 );
        pMks->push_back( M0 );
        pMks->push_back( M1 );
    }
    for( int i{nFns-1}; i>=0; i-- )
    {
        pFns->erase( pFns->begin()+i );
        pMks->erase( pMks->begin()+i );
    }
}

std::vector<double> compute_costs_( std::vector<divisor_t> * pDivs, std::vector<DTT> * pFns, std::vector<DTT> * pMks, std::vector<int> ids )
{
    std::vector<double> costs;
    double nBits = (*pFns)[0].num_bits();
    for( int i{0}; i<ids.size(); ++i )
    {
        double cost{0};
        for( int j{0}; j<pFns->size(); ++j )
        {
            double nOn = kitty::count_ones( (*pMks)[j] );
            if( nOn > 0 )
            {
                //sum += nOn;

                double T00 = kitty::count_ones( (*pMks)[j] & (~(*pFns)[j] & ~(*pDivs)[ids[i]].tt ) );
                double T01 = kitty::count_ones( (*pMks)[j] & (~(*pFns)[j] &  (*pDivs)[ids[i]].tt ) );
                double T10 = kitty::count_ones( (*pMks)[j] & ( (*pFns)[j] & ~(*pDivs)[ids[i]].tt ) );
                double T11 = kitty::count_ones( (*pMks)[j] & ( (*pFns)[j] &  (*pDivs)[ids[i]].tt ) );
                cost += 2*( T00*T01 + T11*T10 );///sum;

            }
        }
        costs.push_back(cost/pFns->size());
        //printf("cost %f\n", 1-max_reward/pFns->size());
    }
    return costs;
}



template<>
std::vector<int> support_generator_t::find_new<supp_selection_t::SUP_DECT>( int nIters )
{
    std::vector<int> support0;
    for( int i{0}; i<divisors.size();++i )
    {
        if( divisors[i].isPo )
        {
            support0.push_back(i);
        }
    }
    bool isEnd = true;
    for( int i{0}; i<targets.size(); ++i )
        isEnd &= targets[i].isDone;
    if( isEnd )
        return support0;

    std::vector<int> support;
    double BETA;
    for( int it{0}; it<nIters; ++it )
    {
        BETA = nIters <= 1 ? ps.BETA0 : ps.BETA0 + it*(ps.BETAZ-ps.BETA0)/(nIters-1);
        support = support0;
        std::vector<DTT> target_graphs;
        std::vector<DTT> target_masks;
        for( int i{0}; i<targets.size(); ++i )
        {
            target_graphs.push_back( targets[i].tt );
            target_masks.push_back( ~targets[i].tt.construct() );
        }
        std::vector<int> divisors_id;    
        for( int i{0}; i<divisors.size(); ++i )
            divisors_id.push_back( i );

        for( int i{support.size()-1}; i>=0; --i )
        {
            int nGraphs = target_graphs.size();
            for( int j{0}; j<nGraphs; ++j )
            {
                target_masks[j] &= divisors[support[i]].tt;
                target_graphs.push_back(target_graphs[j]);
                target_masks.push_back(target_masks[j] & ~divisors[support[i]].tt);
            }
            //printf("%d -> erase init ", target_graphs.size());
            for( int j{0}; j<nGraphs; ++j )
            {
                target_graphs.erase( target_graphs.begin() );
                target_masks.erase( target_masks.begin() );
            }
            
            //printf("%d\n", target_graphs.size());
            divisors_id.erase( divisors_id.begin() + support[i] );
        }
        int iNdPar{0};
        int iNd;
        while( target_graphs.size() > 0 )
        {
            printf("%d\n", support.size());
            std::vector<double> costs = compute_costs_( &divisors, &target_graphs, &target_masks, divisors_id );
            std::vector<double> CDF = compute_cdf( costs, BETA );
            for( int iCost{0}; iCost<costs.size(); ++iCost )
            {
                printf("|%f|", CDF[iCost]);
            }
            printf("\n");

            int iNew = choose_divisor_from_cdf( CDF );
            // update the targets
            int nGraphs = target_graphs.size();
            for( int i{0}; i<nGraphs; ++i )
            {
                target_graphs.push_back(target_graphs[i]);
                target_masks.push_back(target_masks[i] & ~divisors[divisors_id[iNew]].tt);
                target_masks[i] &= divisors[divisors_id[iNew]].tt;

                //printf("%d\n",i);
                //kitty::print_binary( target_masks[i] );
                //printf("\n");
                //kitty::print_binary( target_masks.back() );
                //printf("\n");
            }
            support.push_back( divisors_id[iNew] );
            divisors_id.erase( divisors_id.begin() + iNew );

            for( int i{target_graphs.size()-1}; i>=0; --i )
            {
                bool cond0 = kitty::count_ones( target_masks[i] ) == 0;
                bool cond1 = kitty::equal( target_masks[i], target_masks[i]&target_graphs[i] );
                bool cond2 = kitty::count_ones( target_masks[i]&target_graphs[i] ) == 0;
                if( cond0 || cond1 || cond2 )
                {
                    target_graphs.erase( target_graphs.begin() + i );
                    target_masks.erase( target_masks.begin() + i );
                }
            }
        }
        std::sort( support.begin(), support.end() );
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
                divisors[iDiv].isPo = 1;
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
