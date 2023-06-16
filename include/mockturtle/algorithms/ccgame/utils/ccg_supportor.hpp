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

#include "ccg_rng.hpp"

#include <kitty/print.hpp>
#include <stdio.h>
#include <stack>
#include <set>
#include <math.h>

namespace mockturtle
{

namespace ccgame
{

using DTT = kitty::dynamic_truth_table;

#pragma region INFORMATION GRAPH
/*! \brief converts a truth table to a graph representation */
DTT create_information_graph( DTT tt )
{
    int nBits = tt.num_bits();
    int nVars = tt.num_vars();
    TT  graph( 2*nVars );
    TT  tt2( 2*nVars );
    TT  mk2( 2*nVars );
    for( int iBit{0}; iBit < nBits; ++iBit )
    {
        kitty::set_bit( mk2, iBit );
        if( kitty::get_bit( tt, iBit ) == 1 )
            kitty::set_bit( tt2, iBit );
        else
            kitty::clear_bit( tt2, iBit );
    }

    for( int iBit{nBits-1}; iBit >= 0; --iBit )
    {
        if( kitty::get_bit( tt, iBit ) == 0 )
            graph |= ( tt2 << nBits*iBit );
        else
            graph |= ( ( tt2 ^ mk2 ) << nBits*iBit );
    }
    return graph;
}
#pragma endregion INFORMATION GRAPH

#pragma region DIVISOR
struct divisor_t
{
    int     id;
    DTT     tt;
    DTT     graph;
    double   area;
    double   delay;

    divisor_t(){}

    divisor_t( int id, DTT tt, double area, double delay ):
        id(id), tt(tt), area(area), delay(delay)
    {
        graph = create_information_graph( tt );
    }

    ~divisor_t(){}

    void print()
    {
        printf("[div] id:%3d area:%3.2f delay:%3.2f\n", id, area, delay ); 
        kitty::print_binary(tt);    printf("\n");   kitty::print_binary(graph); printf("\n");
    }
};
#pragma endregion DIVISOR

enum method_t
{
    BASE
};

#pragma region TARGET
struct target_t
{
    int id;
    DTT tt;
    DTT graph;
    
    target_t( int id, DTT tt ) : id(id), tt(tt) 
    {
        graph = create_information_graph( tt );
    };
    target_t(){};
    ~target_t(){};

    void print()
    {
        printf("[trg] id:%3d \n", id ); 
        kitty::print_binary(tt);    printf("\n");   kitty::print_binary(graph); printf("\n");
    }
};
#pragma endregion TARGET

class support_generator_t
{
    public:
        std::vector<divisor_t> divisors;
        std::vector<target_t> targets;
        method_t method;
        int nIdentity;
        std::set<std::vector<int>> history;

        /* CONSTRUCT */
        support_generator_t( std::vector<divisor_t>, std::vector<target_t>, method_t, int );
        support_generator_t(){};
        ~support_generator_t(){};

        /* GROW */
        std::vector<int> find_new( int );
        void store_new( std::vector<int> );
        
        /* VISUALIZE */
        void print_costs_graph();
};

std::vector<double> compute_costs( method_t method, std::vector<divisor_t> * pDivs, std::vector<DTT> * pTrgs, std::vector<int> idDivs )
{
    std::vector<double> costs;

    switch ( method )
    {
        case method_t::BASE:
        {
            for( int i{0}; i < idDivs.size(); ++i )
            {
                DTT Gi = (*pDivs)[idDivs[i]].graph;
                costs.push_back(0);
                for( int j{0}; j<pTrgs->size(); ++j )
                {
                    DTT Gf = (*pTrgs)[j];
                    costs[i] += (double)kitty::count_ones( Gf & ~Gi )/( (double)(kitty::count_ones( Gf )*pTrgs->size() ) );
                }
            }
            break;
        }    
        default:
            break;
    }
    return costs;
}

support_generator_t::support_generator_t( std::vector<divisor_t> divisors, std::vector<target_t> targets, method_t method, int nIdentity ):
divisors(divisors), targets(targets), method(method), nIdentity(nIdentity)
{
    std::vector<int> identity;
    for( int i{0}; i < nIdentity; ++i ) identity.push_back( i );
    history.insert( identity );
    for( int i{nIdentity}; i < divisors.size(); ++i ) identity.push_back( i );

    std::vector<DTT> targets_tt;
    for( auto f : targets )
        targets_tt.push_back( f.graph );
};

std::vector<double> compute_cdf( std::vector<double> H, double B )
{
    std::vector<double> P;
    std::vector<double> CDF;
    double Z=0;
    for( int i{0}; i<H.size(); ++i )
    {
        P.push_back( exp( -B*H[i] ) );
        Z += P[i];
    }
    CDF.push_back(P[0]/Z);
    for( int i{1}; i<H.size(); ++i ) CDF.push_back(CDF[i-1]+P[i]/Z);
    return CDF;
}

int choose_divisor_from_cdf( std::vector<double> CDF )
{
    std::uniform_real_distribution<> distrib(0, 1);
    double rnd = distrib(ccg_gen);
    int res;
    for( int i{0}; i<CDF.size(); ++i )
        if( rnd <= CDF[i] )
        {
            res = i;
            break;
        }
    return res;
}

std::vector<DTT> cover_the_targets( std::vector<DTT> * pGfs, DTT Gx )
{
    std::vector<DTT> res = *pGfs;
    for( int i{0}; i<pGfs->size(); ++i )
        res[i]=res[i] & ~Gx;
    return res;
}

std::vector<int> erase_non_essential( std::vector<divisor_t> * pDivs, std::vector<target_t> * pTrgs, std::vector<int> support )
{
    DTT Gx = (*pDivs)[0].graph.construct();
    DTT Gf = Gx.construct();

    for( int i{0}; i<pTrgs->size(); ++i )
        Gf |= (*pTrgs)[i].graph;

    bool isRed{false};
    while( !isRed )
    {
        std::vector<DTT> Gs;
        for( int i{0}; i<support.size(); ++i )
            Gs.push_back( (*pDivs)[support[i]].graph & Gf );
        DTT G1, G2; 
        if( support.size() > 2 )
        {
            for( int n{support.size()-1}; n >= 2; --n )
            {
                G1 = Gs[n] | Gs[n-1];
                G2 = Gs[n-2] | ( Gs[n] & Gs[n-1] );
                Gs[n-1] = G1;
                Gs[n-2] = G2;
            }
        }
        G1 = Gs[0] ^ Gs[1];

        std::vector<int> candidate_to_erase;

        isRed = true;
        for( int i{support.size()-1}; i >= 0; --i )
        {
            if( kitty::count_ones( G1 & (*pDivs)[support[i]].graph ) == 0 )
            {
                isRed = false;
                candidate_to_erase.push_back(i);
            }
        }
        if( !isRed )
        {
            std::uniform_int_distribution<> distrib(0, candidate_to_erase.size()-1);
            int to_erase = distrib(ccg_gen);
            support.erase( support.begin() + candidate_to_erase[to_erase] );
        }
    }  
    return support;
    
}

std::vector<int> support_generator_t::find_new( int nIters )
{
    std::vector<int> support;

    for( int it{0}; it<nIters; ++it )
    {
        support = {};
        std::vector<DTT> target_graphs;
        for( int i{0}; i<targets.size(); ++i )
            target_graphs.push_back( targets[i].graph );

        std::vector<int> divisors_id;    
        for( int i{0}; i<divisors.size(); ++i )
            divisors_id.push_back( i );

        int iNdPar{0};
        int iNd;

        while( target_graphs.size() > 0 )
        {
            std::vector<double> costs = compute_costs( method, &divisors, &target_graphs, divisors_id );
            std::vector<double> CDF = compute_cdf( costs, 0.00001 );
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
        if(support.size() > 0) support = erase_non_essential( &divisors, &targets, support );
        if( history.find( support ) == history.end() )
            return support;
        else
            support = {};
    }
    return support;
    
}

void support_generator_t::store_new( std::vector<int> support )
{
    std::sort( support.begin(), support.end() );
    history.insert( support );
}


} // namespace ccgame

} // namespace mockturtle
