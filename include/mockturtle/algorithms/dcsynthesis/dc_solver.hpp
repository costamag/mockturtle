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
  \file boolean_chain.hpp
  \brief data structure for storing boolean function representations

  \author Andrea Costamagna
*/
#pragma once
#include "../../networks/xag.hpp"
#include <kitty/partial_truth_table.hpp>
#include <kitty/bit_operations.hpp>
#include <kitty/operators.hpp>
#include <kitty/print.hpp>
#include <stdio.h>
#include <stack>

namespace mockturtle
{
using TT = kitty::partial_truth_table;

enum dc_funcs_t
{
    dcA11_,
    dcA01_,
    dcA10_,
    dcA00_,
    dcXOR_,
    dcBUF_,
    dcNOT_
};

struct dc_divisor_info_t
{
    int xj;
    int xi;
    dc_funcs_t fntype;
    int cost;

    dc_divisor_info_t( int xj, int xi, dc_funcs_t fntype, int cost ):
    xj(xj),
    xi(xi),
    fntype(fntype),
    cost(cost)
    {}
};

struct dc_divisors_t
{
    std::vector<TT> funcs;
    std::vector<dc_divisor_info_t> infos;
};

struct covering_table_t
{
    std::vector<TT> sets_in;
    TT              elemnts;
};

template<class Ntk>
class dc_solver
{

private:
    std::vector<TT> X0_;
    std::vector<TT> Y0_;

public:
  /* creation and destruction */
  dc_solver( std::vector<TT> const&, std::vector<TT> const& );
  ~dc_solver();
  /* information graph covering */
  TT create_information_graph( TT );
  covering_table_t create_covering_table( std::vector<TT>, std::vector<TT> );
  dc_divisors_t init_divisors( std::vector<TT> );
  dc_divisors_t create_candidate_divisors( dc_divisors_t );
  void enumerate_subsets( TT *, std::vector<std::vector<int>> *, std::vector<TT> *, std::vector<TT> *, int );
  std::vector<std::vector<int>> greedy_set_covering( covering_table_t *, dc_divisors_t * );

  void erase_invalid_subsets( std::vector<std::vector<int>> *, std::vector<TT> * );
  std::vector<std::pair<int,int>> compute_subsets_cost( dc_divisors_t *, std::vector<std::vector<int>> * );
  void select_dc_maximizers( dc_divisors_t *, std::vector<std::vector<int>> * );
  void solve_greedy( Ntk * );
  void solve_greedy_multioutput( Ntk * );

  /* visualization */
  void show_specs();
  void show_specs( std::vector<TT> );
  void show_table( covering_table_t );

};

/* creation and destruction */
template<class Ntk>
dc_solver<Ntk>::dc_solver( std::vector<TT> const& X, std::vector<TT> const& Y ): X0_(X), Y0_(Y) {}
template<class Ntk>
dc_solver<Ntk>::~dc_solver(){}

/* information graph covering */
#pragma region covering

template<class Ntk>
TT dc_solver<Ntk>::create_information_graph( TT x )
{
    int nbits = x.num_bits();
    TT igraph( nbits * nbits );
    TT xlarge( nbits * nbits );
    TT mlarge( nbits * nbits );
    assert( kitty::is_const0( igraph ) );
    assert( kitty::is_const0( xlarge ) );
    for( int b{0}; b < nbits; ++b )
    {
        kitty::set_bit( mlarge, b );
        if( kitty::get_bit( x, b ) == 1 )
            kitty::set_bit( xlarge, b );
        else
            kitty::clear_bit( xlarge, b );
    }

    for( int b{nbits-1}; b >= 0; --b )
    {
        if( kitty::get_bit( x, b ) == 0 )
            igraph |= ( xlarge << nbits*b );
        else
            igraph |= ( ( xlarge ^ mlarge ) << nbits*b );
    }
    return igraph;
}

template<class Ntk>
covering_table_t dc_solver<Ntk>::create_covering_table( std::vector<TT> X, std::vector<TT> Y )
{
    covering_table_t table;

    for( int i{0}; i < X.size(); ++i )
        table.sets_in.push_back( create_information_graph( X[i] ) );
    table.elemnts = create_information_graph( Y[0] );
    for( int i{1}; i < Y.size(); ++i )
        table.elemnts |= create_information_graph( Y[i] );
    return table;
}

template<class Ntk>
dc_divisors_t dc_solver<Ntk>::init_divisors( std::vector<TT> X )
{
    dc_divisors_t divisors;
    for( int i{0}; i < X.size() ; ++i )
    {
        divisors.funcs.push_back( X[i] );
        divisors.infos.emplace_back( i, i, dc_funcs_t::dcBUF_, 0 );
    }
    return divisors;
}

template<class Ntk>
dc_divisors_t dc_solver<Ntk>::create_candidate_divisors( dc_divisors_t X )
{
    dc_divisors_t divisors;

    for( int i{0}; i < X.funcs.size() ; ++i )
    {
        divisors.funcs.push_back( X.funcs[i] );
        divisors.infos.emplace_back( i, i, dc_funcs_t::dcBUF_, X.infos[i].cost );
    }

    for( int i{0}; i < X.funcs.size()-1 ; ++i )
    {
        for( int j{i}; j < X.funcs.size(); ++j )
        {
            divisors.funcs.push_back( X.funcs[j] & X.funcs[i] );
            divisors.infos.emplace_back( j, i, dc_funcs_t::dcA11_, 1 + X.infos[i].cost + X.infos[j].cost ); 
            
            divisors.funcs.push_back( ~X.funcs[j] & X.funcs[i] );
            divisors.infos.emplace_back( j, i, dc_funcs_t::dcA01_, 1 + X.infos[i].cost + X.infos[j].cost ); 

            divisors.funcs.push_back( X.funcs[j] & ~X.funcs[i] );
            divisors.infos.emplace_back( j, i, dc_funcs_t::dcA10_, 1 + X.infos[i].cost + X.infos[j].cost ); 

            divisors.funcs.push_back( ~X.funcs[j] & ~X.funcs[i] );
            divisors.infos.emplace_back( j, i, dc_funcs_t::dcA00_, 1 + X.infos[i].cost + X.infos[j].cost ); 

            divisors.funcs.push_back( X.funcs[j] ^ X.funcs[i] );
            divisors.infos.emplace_back( j, i, dc_funcs_t::dcXOR_, 1 + X.infos[i].cost + X.infos[j].cost ); 
        }
    }
    return divisors;
}

template<class Ntk>
void dc_solver<Ntk>::enumerate_subsets( TT * pUni0, std::vector<std::vector<int>> * pSubs, std::vector<TT> * pCGs, std::vector<TT> * pUniv, int idx )
{
    if( idx == 0 )
    {
        if( kitty::count_ones( *pUni0 & (*pCGs)[0] ) > 0 )
        {
            *pSubs = {{}, {0}};
            *pUniv = {*pUni0, *pUni0 & ~(*pCGs)[0] };
        }
        else
        {
            *pSubs = {{}};
            *pUniv = {*pUni0};
        }
    }
    else
    {
        enumerate_subsets( pUni0, pSubs, pCGs, pUniv, idx-1 );
        int nsubs = pSubs->size();
        for( int i{0}; i < nsubs; ++i )
        {
            if( kitty::count_ones( (*pUniv)[i] & (*pCGs)[idx] ) > 0 )
            {
                pSubs->push_back( (*pSubs)[i] );
                (*pSubs)[pSubs->size()-1].push_back(idx);
                pUniv->push_back( (*pUniv)[i] & ( ~(*pCGs)[idx] ) );
            }
        }
    }
}

struct dc_problems_t
{
    TT U;
    std::vector<int> avbs;
    std::vector<int> divs;

    dc_problems_t( TT U, std::vector<int> avbs, std::vector<int> divs ): U(U), avbs(avbs), divs(divs){}
};

/* this routine returns a set of subsets for th eapproximate solution of the set covering problem */
template<class Ntk>
std::vector<std::vector<int>> dc_solver<Ntk>::greedy_set_covering( covering_table_t * pTable, dc_divisors_t * pDivs )
{

    std::random_device rd;  // a seed source for the random number engine
    std::mt19937 gen(5); // mersenne_twister_engine seeded with rd()

    std::vector<std::vector<int>> res;
    std::vector<int> divs0;
    for( int i{0}; i<pDivs->funcs.size(); ++i )
        divs0.push_back(i);
    std::vector<dc_problems_t> problems_old;
    std::vector<dc_problems_t> problems_new;
    dc_problems_t problem0( pTable->elemnts, divs0, {} );
    problems_new.push_back( problem0 );
    int n_left = kitty::count_ones( problem0.U ); // number of elements to be covered
    int min_cost = n_left;
    int cost;

    while( n_left > 0 )
    {
        //printf("%d\n", problems_old.size() );
        problems_old = problems_new;
        problems_new = {};
        for( int iPb{0}; iPb < problems_old.size(); ++iPb ) // consider all problems
        {
            dc_problems_t Pb = problems_old[iPb];
            for( int iDv{0}; iDv < Pb.avbs.size(); ++iDv )
            {
                int           Dv = Pb.avbs[iDv];
                TT            S  = pTable->sets_in[Dv];
                cost = kitty::count_ones( Pb.U & ~S );
                if( cost < min_cost )
                {
                    min_cost = cost;
                    problems_new = { Pb };
                    problems_new[0].U &= ~S;
                    if( problems_new[0].divs.size() == 0 || Dv >= problems_new[0].divs[problems_new[0].divs.size()-1] )
                        problems_new[0].divs.push_back( Pb.avbs[iDv] );
                    else
                    {
                        for( int j{0}; j < problems_new[0].divs.size(); ++j )
                        {
                            if( Dv < problems_new[0].divs[j] )
                            {
                                problems_new[0].divs.insert( problems_new[0].divs.begin() + j, Pb.avbs[iDv] );
                                break;
                            }
                        }
                    }
                    problems_new[0].avbs.erase( problems_new[0].avbs.begin() + iDv );
                }
                else if ( cost == min_cost )
                {
                    int idPb = problems_new.size();
                    problems_new.push_back( Pb );
                    problems_new[idPb].U &= ~S;
                    if( problems_new[idPb].divs.size() == 0 || Dv >= problems_new[idPb].divs[problems_new[idPb].divs.size()-1] )
                        problems_new[idPb].divs.push_back( Pb.avbs[iDv] );
                    else
                    {
                        for( int j{0}; j < problems_new[idPb].divs.size(); ++j )
                        {
                            if( Dv < problems_new[idPb].divs[j] )
                            {
                                problems_new[idPb].divs.insert( problems_new[idPb].divs.begin() + j, Pb.avbs[iDv] );
                                break;
                            }
                        }
                    }
                    problems_new[idPb].avbs.erase( problems_new[idPb].avbs.begin() + iDv );
                }
            }
        }
       n_left = min_cost;
       
       if( problems_new.size() > 5 )
       {
            std::shuffle(problems_new.begin(), problems_new.end(), gen);
            for( int j{problems_new.size()-1}; j>=5; --j )
                problems_new.erase( problems_new.begin() + j );
       }

       for( int i{problems_new.size()-2}; i >= 0; --i )
       {
        for( int j{problems_new.size()-1}; j > i; --j )
        {
            if( problems_new[i].divs == problems_new[j].divs )
                problems_new.erase( problems_new.begin() + j );
        }
       }
       //printf("%d\n", n_left);
    }

    for( auto Pb : problems_new )
        res.push_back( Pb.divs );
    return res;
}

template<class Ntk>
void dc_solver<Ntk>::erase_invalid_subsets( std::vector<std::vector<int>> * pSubs, std::vector<TT> * pUnis )
{
    for( int i{pUnis->size()-1}; i >= 0; --i )
    {
        if( kitty::count_ones( (*pUnis)[i] ) > 0 )
        {
            pUnis->erase( pUnis->begin() + i );
            pSubs->erase( pSubs->begin() + i );
        }
    }
}

template<class Ntk>
std::vector<std::pair<int,int>> dc_solver<Ntk>::compute_subsets_cost( dc_divisors_t * pDivs, std::vector<std::vector<int>> * pSubs )
{
    std::vector<std::pair<int,int>> res;

    for( int i{0}; i < pSubs->size(); ++i )
    {
        int gate_cost = 0;
        for( int j{0}; j < (*pSubs)[i].size(); ++j )
            gate_cost += (*pDivs).infos[ (*pSubs)[i][j] ].cost;
        TT ref = (*pDivs).funcs[0] | ~(*pDivs).funcs[0];
        int bit = 0;
        int nDc = 0;
        while( kitty::count_ones( ref ) > 0 )
        {
            TT ttmp = (*pDivs).funcs[0] | ~(*pDivs).funcs[0];
            if( kitty::get_bit( ref, bit ) == 1 )
            {
                for( int j{0}; j < (*pSubs)[i].size(); ++j )
                {
                    if( kitty::get_bit( (*pDivs).funcs[(*pSubs)[i][j]], bit ) == 1 )
                        ttmp &= (*pDivs).funcs[(*pSubs)[i][j]];
                    else
                        ttmp &= ~(*pDivs).funcs[(*pSubs)[i][j]];
                }
                nDc += (kitty::count_ones(ttmp)-1);
                ref &= ~ttmp;
            }
            
            bit++;
        }
        res.push_back(std::make_pair(nDc, gate_cost ) );
    }
    return res;
}

template<class Ntk>
void dc_solver<Ntk>::select_dc_maximizers( dc_divisors_t * pCandidates, std::vector<std::vector<int>> * pSets )
{
    std::vector<std::pair<int, int>> costs = compute_subsets_cost( pCandidates, pSets );
    int max_dc = 0;
    for( int i{pSets->size()-1}; i>=0; --i )
    {
        if( costs[i].first < max_dc )
        {
            pSets->erase( pSets->begin() + i );
            costs.erase( costs.begin() + i );
        }
        else if( costs[i].first > max_dc )
            max_dc = costs[i].first;
    }

    for( int i{pSets->size()-1}; i>=0; --i )
        if( costs[i].first < max_dc )
        {
            pSets->erase( pSets->begin() + i );
            costs.erase( costs.begin() + i );
        }

    /*for( int i{0}; i < pSets->size(); ++i )
    {
        for( int j{0}; j< (*pSets)[i].size(); ++j )
            printf( "%d ", (*pSets)[i][j] );
        printf(" - %d %d\n", costs[i].first, costs[i].second );
    }*/
}

template<class Ntk>
void dc_solver<Ntk>::solve_greedy( Ntk * pntk )
{
    std::random_device rd;  // a seed source for the random number engine
    std::mt19937 gen(5); // mersenne_twister_engine seeded with rd()

    Ntk ntk_best;
    int size_best = 10000;
    for( int k{0}; k<100; ++k )
    {
        Ntk ntk;
        std::vector<signal<Ntk>> xold;
        std::vector<signal<Ntk>> xnew;
    //show_specs();
        std::vector<TT> Xt = X0_;
        std::vector<TT> Yt = Y0_;
        for( int i{0}; i<Xt.size(); ++i )
        {
            xold.push_back( ntk.create_pi() );
        }

        dc_divisors_t divisors = init_divisors( Xt );
        dc_divisors_t candites = create_candidate_divisors( divisors );

        bool SAT = false;

        while( !SAT )
        {
            //for( int i{0}; i<Xt.size(); ++i )
            //{
            //    printf("%d %d\n", i, xold[i].index);
           //}
            //printf("candidates\n");
/*
            for( uint32_t k{0}; k<candites.funcs.size(); ++k )
            {
                kitty::print_binary( candites.funcs[k] );
                if( candites.infos[k].fntype == dcA00_ )
                    printf(" %d and(!%d,!%d)\n", k ,candites.infos[k].xj, candites.infos[k].xi );
                else if( candites.infos[k].fntype == dcA01_ )
                    printf(" %d and(!%d, %d)\n", k, candites.infos[k].xj, candites.infos[k].xi );
                else if( candites.infos[k].fntype == dcA10_ )
                    printf(" %d and( %d,!%d)\n", k, candites.infos[k].xj, candites.infos[k].xi );
                else if( candites.infos[k].fntype == dcA11_ )
                    printf(" %d and( %d, %d)\n", k, candites.infos[k].xj, candites.infos[k].xi );
                else if( candites.infos[k].fntype == dcXOR_ )
                    printf(" %d xor( %d, %d)\n", k, candites.infos[k].xj, candites.infos[k].xi );
                else if( candites.infos[k].fntype == dcBUF_ )
                    printf(" %d buf(%d)\n", k, candites.infos[k].xi );
                else if( candites.infos[k].fntype == dcNOT_ )
                    printf(" %d not(%d)\n", k, candites.infos[k].xi );
            }
*/
            xnew = {};
            covering_table_t table = create_covering_table( candites.funcs, Yt );
    //show_table( table );
            std::vector<std::vector<int>> subsets = greedy_set_covering( &table, &candites );
            //select_dc_maximizers( &candites, &subsets );
            
            std::uniform_int_distribution<> distrib(0, subsets.size()-1);
            int rnum = distrib(gen);
            std::vector<int> Ssel = subsets[rnum];

            Xt = {};
            divisors = {};
            assert( divisors.funcs.size() == 0 );
            for( int s{0}; s < Ssel.size(); ++s )
            {
                Xt.push_back( candites.funcs[Ssel[s]] );
                divisors.funcs.push_back( candites.funcs[Ssel[s]] );
                divisors.infos.push_back( candites.infos[Ssel[s]] );
                int xi = candites.infos[Ssel[s]].xi;
                int xj = candites.infos[Ssel[s]].xj;

                if( candites.infos[Ssel[s]].fntype == dc_funcs_t::dcA00_ )
                {
                    xnew.push_back( ntk.create_and( ntk.create_not(xold[xj]), ntk.create_not(xold[xi]) ) );
                    //printf("and(%d', %d')\n", xold[xj].index, xold[xi].index );
                }
                else if( candites.infos[Ssel[s]].fntype == dc_funcs_t::dcA01_ )
                {
                    xnew.push_back( ntk.create_and( ntk.create_not( xold[xj] ), xold[xi] ) );
                    //printf("and(%d', %d)\n", xold[xj].index, xold[xi].index );
                }
                else if( candites.infos[Ssel[s]].fntype == dc_funcs_t::dcA10_ )
                {
                    xnew.push_back( ntk.create_and( xold[xj], ntk.create_not( xold[xi] ) ) );
                    //printf("and(%d, %d')\n", xold[xj].index, xold[xi].index );
                }
                else if( candites.infos[Ssel[s]].fntype == dc_funcs_t::dcA11_ )
                {
                    xnew.push_back( ntk.create_and( xold[xj], xold[xi] ) );
                    //printf("and(%d, %d)\n", xold[xj].index, xold[xi].index );
                }
                else if( candites.infos[Ssel[s]].fntype == dc_funcs_t::dcXOR_ )
                {
                    xnew.push_back( ntk.create_xor( xold[xj], xold[xi] ) );
                    //printf("xor(%d, %d)\n", xold[xj].index, xold[xi].index );
                }
                else if( candites.infos[Ssel[s]].fntype == dc_funcs_t::dcBUF_ )
                {
                    xnew.push_back( xold[xi] );
                    //printf("%d\n", xold[xi].index );
                }
                else if( candites.infos[Ssel[s]].fntype == dc_funcs_t::dcNOT_ )
                {
                    xnew.push_back( ntk.create_not(xold[xi]) );
                    //printf("%d'\n", xold[xi].index );
                }
            }

            if( Ssel.size() == 1 )
            {
                SAT = true;
                if( kitty::equal( candites.funcs[Ssel[0]], Y0_[0] ) )
                    ntk.create_po( xnew[0] );
                else if( kitty::equal( ~candites.funcs[Ssel[0]], Y0_[0] ) )
                    ntk.create_po( ntk.create_not( xnew[0] ) );
                else
                    assert(0);

            }
            candites = create_candidate_divisors( divisors );
            //show_specs( Xt );

            xold = xnew;
        }
        printf("ngates = %d\n", ntk.num_gates());
        if( ntk.num_gates() < size_best )
        {
            size_best = ntk.num_gates();
            *pntk     = ntk;
        }
    }
}

template<class Ntk>
void dc_solver<Ntk>::solve_greedy_multioutput( Ntk * pntk )
{
    std::random_device rd;  // a seed source for the random number engine
    std::mt19937 gen(5); // mersenne_twister_engine seeded with rd()

    for( int k{0}; k<1; ++k )
    {
        Ntk ntk;
        std::vector<signal<Ntk>> xold;
        std::vector<signal<Ntk>> xnew;
        std::vector<TT> Xt = X0_;
        std::vector<TT> Yt = Y0_;
        for( int i{0}; i<Xt.size(); ++i )
            xold.push_back( ntk.create_pi() );

        dc_divisors_t divisors = init_divisors( Xt );
        dc_divisors_t candites = create_candidate_divisors( divisors );

        bool SAT = false;
        std::vector<signal<Ntk>> out_signals;
        std::vector<int> out_ids;
        for( uint32_t iOut{0}; iOut < Yt.size() ; ++iOut )
        {
            out_signals.emplace_back();
            out_ids.push_back( iOut );
        }
        std::vector<uint32_t> closure;
        std::vector<bool>     inversi;
        uint32_t mask = 0x0FFFFFFF;
        while( Yt.size() > 0 )
        {
            xnew = {};
            covering_table_t table = create_covering_table( candites.funcs, Yt );
            std::vector<std::vector<int>> subsets = greedy_set_covering( &table, &candites );
            //select_dc_maximizers( &candites, &subsets );
            std::uniform_int_distribution<> distrib(0, subsets.size()-1);
            int rnum = distrib(gen);
            std::vector<int> Ssel = subsets[rnum];
            Xt = {};
            divisors = {};
            assert( divisors.funcs.size() == 0 );

            for( int s{0}; s < Ssel.size(); ++s )
            {
                Xt.push_back( candites.funcs[Ssel[s]] );
                divisors.funcs.push_back( candites.funcs[Ssel[s]] );
                divisors.infos.push_back( candites.infos[Ssel[s]] );
                int xi = candites.infos[Ssel[s]].xi;
                int xj = candites.infos[Ssel[s]].xj;

                if( candites.infos[Ssel[s]].fntype == dcA00_ )
                    xnew.push_back( ntk.create_and( ntk.create_not(xold[xj]), ntk.create_not(xold[xi]) ) );
                else if( candites.infos[Ssel[s]].fntype == dcA01_ )
                    xnew.push_back( ntk.create_and( ntk.create_not( xold[xj] ), xold[xi] ) );
                else if( candites.infos[Ssel[s]].fntype == dcA10_ )
                    xnew.push_back( ntk.create_and( xold[xj], ntk.create_not( xold[xi] ) ) );
                else if( candites.infos[Ssel[s]].fntype == dcA11_ )
                    xnew.push_back( ntk.create_and( xold[xj], xold[xi] ) );
                else if( candites.infos[Ssel[s]].fntype == dcXOR_ )
                    xnew.push_back( ntk.create_xor( xold[xj], xold[xi] ) );
                else if( candites.infos[Ssel[s]].fntype == dcBUF_ )
                    xnew.push_back( xold[xi] );
                else if( candites.infos[Ssel[s]].fntype == dcNOT_ )
                    xnew.push_back( ntk.create_not(xold[xi]) );
            }

            for( int iOut{out_ids.size()-1}; iOut >= 0 ; --iOut )
            {
                closure = {};
                inversi = {};
                for( uint32_t iDiv{0}; iDiv < divisors.funcs.size(); ++iDiv)
                {
                    if( kitty::equal( divisors.funcs[iDiv], ~Yt[iOut] ) )
                    {
                        closure.push_back( iDiv );
                        inversi.push_back( true );
                    }
                    else if( kitty::equal( divisors.funcs[iDiv], Yt[iOut] ) )
                    {
                        closure.push_back( iDiv );
                        inversi.push_back( false );
                    }
                }
                if( closure.size() > 0 )
                {
                    std::uniform_int_distribution<> distrib(0, closure.size()-1);
                    int rnum = distrib(gen);
                    uint32_t iClDiv = closure[rnum];
                    if( inversi[rnum] )
                        out_signals[out_ids[iOut]] = ntk.create_not( xnew[iClDiv] );
                    else
                        out_signals[out_ids[iOut]] = xnew[iClDiv];

                    Yt.erase( Yt.begin() + iOut );
                    out_ids.erase( out_ids.begin() + iOut );
                }
            }
            candites = create_candidate_divisors( divisors );
            xold = xnew;
            printf("NUMGATES=%d NUM OUT=%d\n", ntk.num_gates(), Yt.size());
        }
        
        for( int iOut{0}; iOut < out_signals.size() ; ++iOut )
            ntk.create_po(out_signals[iOut]);
        printf("ngates = %d\n", ntk.num_gates());
        *pntk = ntk;
    }
}

#pragma endregion covering

/* visualization */
#pragma region visualization
template<class Ntk>
void dc_solver<Ntk>::show_specs()
{
    printf( "      " );
    for( int i{X0_.size()-1}; i >= 0 ; --i )
        printf( "%d ", i );
    printf("| ");
    for( int i{0}; i < Y0_.size(); ++i )
        printf( "%d ", i );
    printf("\n");
    
    for( int i{0}; i < Y0_.size()+X0_.size()+6; ++i )
        printf( "==", i );
    printf( "\n" );

    for( int b{0}; b < X0_[0].num_bits(); ++b )
    {
        printf("%4d: ", b );
        for( int i{X0_.size()-1}; i >= 0 ; --i )
            printf( "%d ", kitty::get_bit( X0_[i], b ) );
        printf("| ");
        for( int i{0}; i < Y0_.size(); ++i )
            printf( "%d ", kitty::get_bit( Y0_[i], b ) );
        printf("\n");
    }
}

template<class Ntk>
void dc_solver<Ntk>::show_specs( std::vector<TT> X )
{
    printf( "      " );
    for( int i{X.size()-1}; i >= 0 ; --i )
        printf( "%d ", i );
    printf("| ");
    for( int i{0}; i < Y0_.size(); ++i )
        printf( "%d ", i );
    printf("\n");
    
    for( int i{0}; i < Y0_.size()+X.size()+6; ++i )
        printf( "==", i );
    printf( "\n" );

    for( int b{0}; b < X[0].num_bits(); ++b )
    {
        printf("%4d: ", b );
        for( int i{X.size()-1}; i >= 0 ; --i )
            printf( "%d ", kitty::get_bit( X[i], b ) );
        printf("| ");
        for( int i{0}; i < Y0_.size(); ++i )
            printf( "%d ", kitty::get_bit( Y0_[i], b ) );
        printf("\n");
    }
}

template<class Ntk>
void dc_solver<Ntk>::show_table( covering_table_t table )
{
    int nbits = sqrt( table.elemnts.num_bits() );
    printf( "\n          " );
    for( int i{0}; i < table.sets_in.size() ; ++i )
        printf( "%d ", i );
    printf("| ");
    printf("Y\n");
    for( int i{0}; i <table.sets_in.size()+7; ++i )
        printf( "==" );
    printf("\n");

    for( int b{0}; b < table.sets_in[0].num_bits(); ++b )
    {
        printf( "%3d %3d : ", b/nbits, b%nbits );
        for( int i{0}; i < table.sets_in.size(); ++i )
            printf("%d ", kitty::get_bit( table.sets_in[i], b ) );
        printf("| %d\n", kitty::get_bit( table.elemnts, b ) );
    }
}
#pragma endregion visualization
} // namespace mockturtle