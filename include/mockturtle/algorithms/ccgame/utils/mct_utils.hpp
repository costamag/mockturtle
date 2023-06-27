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
  \brief data structure for storing the cuts for the ccgame

  \author Andrea Costamagna
*/
#pragma once

#include "ccg_rng.hpp"
#include "ccg_supportor.hpp"
#include <stdio.h>
#include <stack>
#include <iostream>

#define PRINTER(name) printer(#name)

void printer(char *name) {
    printf("%s\n", name);
}

namespace mockturtle
{

namespace ccgame
{

using PTT = kitty::partial_truth_table;
using DTT = kitty::dynamic_truth_table;


#pragma region INFORMATION GRAPH
/*! \brief converts a truth table to a graph representation */
DTT create_information_graph( DTT tt )
{
    int nBits = tt.num_bits();
    int nVars = tt.num_vars();
    DTT  graph( 2*nVars );
    DTT  tt2( 2*nVars );
    DTT  mk2( 2*nVars );
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

/*! \brief Gate type in the ccgame namespace. Convention Xl=1100, Xr=1010
 */

#pragma region DIVISOR
struct divisor_t
{
    std::vector<int> fanins;
    int     id;
    int     id2;
    DTT     tt;
    DTT     graph;
    double  area;
    double  delay;
    gate_t  type;
    int     isPo{0};
 
    divisor_t(){}

    divisor_t( int id, DTT tt, double area, double delay, gate_t type ):
        id(id), tt(tt), area(area), delay(delay), type(type)
    {
        graph = create_information_graph( tt );
    }

    divisor_t( int id, DTT tt, double area, double delay, gate_t type, std::vector<int> fanins ):
        id(id), tt(tt), area(area), delay(delay), type(type), fanins(fanins)
    {
        graph = create_information_graph( tt );
    }

    divisor_t( int id, DTT tt, double area, double delay ):
        id(id), tt(tt), area(area), delay(delay)
    {
        graph = create_information_graph( tt );
    }

    ~divisor_t(){}

    void print()
    {
        if( isPo )
            printf("[%3d] id:%3d area:%3.2f delay:%3.2f ", id2, id, area, delay );
        else
            printf("[div] id:%3d area:%3.2f delay:%3.2f ", id, area, delay );
            
        switch (type)
        {
        case gate_t::AI00:      printf("AI00 : ");       break;
        case gate_t::AI01:      printf("AI01 : ");       break;
        case gate_t::AI10:      printf("AI10 : ");       break;
        case gate_t::AI11:      printf("AI11 : ");       break;
        case gate_t::CMPL:      printf("CMPL : ");       break;
        case gate_t::CMPR:      printf("CMPR : ");       break;
        case gate_t::CNTR:      printf("CNTR : ");       break;
        case gate_t::EXOR:      printf("EXOR : ");       break;
        case gate_t::OI00:      printf("OI00 : ");       break;
        case gate_t::OI01:      printf("OI01 : ");       break;
        case gate_t::OI10:      printf("OI10 : ");       break;
        case gate_t::OI11:      printf("OI11 : ");       break;
        case gate_t::PIS :      printf("PI   : ");       break;
        case gate_t::POS :      printf("PO   : ");       break;
        case gate_t::PRJL:      printf("PRJL : ");       break;
        case gate_t::PRJR:      printf("PRJR : ");       break;
        case gate_t::TAUT:      printf("TAUT : ");       break;
        case gate_t::XNOR:      printf("XNOR : ");       break;
        
        default:
            break;
        }
        for( auto fi : fanins )
            printf(" %d ", fi);
        printf("\n");
        kitty::print_binary(tt);    printf("\n");   kitty::print_binary(graph); printf("\n");
    }
};
#pragma endregion DIVISOR

enum gen_method_t
{
    BASE
};

#pragma region TARGET
struct target_t
{
    int id;
    int div;
    DTT tt;
    DTT graph;
    gate_t type;
    bool isDone{false};
    
    target_t( int id, DTT tt ) : id(id), tt(tt) 
    {
        div = -1;
        graph = create_information_graph( tt );
    };
    target_t(){};
    ~target_t(){};

    void print()
    {
        printf("[trg] id:%3d is done? %d\n", id, isDone ); 
        kitty::print_binary(tt);    printf("\n");   kitty::print_binary(graph); printf("\n");
    }
};
#pragma endregion TARGET

std::vector<double> compute_costs( gen_method_t method, std::vector<divisor_t> * pDivs, std::vector<DTT> * pTrgs, std::vector<int> idDivs )
{
    std::vector<double> costs;

    switch ( method )
    {
        case gen_method_t::BASE:
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
    int res=0;

    for( int i{0}; i<CDF.size(); ++i )
    {
        if( rnd <= CDF[i] )
        {
            res = i;
            break;
        }
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
        if( support.size() > 1 )
        {
            for( int n{support.size()-1}; n >= 2; --n )
            {
                G1 = Gs[n] | Gs[n-1];
                G2 = Gs[n-2] | ( Gs[n] & Gs[n-1] );
                Gs[n-1] = G1;
                Gs[n-2] = G2;
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
        else
            isRed = true;
    }  
    return support;   
}


} // namespace ccgame

} // namespace mockturtle