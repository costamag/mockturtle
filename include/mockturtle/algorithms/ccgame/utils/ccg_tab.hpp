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
#include <kitty/partial_truth_table.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/bit_operations.hpp>
#include <kitty/operators.hpp>
#include <kitty/print.hpp>
#include <stdio.h>

#include <stack>
#include "ccg_cut.hpp"
#include "ccg_rng.hpp"


namespace mockturtle
{

namespace ccgame
{

using TT = kitty::dynamic_truth_table;

struct problems_t
{
    TT U;
    std::vector<int> avbs;
    std::vector<int> divs;

    problems_t( TT U, std::vector<int> avbs, std::vector<int> divs ): U(U), avbs(avbs), divs(divs){}
};

class tab_t
{

public:
  cut_t inCut;
  cut_t outCut;
  std::vector<TT> sets;
  TT univ;
  std::vector<std::vector<int>> subsets;

  tab_t();
  tab_t( cut_t, cut_t );
  ~tab_t();

  void init_tab( cut_t, cut_t );
  void init_small_tab( cut_t, cut_t );
  void print();
  void greedy_set_covering( int );
  void implicit_greedy_set_covering( );
  void select_dc_maximizers();
  std::vector<int> compute_subsets_cost();
};

#pragma region constructors
tab_t::tab_t( cut_t sets_c, cut_t univ_c ): inCut( sets_c ), outCut( univ_c )
{
  for( int i{0}; i < sets_c.size(); ++i )
      sets.push_back( sets_c.nodes[i].graph() );

  univ = sets[0] & ~sets[0];

  for( int i{0}; i < univ_c.size(); ++i )
  {
    if( univ_c.nodes[i].gate == gate_t::POS )
      univ |= univ_c.nodes[i].graph();
  }
}
tab_t::tab_t(){}
tab_t::~tab_t(){}

void tab_t::init_tab( cut_t sets_c, cut_t univ_c )
{
  inCut = sets_c;
  outCut = univ_c;
  for( int i{0}; i < sets_c.size(); ++i )
      sets.push_back( sets_c.nodes[i].graph() );

  univ = sets[0] & ~sets[0];

  for( int i{0}; i < univ_c.size(); ++i )
  {
    if( univ_c.nodes[i].gate == gate_t::POS )
      univ |= univ_c.nodes[i].graph();
  }
}

void tab_t::init_small_tab( cut_t sets_c, cut_t univ_c )
{
  inCut = sets_c;
  outCut = univ_c;
  for( int i{0}; i < sets_c.size(); ++i )
      sets.push_back( sets_c.nodes[i].tt );

  univ = sets[0] & ~sets[0];

  for( int i{0}; i < univ_c.size(); ++i )
  {
    if( univ_c.nodes[i].gate == gate_t::POS )
      univ |= univ_c.nodes[i].tt;
  }
}

#pragma endregion

#pragma region set-covering
/*! \brief approximate solution of the set covering problem */
void tab_t::greedy_set_covering( int nCap )
{
    subsets = {};
    std::random_device rd;  // a seed source for the random number engine
    std::mt19937 gen(5); // mersenne_twister_engine seeded with rd()

    std::vector<int> divs0;
    for( int i{0}; i< sets.size(); ++i )
        divs0.push_back(i);
    std::vector<problems_t> problems_old;
    std::vector<problems_t> problems_new;
    problems_t problem0( univ, divs0, {} );
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
            problems_t Pb = problems_old[iPb];
            for( int iDv{0}; iDv < Pb.avbs.size(); ++iDv )
            {
                int           Dv = Pb.avbs[iDv];
                TT            S  = sets[Dv];
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
       
       if( nCap > 0 && problems_new.size() > nCap )
       {
            std::shuffle(problems_new.begin(), problems_new.end(), gen);
            for( int j{problems_new.size()-1}; j>=nCap; --j )
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
        subsets.push_back( Pb.divs );
}

/*! \brief approximate solution of the set covering problem */
void tab_t::implicit_greedy_set_covering()
{
    std::uniform_int_distribution<> distrib(0, sets.size()-1);
    subsets = {};

    std::random_device rd;  // a seed source for the random number engine
    std::mt19937 gen(5); // mersenne_twister_engine seeded with rd()

    bool notFound = true;
    std::vector<int> subset;

    while( notFound )
    {
        std::vector<TT> Fn = {univ};
        std::vector<TT> Mk = {~univ.construct() | univ.construct()};
        std::vector<int> to_use;
        for( int i{0}; i<sets.size(); ++i )
            to_use.push_back(i);
        int rnum = distrib(ccg_gen);
        subset = {};
        int isFirst=1;
        while( Fn.size() > 0 )
        {
            if( to_use.size() == 0 )
                break;
            int bestRef=-1;
            /* choose the next divisor */
            if( isFirst == 0 )
            {
                int bestCorr = 0;
                int corr;

                for( int iRef{0}; iRef < to_use.size(); ++iRef )
                {
                    corr = 0;
                    for( int iFn{0}; iFn < Fn.size(); ++iFn )
                        corr += std::max( kitty::count_ones( Mk[iFn]&(Fn[iFn]^sets[to_use[iRef]]) ), kitty::count_ones( Mk[iFn]&(~Fn[iFn]^sets[to_use[iRef]])));
                    if( corr > bestCorr )
                    {
                        bestCorr = corr;
                        bestRef = iRef;
                    }
                }
            }
            else
            {
                bestRef = rnum;
                isFirst = 0;
            }
            if( bestRef < 0 )
                break;

            /* update the remainders */
            std::vector<int> to_eliminateL;
            std::vector<int> to_eliminateR;
            int nFuncs = Fn.size();
            for( int iFn{0}; iFn < nFuncs; ++iFn )
            {
                Fn.push_back(Fn[iFn]);
                Mk.push_back(Mk[iFn]&sets[to_use[bestRef]]);
                Mk[iFn] &= ~sets[to_use[bestRef]];
                if( kitty::count_ones(Mk[iFn]) == 0 || kitty::count_ones(Mk[iFn]&Fn[iFn]) == 0 || kitty::equal(Mk[iFn]&Fn[iFn], Mk[iFn]) )
                    to_eliminateL.push_back( iFn );
                if( kitty::count_ones(Mk[nFuncs+iFn]) == 0 || kitty::count_ones(Mk[nFuncs+iFn]&Fn[nFuncs+iFn]) == 0 || kitty::equal(Mk[nFuncs+iFn]&Fn[nFuncs+iFn], Mk[nFuncs+iFn]) )
                    to_eliminateR.push_back( nFuncs+iFn );
            }
            if( bestRef < 0  )
                break;
            subset.push_back( to_use[bestRef] );
            to_use.erase( to_use.begin() + bestRef );
            for( int iRem{to_eliminateR.size()-1}; iRem >= 0; --iRem )
                Fn.erase( Fn.begin() + to_eliminateR[iRem] );
            for( int iRem{to_eliminateL.size()-1}; iRem >= 0; --iRem )
                Fn.erase( Fn.begin() + to_eliminateL[iRem] );

            if( Fn.size() == 0 )
                notFound = false;
        }
    }
    subsets = { subset };
    printf("subsets.size() %d  | subset.size() %d\n", subsets.size(), subsets[0].size());
}

std::vector<int> tab_t::compute_subsets_cost()
{
    std::vector<int> res;

    for( int i{0}; i < subsets.size(); ++i )
    {
        TT ref = inCut.nodes[0].tt | ~inCut.nodes[0].tt;
        int bit = 0;
        int nDc = 0;
        while( kitty::count_ones( ref ) > 0 )
        {
            TT ttmp = inCut.nodes[0].tt | ~inCut.nodes[0].tt;
            if( kitty::get_bit( ref, bit ) == 1 )
            {
                for( int j{0}; j < subsets[i].size(); ++j )
                {
                    if( kitty::get_bit( inCut.nodes[subsets[i][j]].tt, bit ) == 1 )
                        ttmp &=  inCut.nodes[subsets[i][j]].tt;
                    else
                        ttmp &= ~inCut.nodes[subsets[i][j]].tt;
                }
                nDc += (kitty::count_ones(ttmp)-1);
                ref &= ~ttmp;
            }
            bit++;
        }
        res.push_back( nDc );
    }
    return res;
}

void tab_t::select_dc_maximizers()
{
    std::vector<int> costs = compute_subsets_cost();
    int max_dc = 0;
    for( int i{subsets.size()-1}; i>=0; --i )
    {
        if( costs[i] < max_dc )
        {
            subsets.erase( subsets.begin() + i );
            costs.erase( costs.begin() + i );
        }
        else if( costs[i] > max_dc )
            max_dc = costs[i];
    }

    for( int i{subsets.size()-1}; i>=0; --i )
        if( costs[i] < max_dc )
        {
            subsets.erase( subsets.begin() + i );
            costs.erase( costs.begin() + i );
        }
    
}
#pragma endregion set-covering

#pragma region visualize
void tab_t::print()
{
  int nbits = sqrt( univ.num_bits() );
  printf( "\n          " );
  for( int i{0}; i < sets.size() ; ++i )
      printf( "%d ", i );
  printf("| ");
  printf("Y\n");
  for( int i{0}; i < sets.size()+7; ++i )
      printf( "==" );
  printf("\n");

  for( int b{0}; b < sets[0].num_bits(); ++b )
  {
      printf( "%3d %3d : ", b/nbits, b%nbits );
      for( int i{0}; i < sets.size(); ++i )
          printf("%d ", kitty::get_bit( sets[i], b ) );
      printf("| %d\n", kitty::get_bit( univ, b ) );
  }
}
#pragma endregion visualize

} // namespace ccgame

} // namespace mockturtle

