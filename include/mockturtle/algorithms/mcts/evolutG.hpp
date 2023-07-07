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
  \file mnist_reader.hpp
  \brief handler of the mnist data

  \author Andrea Costamagna
*/
#pragma once

#include <kitty/partial_truth_table.hpp>
#include <kitty/bit_operations.hpp>
#include <kitty/constructors.hpp>
#include <kitty/print.hpp>
#include "ml_rng.hpp"
#include "genet.hpp"
#include <stdio.h>
#include <stack>
#include <iostream>
#include <string>
#include <fstream>
#include <bit>
#include <bitset>
#include <cstdint>
#include <fmt/format.h>
#include <set>

namespace mockturtle
{

namespace mcts
{

struct evolutG_ps_t
{
    double P0 = 0.5;
    double PZ = 0;
    double frac = 0.5;
    int nInd = 20;
    int nCandElim = 5;
    int nGens = 100;
    double mutation_rate = 0.0001;
};

class evolutG_t
{
public:
    genet_t gen0;
    evolutG_ps_t ps;
    std::vector<genet_t> population;
    std::vector<double> rewards;
    std::vector<double> flipProbs;
    genet_t bestInd;
    int nNodes;

public:
    evolutG_t( genet_t, evolutG_ps_t );
    ~evolutG_t();

    void train();
    void create_generation0();
    std::pair<int, int> binary_tournament_selection();
    void create_new_children( std::pair<int, int> );
    void simulated_annealing();
    void simulate_system( int, double );

};

evolutG_t::evolutG_t( genet_t Gen0, evolutG_ps_t Ps ) : gen0(Gen0), bestInd(Gen0), ps(Ps)
{
    int nLyrs = Gen0.net.size();
    if( nLyrs <= 2 )
        flipProbs = {0.0};
    else
        for( int iLyr{1}; iLyr<nLyrs; ++iLyr )
            flipProbs.push_back( ( ps.P0 + (double)(iLyr-1)*(ps.PZ-ps.P0)/(nLyrs-2) ) );
    
    nNodes=0;
    for( int iLyr{1}; iLyr<nLyrs; ++iLyr )
        for( int iNd{0}; iNd<Gen0.net[iLyr].size()-1; ++iNd  )
            nNodes++;
}

evolutG_t::~evolutG_t()
{
}

void evolutG_t::create_generation0()
{
    int nLyrs = gen0.net.size();
    std::uniform_int_distribution<> distrib_layers(1, nLyrs-2 );
    std::uniform_int_distribution<> distrib_seed( 0, 1000 );
    std::uniform_real_distribution<> distrib_thr( 0, 1 );

    for( int iInd{0}; iInd < ps.nInd; ++iInd )
    {
        genet_t geni = gen0;

        for( int iRnd{0}; iRnd < (int)(ps.frac*nNodes); ++iRnd )
        {
            int iLyr = distrib_layers( ml_gen );
            std::uniform_int_distribution<> distrib_nodes( 0, geni.net[iLyr].size()-1 );
            int iNd = distrib_nodes( ml_gen );
            if( distrib_thr(ml_gen) < flipProbs[iLyr-1] )
            {
                PTT rand = geni.net[iLyr][iNd].y_train.construct();
                kitty::create_random( rand, distrib_seed(ml_gen) );
                geni.net[iLyr][iNd].y_train = geni.net[iLyr][iNd].y_train ^ (  rand );//geni.e_train &
            }
        }
        geni.train_network();
        printf("Atr=%f Ava=%f Ate=%f\n", geni.acc_train(), geni.acc_valid(), geni.acc_test() );
        rewards.push_back( geni.acc_valid() );
        population.push_back( geni );
    }
}

std::pair<int, int> evolutG_t::binary_tournament_selection()
{
    std::vector<int> group0;
    std::vector<int> group1;
    double bestRwd0 = -1;
    double bestRwd1 = -1;
    int best0;
    int best1;
    std::uniform_int_distribution<> distrib(0, 1);

    for( int i{0}; i<rewards.size(); ++i )
    {

        if(distrib(ml_gen)==1)
        {
            group1.push_back(i);
            if( rewards[i] > bestRwd1 )
            {
                best1 = i;
                bestRwd1 = rewards[i];
            }
        }
        else
        {
            group0.push_back(i);
            if( rewards[i] > bestRwd0 )
            {
                best0 = i;
                bestRwd0 = rewards[i];
            }
        }
    }
    printf("GROUP 0: ");
    for( auto x : group0 )
        printf("[%d %f] ", x, rewards[x] );
    printf("\n");
    printf("GROUP 1: ");
    for( auto x : group1 )
        printf("[%d %f] ", x, rewards[x] );
    printf("\n");
    if( group0.size() == 0 )
    {
        std::uniform_int_distribution<> distrib2(0, group1.size()-1);
        group0.push_back( group1[distrib2(ml_gen)]);
        best0 = group0[0];

    }
    if( group1.size() == 0 )
    {
        std::uniform_int_distribution<> distrib2(0, group0.size()-1);
        group1.push_back( group0[distrib2(ml_gen)]);
        best1 = group1[0];
    }

    return std::make_pair( best0, best1 );
}

void evolutG_t::create_new_children( std::pair<int, int> parents )
{
    printf("CROSSOVER WITH %d %d\n", parents.first, parents.second );
    genet_t P1 = population[parents.first];
    genet_t P2 = population[parents.second];
    double  f1 = P1.acc_train();
    double  f2 = P2.acc_train();
    genet_t C  = P1;
    std::uniform_real_distribution<> distrib(0, 1);
    int nLyrs = gen0.net.size();
    
    for( int iLyr{1}; iLyr<nLyrs-1; ++iLyr )
    {
        for( int iNd{0}; iNd<gen0.net[iLyr].size(); ++iNd )
        {
            PTT diff = P1.net[iLyr][iNd].y_train^P2.net[iLyr][iNd].y_train;
            if( distrib(ml_gen) < f2/(f1+f2) )
                C.net[iLyr][iNd].y_train = P2.net[iLyr][iNd].y_train;
        }
    }
    std::uniform_int_distribution<> distrib_layers(1, nLyrs-2 );
    std::uniform_real_distribution<> distrib_seed(0, 1000);
    
    
    int nToMutate = ps.mutation_rate * nNodes;
    for( int i{0}; i < nToMutate; ++i )
    {
        int iLyr = distrib_layers( ml_gen );
        std::uniform_int_distribution<> distrib_nodes( 0, C.net[iLyr].size()-1 );
        int iNd = distrib_nodes( ml_gen );
        PTT rand = C.net[iLyr][iNd].y_train.construct();
        kitty::create_random( rand, distrib_seed(ml_gen) );
        C.net[iLyr][iNd].y_train = C.net[iLyr][iNd].y_train ^ ( C.e_train & rand );
    }



    C.train_network();

    if( C.acc_valid() > bestInd.acc_valid() )
        bestInd = C;
    printf("Atr=%f Ava=%f Ate=%f\n", C.acc_train(), C.acc_valid(), C.acc_test() );

    std::vector<double> rewards_tmp = rewards;
    while( rewards_tmp.size() > ps.nCandElim )
    {
        int idx;
        double bestRwd = -1;
        for( int i{0}; i<rewards_tmp.size(); ++i )
        {
            if( rewards_tmp[i]>bestRwd )
            {
                bestRwd = rewards_tmp[i];
                idx=i;
            }
        }
        rewards_tmp.erase( rewards_tmp.begin() + idx );
    }
    std::uniform_int_distribution<> distrib2(0, rewards_tmp.size()-1);
    double reward_to_remove = distrib2( ml_gen );
    
    int idx;
    for( int i{0}; i<rewards.size(); ++i )
    {
        idx = i;
        if( rewards[i] == reward_to_remove )    break;
    }
    if( C.acc_valid() > rewards[idx] )
    {
        rewards[idx] = C.acc_valid();
        population[idx] = C;
    }
}

void evolutG_t::train()
{
    create_generation0();
    
    for( int iGen{0}; iGen < ps.nGens; ++iGen )
    {
        std::pair<int, int> parents = binary_tournament_selection();
        create_new_children( parents );
    }
}

void evolutG_t::simulate_system( int iIn, double Beta )
{
    genet_t net0 = population[iIn];
    genet_t net1 = net0;
    int nLyrs = gen0.net.size();
    std::uniform_int_distribution<> distrib_layers(1, nLyrs-2 );
    std::uniform_int_distribution<> distrib_nodes( 0, gen0.net[0].size()-1 );
    std::uniform_int_distribution<> distrib_seed( 0, nLyrs*1000 );
    std::uniform_real_distribution<> distrib_mc( 0, 1 );

    std::vector<double> energies;
    for( int it{0}; it < 10; ++it )
    {
        net1 = net0;
        for( int iMod{0}; iMod<10; iMod++ )
        {
            int iLyr = distrib_layers( ml_gen );
            int iNd = distrib_nodes( ml_gen ) % gen0.net[iLyr].size();
            PTT rand = net1.net[iLyr][iNd].y_train.construct();
            kitty::create_random( rand, distrib_seed(ml_gen) );
            net1.net[iLyr][iNd].y_train = net1.net[iLyr][iNd].y_train ^ ( net1.e_train & rand) ;
        }
        net1.train_network();

        double dE = net1.acc_valid()-net0.acc_valid();

        printf( "%f %f %f\n", dE, Beta*dE, exp( -Beta*dE ) );

        if( distrib_mc(ml_gen) < exp( -Beta*dE ) )
        {
            net0 = net1;
            printf("accepted\n");
        }        
        else
            printf("rejected\n");
    }  
    population[iIn] = net1;
    if( net1.acc_valid() > bestInd.acc_valid() )
        bestInd = net1;
}

void evolutG_t::simulated_annealing()
{
    create_generation0();
    double Beta = 0.00001;
    for( int iT{0}; iT < 10; ++iT )
    {
        for( int iIn{0}; iIn<population.size(); ++iIn )
        {
            printf("INIT %d\n", iIn );
            simulate_system( iIn, Beta );
            Beta*=10;
        }
    }
}

} // namespace mcts

} // namespace mockturtle