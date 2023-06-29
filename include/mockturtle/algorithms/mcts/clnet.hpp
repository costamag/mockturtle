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
#include <stdio.h>
#include <stack>
#include <iostream>
#include <string>
#include <fstream>
#include <bit>
#include <bitset>
#include <cstdint>
#include <fmt/format.h>


namespace mockturtle
{

namespace mcts
{

struct clnet_ps
{
    int nFilters{1};
    int nGen0{2};
    int d0 = 28;
    int d1 = 28;
    int d2 = 1;
    int nGenerations = 1;
};

class clnet
{
private:
    std::vector<std::vector<PTT>> XYs;
    std::vector<decision_tree> trees;
    std::vector<double> rewards;
    std::vector<PTT> abstract_inputs;
    std::vector<PTT> x_train;
    std::vector<PTT> y_train;
    std::vector<PTT> m_train;
    std::vector<PTT> x_valid;
    std::vector<PTT> y_valid;
    std::vector<PTT> m_valid;
    std::vector<PTT> x_test;
    std::vector<PTT> y_test;
    std::vector<PTT> m_test;
    PTT train0;
    PTT train1;
    PTT valid0;
    PTT valid1;
    PTT test0;
    PTT test1;
    clnet_ps ps;

public:
    clnet( std::vector<PTT>, std::vector<PTT>, std::vector<PTT>, std::vector<PTT>, std::vector<PTT>, std::vector<PTT>, clnet_ps );
    ~clnet();
    void train();
    void create_gen_0();
    void print_genome();
    double filter_evaluation( decision_tree * );
    std::vector<PTT*> get_train_zone( int ); 
    std::vector<PTT*> get_test_zone( int ); 
    std::vector<PTT*> get_valid_zone( int ); 
    std::pair<int, int> binary_tournament();
    void create_next_gen( std::pair<int, int> );
};

clnet::clnet( std::vector<PTT> xtrain, std::vector<PTT> ytrain, std::vector<PTT> xvalid, std::vector<PTT> yvalid, std::vector<PTT> xtest, std::vector<PTT> ytest, clnet_ps ps ):
    x_train(xtrain), y_train(ytrain), x_valid(xvalid), y_valid(yvalid), x_test(xtest), y_test(ytest), ps(ps)
{
    printf("INIT\n");

    for( int iTrg{0}; iTrg<y_train.size(); ++iTrg )
        m_train.push_back( ~y_train[0].construct() );
    for( int iTrg{0}; iTrg<y_test.size(); ++iTrg )
        m_test.push_back( ~y_test[0].construct() );
    for( int iTrg{0}; iTrg<y_valid.size(); ++iTrg )
        m_valid.push_back( ~y_valid[0].construct() );

    assert( x_test[0].num_bits()==y_test[0].num_bits() );
    assert( x_train[0].num_bits()==y_train[0].num_bits() );
    assert( x_valid[0].num_bits()==y_valid[0].num_bits() );
    printf("|y valid|=%d\n", y_valid[0].num_bits());
    assert( m_valid[0].num_bits()==x_valid[0].num_bits() );
    assert( m_test[0].num_bits()==x_test[0].num_bits() );
    assert( m_train[0].num_bits()==x_train[0].num_bits() );

    train0 = x_train[0].construct();
    train1 = ~x_train[0].construct();
    valid0 =  x_valid[0].construct();
    valid1 = ~x_valid[0].construct();
    test0 = x_test[0].construct();
    test1 = ~x_test[0].construct();

    for( int i{0}; i<9; ++i )
    {
        abstract_inputs.emplace_back(512u);
        kitty::create_nth_var( abstract_inputs[i], i );
    }
}

clnet::~clnet()
{
}


std::vector<PTT*> clnet::get_train_zone( int iFtr )
{
    std::vector<int> ids;
    ids.push_back( iFtr - 28 - 1 );
    ids.push_back( iFtr - 28 );
    ids.push_back( iFtr - 28 + 1 );
    ids.push_back( iFtr + 1 );
    ids.push_back( iFtr + 28 + 1 );
    ids.push_back( iFtr + 28 );
    ids.push_back( iFtr + 28 - 1 );
    ids.push_back( iFtr - 1 );
    ids.push_back( iFtr );

    if( iFtr % 28 == 0 )
    {
        ids[0] = -1;
        ids[6] = -1;
        ids[7] = -1;
    }
    if( iFtr % 28 == 27 )
    {
        ids[2] = -1;
        ids[3] = -1;
        ids[4] = -1;    
    }
    if( iFtr / 28 == 0 )
    {
        ids[0] = -1;
        ids[6] = -1;
        ids[7] = -1;
    }
    if( iFtr / 28 == 27 )
    {
        ids[4] = -1;
        ids[5] = -1;
        ids[6] = -1;    
    }
    std::vector<PTT*> zone;
    zone.push_back(&train0);
    zone.push_back(&train1);
    for( int i{0}; i<9; ++i )
    {
        //printf("%d r=%d c=%d\n", ids[i], ids[i]/28, ids[i]%28 );
        if( ids[i] < 0 )
            zone.push_back( &train0 );
        else
            zone.push_back( &x_train[ids[i]] );
    }
    return zone;
}

std::vector<PTT*> clnet::get_test_zone( int iFtr )
{
    std::vector<int> ids;
    ids.push_back( iFtr - 28 - 1 );
    ids.push_back( iFtr - 28 );
    ids.push_back( iFtr - 28 + 1 );
    ids.push_back( iFtr + 1 );
    ids.push_back( iFtr + 28 + 1 );
    ids.push_back( iFtr + 28 );
    ids.push_back( iFtr + 28 - 1 );
    ids.push_back( iFtr - 1 );
    ids.push_back( iFtr );

    if( iFtr % 28 == 0 )
    {
        ids[0] = -1;
        ids[6] = -1;
        ids[7] = -1;
    }
    if( iFtr % 28 == 27 )
    {
        ids[2] = -1;
        ids[3] = -1;
        ids[4] = -1;    
    }
    if( iFtr / 28 == 0 )
    {
        ids[0] = -1;
        ids[6] = -1;
        ids[7] = -1;
    }
    if( iFtr / 28 == 27 )
    {
        ids[4] = -1;
        ids[5] = -1;
        ids[6] = -1;    
    }
    std::vector<PTT*> zone;
    zone.push_back( &test0);
    zone.push_back( &test1);
    for( int i{0}; i<9; ++i )
    {
        if( ids[i] < 0 )
            zone.push_back( &test0 );
        else
            zone.push_back( &x_test[ids[i]] );
    }
    return zone;
}

std::vector<PTT*> clnet::get_valid_zone( int iFtr )
{
    std::vector<int> ids;
    ids.push_back( iFtr - 28 - 1 );
    ids.push_back( iFtr - 28 );
    ids.push_back( iFtr - 28 + 1 );
    ids.push_back( iFtr + 1 );
    ids.push_back( iFtr + 28 + 1 );
    ids.push_back( iFtr + 28 );
    ids.push_back( iFtr + 28 - 1 );
    ids.push_back( iFtr - 1 );
    ids.push_back( iFtr );

    if( iFtr % 28 == 0 )
    {
        ids[0] = -1;
        ids[6] = -1;
        ids[7] = -1;
    }
    if( iFtr % 28 == 27 )
    {
        ids[2] = -1;
        ids[3] = -1;
        ids[4] = -1;    
    }
    if( iFtr / 28 == 0 )
    {
        ids[0] = -1;
        ids[6] = -1;
        ids[7] = -1;
    }
    if( iFtr / 28 == 27 )
    {
        ids[4] = -1;
        ids[5] = -1;
        ids[6] = -1;    
    }
    std::vector<PTT*> zone;
    zone.push_back( &valid0);
    zone.push_back( &valid1);
    for( int i{0}; i<9; ++i )
    {
        if( ids[i] < 0 )
            zone.push_back( &valid0 );
        else
            zone.push_back( &x_valid[ids[i]] );
    }
    return zone;
}

double clnet::filter_evaluation( decision_tree * pKer )
{
    std::vector<PTT> next_layer_train;
    std::vector<PTT> next_layer_test;
    std::vector<PTT> next_layer_valid;

    for( int iFtr{0}; iFtr < x_train.size(); ++iFtr )
    {
        std::vector<PTT*> zone_train = get_train_zone( iFtr );
        std::vector<PTT> conv_train = pKer->compute( zone_train );
        for( int iKer{0}; iKer<conv_train.size(); ++iKer )
            next_layer_train.push_back(conv_train[iKer] );

        std::vector<PTT*> zone_valid = get_valid_zone( iFtr );
        std::vector<PTT> conv_valid = pKer->compute( zone_valid );
        for( int iKer{0}; iKer<conv_valid.size(); ++iKer )
            next_layer_valid.push_back(conv_valid[iKer] );

        std::vector<PTT*> zone_test = get_test_zone( iFtr );
        std::vector<PTT> conv_test = pKer->compute( zone_test );
        for( int iKer{0}; iKer<conv_test.size(); ++iKer )
            next_layer_test.push_back(conv_test[iKer] );
    }
    //printf("tr %d va %d te %d\n", next_layer_train[0].num_bits(), next_layer_valid[0].num_bits(), next_layer_test[0].num_bits());
    //printf("tr %d va %d te %d\n", y_train[0].num_bits(), y_valid[0].num_bits(), y_test[0].num_bits());
    decision_tree tree_eval( next_layer_train, y_train, next_layer_valid, y_valid, next_layer_test, y_test ); 
    tree_eval.train_impurity<entropy_t::SHAN>();
    double Atrain = tree_eval.train_accuracy();
    printf( "Atrain %f ", Atrain );
    double Avalid = tree_eval.valid_accuracy();
    printf( "Avalid %f ", Avalid );
    double Atest  = tree_eval.test_accuracy();
    printf( "Atest %f\n", Atest );
    return Avalid;
}

void clnet::create_gen_0()
{
    std::uniform_int_distribution<> distrib(0, 1000);
    for( int xy{0}; xy < ps.nGen0; ++xy )
    {
        std::vector<PTT> chromosome;
        for( int i{0}; i<ps.nFilters; ++i )
        {
            chromosome.emplace_back(512u);
            create_random(chromosome[i], distrib(ml_gen) );
        }
        decision_tree tree( abstract_inputs, chromosome, abstract_inputs, chromosome, abstract_inputs, chromosome ); 
        tree.train_impurity<entropy_t::SHAN>();
        //printf( "%f %f\n", tree.train_accuracy(), tree.test_accuracy() );
        trees.push_back( tree );
        XYs.push_back( chromosome );
        rewards.push_back( filter_evaluation( &trees[xy]) );
    }
}

std::pair<int, int> clnet::binary_tournament()
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
            if( rewards[i] > bestRwd1 )
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

    return std::make_pair( best0, best1 );

}

void clnet::create_next_gen( std::pair<int, int> parents )
{
    std::vector<PTT> tt1 = XYs[parents.first]; 
    std::vector<PTT> tt2 = XYs[parents.second]; 
    std::vector<PTT> ttn;
    PTT mk;
    std::uniform_real_distribution<> distrib(0, 1);
    for( int i{0}; i<tt1.size(); ++i )
    {
        mk = ~tt1[i]^tt2[i];
        for( int iBit{0}; iBit < mk.num_bits(); ++iBit )
        {
            if( kitty::get_bit( mk, iBit ) == 0 )
            {
                if( distrib(ml_gen) < 0.1 )
                    kitty::flip_bit( tt1[i], iBit );
                if( distrib(ml_gen) < 0.1 )
                    kitty::flip_bit( tt2[i], iBit );
            }
        }
    } 


    decision_tree tree1( abstract_inputs, tt1, abstract_inputs, tt1, abstract_inputs, tt1 ); 
    tree1.train_impurity<entropy_t::SHAN>();
    printf( "%f %f\n", tree1.train_accuracy(), tree1.test_accuracy() );
    double reward1 = filter_evaluation( &tree1 ) ;


    decision_tree tree2( abstract_inputs, tt2, abstract_inputs, tt2, abstract_inputs, tt2 ); 
    tree2.train_impurity<entropy_t::SHAN>();
    printf( "%f %f\n", tree2.train_accuracy(), tree2.test_accuracy() );
    double reward2 = filter_evaluation( &tree2 ) ;

    printf( "reward1 %f reward2 %f\n", reward1, reward2 );
    int idx;
    double minimum=100000000000;
    for( int i{0}; i<rewards.size(); ++i )
    {
        if( rewards[i]<minimum )
        {
            idx = i;
            minimum = rewards[i];
        }
    }
    if( reward1 > rewards[idx] )
    {
        XYs[idx] = tt1;
        trees[idx] = tree1;
        rewards[idx] = reward1;
    }

    minimum=100000000000;
    for( int i{0}; i<rewards.size(); ++i )
    {
        if( rewards[i]<minimum )
        {
            idx = i;
            minimum = rewards[i];
        }
    }
    if( reward2 > rewards[idx] )
    {
        XYs[idx] = tt2;
        trees[idx] = tree2;
        rewards[idx] = reward2;
    }

}

void clnet::train()
{
    create_gen_0();
    for( int it{1}; it<ps.nGenerations; ++it )
    {
        printf("GENERATION %d\n", it);
        std::pair<int, int> parents = binary_tournament();
        printf("PARENTS: %d %d\n", parents.first, parents.second );
        create_next_gen( parents );
    }
}

void clnet::print_genome()
{
    for( int i{0}; i < XYs.size(); ++i )
    {
        printf("GEN %d \n", i);
        for( int j{0}; j<XYs[i].size(); ++j )
        {
            kitty::print_binary( XYs[i][j] );
            printf("\n");
        }
    }
}

} // namespace mcts

} // namespace mockturtle