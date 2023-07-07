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
#include <set>


namespace mockturtle
{

namespace mcts
{

class genet_node_t
{
    public:
        int idNd;
        std::vector<int> inputs;
        PTT s_train;
        PTT s_valid;
        PTT s_test;
        PTT y_train;
        PTT m_train;
        PTT y_valid;
        PTT m_valid;
        PTT y_test;
        PTT m_test;
        bool isPi;
        bool isPo;
        bool isUsed;

    genet_node_t( int IdNd, PTT S_train, PTT S_valid, PTT S_test ) : 
        idNd(IdNd), s_train(S_train), s_valid(S_valid), s_test(S_test) 
        {
            isPi = true;
            isPo = false;
            isUsed = true;
        };

    genet_node_t( int IdNd, std::vector<int> Inputs, PTT Y_train, PTT M_train ) : 
        idNd(IdNd), inputs(Inputs), y_train(Y_train), m_train(M_train) 
        {
            isPi = false;
            isPo = false;
            isUsed = true;
        };

    genet_node_t( int IdNd, std::vector<int> Inputs, PTT Y_train, PTT M_train, PTT Y_valid, PTT M_valid, PTT Y_test, PTT M_test ) : 
        idNd(IdNd), inputs(Inputs), y_train(Y_train), m_train(M_train), y_valid(Y_valid), m_valid(M_valid), y_test(Y_test), m_test(M_test) 
        {
            isPi = false;
            isPo = true;
            isUsed = true;
        };

};

struct genet_data_t
{
    std::vector<PTT> x;
    PTT y;
    PTT m;
    genet_data_t( std::vector<PTT> X, PTT Y, PTT M ) : x(X), y(Y), m(M){};
};

struct genet_ps_t
{
    std::vector<int> specs;
    int K;
};

enum genet_netcreator_t
{
    CREA_RAND = 0
};

class genet_t
{
public:
    PTT e_train;
    genet_data_t train;
    genet_data_t valid;
    genet_data_t test;
    genet_ps_t   ps;
    PTT train0;
    PTT train1;
    std::vector<std::vector<genet_node_t>> net;

    double accTrain;
    double accValid;
    double accTest;
public:
    genet_t();
    genet_t( genet_data_t, genet_data_t, genet_data_t, genet_ps_t );
    ~genet_t();
    void print();
    template<genet_netcreator_t CREATOR>
    void create_network();

    void train_network();
    void train_node( int, int );

    double acc_train();
    double acc_test();
    double acc_valid();
};

genet_t::genet_t( genet_data_t Train, genet_data_t Valid, genet_data_t Test, genet_ps_t Ps ) :
    train(Train), valid(Valid), test(Test), ps(Ps)
{
    /* INPUT LAYER */
    std::vector<genet_node_t> layer0;
    for( int i{0}; i<train.x.size(); ++i )
        layer0.emplace_back( i, train.x[i], valid.x[i], test.x[i] );
    net.push_back(layer0);

    train0 = train.y.construct();
    train1 = ~train0;
}

genet_t::~genet_t()
{
}
template<>
void genet_t::create_network<genet_netcreator_t::CREA_RAND>()
{
    for( int i{0}; i<ps.specs.size(); ++i )
    {
        int iLyrPrev = i;
        std::uniform_int_distribution<> distrib(0, net[i].size()-1 );
        int iLyr = i+1;
        std::vector<genet_node_t> layer;
        for( int iNd{0}; iNd < ps.specs[i]; iNd++ )
        {
            std::vector<int> inputs;
            while( inputs.size() < ps.K )
            {
                int iNew = distrib(ml_gen);
                if( std::find( inputs.begin(), inputs.end(), iNew ) == inputs.end() )
                    inputs.push_back(iNew);
            }
            std::sort( inputs.begin(), inputs.end() );
            layer.emplace_back( iNd, inputs, train.y, train.m );
        }
        net.push_back( layer );
    }
    /* create output */
    std::uniform_int_distribution<> distrib(0, net.back().size()-1 );
    std::vector<int> inputs;
    while( inputs.size() < ps.K )
    {
        int iNew = distrib(ml_gen);
        if( std::find( inputs.begin(), inputs.end(), iNew ) == inputs.end() )
            inputs.push_back(iNew);
    }
    std::sort( inputs.begin(), inputs.end() );
    genet_node_t output( 0, inputs, train.y, train.m, valid.y, valid.m, test.y, valid.m );
    net.push_back( {output} );

    for( int iLyr{ net.size()-1 }; iLyr > 0; --iLyr )
    {
        int iLyrPrev = iLyr-1;
        int nLyrPrev = net[iLyrPrev].size();
        std::set<int> used;
        for( int iNd{0}; iNd < net[iLyr].size(); ++iNd )
        {
            if( net[iLyr][iNd].isUsed )
            {
                for( auto child: net[iLyr][iNd].inputs )
                    used.insert( child );
            }
        }
        for( int i{0}; i<net[iLyrPrev].size(); ++i )
            if( used.find( i ) == used.end() )
                net[iLyrPrev][i].isUsed = false;
    }
}

void genet_t::train_node( int iLyr, int iNd )
{
    std::uniform_real_distribution<> distrib(0, 1 );
    
    std::vector<PTT*> Xtr;
    std::vector<PTT*> Xva;
    std::vector<PTT*> Xte;
    int iLyrPrev = iLyr-1;
    for( auto idChild : net[iLyr][iNd].inputs )
    {
        Xtr.push_back( &(net[iLyrPrev][idChild].s_train) );
        Xva.push_back( &(net[iLyrPrev][idChild].s_valid) );
        Xte.push_back( &(net[iLyrPrev][idChild].s_test) );
    }
    net[iLyr][iNd].s_train = net[0][0].s_train.construct();
    net[iLyr][iNd].s_valid = net[0][0].s_valid.construct();
    net[iLyr][iNd].s_test  = net[0][0].s_test.construct();
    /* training */
    PTT todos = train1;
    for( int iBit{0}; iBit<todos.num_bits(); ++iBit )
    {
        if( kitty::count_ones( todos ) == 0 )   break;
        PTT find_train = ~train.y.construct();
        PTT find_valid = ~valid.y.construct();
        PTT find_test  = ~test.y.construct();
        if( kitty::get_bit( todos, iBit ) == 1 )
        {
            for( int iX{0}; iX<net[iLyr][iNd].inputs.size(); ++iX )
            {
                if( kitty::get_bit( *(Xtr[iX]), iBit ) == 1 )
                {
                    find_train &= *(Xtr[iX]);
                    find_valid &= *(Xva[iX]);
                    find_test  &= *(Xte[iX]);
                }
                else
                {
                    find_train &= ~*(Xtr[iX]);
                    find_valid &= ~*(Xva[iX]);
                    find_test  &= ~*(Xte[iX]);
                }
            }
            PTT tt1 = find_train & net[iLyr][iNd].y_train;
            PTT tt0 = find_train & ~net[iLyr][iNd].y_train;
            int n0 = kitty::count_ones( tt0 );
            int n1 = kitty::count_ones( tt1 );
            if( n1 > n0 )
            {
                net[iLyr][iNd].s_train |= find_train;
                net[iLyr][iNd].s_valid |= find_valid;
                net[iLyr][iNd].s_test |= find_test;
            }
            else if( n1 == n0 )
            {
                if( distrib(ml_gen) > 0.5 )
                {
                    net[iLyr][iNd].s_train |= find_train;
                    net[iLyr][iNd].s_valid |= find_valid;
                    net[iLyr][iNd].s_test  |= find_test;
                }
            }
            todos &= ~find_train;
        }
    }
}

void genet_t::train_network()
{
    for( int iLyr{1}; iLyr<net.size(); ++iLyr )
        for( int iNd{0}; iNd<net[iLyr].size(); ++iNd )
            if( net[iLyr][iNd].isUsed )
                train_node( iLyr, iNd );

    genet_node_t output_node = net[net.size()-1][0];
    e_train = output_node.y_train ^ output_node.s_train;
    accTrain = (double)kitty::count_ones( ~output_node.y_train ^ output_node.s_train )/ (double)output_node.s_train.num_bits();
    accValid = (double)kitty::count_ones( ~output_node.y_valid ^ output_node.s_valid )/ (double)output_node.s_valid.num_bits();
    accTest  = (double)kitty::count_ones( ~output_node.y_test  ^ output_node.s_test  )/ (double)output_node.s_test.num_bits();

    printf("Atrain %f  Avalid %f Atest %f\n", accTrain, accValid, accTest );
}

double genet_t::acc_test()
{
    return accTest;
}

double genet_t::acc_train()
{
    return accTrain;
}

double genet_t::acc_valid()
{
    return accValid;
}

void genet_t::print()
{
    for( int i{0}; i<net.size(); ++i )
    {
        printf( "LAYER %d  [%d nodes]\n", i, net[i].size() );
        for( int j{0}; j<net[i].size(); ++j )
        {
            if( net[i][j].isUsed )
            {
                printf("[ ");
                for ( int k = 0; k < net[i][j].inputs.size(); ++k )
                    printf("%3d ", net[i][j].inputs[k]);
                printf(" : %3d ] ", net[i][j].idNd );
            } 
            else
            {
                printf("{ ");
                for ( int k = 0; k < net[i][j].inputs.size(); ++k )
                    printf("%3d ", net[i][j].inputs[k]);
                printf(" : %3d } ", net[i][j].idNd );
            }           
        }
        printf("\n");
    }
}


} //namespace mcts

} //namespace mockturtle