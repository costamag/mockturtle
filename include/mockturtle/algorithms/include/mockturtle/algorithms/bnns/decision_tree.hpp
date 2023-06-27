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

namespace ccgame
{

using PTT = kitty::partial_truth_table;

class dt_node
{
private:
    uint32_t idx; // index of the input variables    
    uint32_t ctrl; // index of the input variables    
    uint32_t idx1; // index of the tree: case 1
    uint32_t idx0; // index in the tree: case 0

public:
    dt_node( uint32_t idx, uint32_t ctrl, uint32_t idx1, uint32_t idx0 ) : idx(idx), ctrl(ctrl), idx1(idx1), idx0(idx0) {};
    dt_node( uint32_t idx ) : idx(idx), ctrl(idx), idx1(idx), idx0(idx) {};
    ~dt_node(){};
    
    bool is_input() {   return (idx1 == idx0) && (idx0 == idx);    }
    int get_child0(){ return idx0; }
    int get_child1(){ return idx1; }
    int get_idx(){ return idx; }
    int get_ctrl(){ return ctrl; }

    void print()
    {
        printf("%3d=ITE( %3d, %3d, %3d )\n", idx, ctrl, idx1, idx0 );
    }

};

class decision_tree
{
private:
    std::vector<dt_node> nodes;
    std::vector<PTT> x_train;
    std::vector<PTT> y_train;
    std::vector<PTT> m_train;
    std::vector<PTT> x_test;
    std::vector<PTT> y_test;
    std::vector<PTT> m_test;
    std::vector<int> o_nodes;
public:
    decision_tree( std::vector<PTT> xtrain, std::vector<PTT> ytrain, std::vector<PTT> xtest, std::vector<PTT> ytest ) : x_train(xtrain), y_train(ytrain), x_test(xtest), y_test(ytest)
    {
        for( int i{0}; i<y_train.size(); ++i )
            m_train.push_back( ~y_train[0].construct() );
        for( int i{0}; i<y_test.size(); ++i )
            m_test.push_back( ~y_test[0].construct() );

        PTT train0 = x_train[0].construct();
        PTT test0 = x_test[0].construct();
        x_train.insert( x_train.begin(), ~train0 );
        x_train.insert( x_train.begin(), train0 );
        x_test.insert( x_test.begin(), ~test0 );
        x_test.insert( x_test.begin(),  test0 );
        for( int i{0}; i<x_train.size(); ++i )
            nodes.emplace_back( i );
    };
    decision_tree( std::vector<PTT> xtrain, std::vector<PTT> ytrain, std::vector<PTT> mtrain, std::vector<PTT> xtest, std::vector<PTT> ytest, std::vector<PTT> mtest ) : x_train(xtrain), y_train(ytrain), m_train(mtrain), x_test(xtest), y_test(ytest), m_test(mtest)
    {
        PTT train0 = x_train[0].construct();
        PTT test0 = x_test[0].construct();
        x_train.insert( x_train.begin(), ~train0 );
        x_train.insert( x_train.begin(), train0 );
        x_test.insert( x_test.begin(), ~test0 );
        x_test.insert( x_test.begin(),  test0 );
        for( int i{0}; i<x_train.size(); ++i )
            nodes.emplace_back( i );
    };
    ~decision_tree(){};


    uint32_t recursive_train_entropy( std::vector<int> supp, PTT func, PTT mask )
    {
        
        if( supp.size() == 0 )
        {
            int n0 = kitty::count_ones( mask & ~func );
            int n1 = kitty::count_ones( mask &  func );
            if( n1 > n0 )
                return 1;
            else
                return 0;

        }
        if( (kitty::count_ones( mask & func ) == 0) )
            return 0;
        if( kitty::equal( mask, mask & func ) )
            return 1;
        int idx = 0;
        double MI, bestMI = -1;
        double Nb = (double)kitty::count_ones( mask );
        double P00, P01, P10, P11, PX0, PX1, PY0, PY1;
        for( int i{0}; i<supp.size(); ++i )
        {
            int si = supp[i];
            P00 = (double)kitty::count_ones(  ~x_train[si] & ~func & mask );
            P00 = ( P00 > 0 && Nb > 0 ) ? (P00/Nb)*log2(P00/Nb) : 0 ;
            P01 = (double)kitty::count_ones(  ~x_train[si] &  func & mask );
            P01 = ( P01 > 0 && Nb > 0 ) ? (P01/Nb)*log2(P01/Nb) : 0 ;
            P10 = (double)kitty::count_ones(   x_train[si] & ~func & mask );
            P10 = ( P10 > 0 && Nb > 0 ) ? (P10/Nb)*log2(P10/Nb) : 0 ;
            P11 = (double)kitty::count_ones(   x_train[si] & func & mask );
            P11 = ( P11 > 0 && Nb > 0 ) ? (P11/Nb)*log2(P11/Nb) : 0 ;
            PX0  = (double)kitty::count_ones(  ~x_train[si] & mask );
            PX0 = ( PX0 > 0 && Nb > 0 ) ? (PX0/Nb)*log2(PX0/Nb) : 0 ;
            PX1  = (double)kitty::count_ones(   x_train[si] & mask );
            PX1 = ( PX1 > 0 && Nb > 0 ) ? (PX1/Nb)*log2(PX1/Nb) : 0 ;
            PY0  = (double)kitty::count_ones( ~func & mask );
            PY0 = ( PY0 > 0 && Nb > 0 ) ? (PY0/Nb)*log2(PY0/Nb) : 0 ;
            PY1  = (double)kitty::count_ones(  func & mask );
            PY1 = ( PY1 > 0 && Nb > 0 ) ? (PY1/Nb)*log2(PY1/Nb) : 0 ;


            double MI = P00+P01+P10+P11-PX0-PX1-PY0-PY1;
            if( MI > bestMI )
            {
                bestMI = MI;
                idx = i;
            }
        }
        int ftr = supp[idx];
        supp.erase( supp.begin() + idx );
        PTT func0 = func & ~x_train[ftr];
        PTT mask0 = mask & ~x_train[ftr];
        uint32_t idx0 = recursive_train_entropy( supp, func0, mask0 );
        PTT func1 = func &  x_train[ftr];
        PTT mask1 = mask &  x_train[ftr];
        uint32_t idx1 = recursive_train_entropy( supp, func1, mask1 );
        if( idx0 == idx1 )
            return idx0;
        nodes.emplace_back( nodes.size(), ftr, idx1, idx0 );
        return nodes.back().get_idx();
    } 

    void train_entropy()
    {
        for( int iTrg{0}; iTrg < y_train.size(); ++iTrg )
        {
            PTT mask = ~(y_train[0].construct());
            assert( kitty::is_const0(~mask) );

            std::vector<int> support_ids;
            for( int iFtr{2}; iFtr<x_train.size(); ++iFtr )
                support_ids.push_back( iFtr );
            o_nodes.push_back( recursive_train_entropy( support_ids, y_train[iTrg], mask ) );
        }
    }

   uint32_t recursive_train_random( std::vector<int> supp, PTT func, PTT mask )
    {
        
        if( supp.size() == 0 )
        {
            int n0 = kitty::count_ones( mask & ~func );
            int n1 = kitty::count_ones( mask &  func );
            if( n1 > n0 )
                return 1;
            else
                return 0;

        }
        if( (kitty::count_ones( mask & func ) == 0) )
            return 0;
        if( kitty::equal( mask, mask & func ) )
            return 1;

        std::uniform_int_distribution<> distrib(0, supp.size()-1);
        int idx = distrib(ml_gen);
        int ftr = supp[idx];

        supp.erase( supp.begin() + idx );
        PTT func0 = func & ~x_train[ftr];
        PTT mask0 = mask & ~x_train[ftr];
        uint32_t idx0 = recursive_train_entropy( supp, func0, mask0 );
        PTT func1 = func &  x_train[ftr];
        PTT mask1 = mask &  x_train[ftr];
        uint32_t idx1 = recursive_train_entropy( supp, func1, mask1 );
        if( idx0 == idx1 )
            return idx0;
        nodes.emplace_back( nodes.size(), ftr, idx1, idx0 );
        return nodes.back().get_idx();
    } 

    void train_random()
    {
        for( int iTrg{0}; iTrg < y_train.size(); ++iTrg )
        {
            PTT mask = ~(y_train[0].construct());
            assert( kitty::is_const0(~mask) );

            std::vector<int> support_ids;
            for( int iFtr{2}; iFtr<x_train.size(); ++iFtr )
                support_ids.push_back( iFtr );
            o_nodes.push_back( recursive_train_random( support_ids, y_train[iTrg], mask ) );
        }
    }

    size_t size()
    {
        return nodes.size()-x_train.size();
    }

    PTT compute_recursive( std::vector<PTT> * pX, int idx )
    {
        if( nodes[idx].get_child0() == nodes[idx].get_child1() )
            return (*pX)[nodes[idx].get_child0()];
        PTT TTC = (*pX)[nodes[idx].get_ctrl()];
        PTT TT1 = compute_recursive( pX, nodes[idx].get_child1() );
        PTT TT0 = compute_recursive( pX, nodes[idx].get_child0() );
        return TTC&TT1 | ( ~TTC )&TT0;        
    }

    std::vector<PTT> compute( std::vector<PTT> X )
    {
        assert( kitty::is_const0( X[0]));
        assert( kitty::is_const0(~X[1]));
        assert( X.size() == x_train.size() );
        std::vector<PTT> otts;
        for( int i{0}; i<o_nodes.size(); ++i )
        {
            otts.push_back( compute_recursive( &X, nodes[o_nodes[i]].get_idx() ) );
        }
        return otts;
    }

    double accuracy( std::vector<PTT> * pX, std::vector<PTT> * pY, std::vector<PTT> * pM )
    {
        assert( pM->size() == pY->size() );
        std::vector<PTT> sim = compute( *pX );
        PTT avec = ~(*pY)[0].construct();
        double nData = kitty::count_ones( avec );
        for( int i{0}; i<pY->size(); ++i )
            avec ^= ( (*pM)[i] & ( (*pY)[i]^sim[i] ) ); 
        return (double)kitty::count_ones( avec )/nData;
    } 

    double train_accuracy()
    {
        return accuracy( &x_train, &y_train, &m_train );
    }

    double test_accuracy()
    {
        return accuracy( &x_test, &y_test, &m_test );
    }

    void print()
    {
        for( int i{0}; i<nodes.size(); ++i )
            nodes[i].print();
    }
};


} // namespace ccgame

} // namespace mockturtle