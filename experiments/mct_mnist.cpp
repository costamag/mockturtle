
#pragma once

#include <mockturtle/algorithms/mcts/mnist_manager.hpp>
#include <mockturtle/algorithms/mcts/decision_tree.hpp>
#include <mockturtle/algorithms/mcts/genet.hpp>
#include <mockturtle/algorithms/mcts/evolut.hpp>
#include <mockturtle/algorithms/mcts/evolutG.hpp>
#include <mockturtle/algorithms/mcts/clnet.hpp>
#include <kitty/partial_truth_table.hpp>
#include <kitty/partial_truth_table.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <bit>
#include <bitset>
#include <cstdint>
#include <fmt/format.h>

using namespace mockturtle;
using namespace mcts;

int main( int argc, char ** argv )
{
    int nTrain = 128;
    int nValid = 1000;
    int nTest = 1000;

    std::vector<PTT> Xtrain = read_mnist_image_bin( "../experiments/MNIST/train-images.idx3-ubyte", nTrain+nValid );
    std::vector<PTT> Ytrain = read_mnist_label_04_59( "../experiments/MNIST/train-labels.idx1-ubyte", nTrain+nValid );
    std::vector<PTT> x_test = read_mnist_image_bin( "../experiments/MNIST/t10k-images.idx3-ubyte", nTest );
    std::vector<PTT> y_test  = read_mnist_label_04_59( "../experiments/MNIST/t10k-labels.idx1-ubyte", nTest );
    //std::vector<PTT> y_train10 = read_mnist_label_10( "../experiments/MNIST/train-labels.idx1-ubyte", 1000 );
    //std::vector<PTT> y_test10  = read_mnist_label_10( "../experiments/MNIST/t10k-labels.idx1-ubyte" );

    //std::vector<PTT> x_valid = read_mnist_image_bin( "../experiments/MNIST/train-images.idx3-ubyte", 2000 );
    //std::vector<PTT> y_valid = read_mnist_label_04_59( "../experiments/MNIST/train-labels.idx1-ubyte", 2000 );

    std::vector<PTT> x_train;
    std::vector<PTT> m_train;
    std::vector<PTT> y_train;
    std::vector<PTT> x_valid;
    std::vector<PTT> m_valid;
    std::vector<PTT> y_valid;

    std::vector<PTT> m_test;
    for( auto tt : y_test )
        m_test.push_back( ~tt.construct() );

    for( int iX{0}; iX<Xtrain.size(); ++iX )
    {
        PTT tt_train;
        PTT tt_valid;
        for( int i{0}; i<nTrain+nValid; ++i )
        {
            if( i < nTrain )
                tt_train.add_bit( kitty::get_bit( Xtrain[iX], i ) );
            else
                tt_valid.add_bit( kitty::get_bit( Xtrain[iX], i ) );
        }
        x_train.push_back( tt_train );
        x_valid.push_back( tt_valid );
    }
    for( int iY{0}; iY<Ytrain.size(); ++iY )
    {
        PTT tt_train;
        PTT tt_valid;
        for( int i{0}; i<nTrain+nValid; ++i )
        {
            if( i < nTrain )
                tt_train.add_bit( kitty::get_bit( Ytrain[iY], i ) );
            else
                tt_valid.add_bit( kitty::get_bit( Ytrain[iY], i ) );
        }
        y_train.push_back( tt_train );
        y_valid.push_back( tt_valid );
        m_train.push_back( ~tt_train.construct() );
        m_valid.push_back( ~tt_valid.construct() );
    }
    
printf("GENET\n");
    genet_data_t genet_train(x_train, y_train[0], m_train[0] );
    genet_data_t genet_valid(x_valid, y_valid[0], m_valid[0] );
    genet_data_t genet_test(x_test, y_test[0], m_test[0] );
    genet_ps_t genet_ps;
    genet_ps.specs = { 1024, 1024, 1024, 1024, 1024 };
    genet_ps.K = 4;

    genet_t genet( genet_train, genet_valid, genet_test, genet_ps );
    genet.create_network<genet_netcreator_t::CREA_RAND>();
    //genet.print();

    print_mnist_image( &x_train, &y_train, 0 );
    genet.train_network();

    //cln.print_genome();
printf("EVOLUT\n");
    evolutG_ps_t evolutG_ps;
    evolutG_ps.P0 = 0.9;
    evolutG_ps.PZ = 0;
    evolutG_ps.frac = 1.;
    evolutG_t evolutG( genet, evolutG_ps );

    evolutG.train();

    printf("RESULT::: ");
    printf("Atrain = %f  Avalid = %f Atest = %f \n", evolutG.bestInd.acc_train() , evolutG.bestInd.acc_valid(), evolutG.bestInd.acc_test() );

printf("MUTUAL INFORMATION\n");
    decision_tree dt_h( x_train, y_train, x_test, y_test );
    dt_h.train_impurity<entropy_t::MINF>();
    printf("size = %d\n", dt_h.size());
    printf("train acc = %f\n", dt_h.train_accuracy());
    printf("test  acc = %f\n", dt_h.test_accuracy());
printf("GINI\n");
    decision_tree dt_g( x_train, y_train, x_test, y_test );
    dt_g.train_impurity<entropy_t::GINI>();
    printf("size = %d\n", dt_g.size());
    printf("train acc = %f\n", dt_g.train_accuracy());
    printf("test  acc = %f\n", dt_g.test_accuracy());
printf("SHANNON\n");
    decision_tree dt_s( x_train, y_train, x_test, y_test );
    dt_s.train_impurity<entropy_t::SHAN>();
    printf("size = %d\n", dt_s.size());
    printf("train acc = %f\n", dt_s.train_accuracy());
    printf("test  acc = %f\n", dt_s.test_accuracy());
printf("0-1\n");
    decision_tree dt_01( x_train, y_train, x_test, y_test );
    dt_01.train_impurity<entropy_t::EN01>();
    printf("size = %d\n", dt_01.size());
    printf("train acc = %f\n", dt_01.train_accuracy());
    printf("test  acc = %f\n", dt_01.test_accuracy());

printf("ORDERED\n");
    decision_tree dt_o( x_train, y_train, x_test, y_test );
    dt_o.train_ordered();
    printf("size = %d\n", dt_o.size());
    printf("train acc = %f\n", dt_o.train_accuracy());
    printf("test  acc = %f\n", dt_o.test_accuracy());
printf("RANDOM\n");

    decision_tree dt_r( x_train, y_train, x_test, y_test );
    dt_r.train_random();
    printf("size = %d\n", dt_r.size());
    printf("train acc = %f\n", dt_r.train_accuracy());
    printf("test  acc = %f\n", dt_r.test_accuracy());
//printf("GENETIC\n");
//    clnet_ps clps;
//    clps.nFilters = 10;
//    clps.nGen0 = 10;
//    clps.nGenerations = 10;
//    clnet cln( x_train, y_train, x_valid, y_valid, x_test, y_test, clps );
//    cln.train();
//    PTT train1 = y_train10[1] | ~ y_train10[1];
//    PTT test1 = y_test10[1] | ~ y_test10[1];
//    std::vector<PTT> ytrain10 = {y_train10[1],y_train10[2],y_train10[3],y_train10[4]};
//    std::vector<PTT> mtrain10 = {y_train10[0], train1, train1, train1 };
//    std::vector<PTT> ytest10 = {y_test10[1],y_test10[2],y_test10[3],y_test10[4]};
//    std::vector<PTT> mtest10 = {y_test10[0], test1, test1, test1 };
//
//    decision_tree dt10_h( x_train, ytrain10, mtrain10,x_test, ytest10, mtest10 );
//    dt10_h.train_impurity<entropy_t::MINF>();
//    printf("size = %d\n", dt10_h.size());
//    printf("train acc = %f\n", dt10_h.train_accuracy());
//    printf("test  acc = %f\n", dt10_h.test_accuracy());
//
//    for( int i{0}; i<10; ++i )
//    {
//        for( auto b{0}; b<mtrain10.size(); ++b )
//            printf("%d", kitty::get_bit( mtrain10[b], i ) == 1 ? kitty::get_bit( ytrain10[b], i ) : 8 );
//        printf("\n");
//    }
//
//    decision_tree dt10_r( x_train, ytrain10, mtrain10, x_test, ytest10, mtest10 );
//    dt10_r.train_random();
//    printf("size = %d\n", dt10_r.size());
//    printf("train acc = %f\n", dt10_r.train_accuracy());
//    printf("test  acc = %f\n", dt10_r.test_accuracy());
//printf("GINI\n");
//    decision_tree dt10_g( x_train, ytrain10, mtrain10, x_test, ytest10, mtest10 );
//    dt10_g.train_random();
//    printf("size = %d\n", dt10_g.size());
//    printf("train acc = %f\n", dt10_g.train_accuracy());
//    printf("test  acc = %f\n", dt10_g.test_accuracy());
//printf("SHAN\n");
//    decision_tree dt10_s( x_train, ytrain10, mtrain10, x_test, ytest10, mtest10 );
//    dt10_s.train_random();
//    printf("size = %d\n", dt10_s.size());
//    printf("train acc = %f\n", dt10_s.train_accuracy());
//    printf("test  acc = %f\n", dt10_s.test_accuracy());
//printf("EN01\n");
//    decision_tree dt10_e( x_train, ytrain10, mtrain10, x_test, ytest10, mtest10 );
//    dt10_e.train_random();
//    printf("size = %d\n", dt10_e.size());
//    printf("train acc = %f\n", dt10_e.train_accuracy());
//    printf("test  acc = %f\n", dt10_e.test_accuracy());

    return 0;
}
