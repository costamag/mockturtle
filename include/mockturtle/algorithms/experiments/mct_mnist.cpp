
#pragma once

#include <mockturtle/algorithms/bnns/mnist_manager.hpp>
#include <mockturtle/algorithms/bnns/decision_tree.hpp>
#include <kitty/partial_truth_table.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <bit>
#include <bitset>
#include <cstdint>
#include <fmt/format.h>

using namespace mockturtle;
using namespace ccgame;

int main( int argc, char ** argv )
{
    std::vector<PTT> x_train = read_mnist_image_bin( "../experiments/MNIST/train-images.idx3-ubyte", 1000 );
    std::vector<PTT> y_train = read_mnist_label_04_59( "../experiments/MNIST/train-labels.idx1-ubyte", 1000 );
    std::vector<PTT> y_train10 = read_mnist_label_10( "../experiments/MNIST/train-labels.idx1-ubyte", 1000 );
    std::vector<PTT> x_test = read_mnist_image_bin( "../experiments/MNIST/t10k-images.idx3-ubyte" );
    std::vector<PTT> y_test  = read_mnist_label_04_59( "../experiments/MNIST/t10k-labels.idx1-ubyte" );
    std::vector<PTT> y_test10  = read_mnist_label_10( "../experiments/MNIST/t10k-labels.idx1-ubyte" );
    print_mnist_image( &x_train, &y_train, 0 );

    decision_tree dt_h( x_train, y_train, x_test, y_test );
    dt_h.train_entropy();
    printf("size = %d\n", dt_h.size());
    printf("train acc = %f\n", dt_h.train_accuracy());
    printf("test  acc = %f\n", dt_h.test_accuracy());

    decision_tree dt_r( x_train, y_train, x_test, y_test );
    dt_r.train_random();
    printf("size = %d\n", dt_r.size());
    printf("train acc = %f\n", dt_r.train_accuracy());
    printf("test  acc = %f\n", dt_r.test_accuracy());

    PTT train1 = y_train10[1] | ~ y_train10[1];
    PTT test1 = y_test10[1] | ~ y_test10[1];
    std::vector<PTT> ytrain10 = {y_train10[1],y_train10[2],y_train10[3],y_train10[4]};
    std::vector<PTT> mtrain10 = {y_train10[0], train1, train1, train1 };
    std::vector<PTT> ytest10 = {y_test10[1],y_test10[2],y_test10[3],y_test10[4]};
    std::vector<PTT> mtest10 = {y_test10[0], test1, test1, test1 };

    decision_tree dt10_h( x_train, ytrain10, mtrain10,x_test, ytest10, mtest10 );
    dt10_h.train_entropy();
    printf("size = %d\n", dt10_h.size());
    printf("train acc = %f\n", dt10_h.train_accuracy());
    printf("test  acc = %f\n", dt10_h.test_accuracy());

    for( int i{0}; i<10; ++i )
    {
        for( auto b{0}; b<mtrain10.size(); ++b )
            printf("%d", kitty::get_bit( mtrain10[b], i ) == 1 ? kitty::get_bit( ytrain10[b], i ) : 8 );
        printf("\n");
    }

    decision_tree dt10_r( x_train, ytrain10, mtrain10, x_test, ytest10, mtest10 );
    dt10_r.train_random();
    printf("size = %d\n", dt10_r.size());
    printf("train acc = %f\n", dt10_r.train_accuracy());
    printf("test  acc = %f\n", dt10_r.test_accuracy());

    return 0;
}
