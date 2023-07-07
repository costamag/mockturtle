
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
    int nBits=4;
    std::vector<std::vector<double>> accuracies;
    std::vector<int> vNumData;
    std::vector<double> avgs;
    std::vector<double> stds;
    for( int nData{5}; nData < pow(2,2*nBits); nData=nData+5 )
    {
        std::vector<double> acc;
        for( int it{0}; it < 20; ++it )
        {
            binary_classification_dataset_t data = linearly_separable_dataset( nBits, nData, 10, 10000 );

            genet_data_t genet_train( data.x_train, data.y_train, data.m_train );
            genet_data_t genet_valid( data.x_valid, data.y_valid, data.m_valid );
            genet_data_t genet_test( data.x_test, data.y_test, data.m_test );
            genet_ps_t genet_ps;
            genet_ps.K = 5;
            genet_ps.specs = {128, 128, 128, 128};

            genet_t genet( genet_train, genet_valid, genet_test, genet_ps );
            genet.create_network<genet_netcreator_t::CREA_RAND>();
            //genet.print();

            genet.train_network();
//            evolutG_ps_t evolutSA_ps;
//            evolutSA_ps.PZ = 1.;
//            evolutSA_ps.P0 = 1.;
//            evolutSA_ps.frac = 1.;
//            evolutSA_ps.nInd = 5;
//            evolutG_t SA( genet, evolutSA_ps );
//
//            SA.simulated_annealing();
            acc.push_back( genet.acc_test() );
        }
        accuracies.push_back(acc);
        vNumData.push_back( nData );
    }

    for( int iNumData{0}; iNumData<vNumData.size(); ++iNumData )
    {
        double average = 0;
        double std = 0;
        double num = 1.0*accuracies[iNumData].size();
        for( int iRes{0}; iRes<num; ++iRes )
        {
            average += accuracies[iNumData][iRes];
            std += accuracies[iNumData][iRes]*accuracies[iNumData][iRes];
        }
        average /= num;
        std = sqrt(std/num-average*average);
        printf("%3d : %f pm %f\n", vNumData[iNumData], average, std );
        avgs.push_back( average );
        stds.push_back( std );
    }

    printf("[");
    for( auto m : avgs )
        printf("%f,", m );
    printf("]\n");

    printf("[");
    for( auto s : stds )
        printf("%f,", s );
    printf("]\n");


    printf("[");
    for( auto n : vNumData )
        printf("%d,", n );
    printf("]\n");
    //cln.print_genome();

//printf("SIMULATED ANNEALING\n");
//    evolutG_ps_t evolutSA_ps;
//    evolutSA_ps.PZ = 1.;
//
//    evolutSA_ps.P0 = 1.;
//    evolutSA_ps.frac = 1.;
//    evolutSA_ps.nInd = 1;
//    evolutG_t SA( genet, evolutSA_ps );
//
//    SA.simulated_annealing();
//
//    printf("RESULT::: ");
//    printf("Atrain = %f  Avalid = %f Atest = %f \n", SA.bestInd.acc_train() , SA.bestInd.acc_valid(), SA.bestInd.acc_test() );

    return 0;
}
