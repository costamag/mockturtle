#include "experiments.hpp"
#include <mockturtle/traits.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/resyn_engines/xag_resyn.hpp>
//#include <mockturtle/algorithms/experimental/contest.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/detail/mffc_utils.hpp>
#include <mockturtle/algorithms/lfe/sim_muesli.hpp>
#include <mockturtle/algorithms/lfe/hyperdimensional_computing/model.hpp>
//#include <mockturtle/algorithms/lfe/sim_decomposition.hpp>
//#include <mockturtle/algorithms/lfe/sim_decomposition_dp.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <lorina/aiger.hpp>
#include <kitty/statistics.hpp>
#include <iostream>
#include <string>
#include <set>

using namespace mockturtle;
using namespace hdc;
//using namespace mockturtle::experimental;
using namespace experiments;



int main()
{
  uint32_t N=16u;
  kitty::partial_truth_table X0(N);
  kitty::create_from_binary_string( X0, "1010101010101010" );
  kitty::partial_truth_table X1(N);
  kitty::create_from_binary_string( X1, "1100110011001100" );
  kitty::partial_truth_table X2(N);
  kitty::create_from_binary_string( X2, "1111000011110000" );
  kitty::partial_truth_table X3(N);
  kitty::create_from_binary_string( X3, "1111111100000000" );

  kitty::partial_truth_table F(N);
  F =  X1&( X0 | ( X2 ^ X3 ) );
  std::cout << "F : ";
  kitty::print_binary( F );
  std::cout << std::endl;

  std::cout << "X3: ";
  kitty::print_binary( X3 );
  std::cout << std::endl;

  std::cout << "X2: ";
  kitty::print_binary( X2 );
  std::cout << std::endl;

  std::cout << "X1: ";
  kitty::print_binary( X1 );
  std::cout << std::endl;

  std::cout << "X0: ";
  kitty::print_binary( X0 );
  std::cout << std::endl;


  std::vector< kitty::partial_truth_table * > X = {&X0};
  std::cout << "I(X0;F)=" << kitty::mutual_information(X,&F); std::cout << std::endl;
  X={&X1};
  std::cout << "I(X1;F)=" << kitty::mutual_information(X,&F); std::cout << std::endl;
  X={&X2};
  std::cout << "I(X2;F)=" << kitty::mutual_information(X,&F); std::cout << std::endl;
  X={&X3};
  std::cout << "I(X3;F)=" << kitty::mutual_information(X,&F); std::cout << std::endl;
  
  X={&X3,&X2};
  std::cout << "I(X2,X3;F)=" << kitty::mutual_information(X,&F); std::cout << std::endl;  

  kitty::partial_truth_table U = X2^X3;

  X={&U};
  std::cout << "I(X2^X3;F)=" << kitty::mutual_information(X,&F); std::cout << std::endl;

  X={&U,&X2};
  std::cout << "I(X2,X2^X3;F)=" << kitty::mutual_information(X,&F); std::cout << std::endl;

  X={&U,&X3};
  std::cout << "I(X3,X2^X3;F)=" << kitty::mutual_information(X,&F); std::cout << std::endl;

  X={&U,&X2,&X3};
  std::cout << "I(X2^X3,X2,X3;F)=" << kitty::mutual_information(X,&F); std::cout << std::endl;

  F =  (X1^X0)&(X2^X3) | (X2&X3&~(X1^X0));

  std::cout << "F : ";
  kitty::print_binary( F );
  std::cout << std::endl;

  std::cout << "X3: ";
  kitty::print_binary( X3 );
  std::cout << std::endl;

  std::cout << "X2: ";
  kitty::print_binary( X2 );
  std::cout << std::endl;

  std::cout << "X1: ";
  kitty::print_binary( X1 );
  std::cout << std::endl;

  std::cout << "X0: ";
  kitty::print_binary( X0 );
  std::cout << std::endl;
  
  X = {&X0};
  std::cout << "I(X0;F)=" << kitty::mutual_information(X,&F); std::cout << std::endl;
  X={&X1};
  std::cout << "I(X1;F)=" << kitty::mutual_information(X,&F); std::cout << std::endl;
  X={&X2};
  std::cout << "I(X2;F)=" << kitty::mutual_information(X,&F); std::cout << std::endl;
  X={&X3};
  std::cout << "I(X3;F)=" << kitty::mutual_information(X,&F); std::cout << std::endl;
  
  U = X1^X0;
  X={&U};
  std::cout << "I(X1^X0;F)=" << kitty::mutual_information(X,&F); std::cout << std::endl;

  X={&U,&X1};
  std::cout << "I(X1,X1^X0;F)=" << kitty::mutual_information(X,&F); std::cout << std::endl;

  X={&U,&X0};
  std::cout << "I(X0,X1^X0;F)=" << kitty::mutual_information(X,&F); std::cout << std::endl;

  X={&X1,&X0};
  std::cout << "I(X0,X1;F)=" << kitty::mutual_information(X,&F); std::cout << std::endl;

  X={&U,&X1,&X0};
  std::cout << "I(X0^X1,X0,X1;F)=" << kitty::mutual_information(X,&F); std::cout << std::endl;


  return 0;
}