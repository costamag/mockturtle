#include <iostream>
//#include <catch.hpp>

#include <sstream>
#include <string>
#include <vector>
#include <bitset>
#include <cmath>
#include <random>
#include <boost/dynamic_bitset.hpp>

#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/pla.hpp>

#include <fstream>
#include <string>

using namespace mockturtle;
using namespace kitty;

std::vector<boost::dynamic_bitset<>> prepare_inodes( uint32_t nin )
{
  std::vector<boost::dynamic_bitset<>> input_nodes;
  for ( uint32_t i {0u}; i < pow( 2, nin ); ++i )
  {
    boost::dynamic_bitset<> idata( (nin+1), i );
    input_nodes.push_back( idata );
  }
  return input_nodes;
}

std::vector<boost::dynamic_bitset<>> prepare_onodes( std::vector<uint32_t> Voutput_nodes )
{
  std::vector<boost::dynamic_bitset<>> output_nodes;

  for( uint32_t k {0u}; k < Voutput_nodes.size(); ++k )
  {
    boost::dynamic_bitset<> odata( 1, Voutput_nodes.at(k) );
    output_nodes.push_back( odata );
  }
  return output_nodes;
}

void an_binary_fn(std::vector<uint32_t> vFn )
{
  auto input_nodes = prepare_inodes((uint32_t)log2(vFn.size()));
  auto output_nodes = prepare_onodes( vFn );
  pla_network pla( input_nodes, output_nodes, 3 );

  pla.print_pla();

  std::cout << "mi(x1;f)=" << pla.MI({0},{0}) << std::endl;
  std::cout << "mi(x2;f)=" << pla.MI({1},{0}) << std::endl;
  std::cout << "mi(x1,x2;f)=" << pla.MI({0,1},{0}) << std::endl;
}

void an_ternary_fn(std::vector<uint32_t> vFn )
{
  auto input_nodes = prepare_inodes((uint32_t)log2(vFn.size()));
  auto output_nodes = prepare_onodes( vFn );
  pla_network pla( input_nodes, output_nodes, 3 );

  pla.print_pla();

  std::cout << "mi(x0;f)=" << pla.MI({0},{0}) << std::endl;
  std::cout << "mi(x1;f)=" << pla.MI({1},{0}) << std::endl;
  std::cout << "mi(x2;f)=" << pla.MI({2},{0}) << std::endl;

  std::cout << "mi(x0,x1;f)=" << pla.MI({0,1},{0}) << std::endl;
  std::cout << "mi(x0,x2;f)=" << pla.MI({0,2},{0}) << std::endl;
  std::cout << "mi(x1,x2;f)=" << pla.MI({1,2},{0}) << std::endl;

  std::cout << "mi(x0,x1,x2;f)=" << pla.MI({0,1,2},{0}) << std::endl;

}

void an_quaternary_fn(std::vector<uint32_t> vFn )
{
  auto input_nodes = prepare_inodes((uint32_t)log2(vFn.size()));
  auto output_nodes = prepare_onodes( vFn );
  pla_network pla( input_nodes, output_nodes, 3 );

  pla.print_pla();

  std::cout << "mi(x0;f)=" << pla.MI({0},{0}) << std::endl;
  std::cout << "mi(x1;f)=" << pla.MI({1},{0}) << std::endl;
  std::cout << "mi(x2;f)=" << pla.MI({2},{0}) << std::endl;
  std::cout << "mi(x3;f)=" << pla.MI({3},{0}) << std::endl;


  std::cout << "mi(x0,x1;f)=" << pla.MI({0,1},{0}) << std::endl;
  std::cout << "mi(x0,x2;f)=" << pla.MI({0,2},{0}) << std::endl;
  std::cout << "mi(x0,x3;f)=" << pla.MI({0,3},{0}) << std::endl;
  std::cout << "mi(x1,x2;f)=" << pla.MI({1,2},{0}) << std::endl;
  std::cout << "mi(x1,x3;f)=" << pla.MI({1,3},{0}) << std::endl;
  std::cout << "mi(x2,x3;f)=" << pla.MI({2,3},{0}) << std::endl;

  std::cout << "mi(x1,x2,x3;f)=" << pla.MI({0,1,2},{0}) << std::endl;

}

void an_quinary_fn(std::vector<uint32_t> vFn )
{
  auto input_nodes = prepare_inodes((uint32_t)log2(vFn.size()));
  auto output_nodes = prepare_onodes( vFn );
  pla_network pla( input_nodes, output_nodes, 3 );

  pla.print_pla();

  std::cout << "mi(x0;f)=" << pla.MI({0},{0}) << std::endl;
  std::cout << "mi(x1;f)=" << pla.MI({1},{0}) << std::endl;
  std::cout << "mi(x2;f)=" << pla.MI({2},{0}) << std::endl;
  std::cout << "mi(x3;f)=" << pla.MI({3},{0}) << std::endl;
  std::cout << "mi(x4;f)=" << pla.MI({4},{0}) << std::endl;


  std::cout << "mi(x0,x1;f)=" << pla.MI({0,1},{0}) << std::endl;
  std::cout << "mi(x0,x2;f)=" << pla.MI({0,2},{0}) << std::endl;
  std::cout << "mi(x0,x3;f)=" << pla.MI({0,3},{0}) << std::endl;
  std::cout << "mi(x0,x4;f)=" << pla.MI({0,4},{0}) << std::endl;
  std::cout << "mi(x1,x2;f)=" << pla.MI({1,2},{0}) << std::endl;
  std::cout << "mi(x1,x3;f)=" << pla.MI({1,3},{0}) << std::endl;
  std::cout << "mi(x1,x4;f)=" << pla.MI({1,4},{0}) << std::endl;
  std::cout << "mi(x2,x3;f)=" << pla.MI({2,3},{0}) << std::endl;
  std::cout << "mi(x2,x4;f)=" << pla.MI({2,4},{0}) << std::endl;
  std::cout << "mi(x3,x4;f)=" << pla.MI({3,4},{0}) << std::endl;

  std::cout << "mi(x1,x2,x3;f)=" << pla.MI({0,1,2,3,4},{0}) << std::endl;

}

int main()
{

std::cout << "\nclass 1" << std::endl;
  an_binary_fn( {0,0,0,1} );
  an_binary_fn( {0,0,1,0} );
  an_binary_fn( {0,1,0,0} );
  an_binary_fn( {1,0,0,0} );
  an_binary_fn( {1,1,1,0} );
  an_binary_fn( {1,1,0,1} );
  an_binary_fn( {1,0,1,1} );
  an_binary_fn( {0,1,1,1} );
  std::cout << "\nclass 2" << std::endl;
  an_binary_fn( {0,0,1,1} );
  an_binary_fn( {1,1,0,0} );
  an_binary_fn( {0,1,0,1} );
  an_binary_fn( {1,0,1,0} );
  std::cout << "\nclass 3" << std::endl;
  an_binary_fn( {1,0,0,1} );
  an_binary_fn( {0,1,1,0} );
  std::cout << "\nclass 4" << std::endl;
  an_binary_fn( {0,0,0,0} );
  an_binary_fn( {1,1,1,1} );






  std::cout << "--------------- x1 --------------- " << std::endl;
  an_binary_fn( {0,1,0,1} );

  std::cout << "--------------- (x1)' ---------------" << std::endl;
  an_binary_fn( {1,0,1,0} );

  std::cout << "--------------- 2-AND ---------------" << std::endl;
  an_binary_fn( {0,0,0,1} );

  std::cout << "--------------- 2-OR ---------------" << std::endl;
  an_binary_fn( {0,1,1,1} );

  std::cout << "--------------- 2-XOR ---------------" << std::endl;
  an_binary_fn( {0,1,1,0} );

  std::cout << "--------------- 2-NAND ---------------" << std::endl;
  an_binary_fn( {1,1,1,0} );


  std::cout << "--------------- ab + cde ---------------" << std::endl;
  auto inodes = prepare_inodes(5);
  auto onodes = prepare_onodes({0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1});
  pla_network pla( inodes, onodes, 5 );
  pla.print_pla();
  //pla.muesli(2);
  an_quinary_fn({0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1});
  std::cout << "###############################################\n";
  pla_network pla2( inodes, onodes, 4 );
  pla2.print_pla();
  std::cout << "MI:\n";
  std::cout << "mi(x0;f)=" << pla2.MI({0},{0}) << std::endl;
  std::cout << "mi(x1;f)=" << pla2.MI({1},{0}) << std::endl;
  std::cout << "mi(x2;f)=" << pla2.MI({2},{0}) << std::endl;
  std::cout << "mi(x3;f)=" << pla2.MI({3},{0}) << std::endl;
  std::cout << "mi(x4;f)=" << pla2.MI({4},{0}) << std::endl;
  std::cout << "pi1 = {x0,x1,x2}, pi2 = {x3,x4}\n";
  std::cout << "pi1: create functions by symmetry classes" << std::endl;
  std::cout << "Pi1:" << std::endl;
  auto tt_tmp = pla2.create_fn({0,1,2});
  std::cout << tt_tmp << std::endl;
  pla2.create_klut_node( {0,1,2}, tt_tmp );

  tt_tmp = pla2.create_fn({3,4});
  std::cout << tt_tmp << std::endl;
  pla2.create_klut_node( {3,4}, tt_tmp );


  tt_tmp = pla2.create_fn({5,6});
  std::cout << tt_tmp << std::endl;
  pla2.create_klut_node( {5,6}, tt_tmp );

  std::cout << "--------------- ab + cd ---------------" << std::endl;
   inodes = prepare_inodes(4);
  an_quaternary_fn( {0,0,0,1,0,0,0,1,0,0,0,1,1,1,1,1} );
   onodes =prepare_onodes({0,0,0,1,0,0,0,1,0,0,0,1,1,1,1,1});
  pla_network pla1( inodes, onodes, 4 );
  //pla1.print_pla();
  std::cout << "###############################################\n";
  pla_network pla3( inodes, onodes, 4 );
  pla3.print_pla();
  std::cout << "MI:\n";
  std::cout << "mi(x0;f)=" << pla3.MI({0},{0}) << std::endl;
  std::cout << "mi(x1;f)=" << pla3.MI({1},{0}) << std::endl;
  std::cout << "mi(x2;f)=" << pla3.MI({2},{0}) << std::endl;
  std::cout << "mi(x3;f)=" << pla3.MI({3},{0}) << std::endl;
  std::cout << "pi1 = {x0,x1,x2,x3}\n";
  std::cout << "pi1: support larger than 2: split more" << std::endl;
  std::cout << "Pi1 a: pick x0" << std::endl;
  std::cout << "mi(x0,x1;f)=" << pla3.MI({0,1},{0}) << std::endl;
  std::cout << "mi(x0,x2;f)=" << pla3.MI({0,2},{0}) << std::endl;
  std::cout << "mi(x0,x3;f)=" << pla3.MI({0,3},{0}) << std::endl;
  std::cout << "Pi1 a: x0,x1" << std::endl;
  std::cout << "Pi1 not void, still x2,x3. size is equal to max so put in pi1 b: {x2 x3}" << std::endl;


/*  auto tt_tmp = pla2.create_fn({0,1,2});
  std::cout << tt_tmp << std::endl;
  pla2.create_klut_node( {0,1,2}, tt_tmp );

  tt_tmp = pla2.create_fn({3,4});
  std::cout << tt_tmp << std::endl;
  pla2.create_klut_node( {3,4}, tt_tmp );


  tt_tmp = pla2.create_fn({5,6});
  std::cout << tt_tmp << std::endl;
  pla2.create_klut_node( {5,6}, tt_tmp );*/

  std::cout << "--------------- 3-111 ---------------" << std::endl;
  an_ternary_fn( {0,0,0,0,0,0,0,1} );
  std::cout << "--------------- 3-(111)' ---------------" << std::endl;
  an_ternary_fn( {1,1,1,1,1,1,1,0} );

  std::cout << "--------------- 3-00* ---------------" << std::endl;
  an_ternary_fn( {1,1,0,0,0,0,0,0} );
  std::cout << "--------------- 3-000 101 ---------------" << std::endl;
  an_ternary_fn( {1,0,0,0,0,1,0,0} );
  std::cout << "--------------- 3-000 111 ---------------" << std::endl;
  an_ternary_fn( {1,0,0,0,0,0,0,1} );
  std::cout << "--------------- 3-000 100 001 ---------------" << std::endl;
  an_ternary_fn( {1,1,0,0,1,0,0,0} );
  std::cout << "--------------- 3-000 001 110 ---------------" << std::endl;
  an_ternary_fn( {1,1,0,0,0,0,1,0} );
  std::cout << "--------------- 3-000 011 101 ---------------" << std::endl;
  an_ternary_fn( {1,0,0,1,0,1,0,0} );
  std::cout << "--------------- 3-000 001 010 100 ---------------" << std::endl;
  an_ternary_fn( {1,1,1,0,1,0,0,0} );
  std::cout << "--------------- 3-000 001 010 111 ---------------" << std::endl;
  an_ternary_fn( {0,1,1,0,1,0,0,1} );
  std::cout << "--------------- 3-001 010 100 101 ---------------" << std::endl;
  an_ternary_fn( {0,1,1,0,1,1,0,0} );
  std::cout << "--------------- 3-000 001 010 101 ---------------" << std::endl;
  an_ternary_fn( {1,1,1,0,0,1,0,0} );

  std::cout << "before" << std::endl;
  an_quinary_fn({0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1});
  std::cout << "after" << std::endl;


  std::cout << "--------------- ab + cde ---------------" << std::endl;
  inodes = prepare_inodes(5);
  onodes = prepare_onodes({0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1});
  pla_network pla5( inodes, onodes, 5 );
  pla5.print_pla();
  pla5.preprocess_muesli(3);
  pla5.muesli(5);

    std::cout << "--------------- ab + cd ---------------" << std::endl;
   inodes = prepare_inodes(4);
  an_quaternary_fn( {0,0,0,1,0,0,0,1,0,0,0,1,1,1,1,1} );
   onodes =prepare_onodes({0,0,0,1,0,0,0,1,0,0,0,1,1,1,1,1});
  pla_network pla6( inodes, onodes, 4 );
  pla6.print_pla();
  pla6.preprocess_muesli(3);
  pla6.muesli(5);

    pla_network pla7( inodes, onodes, 4 );
  pla7.print_pla();
  pla7.it_shannon_decomposition(2, 0);


  return 0;
}