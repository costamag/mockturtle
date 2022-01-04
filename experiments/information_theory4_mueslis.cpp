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
#include <mockturtle/networks/plaT.hpp>

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

int main()
{
  std::cout << "STUDY #1: Comparing the different muesli algorithms" << std::endl;
  std::cout << "########################################" << std::endl;
  std::cout << "--------------- ab + cde ---------------" << std::endl;
  std::cout << "########################################" << std::endl;
  std::cout << " MUESLI: " << std::endl;
  auto inodes = prepare_inodes(5);
  auto onodes = prepare_onodes({0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1});
  pla_network pla( inodes, onodes, 5, 2 );
  pla.print_pla();
  pla.muesli();
  std::cout << " MUESLI MODIFIED: " << std::endl;
  pla_network pla_mod( inodes, onodes, 5, 2 );
  pla_mod.print_pla();
  pla_mod.muesli_modified();

  std::cout << " MUESLI PREPROCESSED: " << std::endl;
  inodes = prepare_inodes(5);
  onodes = prepare_onodes({0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1});
  pla_network pla_pre( inodes, onodes, 5, 3 );
  pla_pre.print_pla();
  pla_pre.preprocess_muesli();
  pla_pre.muesli();
/*
  std::cout << "########################################" << std::endl;
  std::cout << "--------------- ab + cd ---------------" << std::endl;
  std::cout << "########################################" << std::endl;
  std::cout << " MUESLI: " << std::endl;
  inodes = prepare_inodes(4);
  onodes = prepare_onodes({0,0,0,1,0,0,0,1,0,0,0,1,1,1,1,1});
  pla_network pla4( inodes, onodes, 4, 2 );
  pla4.print_pla();
  pla4.muesli();
  std::cout << " MUESLI MODIFIED: " << std::endl;
  pla_network pla4_mod( inodes, onodes, 4, 2 );
  pla4_mod.print_pla();
  pla4_mod.muesli_modified();
  
  std::cout << " MUESLI PREPROCESSED: " << std::endl;
  pla_network pla4_pre( inodes, onodes, 4,2 );
  pla4_pre.print_pla();
  pla4_pre.preprocess_muesli();
  pla4_pre.muesli();*/
  return 0;
}