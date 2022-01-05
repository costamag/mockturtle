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
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>



#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>

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
  
  kitty::dynamic_truth_table table( 5u );
  kitty::create_from_expression( table, "{((ab)c)(de)}" );

  std::cout << " SHANNON: " << std::endl;
  auto inodes = prepare_inodes(5);
  auto onodes = prepare_onodes({0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1});

  std::cout << " MUESLI: " << std::endl;

  plaT_network pla( inodes, onodes, 5, 2 );
  pla.print_pla();
  pla.muesli();

  std::cout << "functionality equivalence check" << std::endl;
  aig_network aig = convert_klut_to_graph<aig_network>( pla.klut );

  default_simulator<kitty::dynamic_truth_table> sim( table.num_vars() );
  kitty::print_binary(table);
  std::cout << std::endl;
  std::cout << "ARE EQUAL?: " << (simulate<kitty::dynamic_truth_table>( aig, sim )[0] == table ) << std::endl;

  std::cout << " MUESLI MODIFIED: " << std::endl;
  plaT_network pla_mod( inodes, onodes, 5, 2 );
  pla_mod.print_pla();
  pla_mod.muesli_modified();

  std::cout << " MUESLI PREPROCESSED: " << std::endl;
  inodes = prepare_inodes(5);
  onodes = prepare_onodes({0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1});
  plaT_network pla_pre( inodes, onodes, 5, 3 );
  pla_pre.print_pla();
  pla_pre.preprocess_muesli();
  pla_pre.muesli();

  std::cout << "functionality equivalence check" << std::endl;
  aig_network aig_pre = convert_klut_to_graph<aig_network>( pla_pre.klut );

  default_simulator<kitty::dynamic_truth_table> sim_pre( table.num_vars() );
  kitty::print_binary(table);
  std::cout << std::endl;
  std::cout << "ARE EQUAL?: " << (simulate<kitty::dynamic_truth_table>( aig_pre, sim_pre )[0] == table ) << std::endl;

  std::cout << "simulation time:" << std::endl;

  boost::dynamic_bitset<> inpt( 5, 31 );
  std::cout << pla_pre.simulate_input( inpt ) << std::endl;
  boost::dynamic_bitset<> inpt2( 5, 30 );
  std::cout << pla_pre.simulate_input( inpt2 ) << std::endl;
  boost::dynamic_bitset<> inpt3( 5, 0 );
  std::cout << pla_pre.simulate_input( inpt3 ) << std::endl;

  std::cout << "ACCURACY: " << pla_pre.compute_accuracy( inodes, onodes ) << "%" << std::endl;

  plaT_network pla_sh( inodes, onodes, 5, 3 );
  pla_sh.it_shannon_decomposition();


/*
  std::cout << "########################################" << std::endl;
  std::cout << "--------------- ab + cd ---------------" << std::endl;
  std::cout << "########################################" << std::endl;
  std::cout << " MUESLI: " << std::endl;
  inodes = prepare_inodes(4);
  onodes = prepare_onodes({0,0,0,1,0,0,0,1,0,0,0,1,1,1,1,1});
  plaT_network pla4( inodes, onodes, 4, 2 );
  pla4.print_pla();
  pla4.muesli();
  std::cout << " MUESLI MODIFIED: " << std::endl;
  plaT_network pla4_mod( inodes, onodes, 4, 2 );
  pla4_mod.print_pla();
  pla4_mod.muesli_modified();
  
  std::cout << " MUESLI PREPROCESSED: " << std::endl;
  plaT_network pla4_pre( inodes, onodes, 4,2 );
  pla4_pre.print_pla();
  pla4_pre.preprocess_muesli();
  pla4_pre.muesli();*/
  return 0;
}