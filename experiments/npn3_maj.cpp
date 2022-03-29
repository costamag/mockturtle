
#include <iostream>
//#include <catch.hpp>
//#include <execution>
#include <sstream>
#include <string>
#include <fstream>
#include <string>
#include <vector>
#include <bitset>
#include <cmath>
#include <random>
#include <boost/dynamic_bitset.hpp>

#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/io/truth_reader.hpp>


//#include <mockturtle/networks/cover.hpp>
#include <mockturtle/algorithms/it_decomposition.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/muesli.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>

//#include <mockturtle/networks/pla2.hpp>// plaT for bottom, plaT0 for only greedy
//#include <mockturtle/networks/plaT0.hpp>
#include <lorina/truth.hpp>
#include <mockturtle/views/names_view.hpp>
#include <mockturtle/algorithms/graph_to_lfe.hpp>
#include <typeinfo>
#include <lorina/lorina.hpp>
#include <mockturtle/views/depth_view.hpp>

#include <kitty/constructors.hpp>
#include <kitty/print.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/dbs_statistical_bit_operations.hpp>
#include <fstream>
#include <string>
#include <omp.h>
#include <unistd.h>
using namespace mockturtle;
void print_LFE( lfeNtk<klut_network> LFE, bool only_complete = false )  
{
  std::cout << "complete:" << std::endl;
  for( auto x : LFE.complete.first)
  {
    kitty::print_binary(x);std::cout<<std::endl;
  }
  auto n = LFE.complete.first[0].num_bits();
  for( auto i = 0u; i < n; ++i )
    std::cout << "-";
  std::cout << std::endl;
  for( auto x : LFE.complete.second)
  {
    kitty::print_binary(x);std::cout<<std::endl;
  } 
  if( !only_complete )
  {
    std::cout << "partial:" << std::endl;
    for( auto x : LFE.partial.first)
    {
      std::cout << x <<std::endl;
    }
    n = LFE.partial.first[0].size();
    for( auto i = 0u; i < n; ++i )
      std::cout << "-";
    std::cout << std::endl;

    std::cout << LFE.partial.second <<std::endl; 
  }
}

std::pair<std::vector<kitty::dynamic_truth_table>, uint32_t> load( std::string file_name )
{

  std::string line;
  std::ifstream myfile ( file_name );
  std::vector<kitty::dynamic_truth_table> tts;
  uint64_t n;
  if (myfile.is_open())
  {
    while ( getline (myfile,line) )
    {
        std::cout << line << std::endl;
        std::cout << "2^n: " << line.size() << std::endl;
        n = log2(line.size());
        std::cout << "n: " << n << std::endl;
        kitty::dynamic_truth_table tt(n);
        kitty::create_from_binary_string( tt, line );
        tts.push_back( tt );
        kitty::print_binary(tt);
        std::cout<<std::endl;
    }
    myfile.close();
    for( auto tt : tts )
    {
        print_binary(tt); std::cout << std::endl;
    }
  }
  else std::cout << "Unable to open file";
  return std::make_pair(tts,n);
}

void print_mutual_informations2( lfeNtk<klut_network> LFE )
{
  for( auto i = 0; i<LFE.partial.first.size(); ++i )
  {
    auto x = LFE.partial.first[i];
    std::cout << i << " ";
    auto I = kitty::mutual_information( std::vector{x}, LFE.partial.second );
    std::cout << " " << I << std::endl;
  }

  auto x = LFE.partial.first[0];
  auto y = LFE.partial.first[1];

  std::cout << 0 << " " << 1 << " ";
  auto I = kitty::mutual_information( std::vector{x,y}, LFE.partial.second );
  std::cout << " " << I << std::endl;

}

void print_mutual_informations3( lfeNtk<klut_network> LFE )
{
  for( auto i = 0; i<LFE.partial.first.size(); ++i )
  {
    auto x = LFE.partial.first[i];
    std::cout << i << " ";
    auto I = kitty::mutual_information( std::vector{x}, LFE.partial.second );
    std::cout << " " << I << std::endl;
  }

  for( auto i = 0; i<LFE.partial.first.size(); ++i )
  {
    auto x = LFE.partial.first[i];
    for( auto j = 0; j<i; ++j )
    {
      auto y = LFE.partial.first[j];
      std::cout << i << " " << j << " ";
      auto I = kitty::mutual_information( std::vector{x,y}, LFE.partial.second );
      std::cout << " " << I << std::endl;
    }
  }

  auto x = LFE.partial.first[0];
  auto y = LFE.partial.first[1];
  auto z = LFE.partial.first[2];

  std::cout << 0 << " " << 1 << " " << 2 << " ";
  auto I = kitty::mutual_information( std::vector{x,y,z}, LFE.partial.second );
  std::cout << " " << I << std::endl;

}

int main()
{

 

  std::cout << "#######################################################" <<std::endl;
  std::cout << "                     f = <x y z>                       " <<std::endl;
  std::cout << "#######################################################" <<std::endl;
  std::string str_code = "04";
  std::string path3 = "/home/acostama/projects/EPFL/mockturtle/benchmarks/NPN-representatives/nin3/ex"+str_code+".truth";

  muesli_params ps3;
  ps3.max_sup = 2;
  klut_network klut3;
  if( lorina::read_truth( path3, truth_reader( klut3 ) ) == lorina::return_code::parse_error )
    assert( false );

  auto LFE_pre3 = graph_to_lfe( klut3 );
  print_LFE( LFE_pre3, true );

  print_mutual_informations3( LFE_pre3 );
  std::cout << "H(f)= " << kitty::entropy( std::vector{LFE_pre3.partial.second} ) << std::endl;
  auto f = LFE_pre3.partial.second;

  auto g1 = std::vector{dbitset (8,0x88)};
  auto g2 = std::vector{dbitset (8,0xee)};
  auto g3 = std::vector{dbitset (8,0xf8)};
  auto g4 = std::vector{dbitset (8,0xe0)};
  auto x = LFE_pre3.partial.first[0];
  auto y = LFE_pre3.partial.first[1];
  auto z = LFE_pre3.partial.first[2];

  std::cout << "I(g1=xy;f)=" << kitty::mutual_information( g1, f ) << "<-" << g1[0] << std::endl;
  std::cout << "I(g2=x+y;f)=" << kitty::mutual_information( g2, f ) << "<-" << g2[0] << std::endl;
  std::cout << "\nCan g1 alone give us f?"  << std::endl;
  std::cout << "I(g1=xy;f)=" << kitty::mutual_information( g1, f ) << std::endl;
  std::cout << "I(g1=xy,x;f)=" << kitty::mutual_information( std::vector{g1[0],x}, f ) << std::endl;
  std::cout << "I(g1=xy,y;f)=" << kitty::mutual_information( std::vector{g1[0],y}, f ) << std::endl;
  std::cout << "I(g1=xy,x,y;f)=" << kitty::mutual_information( std::vector{g1[0],x,y}, f ) << std::endl;

  std::cout << "\nCan g2 alone give us f?"  << std::endl;
  std::cout << "I(g2=x+y;f)=" << kitty::mutual_information( g2, f ) << std::endl;
  std::cout << "I(g2=x+y,x;f)=" << kitty::mutual_information( std::vector{g2[0],x}, f ) << std::endl;
  std::cout << "I(g2=x+y,y;f)=" << kitty::mutual_information( std::vector{g2[0],y}, f ) << std::endl;
  std::cout << "I(g2=x+y,x,y;f)=" << kitty::mutual_information( std::vector{g2[0],x,y}, f ) << std::endl;

  std::cout << "\nCan g1 and g2 substitute x and y?"<< std::endl;
  std::cout << "I(g1=xy,g2=x+y;f)=" << kitty::mutual_information( std::vector{g1[0],g2[0]}, f ) << std::endl;
  std::cout << "I(g1=xy,g2=x+y,x;f)=" << kitty::mutual_information( std::vector{g1[0],g2[0],x}, f ) << std::endl;
  std::cout << "I(g1=xy,g2=x+y,y;f)=" << kitty::mutual_information( std::vector{g1[0],g2[0],y}, f ) << std::endl;
  std::cout << "I(g1=xy,g2=x+y,x,y;f)=" << kitty::mutual_information( std::vector{g1[0],g2[0],x,y}, f ) << std::endl;
  std::cout << "\nConsider g2z"<< std::endl;
  std::cout << "I(g4=g2z;f)=" << kitty::mutual_information( g4, f )<< "<-" << g4[0] << std::endl;
  std::cout << "\nConsider g1+z"<< std::endl;
  std::cout << "I(g3=g1+z;f)=" << kitty::mutual_information( g3, f )<< "<-" << g3[0] << std::endl;

 
  return 0;
}