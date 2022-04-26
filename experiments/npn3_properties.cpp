
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

#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/io/truth_reader.hpp>


//#include <mockturtle/networks/cover.hpp>
#include <mockturtle/algorithms/lfe/mi_decomposition.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
//#include <mockturtle/algorithms/muesli.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>

//#include <mockturtle/networks/pla2.hpp>// plaT for bottom, plaT0 for only greedy
//#include <mockturtle/networks/plaT0.hpp>
#include <lorina/truth.hpp>
#include <mockturtle/views/names_view.hpp>
#include <mockturtle/algorithms/lfe/graph_to_lfe.hpp>
#include <typeinfo>
#include <lorina/lorina.hpp>
#include <mockturtle/views/depth_view.hpp>

#include <kitty/constructors.hpp>
#include <kitty/spectral.hpp>
#include <kitty/print.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/statistics.hpp>
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
      kitty::print_binary(x); std::cout <<std::endl;
    }
    n = LFE.partial.first[0].num_bits();
    for( auto i = 0u; i < n; ++i )
      std::cout << "-";
    std::cout << std::endl;

    kitty::print_binary(LFE.partial.second); std::cout <<std::endl; 
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
    auto I = kitty::mutual_information( x, LFE.partial.second );
    std::cout << " " << I << std::endl;
  }

  auto x = LFE.partial.first[0];
  auto y = LFE.partial.first[1];

  std::cout << 0 << " " << 1 << " ";
  auto W = std::vector{x,y};
  auto I = kitty::mutual_information( W, LFE.partial.second );
  std::cout << " " << I << std::endl;

}

void print_mutual_informations3( lfeNtk<klut_network> LFE )
{
  for( auto i = 0; i<LFE.partial.first.size(); ++i )
  {
    auto x = LFE.partial.first[i];
    std::cout << i << " ";
    auto I = kitty::mutual_information( x, LFE.partial.second );
    std::cout << " " << I << std::endl;
  }

  for( auto i = 0; i<LFE.partial.first.size(); ++i )
  {
    auto x = LFE.partial.first[i];
    for( auto j = 0; j<i; ++j )
    {
      auto y = LFE.partial.first[j];
      std::cout << i << " " << j << " ";
      auto W = std::vector{x,y};
      auto I = kitty::mutual_information( W, LFE.partial.second );
      std::cout << " " << I << std::endl;
    }
  }

  auto x = LFE.partial.first[0];
  auto y = LFE.partial.first[1];
  auto z = LFE.partial.first[2];

  std::cout << 0 << " " << 1 << " " << 2 << " ";
  auto W = std::vector{x,y,z};
  auto I = kitty::mutual_information( W, LFE.partial.second );
  std::cout << " " << I << std::endl;

}

struct divisors
{
  std::string tt;
  dbitset fn;
};


int main()
{

  std::cout << "#######################################################" <<std::endl;
  std::cout << "                           f = ab                      " <<std::endl;
  std::cout << "#######################################################" <<std::endl;
  std::string str_code = "00";
  std::string path2 = "/home/acostama/projects/EPFL/mockturtle/benchmarks/NPN-representatives/nin2/ex"+str_code+".truth";

  klut_network klut2_0;
  if( lorina::read_truth( path2, truth_reader( klut2_0 ) ) == lorina::return_code::parse_error )
    assert( false );

  auto LFE_pre = graph_to_lfe( klut2_0 );

  print_mutual_informations2( LFE_pre );
  auto f = LFE_pre.complete.second[0];
  auto walsh = kitty::rademacher_walsh_spectrum( f );
  for( auto w : walsh )
    std::cout << w << " ";
  std::cout << std::endl;
  std::cout << "a+b" << std::endl;
  kitty::dynamic_truth_table aORb(2u);
  kitty::create_from_binary_string(aORb, "1110");
  walsh = kitty::rademacher_walsh_spectrum( aORb );
  for( auto w : walsh )
    std::cout << w << " ";
  std::cout << std::endl;

  std::cout << "a'b" << std::endl;
  kitty::dynamic_truth_table aLTb(2u);
  kitty::create_from_binary_string(aLTb, "0010");
  walsh = kitty::rademacher_walsh_spectrum( aLTb );
  for( auto w : walsh )
    std::cout << w << " ";
  std::cout << std::endl;
  std::cout << "a'b'" << std::endl;
  kitty::dynamic_truth_table aNANDb(2u);
  kitty::create_from_binary_string(aNANDb, "0001");
  walsh = kitty::rademacher_walsh_spectrum( aNANDb );
  for( auto w : walsh )
    std::cout << w << " ";
  std::cout << std::endl;

  std::cout << "#######################################################" <<std::endl;
  std::cout << "                           f = a                       " <<std::endl;
  std::cout << "#######################################################" <<std::endl;
  str_code = "01";
  path2 = "/home/acostama/projects/EPFL/mockturtle/benchmarks/NPN-representatives/nin2/ex"+str_code+".truth";

  klut_network klut2_1;
  if( lorina::read_truth( path2, truth_reader( klut2_1 ) ) == lorina::return_code::parse_error )
    assert( false );

  LFE_pre = graph_to_lfe( klut2_1 );
  print_mutual_informations2( LFE_pre );
  f = LFE_pre.complete.second[0];
  walsh = kitty::rademacher_walsh_spectrum( f );
  for( auto w : walsh )
    std::cout << w << " ";
  std::cout << std::endl;

  std::cout << "#######################################################" <<std::endl;
  std::cout << "                           f = a^b                     " <<std::endl;
  std::cout << "#######################################################" <<std::endl;
  str_code = "02";
  path2 = "/home/acostama/projects/EPFL/mockturtle/benchmarks/NPN-representatives/nin2/ex"+str_code+".truth";

  klut_network klut2_2;
  if( lorina::read_truth( path2, truth_reader( klut2_2 ) ) == lorina::return_code::parse_error )
    assert( false );

  LFE_pre = graph_to_lfe( klut2_2 );
  print_mutual_informations2( LFE_pre );
  f = LFE_pre.complete.second[0];
  walsh = kitty::rademacher_walsh_spectrum( f );
  for( auto w : walsh )
    std::cout << w << " ";
  std::cout << std::endl;

  std::cout << "#######################################################" <<std::endl;
  std::cout << "                     f = <x y z>                       " <<std::endl;
  std::cout << "#######################################################" <<std::endl;
  str_code = "04";
  std::string path_maj = "/home/acostama/projects/EPFL/mockturtle/benchmarks/NPN-representatives/nin3/ex"+str_code+".truth";

  //muesli_params ps3;
  //ps3.max_sup = 2;
  klut_network klut_maj;
  if( lorina::read_truth( path_maj, truth_reader( klut_maj ) ) == lorina::return_code::parse_error )
    assert( false );

  auto LFE_pre_maj = graph_to_lfe( klut_maj );

  print_mutual_informations3( LFE_pre_maj );
  auto Q = std::vector{LFE_pre_maj.partial.second};
  std::cout << "H(f)= " << kitty::entropy( Q ) << std::endl;
  f = LFE_pre_maj.complete.second[0];

  auto x = LFE_pre_maj.partial.first[0];
  auto y = LFE_pre_maj.partial.first[1];
  auto z = LFE_pre_maj.partial.first[2];

  walsh = kitty::rademacher_walsh_spectrum( f );
  for( auto w : walsh )
    std::cout << w << " ";
  std::cout << std::endl;
  
  return 0;
}