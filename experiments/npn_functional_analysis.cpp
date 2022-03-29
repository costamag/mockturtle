
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
#include <kitty/operations.hpp>
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
  std::cout << "                           f = ab                      " <<std::endl;
  std::cout << "#######################################################" <<std::endl;
  std::string str_code = "00";
  std::string path2 = "/home/acostama/projects/EPFL/mockturtle/benchmarks/NPN-representatives/nin2/ex"+str_code+".truth";

  muesli_params ps2;
  ps2.max_sup = 2;
  klut_network klut2_0;
  if( lorina::read_truth( path2, truth_reader( klut2_0 ) ) == lorina::return_code::parse_error )
    assert( false );

  auto LFE_pre = graph_to_lfe( klut2_0 );
  print_LFE( LFE_pre, true );

  print_mutual_informations2( LFE_pre );
  std::cout << "H(f)= " << kitty::entropy( std::vector{LFE_pre.partial.second} ) << std::endl;


  std::cout << "#######################################################" <<std::endl;
  std::cout << "                           f = a                       " <<std::endl;
  std::cout << "#######################################################" <<std::endl;
  str_code = "01";
  path2 = "/home/acostama/projects/EPFL/mockturtle/benchmarks/NPN-representatives/nin2/ex"+str_code+".truth";

  klut_network klut2_1;
  if( lorina::read_truth( path2, truth_reader( klut2_1 ) ) == lorina::return_code::parse_error )
    assert( false );

  LFE_pre = graph_to_lfe( klut2_1 );
  print_LFE( LFE_pre, true );

  print_mutual_informations2( LFE_pre );
  std::cout << "H(f)= " << kitty::entropy( std::vector{LFE_pre.partial.second} ) << std::endl;

  std::cout << "#######################################################" <<std::endl;
  std::cout << "                           f = a^b                     " <<std::endl;
  std::cout << "#######################################################" <<std::endl;
  str_code = "02";
  path2 = "/home/acostama/projects/EPFL/mockturtle/benchmarks/NPN-representatives/nin2/ex"+str_code+".truth";

  klut_network klut2_2;
  if( lorina::read_truth( path2, truth_reader( klut2_2 ) ) == lorina::return_code::parse_error )
    assert( false );

  LFE_pre = graph_to_lfe( klut2_2 );
  print_LFE( LFE_pre, true );

  print_mutual_informations2( LFE_pre );
  std::cout << "H(f)= " << kitty::entropy( std::vector{LFE_pre.partial.second} ) << std::endl;


  std::cout << "#######################################################" <<std::endl;
  std::cout << "                         f = abc                       " <<std::endl;
  std::cout << "#######################################################" <<std::endl;
  str_code = "00";
  std::string path3 = "/home/acostama/projects/EPFL/mockturtle/benchmarks/NPN-representatives/nin3/ex"+str_code+".truth";

  muesli_params ps3;
  ps3.max_sup = 2;
  klut_network klut3_0;
  if( lorina::read_truth( path3, truth_reader( klut3_0 ) ) == lorina::return_code::parse_error )
    assert( false );

  auto LFE_pre3 = graph_to_lfe( klut3_0 );
  print_LFE( LFE_pre3, true );

  print_mutual_informations3( LFE_pre3 );
  auto x = LFE_pre3.partial.first[0];
  auto y = LFE_pre3.partial.first[1];
  auto z = LFE_pre3.partial.first[2];
  dbitset_vector g {dbitset (8,0x88)};
  auto f = LFE_pre3.partial.second;

  std::cout << "I(g=xy;f)=" << kitty::mutual_information( g, f ) << "<-" << g[0] << std::endl;
  std::cout << "I(g,x;f)=" << kitty::mutual_information( std::vector{g[0], x}, f ) << std::endl;
  std::cout << "I(g,y;f)=" << kitty::mutual_information( std::vector{g[0], y}, f ) << std::endl;
  std::cout << "I(g,x,y;f)=" << kitty::mutual_information( std::vector{g[0], x, y}, f ) << std::endl;
  std::cout << "I(x,y;f)=" << kitty::mutual_information( std::vector{ x, y}, f ) << std::endl;


  std::cout << "H(f)= " << kitty::entropy( std::vector{LFE_pre3.partial.second} ) << std::endl;


  std::cout << "#######################################################" <<std::endl;
  std::cout << "                        f = a(b^c)                     " <<std::endl;
  std::cout << "#######################################################" <<std::endl;
  str_code = "01";
  path3 = "/home/acostama/projects/EPFL/mockturtle/benchmarks/NPN-representatives/nin3/ex"+str_code+".truth";

  klut_network klut3_1;
  if( lorina::read_truth( path3, truth_reader( klut3_1 ) ) == lorina::return_code::parse_error )
    assert( false );

  LFE_pre3 = graph_to_lfe( klut3_1 );
  print_LFE( LFE_pre3, true );

  print_mutual_informations3( LFE_pre3 );
  std::cout << "H(f)= " << kitty::entropy( std::vector{LFE_pre3.partial.second} ) << std::endl;

  std::cout << "#######################################################" <<std::endl;
  std::cout << "                        f = a(b+c)                     " <<std::endl;
  std::cout << "#######################################################" <<std::endl;
  str_code = "02";
  path3 = "/home/acostama/projects/EPFL/mockturtle/benchmarks/NPN-representatives/nin3/ex"+str_code+".truth";

  klut_network klut3_2;
  if( lorina::read_truth( path3, truth_reader( klut3_2 ) ) == lorina::return_code::parse_error )
    assert( false );

  LFE_pre3 = graph_to_lfe( klut3_2 );
  print_LFE( LFE_pre3, true );

  print_mutual_informations3( LFE_pre3 );
  std::cout << "H(f)= " << kitty::entropy( std::vector{LFE_pre3.partial.second} ) << std::endl;

  std::cout << "#######################################################" <<std::endl;
  std::cout << "        f = ( ab'c' )^( a'bc' )^( a'b'c )              " <<std::endl;
  std::cout << "#######################################################" <<std::endl;
  str_code = "03";
  path3 = "/home/acostama/projects/EPFL/mockturtle/benchmarks/NPN-representatives/nin3/ex"+str_code+".truth";

  klut_network klut3_3;
  if( lorina::read_truth( path3, truth_reader( klut3_3 ) ) == lorina::return_code::parse_error )
    assert( false );

  LFE_pre3 = graph_to_lfe( klut3_3 );
  print_LFE( LFE_pre3, true );

  print_mutual_informations3( LFE_pre3 );
  std::cout << "H(f)= " << kitty::entropy( std::vector{LFE_pre3.partial.second} ) << std::endl;

  dbitset_vector f0 {dbitset (8,0x11)}; // x'y'
  dbitset_vector f1 {dbitset (8,0x44)}; // 
  dbitset_vector f2 {dbitset (8,0x22)};
  dbitset_vector f3 {dbitset (8,0x66)};
  dbitset_vector f4 {dbitset (8,0x77)};
  dbitset_vector f5 {dbitset (8,0x02)};

  dbitset_vector f6 {dbitset (8,0x02)};
  dbitset_vector f7 {dbitset (8,0xd2)};
  dbitset_vector f8 {dbitset (8,0xdf)};
  dbitset_vector f9 {dbitset (8,0x04)};
  dbitset_vector f10 {dbitset (8,0xb4)};
  dbitset_vector f11 {dbitset (8,0xbf)};
  dbitset_vector f12 {dbitset (8,0x1e)};
  dbitset_vector f13 {dbitset (8,0x10)};

  f = LFE_pre3.partial.second;

  std::cout << "I(x'y';f)=" << kitty::mutual_information( f0, f ) << "<-" << f0[0] << std::endl;
  std::cout << "I(x'y;f)=" << kitty::mutual_information( f1, f ) << "<-" << f1[0] << std::endl;
  std::cout << "I(xy';f)=" << kitty::mutual_information( f2, f ) << "<-" << f2[0] << std::endl;
  std::cout << "X I(x^y;f)=" << kitty::mutual_information( f3, f ) << "<-" << f3[0] << std::endl;
  std::cout << "X I((xy)';f)=" << kitty::mutual_information( f4, f ) << "<-" << f4[0] << std::endl;
  std::cout << "X I(f3z';f)=" << kitty::mutual_information( f6, f ) << "<-" << f6[0] << std::endl;
  std::cout << "X I(f3^z;f)=" << kitty::mutual_information( f7, f ) << "<-" << f7[0] << std::endl;
  std::cout << "X I((f3z)';f)=" << kitty::mutual_information( f8, f ) << "<-" << f8[0] << std::endl;
  std::cout << "X I(f4z';f)=" << kitty::mutual_information( f9, f ) << "<-" << f9[0] << std::endl;
  std::cout << "X I(z^f4;f)=" << kitty::mutual_information( f10, f ) << "<-" << f10[0] << std::endl;
  std::cout << "X I((zf4)';f)=" << kitty::mutual_information( f11, f ) << "<-" << f11[0] << std::endl;
  std::cout << "X I((z^f5)';f)=" << kitty::mutual_information( f12, f ) << "<-" << f12[0] << std::endl;
  std::cout << "X I(zf5;f)=" << kitty::mutual_information( f13, f ) << "<-" << f13[0] << std::endl;


  std::cout << "#######################################################" <<std::endl;
  std::cout << "                    f = xyz^x'y'z'                     " <<std::endl;
  std::cout << "#######################################################" <<std::endl;
  str_code = "05";
  path3 = "/home/acostama/projects/EPFL/mockturtle/benchmarks/NPN-representatives/nin3/ex"+str_code+".truth";

  klut_network klut3_5;
  if( lorina::read_truth( path3, truth_reader( klut3_5 ) ) == lorina::return_code::parse_error )
    assert( false );

  LFE_pre3 = graph_to_lfe( klut3_5 );
  print_LFE( LFE_pre3, true );

  print_mutual_informations3( LFE_pre3 );
  std::cout << "H(f)= " << kitty::entropy( std::vector{LFE_pre3.partial.second} ) << std::endl;

  f = LFE_pre3.partial.second;

  x = LFE_pre3.partial.first[0];
  y = LFE_pre3.partial.first[1];
  z = LFE_pre3.partial.first[2];

  f0 = {dbitset (8,0x88)};
  f1 = {dbitset (8,0x11)};
  f2 = {dbitset (8,0x99)};
  f3 = {dbitset (8,0x80)};
  f4 = {dbitset (8,0x87)};
  f5 = std::vector{dbitset (8,0x01)};

  std::cout << "I(xy;f)=" << kitty::mutual_information( f0, f ) << "<-" << f0[0] << std::endl;
  std::cout << "I(x'y';f)=" << kitty::mutual_information( f1, f ) << "<-" << f1[0] << std::endl;
  std::cout << "I(xy,x'y';f)=" << kitty::mutual_information( std::vector{f0[0],f1[0]}, f ) << std::endl;
  std::cout << "I(xy,x'y',x;f)=" << kitty::mutual_information( std::vector{f0[0],f1[0],x}, f ) << std::endl;
  std::cout << "I(xy,x'y',y;f)=" << kitty::mutual_information( std::vector{f0[0],f1[0],y}, f ) << std::endl;
  std::cout << "I(xy,x'y',x,y;f)=" << kitty::mutual_information( std::vector{f0[0],f1[0],x,y}, f ) << std::endl;

  std::cout << "I((x^y)';f)=" << kitty::mutual_information( f2, f ) << std::endl;
  std::cout << "I((x^y)',x;f)=" << kitty::mutual_information( std::vector{f2[0], x}, f ) << std::endl;
  std::cout << "I((x^y)',y;f)=" << kitty::mutual_information( std::vector{f2[0], y}, f ) << std::endl;
  std::cout << "I((x^y)',x,y;f)=" << kitty::mutual_information( std::vector{f2[0], x,y}, f ) << std::endl;
  
  
  std::cout << "I(z(xy);f)=" << kitty::mutual_information( f3, f ) << "<-" << f3[0] << std::endl;
  std::cout << "I(z^(xy);f)=" << kitty::mutual_information( f4, f ) << "<-" << f4[0] << std::endl;
  std::cout << "I(z'(x'y');f)=" << kitty::mutual_information( f5, f ) << "<-" << f5[0] << std::endl;

  std::cout << "#######################################################" <<std::endl;
  std::cout << "                     f = x^(z+xy)                      " <<std::endl;
  std::cout << "#######################################################" <<std::endl;
  str_code = "06";
  path3 = "/home/acostama/projects/EPFL/mockturtle/benchmarks/NPN-representatives/nin3/ex"+str_code+".truth";

  klut_network klut3_6;
  if( lorina::read_truth( path3, truth_reader( klut3_6 ) ) == lorina::return_code::parse_error )
    assert( false );

  LFE_pre3 = graph_to_lfe( klut3_6 );
  print_LFE( LFE_pre3, true );

  print_mutual_informations3( LFE_pre3 );
  std::cout << "H(f)= " << kitty::entropy( std::vector{LFE_pre3.partial.second} ) << std::endl;

  f = LFE_pre3.partial.second;

  std::cout << "#######################################################" <<std::endl;
  std::cout << "                      f = xy + x'z                     " <<std::endl;
  std::cout << "#######################################################" <<std::endl;
  str_code = "07";
  path3 = "/home/acostama/projects/EPFL/mockturtle/benchmarks/NPN-representatives/nin3/ex"+str_code+".truth";

  klut_network klut3_7;
  if( lorina::read_truth( path3, truth_reader( klut3_7 ) ) == lorina::return_code::parse_error )
    assert( false );

  LFE_pre3 = graph_to_lfe( klut3_7 );
  print_LFE( LFE_pre3, true );

  print_mutual_informations3( LFE_pre3 );
  std::cout << "H(f)= " << kitty::entropy( std::vector{LFE_pre3.partial.second} ) << std::endl;

  f = LFE_pre3.partial.second;

  std::cout << "#######################################################" <<std::endl;
  std::cout << "                        f = x^(yz)                     " <<std::endl;
  std::cout << "#######################################################" <<std::endl;
  str_code = "08";
  path3 = "/home/acostama/projects/EPFL/mockturtle/benchmarks/NPN-representatives/nin3/ex"+str_code+".truth";

  klut_network klut3_8;
  if( lorina::read_truth( path3, truth_reader( klut3_8 ) ) == lorina::return_code::parse_error )
    assert( false );

  LFE_pre3 = graph_to_lfe( klut3_8 );
  print_LFE( LFE_pre3, true );

  print_mutual_informations3( LFE_pre3 );
  std::cout << "H(f)= " << kitty::entropy( std::vector{LFE_pre3.partial.second} ) << std::endl;

  f = LFE_pre3.partial.second;

  std::cout << "#######################################################" <<std::endl;
  std::cout << "                         f = x^y^z                     " <<std::endl;
  std::cout << "#######################################################" <<std::endl;
  str_code = "09";
  path3 = "/home/acostama/projects/EPFL/mockturtle/benchmarks/NPN-representatives/nin3/ex"+str_code+".truth";

  klut_network klut3_9;
  if( lorina::read_truth( path3, truth_reader( klut3_9 ) ) == lorina::return_code::parse_error )
    assert( false );

  LFE_pre3 = graph_to_lfe( klut3_9 );
  print_LFE( LFE_pre3, true );

  print_mutual_informations3( LFE_pre3 );
  std::cout << "H(f)= " << kitty::entropy( std::vector{LFE_pre3.partial.second} ) << std::endl;

  f = LFE_pre3.partial.second;

  std::cout << "#######################################################" <<std::endl;
  std::cout << "                     f = <x y z>                       " <<std::endl;
  std::cout << "#######################################################" <<std::endl;
  str_code = "04";
  path3 = "/home/acostama/projects/EPFL/mockturtle/benchmarks/NPN-representatives/nin3/ex"+str_code+".truth";

  klut_network klut3_4;
  if( lorina::read_truth( path3, truth_reader( klut3_4 ) ) == lorina::return_code::parse_error )
    assert( false );

  LFE_pre3 = graph_to_lfe( klut3_4 );
  print_LFE( LFE_pre3, true );

  print_mutual_informations3( LFE_pre3 );
  std::cout << "H(f)= " << kitty::entropy( std::vector{LFE_pre3.partial.second} ) << std::endl;

  f = LFE_pre3.partial.second;

  f0 = {dbitset (8,0x88)};
  f1 = {dbitset (8,0xa0)};
  f2 = {dbitset (8,0xc0)};
  f3 = {dbitset (8,0xee)};

  std::cout << "I(xy;f)=" << kitty::mutual_information( f0, f ) << "<-" << f0[0] << std::endl;
  std::cout << "I(xz;f)=" << kitty::mutual_information( f1, f ) << "<-" << f1[0] << std::endl;
  std::cout << "I(yz;f)=" << kitty::mutual_information( f2, f ) << "<-" << f2[0] << std::endl;
  std::cout << "I(x+y;f)=" << kitty::mutual_information( f3, f ) << "<-" << f3[0] << std::endl;

  f = LFE_pre3.partial.second;

  std::cout << "#######################################################" <<std::endl;

  kitty::dynamic_truth_table tt_s(3u);
  kitty::create_from_hex_string( tt_s, "16" );
  kitty::print_binary(tt_s); std::cout << std::endl;
  kitty::swap_inplace( tt_s, 0, 1 );
  kitty::print_binary(tt_s);
  kitty::swap_inplace( tt_s, 1, 2 );
  kitty::print_binary(tt_s);
  //swap_adjacent_inplace( tt_s, GetParam().first );


  return 0;
}