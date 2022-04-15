
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

struct divisors
{
  std::string tt;
  dbitset fn;
};

std::pair<std::string, std::vector<dbitset>> compute_divisors( dbitset const& x1, dbitset const& x2, dbitset const& f )
{
  std::vector<divisors> res; 
  std::vector<dbitset> X1;
  std::vector<dbitset> X0;
  std::vector<size_t> C0;
  std::vector<size_t> C1;


  X1.push_back(~x1 & ~x2 & f);
  X0.push_back(~x1 & ~x2 & ~f);  
  X1.push_back(~x1 & x2 & f);
  X0.push_back(~x1 & x2 & ~f);  
  X1.push_back(x1 & ~x2 & f);
  X0.push_back(x1 & ~x2 & ~f);
  X1.push_back(x1 & x2 & f);
  X0.push_back(x1 & x2 & ~f);

  std::string tt_str = "";
  dbitset mask;
  dbitset values;
  for( auto i=0u; i<4; ++i )
  {
    C0.push_back( X0[i].count() );
    C1.push_back( X1[i].count() );
    if( C0[i] == 0 && C1[i] != 0 )
    {
      tt_str = "1" + tt_str;
      mask.push_back(1u);
      values.push_back(1u);
    }
    else if( C1[i] == 0 && C0[i] != 0 )
    {
      tt_str = "0" + tt_str;
      mask.push_back(1u);
      values.push_back(0u);
    }
    else
    {
      tt_str = "x" + tt_str;
      mask.push_back(0u);
      values.push_back(1u);
    }
  }
  std::vector<uint32_t> D = { 1, 2, 4, 6, 7, 8, 9, 11, 13, 14 };
  std::vector<dbitset> divisors;
  for( auto i : D )
  {
    dbitset new_vect(4u,i);
    if((new_vect & mask) == (values & mask))
    {
      std::cout << new_vect << " A " << std::endl;
      dbitset new_divisor ( 8u, 0 ); 
      if( new_vect[0] == 1 )
        new_divisor |= ~x1 & ~x2;
      if( new_vect[1] == 1 )
        new_divisor |= ~x1 & x2;
      if( new_vect[2] == 1 )
        new_divisor |= x1 & ~x2;
      if( new_vect[3] == 1 )
        new_divisor |= x1 & x2;
      divisors.push_back(new_divisor);
      std::cout << "new divisor " << new_divisor << std::endl;
    }
    else
      std::cout << new_vect << " R" << std::endl;


  }
  return std::make_pair(tt_str, divisors);

}

int main()
{

 

  std::cout << "#######################################################" <<std::endl;
  std::cout << "                     f = <x y z>                       " <<std::endl;
  std::cout << "#######################################################" <<std::endl;
  std::string str_code = "04";
  std::string path_maj = "/home/acostama/projects/EPFL/mockturtle/benchmarks/NPN-representatives/nin3/ex"+str_code+".truth";

  //muesli_params ps3;
  //ps3.max_sup = 2;
  klut_network klut_maj;
  if( lorina::read_truth( path_maj, truth_reader( klut_maj ) ) == lorina::return_code::parse_error )
    assert( false );

  auto LFE_pre_maj = graph_to_lfe( klut_maj );
  print_LFE( LFE_pre_maj, true );

  print_mutual_informations3( LFE_pre_maj );
  std::cout << "H(f)= " << kitty::entropy( std::vector{LFE_pre_maj.partial.second} ) << std::endl;
  auto f = LFE_pre_maj.partial.second;

  auto x = LFE_pre_maj.partial.first[0];
  auto y = LFE_pre_maj.partial.first[1];
  auto z = LFE_pre_maj.partial.first[2];

  std::cout << "divisors: " << std::endl;
  auto itt_fin = compute_divisors( x, y, f );
  std::cout << "itt_str = "<< itt_fin.first << std::endl;
  
  auto g1 = itt_fin.second[0];
  auto g2 = itt_fin.second[1];

  std::cout << "I(g1=xy;f)=" << kitty::mutual_information( g1, f ) << "<-" << g1 << std::endl;
  std::cout << "I(g2=x+y;f)=" << kitty::mutual_information( g2, f ) << "<-" << g2 << std::endl;
  std::cout << "\nCan g1 alone give us f?"  << std::endl;
  std::cout << "I(g1=xy;f)=" << kitty::mutual_information( g1, f ) << std::endl;
  std::cout << "I(g1=xy,x;f)=" << kitty::mutual_information( std::vector{g1,x}, f ) << std::endl;
  std::cout << "I(g1=xy,y;f)=" << kitty::mutual_information( std::vector{g1,y}, f ) << std::endl;
  std::cout << "I(g1=xy,x,y;f)=" << kitty::mutual_information( std::vector{g1,x,y}, f ) << std::endl;
  std::cout << " NO " << std::endl;

  std::cout << "\nCan g2 alone give us f?"  << std::endl;
  std::cout << "I(g2=x+y;f)=" << kitty::mutual_information( g2, f ) << std::endl;
  std::cout << "I(g2=x+y,x;f)=" << kitty::mutual_information( std::vector{g2,x}, f ) << std::endl;
  std::cout << "I(g2=x+y,y;f)=" << kitty::mutual_information( std::vector{g2,y}, f ) << std::endl;
  std::cout << "I(g2=x+y,x,y;f)=" << kitty::mutual_information( std::vector{g2,x,y}, f ) << std::endl;
  std::cout << " NO " << std::endl;

  std::cout << "\nCan g1 and g2 substitute x and y?"<< std::endl;
  std::cout << "I(g1=xy,g2=x+y;f)=" << kitty::mutual_information( std::vector{g1,g2}, f ) << std::endl;
  std::cout << "I(g1=xy,g2=x+y,x;f)=" << kitty::mutual_information( std::vector{g1,g2,x}, f ) << std::endl;
  std::cout << "I(g1=xy,g2=x+y,y;f)=" << kitty::mutual_information( std::vector{g1,g2,y}, f ) << std::endl;
  std::cout << "I(g1=xy,g2=x+y,x,y;f)=" << kitty::mutual_information( std::vector{g1,g2,x,y}, f ) << std::endl;
  std::cout << " YES " << std::endl;
  std::cout << " in the list we have : [ z, g1, g2 ] " << std::endl;
  std::cout << " due to the symmetry of the vatriables we first have to try to assemble the functions f(z,g1), f(z,g2) " << std::endl;
  std::cout << " f(z,g1) " << std::endl;
  std::cout << "divisors: " << std::endl;
  auto itt2_fin = compute_divisors( z, g1, f );
  std::cout << "itt2_str = "<< itt2_fin.first << std::endl;
  auto g3 = itt2_fin.second[0];
  std::cout << "g3 = z+g1 = "<< g3 << std::endl;
  std::cout << "I(g3=g1+z;f)=" << kitty::mutual_information( g3, f )<< "<-" << g3 << std::endl;
  std::cout << std::endl;
  std::cout << " f(z,g2) " << std::endl;
  std::cout << "divisors: " << std::endl;
  auto itt3_fin = compute_divisors( z, g2, f );
  std::cout << "itt3_str = "<< itt3_fin.first << std::endl;
  auto g4 = itt3_fin.second[0];
  std::cout << "g4 = zg2 = "<< g4 << std::endl;
  std::cout << "I(g4=g2z;f)=" << kitty::mutual_information( g4, f )<< "<-" << g4 << std::endl;
  std::cout << std::endl;
  std::cout << "list = [z,g1,g2,g3,g4]" << std::endl;
  std::cout << "new functions introduced from z,g1 and z,g2" << std::endl;
  std::cout << "is any info redundant?" << std::endl;
  std::cout << "I(g3,z;f)=" << kitty::mutual_information( std::vector{g3,z}, f )<< "=?="<< "I(g3;f)=" << kitty::mutual_information( g3, f ) << std::endl;
  std::cout << "I(g3,g1;f)=" << kitty::mutual_information( std::vector{g3,g1}, f )<< "=?="<< "I(g1;f)=" << kitty::mutual_information( g1, f ) << std::endl;
  std::cout << "I(g4,z;f)=" << kitty::mutual_information( std::vector{g4,z}, f )<< "=?="<< "I(g4;f)=" << kitty::mutual_information( g4, f ) << std::endl;
  std::cout << "I(g4,g2;f)=" << kitty::mutual_information( std::vector{g4,g2}, f )<< "=?="<< "I(g2;f)=" << kitty::mutual_information( g2, f ) << std::endl;
  std::cout << "no support reduction seems legal. Must consider f(g3,z), f(g3,g1), f(g3,g2), f(g4,z), f(g4,g1), f(g4,g2), f(g3, g4)" << std::endl;
  std::cout << "I(g3,z;f)=" << kitty::mutual_information( std::vector{g3,z}, f )<< std::endl;
  std::cout << "I(g3,g1;f)=" << kitty::mutual_information( std::vector{g3,g1}, f )<< std::endl;
  std::cout << "I(g3,g2;f)=" << kitty::mutual_information( std::vector{g3,g2}, f )<< std::endl;
  std::cout << "I(g4,z;f)=" << kitty::mutual_information( std::vector{g4,z}, f )<< std::endl;
  std::cout << "I(g4,g1;f)=" << kitty::mutual_information( std::vector{g4,g1}, f )<< std::endl;
  std::cout << "I(g4,g2;f)=" << kitty::mutual_information( std::vector{g4,g2}, f )<< std::endl;
  std::cout << "I(g3,g4;f)=" << kitty::mutual_information( std::vector{g3,g4}, f )<< std::endl;
  std::cout << "TWO COMBINATIONS HAVE THE SAME ENTROPY AS THE FUNCTION ITSELF!: g3,g2 and g4,g1"<<std::endl;
  std::cout << "g3,g2"<<std::endl;
  auto itt4_fin = compute_divisors( g2, g3, f );
  std::cout << "itt4_str = "<< itt4_fin.first << std::endl;
  std::cout << "FOUND! f= g2 & g3 = ( x + y ) & ( z + xy) "<<std::endl;
  std::cout << std::endl;
  std::cout << "g4,g1"<<std::endl;
  auto itt5_fin = compute_divisors( g4, g1, f );
  std::cout << "itt5_str = "<< itt5_fin.first << std::endl;
  std::cout << "FOUND! f= g4 + g1 =  z( x + y ) + xy "<<std::endl;
  std::cout << "#######################################################" <<std::endl;
  std::cout << "                    f = xyz^x'y'z'                     " <<std::endl;
  std::cout << "#######################################################" <<std::endl;
  str_code = "05";
  std::string path_gamble = "/home/acostama/projects/EPFL/mockturtle/benchmarks/NPN-representatives/nin3/ex"+str_code+".truth";
 
   klut_network klut_gamble;
  if( lorina::read_truth( path_gamble, truth_reader( klut_gamble ) ) == lorina::return_code::parse_error )
    assert( false );

  auto LFE_pre_gamble = graph_to_lfe( klut_gamble );
  print_LFE( LFE_pre_gamble, true );

  print_mutual_informations3( LFE_pre_gamble );
  std::cout << "H(f)= " << kitty::entropy( std::vector{LFE_pre_gamble.partial.second} ) << std::endl;
  f = LFE_pre_gamble.partial.second;

  x = LFE_pre_gamble.partial.first[0];
  y = LFE_pre_gamble.partial.first[1];
  z = LFE_pre_gamble.partial.first[2];

  std::cout << "Pick x,y -> divisors: " << std::endl;
  auto itt1_fin = compute_divisors( x, y, f );
  std::cout << "itt1_str = "<< itt1_fin.first << std::endl;
  g1 = itt1_fin.second[0];
  g2 = itt1_fin.second[1];
  g3 = itt1_fin.second[2];

  std::cout << "I(g1=x'y';f)=" << kitty::mutual_information( g1, f ) << "<-" << g1 << std::endl;
  std::cout << "I(g2=xy;f)=" << kitty::mutual_information( g2, f ) << "<-" << g2 << std::endl;
  std::cout << "I(g3=(x^y)';f)=" << kitty::mutual_information( g3, f ) << "<-" << g3 << std::endl;
  std::cout << "we know they are symmetric so we can divide in 2 symmetry groups {g1,g2} and {g3}" << std::endl;
  std::cout << "Let us start from the bigger in size: {g1,g2}. We must compute the functions with z due to symmetry, but before we check if we can remove x and y" << std::endl;
  
  std::cout << "\nCan g1 alone give us f?"  << std::endl;
  std::cout << "I(g1;f)=" << kitty::mutual_information( g1, f ) << std::endl;
  std::cout << "I(g1,x;f)=" << kitty::mutual_information( std::vector{g1,x}, f ) << std::endl;
  std::cout << "I(g1,y;f)=" << kitty::mutual_information( std::vector{g1,y}, f ) << std::endl;
  std::cout << "I(g1,x,y;f)=" << kitty::mutual_information( std::vector{g1,x,y}, f ) << std::endl;
  std::cout << " NO " << std::endl;

  std::cout << "\nCan g2 alone give us f?"  << std::endl;
  std::cout << "I(g2;f)=" << kitty::mutual_information( g2, f ) << std::endl;
  std::cout << "I(g2,x;f)=" << kitty::mutual_information( std::vector{g2,x}, f ) << std::endl;
  std::cout << "I(g2,y;f)=" << kitty::mutual_information( std::vector{g2,y}, f ) << std::endl;
  std::cout << "I(g2,x,y;f)=" << kitty::mutual_information( std::vector{g2,x,y}, f ) << std::endl;
  std::cout << " NO " << std::endl;

  std::cout << "\nCan g1 and g2 substitute x and y?"<< std::endl;
  std::cout << "I(g1,g2;f)=" << kitty::mutual_information( std::vector{g1,g2}, f ) << std::endl;
  std::cout << "I(g1,g2,x;f)=" << kitty::mutual_information( std::vector{g1,g2,x}, f ) << std::endl;
  std::cout << "I(g1,g2,y;f)=" << kitty::mutual_information( std::vector{g1,g2,y}, f ) << std::endl;
  std::cout << "I(g1,g2,x,y;f)=" << kitty::mutual_information( std::vector{g1,g2,x,y}, f ) << std::endl;
  std::cout << " YES " << std::endl;
  
  std::cout << " in the list we have : [ z, g1, g2 ] " << std::endl;
  std::cout << " due to the symmetry of the vatriables we first have to try to assemble the functions f(z,g1), f(z,g2) " << std::endl;
  std::cout << " f(z,g1) " << std::endl;
  std::cout << "divisors: " << std::endl;
  itt2_fin = compute_divisors( z, g1, f );
  std::cout << "itt2_str = "<< itt2_fin.first << std::endl;
  g4 = itt2_fin.second[0];
  auto g5 = itt2_fin.second[1];
  std::cout << " f(z,g2) " << std::endl;
  std::cout << "divisors: " << std::endl;
  itt3_fin = compute_divisors( z, g2, f );
  std::cout << "itt3_str = "<< itt3_fin.first << std::endl;
  auto g6 = itt3_fin.second[0];
  auto g7 = itt3_fin.second[1];
  std::cout << "I(g4=z'g1 = z'(x'y');f)=" << kitty::mutual_information( g4, f )<< "<-" << g4 << std::endl;
  std::cout << "I(g5=z^g1 = z^(x'y');f)=" << kitty::mutual_information( g5, f )<< "<-" << g5 << std::endl;
  std::cout << "I(g6=zg2 = z^(x'y');f)=" << kitty::mutual_information( g6, f )<< "<-" << g6 << std::endl;
  std::cout << "I(g7=(z^g2)' = z'^(xy)';f)=" << kitty::mutual_information( g7, f )<< "<-" << g7 << std::endl;
  std::cout << "list [g5,g7,z,g1,g2] [g4,g6,z,g1,g2]" << std::endl;

  std::cout << "\nCan g5 and g7 give us f, removing g1, g2 and z?"  << std::endl;
  std::cout << "I(g5,g7;f)=" << kitty::mutual_information( std::vector{g5,g7}, f ) << std::endl;
  std::cout << "I(g5,g7,g1;f)=" << kitty::mutual_information( std::vector{g5,g7,g1}, f ) << std::endl;
  std::cout << "I(g5,g7,g2;f)=" << kitty::mutual_information( std::vector{g5,g7,g2}, f ) << std::endl;
  std::cout << "I(g5,g7,z;f)=" << kitty::mutual_information( std::vector{g5,g7,z}, f ) << std::endl;
  std::cout << "YES => [g5,g7]" << std::endl;
  std::cout << "divisors: " << std::endl;
  itt4_fin = compute_divisors( g5, g7, f );
  std::cout << "itt4_str = "<< itt4_fin.first << std::endl;
  std::cout << "FOUND! f= g5 & g7 = ( z^g1 ) & ( (z^g2)' ) = ( z^(x'y') ) & ( (z^(xy))' ) "<<std::endl;


  std::cout << "\nCan g4 and g6 give us f, removing g1, g2 and z?"  << std::endl;
  std::cout << "I(g4,g6;f)=" << kitty::mutual_information( std::vector{g4,g6}, f ) << std::endl;
  std::cout << "I(g4,g6,g1;f)=" << kitty::mutual_information( std::vector{g4,g6,g1}, f ) << std::endl;
  std::cout << "I(g4,g6,g2;f)=" << kitty::mutual_information( std::vector{g4,g6,g2}, f ) << std::endl;
  std::cout << "I(g4,g6,z;f)=" << kitty::mutual_information( std::vector{g4,g6,z}, f ) << std::endl;
  std::cout << "YES => [g4,g6]" << std::endl;
  std::cout << "divisors: " << std::endl;
  itt5_fin = compute_divisors( g4, g6, f );
  std::cout << "itt5_str = "<< itt5_fin.first << std::endl;
  std::cout << "FOUND! f= g4 ^ g6 =  z'(x'y')^z(xy) "<<std::endl;
  std::cout << "FOUND! f= g4 + g6 =  z'(x'y')+z(xy) "<<std::endl;
  std::cout << std::endl;
  std::cout << "let's go back to g3 = (x^y)'. Doe it do anything useful with z?" << std::endl;
  std::cout << "divisors: " << std::endl;
  auto itt6_fin = compute_divisors( g3, z, f );
  std::cout << "itt6_str = "<< itt6_fin.first << std::endl;
  std::cout << "UNDERSTAND HOW TO FORMALIZE THE NEGLECTANCE OF THE FIRST 10011001: is this because it is included in g1 g2" <<std::endl;

  std::cout << "#######################################################" <<std::endl;
  std::cout << "        f = ( xy'z' )^( x'yz' )^( x'y'z )              " <<std::endl;
  std::cout << "#######################################################" <<std::endl;
  str_code = "03";
  auto path_1hot = "/home/acostama/projects/EPFL/mockturtle/benchmarks/NPN-representatives/nin3/ex"+str_code+".truth";

  klut_network klut_1hot;
  if( lorina::read_truth( path_1hot, truth_reader( klut_1hot ) ) == lorina::return_code::parse_error )
    assert( false );

  auto LFE_pre_1hot = graph_to_lfe( klut_1hot );
  print_LFE( LFE_pre_1hot, true );

  print_mutual_informations3( LFE_pre_1hot );
  std::cout << "H(f)= " << kitty::entropy( std::vector{LFE_pre_1hot.partial.second} ) << std::endl;
  f = LFE_pre_1hot.partial.second;

  x = LFE_pre_1hot.partial.first[0];
  y = LFE_pre_1hot.partial.first[1];
  z = LFE_pre_1hot.partial.first[2];

  std::cout << "divisors: " << std::endl;
  auto ST1hot = compute_divisors( x, y, f );
  std::cout << "ST1hot = "<< ST1hot.first << std::endl;
  auto w1 = ST1hot.second[0];
  auto w2 = ST1hot.second[1];
  auto w3 = ST1hot.second[2];
  auto w4 = ST1hot.second[3];
  auto w5 = ST1hot.second[4];
  std::cout << "I(w1=x'y';f)=" << kitty::mutual_information( w1, f ) << "<-" << w1 << std::endl;
  std::cout << "I(w2=x'y;f)=" << kitty::mutual_information( w2, f ) << "<-" << w2 << std::endl;
  std::cout << "I(w3=xy';f)=" << kitty::mutual_information( w3, f ) << "<-" << w3 << std::endl;
  std::cout << "I(w4=x^y;f)=" << kitty::mutual_information( w4, f ) << "<-" << w4 << std::endl;
  std::cout << "I(w5=(xy)';f)=" << kitty::mutual_information( w5, f ) << "<-" << w5 << std::endl;
  std::cout << "a basis is w1, w2, w3" << std::endl;
  std::cout << "w1, z" << std::endl;
  std::cout << "divisors: " << std::endl;
  ST1hot = compute_divisors( w1, z, f );
  std::cout << "ST1hot = "<< ST1hot.first << std::endl;
  auto w6 = ST1hot.second[0];
  std::cout << "I(w6=z&w1=z(x'y');f)=" << kitty::mutual_information( w6, f ) << "<-" << w6 << std::endl;
  auto w7 = ST1hot.second[1];
  std::cout << "I(w7=z&w1=(z^(x'y'))';f)=" << kitty::mutual_information( w7, f ) << "<-" << w7 << std::endl;
  std::cout << "w2, z" << std::endl;
  std::cout << "divisors: " << std::endl;
  ST1hot = compute_divisors( w2, z, f );
  std::cout << "ST1hot = "<< ST1hot.first << std::endl;
  auto w8 = ST1hot.second[0];
  std::cout << "I(w8=w2&z';f)=" << kitty::mutual_information( w8, f ) << "<-" << w8 << std::endl;
  auto w9 = ST1hot.second[1];
  std::cout << "I(w9=w2^z;f)=" << kitty::mutual_information( w9, f ) << "<-" << w9 << std::endl;
  auto w10 = ST1hot.second[2];
  std::cout << "I(w10=(w2&z)';f)=" << kitty::mutual_information( w10, f ) << "<-" << w10 << std::endl;

    std::cout << "w3, z" << std::endl;
  std::cout << "divisors: " << std::endl;
  ST1hot = compute_divisors( w3, z, f );
  std::cout << "ST1hot = "<< ST1hot.first << std::endl;
  auto w11 = ST1hot.second[0];
  std::cout << "I(w11=w3z';f)=" << kitty::mutual_information( w11, f ) << "<-" << w11 << std::endl;
  auto w12 = ST1hot.second[1];
  std::cout << "I(w12=w3^z;f)=" << kitty::mutual_information( w12, f ) << "<-" << w12 << std::endl;
  auto w13 = ST1hot.second[2];
  std::cout << "I(w13=(w3z)';f)=" << kitty::mutual_information( w13, f ) << "<-" << w13 << std::endl;
  std::cout << "we have {w6,w8,w11}"<<std::endl;
  
  std::cout << "w6, w8" << std::endl;
  std::cout << "divisors: " << std::endl;
  ST1hot = compute_divisors( w6, w8, f );
  std::cout << "ST1hot = "<< ST1hot.first << std::endl;
  auto w14 = ST1hot.second[0];
  auto w15 = ST1hot.second[1];
  std::cout << "I(w14=w6^w8;f)=" << kitty::mutual_information( w14, f ) << "<-" << w14 << std::endl;
  std::cout << "I(w15=w6+w8;f)=" << kitty::mutual_information( w15, f ) << "<-" << w15 << std::endl;
  std::cout << "w11 to go: w11 w14" << std::endl;
  std::cout << "divisors: " << std::endl;
  ST1hot = compute_divisors( w11, w14, f );
  std::cout << "ST1hot = "<< ST1hot.first << std::endl;
  return 0;
}