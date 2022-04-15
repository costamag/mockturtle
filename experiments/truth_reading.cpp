
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <string>
#include <vector>
#include <bitset>
#include <mockturtle/algorithms/cut_rewriting.hpp>
#include <mockturtle/algorithms/aig_resub.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/aiger_reader.hpp>

#include <cmath>
#include <random>
#include <boost/dynamic_bitset.hpp>

#include <mockturtle/algorithms/akers_synthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/resubstitution.hpp>

#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/io/truth_reader.hpp>

#include <mockturtle/algorithms/it_decomposition.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>

#include <lorina/truth.hpp>
#include <mockturtle/views/names_view.hpp>
#include <mockturtle/algorithms/graph_to_lfe.hpp>
#include <typeinfo>
#include <lorina/lorina.hpp>
#include <mockturtle/views/depth_view.hpp>

#include <kitty/constructors.hpp>
#include <kitty/print.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/statistical_bit_operations.hpp>
#include <fstream>
#include <string>
#include <unistd.h>
#include <omp.h>
using namespace mockturtle;

template<class Ntk>
Ntk abc_opto( Ntk const& ntk, std::string str_code )
{
  write_aiger( ntk, "/tmp/"+str_code+".aig" );
  std::string command = "abc -q \"r /tmp/"+str_code+".aig; resyn2rs; write_aiger /tmp/"+str_code+".aig\"";

  std::array<char, 128> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }

  Ntk res;
  std::string string_path = ("/tmp/"+str_code+".aig");
  lorina::read_aiger( string_path, aiger_reader( res ) );  

  return res;
}

bool cec( std::string path_to_ntk1, std::string path_to_ntk2 )
{
  std::string command = "abc -q \"read "+path_to_ntk1+"; &get; &cec "+path_to_ntk2+" \" ";

  std::array<char, 128> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }

  /* search for one line which says "Networks are equivalent" and ignore all other debug output from ABC */
  std::stringstream ss( result );
  std::string line;
  while ( std::getline( ss, line, '\n' ) )
  {
    if ( line.size() >= 23u && line.substr( 0u, 23u ) == "Networks are equivalent" )
    {
      return true;
    }
  }

  return false;
}

void print_LFE( auto LFE, bool only_complete = false )  
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
int main()
{
std::cout << "NUM THREADS = " << omp_get_max_threads() << std::endl;
omp_set_num_threads( 6 );

std::vector<size_t> bvect = { 65,66,67, 40,45,48};

#pragma omp parallel for 
for (uint32_t i = 39 ; i<100; i++) { // bvect.size()
  uint32_t bsk = i;

  std::string str_code;
  if( bsk < 10 )
    str_code = "0"+std::to_string(bsk);
  else
    str_code = std::to_string(bsk);

  bool is_verbose = true;

  std::string path = "/home/acostama/projects/EPFL/mockturtle/benchmarks/iwls2022/ex"+str_code+".truth";

  it_decomposition_params ps;
  ps.max_sup = 4;
  ps.is_informed = true;
  ps.try_top_decomposition = true;
  ps.try_bottom_decomposition = true;
  ps.try_xor_decomposition = true;
  ps.is_trivial = true;
  ps.is_bottom_exact = true;
  ps.use_cumsum = true;

  klut_network klut;

  if( lorina::read_truth( path, truth_reader( klut ) ) == lorina::return_code::parse_error )
    assert( false );

  if( is_verbose )
  {
    std::cout << "TRUTH Ntk before" << std::endl;
    std::cout << "num gates " << klut.num_gates() << std::endl;
    std::cout << "num outputs " << klut.num_pos() << std::endl;
  }
  auto LFE_pre = graph_to_lfe( klut );


  it_decomposition( klut, ps );

  if( is_verbose )
  {
    std::cout << "TRUTH Ntk after" << std::endl;
    std::cout << "num gates " << klut.num_gates() << std::endl;
    std::cout << "num outputs " << klut.num_pos() << std::endl;
  }

  auto aig = convert_klut_to_graph<aig_network>( klut );
  aig = cleanup_dangling( aig );

  xag_npn_resynthesis<aig_network> resyn;
  cut_rewriting_params ps_cr;
  ps_cr.cut_enumeration_ps.cut_size = 4;
  aig = cut_rewriting( aig, resyn, ps_cr );

  using view_t = depth_view<fanout_view<aig_network>>;
  fanout_view<aig_network> fanout_view{aig};
  view_t resub_view{fanout_view};

  aig_resubstitution( resub_view );
  aig = cleanup_dangling( aig );

  auto LFE_after = graph_to_lfe( klut );
  bool error = false;
  if( LFE_pre.complete.second != LFE_after.complete.second )
  {
    error = true;
    std::cerr << "[e] not equivalent according to simple simulation check" << std::endl;
  }

  aig = cleanup_dangling( aig );
  std::string my_path = "/home/acostama/projects/EPFL/mockturtle/simulations/iwls22/resub/aig/"+str_code+".aig";
  std::string best_path = "/home/acostama/projects/EPFL/mockturtle/simulations/iwls22/contest_results/ex"+str_code+".aig";

  depth_view_params psD;
  psD.count_complements = true;
  depth_view depth_aig_old{aig, {}, psD};
  auto new_depth = depth_aig_old.depth();
  auto old_depth = new_depth;
  auto new_num_gates = aig.num_gates();
  auto old_num_gates = new_num_gates;
  do{
    old_depth = new_depth;
    old_num_gates = new_num_gates;
    aig = abc_opto( aig, str_code );
    write_aiger( aig, my_path );
    new_num_gates = aig.num_gates();
    depth_view depth_aig{aig, {}, psD};
    new_depth = depth_aig.depth();
    std::cout << "ng " << new_num_gates << std::endl;
    std::cout << "dp " << new_depth << std::endl;
  } while( (new_num_gates < old_num_gates) || (new_depth < old_depth) );
  
  depth_view depth_aig{aig, {}, psD};

  if( cec( my_path, best_path ) == false )
    std::cerr << "[e] not equivalent according to abc" << std::endl;
  else
    std::cout << "equivalent according to abc" << std::endl;


  std::cout <<".b " << str_code << std::endl;
  std::cout << ".g " << depth_aig.num_gates() << std::endl;
  std::cout << ".s " << depth_aig.size() << std::endl; 
  std::cout << ".d " << depth_aig.depth() << std::endl;
  
  std::ofstream myfile;
  myfile.open ( ("/home/acostama/projects/EPFL/mockturtle/simulations/iwls22/resub/"+str_code+".txt") ); 

  myfile <<".b " << str_code << std::endl;
  if( error )
    myfile<< ".e 1" << std::endl;
  myfile << ".g " << depth_aig.num_gates() << std::endl;
  myfile << ".s " << depth_aig.size() << std::endl; 
  myfile << ".d " << depth_aig.depth() << std::endl;

  myfile.close();

}
  return 0;
}
