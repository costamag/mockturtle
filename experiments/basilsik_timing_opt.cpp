#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/fanout_view.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <kitty/spectral.hpp>
#include <mockturtle/networks/sequential.hpp>
#include "experiments.hpp"
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/functional_reduction.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/balancing.hpp>
#include <mockturtle/algorithms/balancing/sop_balancing.hpp>
#include <mockturtle/algorithms/cut_rewriting.hpp>
#include <mockturtle/algorithms/mapper.hpp>
#include <mockturtle/algorithms/mig_resub.hpp>
#include <mockturtle/algorithms/sim_resub.hpp>
#include <mockturtle/algorithms/resubstitution.hpp>
#include <mockturtle/algorithms/xag_algebraic_rewriting.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/utils/sequential_converter.hpp>
#include <thread>
#include <mutex>
#include <algorithm>
#include <set>
#include <random>

#include <iostream>

using namespace mockturtle;
using namespace experiments;

uint32_t CALL{0};
//using experiment_t = experiments::experiment<std::string, uint32_t, uint32_t, std::string>;
//experiment_t exp_res( "basilsik", "benchmark", "#gates", "depth", "method" );

bool read_file( sequential<aig_network>&, std::string const& );
void print_commands();

template<class Ntk>
void print_stats( Ntk const& );

template<class Ntk>
bool optimize( Ntk &, std::string const& );

template<class Ntk>
void print_stats( binding_view<Ntk> const& );

template<class Ntk>
Ntk abc_opto( Ntk const&, std::string, std::string );

int main( int argc, char* argv[] )
{

  if ( argc != 2 )
  {
    printf("[e] argc should be 2\n");
    return 1;
  }
  std::string istring( argv[1] );
  std::string benchmark = "../experiments/benchmarks/" + istring + ".aig";

  sequential<aig_network> saig;
  read_file( saig, benchmark );

  network_converters_stats st;
  aig_network aig = sequential_to_combinatorial( saig, st );

  print_commands();

  std::string command;
  while( true )
  {
    getline( std::cin, command );
    if( command == "map" )
      break;
    if( !optimize( aig, command ) )
      printf("wrong command\n");
    
    print_stats( aig );
  }

  sequential<aig_network> saig2 = combinatorial_to_sequential( aig, st );

//  /* library to map to technology */
//  std::vector<gate> gates;
//  std::ifstream in( cell_libraries_path( "mcnc" ) );
//
//  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
//  {
//    return 1;
//  }
//
//  tech_library_params tps;
//  tech_library<5, classification_type::np_configurations> tech_lib( gates, tps );
//  map_params ps2;
//  ps2.cut_enumeration_ps.minimize_truth_table = true;
//  ps2.cut_enumeration_ps.cut_limit = 24;
//  map_stats st2;
//
//  binding_view<sequential<klut_network>> res2 = seq_map( dsaig, tech_lib, ps2, &st2 );
//
//  print_stats( res2 );

  return 0;
}


#pragma region parsing

bool read_file( sequential<aig_network>& saig, std::string const& path )
{
  if( lorina::read_aiger( path, aiger_reader( saig ) ) != lorina::return_code::success )
  {
    std::cerr << "read_aiger failed" << std::endl;
    return false;
  }
  return 1;
}

void print_commands()
{
  printf("===============================\n");
  printf("map         : map to technology\n");
  printf("abc-<script>: abc script\n");
  printf("===============================\n");
}

template<class Ntk>
void print_stats( Ntk const& ntk )
{
  depth_view<Ntk> dntk{ntk};
  printf("#gates = %5d #levels = %5d\n", dntk.num_gates(), dntk.depth() );
}

template<class Ntk>
void print_stats( binding_view<Ntk> const& bntk )
{
  printf("area = %.2f delay = %.2f\n", bntk.compute_area(), bntk.compute_worst_delay() );
}

#pragma endregion parsing

#pragma region optimization

template<class Ntk>
bool optimize( Ntk& ntk, std::string const& cmd )
{
  if( cmd.substr(0, 3) == "abc" )
  {
    std::string script = cmd.substr( 4, cmd.size()-1 );
    auto tmp = abc_opto( ntk, cmd+std::to_string(CALL++), script );
    ntk = tmp;
    return true;
  }
  return false;
}

template<class Ntk>
Ntk abc_opto( Ntk const& ntk, std::string str_code, std::string abc_script )
{
  write_aiger( ntk, "/tmp/"+str_code+".aig" );
  std::string command = "abc -q \"r /tmp/"+str_code+".aig; " + abc_script + "; write_aiger /tmp/" + str_code + ".aig\"";

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
  if( lorina::read_aiger( string_path, aiger_reader( res ) ) != lorina::return_code::success )
    std::cerr << "read_aiger failed" << std::endl;

  return res;
}
#pragma endregion optimization