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
#include <mockturtle/algorithms/mig_algebraic_rewriting.hpp>
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
//experiment_t exp_res( "basilsik", "benchmark", "gates", "depth", "method" );

bool read_file( sequential<aig_network>&, std::string const& );
void print_aig_commands();
void print_mig_commands();
void print_end_commands();

template<class Ntk>
void print_stats( Ntk const&, bool dot = false );

bool optimize_aig( aig_network &, std::string const& );
bool optimize_mig( depth_view<mig_network> &, std::string const& );

template<class Ntk>
void print_stats( binding_view<Ntk> const&, bool dot = false );

template<class Ntk>
mig_network aig_to_mig( Ntk const& );

template<class Ntk>
aig_network mig_to_aig( Ntk const&, double required_time = std::numeric_limits<double>::max() );

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

  while( true )
  {
    print_aig_commands();

    print_stats( aig, true );

    std::string command;
    while( true )
    {
      getline( std::cin, command );
      if( command == "aig->mig" || command == "exit" )
        break;
      if( !optimize_aig( aig, command ) )
        printf("wrong command\n");
      
      print_stats( aig, true );
    }
    if( command == "exit" )
      break;

    print_mig_commands();

    mig_network mig = aig_to_mig( aig );
    print_stats( mig, true );
    depth_view depth_mig{mig};

    while( true )
    {
      getline( std::cin, command );
      if( command == "mig->aig"  )
        break;
      if( !optimize_mig( depth_mig, command ) )
        printf("wrong command\n");

      print_stats( mig, true );
    }

    mig = cleanup_dangling(mig);
    aig = mig_to_aig( mig, 0 );

    print_end_commands();
    getline( std::cin, command );
    if( command == "exit" )
      break;
  }

  printf("=================================\n");
  printf("            AIG                  \n");
  printf("---------------------------------\n");

  print_stats( aig, true );

  sequential<aig_network> saig2 = combinatorial_to_sequential( aig, st );
  //print_stats( saig2 );

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

void print_aig_commands()
{
  printf("=================================\n");
  printf("            AIG                  \n");
  printf("---------------------------------\n");
  printf("resyn2rs    |                    \n");
  printf("it-resyn2rs | (resyn2rs)^infty   \n");
  printf("compress2rs |                    \n");
  printf("a-map       | if -a; fraig;      \n");
  printf("d-map       | if -g; fraig;      \n");
  printf("lazy        | if -y -K 6; fraig; \n");
  printf("---------------------------------\n");
  printf("exit                             \n");
  printf("aig->mig    | map to MIGs        \n");
  printf("---------------------------------\n\n");

}

void print_mig_commands()
{
  printf("=================================\n");
  printf("            MIG                  \n");
  printf("---------------------------------\n");
  printf("dfs                              \n");
  printf("selective                        \n");
  printf("aggressive                       \n");
  printf("---------------------------------\n");
  printf("mig->aig    | map to MIGs        \n");
  printf("---------------------------------\n");
}

void print_end_commands()
{
  printf("=================================\n");
  printf("            END                  \n");
  printf("---------------------------------\n");
  printf("restart                          \n");
  printf("exit                             \n");
  printf("---------------------------------\n");
}

template<class Ntk>
void print_stats( Ntk const& ntk, bool dot )
{
  depth_view<Ntk> dntk{ntk};
  if( dot )
    printf("> gates = %5d levels = %5d.\n", dntk.num_gates(), dntk.depth() );
  else
    printf("  gates = %5d levels = %5d\n", dntk.num_gates(), dntk.depth() );

}

template<class Ntk>
void print_stats( binding_view<Ntk> const& bntk, bool dot )
{
  if( dot )
    printf("> area = %.2f delay = %.2f.\n", bntk.compute_area(), bntk.compute_worst_delay() );
  else
    printf("  area = %.2f delay = %.2f\n", bntk.compute_area(), bntk.compute_worst_delay() );
}

#pragma endregion parsing

#pragma region optimization

bool optimize_mig( depth_view<mig_network>& depth_mig, std::string const& cmd )
{
  if( cmd == "dfs" )
  {
    mig_algebraic_depth_rewriting_params ps;
    ps.strategy = mig_algebraic_depth_rewriting_params::dfs;
    mig_algebraic_depth_rewriting( depth_mig, ps );
    return true;
  }
  if( cmd == "aggressive" )
  {
    mig_algebraic_depth_rewriting_params ps;
    ps.strategy = mig_algebraic_depth_rewriting_params::aggressive;
    mig_algebraic_depth_rewriting( depth_mig, ps );
    return true;
  }
  if( cmd == "selective" )
  {
    mig_algebraic_depth_rewriting_params ps;
    ps.strategy = mig_algebraic_depth_rewriting_params::selective;
    mig_algebraic_depth_rewriting( depth_mig, ps );
    return true;
  }
  return false;
}

bool optimize_aig( aig_network& ntk, std::string const& cmd )
{
  if( cmd == "resyn2rs" )
  {
    auto tmp = abc_opto( ntk, cmd+std::to_string(CALL++), "resyn2rs" );
    ntk = tmp;
    return true;
  }
  if( cmd == "it-resyn2rs" )
  {
    int nT = ntk.num_gates()+1;
    while( nT > ntk.num_gates() )
    {
      nT=ntk.num_gates();
      ntk = abc_opto( ntk, cmd+std::to_string(CALL++), "resyn2rs" );
      print_stats(ntk);
    }
    return true;
  }  
  if( cmd == "compress2rs" )
  {
    auto tmp = abc_opto( ntk, cmd+std::to_string(CALL++), "compress2rs" );
    ntk = tmp;
    return true;
  }
  if( cmd == "d-map" )
  {
    auto tmp = abc_opto( ntk, cmd+std::to_string(CALL++), "if -g; fraig;" );
    ntk = tmp;
    return true;
  }
  if( cmd == "a-map" )
  {
    auto tmp = abc_opto( ntk, cmd+std::to_string(CALL++), "if -a; fraig" );
    ntk = tmp;
    return true;
  }
  if( cmd == "lazy" )
  {
    auto tmp = abc_opto( ntk, cmd+std::to_string(CALL++), "if -y -K 6; fraig;" );
    ntk = tmp;
    return true;
  }
  return false;
}

template<class Ntk>
mig_network aig_to_mig( Ntk const& aig )
{
  /* library to map to MIGs */
  mig_npn_resynthesis resyn{ true };
  exact_library_params eps;
  eps.np_classification = true;
  exact_library<mig_network> exact_lib( resyn, eps );

  map_params ps1;
  ps1.skip_delay_round = false;
  ps1.required_time = std::numeric_limits<double>::max();
  map_stats st1;

  mig_network mig = map( aig, exact_lib, ps1, &st1 );
  return mig;
}

template<class Ntk>
aig_network mig_to_aig( Ntk const& mig, double required_time )
{
  /* library to map to MIGs */
  xag_npn_resynthesis<aig_network, xag_network, xag_npn_db_kind::aig_complete> aig_resyn{};
  exact_library_params eps;
  eps.np_classification = true;
  exact_library<aig_network> exact_lib(aig_resyn, eps);

  map_params ps1;
  ps1.skip_delay_round = false;
  ps1.required_time = required_time;
  map_stats st1;

  aig_network aig = map( mig, exact_lib, ps1, &st1 );
  return aig;
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