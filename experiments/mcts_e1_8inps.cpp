#include <kitty/constructors.hpp>
#include <mockturtle/algorithms/mcts/mct_tree.hpp>
#include <mockturtle/algorithms/mcts/method.hpp>
#include <mockturtle/algorithms/mcts/nodes/nd_size.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/write_dot.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/algorithms/sim_resub.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <string>
#include <bit>
#include <bitset>
#include <cstdint>
#include <fmt/format.h>

using namespace mockturtle;
using namespace mcts;

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

template<class Ntk> Ntk game_on( kitty::dynamic_truth_table * );

template<class Ntk> Ntk abc_deepsyn( DTT );
template<class Ntk> Ntk abc_transduction( DTT );

int main()
{
  printf(ANSI_COLOR_RED     "============================================================="     ANSI_COLOR_RESET "\n");
  printf(ANSI_COLOR_RED     "============================================================="     ANSI_COLOR_RESET "\n\n");

  printf(ANSI_COLOR_RED     "  ####         ####         ####     ####   #      # ########"     ANSI_COLOR_RESET "\n");
  printf(ANSI_COLOR_RED     " ######       ######       ######   ######  ##    ## ########"     ANSI_COLOR_RESET "\n");
  printf(ANSI_COLOR_GREEN   "###  ###     ###  ###     ###  ### ###  ### ###  ### ##      "   ANSI_COLOR_RESET "\n");
  printf(ANSI_COLOR_GREEN   "##    ##     ##    ##     ##    ## ##    ## ######## ##      "   ANSI_COLOR_RESET "\n");
  printf(ANSI_COLOR_YELLOW  "##           ##           ##       ##    ## ## ## ## ##      "  ANSI_COLOR_RESET "\n");
  printf(ANSI_COLOR_YELLOW  "##           ##           ##       ##    ## ## ## ## #####   "  ANSI_COLOR_RESET "\n");
  printf(ANSI_COLOR_BLUE    "##           ##           ##  #### ######## ## ## ## #####   "    ANSI_COLOR_RESET "\n");
  printf(ANSI_COLOR_BLUE    "##           ##           ##  #### ######## ## ## ## ##      "    ANSI_COLOR_RESET "\n");
  printf(ANSI_COLOR_MAGENTA "##    ##     ##    ##     ##    ## ##    ## ##    ## ##      " ANSI_COLOR_RESET "\n");
  printf(ANSI_COLOR_MAGENTA "##   ###     ##   ###     ###  ### ##    ## ##    ## ##      " ANSI_COLOR_RESET "\n");
  printf(ANSI_COLOR_CYAN    " ######  ##   ######  ##   ######  ##    ## ##    ## ########"    ANSI_COLOR_RESET "\n");
  printf(ANSI_COLOR_CYAN    "  ####   ##    ####   ##    ####   ##    ## ##    ## ########"    ANSI_COLOR_RESET "\n");
  printf( "\n\n" );
  printf(ANSI_COLOR_CYAN     "============================================================="     ANSI_COLOR_RESET "\n");
  printf(ANSI_COLOR_CYAN     "============================================================="     ANSI_COLOR_RESET "\n\n");

  printf(ANSI_COLOR_YELLOW   "EXPERIMENT 3: COMPARISON WITH EXACT SYNTHESIS" ANSI_COLOR_RESET "" );
  printf(ANSI_COLOR_YELLOW   " In this experiment we compare the CUSCO heuristic with the" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW   "with the exact synthesis results obtained by Knuth [1]" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_CYAN     "======================= REFERENCES =========================="     ANSI_COLOR_RESET "\n\n");
  printf(ANSI_COLOR_YELLOW " [1] Knuth: 'The art of computer programming' fascicle 1 vol. 4" ANSI_COLOR_RESET "\n\n" );


  kitty::dynamic_truth_table F(3u);

  printf("function | deepsyn  | transd 1 | transd N\n");

  for( uint32_t i{0}; i<20; ++i )
  {
    kitty::create_random( F );
    aig_network aig0 = abc_deepsyn<aig_network>( F );
    aig_network aig1 = abc_transduction<aig_network>( F );
    aig_network aig = game_on<aig_network>(&F);
    kitty::print_hex(F);
    printf(" | %8d | %8d | %8d \n", aig0.num_gates(), aig1.num_gates(), aig.num_gates());

  }
  
  return 0;
}

template<class Ntk>
Ntk game_on( kitty::dynamic_truth_table * pF )
{
  using TT  = kitty::dynamic_truth_table;

  Ntk rep;
  
  kitty::dynamic_truth_table mask = (*pF).construct();
  mask = mask | ~mask;
  //solver.PrintSpecs();
  
  std::vector<double> ts;
  std::vector<kitty::dynamic_truth_table> xs;
  for( int i{0}; i<pF->num_vars(); ++i )
  {
    ts.push_back(0.0);
    xs.emplace_back( pF->num_vars() );
    kitty::create_nth_var( xs[i], i );
  }

  node_ps ndps;
  detailed_gate_t ai00_( gate_t::AI00, 2, 1.0, 1.0, &hpcompute_ai00 );
  detailed_gate_t ai01_( gate_t::AI01, 2, 1.0, 1.0, &hpcompute_ai01 );
  detailed_gate_t ai10_( gate_t::AI10, 2, 1.0, 1.0, &hpcompute_ai10 );
  detailed_gate_t ai11_( gate_t::AI11, 2, 1.0, 1.0, &hpcompute_ai11 );
  detailed_gate_t exor_( gate_t::EXOR, 2, 2.0, 1.0, &hpcompute_exor );
  ndps.lib = { ai00_, ai01_, ai10_, ai11_};//, exor_ };

  mct_ps mctps;
  ndps.sel_type = supp_selection_t::SUP_NORM;
  //mctps.nIters =100;
  //mctps.nSims = 1;
  //mctps.verbose =true;
  //ndps.BETA0 = 100;
  //ndps.nIters = 100;


  mctps.nIters = 20;
  mctps.nSims = 10;
  mctps.verbose = true;
  ndps.BETA0 = 20;
  ndps.BETAZ = 20;
  ndps.nIters = 5;//Ntrails
  ndps.thresh = 15;

  nd_size_t<aig_network> root( xs, ts, {*pF}, ndps );
  mct_method_ps metps;
  metps.sel_type = node_selection_t::NODE_LAY0;

  mct_method_t<nd_size_t<aig_network> > meth( metps );
  mct_tree_t<nd_size_t<aig_network> , mct_method_t> mct( root, meth, mctps );
  int iSol = mct.solve();
  if( iSol == -1 ) printf("ERROR no solution found\n");
  rep = mct.nodes[iSol].ntk;
  std::string nameFile = "ntk.dot";
  write_dot(rep, nameFile );
  //mct.path_print( iSol );

  default_simulator<kitty::dynamic_truth_table> sim( (*pF).num_vars() );
  const auto tt = simulate<kitty::dynamic_truth_table>( rep, sim )[0];
  assert( kitty::equal( tt, *pF ) );

  return rep;

}

template<class Ntk>
Ntk abc_transduction( DTT truth )
{
  Ntk res;

  std::string command = "abc -q \"read_truth -x " + kitty::to_binary(truth) + "; fraig; &get; &transduction -T 8; &put; write_aiger /tmp/pre.aig\"";

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

  std::string string_path = ( "/tmp/pre.aig" );
  if( lorina::read_aiger( string_path, aiger_reader( res ) ) != lorina::return_code::success )
    std::cerr << "read_aiger failed" << std::endl;

  return res;
}

template<class Ntk>
Ntk abc_deepsyn( DTT truth )
{
  Ntk res;

  std::string command = "abc -q \"read_truth -x " + kitty::to_binary(truth) + "; fraig; &get; &deepsyn -I 10 -J 100; &put; write_aiger /tmp/pre.aig\"";

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

  std::string string_path = ( "/tmp/pre.aig" );
  if( lorina::read_aiger( string_path, aiger_reader( res ) ) != lorina::return_code::success )
    std::cerr << "read_aiger failed" << std::endl;

  return res;
}


/*
    ndps.BETA0 = 10000;
    ndps.BETAZ = 10000;
    ndps.nIters = 100;//Ntrails
    ndps.thresh = 20;


    ndps.BETA0 = 10;
    ndps.BETAZ = 10;
    ndps.nIters = 1;//Ntrails
    ndps.thresh = 20;
*/
