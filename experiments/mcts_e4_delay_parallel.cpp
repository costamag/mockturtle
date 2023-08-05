
#include <thread>
#include <kitty/constructors.hpp>
#include <mockturtle/algorithms/mcts/mct_tree.hpp>
#include <mockturtle/algorithms/mcts/method.hpp>
#include <mockturtle/algorithms/mcts/nodes/nd_delay.hpp>
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

#pragma region mutex
std::atomic<uint32_t> KEY;
std::vector<double> DELS;
std::vector<bool> USED;
std::vector<double> SIZS;
#pragma endregion

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

DTT create_from_integer( uint32_t int_tt )
{
  DTT res(4u);
  res ^= res;
  for( int i{0}; i<16; ++i )
  {
    if( ( int_tt >> i ) & 1u == 1u )
      kitty::set_bit( res, i );
    else
      kitty::clear_bit( res, i );
  }
  return res;
}

uint32_t tt_to_key( DTT tt )
{
  uint32_t uint_tt;
  std::string string_tt = kitty::to_hex(tt);
  sscanf( string_tt.c_str(), "%x", &uint_tt ); 
  return uint_tt & 0xFFFF;
}

DTT key_to_tt( uint32_t key )
{
  key = key & 0xFFFF;
  std::string bstring = "";
  for( int iBit{0}; iBit < 16u; ++iBit )
    bstring = ( ((key >> iBit) & 1u) == 1u) ? "1"+bstring : "0"+bstring;
  DTT res(4u);
  kitty::create_from_binary_string(res,bstring);
  return res;
}

struct result_mctsolve
{
  double delay;
  double area;
  bool isValid{false};
};

result_mctsolve mct_solve( kitty::dynamic_truth_table * );

void thread_run()
{
  kitty::dynamic_truth_table F = key_to_tt(KEY);
  kitty::next_inplace( F );
  KEY = tt_to_key(F);
  do
  {
    uint32_t key = tt_to_key( F );
    printf("FUNC %d\n", key );
    result_mctsolve rep = mct_solve( &F );
    if( rep.isValid )
    {
      USED[key]=true;
      DELS[key]=rep.delay;
      SIZS[key]=rep.area;
    }
    kitty::next_inplace( F );
    KEY = tt_to_key(F);
  } while ( !kitty::is_const0( F ) );
}

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

  printf(ANSI_COLOR_YELLOW   "DELAY EXPERIMENT 0: COMPARISON WITH EXACT SYNTHESIS" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_CYAN     "======================= ++++++++++ =========================="     ANSI_COLOR_RESET "\n\n");


  using TT = kitty::dynamic_truth_table;
  using NTK = aig_network;
  TT target(4u);



  for( int i{0}; i<pow(2,pow(2,4)); ++i )
  {
    DELS.push_back(0.0);
    USED.emplace_back(false);
    SIZS.emplace_back(0);
  }


  KEY.store( 0 );
  std::vector<std::thread> threads;

  /* generate threads */
  const auto processor_count = std::thread::hardware_concurrency();

  printf( "[i] Running on %d threads\n", processor_count );
  for ( auto i = 0u; i < processor_count; ++i )
    threads.emplace_back( thread_run );//, ps_contest, run_only_one );

  /* wait threads */
  for ( auto i = 0u; i < processor_count; ++i )
    threads[i].join();

  std::ofstream myfile;
  myfile.open ("MCTS1_0_0_4_4.txt");
  for( uint32_t i{0}; i<DELS.size(); ++i )
  {
    if( USED[i] )
      myfile << i << " " << DELS[i] << " " << SIZS[i] << "\n";
  }
  myfile.close();


  return 0;
}

result_mctsolve mct_solve( kitty::dynamic_truth_table * pF )
{
  result_mctsolve res;
  std::vector<double> T = {0,0,4,4};
  std::vector<kitty::dynamic_truth_table> X;

  for( int i{0}; i<4u; ++i )
  {
    X.emplace_back( 4u );
    kitty::create_nth_var( X[i], i );
  }

  node_ps ndps;
  detailed_gate_t cmpr_( gate_t::CMPR, 1, 0.5, 1.0, &hpcompute_cmpr );
  detailed_gate_t cmpl_( gate_t::CMPR, 1, 0.5, 1.0, &hpcompute_cmpl );
  detailed_gate_t ai00_( gate_t::AI00, 2, 1.0, 1.0, &hpcompute_ai00 );
  detailed_gate_t ai11_( gate_t::AI11, 2, 1.5, 1.0, &hpcompute_ai11 );
  detailed_gate_t exor_( gate_t::EXOR, 2, 2.0, 1.0, &hpcompute_exor );
  ndps.lib = { cmpl_, cmpr_, ai00_, ai11_, exor_ };

  mct_ps mctps;
  ndps.sel_type = supp_selection_t::SUP_NORM;
  mctps.nIters = 5;
  mctps.nSims = 5;
  mctps.verbose =false;
  ndps.BETA0 = 100;
  ndps.BETAZ = 100;
  ndps.nIters = 1;
  ndps.thresh = 5;
  ndps.delay_inv = 0.5;

  nd_delay_t<xag_network> root( X, T, {*pF}, ndps );
  mct_method_ps metps;

  mct_method_t<nd_delay_t<xag_network> > meth( metps );
  mct_tree_t<nd_delay_t<xag_network> , mct_method_t> mct( root, meth, mctps );
  int iSol = mct.solve();
  if( iSol == -1 ) 
  {
    printf("no solution found\n");
    return res;
  }
  else
    res.isValid = true;
  xag_network xag = mct.nodes[iSol].ntk;
  //printf( "size %d || delay %f\n", xag.num_gates(), mct.evaluate(iSol) );
  res.area = xag.num_gates();
  res.delay = mct.evaluate(iSol);
  default_simulator<kitty::dynamic_truth_table> sim( 4u );
  const auto tt = simulate<kitty::dynamic_truth_table>( xag, sim )[0];
  //kitty::print_binary(*pF);
  //printf("\n");
  //kitty::print_binary(tt);
  assert( kitty::equal( tt, *pF ) );
  return res;
} 