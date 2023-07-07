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

  
  std::vector<double> T = {0,0,4,4};
  std::vector<kitty::dynamic_truth_table> X;

  printf( "ENTER INTEGER: " );
  uint32_t INTEGER;
  std::cin >> INTEGER;

  DTT F = create_from_integer( INTEGER );
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
  ndps.sel_type = supp_selection_t::SUP_ENER;
  mctps.nIters =100;
  mctps.nSims = 1;
  mctps.verbose =true;
  ndps.BETA0 = 100;
  ndps.nIters = 100;

  nd_delay_t<xag_network> root( X, T, {F}, ndps );
  mct_method_ps metps;

  mct_method_t<nd_delay_t<xag_network> > meth( metps );
  mct_tree_t<nd_delay_t<xag_network> , mct_method_t> mct( root, meth, mctps );
  int iSol = mct.solve();
  if( iSol == -1 ) printf("no solution found\n");
  xag_network xag = mct.nodes[iSol].ntk;
  printf( "size %d || delay %f\n", xag.num_gates(), mct.evaluate(iSol) );

  default_simulator<kitty::dynamic_truth_table> sim( 4u );
  const auto tt = simulate<kitty::dynamic_truth_table>( xag, sim )[0];
kitty::print_binary(F);
printf("\n");
kitty::print_binary(tt);
  assert( kitty::equal( tt, F ) );

  return 0;
}