#include <kitty/constructors.hpp>
#include <kitty/print.hpp>
#include <kitty/operations.hpp>
#include <kitty/operators.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/write_dot.hpp>
#include <mockturtle/algorithms/sim_resub.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/io/truth_reader.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <lorina/truth.hpp>
#include <mockturtle/algorithms/mcts/mct_tree.hpp>
#include <mockturtle/algorithms/mcts/method.hpp>
#include <mockturtle/algorithms/mcts/nodes/nd_size.hpp>
#include <fmt/format.h>
#include <cstdint>
#include <string>
#include <bit>
#include <bitset>



using namespace mockturtle;
using namespace mcts;

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"


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

  printf(ANSI_COLOR_YELLOW " DO YOU WANT ME TO PROPOSE YOU A GAME [Y/N/K]? " ANSI_COLOR_RESET "" );
  
  std::string info = "";

  std::string benchmark_path = "../experiments/";
  std::string benchmark = "mul4"; 

  klut_network klut;

  auto res0 = lorina::read_truth( benchmark_path + benchmark + ".truth", truth_reader( klut ) );
  if ( res0 != lorina::return_code::success )
  {
    printf(ANSI_COLOR_RED " READ FAILED " ANSI_COLOR_RESET "\n" );
    return 1;
  }
  std::vector<kitty::dynamic_truth_table> xs;
  for( int i{0}; i<8; ++i )
  {
    xs.emplace_back(8u);
    kitty::create_nth_var( xs[i], i );
  }

  std::vector<kitty::dynamic_truth_table> fns;
  klut.foreach_po( [&]( const auto& x, auto index ) {
      fns.emplace_back(8u);
      kitty::create_from_binary_string(fns[index], kitty::to_binary( klut.node_function( x ) ));
    } );

/*
    dc_solver<xag_network> solver( xs, fns );
    xag_network xag;
    solver.solve_greedy_multioutput( &xag );
*/

    std::vector<double> ts;
    for( int i{0}; i<8u; ++i )
    {
      ts.push_back( 0u );
      kitty::create_nth_var( xs[i], i );
    }

    node_ps ndps;
    detailed_gate_t ai00_( gate_t::AI00, 2, 1.0, 1.0, &hpcompute_ai00 );
    detailed_gate_t ai01_( gate_t::AI01, 2, 1.0, 1.0, &hpcompute_ai01 );
    detailed_gate_t ai10_( gate_t::AI10, 2, 1.0, 1.0, &hpcompute_ai10 );
    detailed_gate_t ai11_( gate_t::AI11, 2, 1.0, 1.0, &hpcompute_ai11 );
    detailed_gate_t exor_( gate_t::EXOR, 2, 1.0, 1.0, &hpcompute_exor );
    ndps.lib = { ai00_, ai01_, ai10_, ai11_, exor_ };

    mct_ps mctps;
    ndps.sel_type = supp_selection_t::SUP_ENER;

    mctps.nIters = 1;
    mctps.nSims = 1;
    mctps.verbose = true;
    ndps.BETA0 = 100;
    ndps.thresh = 10;
    ndps.BETAZ = 1;
    ndps.nIters = 10;


    nd_size_t<xag_network> root( xs, ts, fns, ndps );
    mct_method_ps metps;

    mct_method_t<nd_size_t<xag_network> > meth( metps );
    mct_tree_t<nd_size_t<xag_network> , mct_method_t> mct( root, meth, mctps );
    int iSol = mct.solve();
    if( iSol == -1 ) printf("no solution found\n");

  return 0;
}
