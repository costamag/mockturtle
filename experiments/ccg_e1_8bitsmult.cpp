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
#include <mockturtle/algorithms/ccgame/solvers/cusco.hpp>
#include <cstdint>
#include <string>
#include <bit>
#include <bitset>

using namespace mockturtle;
using namespace ccgame;

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


    /* define the parameters */
    solver_t type = solver_t::_COV_RND;
    int nIters = 20;
    cusco_ps ps( type, nIters );
    /* solve */
    cusco<xag_network> solver( xs, fns );
    solver.solve( ps );

  return 0;
}
