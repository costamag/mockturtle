#include <kitty/constructors.hpp>
#include <kitty/print.hpp>
#include <kitty/operations.hpp>
#include <kitty/operators.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/algorithms/decompose/DecSolver.hpp>
#include <mockturtle/algorithms/ccgame/solvers/cusco.hpp>
#include <mockturtle/algorithms/dcsynthesis/dc_solver.hpp>
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

void analysis();
void solve();

template<class Ntk>
Ntk game_on( std::vector<kitty::dynamic_truth_table> );

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

  printf(ANSI_COLOR_YELLOW " ANALYSIS [A] OR SOLVING [S]? " ANSI_COLOR_RESET "" );
  char TODO;
  std::cin >> TODO;
  
  if( TODO == 'S' )
    solve();
  else
    analysis();

  /* let's see what they gave us */


/*
    dc_solver<xag_network> solver( xs, fns );
    xag_network xag;
    solver.solve_greedy_multioutput( &xag );
*/




  return 0;
}

void analysis()
{
  std::string sBench;
  for( int iBench{0}; iBench < 100; iBench++ )
  {
    sBench = iBench < 10 ? fmt::format( "ex0{:d}", iBench ) : fmt::format( "ex{:d}", iBench );
    std::string benchmark = "../experiments/IWLS_2023/" + sBench + ".truth";
    klut_network klut;
    auto res0 = lorina::read_truth( benchmark, truth_reader( klut ) );
    if ( res0 != lorina::return_code::success )
    {
      printf(ANSI_COLOR_RED " READ %s FAILED " ANSI_COLOR_RESET "\n", benchmark.c_str() );
      assert(0);
    }
    int nIns = klut.num_pis();
    int nOuts = klut.num_pos();
    printf( "%2d nIns=%2d nOuts=%2d\n", iBench, nIns, nOuts );

    std::vector<kitty::dynamic_truth_table> xs;
    for( int i{0}; i<nIns; ++i )
    {
      xs.emplace_back(nIns);
      kitty::create_nth_var( xs[i], i );
    }

    std::vector<kitty::dynamic_truth_table> fns;
    klut.foreach_po( [&]( const auto& x, auto index ) {
        fns.emplace_back(nIns);
        kitty::create_from_binary_string(fns[index], kitty::to_binary( klut.node_function( x ) ));
      } );
  }
}

void solve()
{
  printf(ANSI_COLOR_YELLOW " BENCHMARK: " ANSI_COLOR_RESET "" );

  int iBench;
  std::string sBench = "";
  std::cin >> iBench;
  sBench = iBench < 10 ? fmt::format( "ex0{:d}", iBench ) : fmt::format( "ex{:d}", iBench );

  std::string benchmark = "../experiments/IWLS_2023/" + sBench + ".truth";

  klut_network klut;
  auto res0 = lorina::read_truth( benchmark, truth_reader( klut ) );
  if ( res0 != lorina::return_code::success )
  {
    printf(ANSI_COLOR_RED " READ %s FAILED " ANSI_COLOR_RESET "\n", benchmark.c_str() );
    assert(0);
  }
  int nIns = klut.num_pis();
  int nOuts = klut.num_pos();
  printf( "nIns=%2d nOuts=%2d\n", nIns, nOuts );

   std::vector<kitty::dynamic_truth_table> xs;
  for( int i{0}; i<nIns; ++i )
  {
    xs.emplace_back(nIns);
    kitty::create_nth_var( xs[i], i );
  }

   std::vector<kitty::dynamic_truth_table> fns;
  klut.foreach_po( [&]( const auto& x, auto index ) {
      fns.emplace_back(nIns);
      kitty::create_from_binary_string(fns[index], kitty::to_binary( klut.node_function( x ) ));
    } );

    game_on<xag_network>( fns );
}




template<class Ntk>
Ntk game_on( std::vector<kitty::dynamic_truth_table> vF )
{
  using TT  = kitty::dynamic_truth_table;
  std::vector<kitty::dynamic_truth_table> vM;
  typedef DecSolver<TT, Ntk> solver_t;
  Ntk ntk;
  
  kitty::dynamic_truth_table mask = vF[0].construct();
  mask = mask | ~mask;
  for( int iOut = 0; iOut < vF.size(); ++iOut )
    vM.push_back( mask );
  solver_t solver( vF, vM );

  printf(ANSI_COLOR_YELLOW " 0 SYM MANUAL" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " 1 DEC MANUAL" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " 2 SYM AUTOMATIC" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " 3 DEC AUTOMATIC" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " 4 DEC AUTOMATIC WEAK" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " 5 SYM MANUAL RS" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " 6 SYM AUTOMATIC RS" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " 7 SYM AUTOMATIC XOR" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " ===================" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " 8 CGG-RELAX" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " 9 CGG-XOR" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " 10 CGG-SPEC" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " 11 CGG-X" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " ===================" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " =   NEW VERSION   =" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " ===================" ANSI_COLOR_RESET "\n" );

  printf(ANSI_COLOR_YELLOW " CHOOSE YOUR METHOD: " ANSI_COLOR_RESET "" );
  int MET;
  
  std::cin >> MET;
  //MET = 2; printf("CCG-DEC\n");  
  
  if( MET == 0 )
    ntk = solver.man_sym_solve();
  else if( MET == 1 )
    ntk = solver.man_rdec_solve();
  else if( MET == 2 )
  {
    printf(ANSI_COLOR_YELLOW " NUMBER OF ITERATIONS: " ANSI_COLOR_RESET "" );
    
    std::cin >> MET;
    //MET = 33; printf("33\n");  
  
    ntk = solver.aut_sym_solve( MET );
  }
  else if( MET == 3 )
  {
    printf(ANSI_COLOR_YELLOW " NUMBER OF ITERATIONS: " ANSI_COLOR_RESET "" );
    std::cin >> MET;
    ntk = solver.aut_rdec_solve( MET );
  }
  else if( MET == 4 )
  {
    printf(ANSI_COLOR_YELLOW " NUMBER OF ITERATIONS: " ANSI_COLOR_RESET "" );
    std::cin >> MET;
    ntk = solver.aut_symGT_solve( MET );
  }
  else if( MET == 5 )
  {
    ntk = solver.man_sym_solve_rs( );
  }
  else if( MET == 6 )
  {
    printf(ANSI_COLOR_YELLOW " NUMBER OF ITERATIONS: " ANSI_COLOR_RESET "" );
    std::cin >> MET;
    ntk = solver.aut_sym_solve_rs( MET );
  }
  else if( MET == 7 )
  {
    printf(ANSI_COLOR_YELLOW " NUMBER OF ITERATIONS: " ANSI_COLOR_RESET "" );
    std::cin >> MET;
    ntk = solver.aut_sym_solve_xor( MET );
  }
  else if( MET == 8 )
  {
    printf(ANSI_COLOR_YELLOW " NUMBER OF ITERATIONS: " ANSI_COLOR_RESET "" );
    
    //std::cin >> MET;
    MET = 10; printf("10s\n");
    
    ntk = solver.ccg_relax( MET );
  }
  else if( MET == 9 )
  {
    printf(ANSI_COLOR_YELLOW " NUMBER OF ITERATIONS: " ANSI_COLOR_RESET "" );
    
    std::cin >> MET;
    //MET = 10; printf("20s\n");
    
    ntk = solver.ccg_xor( MET );
  }
  else if( MET == 10 )
  {
    printf(ANSI_COLOR_YELLOW " TIME[s]: " ANSI_COLOR_RESET "" );
    
    std::cin >> MET;
    //MET = 20; printf("20s\n");
    int PRC;
    printf(ANSI_COLOR_YELLOW " PERCENTAGE [-1,100]: " ANSI_COLOR_RESET "" );
    std::cin >> PRC;

    ntk = solver.ccg_spectral( MET, PRC );
  }
  else if( MET == 11 )
  {
    printf(ANSI_COLOR_YELLOW " NUMBER OF ITERATIONS: " ANSI_COLOR_RESET "" );
    
    std::cin >> MET;
    ntk = solver.ccgX( MET );
  }
  /*else if( MET == 12 )
  {
    printf(ANSI_COLOR_YELLOW " NUMBER OF ITERATIONS: " ANSI_COLOR_RESET "" );
    MET = 100;
    //std::cin >> MET;
    ntk = sym_solver.ccg_sym( MET );
  }*/
  else
  {
    printf(ANSI_COLOR_RED " CHOICE NOT MATCHING ANY METHOD " ANSI_COLOR_RESET "\n" );   
    assert(0);
  }


  printf("best #nodes: %d\n", ntk.num_gates() );

  default_simulator<kitty::dynamic_truth_table> sim( vF[0].num_vars() );
  const auto tt = simulate<kitty::dynamic_truth_table>( ntk, sim )[0];

  //if( tt.num_vars() < 10 )
  ///{
  //  printf("\n simulation\n");
  //  kitty::print_binary( tt );
  //  std::cout << std::endl;
  //  printf("function returned\n");
  //  kitty::print_binary( vF );
  //  std::cout << std::endl;
  //}

  //std::cout << ( kitty::equal( tt, vF[] ) ? " equal " : " different " ) << std::endl;

  return ntk;

}

