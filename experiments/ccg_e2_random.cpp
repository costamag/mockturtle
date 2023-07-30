#include <kitty/constructors.hpp>
#include <mockturtle/algorithms/decompose/DecSolver.hpp>
#include <mockturtle/algorithms/ccgame/solvers/cusco.hpp>
#include <mockturtle/algorithms/dcsynthesis/dc_solver.hpp>
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
using namespace ccgame;

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

void propose_thresh( kitty::dynamic_truth_table *, int );
void propose_khot( kitty::dynamic_truth_table *, int );

kitty::dynamic_truth_table propose_game( int, int, int );
kitty::dynamic_truth_table userdef_game();
kitty::dynamic_truth_table knuth_game( int );

template<class Ntk> report_t<Ntk> game_on( kitty::dynamic_truth_table *, int, int );

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

  printf(ANSI_COLOR_YELLOW   "EXPERIMENT 2: SUB-OPTIMALITY OF THE HEURISTIC" ANSI_COLOR_RESET "" );
  printf(ANSI_COLOR_YELLOW   " In this experiment we show that multiple not-informad runs of  " ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW   "the heuristic can yield better results." ANSI_COLOR_RESET "\n\n" );


  std::string dot_file  = "";
  std::string blif_file = "";
  std::string aig_file  = "";
  std::vector<std::string> fNames = { "THRESH", "ONEHOT" };
  for( int iFn{0}; iFn<2; ++iFn )
  {
    printf("%s\n", fNames[iFn].c_str() );
    printf("%2s ||%20s|%20s|%20s|%20s|%20s|%20s|%20s|%20s|\n", "n", "S1", "S2", "S3", "S3", "S5", "S6", "S7", "S8" );
    for( uint32_t nVars{2}; nVars < 10; ++nVars )
    {
      std::string info = fmt::format("{:2} ||", nVars );
      for( uint32_t iThr{1}; iThr < ( floor( nVars/2 ) + 2 ); ++iThr )
      {
        kitty::dynamic_truth_table F;
        F = propose_game( nVars, iThr, iFn );

        std::clock_t start;
        double duration;
        start = std::clock();

        report_t<xag_network> rep = game_on<xag_network>(&F, 0, 33 ); // set to random sampling 1 for 33 iterations @.@
        xag_network xag = rep.ntk;
        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

        resubstitution_params ps;
        ps.max_pis = xag.num_pis();
        ps.max_inserts = 20u;
        ps.max_divisors = 1000u;
        ps.odc_levels = -1;
        ps.progress = true;

        auto xag_resub = cleanup_dangling( xag );
        sim_resubstitution( xag_resub, ps );
        xag_resub = cleanup_dangling( xag_resub );

        info += fmt::format("{:2d}.{:3d}.{:2d}>{:3d} {:6.2f}|", (rep.nIt0-rep.nMin), rep.nMin, (rep.nMax-rep.nMin), xag_resub.num_gates(), duration );
        std::string sName = fmt::format("S{:2d}_{:2d}", nVars, iThr );
        dot_file  = "EXPS/EXP2/" + fNames[iFn] +"/dot/"  + sName + ".dot";
        blif_file = "EXPS/EXP2/" + fNames[iFn] +"/blif/" + sName + ".blif";
        aig_file  = "EXPS/EXP2/" + fNames[iFn] +"/aig/"  + sName + ".aig";
        write_dot( xag, dot_file );
        write_dot( xag, blif_file );
        write_aiger( xag, aig_file );
        dot_file  = "EXPS/EXP2/" + fNames[iFn] +"/dot/"  + sName + "rs.dot";
        blif_file = "EXPS/EXP2/" + fNames[iFn] +"/blif/" + sName + "rs.blif";
        aig_file  = "EXPS/EXP2/" + fNames[iFn] +"/aig/"  + sName + "rs.aig";
      }
      printf("%s\n", info.c_str() );
    }
  }
  return 0;
}

void propose_gamble( kitty::dynamic_truth_table * pF, std::string * pInfo )
{
    kitty::dynamic_truth_table gambleP = (*pF).construct();
    kitty::dynamic_truth_table gambleN = (*pF).construct();

    gambleP = gambleP | ~gambleP;
    gambleN = gambleN | ~gambleN;
    int nVars = (*pF).num_vars();
    for( int j{0u}; j < nVars; ++j )
    {
      kitty::dynamic_truth_table x(nVars);
      kitty::create_nth_var( x, j );
      gambleP &= x ;
      gambleN &= ~x ;
    }
    *pF = gambleP | gambleN;

    *pInfo += "gamble/s" + std::to_string(nVars);
}

void propose_parity( kitty::dynamic_truth_table * pF, std::string * pInfo )
{
    kitty::create_parity(*pF);
    int nVars = (*pF).num_vars();
    *pInfo += "parity/s" + std::to_string(nVars);
}

void propose_symmetric( kitty::dynamic_truth_table * pF, std::vector<uint32_t> vals )
{
  uint32_t nVars = (*pF).num_vars();
  *pF = *pF ^ *pF;
  for( auto v : vals )
  {
    assert( v < (nVars + 1) );
    for( uint64_t j{0}; j < pow(2,nVars); ++j )
    {
      unsigned int POP_CNT=0;
      uint64_t mint = j;
      while( mint )
      {
        POP_CNT += mint & 1;
        mint >>= 1;
      }

      if( POP_CNT == v )
        kitty::set_bit( *pF, j );
    }
  }
}

void propose_khot( kitty::dynamic_truth_table * pF, int popcount )
{
  int nVars = (*pF).num_vars();
  *pF = *pF ^ *pF;
  for( uint64_t j{0}; j < pow(2,nVars); ++j )
  {
    unsigned int POP_CNT=0;
    uint64_t mint = j;
    while( mint )
    {
      POP_CNT += mint & 1;
      mint >>= 1;
    }
    if( POP_CNT == popcount )
      kitty::set_bit( *pF, j );
  }
}

void propose_thresh( kitty::dynamic_truth_table * pF, int thresh )
{
  kitty::create_threshold( *pF, thresh );
}

kitty::dynamic_truth_table propose_game( int nVars, int iThr, int Id )
{
    kitty::dynamic_truth_table F( nVars );
 //   printf(ANSI_COLOR_YELLOW " 0 THRESHOLD " ANSI_COLOR_RESET "\n" );
 //   printf(ANSI_COLOR_YELLOW " 1 k-HOT    " ANSI_COLOR_RESET "\n" );

    switch (Id)
    {
    case 0:
        propose_thresh( &F, iThr );
        break;
    case 1:
        propose_khot( &F, iThr );
        break;
    default:
        break;
    }
    return F;
}

kitty::dynamic_truth_table userdef_game()
{
  printf(ANSI_COLOR_YELLOW " 0 CREATE FROM BINARY " ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " 1 CREATE FROM HEX    " ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " CHOOSE THE INPUT ENCODING: " ANSI_COLOR_RESET "" );
  int ENC;
  std::cin >> ENC;
  std::string ISTR;
  
  if( ENC == 0 )
  {
    printf(ANSI_COLOR_YELLOW " ENTER THE BINARY STRING: " ANSI_COLOR_RESET "" );
    std::cin >> ISTR;
    int nBits = ISTR.size();
    if( nBits % 2 != 0 )
    {
      printf(ANSI_COLOR_RED " BAD FUNCTION DEFINITION " ANSI_COLOR_RESET "" );   
      assert(0);
    }
    int nVars = log2( nBits );
    kitty::dynamic_truth_table F( nVars );
    kitty::create_from_binary_string( F, ISTR );
    return F;
  }
  else if( ENC == 1 )
  {
    printf(ANSI_COLOR_YELLOW " ENTER THE HEX STRING: " ANSI_COLOR_RESET "" );
    std::cin >> ISTR;
    int nHexs = ISTR.size();
    if( ( nHexs % 2 != 0 ) && ( nHexs != 1 ) )
    {
      printf(ANSI_COLOR_RED " BAD FUNCTION DEFINITION " ANSI_COLOR_RESET "" );   
      assert(0);
    }
    int nVars = log2( nHexs ) + 2;
    kitty::dynamic_truth_table F( nVars );
    kitty::create_from_hex_string( F, ISTR );
    return F;
  }
  else
  {
    printf(ANSI_COLOR_RED " CHOICE NOT MATCHING ANY ENCODING " ANSI_COLOR_RESET "\n" );   
    assert(0);
  }
}

kitty::dynamic_truth_table knuth_game( int idGame )
{
  int nVars;
  if( idGame < 15 )
    nVars = 4;
  else if( idGame < 31 )
    nVars = 5;
  else
  {
    printf(ANSI_COLOR_RED " PROBLEM NOT DEFINED BY KNUTH" ANSI_COLOR_RESET "\n" );
    assert(0);
  }

  kitty::dynamic_truth_table F( nVars );
  std::vector<uint32_t> VALS;
  switch (idGame)
  {
    case 0 : VALS = {4};        break;
    case 1 : VALS = {3};        break;
    case 2 : VALS = {3,4};      break;
    case 3 : VALS = {2};        break;
    case 4 : VALS = {2,4};      break;
    case 5 : VALS = {2,3};      break;
    case 6 : VALS = {2,3,4};    break;
    case 7 : VALS = {1};        break;
    case 8 : VALS = {1,4};      break;
    case 9 : VALS = {1,3};      break;
    case 10: VALS = {1,3,4};    break;
    case 11: VALS = {1,2};      break;
    case 12: VALS = {1,2,4};    break;
    case 13: VALS = {1,2,3};    break;
    case 14: VALS = {1,2,3,4};  break;
    case 15: VALS = {4};        break;
    case 16: VALS = {4,5};      break;
    case 17: VALS = {3};        break;
    case 18: VALS = {3,5};      break;
    case 19: VALS = {3,4};      break;
    case 20: VALS = {3,4,5};    break;
    case 21: VALS = {2,5};      break;
    case 22: VALS = {2,4};      break;
    case 23: VALS = {2,4,5};    break;
    case 24: VALS = {2,3,5};    break;
    case 25: VALS = {2,3};      break;
    case 26: VALS = {2,3,4};    break;
    case 27: VALS = {1,5};      break;
    case 28: VALS = {1,4};      break;
    case 29: VALS = {1,3,4};    break;
    case 30: VALS = {1,2,5};    break;
    default:
      break;
  }
  propose_symmetric( &F, VALS );
  return F;
}

template<class Ntk>
report_t<Ntk> game_on( kitty::dynamic_truth_table * pF, int MET, int NITERS )
{
  using TT  = kitty::dynamic_truth_table;

  typedef DecSolver<TT, Ntk> solver_t;
  report_t<Ntk> rep;
  
  kitty::dynamic_truth_table mask = (*pF).construct();
  mask = mask | ~mask;
  solver_t solver( {*pF}, {mask} );
  //solver.PrintSpecs();
  
  if( MET == 0 )
  {
    std::vector<kitty::dynamic_truth_table> xs;
    for( int i{0}; i<pF->num_vars(); ++i )
    {
      xs.emplace_back( pF->num_vars() );
      kitty::create_nth_var( xs[i], i );
    }
    /* define the parameters */
    int nIters = 1;
    cusco_ps ps( ccgame::solver_t::_SYM_1SH, nIters );
    /* solve */
    cusco<Ntk> solver( xs, {*pF} );
    rep = solver.solve( ps );
  }
  else if( MET == 1 )
  {
    //ntk = solver.aut_sym_solve( MET );

    std::vector<kitty::dynamic_truth_table> xs;
    for( int i{0}; i<pF->num_vars(); ++i )
    {
      xs.emplace_back( pF->num_vars() );
      kitty::create_nth_var( xs[i], i );
    }
  
    /* define the parameters */
    int nIters = NITERS;
    cusco_ps ps( ccgame::solver_t::_SYM_RND, nIters );
    /* solve */
    cusco<Ntk> solver( xs, {*pF} );
    rep = solver.solve( ps );
  }
  else if( MET == 2 )
  {
    std::vector<kitty::dynamic_truth_table> xs;
    for( int i{0}; i<pF->num_vars(); ++i )
    {
      xs.emplace_back( pF->num_vars() );
      kitty::create_nth_var( xs[i], i );
    }

    /* define the parameters */    
    int nIters = NITERS;
    cusco_ps ps( ccgame::solver_t::_COV_RND, nIters, -1 );
    /* solve */
    cusco<Ntk> solver( xs, {*pF} );
    rep = solver.solve( ps );
  }
  else
  {
    printf(ANSI_COLOR_RED " CHOICE NOT MATCHING ANY METHOD " ANSI_COLOR_RESET "\n" );   
    assert(0);
  }
  default_simulator<kitty::dynamic_truth_table> sim( (*pF).num_vars() );
  const auto tt = simulate<kitty::dynamic_truth_table>( rep.ntk, sim )[0];
  assert( kitty::equal( tt, *pF ) );


  return rep;

}