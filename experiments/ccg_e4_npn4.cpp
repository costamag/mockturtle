#include <kitty/constructors.hpp>
#include <kitty/print.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <mockturtle/algorithms/decompose/DecSolver.hpp>
#include <mockturtle/algorithms/ccgame/solvers/cusco.hpp>
#include <mockturtle/algorithms/dcsynthesis/dc_solver.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/resyn_engines/xag_resyn.hpp>
#include <mockturtle/algorithms/node_resynthesis/exact.hpp>
#include <mockturtle/algorithms/resyn_engines/xag_resyn.hpp>
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

#include <iostream>
#include <fstream>
#include <string>
#include <set>

#include <sys/types.h>
#include <sys/stat.h>

using namespace mockturtle;
using namespace ccgame;

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

void propose_thresh( kitty::dynamic_truth_table *, std::string * );
void propose_gamble( kitty::dynamic_truth_table *, std::string * );
void propose_khot( kitty::dynamic_truth_table *, std::string * );
void propose_parity( kitty::dynamic_truth_table *, std::string * );

kitty::dynamic_truth_table propose_game(std::string*);
kitty::dynamic_truth_table userdef_game();
kitty::dynamic_truth_table knuth_game( int );

template<class Ntk> report_t<Ntk> game_on( kitty::dynamic_truth_table *, int );


std::vector<int> Nsat;
std::vector<int> Ncov;
std::vector<int> Nsym;
uint32_t Nacc;

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


  std::string dot_file  ;
  std::string blif_file ;
  std::string aig_file  ;
  std::vector<int> MET = {1,1};


  printf("%20s| %7s %8s %7s | %7s %8s %7s |\n", "", "", "SYM-CUSCO", "", "", "COV-CUSCO", "" );
  printf("%10s | %6s | %6s | %6s | %6s | %6s | %6s | %6s |\n", "f", "C(f)", "AIG", "AIG*", "T[s]", "AIG", "AIG*", "T[s]" );
  
  using TT = kitty::dynamic_truth_table;
  using Ntk = aig_network;
  
  // prepare the database for lookup
  xag_npn_resynthesis<Ntk> resyn;
  xag_npn_resynthesis<Ntk, Ntk, xag_npn_db_kind::aig_complete> resyn_complete;
  
  std::string info;
  std::vector<std::vector<double>> npn_fractions;
  TT target(4u);
  std::set<TT> reprs;
  int id=0;
  do
  {
    const auto repr = kitty::exact_npn_canonization(target);
   if(kitty::is_const0(target)) 
      reprs.insert( std::get<0>(repr) );
    if( kitty::is_const0(target) || reprs.find( std::get<0>(repr) ) != reprs.end() ) 
    {
      kitty::next_inplace( target );
      continue;
    }
    else
    {
      reprs.insert( std::get<0>(repr) );

      aig_network aig;
      aig_network aig_sat;

      const auto a = aig_sat.create_pi();
      const auto b = aig_sat.create_pi();
      const auto c = aig_sat.create_pi();
      const auto d = aig_sat.create_pi();
      std::vector<aig_network::signal> pis = { a, b, c, d };
      exact_aig_resynthesis<aig_network> resyn( false );
      resyn( aig_sat, target, pis.begin(), pis.end(), [&]( auto const& f ) {
        aig_sat.create_po( f );
      } );

      /* exact */
      info = fmt::format("{:10d} | {:6d} | ", id++, aig_sat.num_gates() );
      Nsat.push_back( aig_sat.num_gates() );
      /* SYM-CUSCO */
      /* COV-CUSCO */
      for( int i{0}; i < 2; ++i )
      {
        report_t<aig_network> report = game_on<aig_network>( &target, MET[i] );
        aig_network aig = report.ntk;

        resubstitution_params ps;
        ps.max_pis = aig.num_pis();
        ps.max_inserts = 20u;
        ps.max_divisors = 1000u;
        ps.odc_levels = -1;
        ps.progress = true;
        auto aig_resub = cleanup_dangling( aig );
        sim_resubstitution( aig_resub, ps );
        aig_resub = cleanup_dangling( aig_resub );
        if( report.nMin < 0 )
          info += fmt::format("{:6s} | {:6s} | {:6s} | ", "-", "-", "-");
        else
          info += fmt::format("{:6d} | {:6d} | {:6.2f} | ", aig.num_gates(), aig_resub.num_gates(), report.time );

        if( i == 0 )
          Nsym.push_back( aig.num_gates() );
        else
          Ncov.push_back( aig.num_gates() );
      }
      printf("%s\n", info.c_str() );
      kitty::next_inplace( target );

    }
  } while ( !kitty::is_const0( target ) );

  printf("\n SAT\n");
  for( int i{0}; i<Nsat.size(); ++i )
    printf("%d, ", Nsat[i] );

  printf("\nSYM\n");
  for( int i{0}; i<Nsym.size(); ++i )
    printf("%d, ", Nsym[i] );

  printf("\nCOV\n");
  for( int i{0}; i<Ncov.size(); ++i )
    printf("%d, ", Ncov[i] );

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

void propose_khot( kitty::dynamic_truth_table * pF, std::string * pInfo )
{
  int nVars = (*pF).num_vars();
  printf( ANSI_COLOR_YELLOW " ENTER THE POPCOUNT [0-%d] " ANSI_COLOR_RESET "", nVars );
  uint32_t POP;
  std::cin >> POP;

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

    if( POP_CNT == POP )
      kitty::set_bit( *pF, j );
  }
  *pInfo += "khot/s" + std::to_string(nVars) + "_" + std::to_string(POP);
}

void propose_thresh( kitty::dynamic_truth_table * pF, std::string * pInfo )
{
  int nVars = (*pF).num_vars();
  printf( ANSI_COLOR_YELLOW " ENTER THE BIAS [0-%d] " ANSI_COLOR_RESET "", nVars+1 );
  int BIAS;
  std::cin >> BIAS;
  kitty::create_threshold( *pF, BIAS );

  *pInfo += "threshold/s" + std::to_string(nVars) + "_" + std::to_string(BIAS);

}

kitty::dynamic_truth_table propose_game( std::string * pInfo )
{
    printf(ANSI_COLOR_YELLOW " ENTER THE NUMBER OF INPUTS: " ANSI_COLOR_RESET "" );
    int INPUT;
    std::cin >> INPUT;
    kitty::dynamic_truth_table F( INPUT );
    printf(ANSI_COLOR_YELLOW " 0 THRESHOLD " ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 1 GAMBLE    " ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 2 k-HOT    " ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 3 PARITY    " ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " CHOOSE THE FUNCTION TYPE: " ANSI_COLOR_RESET "" );
    std::cin >> INPUT;
    //INPUT = 2; printf("2\n");

    switch (INPUT)
    {
    case 0:
        propose_thresh( &F, pInfo );
        break;
    case 1:
        propose_gamble( &F, pInfo );
        break;
    case 2:
        propose_khot( &F, pInfo );
        break;
    case 3:
        propose_parity( &F, pInfo );
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
report_t<Ntk> game_on( kitty::dynamic_truth_table * pF, int MET )
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
    //ntk = solver.aut_sym_solve( MET );
    std::vector<kitty::dynamic_truth_table> xs;
    for( int i{0}; i<pF->num_vars(); ++i )
    {
      xs.emplace_back( pF->num_vars() );
      kitty::create_nth_var( xs[i], i );
    }
    /* define the parameters */
    int nIters = 100;
    cusco_ps ps( ccgame::solver_t::_SYM_ENT, nIters );
    /* solve */
    cusco<Ntk> solver( xs, {*pF} );
    rep = solver.solve( ps );
  }
  else if( MET == 1 )
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
  else if( MET == 2 )
  {
    //ntk = solver.aut_sym_solve( MET );
    std::vector<kitty::dynamic_truth_table> xs;
    for( int i{0}; i<pF->num_vars(); ++i )
    {
      xs.emplace_back( pF->num_vars() );
      kitty::create_nth_var( xs[i], i );
    }
    /* define the parameters */
    int nIters = 100;
    cusco_ps ps( ccgame::solver_t::_SYM_RND, nIters );
    /* solve */
    cusco<Ntk> solver( xs, {*pF} );
    rep = solver.solve( ps );
  }
  else if( MET == 3 )
  {
    std::vector<kitty::dynamic_truth_table> xs;
    for( int i{0}; i<pF->num_vars(); ++i )
    {
      xs.emplace_back( pF->num_vars() );
      kitty::create_nth_var( xs[i], i );
    }

    /* define the parameters */    
    int nIters = 10;
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
  //std::cout << ( kitty::equal( tt, *pF ) ? " equal " : " different " ) << std::endl;


  return rep;

}