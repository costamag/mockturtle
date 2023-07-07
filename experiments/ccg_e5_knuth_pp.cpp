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

void propose_thresh( kitty::dynamic_truth_table *, std::string * );
void propose_gamble( kitty::dynamic_truth_table *, std::string * );
void propose_khot( kitty::dynamic_truth_table *, std::string * );
void propose_parity( kitty::dynamic_truth_table *, std::string * );

kitty::dynamic_truth_table propose_game(std::string*);
kitty::dynamic_truth_table userdef_game();
kitty::dynamic_truth_table knuth_game( int );

template<class Ntk> Ntk game_on( kitty::dynamic_truth_table *, int );

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


  printf(ANSI_COLOR_YELLOW " KNUTH OR TRUTH TABLE OF YOUR CHOICHE [K/C]? " ANSI_COLOR_RESET "" );
  
  char selection;
  
  std::cin >> selection;
  
  kitty::dynamic_truth_table F;
  if( selection == 'C' || selection == 'c' )
  {
    F = userdef_game();
    int iMet = 0;
    xag_network xag = game_on<xag_network>(&F, iMet);
    printf("%d\n", xag.num_gates());
    //printf(ANSI_COLOR_YELLOW " THE FUNCTION IS " ANSI_COLOR_RESET "" );
  }
  else if( selection == 'K' || selection == 'k' )
  {

    std::string dot_file  ;
    std::string blif_file ;
    std::string aig_file  ;

    std::vector<std::string> sF = {  "S(4)", "S(3)", "S(3,4)", "S(2)", "S(2,4)", "S(2,3)", "S(2,3,4)", 
                                          "S(1)", "S(1,4)", "S(1,3)", "S(1,3,4)", "S(1,2)", 
                                          "S(1,2,4)", "S(1,2,3)", "S(1,2,3,4)", 
                                          "S(4)", "S(4,5)", "S(3)", "S(3,5)", "S(3,4)", "S(3,4,5)", "S(2,5)", "S(2,4)", "S(2,4,5)",
                                        "S(2,3,5)", "S(2,3)", "S(2,3,4)", "S(1,5)", "S(1,4)", "S(1,3,4)", "S(1,2,5)" };
    std::vector<int> CC = { 3, 7, 7, 6, 6, 6, 7, 7, 7, 3, 6, 6, 7, 5, 3,
                                        10, 10, 9, 10, 10, 9, 10, 8, 9, 10, 8, 10, 9, 9, 11, 9 };
    std::vector<std::string> sNames = {  "S4_4", "S4_3", "S4_3_4", "S4_2", "S4_2_4", "S4_2_3", "S4_2_3_4", 
                                          "S4_1", "S4_1_4", "S4_1_3", "S4_1_3_4", "S4_1_2", 
                                          "S4_1_2_4", "S4_1_2_3", "S4_1_2_3_4", 
                                          "S5_4", "S5_4_5", "S5_3", "S5_3_5", "S5_3_4", "S5_3_4_5", "S5_2_5", "S5_2_4", "S5_2_4_5",
                                        "S5_2_3_5", "S5_2_3", "S5_2_3_4", "S5_1_5", "S5_1_4", "S5_1_3_4", "S5_1_2_5" };
    std::vector<std::string> sMet = {"UNINF", "REM1", "REM100", "COV100"};
    printf("%20s| %7s %8s %7s | %7s %8s %7s | %7s %8s %7s | %7s %8s %7s |\n", "", "", "UNINF", "", "", "REM-1", "", "", "REM-100", "", "", "CONV-100", "" );
    printf("%10s | %6s | %6s | %6s | %6s | %6s | %6s | %6s | %6s | %6s | %6s | %6s | %6s | %6s |\n", "f", "C(f)", "XAIG", "XAIG*", "T[s]", "XAIG", "XAIG*", "T[s]", "XAIG", "XAIG*", "T[s]", "XAIG", "XAIG*", "T[s]" );
    for( uint32_t i{0}; i < 31; ++i )
    {
      std::string info = fmt::format("{:10s} | {:6d} | ", sF[i].c_str(), CC[i] );

      for( int iMet{0}; iMet < 1; ++iMet )
      {
        kitty::dynamic_truth_table F;
        F = knuth_game( i );
        std::clock_t start;
        double duration;
        start = std::clock();
        xag_network xag = game_on<xag_network>(&F, iMet);
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
        info += fmt::format("{:6d} | {:6d} | {:6.2f} | ", xag.num_gates(), xag_resub.num_gates(), duration );

        dot_file  = "EXPS/EXP3/" + sMet[iMet] + "/dot/"  + sNames[i] + ".dot";
        blif_file = "EXPS/EXP3/" + sMet[iMet] + "/blif/" + sNames[i] + ".blif";
        aig_file  = "EXPS/EXP3/" + sMet[iMet] + "/aig/"  + sNames[i] + ".aig";
        write_dot( xag, dot_file );
        write_dot( xag, blif_file );
        write_aiger( xag, aig_file );
        dot_file  = "EXPS/EXP3/" + sMet[iMet] + "/dot/"  + sNames[i] + "rs.dot";
        blif_file = "EXPS/EXP3/" + sMet[iMet] + "/blif/" + sNames[i] + "rs.blif";
        aig_file  = "EXPS/EXP3/" + sMet[iMet] + "/aig/"  + sNames[i] + "rs.aig";
        write_dot( xag_resub, dot_file );
        write_dot( xag_resub, blif_file );
        write_aiger( xag_resub, aig_file );
      }
      printf("%s\n", info.c_str() );
    }
  }
  else
  {
    printf("CHOICE NOT VALID\n");
    assert(0);
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
Ntk game_on( kitty::dynamic_truth_table * pF, int MET )
{
  using TT  = kitty::dynamic_truth_table;

  Ntk rep;
  
  kitty::dynamic_truth_table mask = (*pF).construct();
  mask = mask | ~mask;
  //solver.PrintSpecs();
  
  if( MET == 0 )
  {
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
    detailed_gate_t exor_( gate_t::EXOR, 2, 1.0, 1.0, &hpcompute_exor );
    ndps.lib = { ai00_, ai01_, ai10_, ai11_, exor_ };

    mct_ps mctps;
    ndps.sel_type = supp_selection_t::SUP_ENER;
    mctps.nIters =100;
    mctps.nSims = 1;
    mctps.verbose =true;
    ndps.BETA0 = 100;
    ndps.nIters = 100;

    nd_size_t<xag_network> root( xs, ts, {*pF}, ndps );
    mct_method_ps metps;

    mct_method_t<nd_size_t<xag_network> > meth( metps );
    mct_tree_t<nd_size_t<xag_network> , mct_method_t> mct( root, meth, mctps );
    int iSol = mct.solve();
    if( iSol == -1 ) printf("no solution found\n");
    rep = mct.nodes[iSol].ntk;
  }
  else
  {
    printf(ANSI_COLOR_RED " CHOICE NOT MATCHING ANY METHOD " ANSI_COLOR_RESET "\n" );   
    assert(0);
  }

  default_simulator<kitty::dynamic_truth_table> sim( (*pF).num_vars() );
  const auto tt = simulate<kitty::dynamic_truth_table>( rep, sim )[0];
  assert( kitty::equal( tt, *pF ) );


//  std::cout << ( kitty::equal( tt, *pF ) ? " equal " : " different " ) << std::endl;


  return rep;

}