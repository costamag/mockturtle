#include <kitty/constructors.hpp>
#include <mockturtle/algorithms/decompose/DecSolver.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/write_dot.hpp>
#include <string>
#include <bit>
#include <bitset>
#include <cstdint>



using namespace mockturtle;

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
kitty::dynamic_truth_table knuth_game(std::string*);

template<class Ntk> Ntk game_on( kitty::dynamic_truth_table * );

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

  char selection;
  std::cin >> selection;

  kitty::dynamic_truth_table F;
  if( selection == 'Y' || selection == 'y' )
  {
    F = propose_game(&info);
    //printf(ANSI_COLOR_YELLOW " THE FUNCTION IS " ANSI_COLOR_RESET "" );
  }
  else if( selection == 'N' || selection == 'n' )
  {
    F = userdef_game();
  }
  else if( selection == 'K' || selection == 'k' )
  {
    F = knuth_game(&info);
  }
  else
  {
    printf(ANSI_COLOR_RED " CHOICE NOT VALID " ANSI_COLOR_RESET "\n" );
    return 1;
  }

  printf(ANSI_COLOR_YELLOW " 0 XAG " ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " 1 AIG " ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " CHOOSE A NETWORK TYPE: " ANSI_COLOR_RESET "" );
  int NTK;
  
  std::cin >> NTK;

  if( NTK == 0 )
  {
    std::string dot_file = "EXPS/XAG/dot/" + info + ".dot";
    std::string aig_file = "EXPS/XAG/aig/" + info + ".aig";
    xag_network xag = game_on<xag_network>(&F);
    write_dot( xag, dot_file );
    write_aiger( xag, aig_file );
  }
  else if( NTK == 1 )
  {
    std::string dot_file = "EXPS/AIG/dot/" + info + ".dot";
    std::string aig_file = "EXPS/AIG/aig/" + info + ".aig";
    aig_network aig = game_on<aig_network>(&F);
    write_dot( aig, dot_file );
    write_aiger( aig, aig_file );
  }
  else
  {
    printf(ANSI_COLOR_RED " NETWORK TYPE NOT VALID " ANSI_COLOR_RESET "\n" );
    assert(0);
  }

  printf(ANSI_COLOR_YELLOW " GAME TIME! " ANSI_COLOR_RESET "\n" );
  
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

void propose_symmetric( kitty::dynamic_truth_table * pF, std::vector<uint32_t> vals, std::string * pInfo )
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
  *pInfo += "sym/s" + std::to_string(nVars) + "_"; 
  for( auto v : vals )
    *pInfo += std::to_string( v ) + "_";
  
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

kitty::dynamic_truth_table knuth_game( std::string * pInfo )
{
    int INPUT;
    int nVars{0};

    printf(ANSI_COLOR_YELLOW " 0  n=4 S(4)        C(f)=3 <= 3" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 1  n=4 S(3)        C(f)=7 <= 7" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 2  n=4 S(3,4)      C(f)=7 <= 7" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 3  n=4 S(2)        C(f)=6 <= 7" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 4  n=4 S(2,4)      C(f)=6 <= 7" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 5  n=4 S(2,3)      C(f)=6 <= 9" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 6  n=4 S(2,3,4)    C(f)=7 <= 7" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 7  n=4 S(1)        C(f)=7 <= 7" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 8  n=4 S(1,4)      C(f)=7 <= 9" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 9  n=4 S(1,3)      C(f)=3 <= 3" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 10 n=4 S(1,3,4)    C(f)=6 <= 7" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 11 n=4 S(1,2)      C(f)=6 <= 9" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 12 n=4 S(1,2,4)    C(f)=7 <= 9" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 13 n=4 S(1,2,3)    C(f)=5 <= 7" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 14 n=4 S(1,2,3,4)  C(f)=3 <= 3" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " =================" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 15 n=5 S(4)        C(f)=10 <= 10" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 16 n=5 S(4,5)      C(f)=10 <= 10" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 17 n=5 S(3)        C(f)=9  <= 12" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 18 n=5 S(3,5)      C(f)=10 <= 10" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 19 n=5 S(3,4)      C(f)=10 <= 13" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 20 n=5 S(3,4,5)    C(f)=9  <= 10" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 21 n=5 S(2,5)      C(f)=10 <= 14" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 22 n=5 S(2,4)      C(f)=8  <= 10" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 23 n=5 S(2,4,5)    C(f)=9  <= 12" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 24 n=5 S(2,3,5)    C(f)=10 <= 15" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 25 n=5 S(2,3)      C(f)=8  <= 15" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 26 n=5 S(2,3,4)    C(f)=10 <= 13" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 27 n=5 S(1,5)      C(f)=9  <= 13" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 28 n=5 S(1,4)      C(f)=9  <= 15" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 29 n=5 S(1,3,4)    C(f)=11 <= 13" ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 30 n=5 S(1,2,5)    C(f)=9  <= 15" ANSI_COLOR_RESET "\n" );


    printf(ANSI_COLOR_YELLOW " CHOOSE THE FUNCTION TYPE: " ANSI_COLOR_RESET "" );
    std::cin >> INPUT;
    if( INPUT <= 14 )
      nVars = 4;
    else if( INPUT < 31 )
      nVars = 5;
    else
    {
      printf(ANSI_COLOR_RED " PROBLEM NOT DEFINED BY KNUTH" ANSI_COLOR_RESET "\n" );
      assert(0);
    }

    kitty::dynamic_truth_table F( nVars );

    std::vector<uint32_t> VALS;
    switch (INPUT)
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

    case 15: VALS = {4}; break;
    case 16: VALS = {4,5}; break;
    case 17: VALS = {3}; break;
    case 18: VALS = {3,5}; break;
    case 19: VALS = {3,4}; break;
    case 20: VALS = {3,4,5}; break;
    case 21: VALS = {2,5}; break;
    case 22: VALS = {2,4}; break;
    case 23: VALS = {2,4,5}; break;
    case 24: VALS = {2,3,5}; break;
    case 25: VALS = {2,3}; break;
    case 26: VALS = {2,3,4}; break;
    case 27: VALS = {1,5}; break;
    case 28: VALS = {1,4}; break;
    case 29: VALS = {1,3,4}; break;
    case 30: VALS = {1,2,5}; break;

    default:

        break;
    }

    propose_symmetric( &F, VALS, pInfo );

    return F;
}

template<class Ntk>
Ntk game_on( kitty::dynamic_truth_table * pF )
{
  using TT  = kitty::dynamic_truth_table;

  typedef DecSolver<TT, Ntk> solver_t;
  Ntk ntk;
  
  kitty::dynamic_truth_table mask = (*pF).construct();
  mask = mask | ~mask;
  solver_t solver( {*pF}, {mask} );
  solver.PrintSpecs();

  printf(ANSI_COLOR_YELLOW " 0 SYM MANUAL" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " 1 DEC MANUAL" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " 2 SYM AUTOMATIC" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " 3 DEC AUTOMATIC" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " 4 DEC AUTOMATIC WEAK" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " 5 SYM MANUAL RS" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " 6 SYM AUTOMATIC RS" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW " CHOOSE YOUR METHOD: " ANSI_COLOR_RESET "" );
  int MET;
  std::cin >> MET;
  
  if( MET == 0 )
    ntk = solver.man_sym_solve();
  else if( MET == 1 )
    ntk = solver.man_rdec_solve();
  else if( MET == 2 )
  {
    printf(ANSI_COLOR_YELLOW " NUMBER OF ITERATIONS: " ANSI_COLOR_RESET "" );
    std::cin >> MET;
    MET = 30;
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
  else
  {
    printf(ANSI_COLOR_RED " CHOICE NOT MATCHING ANY METHOD " ANSI_COLOR_RESET "\n" );   
    assert(0);
  }

  default_simulator<kitty::dynamic_truth_table> sim( (*pF).num_vars() );
  const auto tt = simulate<kitty::dynamic_truth_table>( ntk, sim )[0];
  if( tt.num_vars() < 10 )
  {
    printf("\n simulation\n");
    kitty::print_binary( tt );
    std::cout << std::endl;
    printf("function returned\n");
    kitty::print_binary( *pF );
    std::cout << std::endl;
  }

  std::cout << ( kitty::equal( tt, *pF ) ? " equal " : " different " ) << std::endl;

  return ntk;

}