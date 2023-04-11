#include <kitty/constructors.hpp>
#include <mockturtle/algorithms/decompose/DecSolver.hpp>

using namespace mockturtle;

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

void propose_thresh( kitty::dynamic_truth_table * );
void propose_gamble( kitty::dynamic_truth_table * );
void propose_onehot( kitty::dynamic_truth_table * );
void propose_parity( kitty::dynamic_truth_table * );

kitty::dynamic_truth_table propose_game();
kitty::dynamic_truth_table userdef_game();

void game_on( kitty::dynamic_truth_table * );

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

  printf(ANSI_COLOR_YELLOW " DO YOU WANT ME TO PROPOSE YOU A GAME [Y/N]? " ANSI_COLOR_RESET "" );
  char selection;
  std::cin >> selection;
  kitty::dynamic_truth_table F;
  if( selection == 'Y' || selection == 'y' )
  {
    F = propose_game();
    printf(ANSI_COLOR_YELLOW " THE FUNCTION IS " ANSI_COLOR_RESET "" );
    kitty::print_binary( F ); printf("\n");
  }
  else if( selection == 'N' || selection == 'n' )
  {
    F = userdef_game();
  }
  else
  {
    printf(ANSI_COLOR_RED " CHOICE NOT VALID " ANSI_COLOR_RESET "\n" );
    return 1;
  }

  printf(ANSI_COLOR_YELLOW " GAME TIME! " ANSI_COLOR_RESET "\n" );
  game_on( &F );


  return 0;
}

void propose_gamble( kitty::dynamic_truth_table * pF )
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
}

void propose_parity( kitty::dynamic_truth_table * pF )
{
    kitty::create_parity(*pF);
}

void propose_onehot( kitty::dynamic_truth_table * pF )
{
  int nVars = (*pF).num_vars();
  for( int j{0}; j < nVars; ++j )
  {
    kitty::dynamic_truth_table x(nVars);
    kitty::create_nth_var( x, j );
    kitty::set_bit( *pF, pow(2, j) );
  }
}

void propose_thresh( kitty::dynamic_truth_table * pF )
{
  int nVars = (*pF).num_vars();
  printf( ANSI_COLOR_YELLOW " ENTER THE BIAS [0-%d] " ANSI_COLOR_RESET "", nVars+1 );
  int BIAS;
  std::cin >> BIAS;
  kitty::create_threshold( *pF, BIAS );
}

kitty::dynamic_truth_table propose_game()
{
    printf(ANSI_COLOR_YELLOW " ENTER THE NUMBER OF INPUTS: " ANSI_COLOR_RESET "" );
    int INPUT;
    std::cin >> INPUT;
    kitty::dynamic_truth_table F( INPUT );
    printf(ANSI_COLOR_YELLOW " 0 THRESHOLD " ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 1 GAMBLE    " ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 2 ONEHOT    " ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " 3 PARITY    " ANSI_COLOR_RESET "\n" );
    printf(ANSI_COLOR_YELLOW " CHOOSE THE FUNCTION TYPE: " ANSI_COLOR_RESET "" );
    std::cin >> INPUT;
    switch (INPUT)
    {
    case 0:
        propose_thresh( &F );
        break;
    case 1:
        propose_gamble( &F );
        break;
    case 2:
        propose_onehot( &F );
        break;
    case 3:
        propose_parity( &F );
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
    //printf(ANSI_COLOR_YELLOW " ENTER THE NUMBER OF INPUTS: " ANSI_COLOR_RESET "" );
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

void game_on( kitty::dynamic_truth_table * pF )
{
  using TT  = kitty::dynamic_truth_table;
  using Ntk = aig_network;

  typedef DecSolver<TT, Ntk> solver_t;
  Ntk aig;
  
  kitty::dynamic_truth_table mask = (*pF).construct();
  mask = mask | ~mask;
  solver_t solver( {*pF}, {mask} );
  solver.PrintSpecs();
  solver.man_sym_solve();
}