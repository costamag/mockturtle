#include <kitty/constructors.hpp>
#include <mockturtle/algorithms/mcts/mct_tree.hpp>
#include <mockturtle/algorithms/mcts/method.hpp>
#include <mockturtle/algorithms/mcts/nodes/nd_size.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/write_dot.hpp>

#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/resyn_engines/xag_resyn.hpp>
#include <mockturtle/algorithms/node_resynthesis/exact.hpp>
#include <mockturtle/algorithms/resyn_engines/xag_resyn.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/algorithms/sim_resub.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <string>
#include <bit>
#include <set>
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

  std::vector<double> DELS;
  std::vector<bool> USED;
  std::vector<double> AIGHEU;
  std::vector<double> AIGSAT;

  xag_npn_resynthesis<NTK> resyn;
  xag_npn_resynthesis<NTK, NTK, xag_npn_db_kind::aig_complete> resyn_complete;
  
  std::string info;
  std::vector<std::vector<double>> npn_fractions;
  std::set<TT> reprs;
  int id=0;

  std::vector<int> delta_covering;
  std::vector<int> delta_sym;
  int it{0};
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
      if (kitty::is_const0(target) || kitty::is_const0(~target) )
      {
        continue;
      }
      //kitty::print_binary(target);
      //printf("\n");
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
      AIGSAT.push_back( aig_sat.num_gates() );
        
      result_mctsolve rep = mct_solve( &target );
      AIGHEU.push_back( rep.area );
      printf("%f %f\n", it++, AIGHEU.back(), AIGSAT.back() );
      //printf("======================================\n");

      kitty::next_inplace( target );

    }
  } while ( !kitty::is_const0( target ) );


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
  //detailed_gate_t cmpr_( gate_t::CMPR, 1, 1, 1.0, &hpcompute_cmpr );
  //detailed_gate_t cmpl_( gate_t::CMPR, 1, 1, 1.0, &hpcompute_cmpl );
  detailed_gate_t aig00_( gate_t::AI00, 2, 1, 1.0, &hpcompute_ai00 );
  detailed_gate_t aig01_( gate_t::AI01, 2, 1, 1.0, &hpcompute_ai01 );
  detailed_gate_t aig10_( gate_t::AI10, 2, 1, 1.0, &hpcompute_ai10 );
  detailed_gate_t aig11_( gate_t::AI11, 2, 1, 1.0, &hpcompute_ai11 );
  //detailed_gate_t exor_( gate_t::EXOR, 2, 1, 1.0, &hpcompute_exor );
  ndps.lib = { aig00_, aig01_, aig10_, aig11_ };

  mct_ps mctps;
  ndps.sel_type = supp_selection_t::SUP_BDD;
  mctps.nIters = 10;
  mctps.nSims = 10;
  mctps.verbose =false;
  ndps.BETA0 = 100;
  ndps.BETAZ = 100;
  ndps.nIters = 5;
  ndps.thresh = 6;
  ndps.delay_inv = 0.5;

  nd_size_t<aig_network> root( X, T, {*pF}, ndps );
  mct_method_ps metps;

  mct_method_t<nd_size_t<aig_network> > meth( metps );
  mct_tree_t<nd_size_t<aig_network> , mct_method_t> mct( root, meth, mctps );
  int iSol = mct.solve();
  if( iSol == -1 ) 
  {
    printf("no solution found\n");
    return res;
  }
  else
    res.isValid = true;
  aig_network xag = mct.nodes[iSol].ntk;
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