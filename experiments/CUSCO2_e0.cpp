#include <kitty/constructors.hpp>
#include <mockturtle/algorithms/mcts/supportor.hpp>
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


int main()
{
  printf(ANSI_COLOR_RED     "============================================================="     ANSI_COLOR_RESET "\n");
  printf(ANSI_COLOR_YELLOW   "CUSCO 2     : Generic Set Covering-Based Synthesis " ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW   "EXPERIMENT 0: Effect of Temperature on set covering" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_RED     "============================================================="     ANSI_COLOR_RESET "\n\n");

  using DTT = kitty::dynamic_truth_table;
  std::vector<DTT> xs;  
  std::vector<DTT> fs;

    /* initialize the divisors */
  for( int i{0}; i < 5u; ++i )
  {
      xs.emplace_back( 5u );
      kitty::create_nth_var( xs[i], i );
  }

  std::vector<divisor_t> divisors;
  for( uint32_t i{0}; i < xs.size(); ++i )
  {
      int div_id = i;
      DTT div_tt = xs[i];
      double div_area = 0;
      double div_delay = 0;
      divisor_t div( true, div_id, div_tt, div_area, div_delay );
      divisors.push_back(div);
  }

  kitty::dynamic_truth_table F(5u);
  kitty::create_from_binary_string( F, "01000011101110000110110000100101" );
  std::vector<target_t> targets;
  target_t trg( true, 0, F );
  targets.push_back( trg );

  node_ps ndps;
  detailed_gate_t ai00_( gate_t::AI00, 2, 0.0, 0.0, &hpcompute_ai00 );
  detailed_gate_t ai01_( gate_t::AI01, 2, 0.0, 0.0, &hpcompute_ai01 );
  detailed_gate_t ai10_( gate_t::AI10, 2, 0.0, 0.0, &hpcompute_ai10 );
  detailed_gate_t ai11_( gate_t::AI11, 2, 0.0, 0.0, &hpcompute_ai11 );
  detailed_gate_t exor_( gate_t::EXOR, 2, 0.0, 0.0, &hpcompute_exor );
  ndps.lib = { ai00_, ai10_, ai01_, ai11_, exor_ };
  ndps.nIters = 1;

  ndps.sel_type = supp_selection_t::SUP_ENER;

  for( int order = -10; order < 11; ++order )
  {
    double Beta = pow( 10, order );
    ndps.BETA0 = Beta;
    ndps.BETAZ = Beta;
    ndps.use_inf_graph = true;
    /* support genenrator initialization */
    support_generator_t suppor( &divisors, &targets, ndps );
    int nMax = 30;
        printf("[");
        for( int i{0}; i<nMax; ++i )
        {
            auto sol = suppor.find_new<supp_selection_t::SUP_ENER>(1);
            printf("%d", sol.size());
            if( sol.size() > 0 )
                suppor.store_new(sol);
            if( i < nMax -1 ) printf(", ");
        }
    printf("]\n");
  }

    int maxOrder = 10;
  printf("[");
  for( int order = -10; order < 10; ++order )
  {
    double Beta = pow( 10, order );
    printf("%10.10f", Beta );
    if( order < maxOrder -1 ) printf(", ");
  }
  printf("]");


  return 0;
}
