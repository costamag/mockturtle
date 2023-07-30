#include <kitty/constructors.hpp>
#include <mockturtle/algorithms/mcts/supportor.hpp>
#include <string>
#include <bit>
#include <bitset>
#include <cstdint>
#include <fmt/format.h>

#include <iostream>
#include <fstream>

using namespace mockturtle;
using namespace mcts;

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"


std::vector<int> get_supports( kitty::dynamic_truth_table, double, double, int, int, bool );

int main()
{
  printf(ANSI_COLOR_RED     "============================================================="     ANSI_COLOR_RESET "\n");
  printf(ANSI_COLOR_YELLOW  "             Set Covering For Logic Synthesis                " ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_YELLOW  " EXPERIMENT 0: Effect of Temperature on set covering" ANSI_COLOR_RESET "\n" );
  printf(ANSI_COLOR_RED     "============================================================="     ANSI_COLOR_RESET "\n\n");

  std::vector<std::vector<int>> supports_collection_reduced;
  std::vector<std::vector<int>> supports_collection_vanilla;
  std::vector<double> Betas;

  kitty::dynamic_truth_table F(5u);
  kitty::create_from_binary_string( F, "01000011101110000110110000100101" );

  for( int order = -5; order < 6; ++order )
  {
    double Beta = pow( 10, order );
    Betas.push_back(Beta);
    std::vector<int> supports_reduced = get_supports( F, Beta, Beta, 1, 100, true );
    supports_collection_reduced.push_back(supports_reduced);

    std::vector<int> supports_vanilla = get_supports( F, Beta, Beta, 1, 100, false );
    supports_collection_vanilla.push_back(supports_vanilla);
  }

  std::ofstream myfile;
  myfile.open ("../../EXPS/EXP0/HARD_REDUCED.txt");
  for( auto i{0u}; i<Betas.size(); ++i )
  {
    myfile << Betas[i] << " ";
    for( auto j{0u}; j < supports_collection_reduced[i].size(); ++j )
      myfile << supports_collection_reduced[i][j] << " ";
    myfile << "\n";
  }
  myfile.close();

  myfile.open ("../../EXPS/EXP0/HARD_VANILLA.txt");
  for( auto i{0u}; i<Betas.size(); ++i )
  {
    myfile << Betas[i] << " ";
    for( auto j{0u}; j < supports_collection_vanilla[i].size(); ++j )
      myfile << supports_collection_vanilla[i][j] << " ";
    myfile << "\n";
  }
  myfile.close();

  return 0;
}


std::vector<int> get_supports( kitty::dynamic_truth_table F, double Beta0, double BetaZ, int nAttempts, int nMax, bool erase_non_essential )
{
    std::vector<int> res;
    std::vector<kitty::dynamic_truth_table> xs;
    /* initialize the divisors */
  for( uint32_t i{0}; i < 5u; ++i )
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

  ndps.sel_type = supp_selection_t::SUP_NORM;

  ndps.BETA0 = Beta0;
  ndps.BETAZ = BetaZ;
  ndps.nIters = nAttempts;
  ndps.use_inf_graph = true;
  ndps.erase_not_essentials = erase_non_essential;

  /* support genenrator initialization */
  support_generator_t suppor( &divisors, &targets, ndps );
  for( int i{0}; i<nMax; ++i )
  {
      auto sol = suppor.find_new<supp_selection_t::SUP_NORM>(nAttempts);
      res.push_back(sol.size());
      if( sol.size() > 0 )
          suppor.store_new(sol);
  }
  return res;
}



