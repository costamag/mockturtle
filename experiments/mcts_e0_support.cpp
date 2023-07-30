#include <kitty/constructors.hpp>
#include <kitty/print.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/resyn_engines/xag_resyn.hpp>
#include <mockturtle/algorithms/node_resynthesis/exact.hpp>
#include <mockturtle/algorithms/resyn_engines/xag_resyn.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/algorithms/mcts/mct_utils.hpp>
#include <mockturtle/algorithms/mcts/supportor.hpp>
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

  std::string dot_file  ;
  std::string blif_file ;
  std::string aig_file  ;

  using TT = kitty::dynamic_truth_table;
  using Ntk = aig_network;
  
  // prepare the database for lookup
  TT target(4u);
  std::set<TT> reprs;
  int id=0;

  std::vector<TT> xs;

  std::vector<divisor_t> divs;
  for( int i{0}; i<4u; ++i )
  {
    xs.emplace_back( 4u );
    kitty::create_nth_var( xs[i], i );
    divs.emplace_back( true, i, xs[i], 0.0, 0.0, gate_t::PIS );
  }

    node_ps ndps;
    detailed_gate_t ai00_( gate_t::AI00, 2, 1.0, 1.0, &hpcompute_ai00 );
    detailed_gate_t ai01_( gate_t::AI01, 2, 1.0, 1.0, &hpcompute_ai01 );
    detailed_gate_t ai10_( gate_t::AI10, 2, 1.0, 1.0, &hpcompute_ai10 );
    detailed_gate_t ai11_( gate_t::AI11, 2, 1.0, 1.0, &hpcompute_ai11 );
    detailed_gate_t exor_( gate_t::EXOR, 2, 1.0, 1.0, &hpcompute_exor );
    ndps.lib = { ai00_, ai01_, ai10_, ai11_, exor_ };
    ndps.BETA0 = 1000;
    ndps.BETAZ = 1;
    ndps.nIters = 1;



  do
  {
    const auto repr = kitty::exact_npn_canonization(target);

  std::vector<target_t> trgs;
  trgs.emplace_back( true, 0, target );

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

      support_generator_t supportor( &divs, &trgs, ndps );

      for( int i{0}; i<10; ++i )
      {
        auto S = supportor.find_new<supp_selection_t::SUP_ENER>( 10 );
        printf("%d\n", S.size());
      }

      kitty::next_inplace( target );
      kitty::print_binary( target );
      printf("\n");
    }
  } while ( !kitty::is_const0( target ) );

  return 0;
}
