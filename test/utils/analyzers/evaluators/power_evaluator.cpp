#include <catch.hpp>

#include <cstdint>
#include <vector>

#include <lorina/genlib.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/super_reader.hpp>
#include <mockturtle/networks/mapped/bound_network.hpp>
#include <mockturtle/utils/analyzers/evaluators/power_evaluator.hpp>
#include <mockturtle/utils/tech_library.hpp>

using namespace mockturtle;

std::string const test_library = "GATE   inv1    1 O=!a;            PIN * INV 1 999 0.9 0.3 0.9 0.3\n"
                                 "GATE   inv2    2 O=!a;            PIN * INV 2 999 1.0 0.1 1.0 0.1\n"
                                 "GATE   nand2   2 O=!(a*b);        PIN * INV 1 999 1.0 0.2 1.0 0.2\n"
                                 "GATE   and2    3 O=a*b;           PIN * INV 1 999 1.7 0.2 1.7 0.2\n"
                                 "GATE   xor2    4 O=a^b;           PIN * UNKNOWN 2 999 1.9 0.5 1.9 0.5\n"
                                 "GATE   mig3    3 O=a*b+a*c+b*c;   PIN * INV 1 999 2.0 0.2 2.0 0.2\n"
                                 "GATE   xor3    5 O=a^b^c;         PIN * UNKNOWN 2 999 3.0 0.5 3.0 0.5\n"
                                 "GATE   buf     2 O=a;             PIN * NONINV 1 999 1.0 0.0 1.0 0.0\n"
                                 "GATE   zero    0 O=CONST0;\n"
                                 "GATE   one     0 O=CONST1;\n"
                                 "GATE   ha      5 C=a*b;           PIN * INV 1 999 1.7 0.4 1.7 0.4\n"
                                 "GATE   ha      5 S=!a*b+a*!b;     PIN * INV 1 999 2.1 0.4 2.1 0.4\n"
                                 "GATE   fa      6 C=a*b+a*c+b*c;   PIN * INV 1 999 2.1 0.4 2.1 0.4\n"
                                 "GATE   fa      6 S=a^b^c;         PIN * INV 1 999 3.0 0.4 3.0 0.4";

TEST_CASE( "Power evaluation in Bound networks", "[power_evaluator]" )
{
  using Ntk = mockturtle::bound_network<2>;
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  Ntk ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  /* chain of inverters */
  auto const f1 = ntk.create_node( { a, b }, 4 );
  auto const f2 = ntk.create_node( { a, f1 }, 4 );
  auto const f3 = ntk.create_node( { a, f2 }, 4 );
  ntk.create_po( f3 );

  using TT = kitty::static_truth_table<4u>;
  static constexpr uint32_t num_steps = 10;

  std::vector<TT> tts_init( 2u );
  std::vector<TT> tts_end( 2u );
  kitty::create_nth_var( tts_init[0], 0 );
  kitty::create_nth_var( tts_init[1], 1 );
  kitty::create_nth_var( tts_end[0], 2 );
  kitty::create_nth_var( tts_end[1], 3 );
  workload<TT, num_steps> work( tts_init, tts_end );

  power_evaluator_stats st;
  power_evaluator<Ntk, TT, num_steps> power( ntk, st );

  power.run( work );

  CHECK( st.glitching == 1.5 );
  CHECK( st.switching == 3 );
  CHECK( st.dyn_power == 4 );
  auto sim_eval = power.to_string();

  std::string sim_true =
      "2 0 __________ "
      "-----_____ "
      "__________ "
      "-----_____ "
      "_____----- "
      "---------- "
      "_____----- "
      "---------- "
      "__________ "
      "-----_____ "
      "__________ "
      "-----_____ "
      "_____----- "
      "---------- "
      "_____----- "
      "---------- \n"
      "3 0 __________ "
      "__________ "
      "-----_____ "
      "-----_____ "
      "__________ "
      "__________ "
      "-----_____ "
      "-----_____ "
      "_____----- "
      "_____----- "
      "---------- "
      "---------- "
      "_____----- "
      "_____----- "
      "---------- "
      "---------- \n"
      "4 0 __________ "
      "-------___ "
      "-------___ "
      "__________ "
      "_______--- "
      "---------- "
      "---------- "
      "_______--- "
      "_______--- "
      "---------- "
      "---------- "
      "_______--- "
      "__________ "
      "-------___ "
      "-------___ "
      "__________ \n"
      "5 0 __________ "
      "__________ "
      "-------___ "
      "-------___ "
      "__________ "
      "__________ "
      "-------___ "
      "-------___ "
      "_______--- "
      "_______--- "
      "---------- "
      "---------- "
      "_______--- "
      "_______--- "
      "---------- "
      "---------- \n"
      "6 0 __________ "
      "---------_ "
      "-____---__ "
      "_----___-_ "
      "_________- "
      "---------- "
      "-____---_- "
      "_----___-- "
      "_----___-- "
      "-____---_- "
      "---------- "
      "_________- "
      "_----___-_ "
      "-____---__ "
      "---------_ "
      "__________ \n";

  CHECK( sim_true == sim_eval );
}
