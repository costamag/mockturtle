#include <catch.hpp>

#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>

#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/resynthesis/networks/bound.hpp>
#include <mockturtle/resynthesis/timing/arrival_analyzer.hpp>
#include <mockturtle/resynthesis/timing/sensing_analyzer.hpp>

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

TEST_CASE( "Compute the arrival time of the gates in a network", "[delay_analyzer]" )
{
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bound_network bound( gates );
  /* network's definition */
  auto const x0 = bound.create_pi();
  auto const x1 = bound.create_pi();
  auto const x2 = bound.create_pi();
  auto const f0 = bound.create_node( { x0, x1 }, "and2" )[0];
  auto const f1 = bound.create_node( { f0, x2 }, "and2" )[0];
  bound.create_po( f1 );
  /* compute delay */
  mockturtle::arrival_analyzer analyzer( bound );

  analyzer.compute_arrival_times( { -1.7, -1.7, 3 } );
  double worst_delay = analyzer.worst_delay();
  CHECK( ( ( worst_delay <= 4.8 ) && ( worst_delay >= 4.6 ) ) );

  analyzer.compute_arrival_times();
  worst_delay = analyzer.worst_delay();
  CHECK( ( ( worst_delay <= 3.5 ) && ( worst_delay >= 3.3 ) ) );


  auto const f2 = bound.create_node( { f1, x2 }, "and2" )[0];
  bound.create_po( f2 );
  worst_delay = analyzer.worst_delay();
  CHECK( ( ( worst_delay <= 5.2 ) && ( worst_delay >= 5.0 ) ) );

  auto const f3 = bound.create_node( { x0, x1 }, "xor2" )[0];
  bound.substitute_node( bound.get_node( f0 ), f3 );
  worst_delay = analyzer.worst_delay();
  CHECK( ( ( worst_delay <= 5.4 ) && ( worst_delay >= 5.2 ) ) );
}

TEST_CASE( "Compute the sensing time of the gates in a network", "[delay_analyzer]" )
{
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bound_network bound( gates );
  /* network's definition */
  auto const x0 = bound.create_pi();
  auto const x1 = bound.create_pi();
  auto const x2 = bound.create_pi();
  auto const f0 = bound.create_node( { x0, x1 }, "and2" )[0];
  auto const f1 = bound.create_node( { f0, x2 }, "and2" )[0];
  bound.create_po( f1 );
  /* compute delay */
  mockturtle::sensing_analyzer analyzer( bound );

  analyzer.compute_sensing_times( );
  CHECK( ( analyzer.get_time( x0 ) == 0 ) );
  CHECK( ( analyzer.get_time( x1 ) == 0 ) );
  CHECK( ( analyzer.get_time( x2 ) == 0 ) );
  CHECK( ( analyzer.get_time( f0 ) == 1.7 ) );
  CHECK( ( analyzer.get_time( f1 ) == 1.7 ) );

}
