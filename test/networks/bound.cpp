#include <catch.hpp>

#include <cstdint>
#include <vector>

#include <lorina/genlib.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/super_reader.hpp>
#include <mockturtle/networks/bound.hpp>
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

TEST_CASE( "Initialize bound network", "[bound]" )
{
  using bound_network = mockturtle::bound_network<2>;
  using signal = bound_network::signal;

  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  bound_network ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const f1 = ntk.create_node( { a, b, c }, { 5 } ); // maj3
  ntk.create_po( f1 );
  CHECK( ntk.is_combinational() );

}
//  CHECK( ntk.is_pi( a.index ) );
//  CHECK( ntk.is_pi( b.index ) );
//  CHECK( ntk.is_pi( c.index ) );
//  CHECK( ntk.is_po( f1 ) );
//  CHECK( ntk.is_constant( ntk.get_constant( false ) ) );
//  CHECK( ntk.is_constant( ntk.get_constant( true ) ) );
//  CHECK( !ntk.is_multioutput( f1 ) );
//  CHECK( ntk.is_multioutput( ntk.get_constant( false ) ) );
//  CHECK( ntk.is_multioutput( ntk.get_constant( true ) ) );
//  CHECK( ntk.is_constant( 0 ) );
//  CHECK( ntk.is_constant( 1 ) );
//  CHECK( ntk.constant_value( 0 ) );
//  CHECK( !ntk.constant_value( 1 ) );
//  CHECK( !ntk.is_ci( 0 ) );
//  CHECK( !ntk.is_ci( 1 ) );
//  CHECK( !ntk.is_ci( f1.index ) );
//  CHECK( !ntk.is_pi( 0 ) );
//  CHECK( !ntk.is_pi( 1 ) );
//  CHECK( !ntk.is_pi( f1.index ) );//