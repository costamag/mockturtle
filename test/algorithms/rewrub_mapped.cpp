#include <catch.hpp>

#include <cstdint>
#include <vector>

#include <lorina/genlib.hpp>
#include <lorina/super.hpp>
#include <mockturtle/algorithms/emap.hpp>
#include <mockturtle/algorithms/rewrub_mapped.hpp>
#include <mockturtle/generators/arithmetic.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/super_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/block.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mapped.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/utils/mapped_database.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/views/cell_view.hpp>
#include <mockturtle/views/dont_touch_view.hpp>

using namespace mockturtle;

std::string const test_library = "GATE   inv1    1 O=!a;            PIN * INV 1 999 0.9 0.3 0.9 0.3\n"
                                 "GATE   inv2    2 O=!a;            PIN * INV 2 999 1.0 0.1 1.0 0.1\n"
                                 "GATE   nand2   2 O=!(a*b);        PIN * INV 1 999 1.0 0.2 1.0 0.2\n"
                                 "GATE   and2    3 O=a*b;           PIN * INV 1 999 1.7 0.2 1.7 0.2\n"
                                 "GATE   xor2    4 O=a^b;           PIN * UNKNOWN 2 999 1.9 0.5 1.9 0.5\n"
                                 "GATE   mig3    3 O=a*b+a*c+b*c;   PIN * INV 1 999 2.0 0.2 2.0 0.2\n"
                                 "GATE   xor3    5 O=a^b^c;         PIN * UNKNOWN 2 999 3.0 0.5 3.0 0.5\n"
                                 "GATE   buf     2 O=a;             PIN * NONINV 1 999 1.0 0.0 1.0 0.0\n"
                                 "GATE   nand3   3 O=!((a*b)*c);        PIN * INV 1 999 1.0 0.2 1.0 0.2\n"
                                 "GATE   zero    0 O=CONST0;\n"
                                 "GATE   one     0 O=CONST1;\n"
                                 "GATE   ha      5 C=a*b;           PIN * INV 1 999 1.7 0.4 1.7 0.4\n"
                                 "GATE   ha      5 S=!a*b+a*!b;     PIN * INV 1 999 2.1 0.4 2.1 0.4\n"
                                 "GATE   fa      6 C=a*b+a*c+b*c;   PIN * INV 1 999 2.1 0.4 2.1 0.4\n"
                                 "GATE   fa      6 S=a^b^c;         PIN * INV 1 999 3.0 0.4 3.0 0.4";

TEST_CASE( "Rewrub tree of nand2", "[rewrub_mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<3> lib( gates );

  mapped_database_params ps;
  mapped_database<tech_library<3>> db( "../experiments/databases/mapped/test4.db", lib, ps );

  mapped_network ntk( gates );
  std::vector<mapped_network::signal> pis( 4 );
  std::generate( std::begin( pis ), std::end( pis ),
                 [&]() { return ntk.create_pi(); } );
  /* constants */

  auto const n0 = ntk.create_node( { pis[0], pis[1] }, 2 );
  auto const n1 = ntk.create_node( { pis[2], pis[3] }, 2 );
  auto const n2 = ntk.create_node( { n0, n1 }, 2 );
  ntk.create_po( n2 );

  rewrub_mapped( ntk, db );

}

TEST_CASE( "Rewrub with structural opportunity", "[rewrub_mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<3> lib( gates );

  mapped_database_params ps;
  mapped_database<tech_library<3>> db( "../experiments/databases/mapped/test4.db", lib, ps );

  mapped_network ntk( gates );
  std::vector<mapped_network::signal> pis( 6 );
  std::generate( std::begin( pis ), std::end( pis ),
                 [&]() { return ntk.create_pi(); } );
  /* constants */

  auto const n1 = ntk.create_node( { pis[0], pis[1] }, 3 );
  auto const n2 = ntk.create_node( { pis[2], pis[3] }, 3 );
  auto const n3 = ntk.create_node( { pis[4], pis[5] }, 3 );
  auto const n4 = ntk.create_node( { n2, n3 }, 3 );
  auto const n5 = ntk.create_node( { n1, n4 }, 2 );
  ntk.create_po( n5 );

  CHECK( ntk.compute_area() == 14 );
  rewrub_mapped( ntk, db );
  CHECK( ntk.compute_area() == 10 );

}
