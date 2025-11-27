#include <catch.hpp>

#include <kitty/kitty.hpp>

#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/resynthesis/aig_syn.hpp>
#include <mockturtle/resynthesis/networks/bindings.hpp>
#include <mockturtle/resynthesis/list_simulator.hpp>

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

TEST_CASE( "Initialize the bindings data structure", "[bindings]" )
{
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bindings lib( gates );

  auto const lgates = lib.get_gates();

  using Tt_t = kitty::dynamic_truth_table;
  mockturtle::aig_list_simulator<Tt_t> simulator;
  for( auto g : lgates )
  {
    std::vector<Tt_t> xs;
    std::vector<uint32_t> ids;
    for( uint32_t i = 0; i < g.num_vars; ++i )
    {
      ids.push_back( i );
      xs.emplace_back( g.num_vars );
      kitty::create_nth_var( xs[i], i );
    }
    Tt_t res( g.num_vars );

    simulator( res, g.list, ids.begin(), ids.end(), xs );
    CHECK( kitty::equal( res, g.function ) );
  }
  for ( unsigned int id = 0; id < 10; ++id )
  {
    CHECK( !lib.is_multiple_output( id ) );
  }  
  CHECK( lib.is_multiple_output( 10u ) );
  CHECK( lib.is_multiple_output( 11u ) );
  CHECK( lib.is_multiple_output( 12u ) );
  CHECK( lib.is_multiple_output( 13u ) );
  CHECK( lib.get_group( 10u ) == lib.get_group( 11u ) );
  CHECK( lib.get_group( 12u ) == lib.get_group( 13u ) );

}
