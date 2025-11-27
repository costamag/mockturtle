#include <catch.hpp>

#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>

#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/resynthesis/aig_syn.hpp>
#include <mockturtle/resynthesis/networks/bindings.hpp>
#include <mockturtle/resynthesis/networks/bound.hpp>
#include <mockturtle/resynthesis/list_simulator.hpp>

std::string const test_library = "GATE   inv1    1 O=!a;            PIN * INV 1 999 0.9 0.3 0.9 0.3\n"
                                 "GATE   inv2    2 O=!a;            PIN * INV 2 999 1.0 0.1 1.0 0.1\n"
                                 "GATE   nand2   2 O=!(a*b);        PIN * INV 1 999 1.0 0.2 1.0 0.2\n"
                                 "GATE   and2    3 O=a*b;           PIN * INV 1 999 1.7 0.2 1.7 0.2\n"
                                 "GATE   or2     4 O=a+b;           PIN * UNKNOWN 2 999 1.9 0.5 1.9 0.5\n"
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

// TODO: adjust multi-output
//TEST_CASE( "Storage initialization of a bound network", "[bound_network]" )
//{
//  std::vector<mockturtle::gate> gates;
//
//  std::istringstream in( test_library );
//  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
//  CHECK( result == lorina::return_code::success );
//
//  mockturtle::bound_network bound( gates );
//  /* get the gates with the synthesized list */
//  auto const lib = bound.get_bindings();
//  auto const lgates = lib.get_gates();
//
//  using Tt_t = kitty::dynamic_truth_table;
//  mockturtle::aig_list_simulator<Tt_t> simulator;
//  for( auto g : lgates )
//  {
//    std::vector<Tt_t> xs;
//    std::vector<uint32_t> ids;
//    for( uint32_t i = 0; i < g.num_vars; ++i )
//    {
//      ids.push_back( i );
//      xs.emplace_back( g.num_vars );
//      kitty::create_nth_var( xs[i], i );
//    }
//    Tt_t res( g.num_vars );
//
//    simulator( res, g.list, ids.begin(), ids.end(), xs );
//    CHECK( kitty::equal( res, g.function ) );
//  }
//  for ( unsigned int id = 0; id < 10; ++id )
//  {
//    CHECK( !lib.is_multiple_output( id ) );
//  }  
//  CHECK( lib.is_multiple_output( 10u ) );
//  CHECK( lib.is_multiple_output( 11u ) );
//  CHECK( lib.is_multiple_output( 12u ) );
//  CHECK( lib.is_multiple_output( 13u ) );
//  CHECK( lib.get_group( 10u ) == lib.get_group( 11u ) );
//  CHECK( lib.get_group( 12u ) == lib.get_group( 13u ) );
//}

TEST_CASE( "Inputs, outputs, and constants", "[bound_network]" )
{
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bound_network bound( gates );
  /* constants storage */
  mockturtle::bound_network::signal const0( static_cast<uint32_t>( 0 ) );
  mockturtle::bound_network::signal const1( static_cast<uint32_t>( 1 ) );
  CHECK( ( bound.get_constant( true ) == const1 ) );
  CHECK( ( bound.get_constant( false ) == const0 ) );
  CHECK( ( bound.constant_value( 0 ) == false ) );
  CHECK( ( bound.constant_value( 1 ) == true ) );
  /* PI construction */
  auto const x0 = bound.create_pi();
  auto const x1 = bound.create_pi();
  auto const n0 = bound.get_node( x0 );
  auto const n1 = bound.get_node( x1 );
  CHECK( bound.is_ci( n0 ) );
  CHECK( bound.is_ci( n1 ) );
  CHECK( bound.is_pi( n0 ) );
  CHECK( bound.is_pi( n1 ) );
  /* PO construction */
  bound.create_po( n1 );
  CHECK( bound.fanout_size( n1 ) == 1 );
  CHECK( bound.fanout_size( n0 ) == 0 );
}

TEST_CASE( "Simulate with constant inputs", "[bound_network]" )
{
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bound_network bound( gates );
  /* constants storage */
  mockturtle::bound_network::signal const0( static_cast<uint32_t>( 0 ) );
  mockturtle::bound_network::signal const1( static_cast<uint32_t>( 1 ) );
  /* PI construction */
  auto const x0 = bound.create_pi();
  auto const a = bound.create_node( { x0, const0 }, "or2" )[0];

  std::vector<kitty::static_truth_table<2u>> tts;
  tts.emplace_back();
  tts.emplace_back();
  kitty::create_nth_var( tts[0], 0 );
  kitty::create_nth_var( tts[1], 1 );
  auto res = bound.compute( a.index, tts.begin(), tts.end() );
  CHECK( kitty::equal( res, tts[0] ) );
  kitty::print_binary( res ); std::cout << std::endl;
}

TEST_CASE( "Create nodes in the network", "[bound_network]" )
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
  auto const x3 = bound.create_pi();

  using signal = typename mockturtle::bound_network::signal;
  signal const f0 = bound.create_node( std::vector<signal>( { x0, x1 } ), "and2" )[0];
  signal const f1 = bound.create_node( std::vector<signal>( { x2, x3 } ), "and2" )[0];
  signal const f2 = bound.create_node( std::vector<signal>( { f0, f1 } ), "and2" )[0];
  bound.create_po( f2 );

  std::vector<kitty::dynamic_truth_table> xs;
  for( int i = 0; i < 4; ++i )
  {
    xs.emplace_back( 4 );
    kitty::create_nth_var( xs[i], i );
  }
  auto tt = bound.compute( bound.get_node( f0 ), xs.begin(), xs.begin() + 2 );

  CHECK( kitty::equal( tt, xs[0] & xs[1] ) );
}

TEST_CASE( "Replace a node in the network", "[bound_network]" )
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
  auto const x3 = bound.create_pi();

  using signal = typename mockturtle::bound_network::signal;
  signal const f0 = bound.create_node( std::vector<signal>( { x0, x1 } ), "and2" )[0];
  signal const r0 = bound.create_node( std::vector<signal>( { x0, x1 } ), "xor2" )[0];
  signal const f1 = bound.create_node( std::vector<signal>( { x2, x3 } ), "and2" )[0];
  signal const f2 = bound.create_node( std::vector<signal>( { f0, f1 } ), "and2" )[0];
  bound.create_po( f2 );

  std::vector<kitty::dynamic_truth_table> xs;
  for( int i = 0; i < 4; ++i )
  {
    xs.emplace_back( 4 );
    kitty::create_nth_var( xs[i], i );
  }
  auto tt1 = bound.compute( bound.get_node( f0 ), xs.begin(), xs.begin() + 2 );

  CHECK( kitty::equal( tt1, xs[0] & xs[1] ) );
  auto tt2 = bound.compute( bound.get_node( r0 ), xs.begin(), xs.begin() + 2 );
  bound.substitute_node( bound.get_node( f0 ), r0 );

  CHECK( kitty::equal( tt2, xs[0] ^ xs[1] ) );
}

TEST_CASE( "Replace an equivalent node in the network", "[bound_network]" )
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
  auto const x3 = bound.create_pi();

  using signal = typename mockturtle::bound_network::signal;
  signal const r0 = bound.create_node( std::vector<signal>( { x0, x1 } ), "and2" )[0];
  signal const f0 = bound.create_node( std::vector<signal>( { x0, x1 } ), "nand2" )[0];
  signal const f1 = bound.create_node( std::vector<signal>( { f0 } ), "inv2" )[0];
  signal const f2 = bound.create_node( std::vector<signal>( { x2, x3 } ), "and2" )[0];
  bound.create_po( f2 );
  double const area1 = bound.area();
  const auto tt1 = mockturtle::simulate<kitty::static_truth_table<5u>>( bound )[0];
  bound.substitute_node( bound.get_node( f1 ), r0 );
  double const area2 = bound.area();
  const auto tt2 = mockturtle::simulate<kitty::static_truth_table<5u>>( bound )[0];
  CHECK( kitty::equal( tt1, tt2 ) );
  CHECK( area2 == area1 );
}

TEST_CASE( "Replace an output node", "[bound_network]" )
{
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bound_network bound( gates );
  /* network's definition */
  auto const x0 = bound.create_pi();
  auto const x1 = bound.create_pi();

  using signal = typename mockturtle::bound_network::signal;
  signal const f0 = bound.create_node( std::vector<signal>( { x0 } ), "inv1" )[0];
  signal const f1 = bound.create_node( std::vector<signal>( { x1 } ), "inv1" )[0];
  signal const f2 = bound.create_node( std::vector<signal>( { f0, f1 } ), "and2" )[0];
  signal const f3 = bound.create_node( std::vector<signal>( { f2 } ), "inv2" )[0];
  signal const f4 = bound.create_node( std::vector<signal>( { f0, f1 } ), "nand2" )[0];

  bound.create_po( f3 );
  double const area1 = bound.area();
  const auto tt1 = mockturtle::simulate<kitty::static_truth_table<5u>>( bound )[0];
  bound.substitute_node( bound.get_node( f3 ), f4 );
  double const area2 = bound.area();
  const auto tt2 = mockturtle::simulate<kitty::static_truth_table<5u>>( bound )[0];
  CHECK( kitty::equal( tt1, tt2 ) );
  CHECK( area2 < area1 );
}

TEST_CASE( "Replace with ripple effect", "[bound_network]" )
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

  using signal = typename mockturtle::bound_network::signal;
  signal const f0 = bound.create_node( std::vector<signal>( { x0 } ), "inv1" )[0];
  signal const f1 = bound.create_node( std::vector<signal>( { x1 } ), "inv1" )[0];
  signal const f2 = bound.create_node( std::vector<signal>( { x0, f1 } ), "and2" )[0];
  signal const f3 = bound.create_node( std::vector<signal>( { f0, x1 } ), "and2" )[0];
  signal const f4 = bound.create_node( std::vector<signal>( { f2, f3 } ), "or2" )[0];
  signal const r0 = bound.create_node( std::vector<signal>( { x0, x1 } ), "xor2" )[0];
  signal const f5 = bound.create_node( std::vector<signal>( { x2, f4 } ), "and2" )[0];

  bound.create_po( f5 );
  double const area1 = bound.area();
  const auto tt1 = mockturtle::simulate<kitty::static_truth_table<5u>>( bound )[0];
  bound.substitute_node( bound.get_node( f4 ), r0 );
  double const area2 = bound.area();
  const auto tt2 = mockturtle::simulate<kitty::static_truth_table<5u>>( bound )[0];
  CHECK( kitty::equal( tt1, tt2 ) );
  CHECK( area2 < area1 );
}