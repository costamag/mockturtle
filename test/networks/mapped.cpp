#include <catch.hpp>

#include <cstdint>
#include <vector>

#include <lorina/genlib.hpp>
#include <lorina/super.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/emap.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/generators/arithmetic.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/super_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/block.hpp>
#include <mockturtle/networks/mapped.hpp>
#include <mockturtle/utils/tech_library.hpp>
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
                                 "GATE   zero    0 O=CONST0;\n"
                                 "GATE   one     0 O=CONST1;\n"
                                 "GATE   ha      5 C=a*b;           PIN * INV 1 999 1.7 0.4 1.7 0.4\n"
                                 "GATE   ha      5 S=!a*b+a*!b;     PIN * INV 1 999 2.1 0.4 2.1 0.4\n"
                                 "GATE   fa      6 C=a*b+a*c+b*c;   PIN * INV 1 999 2.1 0.4 2.1 0.4\n"
                                 "GATE   fa      6 S=a^b^c;         PIN * INV 1 999 3.0 0.4 3.0 0.4";

TEST_CASE( "Definition of a mapped network", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );
  CHECK( ntk._storage->data.aig.num_pos() == 14u );
  CHECK( ntk._storage->data.library.size() == 14u );
  
  for( size_t i{0}; i<gates.size(); ++i )
  {
    CHECK( ntk._storage->data.library[i].id == gates[i].id );
    CHECK( ntk._storage->data.library[i].expression == gates[i].expression );
  }
}

TEST_CASE( "create and use constants in a mapped network", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );

  CHECK( has_size_v<mapped_network> );
  CHECK( has_get_constant_v<mapped_network> );
  CHECK( has_is_constant_v<mapped_network> );
  CHECK( has_is_pi_v<mapped_network> );
  CHECK( has_is_constant_v<mapped_network> );
  CHECK( has_get_node_v<mapped_network> );
  CHECK( has_is_complemented_v<mapped_network> );

  CHECK( ntk.size() == 2 );

  auto c0 = ntk.get_constant( false );
  auto c1 = ntk.get_constant( true );

  CHECK( ntk.size() == 2 );
  CHECK( c0 != c1 );
  CHECK( ntk.get_node( c0 ) == 0 );
  CHECK( ntk.get_node( c1 ) == 1 );
  CHECK( !ntk.is_complemented( c0 ) );
  CHECK( !ntk.is_complemented( c1 ) );
  CHECK( ntk.is_constant( c0 ) );
  CHECK( ntk.is_constant( c1 ) );
  CHECK( !ntk.is_pi( c0 ) );
  CHECK( !ntk.is_pi( c1 ) );
}

TEST_CASE( "create and use primary inputs in a mapped network", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );

  CHECK( has_create_pi_v<mapped_network> );
  CHECK( has_is_constant_v<mapped_network> );
  CHECK( has_is_pi_v<mapped_network> );
  CHECK( has_num_pis_v<mapped_network> );

  CHECK( ntk.num_pis() == 0 );

  auto x1 = ntk.create_pi();
  auto x2 = ntk.create_pi();

  CHECK( ntk.size() == 4 );
  CHECK( ntk.num_pis() == 2 );
  CHECK( x1 != x2 );
}

TEST_CASE( "create and use primary outputs in a mapped network", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );

  CHECK( has_create_po_v<mapped_network> );
  CHECK( has_num_pos_v<mapped_network> );

  auto c0 = ntk.get_constant( false );
  auto c1 = ntk.get_constant( true );
  auto x = ntk.create_pi();

  ntk.create_po( c0 );
  ntk.create_po( c1 );
  ntk.create_po( x );

  CHECK( ntk.size() == 3 );
  CHECK( ntk.num_pis() == 1 );
  CHECK( ntk.num_pos() == 3 );
}

TEST_CASE( "create unary operations in a mapped network", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );

  CHECK( has_create_buf_v<mapped_network> );
  CHECK( has_create_not_v<mapped_network> );

  auto x1 = ntk.create_pi();

  CHECK( ntk.size() == 3 );

  auto f1 = ntk.create_buf( x1 );
  auto f2 = ntk.create_not( x1 );

  CHECK( ntk.size() == 5 );
  CHECK( f1 != x1 );
  CHECK( f2 != x1 );
}

TEST_CASE( "create binary operations in a mapped network", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );

  const auto x1 = ntk.create_pi();
  const auto x2 = ntk.create_pi();

  CHECK( ntk.size() == 4 );

  ntk.create_node( { x1, x2 }, 3 );
  CHECK( ntk.size() == 5 );

  ntk.create_node( { x1, x2 }, 3 );
  CHECK( ntk.size() == 5 );

  ntk.create_node( { x2, x1 }, 3 );
  CHECK( ntk.size() == 6 );
}

TEST_CASE( "clone a mapped network", "[mapped]" )
{
  CHECK( has_clone_v<mapped_network> );

  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk1( gates );

  auto a = ntk1.create_pi();
  auto b = ntk1.create_pi();
  auto f1 = ntk1.create_node( { a, b }, 3 );
  ntk1.create_po( f1 );
  CHECK( ntk1.size() == 5 );
  CHECK( ntk1.num_gates() == 1 );
  CHECK( ntk1.num_pos() == 1 );

  auto ntk2 = ntk1;
  auto ntk3 = ntk1.clone();

  auto c = ntk2.create_pi();
  auto f2 = ntk2.create_node( { f1, c }, 2 );
  ntk2.create_po( f2 );
  CHECK( ntk1.size() == 7 );
  CHECK( ntk1.num_gates() == 2 );
  CHECK( ntk1.num_pos() == 2 );

  CHECK( ntk3.size() == 5 );
  CHECK( ntk3.num_gates() == 1 );
  CHECK( ntk3.num_pos() == 1 );
}

TEST_CASE( "clone a node in a mapped network", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk1( gates ), ntk2( gates );

  CHECK( has_clone_node_v<mapped_network> );

  auto a1 = ntk1.create_pi();
  auto b1 = ntk1.create_pi();
  auto f1 = ntk1.create_node( { a1, b1 }, 3 );
  CHECK( ntk1.size() == 5 );

  auto a2 = ntk2.create_pi();
  auto b2 = ntk2.create_pi();
  CHECK( ntk2.size() == 4 );

  auto f2 = ntk2.clone_node( ntk1, ntk1.get_node( f1 ), { a2, b2 } );
  CHECK( ntk2.size() == 5 );

  ntk2.foreach_fanin( ntk2.get_node( f2 ), [&]( auto const& s ) {
    CHECK( !ntk2.is_complemented( s ) );
  } );
}

TEST_CASE( "compute functions from AND and NOT gates in mapped networks", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );

  const auto a = ntk.create_pi();
  const auto b = ntk.create_pi();

  const auto f1 = ntk.create_not( a );
  const auto f2 = ntk.create_node( { a, b }, 3 );

  std::vector<kitty::dynamic_truth_table> xs;
  xs.emplace_back( 3u );
  xs.emplace_back( 3u );
  kitty::create_nth_var( xs[0], 0 );
  kitty::create_nth_var( xs[1], 1 );

  const auto sim_f1 = ntk.compute( ntk.get_node( f1 ), xs.begin(), xs.begin() + 1 );
  const auto sim_f2 = ntk.compute( ntk.get_node( f2 ), xs.begin(), xs.end() );

  CHECK( sim_f1 == ~xs[0] );
  CHECK( sim_f2 == ( xs[0] & xs[1] ) );
}

TEST_CASE( "create nodes and compute a function in a mapped network", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );

  CHECK( has_compute_v<mapped_network, kitty::dynamic_truth_table> );

  const auto a = ntk.create_pi();
  const auto b = ntk.create_pi();
  const auto c = ntk.create_pi();

  kitty::dynamic_truth_table tt_maj( 3u ), tt_xor( 3u ), tt_const0( 0u );
  kitty::create_from_hex_string( tt_maj, "e8" );
  kitty::create_from_hex_string( tt_xor, "96" );

  CHECK( ntk.size() == 5 );

  const auto _const0 = ntk.create_node( {}, 8 );
  const auto _const1 = ntk.create_node( {}, 9 );
  CHECK( _const0 == ntk.get_constant( false ) );
  CHECK( _const1 == ntk.get_constant( true ) );

  const auto _maj = ntk.create_node( { a, b, c }, 5 );
  const auto _xor = ntk.create_node( { a, b, c }, 6 );

  CHECK( ntk.size() == 7 );

  std::vector<kitty::dynamic_truth_table> xs;
  xs.emplace_back( 3u );
  xs.emplace_back( 3u );
  xs.emplace_back( 3u );
  kitty::create_nth_var( xs[0], 0 );
  kitty::create_nth_var( xs[1], 1 );
  kitty::create_nth_var( xs[2], 2 );

  const auto sim_maj = ntk.compute( ntk.get_node( _maj ), xs.begin(), xs.end() );
  const auto sim_xor = ntk.compute( ntk.get_node( _xor ), xs.begin(), xs.end() );

  CHECK( sim_maj == kitty::ternary_majority( xs[0], xs[1], xs[2] ) );
  CHECK( sim_xor == ( xs[0] ^ xs[1] ^ xs[2] ) );

}

TEST_CASE( "hash nodes in mapped network", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );

  const auto a = ntk.create_pi();
  const auto b = ntk.create_pi();
  const auto c = ntk.create_pi();

  kitty::dynamic_truth_table tt_maj( 3u ), tt_xor( 3u );
  kitty::create_from_hex_string( tt_maj, "e8" );
  kitty::create_from_hex_string( tt_xor, "96" );

  ntk.create_node( { a, b, c }, 5 );
  ntk.create_node( { a, b, c }, 6 );

  CHECK( ntk.size() == 7 );

  ntk.create_node( { a, b, c }, 5 );

  CHECK( ntk.size() == 7 );
}

TEST_CASE( "substitute node by another in mapped", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );

  const auto c0 = ntk.get_node( ntk.get_constant( false ) );
  const auto c1 = ntk.get_node( ntk.get_constant( true ) );
  const auto a = ntk.create_pi();
  const auto b = ntk.create_pi();

  kitty::dynamic_truth_table tt_nand( 2u ), tt_le( 2u ), tt_ge( 2u ), tt_or( 2u );
  kitty::create_from_hex_string( tt_nand, "7" );
  kitty::create_from_hex_string( tt_le, "2" );
  kitty::create_from_hex_string( tt_ge, "4" );
  kitty::create_from_hex_string( tt_or, "e" );

  // XOR with NAND
  const auto n1 = ntk.create_node( { a, b }, 2 );
  const auto n2 = ntk.create_node( { a, n1 }, 2 );
  const auto n3 = ntk.create_node( { b, n1 }, 2 );
  const auto n4 = ntk.create_node( { n2, n3 }, 2 );
  ntk.create_po( n4 );

  std::vector<node<mapped_network>> nodes;
  ntk.foreach_node( [&]( auto node ) { nodes.push_back( node ); } );

  CHECK( nodes == std::vector<mapped_network::node>{ c0, c1, a, b, n1, n2, n3, n4 } );
  CHECK( ntk.fanout_size( n4 ) == 1 );
  ntk.foreach_po( [&]( auto f ) {
    CHECK( f == n4 );
    return false;
  } );

  // XOR with AND and OR
  const auto ia = ntk.create_node( { a }, 0 );
  const auto ib = ntk.create_node( { b }, 0 );
  const auto n5 = ntk.create_node( { ia, b }, 2 );
  const auto n6 = ntk.create_node( { a, ib }, 2 );
  const auto n7 = ntk.create_node( { n5, n6 }, 2 );

  nodes.clear();
  ntk.foreach_node( [&]( auto node ) { nodes.push_back( node ); } );

  CHECK( nodes == std::vector{ c0, c1, a, b, n1, n2, n3, n4, ia, ib, n5, n6, n7 } );
  CHECK( ntk.fanout_size( n7 ) == 0 );

  // substitute nodes
  ntk.substitute_node( n4, n7 );

  CHECK( ntk.size() == 13 );
  CHECK( ntk.fanout_size( n4 ) == 0 );
  CHECK( ntk.fanout_size( n7 ) == 1 );
  ntk.foreach_po( [&]( auto f ) {
    CHECK( f == n7 );
    return false;
  } );
}

TEST_CASE( "structural properties of a mapped network", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );

  CHECK( has_size_v<mapped_network> );
  CHECK( has_num_pis_v<mapped_network> );
  CHECK( has_num_pos_v<mapped_network> );
  CHECK( has_num_gates_v<mapped_network> );
  CHECK( has_fanin_size_v<mapped_network> );
  CHECK( has_fanout_size_v<mapped_network> );

  const auto x1 = ntk.create_pi();
  const auto x2 = ntk.create_pi();

  const auto f1 = ntk.create_node( { x1, x2 }, 3 );
  const auto f2 = ntk.create_node( { x2, x1 }, 3 );

  ntk.create_po( f1 );
  ntk.create_po( f2 );

  CHECK( ntk.size() == 6 );
  CHECK( ntk.num_pis() == 2 );
  CHECK( ntk.num_pos() == 2 );
  CHECK( ntk.num_gates() == 2 );
  CHECK( ntk.fanin_size( ntk.get_node( x1 ) ) == 0 );
  CHECK( ntk.fanin_size( ntk.get_node( x2 ) ) == 0 );
  CHECK( ntk.fanin_size( ntk.get_node( f1 ) ) == 2 );
  CHECK( ntk.fanin_size( ntk.get_node( f2 ) ) == 2 );
  CHECK( ntk.fanout_size( ntk.get_node( x1 ) ) == 2 );
  CHECK( ntk.fanout_size( ntk.get_node( x2 ) ) == 2 );
  CHECK( ntk.fanout_size( ntk.get_node( f1 ) ) == 1 );
  CHECK( ntk.fanout_size( ntk.get_node( f2 ) ) == 1 );
}

TEST_CASE( "node and signal iteration in a mapped network", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );

  CHECK( has_foreach_node_v<mapped_network> );
  CHECK( has_foreach_pi_v<mapped_network> );
  CHECK( has_foreach_po_v<mapped_network> );
  CHECK( has_foreach_fanin_v<mapped_network> );

  const auto x1 = ntk.create_pi();
  const auto x2 = ntk.create_pi();
  const auto f1 = ntk.create_node( { x1, x2 }, 3 );
  const auto f2 = ntk.create_node( { x2, x1 }, 3 );
  ntk.create_po( f1 );
  ntk.create_po( f2 );

  CHECK( ntk.size() == 6 );

  /* iterate over nodes */
  uint32_t mask{ 0 }, counter{ 0 };
  ntk.foreach_node( [&]( auto n, auto i ) { mask |= ( 1 << n ); counter += i; } );
  CHECK( mask == 63 );
  CHECK( counter == 15 );

  mask = 0;
  ntk.foreach_node( [&]( auto n ) { mask |= ( 1 << n ); } );
  CHECK( mask == 63 );

  mask = counter = 0;
  ntk.foreach_node( [&]( auto n, auto i ) { mask |= ( 1 << n ); counter += i; return false; } );
  CHECK( mask == 1 );
  CHECK( counter == 0 );

  mask = 0;
  ntk.foreach_node( [&]( auto n ) { mask |= ( 1 << n ); return false; } );
  CHECK( mask == 1 );

  /* iterate over PIs */
  mask = counter = 0;
  ntk.foreach_pi( [&]( auto n, auto i ) { mask |= ( 1 << n ); counter += i; } );
  CHECK( mask == 12 );
  CHECK( counter == 1 );

  mask = 0;
  ntk.foreach_pi( [&]( auto n ) { mask |= ( 1 << n ); } );
  CHECK( mask == 12 );

  mask = counter = 0;
  ntk.foreach_pi( [&]( auto n, auto i ) { mask |= ( 1 << n ); counter += i; return false; } );
  CHECK( mask == 4 );
  CHECK( counter == 0 );

  mask = 0;
  ntk.foreach_pi( [&]( auto n ) { mask |= ( 1 << n ); return false; } );
  CHECK( mask == 4 );

  /* iterate over POs */
  mask = counter = 0;
  ntk.foreach_po( [&]( auto s, auto i ) { mask |= ( 1 << ntk.get_node( s ) ); counter += i; } );
  CHECK( mask == 48 );
  CHECK( counter == 1 );

  mask = 0;
  ntk.foreach_po( [&]( auto s ) { mask |= ( 1 << ntk.get_node( s ) ); } );
  CHECK( mask == 48 );

  mask = counter = 0;
  ntk.foreach_po( [&]( auto s, auto i ) { mask |= ( 1 << ntk.get_node( s ) ); counter += i; return false; } );
  CHECK( mask == 16 );
  CHECK( counter == 0 );

  mask = 0;
  ntk.foreach_po( [&]( auto s ) { mask |= ( 1 << ntk.get_node( s ) ); return false; } );
  CHECK( mask == 16 );
}


TEST_CASE( "visited values in mapped networks", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );

  CHECK( has_clear_visited_v<mapped_network> );
  CHECK( has_visited_v<mapped_network> );
  CHECK( has_set_visited_v<mapped_network> );

  const auto x1 = ntk.create_pi();
  const auto x2 = ntk.create_pi();
  const auto f1 = ntk.create_node( { x1, x2 }, 3 );
  const auto f2 = ntk.create_node( { x2, x1 }, 3 );
  ntk.create_po( f1 );
  ntk.create_po( f2 );

  CHECK( ntk.size() == 6 );

  ntk.clear_visited();
  ntk.foreach_node( [&]( auto n ) {
    CHECK( ntk.visited( n ) == 0 );
    ntk.set_visited( n, static_cast<uint32_t>( n ) );
    CHECK( ntk.visited( n ) == n );
  } );
  ntk.clear_visited();
  ntk.foreach_node( [&]( auto n ) {
    CHECK( ntk.visited( n ) == 0 );
  } );
}

TEST_CASE( "compute values in mapped", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );

  CHECK( has_compute_v<mapped_network, bool> );
  CHECK( has_compute_v<mapped_network, kitty::dynamic_truth_table> );
  CHECK( has_compute_v<mapped_network, kitty::partial_truth_table> );
  CHECK( has_compute_inplace_v<mapped_network, kitty::partial_truth_table> );

  const auto x1 = ntk.create_pi();
  const auto x2 = ntk.create_pi();
  const auto n1 = ntk.create_not( x1 );
  const auto n2 = ntk.create_not( x2 );
  const auto f1 = ntk.create_node( { n1, x2 }, 2 );
  const auto f2 = ntk.create_node( { x1, n2 }, 3 );
  ntk.create_po( f1 );
  ntk.create_po( f2 );

  {
    std::vector<bool> values{ { true, false } };

    CHECK( ntk.compute( ntk.get_node( f1 ), values.begin(), values.end() ) == true );
    CHECK( ntk.compute( ntk.get_node( f2 ), values.begin(), values.end() ) == false );
  }

  {
    std::vector<kitty::dynamic_truth_table> xs{ 2, kitty::dynamic_truth_table( 2 ) };
    kitty::create_nth_var( xs[0], 0 );
    kitty::create_nth_var( xs[1], 1 );

    CHECK( ntk.compute( ntk.get_node( f1 ), xs.begin(), xs.end() ) == ( ~xs[0] | ~xs[1] ) );
    CHECK( ntk.compute( ntk.get_node( f2 ), xs.begin(), xs.end() ) == ( xs[0] & xs[1] ) );
  }

  {
    std::vector<kitty::partial_truth_table> xs{ 2 };

    CHECK( ntk.compute( ntk.get_node( f1 ), xs.begin(), xs.end() ) == ( ~xs[0] | ~xs[1] ) );
    CHECK( ntk.compute( ntk.get_node( f2 ), xs.begin(), xs.end() ) == ( xs[0] & xs[1] ) );

    xs[0].add_bit( 0 );
    xs[1].add_bit( 1 );

    CHECK( ntk.compute( ntk.get_node( f1 ), xs.begin(), xs.end() ) == ( ~xs[0] | ~xs[1] ) );
    CHECK( ntk.compute( ntk.get_node( f2 ), xs.begin(), xs.end() ) == ( xs[0] & xs[1] ) );

    xs[0].add_bit( 1 );
    xs[1].add_bit( 0 );

    CHECK( ntk.compute( ntk.get_node( f1 ), xs.begin(), xs.end() ) == ( ~xs[0] | ~xs[1] ) );
    CHECK( ntk.compute( ntk.get_node( f2 ), xs.begin(), xs.end() ) == ( xs[0] & xs[1] ) );

    xs[0].add_bit( 0 );
    xs[1].add_bit( 0 );

    CHECK( ntk.compute( ntk.get_node( f1 ), xs.begin(), xs.end() ) == ( ~xs[0] | ~xs[1] ) );
    CHECK( ntk.compute( ntk.get_node( f2 ), xs.begin(), xs.end() ) == ( xs[0] & xs[1] ) );

    xs[0].add_bit( 1 );
    xs[1].add_bit( 1 );

    CHECK( ntk.compute( ntk.get_node( f1 ), xs.begin(), xs.end() ) == ( ~xs[0] | ~xs[1] ) );
    CHECK( ntk.compute( ntk.get_node( f2 ), xs.begin(), xs.end() ) == ( xs[0] & xs[1] ) );
  }

  {
    std::vector<kitty::partial_truth_table> xs{ 2 };
    kitty::partial_truth_table result;

    xs[0].add_bit( 0 );
    xs[1].add_bit( 1 );

    ntk.compute( ntk.get_node( f1 ), result, xs.begin(), xs.end() );
    CHECK( result == ( ~xs[0] | ~xs[1] ) );
    ntk.compute( ntk.get_node( f2 ), result, xs.begin(), xs.end() );
    CHECK( result == ( xs[0] & xs[1] ) );

    xs[0].add_bit( 1 );
    xs[1].add_bit( 0 );

    ntk.compute( ntk.get_node( f1 ), result, xs.begin(), xs.end() );
    CHECK( result == ( ~xs[0] | ~xs[1] ) );
    ntk.compute( ntk.get_node( f2 ), result, xs.begin(), xs.end() );
    CHECK( result == ( xs[0] & xs[1] ) );

    xs[0].add_bit( 0 );
    xs[1].add_bit( 0 );

    ntk.compute( ntk.get_node( f1 ), result, xs.begin(), xs.end() );
    CHECK( result == ( ~xs[0] | ~xs[1] ) );
    ntk.compute( ntk.get_node( f2 ), result, xs.begin(), xs.end() );
    CHECK( result == ( xs[0] & xs[1] ) );

    xs[0].add_bit( 1 );
    xs[1].add_bit( 1 );

    ntk.compute( ntk.get_node( f1 ), result, xs.begin(), xs.end() );
    CHECK( result == ( ~xs[0] | ~xs[1] ) );
    ntk.compute( ntk.get_node( f2 ), result, xs.begin(), xs.end() );
    CHECK( result == ( xs[0] & xs[1] ) );
  }
}

TEST_CASE( "Simulate XOR AIG mapped circuit with Booleans", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );

  const auto a = ntk.create_pi();
  const auto b = ntk.create_pi();
  const auto f1 = ntk.create_node( { a, b }, 2 );
  const auto f2 = ntk.create_node( { a, f1 }, 2 );
  const auto f3 = ntk.create_node( { b, f1 }, 2 );
  const auto f4 = ntk.create_node( { f2, f3 }, 2 );
  ntk.create_po( f4 );

  CHECK( !simulate<bool>( ntk, default_simulator<bool>( { false, false } ) )[0] );
  CHECK( simulate<bool>( ntk, default_simulator<bool>( { false, true } ) )[0] );
  CHECK( simulate<bool>( ntk, default_simulator<bool>( { true, false } ) )[0] );
  CHECK( !simulate<bool>( ntk, default_simulator<bool>( { false, false } ) )[0] );
}

TEST_CASE( "Simulate XOR AIG mapped circuit with static truth table", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );

  const auto a = ntk.create_pi();
  const auto b = ntk.create_pi();
  const auto f1 = ntk.create_node( { a, b }, 2 );
  const auto f2 = ntk.create_node( { a, f1 }, 2 );
  const auto f3 = ntk.create_node( { b, f1 }, 2 );
  const auto f4 = ntk.create_node( { f2, f3 }, 2 );
  ntk.create_po( f4 );

  const auto tt = simulate<kitty::static_truth_table<2u>>( ntk )[0];
  CHECK( tt._bits == 0x6 );
}

TEST_CASE( "Simulate XOR AIG mapped circuit with dynamic truth table", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );

  const auto a = ntk.create_pi();
  const auto b = ntk.create_pi();
  const auto f1 = ntk.create_node( { a, b }, 2 );
  const auto f2 = ntk.create_node( { a, f1 }, 2 );
  const auto f3 = ntk.create_node( { b, f1 }, 2 );
  const auto f4 = ntk.create_node( { f2, f3 }, 2 );
  ntk.create_po( f4 );

  default_simulator<kitty::dynamic_truth_table> sim( 2 );
  const auto tt = simulate<kitty::dynamic_truth_table>( ntk, sim )[0];
  CHECK( tt._bits[0] == 0x6 );
}

TEST_CASE( "Simulate XOR AIG mpped circuit with pre-defined values", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );

  const auto a = ntk.create_pi();
  const auto b = ntk.create_pi();
  const auto f1 = ntk.create_node( { a, b }, 2 );
  const auto f2 = ntk.create_node( { a, f1 }, 2 );
  const auto f3 = ntk.create_node( { b, f1 }, 2 );
  const auto f4 = ntk.create_node( { f2, f3 }, 2 );
  ntk.create_po( f4 );

  default_simulator<kitty::dynamic_truth_table> sim( 2 );

  unordered_node_map<kitty::dynamic_truth_table, mapped_network> node_to_value( ntk );
  simulate_nodes<kitty::dynamic_truth_table>( ntk, node_to_value, sim );

  CHECK( node_to_value[f4]._bits[0] == 0x6 );

  node_to_value.reset();

  /* set node f1 to false, such that function f1 becomes true */
  node_to_value[ntk.get_node( f1 )] = ~kitty::dynamic_truth_table( 2 );

  /* re-simulated with the fixed value for f1 */
  simulate_nodes<kitty::dynamic_truth_table>( ntk, node_to_value, sim );
  CHECK( node_to_value[f1]._bits[0] == 0xf );
  CHECK( node_to_value[f2]._bits[0] == 0x5 );
  CHECK( node_to_value[f3]._bits[0] == 0x3 );
  CHECK( node_to_value[f4]._bits[0] == 0xe );
}

TEST_CASE( "Partial simulator for mapped network", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );

  const auto a = ntk.create_pi();
  const auto b = ntk.create_pi();
  const auto f1 = ntk.create_node( { a, b }, 2 );
  const auto f2 = ntk.create_node( { a, f1 }, 2 );
  const auto f3 = ntk.create_node( { b, f1 }, 2 );
  const auto f4 = ntk.create_node( { f2, f3 }, 2 );
  ntk.create_po( f4 );

  std::vector<kitty::partial_truth_table> pats( 2 );
  pats[0].add_bits( 0x0a, 5 ); /* a = 01010 */
  pats[1].add_bits( 0x13, 5 ); /* b = 10011 */
  partial_simulator sim( pats );

  unordered_node_map<kitty::partial_truth_table, mapped_network> node_to_value( ntk );
  simulate_nodes( ntk, node_to_value, sim );

  CHECK( node_to_value[f1]._bits[0] == 0x1d ); /* f1 = 11101 */
  CHECK( node_to_value[f2]._bits[0] == 0x17 ); /* f2 = 10111 */
  CHECK( node_to_value[f3]._bits[0] == 0x0e ); /* f3 = 01110 */
  CHECK( node_to_value[f4]._bits[0] == 0x19 ); /* f4 = 11001 */
}

TEST_CASE( "Add pattern and re-simulate with partial_simulator in mapped network", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );

  const auto a = ntk.create_pi();
  const auto b = ntk.create_pi();
  const auto c = ntk.create_pi();
  const auto f1 = ntk.create_node( { a, b }, 4 );
  const auto f2 = ntk.create_node( { f1, c }, 4 );
  ntk.create_po( f2 ); /* f2 = a ^ b ^ c */

  std::vector<kitty::partial_truth_table> pats( 3 );
  pats[0].add_bits( 0x1, 2 ); /* a = 01 */
  pats[1].add_bits( 0x1, 2 ); /* b = 01 */
  pats[2].add_bits( 0x1, 2 ); /* c = 01 */
  partial_simulator sim( pats );

  unordered_node_map<kitty::partial_truth_table, mapped_network> node_to_value( ntk );
  simulate_nodes( ntk, node_to_value, sim );
  CHECK( sim.num_bits() == 2u );
  CHECK( node_to_value[f2]._bits[0] == 0x1 ); /* f2 = 01 */

  std::vector<bool> pattern( 3 );
  pattern[0] = 0;
  pattern[1] = 1;
  pattern[2] = 0;
  sim.add_pattern( pattern );
  pattern[0] = 1;
  pattern[1] = 0;
  pattern[2] = 0;
  sim.add_pattern( pattern );
  /* a = 1001, b = 0101, c = 0001 */

  simulate_nodes( ntk, node_to_value, sim, false );
  CHECK( sim.num_bits() == 4u );
  CHECK( node_to_value[f2]._bits[0] == 0xd ); /* f2 = 1101 */
}

TEST_CASE( "Incremental simulation with partial_simulator for mapped network", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );

  const auto a = ntk.create_pi();
  const auto b = ntk.create_pi();
  const auto c = ntk.create_pi();
  const auto f1 = ntk.create_node( { a, b }, 4 );
  const auto f2 = ntk.create_node( { f1, c }, 4 );
  ntk.create_po( f2 ); /* f2 = a ^ b ^ c */
  const auto f3 = ntk.create_node( { b, ntk.create_not( c )}, 4 );
  const auto f4 = ntk.create_node( { f3, a }, 4 );
  ntk.create_po( f4 ); /* f4 = !(a ^ b ^ c) = !f2 */

  std::vector<kitty::partial_truth_table> pats( 3 );
  pats[0].add_bits( 0, 64 );
  pats[1].add_bits( 0, 64 );
  pats[2].add_bits( 0, 64 );
  partial_simulator sim( pats );

  unordered_node_map<kitty::partial_truth_table, mapped_network> node_to_value( ntk );
  simulate_nodes( ntk, node_to_value, sim );
  CHECK( node_to_value[f2]._bits[0] == 0 );
  CHECK( node_to_value[f4]._bits[0] == ~uint64_t( 0 ) );

  std::vector<bool> pattern( 3 );
  pattern[0] = 0;
  pattern[1] = 1;
  pattern[2] = 0;
  sim.add_pattern( pattern );

  simulate_node( ntk, ntk.get_node( f2 ), node_to_value, sim );
  CHECK( node_to_value[f2].num_bits() == 65 );
  CHECK( node_to_value[f2]._bits[1] == 1 );
  CHECK( node_to_value[f1].num_bits() == 65 ); /* f1 is in f2's fanin cone hence it is re-simulated. */
  CHECK( node_to_value[f4].num_bits() == 64 ); /* f4 is not in f2's fanin cone hence it is not re-simulated. */

  const auto f5 = ntk.create_node( { f2, f4 }, 3 );
  simulate_node( ntk, ntk.get_node( f5 ), node_to_value, sim );
  CHECK( node_to_value[f5] == kitty::partial_truth_table( 65 ) );
}

TEST_CASE( "Construct from mapped for mapped network", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in_lib( test_library );
  auto result = lorina::read_genlib( in_lib, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<3, classification_type::p_configurations> lib( gates );

  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto c = aig.create_pi();

  const auto n4 = aig.create_and( a, b );
  const auto n5 = aig.create_and( b, c );
  const auto f = aig.create_and( n4, n5 );
  aig.create_po( f );

  emap_params ps;
  emap_stats st;
  binding_view<klut_network> luts = emap_klut( aig, lib, ps, &st );

  mapped_network ntk( gates, luts );

  CHECK( luts.compute_area() == ntk.compute_area() );

}

TEST_CASE( "create nodes and cleanup network", "[mapped]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mapped_network ntk( gates );

  CHECK( has_compute_v<mapped_network, kitty::dynamic_truth_table> );

  const auto a = ntk.create_pi();
  const auto b = ntk.create_pi();
  const auto c = ntk.create_pi();

  const auto n0 = ntk.create_node( { a, b, c }, 5 );
  const auto n1 = ntk.create_node( { a, b, c }, 6 );
  const auto n2 = ntk.create_node( { a, b }, 2 );
  const auto n3 = ntk.create_node( { c, n1 }, 2 );

  ntk.create_po( n0 );

  CHECK( ntk.size() == 9 );

  auto ntk1 = cleanup_mapped( ntk );

  CHECK( ntk1.size() == 6 );
}

//std::string const test_library = "GATE   inv1    1 O=!a;            PIN * INV 1 999 0.9 0.3 0.9 0.3\n"
//                                 "GATE   inv2    2 O=!a;            PIN * INV 2 999 1.0 0.1 1.0 0.1\n"
//                                 "GATE   nand2   2 O=!(a*b);        PIN * INV 1 999 1.0 0.2 1.0 0.2\n"
//                                 "GATE   and2    3 O=a*b;           PIN * INV 1 999 1.7 0.2 1.7 0.2\n"
//                                 "GATE   xor2    4 O=a^b;           PIN * UNKNOWN 2 999 1.9 0.5 1.9 0.5\n"
//                                 "GATE   mig3    3 O=a*b+a*c+b*c;   PIN * INV 1 999 2.0 0.2 2.0 0.2\n"
//                                 "GATE   xor3    5 O=a^b^c;         PIN * UNKNOWN 2 999 3.0 0.5 3.0 0.5\n"
//                                 "GATE   buf     2 O=a;             PIN * NONINV 1 999 1.0 0.0 1.0 0.0\n"
//                                 "GATE   zero    0 O=CONST0;\n"
//                                 "GATE   one     0 O=CONST1;\n"
//                                 "GATE   ha      5 C=a*b;           PIN * INV 1 999 1.7 0.4 1.7 0.4\n"
//                                 "GATE   ha      5 S=!a*b+a*!b;     PIN * INV 1 999 2.1 0.4 2.1 0.4\n"
//                                 "GATE   fa      6 C=a*b+a*c+b*c;   PIN * INV 1 999 2.1 0.4 2.1 0.4\n"
//                                 "GATE   fa      6 S=a^b^c;         PIN * INV 1 999 3.0 0.4 3.0 0.4";