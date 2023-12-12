#include <catch.hpp>

#include <mockturtle/networks/rig.hpp>
#include <mockturtle/traits.hpp>
#include <mockturtle/algorithms/cleanup.hpp>


using namespace mockturtle;
using namespace rils;

TEST_CASE( "create and use constants in an RIG", "[rig]" )
{
  rig_network rig;

  CHECK( has_get_constant_v<rig_network> );
  CHECK( has_is_constant_v<rig_network> );
  CHECK( has_get_node_v<rig_network> );
  CHECK( has_is_complemented_v<rig_network> );
  const auto c0 = rig.get_constant( false );
  CHECK( rig.is_constant( rig.get_node( c0 ) ) );
  CHECK( !rig.is_pi( rig.get_node( c0 ) ) );
  CHECK( rig.size() == 1 );
  CHECK( std::is_same_v<std::decay_t<decltype( c0 )>, rig_network::signal> );
  CHECK( rig.get_node( c0 ) == 0 );
  CHECK( !rig.is_complemented( c0 ) );

  const auto c1 = rig.get_constant( true );

  CHECK( rig.get_node( c1 ) == 0 );
  CHECK( rig.is_complemented( c1 ) );
  CHECK( c0 != c1 );
  CHECK( c0 == !c1 );
  CHECK( ( !c0 ) == c1 );
  CHECK( ( !c0 ) != !c1 );
  CHECK( -c0 == c1 );
  CHECK( -c1 == c1 );
  CHECK( c0 == +c1 );
  CHECK( c0 == +c0 );
}


TEST_CASE( "create and use primary inputs in an RIG", "[rig]" )
{
  rig_network rig;

  CHECK( has_create_pi_v<rig_network> );

  auto a = rig.create_pi();
  auto b = rig.create_pi();

  CHECK( rig.size() == 3 ); // constant + two primary inputs
  CHECK( rig.num_pis() == 2 );
  CHECK( rig.num_gates() == 0 );
  CHECK( rig.is_pi( rig.get_node( a ) ) );
  CHECK( rig.is_pi( rig.get_node( b ) ) );
  CHECK( rig.pi_index( rig.get_node( a ) ) == 0 );
  CHECK( rig.pi_index( rig.get_node( b ) ) == 1 );
  CHECK( std::is_same_v<std::decay_t<decltype( a )>, rig_network::signal> );
  CHECK( a.index == 1 );
  CHECK( a.complement == 0 );

  a = !a;

  CHECK( a.index == 1 );
  CHECK( a.complement == 1 );

  a = +a;

  CHECK( a.index == 1 );
  CHECK( a.complement == 0 );

  a = +a;

  CHECK( a.index == 1 );
  CHECK( a.complement == 0 );

  a = -a;

  CHECK( a.index == 1 );
  CHECK( a.complement == 1 );

  a = -a;

  CHECK( a.index == 1 );
  CHECK( a.complement == 1 );

  a = a ^ true;

  CHECK( a.index == 1 );
  CHECK( a.complement == 0 );

  a = a ^ true;

  CHECK( a.index == 1 );
  CHECK( a.complement == 1 );
}

TEST_CASE( "create and use primary outputs in an RIG", "[rig]" )
{
  rig_network rig;

  CHECK( has_create_po_v<rig_network> );

  const auto c0 = rig.get_constant( false );
  const auto x1 = rig.create_pi();

  CHECK( rig.size() == 2 );
  CHECK( rig.num_pis() == 1 );
  CHECK( rig.num_pos() == 0 );

  rig.create_po( c0 );
  rig.create_po( x1 );
  rig.create_po( !x1 );

  CHECK( rig.size() == 2 );
  CHECK( rig.num_pos() == 3 );

  rig.foreach_po( [&]( auto s, auto i ) {
    switch ( i )
    {
    case 0:
      CHECK( s == c0 );
      break;
    case 1:
      CHECK( s == x1 );
      break;
    case 2:
      CHECK( s == !x1 );
      break;
    }
  } );
}

TEST_CASE( "create unary operations in a RIG network", "[rig]" )
{
  rig_network rig;

  CHECK( has_create_buf_v<rig_network> );
  CHECK( has_create_not_v<rig_network> );

  auto x1 = rig.create_pi();

  CHECK( rig.size() == 2 );

  auto f1 = rig.create_buf( x1 );
  auto f2 = rig.create_not( x1 );
  auto f3 = rig.create_buf( x1 );

  CHECK( rig.size() == 2 );
  CHECK( rig.is_pi( rig.get_node(f1)));
  CHECK( rig.is_pi( rig.get_node(x1)));
  CHECK( rig.is_buf( rig.get_node(f1)));
  CHECK( rig.is_not( rig.get_node(f2)));
  CHECK( f1 == f3 );
}

TEST_CASE( "create binary operations in an RIG", "[rig]" )
{
  rig_network rig;

  CHECK( has_create_and_v<rig_network> );
  CHECK( has_create_nand_v<rig_network> );
  CHECK( has_create_or_v<rig_network> );
  CHECK( has_create_nor_v<rig_network> );
  CHECK( has_create_xor_v<rig_network> );
  CHECK( has_create_xnor_v<rig_network> );

  const auto x1 = rig.create_pi();
  const auto x2 = rig.create_pi();

  CHECK( rig.size() == 3 );

  const auto f1 = rig.create_and( x1, x2 );
  CHECK( rig.size() == 4 );

  const auto s1 = rig.create_and( x1, x2 );
  CHECK( rig.is_and(rig.get_node(s1)) );

  CHECK( rig.size() == 4 );
  CHECK( s1 == f1 );

  const auto f2 = rig.create_nand( x1, x2 );
  CHECK( rig.size() == 5 );
  CHECK( f1 != !f2 );

  const auto s2 = rig.create_nand( x2, x1 );
  CHECK( rig.is_nand(rig.get_node(s2)) );

  CHECK( rig.size() == 5 );
  CHECK( s2 == f2 );

  const auto f3 = rig.create_or( x1, x2 );
  const auto s3 = rig.create_or( x2, x1 );
  CHECK( rig.is_or(rig.get_node(f3)) );

  CHECK( rig.size() == 6 );

  const auto f4 = rig.create_nor( x1, x2 );
  CHECK( rig.is_nor(rig.get_node(f4)) );

  CHECK( rig.size() == 7 );
  CHECK( f3 != !f4 );
  CHECK( s3 == f3 );

  const auto f5 = rig.create_xor( x1, x2 );
  const auto s5 = rig.create_xor( x2, x1 );
  CHECK( rig.is_xor(rig.get_node(s5)) );
  CHECK( rig.size() == 8 );

  const auto f6 = rig.create_xnor( x1, x2 );
  const auto s6 = rig.create_xnor( x2, x1 );
  CHECK( rig.is_xnor(rig.get_node(s6)) );
  CHECK( rig.size() == 9 );
  CHECK( f5 != !f6 );
  CHECK( f5 == s5 );
  CHECK( f6 == s6 );

  const auto f7 = !rig.create_xnor( x2, x1 );
  CHECK( f7 == !s6 );
  CHECK( rig.size() == 9 );

}

TEST_CASE( "hash nodes in RIG network", "[rig]" )
{
  rig_network rig;

  auto a = rig.create_pi();
  auto b = rig.create_pi();

  auto f = rig.create_and( a, b );
  auto g = rig.create_and( a, b );

  CHECK( rig.size() == 4u );
  CHECK( rig.num_gates() == 1u );

  CHECK( rig.get_node( f ) == rig.get_node( g ) );
}

TEST_CASE( "clone a RIG network", "[rig]" )
{
  CHECK( has_clone_v<rig_network> );

  rig_network rig0;
  auto a = rig0.create_pi();
  auto b = rig0.create_pi();
  auto f0 = rig0.create_and( a, b );
  CHECK( rig0.size() == 4 );
  CHECK( rig0.num_gates() == 1 );

  auto rig1 = rig0;
  auto rig_clone = rig0.clone();

  auto c = rig1.create_pi();
  rig1.create_and( f0, c );
  CHECK( rig0.size() == 6 );
  CHECK( rig0.num_gates() == 2 );

  CHECK( rig_clone.size() == 4 );
  CHECK( rig_clone.num_gates() == 1 );
}

TEST_CASE( "clone a node in RIG network", "[rig]" )
{
  rig_network rig1, rig2;

  CHECK( has_clone_node_v<rig_network> );

  auto a1 = rig1.create_pi();
  auto b1 = rig1.create_pi();
  auto f1 = rig1.create_and( a1, b1 );
  CHECK( rig1.size() == 4 );

  auto a2 = rig2.create_pi();
  auto b2 = rig2.create_pi();
  CHECK( rig2.size() == 3 );

  auto f2 = rig2.clone_node( rig1, rig1.get_node( f1 ), { a2, b2 } );
  CHECK( rig2.size() == 4 );

  rig2.foreach_fanin( rig2.get_node( f2 ), [&]( auto const& s, auto ) {
    CHECK( !rig2.is_complemented( s ) );
  } );
}

TEST_CASE( "structural properties of an RIG", "[rig]" )
{
  rig_network rig;

  CHECK( has_size_v<rig_network> );
  CHECK( has_num_pis_v<rig_network> );
  CHECK( has_num_pos_v<rig_network> );
  CHECK( has_num_gates_v<rig_network> );
  CHECK( has_fanin_size_v<rig_network> );
  CHECK( has_fanout_size_v<rig_network> );

  const auto x1 = rig.create_pi();
  const auto x2 = rig.create_pi();

  const auto f1 = rig.create_and( x1, x2 );
  const auto f2 = rig.create_or( x1, x2 );

  rig.create_po( f1 );
  rig.create_po( f2 );

  CHECK( rig.size() == 5 );
  CHECK( rig.num_pis() == 2 );
  CHECK( rig.num_pos() == 2 );
  CHECK( rig.num_gates() == 2 );
  CHECK( rig.fanin_size( rig.get_node( x1 ) ) == 1 );
  CHECK( rig.fanin_size( rig.get_node( x2 ) ) == 1 );
  CHECK( rig.fanin_size( rig.get_node( f1 ) ) == 2 );
  CHECK( rig.fanin_size( rig.get_node( f2 ) ) == 2 );
  CHECK( rig.fanout_size( rig.get_node( x1 ) ) == 2 );
  CHECK( rig.fanout_size( rig.get_node( x2 ) ) == 2 );
  CHECK( rig.fanout_size( rig.get_node( f1 ) ) == 1 );
  CHECK( rig.fanout_size( rig.get_node( f2 ) ) == 1 );
}

TEST_CASE( "hash generic nodes in RIG network", "[rig]" )
{
  rig_network rig;

  const auto a = rig.create_pi();
  const auto b = rig.create_pi();
  const auto c = rig.create_pi();

  std::vector<kitty::dynamic_truth_table> sims;
  for( int i{0}; i<3; ++i )
  {
    sims.emplace_back(3u);
    kitty::create_nth_var( sims[i], i );
  }

  kitty::dynamic_truth_table tt_maj( 3u ), tt_xor( 3u );
  kitty::create_from_hex_string( tt_maj, "e8" );
  kitty::create_from_hex_string( tt_xor, "96" );

  auto s1 = rig.create_node( { a, b, c }, tt_maj );
  auto s2 = rig.create_node( { a, b, c }, tt_xor );

  CHECK( rig.size() == 6 );

  rig.create_node( { a, b, c }, tt_maj );

  CHECK( rig.size() == 6 );

  auto sim_1 = rig.compute( rig.get_node(s1), sims );
  auto sim_2 = rig.compute( rig.get_node(s2), sims );
  
  CHECK( kitty::equal( sim_1, tt_maj ) );
  CHECK( kitty::equal( sim_2, tt_xor ) );
}

TEST_CASE( "node and signal iteration in an RIG", "[rig]" )
{
  rig_network rig;

  CHECK( has_foreach_node_v<rig_network> );
  CHECK( has_foreach_pi_v<rig_network> );
  CHECK( has_foreach_po_v<rig_network> );
  CHECK( has_foreach_gate_v<rig_network> );
  CHECK( has_foreach_fanin_v<rig_network> );

  const auto x1 = rig.create_pi();
  const auto x2 = rig.create_pi();
  const auto f1 = rig.create_and( x1, x2 );
  const auto f2 = rig.create_or( x1, x2 );
  rig.create_po( f1 );
  rig.create_po( f2 );

  CHECK( rig.size() == 5 );

  /* iterate over nodes */
  uint32_t mask{ 0 }, counter{ 0 };
  rig.foreach_node( [&]( auto n, auto i ) { mask |= ( 1 << n ); counter += i; } );
  CHECK( mask == 31 );
  CHECK( counter == 10 );

  mask = 0;
  rig.foreach_node( [&]( auto n ) { mask |= ( 1 << n ); } );
  CHECK( mask == 31 );

  mask = counter = 0;
  rig.foreach_node( [&]( auto n, auto i ) { mask |= ( 1 << n ); counter += i; return false; } );
  CHECK( mask == 1 );
  CHECK( counter == 0 );

  mask = 0;
  rig.foreach_node( [&]( auto n ) { mask |= ( 1 << n ); return false; } );
  CHECK( mask == 1 );

  /* iterate over PIs */
  mask = counter = 0;
  rig.foreach_pi( [&]( auto n, auto i ) { mask |= ( 1 << n ); counter += i; } );
  CHECK( mask == 6 );
  CHECK( counter == 1 );

  mask = 0;
  rig.foreach_pi( [&]( auto n ) { mask |= ( 1 << n ); } );
  CHECK( mask == 6 );

  mask = counter = 0;
  rig.foreach_pi( [&]( auto n, auto i ) { mask |= ( 1 << n ); counter += i; return false; } );
  CHECK( mask == 2 );
  CHECK( counter == 0 );

  mask = 0;
  rig.foreach_pi( [&]( auto n ) { mask |= ( 1 << n ); return false; } );
  CHECK( mask == 2 );

  /* iterate over POs */
  mask = counter = 0;
  rig.foreach_po( [&]( auto s, auto i ) { mask |= ( 1 << rig.get_node( s ) ); counter += i; } );
  CHECK( mask == 24 );
  CHECK( counter == 1 );

  mask = 0;
  rig.foreach_po( [&]( auto s ) { mask |= ( 1 << rig.get_node( s ) ); } );
  CHECK( mask == 24 );

  mask = counter = 0;
  rig.foreach_po( [&]( auto s, auto i ) { mask |= ( 1 << rig.get_node( s ) ); counter += i; return false; } );
  CHECK( mask == 8 );
  CHECK( counter == 0 );

  mask = 0;
  rig.foreach_po( [&]( auto s ) { mask |= ( 1 << rig.get_node( s ) ); return false; } );
  CHECK( mask == 8 );

  /* iterate over gates */
  mask = counter = 0;
  rig.foreach_gate( [&]( auto n, auto i ) { mask |= ( 1 << n ); counter += i; } );
  CHECK( mask == 24 );
  CHECK( counter == 1 );

  mask = 0;
  rig.foreach_gate( [&]( auto n ) { mask |= ( 1 << n ); } );
  CHECK( mask == 24 );

  mask = counter = 0;
  rig.foreach_gate( [&]( auto n, auto i ) { mask |= ( 1 << n ); counter += i; return false; } );
  CHECK( mask == 8 );
  CHECK( counter == 0 );

  mask = 0;
  rig.foreach_gate( [&]( auto n ) { mask |= ( 1 << n ); return false; } );
  CHECK( mask == 8 );

  /* iterate over fanins */
  mask = counter = 0;
  rig.foreach_fanin( rig.get_node( f1 ), [&]( auto s, auto i ) { mask |= ( 1 << rig.get_node( s ) ); counter += i; } );
  CHECK( mask == 6 );
  CHECK( counter == 1 );

  mask = 0;
  rig.foreach_fanin( rig.get_node( f1 ), [&]( auto s ) { mask |= ( 1 << rig.get_node( s ) ); } );
  CHECK( mask == 6 );

  mask = counter = 0;
  rig.foreach_fanin( rig.get_node( f1 ), [&]( auto s, auto i ) { mask |= ( 1 << rig.get_node( s ) ); counter += i; return false; } );
  CHECK( mask == 2 );
  CHECK( counter == 0 );

  mask = 0;
  rig.foreach_fanin( rig.get_node( f1 ), [&]( auto s ) { mask |= ( 1 << rig.get_node( s ) ); return false; } );
  CHECK( mask == 2 );
}

TEST_CASE( "compute values in RIGs", "[rig]" )
{
  rig_network rig;

  CHECK( has_compute_v<rig_network, bool> );
  CHECK( has_compute_v<rig_network, kitty::dynamic_truth_table> );
  CHECK( has_compute_v<rig_network, kitty::partial_truth_table> );
  //CHECK( has_compute_inplace_v<rig_network, kitty::partial_truth_table> );

  const auto x1 = rig.create_pi();
  const auto x2 = rig.create_pi();
  const auto f1 = rig.create_and( !x1, x2 );
  const auto f2 = rig.create_and( x1, !x2 );
  rig.create_po( f1 );
  rig.create_po( f2 );

  {
    std::vector<bool> values{ { true, false } };

    CHECK( rig.compute( rig.get_node( f1 ), values.begin(), values.end() ) == false );
    CHECK( rig.compute( rig.get_node( f2 ), values.begin(), values.end() ) == true );
  }

  {
    std::vector<kitty::dynamic_truth_table> xs{ 2, kitty::dynamic_truth_table( 2 ) };
    kitty::create_nth_var( xs[0], 0 );
    kitty::create_nth_var( xs[1], 1 );

    CHECK( rig.compute( rig.get_node( f1 ), xs.begin(), xs.end() ) == ( ~xs[0] & xs[1] ) );
    CHECK( rig.compute( rig.get_node( f2 ), xs.begin(), xs.end() ) == ( xs[0] & ~xs[1] ) );
  }

  {
    std::vector<kitty::partial_truth_table> xs{ 2 };

    CHECK( rig.compute( rig.get_node( f1 ), xs.begin(), xs.end() ) == ( ~xs[0] & xs[1] ) );
    CHECK( rig.compute( rig.get_node( f2 ), xs.begin(), xs.end() ) == ( xs[0] & ~xs[1] ) );

    xs[0].add_bit( 0 );
    xs[1].add_bit( 1 );

    CHECK( rig.compute( rig.get_node( f1 ), xs.begin(), xs.end() ) == ( ~xs[0] & xs[1] ) );
    CHECK( rig.compute( rig.get_node( f2 ), xs.begin(), xs.end() ) == ( xs[0] & ~xs[1] ) );

    xs[0].add_bit( 1 );
    xs[1].add_bit( 0 );

    CHECK( rig.compute( rig.get_node( f1 ), xs.begin(), xs.end() ) == ( ~xs[0] & xs[1] ) );
    CHECK( rig.compute( rig.get_node( f2 ), xs.begin(), xs.end() ) == ( xs[0] & ~xs[1] ) );

    xs[0].add_bit( 0 );
    xs[1].add_bit( 0 );

    CHECK( rig.compute( rig.get_node( f1 ), xs.begin(), xs.end() ) == ( ~xs[0] & xs[1] ) );
    CHECK( rig.compute( rig.get_node( f2 ), xs.begin(), xs.end() ) == ( xs[0] & ~xs[1] ) );

    xs[0].add_bit( 1 );
    xs[1].add_bit( 1 );

    CHECK( rig.compute( rig.get_node( f1 ), xs.begin(), xs.end() ) == ( ~xs[0] & xs[1] ) );
    CHECK( rig.compute( rig.get_node( f2 ), xs.begin(), xs.end() ) == ( xs[0] & ~xs[1] ) );
  }

  {
    std::vector<kitty::partial_truth_table> xs{ 2 };
    kitty::partial_truth_table result;

    xs[0].add_bit( 0 );
    xs[1].add_bit( 1 );

    rig.compute( rig.get_node( f1 ), result, xs.begin(), xs.end() );
    CHECK( result == ( ~xs[0] & xs[1] ) );
    rig.compute( rig.get_node( f2 ), result, xs.begin(), xs.end() );
    CHECK( result == ( xs[0] & ~xs[1] ) );

    xs[0].add_bit( 1 );
    xs[1].add_bit( 0 );

    rig.compute( rig.get_node( f1 ), result, xs.begin(), xs.end() );
    CHECK( result == ( ~xs[0] & xs[1] ) );
    rig.compute( rig.get_node( f2 ), result, xs.begin(), xs.end() );
    CHECK( result == ( xs[0] & ~xs[1] ) );

    xs[0].add_bit( 0 );
    xs[1].add_bit( 0 );

    rig.compute( rig.get_node( f1 ), result, xs.begin(), xs.end() );
    CHECK( result == ( ~xs[0] & xs[1] ) );
    rig.compute( rig.get_node( f2 ), result, xs.begin(), xs.end() );
    CHECK( result == ( xs[0] & ~xs[1] ) );

    xs[0].add_bit( 1 );
    xs[1].add_bit( 1 );

    rig.compute( rig.get_node( f1 ), result, xs.begin(), xs.end() );
    CHECK( result == ( ~xs[0] & xs[1] ) );
    rig.compute( rig.get_node( f2 ), result, xs.begin(), xs.end() );
    CHECK( result == ( xs[0] & ~xs[1] ) );
  }
}

TEST_CASE( "custom node values in RIGs", "[rig]" )
{
  rig_network rig;

  CHECK( has_clear_values_v<rig_network> );
  CHECK( has_value_v<rig_network> );
  CHECK( has_set_value_v<rig_network> );
  CHECK( has_incr_value_v<rig_network> );
  CHECK( has_decr_value_v<rig_network> );

  const auto x1 = rig.create_pi();
  const auto x2 = rig.create_pi();
  const auto f1 = rig.create_and( x1, x2 );
  const auto f2 = rig.create_or( x1, x2 );
  rig.create_po( f1 );
  rig.create_po( f2 );

  CHECK( rig.size() == 5 );

  rig.clear_values();
  rig.foreach_node( [&]( auto n ) {
    CHECK( rig.value( n ) == 0 );
    rig.set_value( n, static_cast<uint32_t>( n ) );
    CHECK( rig.value( n ) == n );
    CHECK( rig.incr_value( n ) == n );
    CHECK( rig.value( n ) == n + 1 );
    CHECK( rig.decr_value( n ) == n );
    CHECK( rig.value( n ) == n );
  } );
  rig.clear_values();
  rig.foreach_node( [&]( auto n ) {
    CHECK( rig.value( n ) == 0 );
  } );
}

TEST_CASE( "visited values in RIGs", "[rig]" )
{
  rig_network rig;

  CHECK( has_clear_visited_v<rig_network> );
  CHECK( has_visited_v<rig_network> );
  CHECK( has_set_visited_v<rig_network> );

  const auto x1 = rig.create_pi();
  const auto x2 = rig.create_pi();
  const auto f1 = rig.create_and( x1, x2 );
  const auto f2 = rig.create_or( x1, x2 );
  rig.create_po( f1 );
  rig.create_po( f2 );

  CHECK( rig.size() == 5 );

  rig.clear_visited();
  rig.foreach_node( [&]( auto n ) {
    CHECK( rig.visited( n ) == 0 );
    rig.set_visited( n, static_cast<uint32_t>( n ) );
    CHECK( rig.visited( n ) == static_cast<uint32_t>( n ) );
  } );
  rig.clear_visited();
  rig.foreach_node( [&]( auto n ) {
    CHECK( rig.visited( n ) == 0 );
  } );
}

TEST_CASE( "simulate some special functions in RIGs", "[rig]" )
{
  rig_network rig;
  const auto x1 = rig.create_pi();
  const auto x2 = rig.create_pi();
  const auto x3 = rig.create_pi();

  const auto f1 = rig.create_maj( x1, x2, x3 );
  const auto f2 = rig.create_ite( x1, x2, x3 );

  rig.create_po( f1 );
  rig.create_po( f2 );

  CHECK( rig.num_gates() == 2u );

  auto result = simulate<kitty::dynamic_truth_table>( rig, default_simulator<kitty::dynamic_truth_table>( 3 ) );

  CHECK( result[0]._bits[0] == 0xe8u );
  CHECK( result[1]._bits[0] == 0xd8u );
}

TEST_CASE( "substitute nodes with propagation in rigs (test case 1)", "[rig]" )
{
  CHECK( has_substitute_node_v<rig_network> );
  CHECK( has_replace_in_node_v<rig_network> );

  rig_network rig;
  const auto x1 = rig.create_pi();
  const auto x2 = rig.create_pi();
  const auto x3 = rig.create_pi();
  const auto x4 = rig.create_pi();

  const auto f1 = rig.create_and( x1, x2 );
  const auto f2 = rig.create_and( x3, x4 );
  const auto f3 = rig.create_and( x1, x3 );
  const auto f4 = rig.create_and( f1, f2 );
  const auto f5 = rig.create_and( f3, f4 );

  rig.create_po( f5 );

  CHECK( rig.size() == 10u );
  CHECK( rig.num_gates() == 5u );
  CHECK( rig._e_storage->hash.size() == 5u );
  CHECK( rig._e_storage->nodes[f1.index].children[0u].index == x1.index );
  CHECK( rig._e_storage->nodes[f1.index].children[1u].index == x2.index );

  CHECK( rig._e_storage->nodes[f5.index].children[0u].index == f3.index );
  CHECK( rig._e_storage->nodes[f5.index].children[1u].index == f4.index );

  CHECK( rig.fanout_size( rig.get_node( f1 ) ) == 1u );
  CHECK( rig.fanout_size( rig.get_node( f3 ) ) == 1u );
  CHECK( rig.fanout_size( rig.get_node( x1 ) ) == 2u );
  CHECK( rig.fanout_size( rig.get_node( x3 ) ) == 2u );
  CHECK( !rig.is_dead( rig.get_node( f1 ) ) );

  rig.substitute_node( rig.get_node( x2 ), x3 );

  // Node of signal f1 is now relabelled
  CHECK( rig.size() == 10u );
  CHECK( rig.num_gates() == 4u );
  CHECK( rig._e_storage->hash.size() == 4u );
  CHECK( rig._e_storage->nodes[f1.index].children[0u].index == x1.index );
  CHECK( rig._e_storage->nodes[f1.index].children[1u].index == x2.index );

  CHECK( rig._e_storage->nodes[f5.index].children[0u].index == f3.index );
  CHECK( rig._e_storage->nodes[f5.index].children[1u].index == f4.index );

  CHECK( rig.fanout_size( rig.get_node( f1 ) ) == 0u );
  CHECK( rig.fanout_size( rig.get_node( f3 ) ) == 2u );
  CHECK( rig.is_dead( rig.get_node( f1 ) ) );

  rig = cleanup_dangling( rig );

  CHECK( rig.num_gates() == 4u );
}

TEST_CASE( "substitute nodes with propagation in rigs (test case 2)", "[rig]" )
{
  rig_network rig;
  const auto x1 = rig.create_pi();
  const auto x2 = rig.create_pi();
  const auto x3 = rig.create_pi();

  const auto f1 = rig.create_and( x1, x2 );
  const auto f2 = rig.create_and( x1, x3 );
  const auto f3 = rig.create_and( f1, f2 );

  rig.create_po( f3 );

  CHECK( rig.num_gates() == 3u );
  CHECK( rig._e_storage->hash.size() == 3u );
  CHECK( rig._e_storage->nodes[f1.index].children[0u].index == x1.index );
  CHECK( rig._e_storage->nodes[f1.index].children[1u].index == x2.index );
  CHECK( rig._e_storage->nodes[f2.index].children[0u].index == x1.index );
  CHECK( rig._e_storage->nodes[f2.index].children[1u].index == x3.index );
  CHECK( rig._e_storage->nodes[f3.index].children[0u].index == f1.index );
  CHECK( rig._e_storage->nodes[f3.index].children[1u].index == f2.index );
  CHECK( rig._e_storage->outputs[0].index == f3.index );

  CHECK( rig.fanout_size( rig.get_node( f1 ) ) == 1u );
  CHECK( rig.fanout_size( rig.get_node( f2 ) ) == 1u );
  CHECK( rig.fanout_size( rig.get_node( f3 ) ) == 1u );

  rig.substitute_node( rig.get_node( x2 ), x3 );

  // Node of signal f1 is now relabelled
  CHECK( rig.num_gates() == 1u );
  CHECK( rig._e_storage->hash.size() == 1u );
  CHECK( rig._e_storage->nodes[f1.index].children[0u].index == x1.index );
  CHECK( rig._e_storage->nodes[f1.index].children[1u].index == x2.index );
  CHECK( rig._e_storage->nodes[f2.index].children[0u].index == x1.index );
  CHECK( rig._e_storage->nodes[f2.index].children[1u].index == x3.index );
  CHECK( rig._e_storage->nodes[f3.index].children[0u].index == f1.index );
  CHECK( rig._e_storage->nodes[f3.index].children[1u].index == f2.index );
  CHECK( rig._e_storage->outputs[0].index == f2.index );

  CHECK( rig.fanout_size( rig.get_node( f1 ) ) == 0u );
  CHECK( rig.fanout_size( rig.get_node( f2 ) ) == 1u );
  CHECK( rig.fanout_size( rig.get_node( f3 ) ) == 0u );

  rig = cleanup_dangling( rig );

  CHECK( rig.num_gates() == 1u );
}