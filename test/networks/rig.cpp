#include <catch.hpp>

#include <mockturtle/networks/rig.hpp>
#include <mockturtle/traits.hpp>


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