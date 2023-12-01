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

  CHECK( rig.size() == 4 );
  CHECK( !rig.is_pi( rig.get_node(f1)));
  CHECK( rig.is_pi( rig.get_node(x1)));
  CHECK( rig.is_buf( rig.get_node(f1)));
  CHECK( rig.is_not( rig.get_node(f2)));
  CHECK( f1 == f3 );
}