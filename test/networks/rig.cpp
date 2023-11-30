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