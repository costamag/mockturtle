#include <catch.hpp>

#include <mockturtle/networks/rig.hpp>

using namespace mockturtle;

TEST_CASE( "create and use constants in an RIG", "[rig]" )
{
  rig_network rig;

//  CHECK( rig.size() == 1 );
//  CHECK( has_get_constant_v<aig_network> );
//  CHECK( has_is_constant_v<aig_network> );
//  CHECK( has_get_node_v<aig_network> );
//  CHECK( has_is_complemented_v<aig_network> );
//
//  const auto c0 = aig.get_constant( false );
//  CHECK( aig.is_constant( aig.get_node( c0 ) ) );
//  CHECK( !aig.is_pi( aig.get_node( c0 ) ) );
//
//  CHECK( aig.size() == 1 );
//  CHECK( std::is_same_v<std::decay_t<decltype( c0 )>, aig_network::signal> );
//  CHECK( aig.get_node( c0 ) == 0 );
//  CHECK( !aig.is_complemented( c0 ) );
//
//  const auto c1 = aig.get_constant( true );
//
//  CHECK( aig.get_node( c1 ) == 0 );
//  CHECK( aig.is_complemented( c1 ) );
//
//  CHECK( c0 != c1 );
//  CHECK( c0 == !c1 );
//  CHECK( ( !c0 ) == c1 );
//  CHECK( ( !c0 ) != !c1 );
//  CHECK( -c0 == c1 );
//  CHECK( -c1 == c1 );
//  CHECK( c0 == +c1 );
//  CHECK( c0 == +c0 );
}