#include <catch.hpp>

#include <mockturtle/networks/rig.hpp>
#include <mockturtle/traits.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
//#include <mockturtle/io/blif_reader.hpp>
//#include <mockturtle/io/write_blif.hpp>

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
  //CHECK( rig.is_buf( rig.get_node(f1)));
  //CHECK( rig.is_not( rig.get_node(f2)));
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

  {
    // check strashing
    auto q = rig.create_and( x1, x2 );
    CHECK( rig.is_and(rig.get_node(q)) );

    CHECK( rig.size() == 4 );
    CHECK( q == f1 );

    // check permutation
    q = rig.create_and( x2, x1 );
    CHECK( rig.is_and(rig.get_node(q)) );
    CHECK( rig.size() == 4 );
    CHECK( q == f1 );
    // check and with constant 0
    q = rig.create_and( x2, rig.get_constant(0) );
    CHECK( !rig.is_and(rig.get_node(q)) );
    CHECK( rig.size() == 4 );
    CHECK( q == rig.get_constant(0) );
    q = rig.create_and( x1, rig.get_constant(0) );
    CHECK( !rig.is_and(rig.get_node(q)) );
    CHECK( rig.size() == 4 );
    CHECK( q == rig.get_constant(0) );

    // check and with constant 1
    q = rig.create_and( x2, rig.get_constant(1) );
    CHECK( !rig.is_and(rig.get_node(q)) );
    CHECK( rig.size() == 4 );
    CHECK( q == x2 );
    q = rig.create_and( x1, rig.get_constant(1) );
    CHECK( !rig.is_and(rig.get_node(q)) );
    CHECK( rig.size() == 4 );
    CHECK( q == x1 );

    // check and with same input
    q = rig.create_and( x2, x2 );
    CHECK( !rig.is_and(rig.get_node(q)) );
    CHECK( rig.size() == 4 );
    CHECK( q == x2 );
    q = rig.create_and( x1, !x1 );
    CHECK( !rig.is_and(rig.get_node(q)) );
    CHECK( rig.size() == 4 );
    CHECK( q == rig.get_constant(0) );
  }


  const auto f2 = rig.create_nand( x1, x2 );
  CHECK( rig.size() == 5 );
  CHECK( f1 != !f2 );

  {
    // check strashing
    auto q = rig.create_nand( x1, x2 );
    CHECK( rig.is_nand(rig.get_node(q)) );

    CHECK( rig.size() == 5 );
    CHECK( q == f2 );

    // check permutation
    q = rig.create_nand( x2, x1 );
    CHECK( rig.is_nand(rig.get_node(q)) );
    CHECK( rig.size() == 5 );
    CHECK( q == f2 );
    // check and with constant 0
    q = rig.create_nand( x2, rig.get_constant(0) );
    CHECK( !rig.is_nand(rig.get_node(q)) );
    CHECK( rig.size() == 5 );
    CHECK( q == rig.get_constant(1) );

    q = rig.create_nand( x1, rig.get_constant(0) );
    CHECK( !rig.is_nand(rig.get_node(q)) );
    CHECK( rig.size() == 5 );
    CHECK( q == rig.get_constant(1) );

    // check and with constant 1
    q = rig.create_nand( x2, rig.get_constant(1) );
    CHECK( !rig.is_nand(rig.get_node(q)) );
    CHECK( rig.size() == 5 );
    CHECK( q == !x2 );

    q = rig.create_nand( x1, rig.get_constant(1) );
    CHECK( !rig.is_nand(rig.get_node(q)) );
    CHECK( rig.size() == 5 );
    CHECK( q == !x1 );

    // check nand with same input
    q = rig.create_nand( x2, x2 );
    CHECK( !rig.is_nand(rig.get_node(q)) );
    CHECK( rig.size() == 5 );
    CHECK( q == !x2 );
    q = rig.create_nand( x1, !x1 );
    CHECK( !rig.is_nand(rig.get_node(q)) );
    CHECK( rig.size() == 5 );
    CHECK( q == rig.get_constant(1) );    
  }

  const auto f3 = rig.create_or( x1, x2 );
  CHECK( rig.is_or(rig.get_node(f3)) );
  CHECK( rig.size() == 6 );

  {
    // check strashing
    auto q = rig.create_or( x1, x2 );
    CHECK( rig.is_or(rig.get_node(q)) );

    CHECK( rig.size() == 6 );
    CHECK( q == f3 );

    // check permutation
    q = rig.create_or( x2, x1 );
    CHECK( rig.is_or(rig.get_node(q)) );
    CHECK( rig.size() == 6 );
    CHECK( q == f3 );
    // check or with constant 0
    q = rig.create_or( x2, rig.get_constant(0) );
    CHECK( !rig.is_or(rig.get_node(q)) );
    CHECK( rig.size() == 6 );
    CHECK( q == x2 );
    q = rig.create_or( rig.get_constant(0), x1 );
    CHECK( !rig.is_or(rig.get_node(q)) );
    CHECK( rig.size() == 6 );
    CHECK( q == x1 );

    // check and with constant 1
    q = rig.create_or( x2, rig.get_constant(1) );
    CHECK( !rig.is_or(rig.get_node(q)) );
    CHECK( rig.size() == 6 );
    CHECK( q == rig.get_constant(1) );
    q = rig.create_or( x1, rig.get_constant(1) );
    CHECK( !rig.is_or(rig.get_node(q)) );
    CHECK( rig.size() == 6 );
    CHECK( q == rig.get_constant(1) );

    // check nand with same input
    q = rig.create_or( x2, x2 );
    CHECK( !rig.is_or(rig.get_node(q)) );
    CHECK( rig.size() == 6 );
    CHECK( q == x2 );
    q = rig.create_or( x1, !x1 );
    CHECK( !rig.is_or(rig.get_node(q)) );
    CHECK( rig.size() == 6 );
    CHECK( q == rig.get_constant(1) );    
  }

  const auto f4 = rig.create_nor( x1, x2 );
  CHECK( rig.is_nor(rig.get_node(f4)) );

  CHECK( rig.size() == 7 );
  CHECK( f3 != !f4 );

  {
    // check strashing
    auto q = rig.create_nor( x1, x2 );
    CHECK( rig.is_nor(rig.get_node(q)) );

    CHECK( rig.size() == 7 );
    CHECK( q == f4 );

    // check permutation
    q = rig.create_nor( x2, x1 );
    CHECK( rig.is_nor(rig.get_node(q)) );
    CHECK( rig.size() == 7 );
    CHECK( q == f4 );

    // check nor with constant 0
    q = rig.create_nor( x2, rig.get_constant(0) );
    CHECK( !rig.is_nor(rig.get_node(q)) );
    CHECK( rig.size() == 7 );
    CHECK( q == !x2 );
    q = rig.create_nor( rig.get_constant(0), x1 );
    CHECK( !rig.is_nor(rig.get_node(q)) );
    CHECK( rig.size() == 7 );
    CHECK( q == !x1 );

    // check and with constant 1
    q = rig.create_nor( x2, rig.get_constant(1) );
    CHECK( !rig.is_nor(rig.get_node(q)) );
    CHECK( rig.size() == 7 );
    CHECK( q == rig.get_constant(0) );
    q = rig.create_nor( x1, rig.get_constant(1) );
    CHECK( !rig.is_nor(rig.get_node(q)) );
    CHECK( rig.size() == 7 );
    CHECK( q == rig.get_constant(0) );

    // check nand with same input
    q = rig.create_nor( x2, x2 );
    CHECK( !rig.is_nor(rig.get_node(q)) );
    CHECK( rig.size() == 7 );
    CHECK( q == !x2 );
    q = rig.create_nor( x1, !x1 );
    CHECK( !rig.is_nor(rig.get_node(q)) );
    CHECK( rig.size() == 7 );
    CHECK( q == rig.get_constant(0) );    
  }

  const auto f5 = rig.create_lt( x1, x2 );
  CHECK( rig.is_lt(rig.get_node(f5)) );
  CHECK( rig.size() == 8 );
  {
    // check strashing
    auto q = rig.create_lt( x1, x2 );
    CHECK( rig.is_lt(rig.get_node(q)) );

    CHECK( rig.size() == 8 );
    CHECK( q == f5 );

    // check permutation
    q = rig.create_lt( x2, x1 );
    CHECK( rig.is_lt(rig.get_node(q)) );
    CHECK( rig.size() == 9 );
    CHECK( q != f5 );

    // check nor with constant 0
    q = rig.create_lt( rig.get_constant(0), x2 );
    CHECK( !rig.is_lt(rig.get_node(q)) );
    CHECK( rig.size() == 9 );
    CHECK( q == x2 );
    q = rig.create_lt( x2, rig.get_constant(0) );
    CHECK( !rig.is_lt(rig.get_node(q)) );
    CHECK( rig.size() == 9 );
    CHECK( q == rig.get_constant(0) );

    // check and with constant 1
    q = rig.create_lt( rig.get_constant(1), x1 );
    CHECK( !rig.is_lt(rig.get_node(q)) );
    CHECK( rig.size() == 9 );
    CHECK( q == rig.get_constant(0) );
    q = rig.create_lt( x1, rig.get_constant(1) );
    CHECK( !rig.is_lt(rig.get_node(q)) );
    CHECK( rig.size() == 9 );
    CHECK( q == !x1 );

    // check nand with same input
    q = rig.create_lt( x2, x2 );
    CHECK( !rig.is_lt(rig.get_node(q)) );
    CHECK( rig.size() == 9 );
    CHECK( q == rig.get_constant(0) );
    q = rig.create_lt( x1, !x1 );
    CHECK( !rig.is_lt(rig.get_node(q)) );
    CHECK( rig.size() == 9 );
    CHECK( q == !x1 ); 
    q = rig.create_lt( !x1, x1 );
    CHECK( !rig.is_lt(rig.get_node(q)) );
    CHECK( rig.size() == 9 );
    CHECK( q == x1 );    
  }

  const auto f6 = rig.create_ge( x1, x2 );
  CHECK( rig.is_ge(rig.get_node(f6)) );
  CHECK( rig.size() == 10 );
  {
    // check strashing
    auto q = rig.create_ge( x1, x2 );
    CHECK( rig.is_ge(rig.get_node(q)) );

    CHECK( rig.size() == 10 );
    CHECK( q == f6 );

    // check permutation
    q = rig.create_ge( x2, x1 );
    CHECK( rig.is_ge(rig.get_node(q)) );
    CHECK( rig.size() == 11 );
    CHECK( q != f6 );

    // check nor with constant 0
    q = rig.create_ge( rig.get_constant(0), x2 );
    CHECK( !rig.is_ge(rig.get_node(q)) );
    CHECK( rig.size() == 11 );
    CHECK( q == !x2 );
    q = rig.create_ge( x2, rig.get_constant(0) );
    CHECK( !rig.is_ge(rig.get_node(q)) );
    CHECK( rig.size() == 11 );
    CHECK( q == rig.get_constant(1) );

    // check and with constant 1
    q = rig.create_ge( rig.get_constant(1), x1 );
    CHECK( !rig.is_ge(rig.get_node(q)) );
    CHECK( rig.size() == 11 );
    CHECK( q == rig.get_constant(1) );
    q = rig.create_ge( x1, rig.get_constant(1) );
    CHECK( !rig.is_ge(rig.get_node(q)) );
    CHECK( rig.size() == 11 );
    CHECK( q == x1 );

    // check nand with same input
    q = rig.create_ge( x2, x2 );
    CHECK( !rig.is_ge(rig.get_node(q)) );
    CHECK( rig.size() == 11 );
    CHECK( q == rig.get_constant(1) );
    q = rig.create_ge( x1, !x1 );
    CHECK( !rig.is_ge(rig.get_node(q)) );
    CHECK( rig.size() == 11 );
    CHECK( q == x1 ); 
    q = rig.create_ge( !x1, x1 );
    CHECK( !rig.is_ge(rig.get_node(q)) );
    CHECK( rig.size() == 11 );
    CHECK( q == !x1 );    
  }

  const auto f7 = rig.create_gt( x1, x2 );
  CHECK( rig.is_gt(rig.get_node(f7)) );
  CHECK( rig.size() == 12 );
  {
    // check strashing
    auto q = rig.create_gt( x1, x2 );
    CHECK( rig.is_gt(rig.get_node(q)) );

    CHECK( rig.size() == 12 );
    CHECK( q == f7 );

    // check permutation
    q = rig.create_gt( x2, x1 );
    CHECK( rig.is_gt(rig.get_node(q)) );
    CHECK( rig.size() == 13 );
    CHECK( q != f7 );

    // check nor with constant 0
    q = rig.create_gt( rig.get_constant(0), x2 );
    CHECK( !rig.is_gt(rig.get_node(q)) );
    CHECK( rig.size() == 13 );
    CHECK( q == rig.get_constant(0) );
    q = rig.create_gt( x2, rig.get_constant(0) );
    CHECK( !rig.is_gt(rig.get_node(q)) );
    CHECK( rig.size() == 13 );
    CHECK( q == x2 );

    // check and with constant 1
    q = rig.create_gt( rig.get_constant(1), x1 );
    CHECK( !rig.is_gt(rig.get_node(q)) );
    CHECK( rig.size() == 13 );
    CHECK( q == !x1 );
    q = rig.create_gt( x1, rig.get_constant(1) );
    CHECK( !rig.is_gt( rig.get_node(q) ) );
    CHECK( rig.size() == 13 );
    CHECK( q == rig.get_constant(false) );

    // check nand with same input
    q = rig.create_gt( x2, x2 );
    CHECK( !rig.is_gt(rig.get_node(q)) );
    CHECK( rig.size() == 13 );
    CHECK( q == rig.get_constant(0) );
    q = rig.create_gt( x1, !x1 );
    CHECK( !rig.is_gt( rig.get_node(q) ) );
    CHECK( rig.size() == 13 );
    CHECK( q == x1 ); 
    q = rig.create_gt( !x1, x1 );
    CHECK( !rig.is_gt(rig.get_node(q)) );
    CHECK( rig.size() == 13 );
    CHECK( q == !x1 );    
  }

  const auto f8 = rig.create_le( x1, x2 );
  CHECK( rig.is_le(rig.get_node(f8)) );
  CHECK( rig.size() == 14 );
  {
    // check strashing
    auto q = rig.create_le( x1, x2 );
    CHECK( rig.is_le(rig.get_node(q)) );

    CHECK( rig.size() == 14 );
    CHECK( q == f8 );

    // check permutation
    q = rig.create_le( x2, x1 );
    CHECK( rig.is_le(rig.get_node(q)) );
    CHECK( rig.size() == 15 );
    CHECK( q != f8 );

    // check nor with constant 0
    q = rig.create_le( rig.get_constant(0), x2 );
    CHECK( !rig.is_le(rig.get_node(q)) );
    CHECK( rig.size() == 15 );
    CHECK( q == rig.get_constant(1) );
    q = rig.create_le( x2, rig.get_constant(0) );
    CHECK( !rig.is_le(rig.get_node(q)) );
    CHECK( rig.size() == 15 );
    CHECK( q == !x2 );

    // check and with constant 1
    q = rig.create_le( rig.get_constant(1), x1 );
    CHECK( !rig.is_le( rig.get_node(q) ) );
    CHECK( rig.size() == 15 );
    CHECK( q == x1 );
    q = rig.create_le( x1, rig.get_constant(1) );
    CHECK( !rig.is_le( rig.get_node(q) ) );
    CHECK( rig.size() == 15 );
    CHECK( q == rig.get_constant(true) );

    // check nand with same input
    q = rig.create_le( x2, x2 );
    CHECK( !rig.is_le(rig.get_node(q)) );
    CHECK( rig.size() == 15 );
    CHECK( q == rig.get_constant(true) );
    q = rig.create_le( x1, !x1 );
    CHECK( !rig.is_le( rig.get_node(q) ) );
    CHECK( rig.size() == 15 );
    CHECK( q == !x1 ); 
    q = rig.create_le( !x1, x1 );
    CHECK( !rig.is_le(rig.get_node(q)) );
    CHECK( rig.size() == 15 );
    CHECK( q == x1 );    
  }

  const auto f9 = rig.create_xor( x1, x2 );
  CHECK( rig.is_xor(rig.get_node(f9)) );
  CHECK( rig.size() == 16 );
  {
    // check strashing
    auto q = rig.create_xor( x1, x2 );
    CHECK( rig.is_xor(rig.get_node(q)) );

    CHECK( rig.size() == 16 );
    CHECK( q == f9 );

    // check permutation
    q = rig.create_xor( x2, x1 );
    CHECK( rig.is_xor(rig.get_node(q)) );
    CHECK( q == f9 );

    // check nor with constant 0
    q = rig.create_xor( rig.get_constant(0), x2 );
    CHECK( !rig.is_xor( rig.get_node(q) ) );
    CHECK( q == x2 );

    // check and with constant 1
    q = rig.create_xor( rig.get_constant(1), x1 );
    CHECK( !rig.is_xor( rig.get_node(q) ) );
    CHECK( q == !x1 );

    // check nand with same input
    q = rig.create_xor( x2, x2 );
    CHECK( !rig.is_xor( rig.get_node(q) ) );
    CHECK( q == rig.get_constant( false ) );
    q = rig.create_xor( x1, !x1 );
    CHECK( !rig.is_xor( rig.get_node(q) ) );
    CHECK( q == rig.get_constant( true ) );    
  }

  const auto f10 = rig.create_xnor( x1, x2 );
  CHECK( rig.is_xnor(rig.get_node(f10)) );
  CHECK( rig.size() == 17 );
  {
    // check strashing
    auto q = rig.create_xnor( x1, x2 );
    CHECK( rig.is_xnor(rig.get_node(q)) );

    CHECK( rig.size() == 17 );
    CHECK( q == f10 );

    // check permutation
    q = rig.create_xnor( x2, x1 );
    CHECK( rig.is_xnor(rig.get_node(q)) );
    CHECK( q == f10 );

    // check nor with constant 0
    q = rig.create_xnor( rig.get_constant(0), x2 );
    CHECK( !rig.is_xnor( rig.get_node(q) ) );
    CHECK( q == !x2 );

    // check and with constant 1
    q = rig.create_xnor( rig.get_constant(1), x1 );
    CHECK( !rig.is_xnor( rig.get_node(q) ) );
    CHECK( q == x1 );

    // check nand with same input
    q = rig.create_xnor( x2, x2 );
    CHECK( !rig.is_xnor( rig.get_node(q) ) );
    CHECK( q == rig.get_constant( true ) );
    q = rig.create_xnor( x1, !x1 );
    CHECK( !rig.is_xnor( rig.get_node(q) ) );
    CHECK( q == rig.get_constant( false ) );    
  }

}

TEST_CASE( "create ternary operations in an RIG", "[rig]" )
{
  rig_network rig;

  CHECK( has_create_maj_v<rig_network> );
  CHECK( has_create_xor3_v<rig_network> );
  CHECK( has_create_ite_v<rig_network> );

  const auto x1 = rig.create_pi();
  const auto x2 = rig.create_pi();
  const auto x3 = rig.create_pi();

  CHECK( rig.size() == 4 );

  const auto f1 = rig.create_maj( x1, x2, x3 );
  const auto f2 = rig.create_maj( !x1, x2, !x3 );
  const auto f3 = rig.create_maj( x1, !x2, x3 );
  CHECK( rig.size() == 6 );
  CHECK( f2 == !f3 );
  {
    auto q = rig.create_maj( x3, x1, x2 );
    CHECK( q == f1 );
    q = rig.create_maj( x2, x3, x1 );
    CHECK( q == f1 );
    q = rig.create_maj( x1, x3, x2 );
    CHECK( q == f1 );
    q = rig.create_maj( x2, x1, x3 );
    CHECK( q == f1 );
    q = rig.create_maj( x3, x2, x1 );
    CHECK( q == f1 );
    q = rig.create_maj( x3, x2, x1 );
    CHECK( q == f1 );
    q = rig.create_maj( x1, x3, x2 );
    CHECK( q == f1 );
    q = rig.create_maj( x2, x1, x3 );
    CHECK( q == f1 );

    // two inputs are the same
    q = rig.create_maj( x1, x1, x2 );
    CHECK( q == x1 );
    q = rig.create_maj( x1, !x1, x2 );
    CHECK( q == x2 );
  }

  // TODO: test for xor3 and ite.

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

TEST_CASE( "simulate some 2-inputs functions in RIGs", "[rig]" )
{
  rig_network rig;
  const auto x1 = rig.create_pi();
  const auto x2 = rig.create_pi();

  const auto f1 = rig.create_and( x1, x2 );
  const auto f2 = rig.create_nand( x1, x2 );
  const auto f3 = rig.create_or( x1, x2 );

  rig.create_po( f1 );
  rig.create_po( f2 );
  rig.create_po( f3 );

  CHECK( rig.num_gates() == 3u );

  auto result = simulate<kitty::dynamic_truth_table>( rig, default_simulator<kitty::dynamic_truth_table>( 2 ) );

  CHECK( result[0]._bits[0] == 0x8 );
  CHECK( result[1]._bits[0] == 0x7 );
  CHECK( result[2]._bits[0] == 0xe );
}

TEST_CASE( "substitute input by constant in NAND-based XOR RIG", "[rig]" )
{
  rig_network rig;
  const auto x1 = rig.create_pi();
  const auto x2 = rig.create_pi();

  const auto f1 = rig.create_nand( x1, x2 );
  const auto f2 = rig.create_nand( x1, f1 );
  const auto f3 = rig.create_nand( x2, f1 );
  const auto f4 = rig.create_nand( f2, f3 );

  rig.create_po( f4 );

  CHECK( rig.num_gates() == 4u );
  auto sims = simulate<kitty::dynamic_truth_table>( rig, default_simulator<kitty::dynamic_truth_table>( 2 ) );

  CHECK( sims[0]._bits[0] == 0x6 );
  rig.substitute_node( rig.get_node( x1 ), rig.get_constant( true ) );
  CHECK( simulate<kitty::static_truth_table<2u>>( rig )[0]._bits == 0x3 );

  CHECK( rig.fanout_size( rig.get_node( f1 ) ) == 0u );
  CHECK( rig.fanout_size( rig.get_node( f2 ) ) == 0u );
  CHECK( rig.fanout_size( rig.get_node( f3 ) ) == 0u );
  CHECK( rig.fanout_size( rig.get_node( f4 ) ) == 0u );
}

TEST_CASE( "substitute node by constant in NAND-based XOR RIG", "[rig]" )
{
  rig_network rig;
  const auto x1 = rig.create_pi();
  const auto x2 = rig.create_pi();

  const auto f1 = rig.create_nand( x1, x2 );
  const auto f2 = rig.create_nand( x1, f1 );
  const auto f3 = rig.create_nand( x2, f1 );
  const auto f4 = rig.create_nand( f2, f3 );
  rig.create_po( f4 );

  CHECK( rig.num_gates() == 4u );
  CHECK( simulate<kitty::static_truth_table<2u>>( rig )[0]._bits == 0x6 );

  rig.substitute_node( rig.get_node( f3 ), rig.get_constant( true ) );

  CHECK( simulate<kitty::static_truth_table<2u>>( rig )[0]._bits == 0x2 );

  CHECK( rig.num_gates() == 2u );
  CHECK( rig.fanout_size( rig.get_node( f1 ) ) == 1u );
  CHECK( rig.fanout_size( rig.get_node( f2 ) ) == 1u );
  CHECK( rig.fanout_size( rig.get_node( f3 ) ) == 0u );
  CHECK( rig.fanout_size( rig.get_node( f4 ) ) == 0u );
  CHECK( !rig.is_dead( rig.get_node( f1 ) ) );
  CHECK( !rig.is_dead( rig.get_node( f2 ) ) );
  CHECK( rig.is_dead( rig.get_node( f3 ) ) );
  CHECK( rig.is_dead( rig.get_node( f4 ) ) );
}

TEST_CASE( "substitute node by constant in NAND-based XOR RIG (test case 2)", "[rig]" )
{
  rig_network rig;
  const auto x1 = rig.create_pi();
  const auto x2 = rig.create_pi();

  const auto f1 = rig.create_nand( x1, x2 );
  const auto f2 = rig.create_nand( x1, f1 );
  const auto f3 = rig.create_nand( x2, f1 );
  const auto f4 = rig.create_nand( f2, f3 );
  rig.create_po( f4 );

  CHECK( rig.num_gates() == 4u );
  CHECK( simulate<kitty::static_truth_table<2u>>( rig )[0]._bits == 0x6 );

  rig.substitute_node( rig.get_node( f1 ), rig.get_constant( true ) );

  CHECK( simulate<kitty::static_truth_table<2u>>( rig )[0]._bits == 0xe );

  CHECK( rig.fanout_size( rig.get_node( f1 ) ) == 0u );
  CHECK( rig.fanout_size( rig.get_node( f2 ) ) == 0u );
  CHECK( rig.fanout_size( rig.get_node( f3 ) ) == 0u );
  CHECK( rig.fanout_size( rig.get_node( f4 ) ) == 1u );
}

TEST_CASE( "invoke take_out_node two times on the same node RIG", "[rig]" )
{
  rig_network rig;
  const auto x1 = rig.create_pi();
  const auto x2 = rig.create_pi();

  const auto f1 = rig.create_and( x1, x2 );
  const auto f2 = rig.create_or( x1, x2 );
  (void)f2;

  CHECK( rig.fanout_size( rig.get_node( x1 ) ) == 2u );
  CHECK( rig.fanout_size( rig.get_node( x2 ) ) == 2u );

  /* delete node */
  CHECK( !rig.is_dead( rig.get_node( f1 ) ) );
  rig.take_out_node( rig.get_node( f1 ) );
  CHECK( rig.is_dead( rig.get_node( f1 ) ) );
  CHECK( rig.fanout_size( rig.get_node( x1 ) ) == 1u );
  CHECK( rig.fanout_size( rig.get_node( x2 ) ) == 1u );

  /* ensure that double-deletion has no effect on the fanout-size of x1 and x2 */
  CHECK( rig.is_dead( rig.get_node( f1 ) ) );
  rig.take_out_node( rig.get_node( f1 ) );
  CHECK( rig.is_dead( rig.get_node( f1 ) ) );
  CHECK( rig.fanout_size( rig.get_node( x1 ) ) == 1u );
  CHECK( rig.fanout_size( rig.get_node( x2 ) ) == 1u );
}

TEST_CASE( "substitute node and restrash RIG", "[rig]" )
{
  rig_network rig;
  auto const x1 = rig.create_pi();
  auto const x2 = rig.create_pi();

  auto const f1 = rig.create_and( x1, x2 );
  auto const f2 = rig.create_and( f1, x2 );
  rig.create_po( f2 );

  CHECK( rig.fanout_size( rig.get_node( x1 ) ) == 1 );
  CHECK( rig.fanout_size( rig.get_node( x2 ) ) == 2 );
  CHECK( rig.fanout_size( rig.get_node( f1 ) ) == 1 );
  CHECK( rig.fanout_size( rig.get_node( f2 ) ) == 1 );

  CHECK( simulate<kitty::static_truth_table<2u>>( rig )[0]._bits == 0x8 );

  /* substitute f1 with x1
   *
   * this is a very interesting test case because replacing f1 with x1
   * in f2 makes f2 and f1 equal.  a correct implementation will
   * create a new entry in the hash, although (x1, x2) is already
   * there, because (x1, x2) will be deleted in the next step.
   */
  rig.substitute_node( rig.get_node( f1 ), x1 );
  CHECK( simulate<kitty::static_truth_table<2u>>( rig )[0]._bits == 0x8 );

  CHECK( rig.fanout_size( rig.get_node( x1 ) ) == 1 );
  CHECK( rig.fanout_size( rig.get_node( x2 ) ) == 1 );
  CHECK( rig.fanout_size( rig.get_node( f1 ) ) == 0 );
  CHECK( rig.fanout_size( rig.get_node( f2 ) ) == 1 );
}

TEST_CASE( "substitute node with complemented node in rig_network", "[rig]" )
{
  rig_network rig;
  auto const x1 = rig.create_pi();
  auto const x2 = rig.create_pi();

  auto const f1 = rig.create_and( x1, x2 );
  auto const f2 = rig.create_and( x1, f1 );
  rig.create_po( f2 );

  CHECK( rig.fanout_size( rig.get_node( x1 ) ) == 2 );
  CHECK( rig.fanout_size( rig.get_node( x2 ) ) == 1 );
  CHECK( rig.fanout_size( rig.get_node( f1 ) ) == 1 );
  CHECK( rig.fanout_size( rig.get_node( f2 ) ) == 1 );

  CHECK( simulate<kitty::static_truth_table<2u>>( rig )[0]._bits == 0x8 );

  rig.substitute_node( rig.get_node( f2 ), !f2 );

  CHECK( rig.fanout_size( rig.get_node( x1 ) ) == 2 );
  CHECK( rig.fanout_size( rig.get_node( x2 ) ) == 1 );
  CHECK( rig.fanout_size( rig.get_node( f1 ) ) == 1 );
  CHECK( rig.fanout_size( rig.get_node( f2 ) ) == 1 );

  CHECK( simulate<kitty::static_truth_table<2u>>( rig )[0]._bits == 0x7 );
}

TEST_CASE( "substitute multiple RIG nodes", "[rig]" )
{
  using node = rig_network::node;
  using signal = rig_network::signal;

  rig_network rig;
  auto const x1 = rig.create_pi();
  auto const x2 = rig.create_pi();
  auto const x3 = rig.create_pi();

  auto const n4 = rig.create_and( !x1, x2 );
  auto const n5 = rig.create_and( x1, n4 );
  auto const n6 = rig.create_and( x3, n5 );
  auto const n7 = rig.create_and( n4, x2 );
  auto const n8 = rig.create_and( !n5, !n7 );
  auto const n9 = rig.create_and( !n8, n4 );

  rig.create_po( n6 );
  rig.create_po( n9 );

  rig.substitute_nodes( std::list<std::pair<node, signal>>{
      { rig.get_node( n5 ), rig.get_constant( false ) },
      { rig.get_node( n9 ), n4 } } );

  CHECK( !rig.is_dead( rig.get_node( rig.get_constant( false ) ) ) );
  CHECK( !rig.is_dead( rig.get_node( x1 ) ) );
  CHECK( !rig.is_dead( rig.get_node( x2 ) ) );
  CHECK( !rig.is_dead( rig.get_node( x3 ) ) );
  CHECK( !rig.is_dead( rig.get_node( n4 ) ) );
  CHECK( rig.is_dead( rig.get_node( n5 ) ) );
  CHECK( rig.is_dead( rig.get_node( n6 ) ) );
  CHECK( rig.is_dead( rig.get_node( n7 ) ) );
  CHECK( rig.is_dead( rig.get_node( n8 ) ) );
  CHECK( rig.is_dead( rig.get_node( n9 ) ) );

  CHECK( rig.fanout_size( rig.get_node( rig.get_constant( false ) ) ) == 1u );
  CHECK( rig.fanout_size( rig.get_node( x1 ) ) == 1u );
  CHECK( rig.fanout_size( rig.get_node( x2 ) ) == 1u );
  CHECK( rig.fanout_size( rig.get_node( x3 ) ) == 0u );
  CHECK( rig.fanout_size( rig.get_node( n4 ) ) == 1u );
  CHECK( rig.fanout_size( rig.get_node( n5 ) ) == 0u );
  CHECK( rig.fanout_size( rig.get_node( n6 ) ) == 0u );
  CHECK( rig.fanout_size( rig.get_node( n7 ) ) == 0u );
  CHECK( rig.fanout_size( rig.get_node( n8 ) ) == 0u );
  CHECK( rig.fanout_size( rig.get_node( n9 ) ) == 0u );

  rig.foreach_po( [&]( signal const o, uint32_t index ) {
    switch ( index )
    {
    case 0:
      CHECK( o == rig.get_constant( false ) );
      break;
    case 1:
      CHECK( o == n4 );
      break;
    default:
      CHECK( false );
    }
  } );
}

TEST_CASE( "substitute node with dependency in rig_network", "[rig]" )
{
  rig_network rig{};

  auto const a = rig.create_pi();
  auto const b = rig.create_pi();
  auto const c = rig.create_pi();          /* place holder */
  auto const tmp = rig.create_and( b, c ); /* place holder */
  auto const f1 = rig.create_and( a, b );
  auto const f2 = rig.create_and( f1, tmp );
  auto const f3 = rig.create_and( f1, a );
  rig.create_po( f2 );
  rig.substitute_node( rig.get_node( tmp ), f3 );

  /**
   * issue #545
   *
   *      f2
   *     /  \
   *    /   f3
   *    \  /  \
   *  1->f1    a
   *
   * stack:
   * 1. push (f2->f3)
   * 2. push (f3->a)
   * 3. pop (f3->a)
   * 4. pop (f2->f3) but, f3 is dead !!!
   */

  rig.substitute_node( rig.get_node( f1 ), rig.get_constant( 1 ) /* constant 1 */ );

  CHECK( rig.is_dead( rig.get_node( f1 ) ) );
  CHECK( rig.is_dead( rig.get_node( f2 ) ) );
  CHECK( rig.is_dead( rig.get_node( f3 ) ) );
  rig.foreach_po( [&]( auto s ) {
    CHECK( rig.is_dead( rig.get_node( s ) ) == false );
  } );
}

TEST_CASE( "substitute node and re-strash case 2 RIG", "[rig]" )
{
  rig_network rig;

  auto const x1 = rig.create_pi();
  auto const x2 = rig.create_pi();
  auto const x3 = rig.create_pi();
  auto const n4 = rig.create_and( x2, x3 );
  auto const n5 = rig.create_and( x1, n4 );
  auto const n6 = rig.create_and( n5, x3 );
  auto const n7 = rig.create_and( x1, n6 );
  rig.create_po( n7 );

  rig.substitute_node( rig.get_node( n6 ), n4 );
  /* replace in node n7: n6 <- n4 => re-strash with fanins (x1, n4) => n7 <- n5
   * take out node n6 => take out node n5 => take out node n4 (MFFC)
   * execute n7 <- n5, but n5 is dead => revive n5 and n4 */

  CHECK( !rig.is_dead( rig.get_node( n4 ) ) );
  CHECK( !rig.is_dead( rig.get_node( n5 ) ) );
  CHECK( rig.is_dead( rig.get_node( n6 ) ) );
  CHECK( rig.is_dead( rig.get_node( n7 ) ) );
  rig.foreach_fanin( rig.get_node( rig.po_at( 0 ) ), [&]( auto f, auto i ){
    switch ( i )
    {
    case 0:
      CHECK( f == x1 );
      break;
    case 1:
      CHECK( f == n4 );
      break;
    default:
      CHECK( false );
    }
  } );
  CHECK( rig.fanout_size( rig.get_node( n4 ) ) == 1 );
}

TEST_CASE( "substitute node without re-strashing case 1 RIG", "[rig]" )
{
  rig_network rig;
  auto const x1 = rig.create_pi();
  auto const x2 = rig.create_pi();
  auto const f1 = rig.create_and( x1, x2 );
  auto const f2 = rig.create_and( f1, x2 );
  rig.create_po( f2 );

  rig.substitute_node_no_restrash( rig.get_node( f1 ), x1 );
  rig = cleanup_dangling( rig );

  CHECK( rig.num_gates() == 1 );
  CHECK( simulate<kitty::static_truth_table<2u>>( rig )[0]._bits == 0x8 );
}

TEST_CASE( "substitute node with re-strashing case 2 RIG", "[rig]" )
{
  rig_network rig;

  auto const a = rig.create_pi();
  auto const b = rig.create_pi();
  auto const c = rig.create_pi();
  auto const tmp = rig.create_and( b, c );
  auto const f1 = rig.create_and( a, b );
  auto const f2 = rig.create_and( f1, tmp );
  auto const f3 = rig.create_and( f1, a );
  rig.create_po( f2 );
  rig.substitute_node( rig.get_node( tmp ), f3 );
  rig.substitute_node( rig.get_node( f1 ), rig.get_constant( 1 ) );
  rig = cleanup_dangling( rig );

  CHECK( rig.num_gates() == 0 );
  CHECK( !rig.is_dead( rig.get_node( rig.po_at( 0 ) ) ) );
  CHECK( rig.get_node( rig.po_at( 0 ) ) == rig.pi_at( 0 ) );
}

TEST_CASE( "substitute node without re-strashing case 2 RIG", "[rig]" )
{
  rig_network rig;

  auto const a = rig.create_pi();
  auto const b = rig.create_pi();
  auto const c = rig.create_pi();
  auto const tmp = rig.create_and( b, c );
  auto const f1 = rig.create_and( a, b );
  auto const f2 = rig.create_and( f1, tmp );
  auto const f3 = rig.create_and( f1, a );
  rig.create_po( f2 );
  rig.substitute_node_no_restrash( rig.get_node( tmp ), f3 );
  rig.substitute_node_no_restrash( rig.get_node( f1 ), rig.get_constant( 1 ) );
  rig = cleanup_rigs( rig );

  CHECK( rig.num_gates() == 0 );
  CHECK( !rig.is_dead( rig.get_node( rig.po_at( 0 ) ) ) );
  CHECK( rig.get_node( rig.po_at( 0 ) ) == rig.pi_at( 0 ) );
}

TEST_CASE( "substitute node without re-strashing case 3 RIG", "[rig]" )
{
  rig_network rig;

  auto const x1 = rig.create_pi();
  auto const x2 = rig.create_pi();
  auto const x3 = rig.create_pi();
  auto const n4 = rig.create_and( x2, x3 );
  auto const n5 = rig.create_and( x1, n4 );
  auto const n6 = rig.create_and( n5, x3 );
  auto const n7 = rig.create_and( x1, n6 );
  rig.create_po( n7 );

  rig.substitute_node_no_restrash( rig.get_node( n6 ), n4 );
  rig = cleanup_dangling( rig );
  CHECK( rig.num_gates() == 2 );
  CHECK( simulate<kitty::static_truth_table<3u>>( rig )[0]._bits == 0x80 );
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

TEST_CASE( "create a node in a RIG network (test permuting inputs)", "[rig]" )
{
  rig_network rig;

  CHECK( has_create_node_v<rig_network> );
  CHECK( has_compute_v<rig_network, kitty::dynamic_truth_table> );

  const auto x1 = rig.create_pi();
  const auto x2 = rig.create_pi();

  kitty::dynamic_truth_table tt1( 2u ), tt2( 2u ), tt_const0( 0u );
  kitty::create_from_hex_string( tt1, "2" );
  kitty::create_from_hex_string( tt2, "4" );

  CHECK( rig.size() == 3 );

  const auto _const0 = rig.create_node( {}, tt_const0 );
  const auto _const1 = rig.create_node( {}, ~tt_const0 );
  CHECK( _const0 == rig.get_constant( false ) );
  CHECK( _const1 == rig.get_constant( true ) );

  const auto f1 = rig.create_node( { x1, x2 }, tt1 );
  const auto f2 = rig.create_node( { x2, x1 }, tt2 );

  CHECK( rig.size() == 4 );

  std::vector<kitty::dynamic_truth_table> xs;
  xs.emplace_back( 2u );
  xs.emplace_back( 2u );
  kitty::create_nth_var( xs[0], 0 );
  kitty::create_nth_var( xs[1], 1 );

  const auto sim1 = rig.compute( rig.get_node( f1 ), xs.begin(), xs.end() );
  const auto sim2 = rig.compute( rig.get_node( f2 ), xs.begin(), xs.end() );

  CHECK( kitty::equal( sim2, sim1 ) );
  CHECK( kitty::equal( sim1, ~xs[1] & xs[0] ) );
}

TEST_CASE( "create a node in a RIG network (test cosnstant propagation)", "[rig]" )
{
  rig_network rig;

  CHECK( has_create_node_v<rig_network> );
  CHECK( has_compute_v<rig_network, kitty::dynamic_truth_table> );

  const auto x1 = rig.create_pi();
  const auto x2 = rig.create_pi();
  const auto x3 = rig.create_pi();

  kitty::dynamic_truth_table tt1( 3u ), tt2( 2u );
  kitty::create_from_hex_string( tt1, "a0" );
  kitty::create_from_hex_string( tt2, "8" );

  CHECK( rig.size() == 4 );

  const auto f1 = rig.create_node( { x1, x2, x3 }, tt1 );
  const auto f2 = rig.create_node( { x1, x3 }, tt2 );

  CHECK( rig.size() == 5 );

  std::vector<kitty::dynamic_truth_table> xs;
  xs.emplace_back( 3u );
  xs.emplace_back( 3u );
  kitty::create_nth_var( xs[0], 0 );
  kitty::create_nth_var( xs[1], 1 );

  const auto sim1 = rig.compute( rig.get_node( f1 ), xs.begin(), xs.end() );
  const auto sim2 = rig.compute( rig.get_node( f2 ), xs.begin(), xs.end() );

  CHECK( kitty::equal( sim2, sim1 ) );
  CHECK( kitty::equal( sim1, xs[1] & xs[0] ) );
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

//TEST_CASE( "read a combinational BLIF file into RIG network", "[rig]" )
//{
//  rig_network rig;
//
//  std::string file{
//      ".model top\n"
//      ".inputs a b c\n"
//      ".outputs y1 y2 y3 y4 y5\n"
//      ".names n0\n"
//      "0\n"
//      ".names n1\n"
//      "1\n"
//      ".names a b n1\n"
//      "11 1\n"
//      ".names c n1 n2\n"
//      "1- 1\n"
//      "-1 1\n"
//      ".names n2 y1\n"
//      "0 1\n"
//      ".names n0 y2\n"
//      "1 1\n"
//      ".names n0 y3\n"
//      "0 1\n"
//      ".names n1 y4\n"
//      "1 1\n"
//      ".names n1 y5\n"
//      "0 1\n"
//      ".end\n" };
//
//  std::istringstream in( file );
//  auto result = lorina::read_blif( in, blif_reader( rig ) );
//
//  /* structural checks */
//  CHECK( result == lorina::return_code::success );
//  CHECK( rig.size() == 6 );
//  CHECK( rig.num_pis() == 3 );
//  CHECK( rig.num_pos() == 5 );
//  CHECK( rig.num_gates() == 2 );
//
//  /* functional checks */
//  default_simulator<kitty::dynamic_truth_table> sim( rig.num_pis() );
//  const auto tts = simulate<kitty::dynamic_truth_table>( rig, sim );
//  rig.foreach_po( [&]( auto const&, auto i ) {
//    switch ( i )
//    {
//    case 0:
//      CHECK( kitty::to_hex( tts[i] ) == "07" );
//      break;
//    case 1:
//      CHECK( kitty::to_hex( tts[i] ) == "00" );
//      break;
//    }
//  } );
//}
//
//TEST_CASE( "read a combinational BLIF file into RIG network case 2", "[rig]" )
//{
//  rig_network rig;
//
//  std::string file{
//      ".model top\n"
//      ".inputs a b c\n"
//      ".outputs y1 y2\n"
//      ".names a b n1\n"
//      "10 1\n"
//      ".names a b n2\n"
//      "00 1\n"
//      "01 1\n"
//      "11 1\n"
//      ".names n1 y1\n"
//      "1 1\n"
//      ".names n2 y2\n"
//      "1 1\n"
//      ".end\n" };
//
//  std::istringstream in( file );
//  auto result = lorina::read_blif( in, blif_reader( rig ) );
//
//  /* structural checks */
//  CHECK( result == lorina::return_code::success );
//  CHECK( rig.size() == 6 );
//  CHECK( rig.num_pis() == 3 );
//  CHECK( rig.num_pos() == 2 );
//  CHECK( rig.num_gates() == 2 );
//
//  /* functional checks */
//  default_simulator<kitty::dynamic_truth_table> sim( rig.num_pis() );
//  const auto tts = simulate<kitty::dynamic_truth_table>( rig, sim );
//  rig.print();
//  rig.foreach_po( [&]( auto const&, auto i ) {
//    rig.print_aig( rig.po_at( i ) );
//    switch ( i )
//    {
//    case 0:
//      CHECK( kitty::to_hex( tts[i] ) == "22" );
//      break;
//    case 1:
//      CHECK( kitty::to_hex( tts[i] ) == "dd" );
//      break;
//    }
//  } );
//}