#include <catch.hpp>

#include <mockturtle/generators/arithmetic.hpp>
#include <mockturtle/networks/mapped/bound_network.hpp>
#include <mockturtle/utils/signal_map.hpp>

#include <cstdint>
#include <vector>

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

TEST_CASE( "create incomplete signal map for full adder", "[signal_map]" )
{
  using Ntk = mockturtle::bound_network<2>;
  using signal = Ntk::signal;

  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  Ntk ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  /* chain of inverters */
  auto const f = ntk.create_node( { a, b }, { 12, 13 } );
  auto carry = signal{ f.index, 0 };
  auto sum = signal{ f.index, 1 };
  ntk.create_po( sum );
  ntk.create_po( carry );

  /* create incomplete node map */
  incomplete_signal_map<uint32_t, Ntk> map{ ntk };
  ntk.foreach_node( [&]( auto n ) {
    ntk.foreach_output_pin( n, [&]( auto pin, auto i ) {
      (void)pin;
      CHECK( !map.has( signal{ n, i } ) );
    } );
  } );

  int index = 1;
  ntk.foreach_node( [&]( auto n ) {
    ntk.foreach_output_pin( n, [&]( auto pin, auto j ) {
      (void)pin;
      map[signal{ n, j }] = index++;
    } );
  } );

  ntk.foreach_node( [&]( auto n ) {
    ntk.foreach_output_pin( n, [&]( auto pin, auto i ) {
      (void)pin;
      CHECK( map.has( signal{ n, i } ) );
    } );
  } );

  int total{ 0 };
  ntk.foreach_node( [&]( auto n ) {
    ntk.foreach_output_pin( n, [&]( auto pin, auto j ) {
      (void)pin;
      total += map[signal{ n, j }];
    } );
  } );

  CHECK( total == ( index * ( index - 1 ) ) / 2 );

  /* reset all values to 1 */
  map.reset();
  ntk.foreach_node( [&]( auto n ) {
    ntk.foreach_output_pin( n, [&]( auto pin, auto i ) {
      (void)pin;
      CHECK( !map.has( signal{ n, i } ) );
    } );
  } );

  ntk.foreach_node( [&]( auto n ) {
    ntk.foreach_output_pin( n, [&]( auto pin, auto j ) {
      (void)pin;
      map[signal{ n, j }] = 1;
    } );
  } );

  total = 0;
  ntk.foreach_node( [&]( auto n ) {
    ntk.foreach_output_pin( n, [&]( auto pin, auto j ) {
      (void)pin;
      total += map[signal{ n, j }];
    } );
  } );

  CHECK( total == 4 );

  /* test erase */
  map.erase( a );
  CHECK( !map.has( a ) );

  /* test resize */
  const auto d = ntk.create_pi();
  map.resize();
  CHECK( !map.has( d ) );

  map[d] = map[a] = 1;
  total = 0;
  ntk.foreach_node( [&]( auto n ) {
    ntk.foreach_output_pin( n, [&]( auto pin, auto j ) {
      (void)pin;
      total += map[signal{ n, j }];
    } );
  } );
  CHECK( total == 5 );

  /* reset with initial value 10 */
  map.reset( 10 );
  /* create with initial value 10 */
  incomplete_signal_map<uint32_t, Ntk> map2( ntk, 10 );

  ntk.foreach_node( [&]( auto n ) {
    ntk.foreach_output_pin( n, [&]( auto pin, auto j ) {
      (void)pin;
      signal s{ n, j };
      CHECK( map[s] == map2[s] );
    } );
  } );
}