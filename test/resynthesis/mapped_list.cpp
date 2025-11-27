#include <catch.hpp>

#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>

#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/resynthesis/aig_syn.hpp>
#include <mockturtle/resynthesis/database/homo_database.hpp>
#include <mockturtle/resynthesis/database/mapped_list.hpp>
#include <mockturtle/resynthesis/networks/bindings.hpp>
#include <mockturtle/resynthesis/networks/bound.hpp>
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

TEST_CASE( "Construct mapped list from network", "[mapped_list]" )
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
  auto const x4 = bound.create_pi();
  auto const x5 = bound.create_pi();

  using signal = typename mockturtle::bound_network::signal;
  using node = typename mockturtle::bound_network::node;
  signal const f0 = bound.create_node( std::vector<signal>( { x0, x1 } ), "and2" )[0];
  signal const f1 = bound.create_node( std::vector<signal>( { x2, x3 } ), "and2" )[0];
  signal const f2 = bound.create_node( std::vector<signal>( { f0, f1 } ), "and2" )[0];
  bound.create_po( f2 );

  std::vector<kitty::dynamic_truth_table> xs;
  for( int i = 0; i < 6; ++i )
  {
    xs.emplace_back( 6 );
    kitty::create_nth_var( xs[i], i );
  }
  auto tt = bound.compute( bound.get_node( f0 ), xs.begin(), xs.begin() + 2 );

  CHECK( kitty::equal( tt, xs[0] & xs[1] ) );

  std::array<node, 6u> cut;
  cut[0] = bound.get_node( x0 );
  cut[1] = bound.get_node( x1 );
  cut[2] = bound.get_node( x2 );
  cut[3] = bound.get_node( x3 );
  cut[4] = bound.get_node( x4 );
  cut[5] = bound.get_node( x5 );

  mockturtle::mapped_list mlist( bound, gates, cut, f2 );

  auto const bs = bound.get_bindings();
  mockturtle::homo_database<6> db( bs );
  db.add_entry( bound, cut, f2 );

  signal const n0 = bound.create_node( std::vector<signal>( { x2, x5 } ), "and2" )[0];
  signal const n1 = bound.create_node( std::vector<signal>( { n0, x0 } ), "and2" )[0];
  signal const n2 = bound.create_node( std::vector<signal>( { x3, x4 } ), "xor2" )[0];
  signal const n3 = bound.create_node( std::vector<signal>( { n2, x1 } ), "nand2" )[0];
  signal const n4 = bound.create_node( std::vector<signal>( { n3, n1 } ), "and2" )[0];
  /* x0 -> x2  x5 -> x1 */

  db.add_entry( bound, cut, n4 );
  


}

