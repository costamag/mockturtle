#include <catch.hpp>

#include <cstdint>
#include <vector>

#include <lorina/genlib.hpp>
#include <mockturtle/algorithms/mapped/evaluators/area_resyn_evaluator.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/super_reader.hpp>
#include <mockturtle/networks/mapped/bound_network.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>

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

TEST_CASE( "Area evaluator for resynthesis of mapped networks", "[area_resyn_evaluator]" )
{
  using bound_network = mockturtle::bound_network<bound::design_type_t::CELL_BASED, 2>;
  using signal = bound_network::signal;
  using node = bound_network::node;

  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  bound_network ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const f1 = ntk.create_node( { a, b }, 2 );
  auto const f2 = ntk.create_node( { b, c }, 2 );
  auto const f3 = ntk.create_node( { c, d }, 2 );
  auto const f4 = ntk.create_node( { f1, f2 }, 2 );
  auto const f5 = ntk.create_node( { f3, f4 }, 2 );
  auto const f6 = ntk.create_node( { f4, f5 }, 2 );
  auto const f7 = ntk.create_node( { f2, f3 }, 2 );
  ntk.create_po( f6 );
  ntk.create_po( f7 );

  evaluator_params ps;
  ps.max_num_roots = 7;
  area_resyn_evaluator evaluator( ntk, ps );

  std::vector<node> sorted_nodes;
  evaluator.foreach_gate( [&]( auto n ) {
    sorted_nodes.push_back( n );
  } );

  CHECK( ntk.area() == 14 );
  CHECK( sorted_nodes[0] == 11 );
  CHECK( sorted_nodes[1] == 9 );
  CHECK( evaluator.evaluate( f6.index, std::vector<signal>( { a, b, c, d } ) ) == 8 );
  CHECK( evaluator.evaluate( f6.index, std::vector<signal>( { f1, f2, f3 } ) ) == 6 );
  CHECK( evaluator.evaluate( f6.index, std::vector<signal>( { f4, f5 } ) ) == 2 );

  bound_list<bound::design_type_t::CELL_BASED> list;
  list.add_inputs( 3 );
  auto const la = list.pi_at( 0 );
  auto const lb = list.pi_at( 1 );
  auto const lc = list.pi_at( 2 );
  auto const l1 = list.add_gate( { la, lb }, 2 );
  auto const l2 = list.add_gate( { lb, lc }, 2 );
  auto const l3 = list.add_gate( { l1, l2 }, 4 );
  list.add_output( l3 );
  auto cost = evaluator.evaluate( list, std::vector<signal>( { a, b, c } ) );
  CHECK( cost == 4 );
  CHECK( ntk.area() == 14 );
}