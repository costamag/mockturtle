#include <catch.hpp>

#include <kitty/kitty.hpp>

#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/resynthesis/aig_syn.hpp>
#include <mockturtle/resynthesis/networks/bindings.hpp>
#include <mockturtle/resynthesis/placement/placement_analyzer.hpp>

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

TEST_CASE( "Basic placer functionalities", "[placement_analyzer]" )
{
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bindings lib( gates );
  using Ntk = mockturtle::bound_network;
  using signal = typename Ntk::signal;
  Ntk ntk( gates );
  std::vector<signal> xs;
  std::vector<signal> fs;
  for ( int i = 0; i < 3; ++i )
  {
    xs.push_back( ntk.create_pi() );
  }
  fs.push_back( ntk.create_node( { xs[0], xs[1] }, "and2" )[0] );
  fs.push_back( ntk.create_node( { xs[1], xs[2] }, "and2" )[0] );
  fs.push_back( ntk.create_node( { fs[0], fs[1] }, "and2" )[0] );
  fs.push_back( ntk.create_node( { fs[1], fs[2] }, "and2" )[0] );
  fs.push_back( ntk.create_node( { fs[2], fs[3] }, "and2" )[0] );
  fs.push_back( ntk.create_node( { fs[4], fs[3] }, "and2" )[0] );
  fs.push_back( ntk.create_node( { fs[0], fs[5] }, "and2" )[0] );
  ntk.create_po( fs.back() );

  mockturtle::placement_params ps;
  ps.area_density = 5.0/25.0;
  mockturtle::placement_analyzer<Ntk, 6, 6> analyzer( ntk, ps );
  analyzer.print_params_placer();
  std::cout << analyzer.compute_cost() << std::endl;
  analyzer.save_result( "/home/costamag/POST_MAPPING/mockturtle/build/graphs/test0_0" );
  analyzer.fire();
  analyzer.save_result( "/home/costamag/POST_MAPPING/mockturtle/build/graphs/test0_1" );
  std::cout << analyzer.compute_cost() << std::endl;
  //CHECK( analyzer.compute_cost_diffn() == 0 );
}
