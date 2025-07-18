#include <catch.hpp>

#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>

#include <mockturtle/algorithms/mapped/windowing/window_manager.hpp>
#include <mockturtle/algorithms/mapped/windowing/window_simulator.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/mapped/bound_network.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>

std::string const test_library = "GATE   inv1    1.0 O=!a ;         PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   and2    1.0 O=a*b;         PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   or2     1.0 O=a+b;         PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   xor2    1.0 O=a^b;         PIN * INV 1   999 3.0 0.0 3.0 0.0";

TEST_CASE( "Simulate a small window", "[window_simulator]" )
{
  using Ntk = mockturtle::bound_network<mockturtle::bound::design_type_t::CELL_BASED, 2>;
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  Ntk ntk( gates );

  using signal = typename Ntk::signal;
  std::vector<signal> fs;
  /* Database construction from network */
  signal const a = ntk.create_pi();
  signal const b = ntk.create_pi();
  signal const c = ntk.create_pi();
  signal const d = ntk.create_pi();

  fs.push_back( ntk.create_node( { a, b }, 1 ) );         // 0
  fs.push_back( ntk.create_node( { c, d }, 1 ) );         // 1
  fs.push_back( ntk.create_node( { fs[0], fs[1] }, 1 ) ); // 2
  fs.push_back( ntk.create_node( { d }, 0 ) );            // 3
  fs.push_back( ntk.create_node( { a }, 0 ) );            // 4
  fs.push_back( ntk.create_node( { fs[4], fs[2] }, 2 ) ); // 5
  fs.push_back( ntk.create_node( { fs[3], fs[5] }, 2 ) ); // 6
  fs.push_back( ntk.create_node( { c }, 0 ) );            // 7
  fs.push_back( ntk.create_node( { fs[7], fs[6] }, 1 ) ); // 8
  fs.push_back( ntk.create_node( { fs[8] }, 0 ) );        // 9

  ntk.create_po( fs[9] );

  using DNtk = mockturtle::depth_view<Ntk>;
  mockturtle::window_manager_stats st;
  DNtk dntk( ntk );

  mockturtle::window_manager_params ps;
  ps.odc_levels = 4;
  ps.cut_limit = 8;
  mockturtle::window_manager<DNtk> window( dntk, ps, st );
  CHECK( window.run( dntk.get_node( fs[2] ) ) );
  mockturtle::window_simulator sim( dntk );
  sim.run( window );
  auto tt9 = sim.get( fs[9] );
  auto tta = sim.get( a );
  auto ttc = sim.get( c );
  auto ttd = sim.get( d );
  CHECK( kitty::equal( tt9, ttc | ( tta & ttd ) ) );
  auto care = sim.compute_observability_careset( window );
  CHECK( kitty::equal( care, ~ttc & ( tta & ttd ) ) );
}
