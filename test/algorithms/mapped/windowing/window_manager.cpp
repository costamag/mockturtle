#include <catch.hpp>

#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>

#include <mockturtle/algorithms/mapped/windowing/window_manager.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/mapped/bound_network.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>

std::string const test_library = "GATE   inv1    1.0 O=!a ;         PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   and2    1.0 O=a*b;         PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   xor2    1.0 O=a^b;         PIN * INV 1   999 3.0 0.0 3.0 0.0";

TEST_CASE( "Window construction simple inverter chain", "[window_manager]" )
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
  fs.push_back( ntk.create_node( { a }, 0 ) );
  fs.push_back( ntk.create_node( { fs[0] }, 0 ) );
  fs.push_back( ntk.create_node( { fs[1] }, 0 ) );

  ntk.create_po( fs.back() );

  mockturtle::window_manager_params wps;
  using DNtk = mockturtle::depth_view<Ntk>;
  mockturtle::window_manager_stats<DNtk> wst;
  DNtk dntk( ntk );
  mockturtle::window_manager<DNtk, 5> win_manager( dntk, wps, wst );
  win_manager.run( dntk.get_node( fs.back() ) );
  CHECK( win_manager.is_valid() );
  auto const& leaves = win_manager.get_leaves();
  auto const& divs = win_manager.get_divs();
  auto const& mffc = win_manager.get_mffc();
  CHECK( leaves.size() == 1 );
  CHECK( leaves[0] == ntk.get_node( a ) );
  CHECK( mffc.size() == 3 );
  CHECK( mffc[0] == ntk.get_node( fs[0] ) );
  CHECK( mffc[1] == ntk.get_node( fs[1] ) );
  CHECK( mffc[2] == ntk.get_node( fs[2] ) );
  CHECK( divs.size() == 1 );
  CHECK( divs[0] == leaves[0] );
}