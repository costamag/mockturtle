#include <catch.hpp>

#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/traits.hpp>
#include <mockturtle/views/simulation_view.hpp>

#include <kitty/static_truth_table.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/partial_truth_table.hpp>
#include <kitty/constructors.hpp>

using namespace mockturtle;

template<typename Ntk, typename TT>
void test_simulation_view()
{
  CHECK( is_network_type_v<Ntk> );
  CHECK( !has_simulation_v<Ntk> );
  CHECK( !has_level_v<Ntk> );
  CHECK( !has_simulation_v<Ntk> );
  CHECK( !has_faninsize_v<Ntk> );

  using simulation_ntk = simulation_view<Ntk,TT>;

  CHECK( is_network_type_v<simulation_ntk> );
  CHECK( has_simulation_v<simulation_ntk> );
  CHECK( has_level_v<simulation_ntk> );
  CHECK( has_simulation_v<simulation_ntk> );
  CHECK( has_faninsize_v<simulation_ntk> );

  using simulationsimulation_ntk = simulation_view<simulation_ntk, TT>;

  CHECK( is_network_type_v<simulationsimulation_ntk> );
  CHECK( has_simulation_v<simulationsimulation_ntk> );
  CHECK( has_level_v<simulationsimulation_ntk> );
  CHECK( has_simulation_v<simulationsimulation_ntk> );
  CHECK( has_faninsize_v<simulationsimulation_ntk> );
};

TEST_CASE( "create different simulation views", "[simulation_view]" )
{

  using TT0 = kitty::static_truth_table<2u>;

  test_simulation_view<aig_network,TT0>();
  test_simulation_view<mig_network,TT0>();
  test_simulation_view<klut_network,TT0>();

  using TT1 = kitty::dynamic_truth_table;

  test_simulation_view<aig_network,TT1>();
  test_simulation_view<mig_network,TT1>();
  test_simulation_view<klut_network,TT1>();

  using TT2 = kitty::partial_truth_table;

  test_simulation_view<aig_network,TT2>();
  test_simulation_view<mig_network,TT2>();
  test_simulation_view<klut_network,TT2>();
}

TEST_CASE( "compute depth, levels, fanin size, and simulations for AIG", "[simulation_view]" )
{
  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto f1 = aig.create_nand( a, b );
  const auto f2 = aig.create_nand( a, f1 );
  const auto f3 = aig.create_nand( b, f1 );
  const auto f4 = aig.create_nand( f2, f3 );
  aig.create_po( f4 );

  using TT = kitty::dynamic_truth_table;
  using Ntk = aig_network;

  std::vector<kitty::dynamic_truth_table> xs{ 2, kitty::dynamic_truth_table( 2u ) };
  kitty::create_nth_var( xs[0], 0 );
  kitty::create_nth_var( xs[1], 1 );
  kitty::dynamic_truth_table const0 = xs[0].construct();

  simulation_view<Ntk,TT> simulation_aig{ aig };
  simulation_aig.set_input_simulations( xs );
  simulation_aig.update_simulations();
  simulation_aig.update_faninsizes();

  CHECK( simulation_aig.simulation( simulation_aig.get_node(simulation_aig.get_constant( false )) ) == const0 );

  CHECK( simulation_aig.depth() == 3 );
  CHECK( simulation_aig.level( aig.get_node( a ) ) == 0 );
  CHECK( simulation_aig.level( aig.get_node( b ) ) == 0 );
  CHECK( simulation_aig.level( aig.get_node( f1 ) ) == 1 );
  CHECK( simulation_aig.level( aig.get_node( f2 ) ) == 2 );
  CHECK( simulation_aig.level( aig.get_node( f3 ) ) == 2 );
  CHECK( simulation_aig.level( aig.get_node( f4 ) ) == 3 );

  CHECK( simulation_aig.faninsize( aig.get_node( a ) ) == 0 );
  CHECK( simulation_aig.faninsize( aig.get_node( b ) ) == 0 );
  CHECK( simulation_aig.faninsize( aig.get_node( f1 ) ) == 2 );
  CHECK( simulation_aig.faninsize( aig.get_node( f2 ) ) == 3 );
  CHECK( simulation_aig.faninsize( aig.get_node( f3 ) ) == 3 );
  CHECK( simulation_aig.faninsize( aig.get_node( f4 ) ) == 5 );

  kitty::dynamic_truth_table tt(2u);
  kitty::create_from_binary_string( tt, "1000" );
  CHECK( simulation_aig.simulation( simulation_aig.get_node(f1) ) == tt );
  kitty::create_from_binary_string( tt, "0010" );
  CHECK( simulation_aig.simulation( simulation_aig.get_node(f2) ) == tt );
  kitty::create_from_binary_string( tt, "0100" );
  CHECK( simulation_aig.simulation( simulation_aig.get_node(f3) ) == tt );
  kitty::create_from_binary_string( tt, "1001" );
  CHECK( simulation_aig.simulation( simulation_aig.get_node(f4) ) == tt );

  simulation_view simulation_aig2{ aig, xs };
  simulation_aig2.update_faninsizes();

  CHECK( simulation_aig2.depth() == 3 );
  CHECK( simulation_aig2.level( aig.get_node( a ) ) == 0 );
  CHECK( simulation_aig2.level( aig.get_node( b ) ) == 0 );
  CHECK( simulation_aig2.level( aig.get_node( f1 ) ) == 1 );
  CHECK( simulation_aig2.level( aig.get_node( f2 ) ) == 2 );
  CHECK( simulation_aig2.level( aig.get_node( f3 ) ) == 2 );
  CHECK( simulation_aig2.level( aig.get_node( f4 ) ) == 3 );

  CHECK( simulation_aig2.faninsize( aig.get_node( a ) ) == 0 );
  CHECK( simulation_aig2.faninsize( aig.get_node( b ) ) == 0 );
  CHECK( simulation_aig2.faninsize( aig.get_node( f1 ) ) == 2 );
  CHECK( simulation_aig2.faninsize( aig.get_node( f2 ) ) == 3 );
  CHECK( simulation_aig2.faninsize( aig.get_node( f3 ) ) == 3 );
  CHECK( simulation_aig2.faninsize( aig.get_node( f4 ) ) == 5 );

  kitty::create_from_binary_string( tt, "1000" );
  CHECK( simulation_aig2.simulation( simulation_aig2.get_node(f1) ) == tt );
  kitty::create_from_binary_string( tt, "0010" );
  CHECK( simulation_aig2.simulation( simulation_aig2.get_node(f2) ) == tt );
  kitty::create_from_binary_string( tt, "0100" );
  CHECK( simulation_aig2.simulation( simulation_aig2.get_node(f3) ) == tt );
  kitty::create_from_binary_string( tt, "1001" );
  CHECK( simulation_aig2.simulation( simulation_aig2.get_node(f4) ) == tt );

  CHECK( simulation_aig2.simulation( simulation_aig2.get_node(simulation_aig2.get_constant( false )) ) == const0 );


}

TEST_CASE( "compute depth, levels, fanin sizes, and simulations for AIG with inverter cost", "[simulation_view]" )
{
  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto f1 = aig.create_nand( a, b );
  const auto f2 = aig.create_nand( a, f1 );
  const auto f3 = aig.create_nand( b, f1 );
  const auto f4 = aig.create_nand( f2, f3 );
  aig.create_po( f4 );

  using TT = kitty::dynamic_truth_table;
  using Ntk = aig_network;

  simulation_view_params ps;
  ps.count_complements = true;
  simulation_view<Ntk,TT> simulation_aig{ aig, {}, ps };

  std::vector<kitty::dynamic_truth_table> xs{ 2, kitty::dynamic_truth_table( 2u ) };
  kitty::create_nth_var( xs[0], 0 );
  kitty::create_nth_var( xs[1], 1 );
  
  simulation_aig.set_input_simulations( xs );
  simulation_aig.update_simulations( );

  TT const0 = xs[0].construct();
  CHECK( simulation_aig.simulation( simulation_aig.get_node(simulation_aig.get_constant( false )) ) == const0 );

  CHECK( simulation_aig.depth() == 6 );
  CHECK( simulation_aig.level( aig.get_node( a ) ) == 0 );
  CHECK( simulation_aig.level( aig.get_node( b ) ) == 0 );
  CHECK( simulation_aig.level( aig.get_node( f1 ) ) == 1 );
  CHECK( simulation_aig.level( aig.get_node( f2 ) ) == 3 );
  CHECK( simulation_aig.level( aig.get_node( f3 ) ) == 3 );
  CHECK( simulation_aig.level( aig.get_node( f4 ) ) == 5 );

  CHECK( simulation_aig.faninsize( aig.get_node( a ) ) == 0 );
  CHECK( simulation_aig.faninsize( aig.get_node( b ) ) == 0 );
  CHECK( simulation_aig.faninsize( aig.get_node( f1 ) ) == 2 );
  CHECK( simulation_aig.faninsize( aig.get_node( f2 ) ) == 3 );
  CHECK( simulation_aig.faninsize( aig.get_node( f3 ) ) == 3 );
  CHECK( simulation_aig.faninsize( aig.get_node( f4 ) ) == 5 );

  kitty::dynamic_truth_table tt(2u);
  kitty::create_from_binary_string( tt, "1000" );
  CHECK( simulation_aig.simulation( simulation_aig.get_node(f1) ) == tt );
  kitty::create_from_binary_string( tt, "0010" );
  CHECK( simulation_aig.simulation( simulation_aig.get_node(f2) ) == tt );
  kitty::create_from_binary_string( tt, "0100" );
  CHECK( simulation_aig.simulation( simulation_aig.get_node(f3) ) == tt );
  kitty::create_from_binary_string( tt, "1001" );
  CHECK( simulation_aig.simulation( simulation_aig.get_node(f4) ) == tt );
}

TEST_CASE( "compute critical path information for simulation view", "[simulation_view]" )
{
  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto c = aig.create_pi();
  const auto d = aig.create_pi();
  const auto e = aig.create_pi();

  const auto f1 = aig.create_and( a, b );
  const auto f2 = aig.create_and( c, f1 );
  const auto f3 = aig.create_and( d, e );
  const auto f = aig.create_and( f2, f3 );
  aig.create_po( f );

  using TT = kitty::dynamic_truth_table;
  using Ntk = aig_network;

  simulation_view<Ntk,TT> simulation_aig{ aig };
  CHECK( !has_is_on_critical_path_v<decltype( aig )> );
  CHECK( has_is_on_critical_path_v<decltype( simulation_aig )> );
  CHECK( simulation_aig.is_on_critical_path( aig.get_node( a ) ) );
  CHECK( simulation_aig.is_on_critical_path( aig.get_node( b ) ) );
  CHECK( !simulation_aig.is_on_critical_path( aig.get_node( c ) ) );
  CHECK( !simulation_aig.is_on_critical_path( aig.get_node( d ) ) );
  CHECK( !simulation_aig.is_on_critical_path( aig.get_node( e ) ) );
  CHECK( simulation_aig.is_on_critical_path( aig.get_node( f1 ) ) );
  CHECK( simulation_aig.is_on_critical_path( aig.get_node( f2 ) ) );
  CHECK( !simulation_aig.is_on_critical_path( aig.get_node( f3 ) ) );
  CHECK( simulation_aig.is_on_critical_path( aig.get_node( f ) ) );
}

TEST_CASE( "compute levels, fanin sizes, and simulations during node construction", "[simulation_view]" )
{
  using TT = kitty::partial_truth_table;
  simulation_view<xag_network,TT> sxag;

  const auto a = sxag.create_pi();
  const auto b = sxag.create_pi();
  const auto c = sxag.create_pi();

  std::vector<TT> xs{ 3, kitty::partial_truth_table( 8u ) };
  kitty::create_nth_var( xs[0], 0 );
  kitty::create_nth_var( xs[1], 1 );
  kitty::create_nth_var( xs[2], 2 );
  
  sxag.set_input_simulations( xs );

  auto fo = sxag.create_xor( b, sxag.create_and( sxag.create_xor( a, b ), sxag.create_xor( b, c ) ) );
  sxag.create_po( fo );

  CHECK( sxag.depth() == 3u );
  TT tt(8u);
  kitty::create_from_binary_string( tt, "11101000" );
  CHECK( sxag.simulation( sxag.get_node(fo) ) == tt );
  CHECK( sxag.faninsize( sxag.get_node(fo) ) == 6 );

  TT const0 = xs[0].construct();
  CHECK( sxag.simulation( sxag.get_node(sxag.get_constant( false )) ) == const0 );

  simulation_view<xag_network,TT> sxag2;

  kitty::create_nth_var( xs[0], 0 );
  kitty::create_nth_var( xs[1], 1 );
  kitty::create_nth_var( xs[2], 2 );

  const auto a2 = sxag2.create_pi( xs[0] );
  const auto b2 = sxag2.create_pi( xs[1] );
  const auto c2 = sxag2.create_pi( xs[2] );
  
  auto fo2 = sxag2.create_xor( b2, sxag2.create_and( sxag2.create_xor( a2, b2 ), sxag2.create_xor( b2, c2 ) ) );
  sxag2.create_po( fo2 );

  CHECK( sxag2.depth() == 3u );
  CHECK( sxag2.simulation( sxag2.get_node(fo2) ) == tt );
  CHECK( sxag.faninsize( sxag2.get_node(fo) ) == 6 );
}