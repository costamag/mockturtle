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

  using simulation_ntk = simulation_view<Ntk,TT>;

  CHECK( is_network_type_v<simulation_ntk> );
  CHECK( has_simulation_v<simulation_ntk> );
  CHECK( has_level_v<simulation_ntk> );
  CHECK( has_simulation_v<simulation_ntk> );

  using simulation_simulation_ntk = simulation_view<simulation_ntk, TT>;

  CHECK( is_network_type_v<simulation_ntk> );
  CHECK( has_simulation_v<simulation_ntk> );
  CHECK( has_level_v<simulation_ntk> );
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

TEST_CASE( "compute depth, levels, and simulations for AIG", "[simulation_view]" )
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
  
  simulation_view<Ntk,TT> simulation_aig{ aig };
  simulation_aig.set_input_simulations( xs );
  simulation_aig.update_simulations( );

  CHECK( simulation_aig.depth() == 3 );
  CHECK( simulation_aig.level( aig.get_node( a ) ) == 0 );
  CHECK( simulation_aig.level( aig.get_node( b ) ) == 0 );
  CHECK( simulation_aig.level( aig.get_node( f1 ) ) == 1 );
  CHECK( simulation_aig.level( aig.get_node( f2 ) ) == 2 );
  CHECK( simulation_aig.level( aig.get_node( f3 ) ) == 2 );
  CHECK( simulation_aig.level( aig.get_node( f4 ) ) == 3 );

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