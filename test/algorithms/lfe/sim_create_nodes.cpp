#include <catch.hpp>

#include <mockturtle/algorithms/lfe/sim_create_nodes.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <kitty/partial_truth_table.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <kitty/print.hpp>

using namespace mockturtle;


TEST_CASE( "sim_create_nodes to learn AND from partial truth tables ", "[sim_create_nodes]" )
{
  kitty::partial_truth_table tt1(4u);
  kitty::partial_truth_table tt2(4u);
  kitty::partial_truth_table tto(4u);
  kitty::dynamic_truth_table dtto(2u);
  kitty::create_from_binary_string( tt1, "1010" );
  kitty::create_from_binary_string( tt2, "1100" );
  std::vector<kitty::partial_truth_table * > ttv = { &tt1, &tt2};
  kitty::create_from_binary_string( tto, "1000" );

  kitty::create_from_binary_string( dtto, "1000" );
  sim_create_nodes_result res = sim_create_nodes_method( ttv, &tto );
  CHECK( res.tt_v.size() == 1 );
  CHECK( res.tt_v[0] == "1000" );
  CHECK( res.pat_v[0] == tto );
  CHECK( res.dtt_v[0] == dtto );

}

TEST_CASE( "sim_create_nodes to learn AND from dymanic truth tables ", "[sim_create_nodes]" )
{
  kitty::dynamic_truth_table tt1(2u);
  kitty::dynamic_truth_table tt2(2u);
  kitty::dynamic_truth_table tto(2u);
  kitty::create_from_binary_string( tt1, "1010" );
  kitty::create_from_binary_string( tt2, "1100" );
  std::vector<kitty::dynamic_truth_table * > ttv = { &tt1, &tt2 };
  kitty::create_from_binary_string( tto, "1000" );
  sim_create_nodes_result res = sim_create_nodes_method( ttv, &tto );
  CHECK( res.tt_v.size() == 1 );
  CHECK( res.tt_v[0] == "1000" );
  CHECK( res.dtt_v[0] == tto );
  CHECK( res.pat_v[0] == tto );
}

TEST_CASE( "sim_create_nodes to learn XOR from dymanic truth tables ", "[sim_create_nodes]" )
{
  kitty::partial_truth_table tt1(6u);
  kitty::partial_truth_table tt2(6u);
  kitty::partial_truth_table tto(6u);
  kitty::dynamic_truth_table dtto(2u);
  kitty::create_from_binary_string( tt1, "101010" );
  kitty::create_from_binary_string( tt2, "110011" );
  std::vector<kitty::partial_truth_table * > ttv = { &tt1, &tt2 };
  kitty::create_from_binary_string( tto, "011001" );
  sim_create_nodes_result res = sim_create_nodes_method( ttv, &tto );
  CHECK( res.tt_v.size() == 1 );
  CHECK( res.tt_v[0] == "0110" );
  kitty::create_from_binary_string( dtto, "0110" );
  CHECK( res.pat_v[0] == tto );
  CHECK( res.dtt_v[0] == dtto );
}

TEST_CASE( "sim create nodes with uncertainty ", "[sim_create_nodes]" )
{
  kitty::partial_truth_table tt1(5u);
  kitty::partial_truth_table tt2(5u);
  kitty::partial_truth_table tto(5u);
  kitty::partial_truth_table tto1(5u);
  kitty::partial_truth_table tto2(5u);
  kitty::dynamic_truth_table dtto1(2u);
  kitty::dynamic_truth_table dtto2(2u);
  kitty::create_from_binary_string( tt1, "11100" );
  kitty::create_from_binary_string( tt2, "11010" );
  std::vector<kitty::partial_truth_table * > ttv = { &tt1, &tt2 };
  kitty::create_from_binary_string( tto, "10001" );
  kitty::create_from_binary_string( tto1, "11001" );
  kitty::create_from_binary_string( tto2, "00001" );
  sim_create_nodes_result res = sim_create_nodes_method( ttv, &tto );
  CHECK( res.tt_v.size() == 2 );
  CHECK( res.tt_v[0] == "1001" );
  CHECK( res.tt_v[1] == "0001" );
  kitty::create_from_binary_string( dtto1, "1001" );
  kitty::create_from_binary_string( dtto2, "0001" );
  CHECK( res.pat_v[0] == tto1 );
  CHECK( res.pat_v[1] == tto2 );
  CHECK( res.dtt_v[0] == dtto1 );
  CHECK( res.dtt_v[1] == dtto2 );
}
