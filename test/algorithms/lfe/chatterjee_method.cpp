#include <catch.hpp>

#include <mockturtle/algorithms/lfe/chatterjee_method.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <kitty/partial_truth_table.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <kitty/print.hpp>

using namespace mockturtle;


TEST_CASE( "chatterjee to learn AND from partial truth tables ", "[chatterjee]" )
{
  kitty::partial_truth_table tt1(4u);
  kitty::partial_truth_table tt2(4u);
  kitty::partial_truth_table tto(4u);
  kitty::create_from_binary_string( tt1, "1010" );
  kitty::create_from_binary_string( tt2, "1100" );
  std::vector<kitty::partial_truth_table * > ttv = { &tt1, &tt2};
  kitty::create_from_binary_string( tto, "1000" );
  chj_result res = chatterjee_method( ttv, &tto );
  CHECK( res.tt == "1000" );
}

TEST_CASE( "chatterjee to learn AND from dymanic truth tables ", "[chatterjee]" )
{
  kitty::dynamic_truth_table tt1(2u);
  kitty::dynamic_truth_table tt2(2u);
  kitty::dynamic_truth_table tto(2u);
  kitty::create_from_binary_string( tt1, "1010" );
  kitty::create_from_binary_string( tt2, "1100" );
  std::vector<kitty::dynamic_truth_table * > ttv = { &tt1, &tt2 };
  kitty::create_from_binary_string( tto, "1000" );
  chj_result res = chatterjee_method( ttv, &tto );
  CHECK( res.tt == "1000" );
}

TEST_CASE( "chatterjee to learn XOR from dymanic truth tables ", "[chatterjee]" )
{
  kitty::partial_truth_table tt1(6u);
  kitty::partial_truth_table tt2(6u);
  kitty::partial_truth_table tto(6u);
  kitty::create_from_binary_string( tt1, "101010" );
  kitty::create_from_binary_string( tt2, "110011" );
  std::vector<kitty::partial_truth_table * > ttv = { &tt1, &tt2 };
  kitty::create_from_binary_string( tto, "011001" );
  chj_result res = chatterjee_method( ttv, &tto );
  CHECK( res.tt == "0110" );
  CHECK( res.pat == tto );
}
