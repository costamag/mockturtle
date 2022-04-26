#include <catch.hpp>

#include <kitty/partial_truth_table.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/static_truth_table.hpp>
#include <kitty/statistics.hpp>
#include <kitty/constructors.hpp>

using namespace kitty;

TEST_CASE( " Probabilities ", "[statistics]" )
{
  partial_truth_table tt1 ( 5 );
  create_from_binary_string( tt1, "00101" );
  auto P1 = probability( tt1 );
  CHECK( (( P1[0] == (double)3/5 ) && ( P1[1] == (double)2/5 )) == true );
}

TEST_CASE( "Probabilities from static truth table", "[statistics]" )
{
  kitty::static_truth_table<2u> tt2;
  kitty::create_from_binary_string(tt2, "1110");
  auto probs = kitty::probability(tt2);
  CHECK( probs[0] == 0.25 );
  CHECK( probs[1] == 0.75 );
}

TEST_CASE( "Probabilities from dynamic dynamic truth table", "[statistics]" )
{
  kitty::dynamic_truth_table tt2(2u);
  kitty::create_from_binary_string(tt2, "1110");
  auto probs = kitty::probability(tt2);
  CHECK( probs[0] == 0.25 );
  CHECK( probs[1] == 0.75 );

  kitty::dynamic_truth_table tt3(3u);
  kitty::create_from_binary_string(tt3, "10010110");
  probs = kitty::probability(tt3);
  CHECK( probs[0] == 0.5 );
  CHECK( probs[1] == 0.5 );
}

TEST_CASE( "Probabilities from partial truth table", "[statistics]" )
{
  kitty::partial_truth_table tt1(5);
  kitty::create_from_binary_string(tt1, "11110");
  auto probs = kitty::probability(tt1);
  CHECK( std::ceil(100*probs[1]) == 80 );
  CHECK( std::ceil(100*probs[0]) == 20 );

  kitty::partial_truth_table tt2(6);
  kitty::create_from_binary_string(tt2, "100101");
  probs = kitty::probability(tt2);
  CHECK( std::ceil(100*probs[1]) == 50 );
  CHECK( std::ceil(100*probs[0]) == 50 );

}

TEST_CASE( "Probabilities from vectors of dynamic truth tables", "[statistics]" )
{
  kitty::dynamic_truth_table tt2_0(2u);
  kitty::dynamic_truth_table tt2_1(2u);
  kitty::create_from_binary_string(tt2_0, "1110");
  kitty::create_from_binary_string(tt2_1, "1010");
  std::vector<kitty::dynamic_truth_table> tts = { tt2_0, tt2_1 };
  auto probs = kitty::probability( tts );
  CHECK( probs[0] == 0.25 );
  CHECK( probs[1] == 0.25 );
  CHECK( probs[2] == 0 );
  CHECK( probs[3] == 0.5 );

}


TEST_CASE( "Entropy from dynamic truth table", "[statistics]" )
{
  kitty::dynamic_truth_table tt2(2u);
  kitty::create_from_binary_string(tt2, "0110");
  auto H = kitty::entropy(tt2);
  CHECK( H == 1 );

  kitty::dynamic_truth_table tt3(2u);
  kitty::create_from_binary_string(tt3, "0111");
  H = kitty::entropy(tt3);
  CHECK( H == -((double)1/4)*log2((double)1/4)-((double)3/4)*log2((double)3/4) );
}
