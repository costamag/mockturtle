#include <catch.hpp>

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/statistics.hpp>

using namespace kitty;
/*
TEST_CASE( "Probabilities from truth table", "[statistics_on_tt]" )
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

TEST_CASE( "Probabilities from vectors of truth tables", "[statistics_on_tt]" )
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

TEST_CASE( "Probabilities from truth table and some variables", "[statistics_on_tt]" )
{
  kitty::dynamic_truth_table f(2u);
  kitty::create_from_binary_string(f, "1100");
  std::vector<uint64_t> indeces = { 0 };
  auto probs = kitty::probability( f, indeces ); // f x0 = 00 01 10 11
  CHECK( probs[0] == 0.25 );  
  CHECK( probs[1] == 0.25 );
  CHECK( probs[2] == 0.25 );
  CHECK( probs[3] == 0.25 );  

  indeces = { 1 };
  probs = kitty::probability( f, indeces ); // f x0 = 00 00 11 11
  CHECK( probs[0] == 0.5 );  
  CHECK( probs[1] == 0 );
  CHECK( probs[2] == 0 );
  CHECK( probs[3] == 0.5 );  
}

TEST_CASE( "Probabilities from vector of truth tables and some variables", "[statistics_on_tt]" )
{
  kitty::dynamic_truth_table f1(3u);
  kitty::create_from_binary_string(f1, "11001100");
  kitty::dynamic_truth_table f2(3u);
  kitty::create_from_binary_string(f2, "01010101");
  std::vector<uint64_t> indeces = { 1,2 };
  std::vector<kitty::dynamic_truth_table> fs { f1, f2 };
  auto probs = kitty::probability( fs, indeces ); 
  // f2| 01010101 <-
  // f1| 11001100 <-
  // --------
  // x2| 11110000 <-
  // x1| 11001100 <-  
  // x0| 10101010

  CHECK( probs[0] == 0.125 );  
  CHECK( probs[1] == 0 );
  CHECK( probs[2] == 0.125 );
  CHECK( probs[3] == 0 );  
  CHECK( probs[4] == 0 );  
  CHECK( probs[5] == 0.125 );  
  CHECK( probs[6] == 0 );  
  CHECK( probs[7] == 0.125 );  
  CHECK( probs[8] == 0.125 );  
}

TEST_CASE( "Entropy from truth table", "[statistics_on_tt]" )
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

TEST_CASE( "Mutual Information from truth table", "[statistics_on_tt]" )
{
  kitty::dynamic_truth_table tt2(2u);
  kitty::create_from_binary_string(tt2, "0110");
  std::vector<uint64_t> indeces { 1 };
  auto I = kitty::mutual_information( tt2, indeces );
  CHECK( I == 0 );
  indeces = {0};
  I = kitty::mutual_information( tt2, indeces );
  CHECK( I == 0 );

  indeces = { 0, 1 };
  I = kitty::mutual_information( tt2, indeces );
  CHECK( I == 1 );

  kitty::create_from_binary_string(tt2, "1000");
  indeces = { 1 };
   I = kitty::mutual_information( tt2, indeces );
  CHECK( round(100*I) == 31 );
  indeces = {0};
  I = kitty::mutual_information( tt2, indeces );
  CHECK( round(100*I) == 31 );

  indeces = { 0, 1 };
  I = kitty::mutual_information( tt2, indeces );
  CHECK( round(100*I) == 81 );
}*/