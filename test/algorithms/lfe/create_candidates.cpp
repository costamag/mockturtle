#include <catch.hpp>

#include <kitty/partial_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <mockturtle/algorithms/lfe/create_candidates.hpp>

using namespace kitty;
using namespace mockturtle;

TEST_CASE( "muesli f=abcd", "[candidates]" )
{
  std::vector<partial_truth_table> X;
  std::vector<partial_truth_table * > X_ptr;
  partial_truth_table Y(16u);
  partial_truth_table tt(16u);
  partial_truth_table * Y_ptr;


  for( uint32_t i = 0; i < 4; ++i )
  {
    create_nth_var(tt, i);
    X.push_back(tt);
  }

  for( uint32_t i = 0; i < 2; ++i )
    X_ptr.push_back(&X[i]);

  Y = ( X[0] & X[1] ) & ( X[2] & X[3] );
  Y_ptr = &Y;

  create_candidates_result candidates = create_candidates_method( X_ptr, Y_ptr );

  CHECK( candidates.sC_v.size() == 1 );
  CHECK( candidates.sC_v[0] == 13 );
  CHECK( candidates.tt_v[0] == "1000" );

}

TEST_CASE( "muesli f=ab^cd", "[candidates]" )
{
  std::vector<partial_truth_table> X;
  std::vector<partial_truth_table * > X_ptr;
  partial_truth_table Y(16u);
  partial_truth_table tt(16u);
  partial_truth_table * Y_ptr;


  for( uint32_t i = 0; i < 4; ++i )
  {
    create_nth_var(tt, i);
    X.push_back(tt);
  }

  for( uint32_t i = 0; i < 2; ++i )
    X_ptr.push_back(&X[i]);

  Y = ( X[0] & X[1] ) ^ ( X[2] & X[3] );
  Y_ptr = &Y;

  create_candidates_result candidates = create_candidates_method( X_ptr, Y_ptr );

  CHECK( candidates.sC_v.size() == 10 );
  CHECK( candidates.sC_v[0] == 12 );
  CHECK( candidates.tt_v[0] == "1000" );

}