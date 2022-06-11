#include <catch.hpp>
#include <vector>
#include <kitty/partial_truth_table.hpp>
#include <kitty/print.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/algorithms/lfe/graph_to_lfe.hpp>

#include <mockturtle/algorithms/lfe/chatterjee_method.hpp>

using namespace mockturtle;
TEST_CASE( "Create from cover: given indeces, single output, exact", "[chj]" )
{
  using TT = kitty::partial_truth_table;
  std::vector<TT> X;
  TT Y(8u);
  X.emplace_back( TT(8u) );
  kitty::create_from_binary_string( X[0], "10101010" );
  X.emplace_back( TT(8u) );
  kitty::create_from_binary_string( X[1], "11001100" );
  X.emplace_back( TT(8u) );
  kitty::create_from_binary_string( X[2], "11110000" );

  kitty::create_from_binary_string( Y, "10001000" );

  std::vector<uint64_t> indeces = { 0, 1 };
  auto chj_res = chatterjee_method( X, Y, indeces );

  CHECK( X.size() == 4 );
  CHECK( chj_res.first == "1000" );
  CHECK( kitty::to_binary(X[3]) == "10001000" );
  CHECK( chj_res.second == true );
}


TEST_CASE( "Create from cover: given indeces, single output, not exact", "[chj]" )
{
  using TT = kitty::partial_truth_table;
  std::vector<TT> X;
  TT Y(5u);
  X.emplace_back( TT(5u) );
  kitty::create_from_binary_string( X[0], "10101" );
  X.emplace_back( TT(5u) );
  kitty::create_from_binary_string( X[1], "11001" );


  kitty::create_from_binary_string( Y, "10010" );

  std::vector<uint64_t> indeces = { 0, 1 };
  auto chj_res = chatterjee_method( X, Y, indeces );


  CHECK( X.size() == 3 );
  CHECK( chj_res.second == false );
}

TEST_CASE( "Create from cover: given indeces, multi-output, not exact", "[chj]" )
{
  using TT = kitty::partial_truth_table;
  std::vector<TT> X;
  std::vector<TT> Y;
  X.emplace_back( TT(10u) );
  kitty::create_from_binary_string( X[0], "1010101001" );
  X.emplace_back( TT(10u) );
  kitty::create_from_binary_string( X[1], "1100110001" );
  X.emplace_back( TT(10u) );
  kitty::create_from_binary_string( X[2], "1111000000" );
  
  Y.emplace_back( TT(10u) );
  kitty::create_from_binary_string( Y[0], "1000100001" );
  Y.emplace_back( TT(10u) );
  kitty::create_from_binary_string( Y[1], "1000110001" );

  std::vector<uint64_t> indeces = { 0, 1 };
  auto chj_res1 = chatterjee_method( X, Y, 0, indeces );

  CHECK( X.size() == 4 );
  CHECK( chj_res1.first == "1000" );
  CHECK( kitty::to_binary(X[3]) == "1000100001" );
  CHECK( chj_res1.second == true );

  auto chj_res2 = chatterjee_method( X, Y, 1, indeces );
  CHECK( ( chj_res2.first == "1000" || chj_res2.first == "1100" ) );
  CHECK( chj_res2.second == false );
}

TEST_CASE( "apply chatterjee on network", "[chj]" )
{
  using TT = kitty::partial_truth_table;
  std::vector<TT> X;
  std::vector<TT> Y;
  klut_network klut;

  X.emplace_back( TT(10u) );
  kitty::create_from_binary_string( X[0], "1010101001" );
  
  std::vector<signal<klut_network>> xs;
  for( size_t i{0}; i < 3; ++i )
    xs.push_back( klut.create_pi() );

  X.emplace_back( TT(10u) );
  kitty::create_from_binary_string( X[1], "1100110001" );
  X.emplace_back( TT(10u) );
  kitty::create_from_binary_string( X[2], "1111000000" );
  
  Y.emplace_back( TT(10u) );
  kitty::create_from_binary_string( Y[0], "1000100001" );
  Y.emplace_back( TT(10u) );
  kitty::create_from_binary_string( Y[1], "1000110001" );

  std::vector<uint64_t> indeces = { 0, 1 };
  std::vector<signal<klut_network>> support = { xs[0], xs[1] };

  auto f0 = apply_chatterjee( klut, support, X, Y, 0, indeces );
  klut.create_po(f0);
  auto f1 = apply_chatterjee( klut, support, X, Y, 1, indeces );
  klut.create_po(f1);

  auto lfe = graph_to_lfe( klut );

  CHECK( X.size() == 5 );
  CHECK( ( kitty::to_binary(lfe.partial_MO.second[1]) == "10001000" || kitty::to_binary(lfe.partial_MO.second[1]) == "11001100" ) );
  CHECK( kitty::to_binary(lfe.partial_MO.second[0]) == "10001000"  );

}