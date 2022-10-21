#include <catch.hpp>
#include <vector>

#include <kitty/static_truth_table.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/partial_truth_table.hpp>

#include <kitty/print.hpp>
#include <mockturtle/networks/klut.hpp>

#include <mockturtle/algorithms/sfps/nodes_creation.hpp>

using namespace mockturtle;

TEST_CASE( "Create from cover: static truth table", "[node_creation]" )
{
  using TT = kitty::static_truth_table<3u>;
  std::vector<TT> X;
  TT Y;
  X.emplace_back();
  kitty::create_nth_var( X[0], 0 );
  X.emplace_back();
  kitty::create_nth_var( X[1], 1 );
  X.emplace_back();
  kitty::create_nth_var( X[2], 2 );

  kitty::create_from_binary_string( Y, "10001000" );

  std::vector<TT*> Xptr = { &X[0], &X[1] };
  auto chj_res = chatterjee_method<TT>( Xptr, &Y );

  CHECK( X.size() == 3 );
  CHECK( chj_res.tt == "1000" );
  CHECK( kitty::to_binary(chj_res.pat) == "10001000" );
}

TEST_CASE( "Create from cover: dynamic truth table", "[node_creation]" )
{
  using TT = kitty::dynamic_truth_table;
  std::vector<TT> X;
  TT Y(3u);
  X.emplace_back( TT(3u) );
  kitty::create_from_binary_string( X[0], "10101010" );
  X.emplace_back( TT(3u) );
  kitty::create_from_binary_string( X[1], "11001100" );
  X.emplace_back( TT(3u) );
  kitty::create_from_binary_string( X[2], "11110000" );

  kitty::create_from_binary_string( Y, "10001000" );

  std::vector<TT*> Xptr = { &X[0], &X[1] };
  auto chj_res = chatterjee_method<TT>( Xptr, &Y );

  CHECK( X.size() == 3 );
  CHECK( chj_res.tt == "1000" );
  CHECK( kitty::to_binary(chj_res.pat) == "10001000" );
}

TEST_CASE( "Create from cover: partial truth table", "[node_creation]" )
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

  std::vector<TT*> Xptr = { &X[0], &X[1] };
  auto chj_res = chatterjee_method<TT>( Xptr, &Y );

  CHECK( X.size() == 3 );
  CHECK( chj_res.tt == "1000" );
  CHECK( kitty::to_binary(chj_res.pat) == "10001000" );
}


TEST_CASE( "Create from cover: given indeces, single output, not exact", "[node_creation]" )
{
  using TT = kitty::partial_truth_table;
  std::vector<TT> X;
  TT Y(5u);
  X.emplace_back( TT(5u) );
  kitty::create_from_binary_string( X[0], "10101" );
  X.emplace_back( TT(5u) );
  kitty::create_from_binary_string( X[1], "11001" );

  kitty::create_from_binary_string( Y, "10010" );

  std::vector<TT*> Xptr = { &X[0], &X[1] };
  auto chj_res = chatterjee_method<TT>( Xptr, &Y );

  CHECK( (chj_res.tt == "0001" || chj_res.tt == "1001") );
}

TEST_CASE( "Create from dynamic truth table cover", "[node_creation]" )
{
  using TT = kitty::dynamic_truth_table;
  std::vector<TT> X;
  TT Y(5u);
  X.emplace_back( TT(5u) );
  kitty::create_nth_var( X[0], 0 );
  X.emplace_back( TT(5u) );
  kitty::create_nth_var( X[1], 1 );
  X.emplace_back( TT(5u) );
  kitty::create_nth_var( X[2], 2 );
  X.emplace_back( TT(5u) );
  kitty::create_nth_var( X[3], 3 );
  X.emplace_back( TT(5u) );
  kitty::create_nth_var( X[4], 4 );

  Y = X[0] ^ X[1] ^ X[2] ^ X[3] ^ X[4];

  chatterjee_method_params ps;
  ps.seed = 42;
  std::vector<TT*> Xptr = { &X[0], &X[1] };
  auto chj_res1 = chatterjee_method<TT>( Xptr, &Y, ps );
  
  ps.seed = 43;
  auto chj_res2 = chatterjee_method<TT>( Xptr, &Y, ps );
  
  ps.seed = 44;
  auto chj_res3 = chatterjee_method<TT>( Xptr, &Y, ps );

  CHECK( ( ( chj_res1.tt != chj_res2.tt ) && ( chj_res1.tt != chj_res3.tt ) && ( chj_res2.tt != chj_res3.tt ) ) );
}

TEST_CASE( "Create static truth table from cover using nodes enumeration", "[node_creation]" )
{
  using TT = kitty::static_truth_table<3u>;
  std::vector<TT> X;
  TT Y;
  X.emplace_back();
  kitty::create_nth_var( X[0], 0 );
  X.emplace_back();
  kitty::create_nth_var( X[1], 1 );
  X.emplace_back();
  kitty::create_nth_var( X[2], 2 );

  kitty::create_from_binary_string( Y, "10001000" );

  std::vector<TT*> Xptr = { &X[0], &X[1] };
  auto nen_res = nodes_enumeration<TT>( Xptr, &Y );

  CHECK( nen_res.tt_v[0] == "1000" );
  CHECK( kitty::to_binary(nen_res.pat_v[0]) == "10001000" );
}

TEST_CASE( "Create dynamic truth table from cover using nodes enumeration", "[node_creation]" )
{
  using TT = kitty::dynamic_truth_table;
  std::vector<TT> X;
  TT Y(3u);
  X.emplace_back( TT(3u) );
  kitty::create_from_binary_string( X[0], "10101010" );
  X.emplace_back( TT(3u) );
  kitty::create_from_binary_string( X[1], "11001100" );
  X.emplace_back( TT(3u) );
  kitty::create_from_binary_string( X[2], "11110000" );

  kitty::create_from_binary_string( Y, "10001000" );

  std::vector<TT*> Xptr = { &X[0], &X[1] };
  auto nen_res = nodes_enumeration<TT>( Xptr, &Y );

  CHECK( nen_res.tt_v[0] == "1000" );
  CHECK( kitty::to_binary(nen_res.pat_v[0]) == "10001000" );
}

TEST_CASE( "Create partial truth table from cover using nodes enumeration", "[node_creation]" )
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

  std::vector<TT*> Xptr = { &X[0], &X[1] };
  auto nen_res = nodes_enumeration<TT>( Xptr, &Y );

  CHECK( nen_res.tt_v[0] == "1000" );
  CHECK( kitty::to_binary(nen_res.pat_v[0]) == "10001000" );
}

TEST_CASE( "Create partial truth table with nodes enumeration in presence of alternative", "[node_creation]" )
{
  using TT = kitty::partial_truth_table;
  std::vector<TT> X;
  TT Y(6u);
  X.emplace_back( TT(6u) );
  kitty::create_from_binary_string( X[0], "101011" );
  X.emplace_back( TT(6u) );
  kitty::create_from_binary_string( X[1], "110011" );

  kitty::create_from_binary_string( Y, "100101" );

  std::vector<TT*> Xptr = { &X[0], &X[1] };
  auto nen_res = nodes_enumeration<TT>( Xptr, &Y );

  CHECK( (nen_res.tt_v[0] == "1001") );
  CHECK( ( kitty::to_binary(nen_res.pat_v[0]) == "100111") );
  CHECK( (nen_res.tt_v[1] == "0001") );
  CHECK( ( kitty::to_binary(nen_res.pat_v[1]) == "000100") );

}