#include <catch.hpp>

#include <kitty/kitty.hpp>

#include <mockturtle/resynthesis/aig_syn.hpp>
#include <mockturtle/resynthesis/list_simulator.hpp>
#include <mockturtle/utils/index_list.hpp>

TEST_CASE( "Synthesize two-input index lists from dynamic truth tables", "[aig_syn]" )
{
  using Tt_t = kitty::dynamic_truth_table;
  mockturtle::aig_syn_decompose<Tt_t> engine;
  mockturtle::aig_list_simulator<Tt_t> simulator;

  /* define the inputs */
  std::vector<int> ids { 0, 1 };
  std::vector<Tt_t> xs;
  for( int i = 0; i < 2; ++i )
  {
    xs.emplace_back( 2u );
    kitty::create_nth_var( xs[i], i );
  }
  // F = x
  Tt_t tt1(2u);
  kitty::create_from_binary_string( tt1, "1010" );
  auto list1 = engine( tt1 );
  CHECK( to_index_list_string( list1 ) == "{2, 1, 0, 2}" );
  Tt_t res1;
  simulator( res1, list1, ids.begin(), ids.end(), xs );
  CHECK( kitty::equal( tt1, res1 ) );

  // F = x'
  Tt_t tt2(2u);
  kitty::create_from_binary_string( tt2, "0101" );
  auto list2 = engine( tt2 );
  CHECK( to_index_list_string( list2 ) == "{2, 1, 0, 3}" );
  Tt_t res2;
  simulator( res2, list2, ids.begin(), ids.end(), xs );
  CHECK( kitty::equal( tt2, res2 ) );
  // F = xF1
  Tt_t tt3(2u);
  kitty::create_from_binary_string( tt3, "1000" );
  auto list3 = engine( tt3 );
  CHECK( to_index_list_string( list3 ) == "{2, 1, 1, 2, 4, 6}" );
  Tt_t res3;
  simulator( res3, list3, ids.begin(), ids.end(), xs );
  CHECK( kitty::equal( tt3, res3 ) );
  // F = x'F0
  Tt_t tt4(2u);
  kitty::create_from_binary_string( tt4, "0100" );
  auto list4 = engine( tt4 );
  CHECK( to_index_list_string( list4 ) == "{2, 1, 1, 3, 4, 6}" );
  Tt_t res4;
  simulator( res4, list4, ids.begin(), ids.end(), xs );
  CHECK( kitty::equal( tt4, res4 ) );
  // F = x + F0
  Tt_t tt5(2u);
  kitty::create_from_binary_string( tt5, "1011" );
  auto list5 = engine( tt5 );
  CHECK( to_index_list_string( list5 ) == "{2, 1, 1, 3, 4, 7}" );
  Tt_t res5;
  simulator( res5, list5, ids.begin(), ids.end(), xs );
  CHECK( kitty::equal( tt5, res5 ) );
  // F = x' + F1
  Tt_t tt6(2u);
  kitty::create_from_binary_string( tt6, "1101" );
  auto list6 = engine( tt6 );
  CHECK( to_index_list_string( list6 ) == "{2, 1, 1, 2, 5, 7}" );
  Tt_t res6;
  simulator( res6, list6, ids.begin(), ids.end(), xs );
  CHECK( kitty::equal( tt6, res6 ) );
}

TEST_CASE( "Synthesize four-input index lists from dynamic truth tables", "[aig_syn]" )
{
  using Tt_t = kitty::dynamic_truth_table;
  mockturtle::aig_syn_decompose<Tt_t> engine;
  mockturtle::aig_list_simulator<Tt_t> simulator;

  /* define the inputs */
  std::vector<int> ids { 0, 1, 2, 3 };
  std::vector<Tt_t> xs;
  for( int i = 0; i < 4; ++i )
  {
    xs.emplace_back( 4u );
    kitty::create_nth_var( xs[i], i );
  }

  Tt_t tt(4u);
  Tt_t res(4u);
  
  do {
    auto list = engine( tt );
    simulator( res, list, ids.begin(), ids.end(), xs );
    CHECK( kitty::equal( tt, res ) );
    kitty::next_inplace( tt );
  } while( !kitty::is_const0( tt ) );
}

TEST_CASE( "Synthesize four-input index lists from static truth tables", "[aig_syn]" )
{
  using Tt_t = kitty::static_truth_table<4u>;
  mockturtle::aig_syn_decompose<Tt_t> engine;
  mockturtle::aig_list_simulator<Tt_t> simulator;

  /* define the inputs */
  std::vector<int> ids { 0, 1, 2, 3 };
  std::vector<Tt_t> xs;
  for( int i = 0; i < 4; ++i )
  {
    xs.emplace_back();
    kitty::create_nth_var( xs[i], i );
  }

  Tt_t tt;
  Tt_t res;
  
  do {
    auto list = engine( tt );
    simulator( res, list, ids.begin(), ids.end(), xs );
    CHECK( kitty::equal( tt, res ) );
    kitty::next_inplace( tt );
  } while( !kitty::is_const0( tt ) );
}

