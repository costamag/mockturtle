#include <catch.hpp>

#include <kitty/kitty.hpp>

#include <mockturtle/algorithms/synth_engines/xag_synth.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/utils/index_list/index_list.hpp>

using namespace mockturtle;

TEST_CASE( "XAIG synththesizer - constants", "[xag_synth]" )
{
  xag_synth_stats st;
  constexpr uint32_t NumVars = 5u;
  xag_synth_decompose<NumVars> engine( st );
  using TT = kitty::static_truth_table<NumVars>;
  TT onset, careset;
  std::vector<uint32_t> raw;
  large_xag_index_list list;
  kitty::create_from_hex_string( onset, "0000000A" );
  kitty::create_from_hex_string( careset, "FFFFFFF0" );
  kitty::ternary_truth_table<TT> const0( onset, careset );

  engine( const0 );
  list = engine.get_list();
  raw = { 5, 1, 0, 0 };
  CHECK( raw == list.raw() );

  kitty::create_from_hex_string( onset, "FFFFFFFA" );
  kitty::ternary_truth_table<TT> const1( onset, careset );
  engine( const1 );
  list = engine.get_list();
  raw = { 5, 1, 0, 1 };
  CHECK( raw == list.raw() );
}

TEST_CASE( "XAIG synththesizer - projections", "[xag_synth]" )
{
  xag_synth_stats st;
  constexpr uint32_t NumVars = 7u;
  xag_synth_decompose<NumVars> engine( st );
  using TT = kitty::static_truth_table<NumVars>;
  TT onset, careset;
  std::vector<uint32_t> raw;
  large_xag_index_list list;
  std::array<kitty::ternary_truth_table<TT>, NumVars> proj_fns;
  for ( auto i = 0u; i < NumVars; ++i )
  {
    kitty::create_nth_var( proj_fns[i], i );
  }
  
  kitty::create_from_hex_string( onset, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAFF" );
  kitty::create_from_hex_string( careset, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00" );
  kitty::ternary_truth_table<TT> func0( onset, careset );
  engine( func0 );
  list = engine.get_list();
  raw = { 7, 1, 0, 2 };
  CHECK( raw == list.raw() );

  kitty::create_from_hex_string( onset, "555555555555555555555555555555FF" );
  kitty::ternary_truth_table<TT> func1( onset, careset );
  engine( func1 );
  list = engine.get_list();
  raw = { 7, 1, 0, 3 };
  CHECK( raw == list.raw() );
}
