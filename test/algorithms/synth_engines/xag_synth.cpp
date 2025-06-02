#include <catch.hpp>

#include <kitty/kitty.hpp>

#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/algorithms/synth_engines/xag_synth.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/utils/index_lists/index_list.hpp>
#include <mockturtle/utils/index_lists/simulators/list_simulator.hpp>

using namespace mockturtle;

/* remove Boolean matching with don't cares from testing */
bool constexpr xag_synth_tests_with_dcs = false;

#if xag_synth_tests_with_dcs
TEST_CASE( "XAIG synththesizer - constants", "[xag_synth]" )
{
  xag_synth_stats st;
  constexpr uint32_t NumVars = 5u;
  using TT = kitty::static_truth_table<NumVars>;
  constexpr bool UseDCs = true;
  xag_synth_decompose<UseDCs> engine( st );
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
  constexpr bool UseDCs = true;
  xag_synth_decompose<UseDCs> engine( st );
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
#endif

template<uint32_t NumVars>
void test_xag_n_input_functions()
{
  xag_synth_stats st;
  xag_synth_decompose engine( st );
  using TT = kitty::static_truth_table<NumVars>;
  kitty::static_truth_table<NumVars> onset;
  large_xag_index_list list;

  list_simulator<xag_index_list<true>, TT> sim;
  std::vector<TT> divisor_functions;
  TT tmp;
  for ( auto i = 0u; i < NumVars; ++i )
  {
    kitty::create_nth_var( tmp, i );
    divisor_functions.emplace_back( tmp );
  }
  std::vector<TT const*> xs_r;
  for ( auto i = 0u; i < NumVars; ++i )
  {
    xs_r.emplace_back( &divisor_functions[i] );
  }
  do
  {
    kitty::ternary_truth_table<TT> tt( onset );
    engine( tt );
    auto index_list = engine.get_list();

    sim( index_list, xs_r );
    TT res;
    sim.get_simulation_inline( res, index_list, xs_r, index_list.po_at( 0 ) );
    CHECK( onset == res );

    kitty::next_inplace( onset );
  } while ( !kitty::is_const0( onset ) );
}

TEST_CASE( "XAIG synththesizer - 3 input functions", "[xag_synth]" )
{
  test_xag_n_input_functions<3>();
}

template<uint32_t NumVars>
void test_xag_n_input_functions_random()
{
  xag_synth_stats st;
  xag_synth_decompose engine( st );
  using TT = kitty::static_truth_table<NumVars>;
  kitty::static_truth_table<NumVars> onset;
  large_xag_index_list list;

  list_simulator<xag_index_list<true>, TT> sim;
  std::vector<TT> divisor_functions;
  TT tmp;
  for ( auto i = 0u; i < NumVars; ++i )
  {
    kitty::create_nth_var( tmp, i );
    divisor_functions.emplace_back( tmp );
  }
  std::vector<TT const*> xs_r;
  for ( auto i = 0u; i < NumVars; ++i )
  {
    xs_r.emplace_back( &divisor_functions[i] );
  }
  int i = 1;
  do
  {
    kitty::ternary_truth_table<TT> tt( onset );
    engine( tt );
    auto index_list = engine.get_list();

    sim( index_list, xs_r );
    TT res;
    sim.get_simulation_inline( res, index_list, xs_r, index_list.po_at( 0 ) );
    CHECK( onset == res );
    kitty::create_random( onset, 2 * i++ );
  } while ( i < 1000 );
}

TEST_CASE( "XAIG synththesizer - random 10 input functions", "[xag_synth]" )
{
  test_xag_n_input_functions_random<10>();
}