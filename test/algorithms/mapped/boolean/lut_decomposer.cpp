#include <catch.hpp>

#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>

#include <mockturtle/algorithms/mapped/boolean/lut_decomposer.hpp>

TEST_CASE( "Termination condition for LUT decomposition", "[lut_decomposer]" )
{
  static constexpr uint32_t MaxNumVars = 2u;
  static constexpr uint32_t MaxCutSize = 3u;
  using CSTT = kitty::static_truth_table<MaxCutSize>;
  using ISTT = kitty::ternary_truth_table<CSTT>;
  mockturtle::lut_decomposer<MaxCutSize, MaxNumVars> decomposer;

  CSTT care1, mask1;
  kitty::create_from_binary_string( care1, "10001100" );
  kitty::create_from_binary_string( mask1, "11111011" );
  ISTT func1( care1, mask1 );
  std::vector<double> times{ 0.0, 0.0, 0.0 };
  CHECK( decomposer.run( func1, times ) );
  decomposer.foreach_spec( [&]( auto const& tt ) {
    kitty::static_truth_table<2u> expected;
    kitty::create_from_binary_string( expected, "1000" );
    CHECK( kitty::equal( expected, tt._bits ) );
    CHECK( kitty::is_const0( ~tt._care ) );\
    return true; },
                           [&]( auto const& sim_ptrs ) {
                             return *sim_ptrs[0] & *sim_ptrs[1];
                           } );
}