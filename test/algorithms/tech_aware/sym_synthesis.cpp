#include <catch.hpp>

#include <algorithm>
#include <vector>

#include <kitty/constructors.hpp>
#include <mockturtle/algorithms/techaware/sym_synthesis.hpp>

using namespace mockturtle;
using namespace techaware;

TEST_CASE( "Input matching", "[techaware]" )
{
  techaware::TT F(3u);
  kitty::create_nth_var(F, 0);
  std::vector<uint32_t> T {1, 5, 2};

  sym_synthesis<xag_network> synt( F, T );
}

TEST_CASE( "Majority of 3", "[techaware]" )
{
  techaware::TT F(3u);
  kitty::create_majority(F);
  std::vector<uint32_t> T {1,5,2};

  sym_synthesis<xag_network> synt( F, T );
}