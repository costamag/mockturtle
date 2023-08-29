#include <catch.hpp>

#include <algorithm>
#include <vector>

#include <kitty/constructors.hpp>
#include <mockturtle/algorithms/techaware/sym_synthesis.hpp>

using namespace mockturtle;
using namespace techaware;

TEST_CASE( "Preliminaries", "[techaware]" )
{
  techaware::TT F(3u);
  kitty::create_majority(F);
  std::vector<uint32_t> T {1, 5, 2};

  sym_synthesis<xag_network> synt( F, T );

}