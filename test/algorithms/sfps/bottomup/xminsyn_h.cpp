#include <catch.hpp>

#include <kitty/constructors.hpp>
#include <mockturtle/algorithms/sfps/bottomup/xminsyn_h.hpp>
#include <mockturtle/networks/aig.hpp>

using namespace mockturtle;

TEST_CASE( "Symmetries of the majority3 AIGs", "[xminsyn_h]" )
{
    kitty::dynamic_truth_table table( 3u );
    kitty::create_majority( table );

    aig_network aig;
    const auto x1 = aig.create_pi();
    const auto x2 = aig.create_pi();
    const auto x3 = aig.create_pi();

    xminsyn_h( aig, table, { x1, x2, x3 } );
}

