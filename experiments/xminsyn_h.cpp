#include <kitty/constructors.hpp>
#include <mockturtle/algorithms/sfps/bottomup/xminsyn_h.hpp>
#include <mockturtle/networks/aig.hpp>

using namespace mockturtle;

int main()
{
    kitty::dynamic_truth_table table( 3u );

    aig_network aig;
    const auto s1 = aig.create_pi();
    const auto s2 = aig.create_pi();
    const auto s3 = aig.create_pi();

    kitty::dynamic_truth_table x1(3u);
    kitty::dynamic_truth_table x2(3u);
    kitty::dynamic_truth_table x3(3u);
    kitty::create_nth_var( x1, 0 ); 
    kitty::create_nth_var( x2, 1 ); 
    kitty::create_nth_var( x3, 2 ); 

    table = (~x1 & ~x2 ) | ( x1 & x2 ) | ( x1 & ~x3 ) | ( x2 & x3 );

    xminsyn_h( aig, table, { s1, s2, s3 } );
}

