#include <catch.hpp>

#include <mockturtle/algorithms/bnns/decision_tree.hpp>
#include <kitty/partial_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <kitty/print.hpp>


using namespace mockturtle;
using namespace ccgame;


TEST_CASE( "initialization", "[DECTREE]" )
{
    using PTT = kitty::partial_truth_table;
    std::vector<PTT> X;
    std::vector<PTT> Y;
    for( uint32_t i{0}; i<3; ++i )
    {
        X.emplace_back( 8u );
        kitty::create_nth_var( X[i], i );
    }
    Y.push_back( X[0] & X[1] & X[2] );
    Y.push_back( X[0] ^ X[1] & X[2] );

    decision_tree dt( X, Y, X, Y );

    dt.train_random();
    dt.print();

}

