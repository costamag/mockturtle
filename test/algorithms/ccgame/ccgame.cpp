#include <catch.hpp>

#include <mockturtle/algorithms/ccgame/solvers/cusco.hpp>
#include <kitty/partial_truth_table.hpp>
#include <kitty/constructors.hpp>

using namespace mockturtle;
using namespace ccgame;

TEST_CASE( "cusco remapping", "[CCG]" )
{
    using TTP = kitty::partial_truth_table;
    using TTD = kitty::dynamic_truth_table;
    std::vector<TTP> xs;
    TTD fun_dyn( 3u );
    TTP fun_par( 8u );
    for( int i{0}; i < 3u; ++i )
    {
        xs.emplace_back( 8u );
        kitty::create_nth_var( xs[i], i );
    }
    kitty::create_majority( fun_dyn );
    kitty::create_from_binary_string( fun_par, kitty::to_binary( fun_dyn ) );

    /* define the parameters */
    solver_t type = solver_t::_SYM_RND;
    int nIters = 1;
    cusco_ps ps( type, nIters );
    /* solve */
    cusco<xag_network> solver( xs, {fun_par} );
    solver.solve( ps );
}

TEST_CASE( "cusco set covering", "[CCG]" )
{
    using TTP = kitty::partial_truth_table;
    using TTD = kitty::dynamic_truth_table;
    std::vector<TTP> xs;
    TTD fun_dyn( 3u );
    TTP fun_par( 8u );
    for( int i{0}; i < 3u; ++i )
    {
        xs.emplace_back( 8u );
        kitty::create_nth_var( xs[i], i );
    }
    kitty::create_majority( fun_dyn );
    kitty::create_from_binary_string( fun_par, kitty::to_binary( fun_dyn ) );

    /* define the parameters */
    solver_t type = solver_t::_COV_RND;
    int nIters = 20;
    cusco_ps ps( type, nIters );
    /* solve */
    cusco<xag_network> solver( xs, {fun_par} );
    solver.solve( ps );
}