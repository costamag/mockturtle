#include <catch.hpp>

#include <mockturtle/algorithms/ccgame/solvers/cusco.hpp>
#include <mockturtle/algorithms/ccgame/utils/ccg_analyzer.hpp>
#include <mockturtle/algorithms/ccgame/utils/ccg_cut.hpp>
#include <kitty/partial_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <mockturtle/algorithms/simulation.hpp>

using namespace mockturtle;
using namespace ccgame;

TEST_CASE( "cusco remapping", "[CCG]" )
{
    using TTP = kitty::partial_truth_table;
    using TTD = kitty::dynamic_truth_table;
    std::vector<TTP> xs;
    std::vector<TTP> fs;
    TTD fun_dyn( 3u );
    TTP fun_par( 8u );
    for( int i{0}; i < 3; ++i )
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
    fs = { fun_par };
    cusco<xag_network> solver( xs, fs );
    xag_network xag = solver.solve( ps );


    partial_simulator sim( xs );
    unordered_node_map<TTP, xag_network> node_to_value( xag );
    simulate_nodes( xag, node_to_value, sim );

    xag.foreach_po( [&]( const auto& x, auto index ) {
      auto sim_res = xag.is_complemented( x ) ? ~node_to_value[x] : node_to_value[x];
      CHECK( kitty::equal( sim_res, fs[index] ) );
    } );

}

TEST_CASE( "symmetry analyzer", "[CCG]" )
{
    using TTD = kitty::dynamic_truth_table;
    std::vector<TTD> xs;
    std::vector<TTD> fs;
    for( int i{0}; i < 2; ++i )
    {
        xs.emplace_back( 2u );
        kitty::create_nth_var( xs[i], i );
    }
    for( uint j{0}; j < 16u; ++j )
        fs.emplace_back(2u);

    kitty::create_from_binary_string( fs[0], "0100" ); // ES
    TTD mk = ~fs[0].construct();
    analyzer_t analyzer;
    auto syms = analyzer.find_symmetries( xs, fs[0], mk, {0, 1} );
    for( int i{0}; i < syms.size(); ++i )
        printf("%d%d->%d%d %d%d->%d%d\n", 1 & syms[i].type >> 7, 1 & syms[i].type >> 6, 1 & syms[i].type >> 5, 1 & syms[i].type >> 4, 1 & syms[i].type >> 3, 1 & syms[i].type >> 2, 1 & syms[i].type >> 1, 1 & syms[i].type );
}

TEST_CASE( "cusco set covering", "[CCG]" )
{
    using TTP = kitty::partial_truth_table;
    using TTD = kitty::dynamic_truth_table;
    std::vector<TTP> xs;
    TTD fun_dyn( 3u );
    TTP fun_par( 8u );
    for( int i{0}; i < 3; ++i )
    {
        xs.emplace_back( 8u );
        kitty::create_nth_var( xs[i], i );
    }
    kitty::create_majority( fun_dyn );
    kitty::create_from_binary_string( fun_par, kitty::to_binary( fun_dyn ) );
    std::vector<TTP> fns {fun_par}; 

    /* define the parameters */
    solver_t type = solver_t::_COV_RND;
    int nIters = 20;
    cusco_ps ps( type, nIters );
    /* solve */
    cusco<xag_network> solver( xs, fns );
    xag_network xag;
    xag = solver.solve( ps );

    partial_simulator sim( xs );
    unordered_node_map<TTP, xag_network> node_to_value( xag );
    simulate_nodes( xag, node_to_value, sim );

    xag.foreach_po( [&]( const auto& x, auto index ) {
      auto sim_res = xag.is_complemented( x ) ? ~node_to_value[x] : node_to_value[x];
      CHECK( kitty::equal( sim_res, fns[index] ) );
    } );

}