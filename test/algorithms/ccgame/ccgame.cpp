#include <catch.hpp>

#include <mockturtle/algorithms/ccgame/solvers/cusco.hpp>
#include <mockturtle/algorithms/ccgame/utils/ccg_analyzer.hpp>
#include <mockturtle/algorithms/ccgame/utils/ccg_supportor.hpp>
#include <mockturtle/algorithms/ccgame/utils/ccg_cut.hpp>
#include <kitty/partial_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <mockturtle/algorithms/simulation.hpp>

using namespace mockturtle;
using namespace ccgame;


TEST_CASE( "cusco remapping", "[CCG]" )
{
    using TT = kitty::dynamic_truth_table;
    std::vector<TT> xs;
    std::vector<TT> fs;
    TT fun( 3u );
    for( int i{0}; i < 3; ++i )
    {
        xs.emplace_back( 3u );
        kitty::create_nth_var( xs[i], i );
    }
    kitty::create_majority( fun );

    /* define the parameters */
    solver_t type = solver_t::_SYM_RND;
    int nIters = 1;
    cusco_ps ps( type, nIters );
    /* solve */
    fs = { fun };
    cusco<xag_network> solver( xs, fs );
    xag_network xag = solver.solve( ps ).ntk;

    default_simulator<kitty::dynamic_truth_table> sim( 3u );
    const auto sim_r = simulate<kitty::dynamic_truth_table>( xag, sim )[0];
    CHECK( kitty::equal( sim_r, fs[0] ) );
}

TEST_CASE( "symmetry analyzer", "[CCG]" )
{
    using TT = kitty::dynamic_truth_table;
    std::vector<TT> xs;
    std::vector<TT> fs;
    for( int i{0}; i < 2; ++i )
    {
        xs.emplace_back( 2u );
        kitty::create_nth_var( xs[i], i );
    }
    for( uint j{0}; j < 16u; ++j )
        fs.emplace_back(2u);

    kitty::create_from_binary_string( fs[0], "0100" ); // ES
    TT mk = ~fs[0].construct();
    analyzer_t analyzer;
    std::vector<int> vNext {0, 1};
    auto syms = analyzer.find_symmetries( &xs, &fs[0], &mk, &vNext );
    for( uint32_t i{0}; i < syms.size(); ++i )
        printf("%d%d->%d%d %d%d->%d%d\n", 1 & syms[i].type >> 7, 1 & syms[i].type >> 6, 1 & syms[i].type >> 5, 1 & syms[i].type >> 4, 1 & syms[i].type >> 3, 1 & syms[i].type >> 2, 1 & syms[i].type >> 1, 1 & syms[i].type );
}

TEST_CASE( "cusco set covering", "[CCG]" )
{
    using TT = kitty::dynamic_truth_table;
    std::vector<TT> xs;
    TT fun_dyn( 3u );
    for( int i{0}; i < 3; ++i )
    {
        xs.emplace_back( 3u );
        kitty::create_nth_var( xs[i], i );
    }
    kitty::create_majority( fun_dyn );
    std::vector<TT> fns {fun_dyn}; 

    /* define the parameters */
    solver_t type = solver_t::_COV_RND;
    int nIters = 20;
    cusco_ps ps( type, nIters );
    /* solve */
    cusco<xag_network> solver( xs, fns );
    xag_network xag;
    xag = solver.solve( ps ).ntk;


    default_simulator<kitty::dynamic_truth_table> sim( 3u );
    const auto sim_r = simulate<kitty::dynamic_truth_table>( xag, sim )[0];
    CHECK( kitty::equal( sim_r, fns[0] ) );

    //partial_simulator sim( xs );
    //unordered_node_map<TT, xag_network> node_to_value( xag );
    //simulate_nodes( xag, node_to_value, sim );

    //xag.foreach_po( [&]( const auto& x, auto index ) {
    //  auto sim_res = xag.is_complemented( x ) ? ~node_to_value[x] : node_to_value[x];
    //  CHECK( kitty::equal( sim_res, fns[index] ) );
    //} );

}

TEST_CASE( "cusco genetic set covering", "[CCG]" )
{
    using TT = kitty::dynamic_truth_table;
    std::vector<TT> xs;
    TT fun_dyn( 3u );
    for( int i{0}; i < 3; ++i )
    {
        xs.emplace_back( 3u );
        kitty::create_nth_var( xs[i], i );
    }
    kitty::create_majority( fun_dyn );
    std::vector<TT> fns {fun_dyn}; 

    /* define the parameters */
    solver_t type = solver_t::_COV_GEN;
    int nIters = 100;
    cusco_ps ps( type, nIters );
    /* solve */
    cusco<xag_network> solver( xs, fns );
    xag_network xag;
    xag = solver.solve( ps ).ntk;
    printf("nNodes=%d\n", xag.num_gates());

}

TEST_CASE( "cusco MCTS set covering", "[CCG]" )
{
    using TT = kitty::dynamic_truth_table;
    std::vector<TT> xs;
    TT fun_dyn( 3u );
    for( int i{0}; i < 3; ++i )
    {
        xs.emplace_back( 3u );
        kitty::create_nth_var( xs[i], i );
    }
    kitty::create_majority( fun_dyn );
    std::vector<TT> fns {fun_dyn}; 

    /* define the parameters */
    solver_t type = solver_t::_COV_MCTS;
    int nIters = 100;
    cusco_ps ps( type, nIters );
    /* solve */
    cusco<xag_network> solver( xs, fns );
    xag_network xag;
    xag = solver.solve( ps ).ntk;
    printf("nNodes=%d\n", xag.num_gates());

}