#include <catch.hpp>

#include <mockturtle/algorithms/ccgame/solvers/cusco.hpp>
#include <mockturtle/algorithms/ccgame/utils/ccg_analyzer.hpp>
#include <mockturtle/algorithms/ccgame/utils/ccg_supportor.hpp>
#include <mockturtle/algorithms/ccgame/utils/ccg_cut.hpp>
#include <mockturtle/algorithms/ccgame/utils/mct_node.hpp>
#include <mockturtle/algorithms/ccgame/utils/mct_method.hpp>
#include <mockturtle/algorithms/ccgame/utils/mct_tree.hpp>
#include <kitty/partial_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <mockturtle/algorithms/simulation.hpp>

using namespace mockturtle;
using namespace ccgame;

//TEST_CASE( "supportor", "[CCG]" )
//{
//    using TT = kitty::dynamic_truth_table;
//    std::vector<TT> xs;
//    std::vector<TT> fs;
//    for( int i{0}; i < 2; ++i )
//    {
//        xs.emplace_back( 2u );
//        kitty::create_nth_var( xs[i], i );
//    }
//    xs.push_back( ~xs[0] & ~xs[1] );
//    xs.push_back( ~xs[0] &  xs[1] );
//    xs.push_back(  xs[0] & ~xs[1] );
//    xs.push_back(  xs[0] &  xs[1] );
//
//    std::vector<int> costs = {1,1,1,1,1,1};
//
//    for( int i{0}; i<2; ++i )
//        fs.emplace_back( 2u );
//    fs[0] = xs[0]^xs[1];
//    fs[1] = ~xs[0]^xs[1];
//
//    /* SUPPORTER INITIALIZATION: 
//     * the only active nodes are the root node and the 
//     * identity remapping.
//     * check the correct initialization.
//    */
//    supportor_t supp( 2u, xs, costs, fs, true ); // 2u # first identity remapped entries
//    supp.print(); // remove at later stage
//    CHECK( supp.costs == std::vector<int>{1,1,1,1,1,1} );
//
//    /* information graph */
//    std::vector<TT> xxs;
//    std::vector<TT> ffs;
//    for( auto x : xs )
//        xxs.emplace_back( 4u );
//    for( auto f : fs )
//        ffs.emplace_back( 4u ); 
//    kitty::create_from_binary_string( xxs[0], "0101101001011010" );
//    kitty::create_from_binary_string( xxs[1], "0011001111001100" );
//    kitty::create_from_binary_string( xxs[2], "0001000100011110" );
//    kitty::create_from_binary_string( xxs[3], "0100101101000100" );
//    kitty::create_from_binary_string( xxs[4], "0010001011010010" );
//    kitty::create_from_binary_string( xxs[5], "0111100010001000" );
//    kitty::create_from_binary_string( ffs[0], "0110100110010110" );
//    kitty::create_from_binary_string( ffs[1], "0110100110010110" );
//
//    for( int i{0}; i<5; ++i )
//        CHECK( kitty::equal( supp.G_H[i], xxs[i] ) );
//    for( int i{0}; i<2; ++i )
//        CHECK( kitty::equal( supp.G_F[i], ffs[i] ) );
//    /* simulation patterns */
//    for( uint32_t i{0}; i<xs.size(); ++i )
//        CHECK( kitty::equal( supp.H[i], xs[i] ) );
//    for( uint32_t i{0}; i<fs.size(); ++i )
//        CHECK( kitty::equal( supp.F[i], fs[i] ) );
//
//    /* root node */
//    supportor_nd_t<TT> * pRootNd = supp.get_root_nodeP();
//    int idRoot = supp.get_root_nd_id();
//    CHECK( pRootNd->idNd == idRoot );
//    CHECK( pRootNd->idNd == 0 );
//    CHECK( pRootNd->idDv == -1 );
//    CHECK( pRootNd->nDvs == 2u+4u );
//    CHECK( pRootNd->children.size() == 2u+4u );
//    CHECK( pRootNd->pTT == nullptr );
//
//    /* identity remapped legal cut */
//    CHECK( supp.legalCuts.size() == 1 );
//    int idLegCut = 0; // identifier of the legal cut
//    std::vector<int> nds_cut = supp.get_cut( idLegCut );
//    CHECK( nds_cut == std::vector<int>{1,2} );
//    int n0 = nds_cut[0];
//    int n1 = nds_cut[1];
//    CHECK( supp.legalCuts[idLegCut].cost == -1 );
//    CHECK( supp.legalCuts[idLegCut].chromosome.sequence[0] == 3 );
//    supp.set_cost( idLegCut, 2 );
//    CHECK( supp.legalCuts[idLegCut].cost == 2 );
//    CHECK( supp.legalCuts[idLegCut].chromosome.sequence[0] == 3 );
//
//    CHECK( supp.nodes[n0].idNd == 1 ); 
//    CHECK( supp.nodes[n1].idNd == 2 );
//    CHECK( supp.nodes[n0].idNd == n0 ); 
//    CHECK( supp.nodes[n1].idNd == n1 ); 
//    supportor_nd_t<TT> nd = supp.nodes[supp.legalCuts[idLegCut].idNd];
//    CHECK( nd.is_leaf() == true );
//    CHECK( kitty::equal( *supp.nodes[n0].pTT, xs[0] ) );
//    CHECK( kitty::equal( *supp.nodes[n1].pTT, xs[1] ) );
//    CHECK( supp.nodes[n1].idNdPar == n0 );
//    CHECK( supp.nodes[n0].idNdPar == 0 );
//
//    /* Now it is initialized: generate a new solution and return the id in the solutions space */
//    chromosome_t Xn = supp.find_new( mutation_t::RND_DEL );
//}

TEST_CASE( "support generator initialization", "[CCG]" )
{
    using DTT = kitty::dynamic_truth_table;
    std::vector<DTT> xs;
    std::vector<DTT> fs;

    /* initialize the divisors */
    for( int i{0}; i < 2; ++i )
    {
        xs.emplace_back( 2u );
        kitty::create_nth_var( xs[i], i );
    }

    std::vector<divisor_t> divisors;
    for( uint32_t i{0}; i < xs.size(); ++i )
    {
        int div_id = i;
        DTT div_tt = xs[i];
        double div_area = 0;
        double div_delay = 0;
        divisor_t div( div_id, div_tt, div_area, div_delay );
        divisors.push_back(div);
    }

    /* initialize the targets */
    for( int i{0}; i<2; ++i )
        fs.emplace_back( 2u );
    fs[0] = xs[0]^xs[1];
    fs[1] = ~xs[0]^xs[1];

    std::vector<target_t> targets;
    for( uint32_t i{0}; i<fs.size(); ++i )
    {
        DTT trg_tt = fs[i];
        target_t trg( i, trg_tt );
        targets.push_back( trg );
        //targets[i].print();
    }

    /* support genenrator initialization */
    support_generator_t suppor( &divisors, &targets, gen_method_t::BASE );
   
    /* CHECKS */
    std::vector<double> area {0, 0, 1, 1, 1, 1}; 
    std::vector<double> delay {0, 0, 1, 1, 1, 1}; 
    std::vector<TT> xxs;
    std::vector<TT> ffs;
    for( int i{0}; i<6; ++i )
        xxs.emplace_back( 4u );
    for( auto f : fs )
        ffs.emplace_back( 4u ); 
    kitty::create_from_binary_string( xxs[0], "0101101001011010" );
    kitty::create_from_binary_string( xxs[1], "0011001111001100" );
    kitty::create_from_binary_string( xxs[2], "0001000100011110" );
    kitty::create_from_binary_string( xxs[3], "0010001011010010" );
    kitty::create_from_binary_string( xxs[4], "0100101101000100" );
    kitty::create_from_binary_string( xxs[5], "0111100010001000" );
    kitty::create_from_binary_string( ffs[0], "0110100110010110" );
    kitty::create_from_binary_string( ffs[1], "0110100110010110" );

    for( uint32_t i{0}; i<suppor.divisors.size(); ++i )
    {
        divisor_t div = suppor.divisors[i];
        //div.print();
        CHECK( div.area == area[i] );
        CHECK( div.delay == delay[i] );
        CHECK( div.id == (int)i );
        CHECK( kitty::equal( div.graph, xxs[i] ) );
    }

    for( uint32_t i{0u}; i<fs.size(); ++i )
    {
        target_t trg = suppor.targets[i];
        //trg.print();
        CHECK( kitty::equal( trg.tt, fs[i] ) );
        CHECK( kitty::equal( trg.graph, ffs[i] ) );
        CHECK( trg.id == (int)i );
    }
    CHECK( suppor.history.find( std::vector<int>{0,1} ) != suppor.history.end() );

    for( int i{0}; i<10; ++i )
    {
        auto sol = suppor.find_new(10);
        if( sol.size() > 0 )
        {
            suppor.store_new(sol);
            printf("size %d\n", sol.size());
            for( auto l : sol )
                printf("%d ", l);
            printf("\n");
        }
    }
}

TEST_CASE( "node of the mcts", "[CCG]" )
{
    using DTT = kitty::dynamic_truth_table;
    std::vector<double> ts;
    std::vector<DTT> xs;
    std::vector<DTT> fs;

    /* initialize the divisors */
    for( int i{0}; i < 2; ++i )
    {
        ts.push_back(0.0);
        xs.emplace_back( 2u );
        kitty::create_nth_var( xs[i], i );
    }

    /* initialize the targets */
    for( int i{0}; i<2; ++i )
        fs.emplace_back( 2u );
    fs[0] = xs[0]^xs[1];
    fs[1] = ~xs[0]^xs[1];

    mct_node_t<xag_network> root( xs, ts, fs );
    root.print();

    mct_method_ps ps( node_selection_t::RAND, 10, 1 );
    mct_method_t<mct_node_t<xag_network> > meth( ps );

    mct_tree_t<mct_node_t<xag_network> , mct_method_t> mct( root, meth );
    for( int it{0}; it<10; ++it )
        CHECK( 0 == mct.select() );
    
    /* expansion */
    printf("new node\n");
    int id = mct.expand( 0 );
    mct_node_t<xag_network>  * pNd = &mct.nodes[id];
    pNd->print();
    printf("new node\n");
    mct.simulate( pNd->id );
}


TEST_CASE( "node of the mcts: network synthesized at the root", "[CCG]" )
{
    using DTT = kitty::dynamic_truth_table;
    std::vector<double> ts;
    std::vector<DTT> xs;
    std::vector<DTT> fs;

    /* initialize the divisors */
    for( int i{0}; i < 2; ++i )
    {
        xs.emplace_back( 2u );
        ts.push_back(0.5);
        kitty::create_nth_var( xs[i], i );
    }

    /* initialize the targets */
    for( int i{0}; i<2; ++i )
        fs.emplace_back( 2u );
    fs[0] = xs[1];
    fs[1] = xs[0];

    mct_node_t<xag_network> root( xs, ts, fs );
    root.print();
    CHECK( kitty::equal( root.targets[0].tt, root.divisors[1].tt ) );
    CHECK( kitty::equal( root.targets[1].tt, root.divisors[0].tt ) );

    CHECK( root.is_leaf() == true );
    
    mct_method_ps ps( node_selection_t::RAND, 10, 1 );
    mct_method_t<mct_node_t<xag_network> > meth( ps );

    mct_tree_t<mct_node_t<xag_network> , mct_method_t> mct( root, meth );
    for( int it{0}; it<10; ++it )
        CHECK( 0 == mct.select() );
    CHECK( mct.nodes[0].TargetsDoneHere.size()==2 );
    CHECK( mct.nodes[0].TargetsDoneHere[0]==0 );
    CHECK( mct.nodes[0].TargetsDoneHere[1]==1 );
    CHECK( mct.nodes[0].TargetsDoneHere.size()==2 );
    CHECK( mct.evaluate(0) == 0 );
    CHECK( mct.nodes[0].is_leaf() );
    
    CHECK( kitty::equal(mct.nodes[0].supportor.divisors[0].tt, xs[0]) );
    CHECK( kitty::equal(mct.nodes[0].supportor.divisors[1].tt, xs[1]) );
    CHECK( kitty::equal(mct.nodes[0].supportor.divisors[2].tt, ~xs[1]&~xs[0]) );
    CHECK( kitty::equal(mct.nodes[0].supportor.divisors[3].tt, ~xs[1]& xs[0]) );
    CHECK( kitty::equal(mct.nodes[0].supportor.divisors[4].tt,  xs[1]&~xs[0]) );
    CHECK( kitty::equal(mct.nodes[0].supportor.divisors[5].tt,  xs[1]& xs[0]) );

    int iSol = mct.solve();
    CHECK( iSol == 0 );
}

TEST_CASE( "node of the mcts: network synthesized after one expansion", "[CCG]" )
{
    using DTT = kitty::dynamic_truth_table;
    std::vector<double> ts;
    std::vector<DTT> xs;
    std::vector<DTT> fs;

    /* initialize the divisors */
    for( int i{0}; i < 2; ++i )
    {
        xs.emplace_back( 2u );
        ts.push_back(0.5);
        kitty::create_nth_var( xs[i], i );
    }

    /* initialize the targets */
    for( int i{0}; i<2; ++i )
        fs.emplace_back( 2u );
    fs[0] = xs[1]&xs[0];
    fs[1] = xs[0]|xs[1];

    mct_node_t<xag_network> root( xs, ts, fs );
    root.print();

    CHECK( root.is_leaf() == false );
    
    mct_method_ps ps( node_selection_t::RAND, 10, 1 );
    mct_method_t<mct_node_t<xag_network> > meth( ps );

    mct_tree_t<mct_node_t<xag_network> , mct_method_t> mct( root, meth );
    
    CHECK( kitty::equal(mct.nodes[0].supportor.divisors[0].tt, xs[0]) );
    CHECK( kitty::equal(mct.nodes[0].supportor.divisors[1].tt, xs[1]) );
    CHECK( kitty::equal(mct.nodes[0].supportor.divisors[2].tt, ~xs[1]&~xs[0]) );
    CHECK( kitty::equal(mct.nodes[0].supportor.divisors[3].tt, ~xs[1]& xs[0]) );
    CHECK( kitty::equal(mct.nodes[0].supportor.divisors[4].tt,  xs[1]&~xs[0]) );
    CHECK( kitty::equal(mct.nodes[0].supportor.divisors[5].tt,  xs[1]& xs[0]) );
    CHECK( mct.nodes[0].supportor.divisors[2].isPo );
    CHECK( mct.nodes[0].supportor.divisors[5].isPo );

    int iSol = mct.solve();
    CHECK( iSol == 1 );
    CHECK( mct.nodes[1].TargetsDoneHere.size()==2 );
    CHECK( mct.nodes[1].TargetsDoneHere[0]==0 );
    CHECK( mct.nodes[1].TargetsDoneHere[1]==1 );
    CHECK( mct.nodes[1].TargetsDoneHere.size()==2 );
    CHECK( mct.nodes[1].is_leaf() );
    CHECK( mct.evaluate(0) == -1 );
    CHECK( mct.evaluate(1) == 2 );
}

TEST_CASE( "node of the mcts: network synthesized in the first two steps", "[CCG]" )
{
    using DTT = kitty::dynamic_truth_table;
    std::vector<double> ts;
    std::vector<DTT> xs;
    std::vector<DTT> fs;

    /* initialize the divisors */
    for( int i{0}; i < 2; ++i )
    {
        xs.emplace_back( 2u );
        ts.push_back(0.5);
        kitty::create_nth_var( xs[i], i );
    }

    /* initialize the targets */
    for( int i{0}; i<2; ++i )
        fs.emplace_back( 2u );
    fs[0] = xs[1];
    fs[1] = xs[0]|xs[1];

    mct_node_t<xag_network> root( xs, ts, fs );
    root.print();

    CHECK( root.is_leaf() == false );
    
    mct_method_ps ps( node_selection_t::RAND, 10, 1 );
    mct_method_t<mct_node_t<xag_network> > meth( ps );

    mct_tree_t<mct_node_t<xag_network> , mct_method_t> mct( root, meth );
    
    int iSol = mct.solve();
    CHECK( iSol == 1 );
    CHECK( mct.nodes[0].TargetsDoneHere.size()==1 );
    CHECK( mct.nodes[0].divisors[1].isPo );
    CHECK( mct.nodes[0].supportor.divisors[2].isPo );
    CHECK( mct.nodes[1].TargetsDoneHere.size()==1 );
    CHECK( mct.nodes[1].is_leaf() );
    CHECK( mct.evaluate(0) == -1 );
    CHECK( mct.evaluate(1) == 1 );
}

TEST_CASE( "node of the mcts: network synthesized in the second level", "[CCG]" )
{
    using DTT = kitty::dynamic_truth_table;
    std::vector<double> ts;
    std::vector<DTT> xs;
    std::vector<DTT> fs;

    /* initialize the divisors */
    for( int i{0}; i < 2; ++i )
    {
        xs.emplace_back( 2u );
        ts.push_back(0.5);
        kitty::create_nth_var( xs[i], i );
    }

    /* initialize the targets */
    for( int i{0}; i<2; ++i )
        fs.emplace_back( 2u );
    fs[0] = xs[1]^xs[0];
    fs[1] = ~xs[1]^xs[0];

    mct_node_t<xag_network> root( xs, ts, fs );
    root.print();

    mct_method_ps ps( node_selection_t::RAND, 10, 1 );
    mct_method_t<mct_node_t<xag_network> > meth( ps );

    mct_tree_t<mct_node_t<xag_network> , mct_method_t> mct( root, meth );
    
    //int iSol = mct.solve();
    int iSel = mct.select();
    CHECK( iSel == 0 );
    int iExp = mct.expand( iSel );
    CHECK( iExp == 1 );
    int iSim = mct.simulate( iSel );
    mct.path_print( iSim );

    CHECK( mct.evaluate(iSim) == 3 );

    int iSol = mct.solve();
    CHECK( mct.evaluate(iSol) == 3 );
    mct.path_print( iSol );


}

TEST_CASE( "node of the mcts: 3 inputs", "[CCG]" )
{
    using DTT = kitty::dynamic_truth_table;
    std::vector<double> ts;
    std::vector<DTT> xs;
    std::vector<DTT> fs;

    /* initialize the divisors */
    for( int i{0}; i < 3; ++i )
    {
        ts.push_back(0.0);
        xs.emplace_back( 3u );
        kitty::create_nth_var( xs[i], i );
    }

    /* initialize the targets */
    for( int i{0}; i<2; ++i )
        fs.emplace_back( 3u );
    fs[0] = xs[0]^xs[1]&xs[2];
    fs[1] = ~xs[0]^(xs[1]|xs[2]);

    mct_node_t<xag_network>  root( xs, ts, fs );
    for( int i{0}; i<xs.size(); ++i )
        CHECK( kitty::equal(root.divisors[i].tt, xs[i]) );

    mct_method_ps ps( node_selection_t::RAND, 10, 1 );
    mct_method_t<mct_node_t<xag_network> > meth( ps );
    mct_tree_t<mct_node_t<xag_network> , mct_method_t> mct( root, meth );
    int id = mct.select();
    id = mct.expand(id);
    id = mct.simulate( id );
    mct.path_print(id);
}

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