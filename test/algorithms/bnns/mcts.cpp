#include <catch.hpp>

#include <mockturtle/algorithms/mcts/decision_tree.hpp>
#include <mockturtle/algorithms/mcts/supportor.hpp>
#include <mockturtle/algorithms/mcts/nodes/nd_size.hpp>
#include <mockturtle/algorithms/mcts/method.hpp>
#include <mockturtle/algorithms/mcts/mct_tree.hpp>
#include <kitty/partial_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <mockturtle/algorithms/simulation.hpp>

using namespace mockturtle;
using namespace mcts;


TEST_CASE( "support generator initialization", "[MCTS]" )
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

    node_ps ndps;
    /* support genenrator initialization */
    support_generator_t suppor( &divisors, &targets, ndps );
   
    /* CHECKS */
    std::vector<double> area {0, 0, 1, 1, 1, 1}; 
    std::vector<double> delay {0, 0, 1, 1, 1, 1}; 
    std::vector<DTT> xxs;
    std::vector<DTT> ffs;
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

TEST_CASE( "node of the mcts", "[MCTS]" )
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

    node_ps ndps;
    nd_size_t<xag_network> root( xs, ts, fs, ndps );
    root.print();

    mct_method_ps met_ps;
    mct_method_t<nd_size_t<xag_network> > meth( met_ps );

    mct_ps mctps;
    mct_tree_t<nd_size_t<xag_network> , mct_method_t> mct( root, meth, mctps );
    for( int it{0}; it<10; ++it )
        CHECK( 0 == mct.select() );
    
    /* expansion */
    printf("new node\n");
    int id = mct.expand( 0 );
    printf( "a\n" );
    nd_size_t<xag_network>  * pNd = &mct.nodes[id];
    printf( "b\n" );
    pNd->print();
    printf("new node\n");
    mct.simulate( pNd->id );
    printf( "c\n" );

}


TEST_CASE( "node of the mcts: network synthesized at the root", "[MCTS]" )
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

    node_ps ndps;
    nd_size_t<xag_network> root( xs, ts, fs, ndps );
    root.print();
    CHECK( kitty::equal( root.targets[0].tt, root.divisors[1].tt ) );
    CHECK( kitty::equal( root.targets[1].tt, root.divisors[0].tt ) );

    CHECK( root.is_leaf() == true );
    
    mct_method_ps metps;
    mct_method_t<nd_size_t<xag_network> > meth( metps );

    mct_ps mctps;
    mct_tree_t<nd_size_t<xag_network> , mct_method_t> mct( root, meth, mctps );
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

TEST_CASE( "node of the mcts: network synthesized after one expansion", "[MCTS]" )
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

    node_ps ndps;
    nd_size_t<xag_network> root( xs, ts, fs, ndps );
    root.print();

    CHECK( root.is_leaf() == false );
    
    mct_method_ps metps;
    mct_method_t<nd_size_t<xag_network> > meth( metps );

    mct_ps mctps;
    mct_tree_t<nd_size_t<xag_network> , mct_method_t> mct( root, meth, mctps );
    
    CHECK( kitty::equal(mct.nodes[0].supportor.divisors[0].tt, xs[0]) );
    CHECK( kitty::equal(mct.nodes[0].supportor.divisors[1].tt, xs[1]) );
    CHECK( kitty::equal(mct.nodes[0].supportor.divisors[2].tt, ~xs[1]&~xs[0]) );
    CHECK( kitty::equal(mct.nodes[0].supportor.divisors[3].tt, ~xs[1]& xs[0]) );
    CHECK( kitty::equal(mct.nodes[0].supportor.divisors[4].tt,  xs[1]&~xs[0]) );
    CHECK( kitty::equal(mct.nodes[0].supportor.divisors[5].tt,  xs[1]& xs[0]) );
    CHECK( mct.nodes[0].supportor.divisors[2].isPo );
    CHECK( mct.nodes[0].supportor.divisors[5].isPo );

    int iSol = mct.solve();
    CHECK( mct.nodes[iSol].TargetsDoneHere.size()==2 );
    CHECK( mct.nodes[iSol].TargetsDoneHere[0]==0 );
    CHECK( mct.nodes[iSol].TargetsDoneHere[1]==1 );
    CHECK( mct.nodes[iSol].TargetsDoneHere.size()==2 );
    CHECK( mct.nodes[iSol].is_leaf() );
    CHECK( mct.evaluate(0) == -1 );
    CHECK( mct.evaluate(iSol) == 2 );
}

TEST_CASE( "node of the mcts: network synthesized in the first two steps", "[MCTS]" )
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

    node_ps ndps;
    nd_size_t<xag_network> root( xs, ts, fs, ndps );
    root.print();

    CHECK( root.is_leaf() == false );
    
    mct_method_ps metps;
    mct_method_t<nd_size_t<xag_network> > meth( metps );

    mct_ps mctps;
    mct_tree_t<nd_size_t<xag_network> , mct_method_t> mct( root, meth, mctps );
    
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

TEST_CASE( "node of the mcts: network synthesized in the second level", "[MCTS]" )
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

    node_ps ndps;
    nd_size_t<xag_network> root( xs, ts, fs, ndps );
    root.print();

    mct_method_ps metps;
    mct_method_t<nd_size_t<xag_network> > meth( metps );

    mct_ps mctps;
    mct_tree_t<nd_size_t<xag_network> , mct_method_t> mct( root, meth, mctps );
    
    int iSol = mct.solve();
    CHECK( mct.evaluate(iSol) == 3 );
    xag_network xag = mct.nodes[iSol].ntk;
    printf("%d\n", xag.num_gates());

}

TEST_CASE( "node of the mcts: 3 inputs", "[MCTS]" )
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

    node_ps ndps;
    nd_size_t<xag_network>  root( xs, ts, fs, ndps );
    for( int i{0}; i<xs.size(); ++i )
        CHECK( kitty::equal(root.divisors[i].tt, xs[i]) );

    mct_method_ps metps;
    mct_method_t<nd_size_t<xag_network> > meth( metps );
    
    mct_ps mctps;
    mct_tree_t<nd_size_t<xag_network> , mct_method_t> mct( root, meth, mctps );
    int id = mct.select();
    id = mct.expand(id);
    id = mct.simulate( id );
    mct.path_print(id);
}


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