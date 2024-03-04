#include <catch.hpp>

#include <mockturtle/utils/spfd_utils.hpp>
#include <kitty/static_truth_table.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/partial_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <cstdint>
#include <vector>

using namespace mockturtle;


TEST_CASE( "spfd manager with static truth tables", "[spfd_utils]" )
{
    using truth_table_t = kitty::static_truth_table<2u>;
    truth_table_t func;
    truth_table_t care;
    kitty::create_from_binary_string( func, "0110" );
    kitty::create_from_binary_string( care, "0111" );
    truth_table_t a, b, c;
    kitty::create_from_binary_string( a, "1010" );
    kitty::create_from_binary_string( b, "1100" );
    kitty::create_from_binary_string( c, "1110" );

    spfd_covering_manager_t<truth_table_t, 4u> manager;
    manager.init( func, care );
    CHECK( kitty::equal( func&care, manager.func[1] ) );
    CHECK( kitty::equal( ~func&care, manager.func[0] ) );

    CHECK( manager.nEdges == 2u );
    manager.update( c );
    CHECK( manager.is_covered() );
    CHECK( !manager.is_saturated() );
    manager.reset();
    CHECK( manager.nEdges == 2u );
    manager.update( a );
    manager.update( b );
    CHECK( manager.is_covered() );
    CHECK( manager.is_saturated() );
}

TEST_CASE( "spfd manager with dynamic truth tables", "[spfd_utils]" )
{
    using truth_table_t = kitty::dynamic_truth_table;
    truth_table_t func(2u);
    truth_table_t care(2u);
    kitty::create_from_binary_string( func, "0110" );
    kitty::create_from_binary_string( care, "0111" );
    truth_table_t a(2u), b(2u), c(2u);
    kitty::create_from_binary_string( a, "1010" );
    kitty::create_from_binary_string( b, "1100" );
    kitty::create_from_binary_string( c, "1110" );

    spfd_covering_manager_t<truth_table_t, 4u> manager;
    manager.init( func, care );
    CHECK( kitty::equal( func&care, manager.func[1] ) );
    CHECK( kitty::equal( ~func&care, manager.func[0] ) );

    CHECK( manager.nEdges == 2u );
    manager.update( c );
    CHECK( manager.is_covered() );
    CHECK( !manager.is_saturated() );
    manager.reset();
    CHECK( manager.nEdges == 2u );
    manager.update( a );
    manager.update( b );
    CHECK( manager.is_covered() );
    CHECK( manager.is_saturated() );
}

TEST_CASE( "spfd manager with partial truth tables", "[spfd_utils]" )
{
    using truth_table_t = kitty::partial_truth_table;
    truth_table_t func(4u);
    truth_table_t care(4u);
    kitty::create_from_binary_string( func, "0110" );
    kitty::create_from_binary_string( care, "0111" );
    truth_table_t a(4u), b(4u), c(4u);
    kitty::create_from_binary_string( a, "1010" );
    kitty::create_from_binary_string( b, "1100" );
    kitty::create_from_binary_string( c, "1110" );

    spfd_covering_manager_t<truth_table_t, 4u> manager;
    manager.init( func, care );
    CHECK( kitty::equal( func&care, manager.func[1] ) );
    CHECK( kitty::equal( ~func&care, manager.func[0] ) );

    CHECK( manager.nEdges == 2u );
    manager.update( c );
    CHECK( manager.is_covered() );
    CHECK( !manager.is_saturated() );
    manager.reset();
    CHECK( manager.nEdges == 2u );
    manager.update( a );
    manager.update( b );
    CHECK( manager.is_covered() );
    CHECK( manager.is_saturated() );
}

TEST_CASE( "spfd manager decompose input function", "[spfd_utils]" )
{
    kitty::dynamic_truth_table tt(4u);
    kitty::create_from_binary_string( tt, "1111111100000000" );
    lut_resynthesis_t<2u, 4u> resyn;
    resyn.print();
    auto lit_out = resyn.decompose( tt, 4 );
    CHECK( *lit_out == 3u );
    CHECK( resyn.funcs[*lit_out].num_vars() == 1u );
    CHECK( resyn.supps[*lit_out].size() == 1u );
}

TEST_CASE( "spfd manager decompose small function", "[spfd_utils]" )
{
    kitty::dynamic_truth_table tt(4u);
    kitty::create_from_binary_string( tt, "1111000000000000" );
    lut_resynthesis_t<2u, 4u> resyn;
    resyn.print();
    auto lit_out = resyn.decompose( tt, 4 );
    CHECK( *lit_out == 4u );
    CHECK( resyn.funcs[*lit_out].num_vars() == 2u );
    CHECK( resyn.supps[*lit_out].size() == 2u );
}

TEST_CASE( "spfd manager decompose complex function", "[spfd_utils]" )
{
    kitty::dynamic_truth_table tt(6u);
    std::vector<kitty::dynamic_truth_table> xs;
    for( int i{0}; i<6; ++i )
    {
        xs.emplace_back(6u);
        kitty::create_nth_var( xs[i], i );
    }
    kitty::create_random( tt, 5 );

    tt = tt & ~xs[3] | ~xs[2] & xs[0] | xs[1];
    lut_resynthesis_t<5u, 9u> resyn;
    resyn.print();
    auto lit_out = resyn.decompose( tt, 10 );
    resyn.print();
    printf("lit out %d\n", *lit_out);
}

TEST_CASE( "spfd manager decompose random K=3 S=4", "[spfd_utils]" )
{
    static constexpr uint32_t K = 3;
    static constexpr uint32_t S = 4;
    std::vector<uint32_t> counter = {0,0,0,0};
    for( int j{0}; j<100; ++j )
    {
        kitty::dynamic_truth_table tt(S);
        std::vector<kitty::dynamic_truth_table> xs;
        for( int i{0}; i<S; ++i )
        {
            xs.emplace_back(S);
            kitty::create_nth_var( xs[i], i );
        }

        kitty::create_random( tt, j );
        lut_resynthesis_t<K, 10u> resyn;
        auto lit_out = resyn.decompose( tt, 20 );
        if( lit_out )
        {
            CHECK( kitty::equal( resyn.sims[*lit_out], tt ) );
            
            counter[resyn.num_luts()]++;
        }
        else
        {
            printf("NOT FOUND \n");
        }
    }
    for( int i{0}; i<counter.size(); ++i )
    {
        printf("[%2d %3d]", i, counter[i] );
    } 
    printf("\n");
}

TEST_CASE( "spfd manager decompose random K=4 S=5", "[spfd_utils]" )
{
    static constexpr uint32_t K = 4;
    static constexpr uint32_t S = 5;
    std::vector<uint32_t> counter = {0,0,0,0,0,0,0,0,0};
    for( int j{0}; j<100; ++j )
    {
        kitty::dynamic_truth_table tt(S);
        std::vector<kitty::dynamic_truth_table> xs;
        for( int i{0}; i<S; ++i )
        {
            xs.emplace_back(S);
            kitty::create_nth_var( xs[i], i );
        }

        kitty::create_random( tt, j );
        lut_resynthesis_t<K, 10u> resyn;
        auto lit_out = resyn.decompose( tt, 20 );
        if( lit_out )
        {
            printf("%d\n", resyn.num_luts());
            kitty::print_binary( tt );
            printf("\n");
            kitty::print_binary( resyn.sims[*lit_out] );
            printf("\n\n");
            CHECK( kitty::equal( resyn.sims[*lit_out], tt ) );
            
            counter[resyn.num_luts()]++;
        }
        else
        {
            printf("NOT FOUND \n");
        }
    }
    for( int i{0}; i<counter.size(); ++i )
    {
        printf("[%2d %3d]", i, counter[i] );
    } 
    printf("\n");
}

TEST_CASE( "spfd manager decompose random K=4 S=6", "[spfd_utils]" )
{
    static constexpr uint32_t K = 4;
    static constexpr uint32_t S = 6;
    std::vector<uint32_t> counter = {0,0,0,0,0,0,0,0,0};
    for( int j{0}; j<100; ++j )
    {
        kitty::dynamic_truth_table tt(S);
        std::vector<kitty::dynamic_truth_table> xs;
        for( int i{0}; i<S; ++i )
        {
            xs.emplace_back(S);
            kitty::create_nth_var( xs[i], i );
        }

        kitty::create_random( tt, j );
        lut_resynthesis_t<K, 10u> resyn;
        //tt = (xs[0] & xs[1]) | kitty::cofactor1(kitty::cofactor0(tt,0),1);
        auto lit_out = resyn.decompose( tt, 20 );
        if( lit_out )
        {
            CHECK( kitty::equal( resyn.sims[*lit_out], tt ) );
            
            counter[resyn.num_luts()]++;
        }
        else
        {
            printf("NOT FOUND \n");
        }
    }
    for( int i{0}; i<counter.size(); ++i )
    {
        printf("[%2d %3d]", i, counter[i] );
    } 
    printf("\n");
}

TEST_CASE( "spfd manager decompose random K=6 S=7", "[spfd_utils]" )
{
    static constexpr uint32_t K = 6;
    static constexpr uint32_t S = 7;
    std::vector<uint32_t> counter = {0,0,0,0,0,0,0,0};
    for( int j{0}; j<100; ++j )
    {
        kitty::dynamic_truth_table tt(S);
        std::vector<kitty::dynamic_truth_table> xs;
        for( int i{0}; i<S; ++i )
        {
            xs.emplace_back(S);
            kitty::create_nth_var( xs[i], i );
        }

        kitty::create_random( tt, j );
        lut_resynthesis_t<K, 10u> resyn;
        auto lit_out = resyn.decompose( tt, 20 );
        if( lit_out )
        {
            CHECK( kitty::equal( resyn.sims[*lit_out], tt ) );
            
            counter[resyn.num_luts()]++;
        }
        else
        {
            printf("NOT FOUND \n");
        }
    }
    for( int i{0}; i<counter.size(); ++i )
    {
        printf("[%2d %3d]", i, counter[i] );
    } 
    printf("\n");
}

TEST_CASE( "spfd manager decompose problematic", "[spfd_utils]" )
{
    static constexpr uint32_t K = 4;
    static constexpr uint32_t S = 5;

        kitty::dynamic_truth_table tt(S);
        std::vector<kitty::dynamic_truth_table> xs;
        for( int i{0}; i<S; ++i )
        {
            xs.emplace_back(S);
            kitty::create_nth_var( xs[i], i );
        }

        kitty::create_from_binary_string( tt, "01100011001111111101001111001001" );
//        kitty::print_binary( tt );
//        printf("\n");
        //tt = ( ( xs[3]^xs[2] | xs[1] ) ) ^ xs[0];
        lut_resynthesis_t<K, 10u> resyn;
        //resyn.print();
        auto lit_out = resyn.decompose( tt, 20 );
        //resyn.print();
        if( lit_out )
        {
            printf("%d\n", resyn.num_luts());
            kitty::print_binary( tt );
            printf("\n");
            kitty::print_binary( resyn.sims[*lit_out] );
            printf("\n");
            CHECK( kitty::equal( resyn.sims[*lit_out], tt ) );
            if( kitty::equal( resyn.sims[*lit_out], tt ) )
            {
                printf(":)\n");
            }
            else
            {
                printf(":(((((((((\n");
            }
        }
        else
        {
            printf("NOT FOUND \n");
        }
}
