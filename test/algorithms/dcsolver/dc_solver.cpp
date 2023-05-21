#include <catch.hpp>

#include <mockturtle/algorithms/dcsynthesis/dc_solver.hpp>
#include <mockturtle/networks/xag.hpp>
#include <kitty/partial_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <mockturtle/algorithms/simulation.hpp>

using namespace mockturtle;

TEST_CASE( "dc solver initialization", "[DCS]" )
{
    using TT = kitty::partial_truth_table;
    std::vector<TT> xs;
    for( int i{0}; i<3; ++i )
    {
        xs.emplace_back(8u);
        kitty::create_nth_var( xs[i], i );
    }
    TT maj3(8u);
    
    kitty::create_from_binary_string( maj3, "11101000");
    std::vector<TT> fns { maj3 };
    dc_solver<xag_network> solver( xs, fns);
    xag_network xag;
    solver.solve_greedy( &xag );

  partial_simulator sim( xs );
  unordered_node_map<kitty::partial_truth_table, xag_network> node_to_value( xag );
  simulate_nodes( xag, node_to_value, sim );
  auto tts = sim.get_patterns();

  printf("J\n");
  kitty::print_binary( maj3 );
  printf("\n");
  
  std::vector<xag_network::signal> out_sigs;
  xag.foreach_po( [&]( const auto& x ) {
      out_sigs.push_back(x);
} );

TT sim_res = xag.is_complemented( out_sigs[0] ) ? ~node_to_value[out_sigs[0]] : node_to_value[out_sigs[0]];
CHECK( kitty::equal( sim_res, maj3 ) );

}

TEST_CASE( "dc solver hard function", "[DCS]" )
{
    using TT = kitty::partial_truth_table;
    std::vector<TT> xs;
    for( int i{0}; i<5; ++i )
    {
        xs.emplace_back(32u);
        kitty::create_nth_var( xs[i], i );
    }
    TT hard(32u);
    
    kitty::create_from_binary_string( hard, "01000011101110000110110000100101");
    std::vector<TT> fns { hard };
    dc_solver<xag_network> solver( xs, fns);
    xag_network xag;
    solver.solve_greedy(&xag);

  partial_simulator sim( xs );
  unordered_node_map<kitty::partial_truth_table, xag_network> node_to_value( xag );
  simulate_nodes( xag, node_to_value, sim );
  auto tts = sim.get_patterns();

  kitty::print_binary( hard );
  std::vector<xag_network::signal> out_sigs;
  xag.foreach_po( [&]( const auto& x ) {
      out_sigs.push_back(x);
} );

TT sim_res = xag.is_complemented( out_sigs[0] ) ? ~node_to_value[out_sigs[0]] : node_to_value[out_sigs[0]];
CHECK( kitty::equal( sim_res, hard ) );

}

TEST_CASE( "dc solver multi output", "[DCS]" )
{
    using TT = kitty::partial_truth_table;
    std::vector<TT> xs;
    for( int i{0}; i<3; ++i )
    {
        xs.emplace_back(8u);
        kitty::create_nth_var( xs[i], i );
    }
    TT maj3(8u);
    TT xor3(8u);
    
    kitty::create_from_binary_string( maj3, "11101000");
    kitty::create_from_binary_string( xor3, "10010110");
    std::vector<TT> fns { maj3, xor3 };
    dc_solver<xag_network> solver( xs, fns);
    xag_network xag;
    solver.solve_greedy_multioutput( &xag );

    partial_simulator sim( xs );
    unordered_node_map<kitty::partial_truth_table, xag_network> node_to_value( xag );
    simulate_nodes( xag, node_to_value, sim );
    auto tts = sim.get_patterns();
    
    std::vector<xag_network::signal> out_sigs;
    CHECK( xag.num_pos() == 2 );

    xag.foreach_po( [&]( const auto& x, auto index ) {
        TT sim_res = xag.is_complemented( x ) ? ~node_to_value[x] : node_to_value[x];
        kitty::print_binary( sim_res );
        printf(" sim %d\n", index );
        kitty::print_binary( fns[index] );
        printf(" fun %d\n", index );

        CHECK( kitty::equal( sim_res, fns[index] ) );
    } );




}
