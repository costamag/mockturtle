#include <catch.hpp>

#include <mockturtle/networks/klut.hpp>
#include <mockturtle/algorithms/lfe/simulation_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <kitty/partial_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <kitty/print.hpp>

using namespace mockturtle;

TEST_CASE( "create network", "[sim_view]" )
{
  klut_network klut;
  simulation_view klut_sim{klut};

  std::vector<kitty::partial_truth_table> pats( 3 );
  kitty::partial_truth_table tta(8u);
  kitty::create_from_binary_string( tta, "10101010" );
  kitty::partial_truth_table ttb(8u);
  kitty::create_from_binary_string( ttb, "11001100" );
  kitty::partial_truth_table ttc(8u);
  kitty::create_from_binary_string( ttc, "11110000" );

  const auto a = klut_sim.create_pi( tta );
  const auto b = klut_sim.create_pi( ttb );
  const auto c = klut_sim.create_pi( ttc );

  CHECK( klut_sim.sim_patterns.size() == (2+3) );  
  CHECK( klut_sim.sim_patterns[2+0].pat == tta );
  CHECK( klut_sim.sim_patterns[2+1].pat == ttb );
  CHECK( klut_sim.sim_patterns[2+2].pat == ttc );

  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );

  CHECK( ipatterns[0].pat == tta );
  CHECK( ipatterns[1].pat == ttb );
  CHECK( ipatterns[2].pat == ttc );

  /* create unary function */
  auto f1 = klut_sim.create_not( a );
  kitty::partial_truth_table target(8u);
  kitty::create_from_binary_string( target, "01010101" );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f1)].pat == target );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f1)].sig == f1 );

  /* create binary function */
  auto f2 = klut_sim.create_and( a, b );
  kitty::create_from_binary_string( target, "10001000" );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f2)].pat == target );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f2)].sig == f2 );

  auto f3 = klut_sim.create_nand( f2, c );
  kitty::create_from_binary_string( target, "01111111" );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f3)].pat == target );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f3)].sig == f3 );

  auto f4 = klut_sim.create_or( f2, f1 );
  kitty::create_from_binary_string( target, "11011101" );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f4)].pat == target );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f4)].sig == f4 );

  auto f5 = klut_sim.create_lt( f2, f3 );
  kitty::create_from_binary_string( target, "01110111" );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f5)].pat == target );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f5)].sig == f5 );

  auto f6 = klut_sim.create_le( f2, f3 );
  kitty::create_from_binary_string( target, "01111111" );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f6)].pat == target );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f6)].sig == f6 );

  auto f7 = klut_sim.create_xor( f2, c );
  kitty::create_from_binary_string( target, "01111000" );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f7)].pat == target );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f7)].sig == f7 );

  /* create ternary function */
  auto f8 = klut_sim.create_maj( f5, f6, f7 );
  kitty::create_from_binary_string( target, "01111111" );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f8)].pat == target );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f8)].sig == f8 );  

  auto f9 = klut_sim.create_ite( f5, f6, f7 );
  kitty::create_from_binary_string( target, "01111111" );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f9)].pat == target );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f9)].sig == f9 ); 

  /* create arbitrary function */ 
  kitty::dynamic_truth_table new_tt(3u);
  kitty::create_from_binary_string( new_tt, "10000000" );
  auto f10 = klut_sim.create_node( std::vector{ f1, b, c }, new_tt );
  kitty::create_from_binary_string( target, "01000000" );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f10)].pat == target );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f10)].sig == f10 ); 

  CHECK( klut_sim.sim_patterns.size() == (2+13) );
  CHECK( ipatterns.size() == 3 );
}

TEST_CASE( "initialization", "[sim_view]" )
{
  klut_network klut1;
  auto x1 = klut1.create_pi();
  auto x2 = klut1.create_pi();
  auto x3 = klut1.create_pi();
  kitty::partial_truth_table tt1(8u);
  kitty::partial_truth_table tt2(8u);
  kitty::partial_truth_table tt3(8u);
  kitty::create_from_binary_string( tt1, "10101010" );
  kitty::create_from_binary_string( tt2, "11001100" );
  kitty::create_from_binary_string( tt3, "11110000" );
  std::vector<kitty::partial_truth_table> tts = { tt1, tt2, tt3 };

  simulation_view klut1_sim{ klut1 };
  klut1_sim.initialize_network( tts );

  CHECK( klut1_sim.sim_patterns[klut1_sim.get_node_pattern(x1)].pat == tt1 );
  CHECK( klut1_sim.sim_patterns[klut1_sim.get_node_pattern(x1)].sig == x1 ); 
  CHECK( klut1_sim.sim_patterns[klut1_sim.get_node_pattern(x2)].pat == tt2 );
  CHECK( klut1_sim.sim_patterns[klut1_sim.get_node_pattern(x2)].sig == x2 ); 
  CHECK( klut1_sim.sim_patterns[klut1_sim.get_node_pattern(x3)].pat == tt3 );
  CHECK( klut1_sim.sim_patterns[klut1_sim.get_node_pattern(x3)].sig == x3 );
  auto ipatterns1 = klut1_sim.get_input_patterns();
  CHECK( ipatterns1[klut1_sim.get_input_pattern(x1)].pat == tt1 );
  CHECK( ipatterns1[klut1_sim.get_input_pattern(x1)].sig == x1 ); 
  CHECK( ipatterns1[klut1_sim.get_input_pattern(x2)].pat == tt2 );
  CHECK( ipatterns1[klut1_sim.get_input_pattern(x2)].sig == x2 ); 
  CHECK( ipatterns1[klut1_sim.get_input_pattern(x3)].pat == tt3 );
  CHECK( ipatterns1[klut1_sim.get_input_pattern(x3)].sig == x3 );
  CHECK( ipatterns1.size() == 3 );
  CHECK( klut1_sim.sim_patterns.size() == 2+3 ) ;

  klut_network klut2;
  simulation_view klut2_sim{ klut2 };
  klut2_sim.initialize_network( tts );
  
  CHECK( klut2_sim.sim_patterns[klut2_sim.get_node_pattern(x1)].pat == tt1 );
  CHECK( klut2_sim.sim_patterns[klut2_sim.get_node_pattern(x1)].sig == x1 ); 
  CHECK( klut2_sim.sim_patterns[klut2_sim.get_node_pattern(x2)].pat == tt2 );
  CHECK( klut2_sim.sim_patterns[klut2_sim.get_node_pattern(x2)].sig == x2 ); 
  CHECK( klut2_sim.sim_patterns[klut2_sim.get_node_pattern(x3)].pat == tt3 );
  CHECK( klut2_sim.sim_patterns[klut2_sim.get_node_pattern(x3)].sig == x3 );
  auto ipatterns2 = klut2_sim.get_input_patterns();
  CHECK( ipatterns2[klut2_sim.get_input_pattern(x1)].pat == tt1 );
  CHECK( ipatterns2[klut2_sim.get_input_pattern(x1)].sig == x1 ); 
  CHECK( ipatterns2[klut2_sim.get_input_pattern(x2)].pat == tt2 );
  CHECK( ipatterns2[klut2_sim.get_input_pattern(x2)].sig == x2 ); 
  CHECK( ipatterns2[klut2_sim.get_input_pattern(x3)].pat == tt3 );
  CHECK( ipatterns2[klut2_sim.get_input_pattern(x3)].sig == x3 );
  CHECK( ipatterns2.size() == 3 );
  CHECK( klut2_sim.sim_patterns.size() == 2+3 ) ;
}

TEST_CASE( "initial simulation", "[sim_view]" )
{
  klut_network klut;
  auto x1 = klut.create_pi();
  auto x2 = klut.create_pi();
  auto x3 = klut.create_pi();
  auto f1 = klut.create_and( x1, x2 );
  auto f2 = klut.create_and( x1, x3 );
  auto f3 = klut.create_and( f1, f2 );
  klut.create_po( f3 );

  kitty::partial_truth_table tt1(8u);
  kitty::partial_truth_table tt2(8u);
  kitty::partial_truth_table tt3(8u);
  kitty::create_from_binary_string( tt1, "10101010" );
  kitty::create_from_binary_string( tt2, "11001100" );
  kitty::create_from_binary_string( tt3, "11110000" );
  std::vector<kitty::partial_truth_table> tts = { tt1, tt2, tt3 };

  simulation_view klut_sim{ klut };
  klut_sim.initialize_network( tts );

  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(x1)].pat == tt1 );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(x1)].sig == x1 ); 
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(x2)].pat == tt2 );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(x2)].sig == x2 ); 
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(x3)].pat == tt3 );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(x3)].sig == x3 );
  auto ipatterns1 = klut_sim.get_input_patterns();
  CHECK( ipatterns1[klut_sim.get_input_pattern(x1)].pat == tt1 );
  CHECK( ipatterns1[klut_sim.get_input_pattern(x1)].sig == x1 ); 
  CHECK( ipatterns1[klut_sim.get_input_pattern(x2)].pat == tt2 );
  CHECK( ipatterns1[klut_sim.get_input_pattern(x2)].sig == x2 ); 
  CHECK( ipatterns1[klut_sim.get_input_pattern(x3)].pat == tt3 );
  CHECK( ipatterns1[klut_sim.get_input_pattern(x3)].sig == x3 );
  CHECK( ipatterns1.size() == 3 );
  CHECK( klut_sim.sim_patterns.size() == 2+6 ) ;

  CHECK( klut_sim.num_gates() == 3 );
  
  CHECK( klut_sim.fanin_size(f1) == 2 );
  CHECK( klut_sim.fanin_size(f2) == 2 );
  CHECK( klut_sim.fanin_size(f3) == 2 );

}

TEST_CASE( "simulate fanin cone", "[sim_view]" )
{
  klut_network klut;
  auto x1 = klut.create_pi();
  auto x2 = klut.create_pi();
  auto x3 = klut.create_pi();
  auto f1 = klut.create_and( x1, x2 );
  auto f2 = klut.create_and( x1, x3 );
  auto f3 = klut.create_and( f1, f2 );
  klut.create_po( f3 );

  kitty::partial_truth_table tt1(8u);
  kitty::partial_truth_table tt2(8u);
  kitty::partial_truth_table tt3(8u);
  kitty::create_from_binary_string( tt1, "10101010" );
  kitty::create_from_binary_string( tt2, "11001100" );
  kitty::create_from_binary_string( tt3, "11110000" );
  std::vector<kitty::partial_truth_table> tts = { tt1, tt2, tt3 };

  simulation_view klut_sim{ klut };
  klut_sim.initialize_network( tts );
  klut_sim.simulate_fanin_cone( klut_sim.get_node(f3) );

  kitty::partial_truth_table target( 8u );
  kitty::create_from_binary_string( target, "10001000" );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f1)].pat == target );
  kitty::create_from_binary_string( target, "10100000" );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f2)].pat == target );
  kitty::create_from_binary_string( target, "10000000" );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f3)].pat == target );

  klut_sim.simulate_network( ); /* if you re-simulate you don't want the ntk to be rewritten */
  kitty::create_from_binary_string( target, "10001000" );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f1)].pat == target );
  kitty::create_from_binary_string( target, "10100000" );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f2)].pat == target );
  kitty::create_from_binary_string( target, "10000000" );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f3)].pat == target );

}

TEST_CASE( "simulate ntk", "[sim_view]" )
{
  klut_network klut;
  auto x1 = klut.create_pi();
  auto x2 = klut.create_pi();
  auto x3 = klut.create_pi();
  auto f1 = klut.create_and( x1, x2 );
  auto f2 = klut.create_and( x1, x3 );
  auto f3 = klut.create_and( f1, f2 );
  auto f4 = klut.create_or( f1, f2 );
  klut.create_po( f3 );
  klut.create_po( f4 );

  kitty::partial_truth_table tt1(8u);
  kitty::partial_truth_table tt2(8u);
  kitty::partial_truth_table tt3(8u);
  kitty::create_from_binary_string( tt1, "10101010" );
  kitty::create_from_binary_string( tt2, "11001100" );
  kitty::create_from_binary_string( tt3, "11110000" );
  std::vector<kitty::partial_truth_table> tts = { tt1, tt2, tt3 };

  simulation_view klut_sim{ klut };
  klut_sim.initialize_network( tts );
  klut_sim.simulate_fanin_cone( klut_sim.get_node(f3) );
  klut_sim.simulate_fanin_cone( klut_sim.get_node(f4) );
  
  kitty::partial_truth_table target( 8u );
  kitty::create_from_binary_string( target, "10001000" );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f1)].pat == target );
  kitty::create_from_binary_string( target, "10100000" );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f2)].pat == target );
  kitty::create_from_binary_string( target, "10000000" );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f3)].pat == target );
  kitty::create_from_binary_string( target, "10101000" );
  CHECK( klut_sim.sim_patterns[klut_sim.get_node_pattern(f4)].pat == target );
}

TEST_CASE( "simulate ntk with depth view", "[sim_view]" )
{
  klut_network klut;
  auto x1 = klut.create_pi();
  auto x2 = klut.create_pi();
  auto x3 = klut.create_pi();
  auto f1 = klut.create_and( x1, x2 );
  auto f2 = klut.create_and( x1, x3 );
  auto f3 = klut.create_and( f1, f2 );
  auto f4 = klut.create_or( f1, f2 );
  klut.create_po( f3 );
  klut.create_po( f4 );

  kitty::partial_truth_table tt1(8u);
  kitty::partial_truth_table tt2(8u);
  kitty::partial_truth_table tt3(8u);
  kitty::create_from_binary_string( tt1, "10101010" );
  kitty::create_from_binary_string( tt2, "11001100" );
  kitty::create_from_binary_string( tt3, "11110000" );
  std::vector<kitty::partial_truth_table> tts = { tt1, tt2, tt3 };

  simulation_view klut_sim{ klut };
  depth_view klut_sd{ klut_sim };

  klut_sd.initialize_network( tts );
  klut_sd.simulate_fanin_cone( klut_sd.get_node(f3) );
  klut_sd.simulate_fanin_cone( klut_sd.get_node(f4) );
  
  kitty::partial_truth_table target( 8u );
  kitty::create_from_binary_string( target, "10001000" );
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(f1)].pat == target );
  kitty::create_from_binary_string( target, "10100000" );
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(f2)].pat == target );
  kitty::create_from_binary_string( target, "10000000" );
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(f3)].pat == target );
  kitty::create_from_binary_string( target, "10101000" );
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(f4)].pat == target );
  CHECK( klut_sd.depth() == 2u );
  CHECK( klut_sd.level( klut_sd.get_node( f1 ) ) == 1u );
  CHECK( klut_sd.level( klut_sd.get_node( f2 ) ) == 1u );
  CHECK( klut_sd.level( klut_sd.get_node( f3 ) ) == 2u );
  CHECK( klut_sd.level( klut_sd.get_node( f4 ) ) == 2u );
}

TEST_CASE( "clear flags", "[sim_view]" )
{
  klut_network klut;
  auto x1 = klut.create_pi();
  auto x2 = klut.create_pi();
  auto x3 = klut.create_pi();
  auto f1 = klut.create_and( x1, x2 );
  auto f2 = klut.create_and( x1, x3 );
  auto f3 = klut.create_and( f1, f2 );
  auto f4 = klut.create_or( f1, f2 );
  klut.create_po( f3 );
  klut.create_po( f4 );

  kitty::partial_truth_table tt1(8u);
  kitty::partial_truth_table tt2(8u);
  kitty::partial_truth_table tt3(8u);
  kitty::create_from_binary_string( tt1, "10101010" );
  kitty::create_from_binary_string( tt2, "11001100" );
  kitty::create_from_binary_string( tt3, "11110000" );
  std::vector<kitty::partial_truth_table> tts = { tt1, tt2, tt3 };

  simulation_view klut_sim{ klut };
  depth_view klut_sd{ klut_sim };

  klut_sd.initialize_network( tts );
  klut_sd.simulate_fanin_cone( klut_sd.get_node(f3) );
  klut_sd.simulate_fanin_cone( klut_sd.get_node(f4) );
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(x1)].simulated == true );
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(x2)].simulated == true );
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(x3)].simulated == true );
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(f1)].simulated == true );
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(f2)].simulated == true );
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(f3)].simulated == true );
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(f4)].simulated == true );
  klut_sd.clear_simulated();
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(x1)].simulated == false );
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(x2)].simulated == false );
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(x3)].simulated == false );
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(f1)].simulated == false );
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(f2)].simulated == false );
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(f3)].simulated == false );
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(f1)].simulated == false );
  klut_sd.sim_patterns[klut_sd.nodes_to_patterns[x1]].flag = true;
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(x1)].flag == true );
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(f2)].flag == false );
  klut_sd.clear_flag();
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(x1)].flag == false );
  CHECK( klut_sd.sim_patterns[klut_sd.get_node_pattern(f2)].flag == false );

}