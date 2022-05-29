#include <catch.hpp>

#include <mockturtle/algorithms/lfe/sim_decomposition.hpp>
#include <mockturtle/algorithms/lfe/simulation_view.hpp>
#include <mockturtle/algorithms/lfe/sim_muesli.hpp>
#include <mockturtle/algorithms/lfe/muesli.hpp>
#include <mockturtle/algorithms/lfe/chatterjee_method.hpp>
#include <mockturtle/networks/klut.hpp>
#include <kitty/partial_truth_table.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <kitty/print.hpp>

using namespace mockturtle;


TEST_CASE( "decompose f = ab+cde ", "[sim_dec]" )
{
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 5; ++i )
  {
    ex.push_back( kitty::partial_truth_table(32u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(32u);
  y = ( ex[0] & ex[1] ) | ( ( ex[2] & ex[3] ) & ex[4] );

  sim_decomposition_params ps;
  ps.verbose = true;

  auto f0 = sim_decomposition( klut_sim, ex, y, ps );
  klut_sim.create_po( f0 );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 5 );
  for( size_t i = 0; i < 5; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[klut_sim.nodes_to_patterns[klut_sim.get_node(f0)]].pat == y );

}

TEST_CASE( "decompose F = {ab+cde, ab+cd} ", "[sim_dec]" )
{
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 5; ++i )
  {
    ex.push_back( kitty::partial_truth_table(32u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y1(32u);
  kitty::partial_truth_table y2(32u);
  y1 = ( ex[0] & ex[1] ) | ( ( ex[2] & ex[3] ) & ex[4] );
  y2 = ( ex[0] & ex[1] ) | ( ex[2] & ex[3] );
  std::vector<kitty::partial_truth_table> Y = { y1, y2 };
  sim_decomposition_params ps;
  ps.verbose = true;

  auto osignals = sim_decomposition( klut_sim, ex, Y, ps );
  klut_sim.create_po( osignals[0] );
  klut_sim.create_po( osignals[1] );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 5 );
  for( size_t i = 0; i < 5; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  for( size_t i = 0; i < osignals.size(); ++i )
    CHECK( klut_sim.sim_patterns[klut_sim.nodes_to_patterns[klut_sim.get_node( osignals[i] )]].pat == Y[i] );
  
}


TEST_CASE( "sim muesli create network f = ab+cde ", "[sim_muesli]" )
{
  std::cout << "sim muesli : f = ab+cde" << std::endl;
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 5; ++i )
  {
    ex.push_back( kitty::partial_truth_table(32u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(32u);
  y = ( ex[0] & ex[1] ) | ( ( ex[2] & ex[3] ) & ex[4] );
  
  sim_muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 2;
  ps.eps_th = 0.99;
  ps.verbose = true;

  auto f0 = sim_muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 5 );
  for( size_t i = 0; i < 5; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  CHECK( klut_sim.num_gates() == 4 );

}

TEST_CASE( "sim muesli create network f = ab+cd ", "[sim_muesli]" )
{
  std::cout << "sim muesli : f = ab+cd" << std::endl;

  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 4; ++i )
  {
    ex.push_back( kitty::partial_truth_table(32u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(32u);
  y = ( ex[0] & ex[1] ) | ( ex[2] & ex[3] );
  
  sim_muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 2;
  ps.eps_th = 0.99;
  ps.verbose = true;

  auto f0 = sim_muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 4 );
  for( size_t i = 0; i < 4; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  CHECK( klut_sim.num_gates() == 4 );


}

TEST_CASE( "sim muesli create network f = abcd |A|=3", "[sim_muesli]" )
{
  std::cout << "sim muesli : f = abcd" << std::endl;

  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 4; ++i )
  {
    ex.push_back( kitty::partial_truth_table(16u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(16u);
  y = ( ex[0] & ex[1] ) & ( ex[2] & ex[3] );
  
  sim_muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 2;
  ps.eps_th = 0.99;
  ps.verbose = true;

  auto f0 = sim_muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 4 );
  for( size_t i = 0; i < 4; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 3 );

}


TEST_CASE( "sim muesli create network f = a^(bcd) |A|=3", "[sim_muesli]" )
{
  std::cout << "sim muesli : f = a^(bcd)" << std::endl;

  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 4; ++i )
  {
    ex.push_back( kitty::partial_truth_table(16u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(16u);
  y =  ex[0] ^ (ex[1]  & ( ex[2] & ex[3] ));
  
  sim_muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 2;
  ps.eps_th = 0.99;
  ps.verbose = true;

  auto f0 = sim_muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 4 );
  for( size_t i = 0; i < 4; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 3 );

}

TEST_CASE( "sim muesli create network f = a+b^c+d |A|=3", "[sim_muesli]" )
{
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 4; ++i )
  {
    ex.push_back( kitty::partial_truth_table(16u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(16u);
  y =  ex[0] | ( (ex[1] ^ ex[2]) | ex[3] );
  
  sim_muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = sim_muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 4 );
  for( size_t i = 0; i < 4; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 3 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 3 );


}

TEST_CASE( "sim muesli create network f = a^(b^c)+d |A|=3", "[sim_muesli]" )
{
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 4; ++i )
  {
    ex.push_back( kitty::partial_truth_table(16u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(16u);
  y =  (ex[0] ^  (ex[1] | ex[2])) | ex[3];
  
  sim_muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = sim_muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 4 );
  for( size_t i = 0; i < 4; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 5 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 3 );
 
}

TEST_CASE( "sim muesli create network and3 = xyz", "[sim_muesli]" )
{
  std::cout << "AND3" << std::endl;

  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 3; ++i )
  {
    ex.push_back( kitty::partial_truth_table(8u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(8u);
  y = ex[0] & ( ex[1] & ex[2] );
  
  sim_muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = sim_muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );
  for( size_t i = 0; i < 3; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 2 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 2 );

}

TEST_CASE( "sim muesli create network XorAnd = x(y^z)", "[sim_muesli]" )
{
  std::cout << "XORAND" << std::endl;

  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 3; ++i )
  {
    ex.push_back( kitty::partial_truth_table(8u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(8u);
  y = ex[0] & ( ex[1] ^ ex[2] );
  
  sim_muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = sim_muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );
  for( size_t i = 0; i < 3; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 2 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 2 );

}

TEST_CASE( "sim muesli create network orAnd = x(y|z)", "[sim_muesli]" )
{
  std::cout << "ORAND" << std::endl;
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 3; ++i )
  {
    ex.push_back( kitty::partial_truth_table(8u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(8u);
  y = ex[0] & ( ex[1] | ex[2] );
  
  sim_muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = sim_muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );
  for( size_t i = 0; i < 3; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 2 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 2 );

}


TEST_CASE( "sim muesli create network OneHot ", "[sim_muesli]" )
{
  std::cout << " sim muesli one hot" << std::endl;
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 3; ++i )
  {
    ex.push_back( kitty::partial_truth_table(8u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(8u);
  y = ( ( ex[0] &( ~ex[1] & ~ ex[2] ) ) ^ ( ~ex[0] &( ex[1] & ~ ex[2] ) ) )^ ( ~ex[0] &( ~ex[1] & ex[2] ) ) ;
  
  sim_muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = sim_muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );
  for( size_t i = 0; i < 3; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 7 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 6 );
}

TEST_CASE( "sim muesli create network majority ", "[sim_muesli]" )
{
  std::cout << " sim muesli maj" << std::endl;

  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 3; ++i )
  {
    ex.push_back( kitty::partial_truth_table(8u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(8u);
  y = ( ( ex[0] & ex[1] ) ^ ( ex[1] & ex[2] ) ) ^ ( ex[0] & ex[2] );
  
  sim_muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = sim_muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );
  for( size_t i = 0; i < 3; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 4 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 6 );
}

TEST_CASE( "sim muesli create network gamble ", "[sim_muesli]" )
{
  std::cout << " sim muesli gamble" << std::endl;

  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 3; ++i )
  {
    ex.push_back( kitty::partial_truth_table(8u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(8u);
  y = ( ( ex[0] & ex[1] ) & ex[2] ) ^ ( ( ~ex[0] & ~ex[1] ) & ~ex[2] );
  
  sim_muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = sim_muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );
  for( size_t i = 0; i < 3; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 3 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 6 );
}

TEST_CASE( "sim muesli create network mux ", "[sim_muesli]" )
{
  std::cout << " sim muesli mux" << std::endl;

  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 3; ++i )
  {
    ex.push_back( kitty::partial_truth_table(8u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(8u);
  y = ( ex[0] & ex[1] ) ^ ( ~ex[0] & ex[2] );
  
  sim_muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = sim_muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );
  for( size_t i = 0; i < 3; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 4 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 6 );
}

TEST_CASE( "sim muesli create network andxor ", "[sim_muesli]" )
{
  std::cout << " sim muesli andxor" << std::endl;

  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 3; ++i )
  {
    ex.push_back( kitty::partial_truth_table(8u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(8u);
  y = ex[0] ^  ( ex[1] & ex[2] );
  
  sim_muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = sim_muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );
  for( size_t i = 0; i < 3; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 5 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 2 );
}

TEST_CASE( "sim muesli create network xor3 ", "[sim_muesli]" )
{
  std::cout << " sim muesli xor3" << std::endl;

  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 3; ++i )
  {
    ex.push_back( kitty::partial_truth_table(8u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(8u);
  y = ex[0] ^ ( ex[1] ^ ex[2] );
  
  sim_muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = sim_muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );
  for( size_t i = 0; i < 3; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 2 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 2 );
}


TEST_CASE( "create network f = ab+cde ", "[muesli]" )
{
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 5; ++i )
  {
    ex.push_back( kitty::partial_truth_table(32u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(32u);
  y = ( ex[0] & ex[1] ) | ( ( ex[2] & ex[3] ) & ex[4] );
  
  muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 2;
  ps.eps_th = 0.99;
  ps.verbose = true;

  auto f0 = muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 5 );
  for( size_t i = 0; i < 5; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  CHECK( klut_sim.num_gates() == 4 );

}

TEST_CASE( "create network f = ab+cd ", "[muesli]" )
{
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 4; ++i )
  {
    ex.push_back( kitty::partial_truth_table(32u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(32u);
  y = ( ex[0] & ex[1] ) | ( ex[2] & ex[3] );
  
  muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 2;
  ps.eps_th = 0.99;
  ps.verbose = true;

  auto f0 = muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 4 );
  for( size_t i = 0; i < 4; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );


}

TEST_CASE( "create network f = ab+cde |A|=3", "[muesli]" )
{
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 4; ++i )
  {
    ex.push_back( kitty::partial_truth_table(16u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(16u);
  y = ( ex[0] & ex[1] ) & ( ex[2] & ex[3] );
  
  muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 2;
  ps.eps_th = 0.99;
  ps.verbose = true;

  auto f0 = muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 4 );
  for( size_t i = 0; i < 4; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 0 );

}


TEST_CASE( "muesli create network dot = x^(z | xy)", "[muesli]" )
{
  std::cout << " #######################################################" << std::endl;
  std::cout << " #######################################################" << std::endl;
  std::cout << " #######################################################" << std::endl;
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 3; ++i )
  {
    ex.push_back( kitty::partial_truth_table(8u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(8u);
  y = ex[0] ^ ( ex[2] | ( ex[0] & ex[1] ) ) ;
  
  muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );
  for( size_t i = 0; i < 3; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 4 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 6 );

}

TEST_CASE( "muesli create network and3 = xyz", "[muesli]" )
{
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 3; ++i )
  {
    ex.push_back( kitty::partial_truth_table(8u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(8u);
  y = ex[0] & ( ex[1] & ex[2] );
  
  muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );
  for( size_t i = 0; i < 3; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 2 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 2 );

}

TEST_CASE( "muesli create network XorAnd = x(y^z)", "[muesli]" )
{
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 3; ++i )
  {
    ex.push_back( kitty::partial_truth_table(8u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(8u);
  y = ex[0] & ( ex[1] ^ ex[2] );
  
  muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );
  for( size_t i = 0; i < 3; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 2 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 2 );

}

TEST_CASE( "muesli create network XorAnd = x(y|z)", "[muesli]" )
{
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 3; ++i )
  {
    ex.push_back( kitty::partial_truth_table(8u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(8u);
  y = ex[0] & ( ex[1] | ex[2] );
  
  muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );
  for( size_t i = 0; i < 3; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 2 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 2 );

}


TEST_CASE( "muesli create network OneHot ", "[muesli]" )
{
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 3; ++i )
  {
    ex.push_back( kitty::partial_truth_table(8u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(8u);
  y = ( ( ex[0] &( ~ex[1] & ~ ex[2] ) ) ^ ( ~ex[0] &( ex[1] & ~ ex[2] ) ) )^ ( ~ex[0] &( ~ex[1] & ex[2] ) ) ;
  
  muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );
  for( size_t i = 0; i < 3; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 6 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 6 );
}

TEST_CASE( "muesli create network majority ", "[muesli]" )
{
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 3; ++i )
  {
    ex.push_back( kitty::partial_truth_table(8u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(8u);
  y = ( ( ex[0] & ex[1] ) ^ ( ex[1] & ex[2] ) ) ^ ( ex[0] & ex[2] );
  
  muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );
  for( size_t i = 0; i < 3; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 5 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 6 );
}

TEST_CASE( "muesli create network gamble ", "[muesli]" )
{
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 3; ++i )
  {
    ex.push_back( kitty::partial_truth_table(8u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(8u);
  y = ( ( ex[0] & ex[1] ) & ex[2] ) ^ ( ( ~ex[0] & ~ex[1] ) & ~ex[2] );
  
  muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );
  for( size_t i = 0; i < 3; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 5 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 6 );
}

TEST_CASE( "muesli create network mux ", "[muesli]" )
{
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 3; ++i )
  {
    ex.push_back( kitty::partial_truth_table(8u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(8u);
  y = ( ex[0] & ex[1] ) ^ ( ~ex[0] & ex[2] );
  
  muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );
  for( size_t i = 0; i < 3; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 6 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 6 );
}

TEST_CASE( "muesli create network andxor ", "[muesli]" )
{
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 3; ++i )
  {
    ex.push_back( kitty::partial_truth_table(8u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(8u);
  y = ex[0] ^  ( ex[1] & ex[2] );
  
  muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );
  for( size_t i = 0; i < 3; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 2 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 2 );
}

TEST_CASE( "muesli create network xor3 ", "[muesli]" )
{
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 3; ++i )
  {
    ex.push_back( kitty::partial_truth_table(8u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(8u);
  y = ex[0] ^ ( ex[1] ^ ex[2] );
  
  muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );
  for( size_t i = 0; i < 3; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 6 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 2 );
}

TEST_CASE( "sim muesli create network dot = x^(z | xy)", "[sim_muesli]" )
{
  std::cout << "DOT" << std::endl;
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 3; ++i )
  {
    ex.push_back( kitty::partial_truth_table(8u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(8u);
  y = ex[0] ^ ( ex[2] | ( ex[0] & ex[1] ) ) ;
  
  sim_muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = sim_muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );
  for( size_t i = 0; i < 3; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 3 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 6 );

}

TEST_CASE( "sim muesli create network OneHot ", "[sim_muesli]" )
{
  std::cout << " sim muesli one hot" << std::endl;
  klut_network klut;
  simulation_view klut_sim{ klut };

  std::vector<kitty::partial_truth_table> ex;
  for( size_t i = 0; i < 3; ++i )
  {
    ex.push_back( kitty::partial_truth_table(8u) );
    create_nth_var( ex[i], i );
  }
  kitty::partial_truth_table y(8u);
  y = ( ( ex[0] &( ~ex[1] & ~ ex[2] ) ) ^ ( ~ex[0] &( ex[1] & ~ ex[2] ) ) )^ ( ~ex[0] &( ~ex[1] & ex[2] ) ) ;
  
  sim_muesli_params ps;
  ps.init_sup = 2;
  ps.max_sup = 3;
  ps.max_act = 3;
  ps.eps_th = 0.99;
  ps.verbose = true;
  ps.try_accuracy_recovery = true;

  auto f0 = sim_muesli( klut_sim, ex, y, ps );
  auto ipatterns = klut_sim.get_input_patterns();
  CHECK( ipatterns.size() == 3 );
  for( size_t i = 0; i < 3; ++i )
  {
    CHECK( ipatterns[i].pat == ex[i] );
    CHECK( klut_sim.sim_patterns[2+i].pat == ex[i] );
  }
  CHECK( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ] .pat == y );
  kitty::print_binary( klut_sim.sim_patterns[ klut_sim.get_node_pattern(f0) ].pat ); std::cout << std::endl;
  CHECK( klut_sim.num_gates() == 7 );

  klut_network klut_dec;
  simulation_view klut_dec_sim{ klut_dec };
  sim_decomposition_params decps;
  decps.verbose = true;

  auto f0_dec = sim_decomposition( klut_dec_sim, ex, y, decps );
  klut_dec_sim.create_po( f0_dec );
  CHECK( klut_dec_sim.num_gates() == 6 );
}