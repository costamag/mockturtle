#include <catch.hpp>

#include <mockturtle/algorithms/aig_algebraic_rewriting.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/algorithms/equivalence_checking.hpp>
#include <mockturtle/algorithms/miter.hpp>
#include <mockturtle/algorithms/functional_reduction.hpp>
#include <kitty/static_truth_table.hpp>
#include <lorina/aiger.hpp>

using namespace mockturtle;

TEST_CASE( "Simple associativity (AND)", "[aig_algebraic_rewriting]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_and( pis[0], pis[1] );
  const auto f2 = aig.create_and( f1, pis[2] );
  const auto f3 = aig.create_and( f2, pis[3] );
  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};
  CHECK( depth_aig.depth() == 2 );

  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

TEST_CASE( "Simple associativity (OR)", "[aig_algebraic_rewriting]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_or( pis[0], pis[1] );
  const auto f2 = aig.create_or( f1, pis[2] );
  const auto f3 = aig.create_or( f2, pis[3] );
  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};
  CHECK( depth_aig.depth() == 2 );

  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

TEST_CASE( "Simple distributivity (OR on top)", "[aig_algebraic_rewriting]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto g = aig.create_xor( pis[0], pis[1] );
  const auto f1 = aig.create_and( g, pis[2] );
  const auto f2 = aig.create_and( g, pis[3] );
  const auto f3 = aig.create_or( f1, f2 );
  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};
  CHECK( depth_aig.depth() == 3 );

  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

TEST_CASE( "Simple distributivity (AND on top)", "[aig_algebraic_rewriting]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto g = aig.create_xor( pis[0], pis[1] );
  const auto f1 = aig.create_or( g, pis[2] );
  const auto f2 = aig.create_or( g, pis[3] );
  const auto f3 = aig.create_and( f1, f2 );
  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};
  CHECK( depth_aig.depth() == 3 );

  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

TEST_CASE( "Three-layer distributivity", "[aig_algebraic_rewriting]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{5};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  /* This rule is not in the pdf, but also simple: 
     ((g x2) + x3 ) x4 = (g x2 x4) + (x3 x4) = (g (x2 x4)) + (x3 x4) */
  const auto g = aig.create_xor( pis[0], pis[1] );
  const auto f1 = aig.create_and( g, pis[2] );
  const auto f2 = aig.create_or( f1, pis[3] );
  const auto f3 = aig.create_and( f2, pis[4] );
  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};
  CHECK( depth_aig.depth() == 4 );

  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

TEST_CASE( "Depth optimization on ISCAS benchmarks", "[aig_algebraic_rewriting]" )
{
  uint32_t benchmark_ids[11] = {17, 432, 499, 880, 1355, 1908, 2670, 3540, 5315, 6288, 7552};
  uint32_t expected_depths[11] = {3, 26, 19, 19, 25, 26, 18, 35, 34, 120, 25};

  for ( uint32_t i = 0u; i < 11; ++i )
  {
    aig_network ntk, ntk_ori;
    auto const result = lorina::read_aiger( fmt::format( "{}/c{}.aig", BENCHMARKS_PATH, benchmark_ids[i] ), aiger_reader( ntk ) );
    if ( result != lorina::return_code::success )
    {
      continue;
    }
    ntk_ori = cleanup_dangling( ntk );

    /* call the algorithm */
    aig_algebraic_rewriting( ntk );

    /* check the resulting depth */
    /* (You should already pass by implementing the rules introduced in the pdf,
        but if you have implemented more rules, better results are possible.) */
    depth_view depth_aig{ntk};
    fmt::print( "[i] On benchmark c{}.aig: Optimized depth = {} (expected at most {})\n", 
                benchmark_ids[i], depth_aig.depth(), expected_depths[i] );
    CHECK( depth_aig.depth() <= expected_depths[i] );

    /* equivalence checking */
    aig_network miter_aig = *miter<aig_network>( ntk_ori, ntk );
    functional_reduction( miter_aig );
    bool cec = *equivalence_checking( miter_aig );
    CHECK( cec == true );
  }
}

/* #############################################################
 *  Two levels two nodes
 * #############################################################
*/
/* ############## S2 ##################*/
// S2 phi1 = phi3
TEST_CASE( "Two levels two nodes - S2 phi2(s) = not(s) phi1=phi3  - 1/2 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{2};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_and( pis[0], pis[1] ); 

  const auto f2 = aig.create_and( !f1, pis[0] );

  aig.create_po( f2 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 1 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

TEST_CASE( "Two levels two nodes - S2 phi2(s) = not(s) phi1=phi4  - 2/2 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{2};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_and( pis[0], pis[1] ); 

  const auto f2 = aig.create_and( !f1, pis[1] );

  aig.create_po( f2 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 1 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

TEST_CASE( "Two levels two nodes - S2 phi2(s) = not(s) phi1=not(phi3)  - 1/2 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{2};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_and( pis[0], pis[1] ); 

  const auto f2 = aig.create_and( !f1, !pis[0] );

  aig.create_po( f2 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 0 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

TEST_CASE( "Two levels two nodes - S2 phi2(s) = not(s) phi1=not(phi4)  - 2/2 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{2};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_and( pis[0], pis[1] ); 

  const auto f2 = aig.create_and( !f1, !pis[1] );

  aig.create_po( f2 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 0 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

/* ############## S3 ##################*/
// S3 phi1 = phi3
TEST_CASE( "Two levels two nodes - S3 phi2(s) = s phi1=phi3  - 1/2 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{2};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_and( pis[0], pis[1] ); 

  const auto f2 = aig.create_and( f1, pis[0] );

  aig.create_po( f2 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 1 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

TEST_CASE( "Two levels two nodes - S3 phi2(s) = s phi1=phi4  - 2/2 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{2};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_and( pis[0], pis[1] ); 

  const auto f2 = aig.create_and( f1, pis[1] );

  aig.create_po( f2 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 1 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

TEST_CASE( "Two levels two nodes - S3 phi2(s) = s phi1=not(phi3)  - 1/2 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{2};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_and( pis[0], pis[1] ); 

  const auto f2 = aig.create_and( f1, !pis[0] );

  aig.create_po( f2 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 0 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

TEST_CASE( "Two levels two nodes - S3 phi2(s) = s phi1=not(phi4)  - 2/2 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{2};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_and( pis[0], pis[1] ); 

  const auto f2 = aig.create_and( f1, !pis[1] );

  aig.create_po( f2 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 0 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}
// s3b has already been tested in aig_algebraic_rewriting
/* #############################################################
 *  Two levels three nodes
 * #############################################################
*/
// S1a
TEST_CASE( "Two levels two nodes - S1a phi2(s) = phi1(s) = s phi4=phi5  - 1/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f0 = aig.create_or( pis[0], pis[1] );
  const auto f1 = aig.create_and( pis[2], f0 ); 
  const auto f2 = aig.create_and( f0, pis[3] );
  const auto f3 = aig.create_and( f1, f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 2 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

TEST_CASE( "Two levels two nodes - S1a phi2(s) = phi1(s) = s phi3=phi5  - 2/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f0 = aig.create_or( pis[0], pis[1] );
  const auto f1 = aig.create_and( f0, pis[2] ); 
  const auto f2 = aig.create_and( f0, pis[3] );
  const auto f3 = aig.create_and( f1, f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 2 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

TEST_CASE( "Two levels two nodes - S1a phi2(s) = phi1(s) = s phi3=phi6  - 3/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f0 = aig.create_or( pis[0], pis[1] );
  const auto f1 = aig.create_and( f0, pis[2] ); 
  const auto f2 = aig.create_and( pis[3], f0 );
  const auto f3 = aig.create_and( f1, f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 2 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

TEST_CASE( "Two levels two nodes - S1a phi2(s) = phi1(s) = s phi4=phi6  - 4/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f0 = aig.create_or( pis[0], pis[1] );
  const auto f1 = aig.create_and( pis[2], f0 ); 
  const auto f2 = aig.create_and( pis[3], f0 );
  const auto f3 = aig.create_and( f1, f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 2 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

// S1b
TEST_CASE( "Two levels two nodes - S1b phi2(s) = phi1(s) = s phi4=not(phi5)  - 1/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f0 = aig.create_or( pis[0], pis[1] );
  const auto f1 = aig.create_and( pis[2], f0 ); 
  const auto f2 = aig.create_and( !f0, pis[3] );
  const auto f3 = aig.create_and( f1, f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 0 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

TEST_CASE( "Two levels two nodes - S1b phi2(s) = phi1(s) = s phi3=not(phi5)  - 2/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f0 = aig.create_or( pis[0], pis[1] );
  const auto f1 = aig.create_and( f0, pis[2] ); 
  const auto f2 = aig.create_and( !f0, pis[3] );
  const auto f3 = aig.create_and( f1, f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 0 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

TEST_CASE( "Two levels two nodes - S1b phi2(s) = phi1(s) = s phi3=not(phi6)  - 3/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f0 = aig.create_or( pis[0], pis[1] );
  const auto f1 = aig.create_and( f0, pis[2] ); 
  const auto f2 = aig.create_and( pis[3], !f0 );
  const auto f3 = aig.create_and( f1, f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 0 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

TEST_CASE( "Two levels two nodes - S1b phi2(s) = phi1(s) = s phi4=not(phi6)  - 4/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f0 = aig.create_or( pis[0], pis[1] );
  const auto f1 = aig.create_and( pis[2], f0 ); 
  const auto f2 = aig.create_and( pis[3], !f0 );
  const auto f3 = aig.create_and( f1, f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 0 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}


// S2a and S2b
TEST_CASE( "Two levels two nodes - S2a phi2(s) = phi1(s) = not(s) phi4=phi5 phi3=phi6'  - 1/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{2};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_and( pis[0], pis[1] ); 
  const auto f2 = aig.create_and( pis[1], !pis[0] );
  const auto f3 = aig.create_and( !f1, !f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 0 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

// S2a and S2b
TEST_CASE( "Two levels two nodes - S2a phi2(s) = phi1(s) = not(s) phi4=phi6 phi3=phi5'  - 2/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{2};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_and( pis[0], pis[1] ); 
  const auto f2 = aig.create_and( !pis[0], pis[1] );
  const auto f3 = aig.create_and( !f1, !f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 0 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

// S2a and S2b
TEST_CASE( "Two levels two nodes - S2a phi2(s) = phi1(s) = not(s) phi4=phi5' phi3=phi6  - 3/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{2};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_and( pis[0], pis[1] ); 
  const auto f2 = aig.create_and( !pis[1], pis[0] );
  const auto f3 = aig.create_and( !f1, !f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 0 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

// S2a and S2b
TEST_CASE( "Two levels two nodes - S2a phi2(s) = phi1(s) = not(s) phi4=phi6' phi3=phi5  - 4/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{2};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_and( pis[0], pis[1] ); 
  const auto f2 = aig.create_and( pis[0], !pis[1] );
  const auto f3 = aig.create_and( !f1, !f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 0 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

// S2c
TEST_CASE( "Two levels two nodes - S2c phi2(s) = phi1(s) = not(s) phi4=phi5  - 1/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f0 = aig.create_or( pis[0], pis[1] );
  const auto f1 = aig.create_and( pis[2], f0 ); 
  const auto f2 = aig.create_and( f0, pis[3] );
  const auto f3 = aig.create_and( !f1, !f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 2 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

TEST_CASE( "Two levels two nodes - S2c phi2(s) = phi1(s) = not(s) phi4=phi6  - 2/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f0 = aig.create_or( pis[0], pis[1] );
  const auto f1 = aig.create_and( pis[2], f0 ); 
  const auto f2 = aig.create_and( pis[3], f0 );
  const auto f3 = aig.create_and( !f1, !f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 2 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

TEST_CASE( "Two levels two nodes - S2c phi2(s) = phi1(s) = not(s) phi3=phi5  - 3/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f0 = aig.create_or( pis[0], pis[1] );
  const auto f1 = aig.create_and( f0, pis[2] ); 
  const auto f2 = aig.create_and( f0, pis[3] );
  const auto f3 = aig.create_and( !f1, !f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 2 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

TEST_CASE( "Two levels two nodes - S2c phi2(s) = phi1(s) = not(s) phi3=phi6  - 4/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f0 = aig.create_or( pis[0], pis[1] );
  const auto f1 = aig.create_and( f0, pis[2] ); 
  const auto f2 = aig.create_and( pis[3], f0 );
  const auto f3 = aig.create_and( !f1, !f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 2 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

// S3a
TEST_CASE( "Two levels two nodes - S3a phi1(s) = phi2(s)' = s' phi4=phi5  - 1/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f0 = aig.create_or( pis[0], pis[1] );
  const auto f1 = aig.create_and( pis[2], f0 ); 
  const auto f2 = aig.create_and( f0, pis[3] );
  const auto f3 = aig.create_and( !f1, f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 2 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

// S3a
TEST_CASE( "Two levels two nodes - S3a phi1(s) = phi2(s)' = s' phi4=phi6  - 2/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f0 = aig.create_or( pis[0], pis[1] );
  const auto f1 = aig.create_and( pis[2], f0 ); 
  const auto f2 = aig.create_and( pis[3], f0 );
  const auto f3 = aig.create_and( !f1, f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 2 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

// S3a
TEST_CASE( "Two levels two nodes - S3a phi1(s) = phi2(s)' = s' phi3=phi5  - 3/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f0 = aig.create_or( pis[0], pis[1] );
  const auto f1 = aig.create_and( f0, pis[2] ); 
  const auto f2 = aig.create_and( f0, pis[3] );
  const auto f3 = aig.create_and( !f1, f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 2 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

// S3a
TEST_CASE( "Two levels two nodes - S3a phi1(s) = phi2(s)' = s' phi3=phi6  - 4/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f0 = aig.create_or( pis[0], pis[1] );
  const auto f1 = aig.create_and( f0, pis[2] ); 
  const auto f2 = aig.create_and( pis[3], f0 );
  const auto f3 = aig.create_and( !f1, f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 2 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

// S3a
TEST_CASE( "Two levels two nodes - S3a phi1(s) = phi2(s)' = s phi4=phi5  - 1/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f0 = aig.create_or( pis[0], pis[1] );
  const auto f1 = aig.create_and( pis[2], f0 ); 
  const auto f2 = aig.create_and( f0, pis[3] );
  const auto f3 = aig.create_and( f1, !f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 2 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

// S3a
TEST_CASE( "Two levels two nodes - S3a phi1(s) = phi2(s)' = s phi4=phi6  - 2/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f0 = aig.create_or( pis[0], pis[1] );
  const auto f1 = aig.create_and( pis[2], f0 ); 
  const auto f2 = aig.create_and( pis[3], f0 );
  const auto f3 = aig.create_and( f1, !f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 2 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

// S3a
TEST_CASE( "Two levels two nodes - S3a phi1(s) = phi2(s)' = s phi3=phi5  - 3/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f0 = aig.create_or( pis[0], pis[1] );
  const auto f1 = aig.create_and( f0, pis[2] ); 
  const auto f2 = aig.create_and( f0, pis[3] );
  const auto f3 = aig.create_and( f1, !f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 2 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

// S3a
TEST_CASE( "Two levels two nodes - S3a phi1(s) = phi2(s)' = s phi3=phi6  - 4/4 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f0 = aig.create_or( pis[0], pis[1] );
  const auto f1 = aig.create_and( f0, pis[2] ); 
  const auto f2 = aig.create_and( pis[3], f0 );
  const auto f3 = aig.create_and( f1, !f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 2 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

TEST_CASE( "Two levels two nodes - S3b phi1(s) = phi2(s)' = s' phi4=phi5  - 1/1 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{3};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_and( pis[0], pis[1] ); 
  const auto f2 = aig.create_and( !pis[1], pis[2] );
  const auto f3 = aig.create_and( !f1, f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 1 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

/* #############################################################
 *  Three levels three nodes
 * #############################################################
*/

/*      S1a1      */

TEST_CASE( "Three levels three nodes - S1a1 phi2(s) = phi4(s) = s' phi4=phi5 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_and( pis[0], pis[1] ); 
  const auto f2 = aig.create_and( !f1 , pis[2] );
  const auto f3 = aig.create_and( pis[0], !f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 2 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

/*      S1a2      */

TEST_CASE( "Three levels three nodes - S1a1 phi2(s) = phi4(s) = s' phi4=phi5' ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_and( pis[0], pis[1] ); 
  const auto f2 = aig.create_and( !f1 , pis[2] );
  const auto f3 = aig.create_and( !pis[0], !f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 1 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}


/*      S1b1     */

TEST_CASE( "Three levels three nodes - S1b1 phi2(s) = s' phi4(s) = s phi1=phi5 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_and( pis[0], pis[1] ); 
  const auto f2 = aig.create_and( f1 , pis[2] );
  const auto f3 = aig.create_and( pis[0], !f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 2 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

/*      S1b2     */

TEST_CASE( "Three levels three nodes - S1b2 phi2(s) = s' phi4(s) = s phi1=phi5' ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_and( pis[0], pis[1] ); 
  const auto f2 = aig.create_and( f1 , pis[2] );
  const auto f3 = aig.create_and( !pis[0], !f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 0 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

/*      S2a1     */

TEST_CASE( "Three levels three nodes - S2a1 phi2(s) = s phi4(s) = s' phi1=phi5 ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_and( pis[0], pis[1] ); 
  const auto f2 = aig.create_and( !f1 , pis[2] );
  const auto f3 = aig.create_and( pis[0], f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 2 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

/*      S2a2     */

TEST_CASE( "Three levels three nodes - S2a2 phi2(s) = s phi4(s) = s' phi1=phi5' ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_and( pis[0], pis[1] ); 
  const auto f2 = aig.create_and( !f1 , pis[2] );
  const auto f3 = aig.create_and( !pis[0], f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 1 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}


/*      S2a2     */

TEST_CASE( "Three levels three nodes - S2b phi2(s) = phi4(s) = s phi1=phi5' ", "[aig_algebraic_rewriting2]" )
{
  /* create the network */
  aig_network aig;
  static const uint32_t num_pis{4};
  std::vector<typename aig_network::signal> pis;
  for ( uint32_t i = 0; i < num_pis; ++i )
    pis.emplace_back( aig.create_pi() );

  const auto f1 = aig.create_and( pis[0], pis[1] ); 
  const auto f2 = aig.create_and( !f1 , pis[2] );
  const auto f3 = aig.create_and( !pis[0], f2 );

  aig.create_po( f3 );

  /* simulate to get the output truth table(s) */
  auto tts = simulate<kitty::static_truth_table<num_pis>>( aig );

  /* call the algorithm */
  aig_algebraic_rewriting( aig );

  /* check the resulting depth */
  depth_view depth_aig{aig};

  CHECK( depth_aig.depth() == 1 );
 
  /* check that the output functions remain the same */
  CHECK( tts == simulate<kitty::static_truth_table<num_pis>>( aig ) );
}

