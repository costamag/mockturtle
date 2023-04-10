#include <catch.hpp>

#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/algorithms/decompose/DecNet.hpp>
#include <mockturtle/algorithms/decompose/DecAnalyzer.hpp>
#include <mockturtle/algorithms/decompose/DecChsToGraph.hpp>
#include <mockturtle/algorithms/decompose/DecSolver.hpp>
#include <mockturtle/networks/aig.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <kitty/operations.hpp>
#include <kitty/print.hpp>

using namespace mockturtle;

TEST_CASE( "Simulations storage", "[DEC]" )
{
  using TT  = kitty::dynamic_truth_table;
  using Ntk = aig_network;

  std::vector<TT> vFuncs;
  std::vector<TT> vMasks;
  for( int i=0; i<3; ++i )
  {
    vFuncs.emplace_back(3u);
    vMasks.emplace_back(3u);
    kitty::create_random( vFuncs[i] );
    kitty::create_random( vMasks[i] );
  }

  /* create a targets data struct */
  typedef DecSims<TT> targets_t;
  targets_t tars;
  CHECK( tars.size() == 0 );
  sim_t t0 = tars.addSim( vFuncs[0], vMasks[0] );
  CHECK( tars.isUsed(0) );
  CHECK( t0 == 0 );
  CHECK( tars.size() == 1 );
  CHECK( * tars.getFuncP(0) == vFuncs[0] );
  CHECK( * tars.getMaskP(0) == vMasks[0] );
  sim_t t1 = tars.addSim( vFuncs[1], vMasks[1] );
  CHECK( tars.isUsed(1) );
  CHECK( t1 == 1 );
  CHECK( tars.size() == 2 );
  CHECK( * tars.getFuncP(1) == vFuncs[1] );
  CHECK( * tars.getMaskP(1) == vMasks[1] );
  tars.remove(0);
  CHECK( !tars.isUsed(0) );
  CHECK( tars.size() == 1 );
  CHECK( kitty::is_const0(* tars.getFuncP(0)) );
  CHECK( kitty::is_const0( ~(* tars.getMaskP(0))) );
  sim_t t2 = tars.addSim( vFuncs[2], vMasks[2] );
  CHECK( tars.isUsed(0) );
  CHECK( t2 == 0 );
  CHECK( * tars.getFuncP(0) == vFuncs[2] );
  CHECK( * tars.getMaskP(0) == vMasks[2] );
}

TEST_CASE( "linking nodes to the simulation storage", "[DEC]" )
{
  using TT = kitty::dynamic_truth_table;
  using Ntk = aig_network;

  std::vector<TT> vFuncs;
  std::vector<TT> vMasks;
  for( int i=0; i<3; ++i )
  {
    vFuncs.emplace_back(3u);
    vMasks.emplace_back(3u);
    kitty::create_random( vFuncs[i] );
    vMasks[i] |= ~vMasks[i];
  }
  typedef DecSims<TT>    targets_t;
  typedef DecNodes<Ntk>  graph_t;

  targets_t tars;
  graph_t   graph;

  std::vector<sim_t> t;
  std::vector<node_t> n;
  for( int i = 0; i<3; ++i )
  {
    t.push_back( tars.addSim( vFuncs[i], vMasks[i] ) );
    n.push_back( graph.addNode( {}, t[i], DecFunc_t::PI_ ) );
  }

  t.push_back(tars.addSim( vFuncs[0] & vFuncs[1], vMasks[0] ));
  n.push_back(graph.addNode( { n[0], n[1] }, t[3], DecFunc_t::AND_ ));
  std::vector<node_t> fanins = {n[0], n[1]};
  CHECK( * graph.getFanInsP(n[3]) == fanins );
  CHECK( graph.getFunc(n[3]) == DecFunc_t::AND_ );
  CHECK( graph.getSim(n[3]) == t[3] );
}
/*
TEST_CASE( "the choicesim network", "[DEC]" )
{
  using TT = kitty::dynamic_truth_table;
  using Ntk = aig_network;

  typedef DecNet<TT, Ntk> net_t;
  net_t net;
  std::vector<signal_t> pis;
  std::vector<TT> xs;
  for( int i{0}; i<3; ++i )
  {
    xs.emplace_back(3u);
    kitty::create_nth_var( xs[i], i );
    pis.push_back( net.create_PI(xs[i]) );
  }
  CHECK( net.numPIs() == 3 );
 for( int i{0}; i<3; ++i )
  {
    CHECK( kitty::equal( xs[i], * net.getFuncP( pis[i] ) ) );
    CHECK( kitty::is_const0( ~ * net.getMaskP( pis[i] ) ) ); 
  }
  typedef DecAnalyzer<TT> analize_t;
  analize_t ana;
  std::vector<bool> isDec = { false, false, false };
  // top and decomposability 
  signal_t x0 = pis[0];
  signal_t x1 = pis[1];
  signal_t x2 = pis[2];
  TT f1 = kitty::dynamic_truth_table(3u);
  f1 = xs[0] & ( xs[1]^xs[2] );
  std::vector<signal_t> targets;
  TT m = f1 | ~f1;
  targets.push_back( net.create_target( f1, m ) );
  signal_t x3 = net.create_xor( x1, x2 );
  signal_t x4 = net.create_and( x0, x3 );
  net.create_PO( x4 );
  for( int i = 0; i < 3; ++i )
    isDec[i] = ana.IsTopAndDec( net.getFuncP( x4 ), net.getMaskP( x4 ), net.getFuncP( pis[i] ), net.getMaskP( pis[i] ) );
  CHECK( isDec[0]  );
  CHECK( !isDec[1] );
  CHECK( !isDec[2] );
  for( int i = 0; i < 3; ++i )
    isDec[i] = ana.IsTopAndDec( net.getFuncP( targets[0] ), net.getMaskP( targets[0] ), net.getFuncP( pis[i] ), net.getMaskP( pis[i] ) );
  CHECK( isDec[0]  );
  CHECK( !isDec[1] );
  CHECK( !isDec[2] );
  CHECK( net.numPOs() == 1 );
  CHECK( kitty::equal( f1, * net.getFuncP( targets[0] ) ) ); 
  CHECK( kitty::equal( m, * net.getMaskP( targets[0] ) ) ); 
}
*/

TEST_CASE( "converting choicesim network", "[DEC]" )
{
  using TT = kitty::dynamic_truth_table;
  using Ntk = aig_network;
  typedef DecNet<TT, Ntk> net_t;
  typedef DecChsToGraph<TT, Ntk> cnv_t;
  net_t net;
  Ntk aig;

  std::vector<signal_t> xs;
  std::vector<TT> tts;
  for( int i{0}; i<4; ++i )
  {
    tts.emplace_back(4u);
    kitty::create_nth_var( tts[i], i );
    xs.push_back( net.create_PI(tts[i]) );
  }
  signal_t x4 = net.create_xor( xs[1], xs[2] );
  signal_t x5 = net.create_and( xs[0], xs[3] );
  signal_t x6 = net.create_or( x4, x5 );
  signal_t x7 = net.create_lt( x6, x5 );
  signal_t x8 = net.create_le( xs[0], x7 );
  signal_t x9 = net.create_ge( xs[1], x8 );
  signal_t x10 = net.create_gt( xs[2], x9 );
  net.create_PO( x10 );

  cnv_t conv( net );  
  aig = conv.convert();
  
  CHECK( aig.num_pis() == 4 );
  CHECK( aig.num_pos() == 1 );
  CHECK( aig.num_gates() == 9 );
}

TEST_CASE( "checking the simulation patterns in choicesim 0", "[DEC]" )
{
  using TT = kitty::dynamic_truth_table;
  using Ntk = aig_network;
  typedef DecNet<TT, Ntk> net_t;
  typedef DecChsToGraph<TT, Ntk> cnv_t;
  net_t net;
  Ntk aig;

  std::vector<signal_t> xs;
  std::vector<TT> tts;
  for( int i{0}; i<2; ++i )
  {
    tts.emplace_back(2u);
    kitty::create_nth_var( tts[i], i );
    xs.push_back( net.create_PI(tts[i]) );
  }
  std::vector<signal_t> vSigs;
  vSigs.push_back(net.create_xor( xs[0], xs[1] ));
  vSigs.push_back(net.create_xnor( xs[0], xs[1] ));
  vSigs.push_back(net.create_and( xs[0], xs[1] ));
  vSigs.push_back(net.create_nand( xs[0], xs[1] ));
  vSigs.push_back(net.create_or( xs[0], xs[1] ));
  vSigs.push_back(net.create_nor( xs[0], xs[1] ));
  vSigs.push_back(net.create_le( xs[0], xs[1] ));
  vSigs.push_back(net.create_gt( xs[0], xs[1] ));
  vSigs.push_back(net.create_lt( xs[0], xs[1] ));
  vSigs.push_back(net.create_ge( xs[0], xs[1] ));
  vSigs.push_back(net.create_not( xs[0] ));
  vSigs.push_back(net.create_buf( xs[1] ));
  for( auto x : vSigs )
    net.create_PO( x );

  cnv_t conv( net );  
  aig = conv.convert();
  default_simulator<kitty::dynamic_truth_table> sim( 2 );
  const auto sims = simulate<kitty::dynamic_truth_table>( aig, sim );
  for( uint32_t i=0; i<vSigs.size(); ++i )
    CHECK( sims[i] == * net.getFuncP( vSigs[i] ) );
}

TEST_CASE( "checking the simulationn patterns in choicesim 1", "[DEC]" )
{
  using TT = kitty::dynamic_truth_table;
  using Ntk = aig_network;
  typedef DecNet<TT, Ntk> net_t;
  typedef DecChsToGraph<TT, Ntk> cnv_t;
  net_t net;
  Ntk aig;

  std::vector<signal_t> xs;
  std::vector<TT> tts;
  for( int i{0}; i<4; ++i )
  {
    tts.emplace_back(4u);
    kitty::create_nth_var( tts[i], i );
    xs.push_back( net.create_PI(tts[i]) );
  }
  signal_t x4 = net.create_xor( xs[1], xs[2] );
  signal_t x5 = net.create_and( xs[0], xs[3] );
  signal_t x6 = net.create_or( x4, x5 );
  signal_t x7 = net.create_lt( x6, x5 );
  signal_t x8 = net.create_le( xs[0], x7 );
  signal_t x9 = net.create_ge( xs[1], x8 );
  signal_t x10 = net.create_gt( xs[2], x9 );
  net.create_PO( x10 );
  cnv_t conv( net );  
  aig = conv.convert();
  default_simulator<kitty::dynamic_truth_table> sim( 4 );
  const auto sims = simulate<kitty::dynamic_truth_table>( aig, sim );
  CHECK( sims[0] == * net.getFuncP( x10 ) );
}

TEST_CASE( "checking the simulation patterns in choicesim 2", "[DEC]" )
{
  using TT = kitty::dynamic_truth_table;
  using Ntk = aig_network;
  typedef DecNet<TT, Ntk> net_t;
  typedef DecChsToGraph<TT, Ntk> cnv_t;
  net_t net;
  Ntk aig;

  std::vector<signal_t> xs;
  std::vector<TT> tts;
  for( int i{0}; i<4; ++i )
  {
    tts.emplace_back(4u);
    kitty::create_nth_var( tts[i], i );
    xs.push_back( net.create_PI(tts[i]) );
  }
  signal_t x4 = net.create_xnor( xs[1], xs[2] );
  signal_t x5 = net.create_nand( xs[0], xs[3] );
  signal_t x6 = net.create_nor( x4, x5 );
  signal_t x7 = net.create_le( x6, x5 );
  signal_t x8 = net.create_gt( xs[0], x7 );
  signal_t x9 = net.create_not( x8 );
  signal_t x10 = net.create_buf( x9 );
  net.create_PO( x10 );
  cnv_t conv( net );  
  aig = conv.convert();
  default_simulator<kitty::dynamic_truth_table> sim( 4 );
  const auto sims = simulate<kitty::dynamic_truth_table>( aig, sim );
  CHECK( sims[0] == * net.getFuncP( x10 ) );
}

TEST_CASE( "solver: top decompositions ", "[DEC]" )
{
  using TT = kitty::dynamic_truth_table;
  using Ntk = aig_network;
  typedef DecSolver<TT, Ntk> solver_t;
  Ntk aig;

  std::vector<TT> xs;
  for( int i{0}; i<5; ++i )
  {
    xs.emplace_back(5u);
    kitty::create_nth_var( xs[i], i );
  }

  std::vector<TT> vTruths;
  std::vector<TT> vMasks;
  vTruths.emplace_back(5u);
  vTruths.emplace_back(5u);
  vTruths.emplace_back(5u);
  vTruths[0] = xs[4] & ( ~xs[3] & ( xs[2]^ xs[0] ) );
  vTruths[1] = (xs[4] & xs[2]) ^ xs[2];
  vTruths[2] = xs[3] & xs[2] &xs[1] ;
  
  vMasks.emplace_back(5u);
  vMasks.emplace_back(5u);
  vMasks.emplace_back(5u);
  vMasks[0] |= ~vMasks[0];
  vMasks[1] = vMasks[0];
  vMasks[2] = vMasks[0];
  solver_t solver( vTruths, vMasks );
  solver.PrintSpecs();
  solver.solve();
}