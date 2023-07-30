

#include <kitty/constructors.hpp>
#include <kitty/print.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/resyn_engines/xag_resyn.hpp>
#include <mockturtle/algorithms/node_resynthesis/exact.hpp>
#include <mockturtle/algorithms/resyn_engines/xag_resyn.hpp>
#include <mockturtle/algorithms/balancing/sop_balancing.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/write_dot.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/algorithms/sim_resub.hpp>
#include <mockturtle/algorithms/cleanup.hpp>

#include <mockturtle/algorithms/decompose/DecSolver.hpp>
#include <mockturtle/algorithms/ccgame/solvers/cusco.hpp>
#include <mockturtle/algorithms/mcts/mct_tree.hpp>
#include <mockturtle/algorithms/mcts/method.hpp>
#include <mockturtle/algorithms/mcts/nodes/nd_delay.hpp>

#include <string>
#include <bit>
#include <bitset>
#include <cstdint>

#include <iostream>
#include <fstream>
#include <string>
#include <set>

#include <sys/types.h>
#include <sys/stat.h>

using namespace mockturtle;
using namespace ccgame;
using namespace mcts;

template<class Ntk> report_t<Ntk> symm_opt( kitty::dynamic_truth_table *, int, int, std::vector<double> );


int main()
{
  
  using TT = kitty::dynamic_truth_table;
  using Ntk = aig_network;
  
  TT target(4u);
  int id=0;

  std::vector<double> T = {0.0, 0.0, 4.0, 4.0 };
  do
  {
    kitty::next_inplace( target );

    kitty::print_binary( target ); 
    report_t<aig_network> rep1 = symm_opt<aig_network>( &target, 3, 33, T );
    report_t<aig_network> rep2 = symm_opt<aig_network>( &target, 4, 33, T );
    printf("\n");

  } while ( !kitty::is_const0( target ) );

  return 0;
}

template<class Ntk>
report_t<Ntk> symm_opt( kitty::dynamic_truth_table * pF, int MET, int NITERS, std::vector<double> T )
{
  using TT  = kitty::dynamic_truth_table;

  report_t<Ntk> rep;
  typedef DecSolver<TT, Ntk> solver_t;

  kitty::dynamic_truth_table mask = (*pF).construct();
  mask = mask | ~mask;
  solver_t solver( {*pF}, {mask} );
  //solver.PrintSpecs();
  
  if( MET == 0 )
  {
    std::vector<kitty::dynamic_truth_table> xs;
    for( int i{0}; i<pF->num_vars(); ++i )
    {
      xs.emplace_back( pF->num_vars() );
      kitty::create_nth_var( xs[i], i );
    }
    /* define the parameters */
    int nIters = 1;
    cusco_ps ps( ccgame::solver_t::_SYM_1SH, nIters );
    /* solve */
    cusco<Ntk> solver( xs, {*pF} );
    rep = solver.solve( ps );
  }
  else if( MET == 1 )
  {
    //ntk = solver.aut_sym_solve( MET );

    std::vector<kitty::dynamic_truth_table> xs;
    for( int i{0}; i<pF->num_vars(); ++i )
    {
      xs.emplace_back( pF->num_vars() );
      kitty::create_nth_var( xs[i], i );
    }
  
    /* define the parameters */
    int nIters = NITERS;
    cusco_ps ps( ccgame::solver_t::_SYM_RND, nIters );
    /* solve */
    cusco<Ntk> solver( xs, {*pF} );
    rep = solver.solve( ps );
  }
  else if( MET == 2 )
  {
    std::vector<kitty::dynamic_truth_table> xs;
    for( int i{0}; i<pF->num_vars(); ++i )
    {
      xs.emplace_back( pF->num_vars() );
      kitty::create_nth_var( xs[i], i );
    }

    /* define the parameters */    
    int nIters = NITERS;
    cusco_ps ps( ccgame::solver_t::_COV_RND, nIters, -1 );
    /* solve */
    cusco<Ntk> solver( xs, {*pF} );
    rep = solver.solve( ps );
  }
  else if( MET == 3 )
  {
    std::vector<kitty::dynamic_truth_table> xs;
    for( int i{0}; i<pF->num_vars(); ++i )
    {
      xs.emplace_back( pF->num_vars() );
      kitty::create_nth_var( xs[i], i );
    }
    /* define the parameters */
    int nIters = 1;
    cusco_ps ps( ccgame::solver_t::_SYM_1DE, nIters );
    ps.T = {0,0,4,4};
    /* solve */
    cusco<Ntk> solver( xs, {*pF} );
    rep = solver.solve( ps );
  }
  else if( MET == 4 )
  {
    std::vector<kitty::dynamic_truth_table> xs;
    for( int i{0}; i<pF->num_vars(); ++i )
    {
      xs.emplace_back( pF->num_vars() );
      kitty::create_nth_var( xs[i], i );
    }
    /* define the parameters */
    int nIters = 10;
    cusco_ps ps( ccgame::solver_t::_SYM_RDE, nIters );
    ps.T = {0,0,4,4};
    /* solve */
    cusco<Ntk> solver( xs, {*pF} );
    rep = solver.solve( ps );
  }
  else
  {
    assert(0);
  }

  if( rep.Esl )
  {
    default_simulator<kitty::dynamic_truth_table> sim( (*pF).num_vars() );
    const auto tt = simulate<kitty::dynamic_truth_table>( rep.ntk, sim )[0];
    assert( kitty::equal( tt, *pF ) );    
      printf(" [symm %d %d] ", rep.nMin, rep.levels );
    return rep;
  }
  else
    printf(" [symm X X ] " );
  return rep;
}

//template<class Ntk>
//Ntk mcts_opt( kitty::dynamic_truth_table * pF, int MET, int NITERS, std::vector<double> T )
//{
//  using TT  = kitty::dynamic_truth_table;
//
//  std::vector<kitty::dynamic_truth_table> X;
//
//  node_ps ndps;
//  detailed_gate_t cmpr_( mcts::gate_t::CMPR, 1, 1.0, 1.0, &hpcompute_cmpr );
//  detailed_gate_t cmpl_( mcts::gate_t::CMPR, 1, 1.0, 1.0, &hpcompute_cmpl );
//  detailed_gate_t ai00_( mcts::gate_t::AI00, 2, 1.0, 1.0, &hpcompute_ai00 );
//  detailed_gate_t ai11_( mcts::gate_t::AI11, 2, 1.0, 1.0, &hpcompute_ai11 );
//  detailed_gate_t exor_( mcts::gate_t::EXOR, 2, 1.0, 1.0, &hpcompute_exor );
//  ndps.lib = { cmpl_, cmpr_, ai00_, ai11_ };
//
//  mct_ps mctps;
//  ndps.sel_type = supp_selection_t::SUP_ENER;
//  mctps.nIters =300;
//  mctps.nSims = 1;
//  mctps.verbose =false;
//  ndps.BETA0 = 100;
//  ndps.BETAZ = 100;
//  ndps.nIters = 1;
//
//  for( int i{0}; i<pF->num_vars(); ++i )
//  {
//    X.emplace_back( pF->num_vars() );
//    kitty::create_nth_var( X[i], i );
//  }
//
//  nd_delay_t<Ntk> root( X, T, {* pF}, ndps );
//  mct_method_ps metps;
//
//  mct_method_t<nd_delay_t<Ntk> > meth( metps );
//  mct_tree_t<nd_delay_t<Ntk> , mct_method_t> mct( root, meth, mctps );
//  int iSol = mct.solve();
//
//  if( iSol == -1 )
//  {
//    printf(" [mcts X X ] " );
//    Ntk xag;
//    return xag;
//  }
//  else
//  {
//      Ntk xag = cleanup_dangling( mct.nodes[iSol].ntk );
//      printf( "[mcts %d %f] ", xag.num_gates(), mct.evaluate(iSol) );
//      default_simulator<kitty::dynamic_truth_table> sim( 4u );
//      const auto tt = simulate<kitty::dynamic_truth_table>( xag, sim )[0];
//      assert( kitty::equal( tt, * pF ) );
//      return xag;
//  }
//
//}