

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

template<class Ntk> report_t<Ntk> game_on( kitty::dynamic_truth_table *, int, int, std::vector<double> );


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
    report_t<aig_network> rep = game_on<aig_network>( &target, 3, 33, T );
    if(!rep.Esl)  continue;
    kitty::print_binary( target ); 
    printf(" -> %d %d\n", rep.nMin, rep.levels );

  } while ( !kitty::is_const0( target ) );

  return 0;
}

template<class Ntk>
report_t<Ntk> game_on( kitty::dynamic_truth_table * pF, int MET, int NITERS, std::vector<double> T )
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
  else
  {
    assert(0);
  }

  if( rep.Esl )
  {
    default_simulator<kitty::dynamic_truth_table> sim( (*pF).num_vars() );
    const auto tt = simulate<kitty::dynamic_truth_table>( rep.ntk, sim )[0];
    assert( kitty::equal( tt, *pF ) );    
    return rep;
  }
  else
    printf("NO SOL FOUND\n");
  return rep;
}