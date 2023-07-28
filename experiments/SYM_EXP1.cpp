#pragma once

#include <kitty/constructors.hpp>
#include <kitty/print.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/networks/aig.hpp>
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

#include <mockturtle/algorithms/ccgame/solvers/cusco.hpp>

using namespace mockturtle;
using namespace ccgame;

template<typename Ntk>
report_t<Ntk> sym_solve( kitty::dynamic_truth_table *, library__t );

uint32_t tt_to_key( DTT tt )
{
  uint32_t uint_tt;
  std::string string_tt = kitty::to_hex(tt);
  sscanf( string_tt.c_str(), "%x", &uint_tt ); 
  return uint_tt & 0xFFFF;
}

DTT key_to_tt( uint32_t key )
{
  key = key & 0xFFFF;
  std::string bstring = "";
  for( int iBit{0}; iBit < 16u; ++iBit )
    bstring = ( ((key >> iBit) & 1u) == 1u) ? "1"+bstring : "0"+bstring;
  DTT res(4u);
  kitty::create_from_binary_string(res,bstring);
  return res;
}


int main()
{
  using TT = kitty::dynamic_truth_table;
  using NTK = aig_network;
  TT target(4u);

  library__t lib;
  lib.AI00 = {1., 1.  , 1.0 };
  lib.AI01 = {1.5, 1.0, 2.0 };
  lib.AI10 = {1.0, 1.5, 2.0 };
  lib.AI11 = {1.5, 1.5, 1.0 };
  lib.CMPL = {0.5, 0.5, 1.0 };
  lib.CMPR = {0.5, 0.5, 1.0 };
  lib.CNTR = {0., 0.  , 0.0 };
  lib.EXOR = {2., 2.  , 1.0 };
  lib.OI00 = {1.5, 1.5, 2.0 };
  lib.OI01 = {2.0, 1.5, 2.0 };
  lib.OI10 = {1.5, 2.0, 2.0 };
  lib.OI11 = {2., 2.  , 1.0 };
  lib.PIS  = {0., 0.  , 0.0 };
  lib.POS  = {0., 0.  , 0.0 };
  lib.PRJL = {0., 0.  , 0.0 };
  lib.PRJR = {0., 0.  , 0.0 };
  lib.TAUT = {0., 0.  , 0.0 };
  lib.XNOR = {2.5, 2.5, 2.0 };

  std::vector<double> DELS;
  std::vector<bool> USED;
  std::vector<double> SIZS;

  for( int i{0}; i<pow(2,pow(2,4)); ++i )
  {
    DELS.push_back(0.0);
    USED.emplace_back(false);
    SIZS.emplace_back(0);
  }

  do
  {
    uint32_t key = tt_to_key( target );

    //target = key_to_tt(10300);
    //printf("%d\n", tt_to_key(target));
    //assert( 10300 == tt_to_key(target) );

    printf("FUNC %d\n", key );

    report_t<NTK> rep = sym_solve<NTK>( &target, lib );
    //kitty::print_binary( target );
    if( rep.Esl )
    {
      USED[key]=true;
      DELS[key]=rep.levels;
      SIZS[key]=rep.area;
    }
    //printf("\n\n");
    kitty::next_inplace( target );
    //if( key > 10302 )
    //  break;
  } while ( !kitty::is_const0( target ) );

  std::ofstream myfile;
  myfile.open ("SYM10_SYN_0_0_0_0.txt");
  for( uint32_t i{0}; i<DELS.size(); ++i )
  {
    if( USED[i] )
      myfile << i << " " << DELS[i] << " " << SIZS[i] << "\n";
  }
  myfile.close();

  return 0;
}

  template<typename Ntk>
  report_t<Ntk> sym_solve( kitty::dynamic_truth_table * pF, library__t lib )
  {
    report_t<Ntk> rep;

    std::vector<kitty::dynamic_truth_table> xs;
    for( uint32_t i{0}; i<pF->num_vars(); ++i )
    {
      xs.emplace_back( pF->num_vars() );
      kitty::create_nth_var( xs[i], i );
    }
    /* define the parameters */
    int nIters = 10;
    cusco_ps ps( ccgame::solver_t::_SYM_RDE, nIters, lib );
    ps.T = std::vector<double>{0.,0.,0.,0.};
    /* solve */
    cusco<Ntk> solver( xs, {*pF} );
    rep = solver.solve( ps );


    if( rep.Esl )
    {
      default_simulator<kitty::dynamic_truth_table> sim( (*pF).num_vars() );
      const auto tt = simulate<kitty::dynamic_truth_table>( rep.ntk, sim )[0];
      assert( kitty::equal( tt, *pF ) );    
      return rep;
    }
    else
      rep.Esl = false;
    return rep;
}