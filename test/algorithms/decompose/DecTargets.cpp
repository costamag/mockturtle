#include <catch.hpp>

#include <mockturtle/algorithms/decompose/DecTargets.hpp>
#include <mockturtle/networks/aig.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <kitty/operations.hpp>
#include <kitty/print.hpp>

using namespace mockturtle;

TEST_CASE( "Targets handling dynamic truth tables", "[DecTargets]" )
{
  using TT  = kitty::dynamic_truth_table;
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
  typedef DecTargets<TT> targets_t;
  targets_t tars;
  CHECK( tars.size() == 0 );
  int t0 = tars.insert( vFuncs[0], vMasks[0] );
  CHECK( t0 == 0 );
  CHECK( tars.size() == 1 );
  CHECK( * tars.getFuncP(0) == vFuncs[0] );
  CHECK( * tars.getMaskP(0) == vMasks[0] );
  int t1 = tars.insert( vFuncs[1], vMasks[1] );
  CHECK( t1 == 1 );
  CHECK( tars.size() == 2 );
  CHECK( * tars.getFuncP(1) == vFuncs[1] );
  CHECK( * tars.getMaskP(1) == vMasks[1] );
  tars.remove(0);
  CHECK( tars.size() == 1 );
  CHECK( kitty::is_const0(* tars.getFuncP(0)) );
  CHECK( kitty::is_const0( ~(* tars.getMaskP(0))) );
  int t2 = tars.insert( vFuncs[2], vMasks[2] );
  CHECK( t2 == 0 );
  CHECK( * tars.getFuncP(0) == vFuncs[2] );
  CHECK( * tars.getMaskP(0) == vMasks[2] );
}
