/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file boolean_chain.hpp
  \brief data structure for storing boolean function representations

  \author Andrea Costamagna
*/
#pragma once
#include "../../../networks/xag.hpp"
#include "../utils/ccg_net.hpp"
#include "../utils/ccg_rng.hpp"
#include "../utils/ccg_mcnodes.hpp"
#include "../utils/ccg_mctree.hpp"
#include <kitty/partial_truth_table.hpp>
#include <kitty/bit_operations.hpp>
#include <kitty/operators.hpp>
#include <kitty/print.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <stdio.h>
#include <stack>
#include <set>

namespace mockturtle
{

namespace ccgame
{

using TT = kitty::dynamic_truth_table;


template<class Ntk>
struct report_cov_t
{
  bool isDone{false};
  bool isFound{false};
  int nIt0;
  int nMin;
  int nMax;
  Ntk ntk;
};

struct cusco_cov_ps
{
  /* method */
  /*! \brief number of iterations */
  int nIters;
  /*! \brief capacity: #candidates considered */
  int nCap;
  /*! \brief capacity: #candidates considered */
  bool dcMax;
  cusco_cov_ps( int nIters, int nCap ) : nIters( nIters ), nCap(nCap) { dcMax=false; }
  cusco_cov_ps( int nIters, int nCap, bool dc_maximize ) : nIters( nIters ), nCap(nCap), dcMax(dc_maximize) {}
};

template<class Ntk>
class cusco_cov
{
public:
    /* problem definition */
    std::vector<TT> X; // input simulations
    std::vector<TT> Y; // output simulations

public:
  /* creation and destruction */
  cusco_cov( std::vector<TT> const&, std::vector<TT> const& );
  ~cusco_cov();
  /* solve */
  report_cov_t<Ntk> solve_random( cusco_cov_ps const& );
  report_cov_t<Ntk> solve_mcts( cusco_cov_ps const& );
  report_cov_t<Ntk> solve_genetic( cusco_cov_ps const& );
  report_cov_t<Ntk> recut( net_t, cusco_cov_ps const&, int );

};

/* creation and destruction */
template<class Ntk>
cusco_cov<Ntk>::cusco_cov( std::vector<TT> const& X, std::vector<TT> const& Y ): X(X), Y(Y) {}
template<class Ntk>
cusco_cov<Ntk>::~cusco_cov(){}

template<class Ntk>
report_cov_t<Ntk> cusco_cov<Ntk>::solve_random( cusco_cov_ps const& ps )
{
  report_cov_t<Ntk> rep;
  rep.nMax = -1;
  rep.nMin = 10000;
  std::random_device rd;  // a seed source for the random number engine
  std::mt19937 gen(5); // mersenne_twister_engine seeded with rd()
  for( int iIt{0}; iIt < ps.nIters; ++iIt )
  {
    Ntk ntk;
    net_t net( X, Y );
    while( net.nHunging > 0 )
    {
      cut_t candidates = net.enumerate_divs();
      cut_t closed_c = net.check_closure( candidates );

      net.add_cut( closed_c );

      if( net.nHunging == 0 )
        break;

      tab_t table( candidates, net.outCut );
      table.greedy_set_covering( ps.nCap );
      //if( ps.dcMax )  
        //table.select_dc_maximizers( );

      std::uniform_int_distribution<> distrib(0, table.subsets.size()-1);
      int rnum = distrib(gen);
      std::vector<int> SelIds = table.subsets[rnum];
      
      cut_t new_c;
      for( int i{0}; i < SelIds.size(); ++i )
        new_c.add_node( candidates.nodes[SelIds[i]] );

      net.complete_cut( new_c );
       
    }
    ntk = net.convert<Ntk>();
    ntk = cleanup_dangling( ntk );
    if( iIt == 0 )  rep.nIt0 = ntk.num_gates();
    if( ntk.num_gates() < rep.nMin )
    {
      rep.ntk = ntk;
      rep.nMin = ntk.num_gates();
    }
    else if ( ntk.num_gates() > rep.nMax )  rep.nMax = ntk.num_gates();
    //printf("|gates*|= %d\n", ntk.num_gates());
  }
  //printf("END\n", rep.ntk.num_gates());
  //printf("|gates*|= %d\n", rep.ntk.num_gates());

  return rep;
}


template<class Ntk>
report_cov_t<Ntk> cusco_cov<Ntk>::solve_mcts( cusco_cov_ps const& ps )
{
  report_cov_t<Ntk> rep;
  rep.nMax = -1;
  rep.nMin = 10000;
  //net_t net(X, Y);
  mcnode_cut_t mcRootNd( X, Y );

  mcnode_cut_t * pNdSel;
  mcnode_cut_t * pNdExp;
  mcnode_cut_t * pNdSim;

  mctree_t mcTree( mcRootNd );

  int idSel;
  int idEnd;

  for( int iIt{0}; iIt < ps.nIters; ++iIt )
  {
    printf("_d\n");

    /* SELECT a node never marked as exhausted ( it could be ) */
    idSel = mcTree.select_random();
    if( idSel < 0 ) continue;
    printf("_a\n");

    idEnd = mcTree.check_closure( idSel );
    printf("_a1\n");

    if( idEnd >= 0 ) continue;
    printf("_a2\n");

    /* EXPAND with a new node out of the current one */
    int idExp = mcTree.expand_random( idSel );
    printf("_a3\n");
    printf("=idExp = %d\n", idExp );
    if( idExp < 0 ) continue;
    printf("_a4\n");

    idEnd = mcTree.check_closure( idExp );
    printf("_a5\n");

    if( idEnd >= 0 ) continue;
    printf("_a61\n");

    /* SIMULATE only if there is a new node and the leaf does not terminate the game */
    idEnd = mcTree.simulate_random( idExp );
    printf("_a7\n");

    if( idEnd < 0 ) continue;
    printf("_b\n");
    
    /* Here we found a solution if is valid */
    Ntk ntk = mcTree.nodes[idEnd].net.convert<Ntk>();
    printf("_c\n");

    //printf("|ntk|=%d\n", ntk.num_gates());
    ntk = cleanup_dangling( ntk );
    if( ntk.num_gates() < rep.nMin )
    {
      rep.ntk = ntk;
      rep.nMin = ntk.num_gates();
    }
    else if ( ntk.num_gates() > rep.nMax )  rep.nMax = ntk.num_gates();
  }
  printf("-\n");
  return rep;
}

template<class Ntk>
report_cov_t<Ntk> cusco_cov<Ntk>::recut( net_t net, cusco_cov_ps const& ps, int nGlbMin )
{ 
  //std::random_device rd;  // a seed source for the random number engine
  //std::mt19937 gen(rSeed); // mersenne_twister_engine seeded with rd()

  report_cov_t<Ntk> rep;

  if( net.nHunging == 0 )
  {
    Ntk ntk = net.convert<Ntk>();
    ntk = cleanup_dangling( ntk );
    rep.ntk = ntk;
    rep.nMin = ntk.num_gates();
    return rep;
  }

  cut_t candidates = net.enumerate_divs();
  cut_t closed_c   = net.check_closure( candidates );
  net.add_cut( closed_c );

  tab_t table( candidates, net.outCut );
  table.greedy_set_covering( ps.nCap );

  int nBest = 100000;
  
  if( table.subsets.size() < 3 )
  {
    for( int iSet{0}; iSet < table.subsets.size(); ++iSet )//setUsed.size() < std::min( (int)table.subsets.size(), 5 )
    {
      report_cov_t<Ntk> repLoc;
      net_t netLoc = net;
      cut_t new_c;
      for( int iId{0}; iId < table.subsets[iSet].size(); ++iId )
        new_c.add_node( candidates.nodes[table.subsets[iSet][iId]] );

      netLoc.complete_cut( new_c );
      repLoc = recut( netLoc, ps, nGlbMin );
      if( repLoc.ntk.num_gates() < nBest )
      {
        nBest = repLoc.ntk.num_gates();
        rep = repLoc;
      }
    }    
  }
  else
  {
    std::uniform_int_distribution<> distrib(0, table.subsets.size()-1);
    std::set<int> setUsed;
    while(setUsed.size() < 3 )
    {
      report_cov_t<Ntk> repLoc;
      net_t netLoc = net;
      cut_t new_c;

      int rnum = distrib(ccg_gen);
      if( setUsed.find( rnum ) != setUsed.end()  )
        continue;
      else
        setUsed.insert( rnum );

      std::vector<int> SelIds = table.subsets[rnum];
      for( int iId{0}; iId < SelIds.size(); ++iId )
        new_c.add_node( candidates.nodes[SelIds[iId]] );
      netLoc.complete_cut( new_c );
      repLoc = recut( netLoc, ps, nGlbMin );
      if( repLoc.ntk.num_gates() < nBest )
      {
        nBest = repLoc.ntk.num_gates();
        rep = repLoc;
      }
    }
  }
      
  return rep;

}

template<class Ntk>
report_cov_t<Ntk> cusco_cov<Ntk>::solve_genetic( cusco_cov_ps const& ps )
{  
  int it{0};
  
  report_cov_t<Ntk> rep;
  rep.nMin = 1000000;
  while( it < ps.nIters )
  {
    net_t net( X, Y );
    report_cov_t<Ntk> repLoc = recut( net, ps, rep.nMin );
    //printf("it %d nBest = %d\n", it, repLoc.nMin );
    if( repLoc.nMin < rep.nMin )
      rep = repLoc;
    it++;
  }
  return rep;
}

} // namespace ccgame

} // namespace mockturtle