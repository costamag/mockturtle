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
  \file ccg_cut.hpp
  \brief data structure for storing the cuts for the ccgame

  \author Andrea Costamagna
*/
#pragma once

#include "ccg_net.hpp"
#include "ccg_cut.hpp"
#include "ccg_rng.hpp"
#include "ccg_supportor.hpp"

#include <stdio.h>
#include <stack>
#include <set>
#include <math.h>

namespace mockturtle
{

namespace ccgame
{

using DTT = kitty::dynamic_truth_table;

class mcnode_cut_t
{
  public:
    bool isExhausted{false};
    bool isValid{true};
    /*! \brief vector of Xi rewards sprouting from this node */
    std::vector<float> rwdPlayouts;
    std::vector<int> children;
    net_t net;
    cut_t candidates;
    int id{0};
    std::set<std::vector<int>> usedSets;
    tab_t table;
    support_generator_t supportor;
    

    mcnode_cut_t();
    mcnode_cut_t( std::vector<DTT> const&, std::vector<DTT> const& );
    mcnode_cut_t( net_t );
    ~mcnode_cut_t();

    void setId(int);
    cut_t check_closure();
  
};

mcnode_cut_t::mcnode_cut_t(){};

mcnode_cut_t::mcnode_cut_t( std::vector<DTT> const& X, std::vector<DTT> const& Y ): net( X, Y ) //initializing the network creates the first layer
{
  // there is one cut at the beginning
  assert( net.cuts.size() == 1 );

  cut_t rootCut = net.cuts[0];

  /* set the previous cut as used nodes */
  std::vector<int> vecUsed;
  for( int iDiv{0}; iDiv < rootCut.nodes.size(); ++iDiv )
    vecUsed.push_back(iDiv);
  usedSets.insert( vecUsed );

  /* find the candidates */
  candidates = net.enumerate_divs();
  table.init_tab( candidates, net.outCut );


  std::vector<divisor_t> divisors;
  for( uint32_t i{0}; i < candidates.size(); ++i )
  {
      int div_id = i;
      DTT div_tt = candidates.nodes[i].tt;
      double div_area = (double)(i>X.size());
      double div_delay = (double)(i>X.size());
      divisor_t div( div_id, div_tt, div_area, div_delay );
      divisors.push_back(div);
  }

  std::vector<target_t> targets;
  for( uint32_t i{0}; i<Y.size(); ++i )
  {
      DTT trg_tt = Y[i];
      target_t trg( i, trg_tt );
      targets.push_back( trg );
  }

    /* support genenrator initialization */
    support_generator_t suppor( divisors, targets, method_t::BASE, 2 );

};

mcnode_cut_t::mcnode_cut_t( net_t eNet ): net(eNet)
{
  cut_t lastCut = net.get_last_cut();
  /* set the previous cut as used nodes */
  std::vector<int> vecUsed;
  for( int iDiv{0}; iDiv < lastCut.nodes.size(); ++iDiv )
    vecUsed.push_back(iDiv);
  usedSets.insert( vecUsed );

  /* find the candidates */
  candidates = net.enumerate_divs();
  table.init_tab( candidates, net.outCut );
};


mcnode_cut_t::~mcnode_cut_t(){};

void mcnode_cut_t::setId( int Id ){ id = Id; }


} // namespace ccgame

} // namespace mockturtle


