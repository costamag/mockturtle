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

#include <stdio.h>
#include <stack>
#include <set>
#include <math.h>

namespace mockturtle
{

namespace mcts
{

template<class TT>
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
    tab_t small_table;

    mcnode_cut_t();
    mcnode_cut_t( std::vector<TT> const&, std::vector<TT> const& );
    mcnode_cut_t( net_t );
    ~mcnode_cut_t();

    void setId(int);
    cut_t check_closure();
  
};

mcnode_cut_t::mcnode_cut_t(){};

mcnode_cut_t::mcnode_cut_t( std::vector<TT> const& X, std::vector<TT> const& Y ): net( X, Y ) //initializing the network creates the first layer
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
  small_table.init_small_tab( candidates, net.outCut );
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
  small_table.init_tab( candidates, net.outCut );
};


mcnode_cut_t::~mcnode_cut_t(){};

void mcnode_cut_t::setId( int Id ){ id = Id; }


} // namespace ccgame

} // namespace mockturtle


