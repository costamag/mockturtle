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
#include "ccg_mcnodes.hpp"
#include "ccg_rng.hpp"
#include <stdio.h>
#include <stack>
#include <iostream>

namespace mockturtle
{

namespace ccgame
{

using TT = kitty::dynamic_truth_table;


template<class Node>
class mctree_t
{
  public:
    std::vector<Node> nodes;
    Node root;
  public:
    mctree_t( Node );
    ~mctree_t();

    /*! \brief select a node at random */
    int select_random();
    int expand_random( int );
    int simulate_random( int );
    int check_closure( int );
    int check_closure( Node );
};

template<class Node> mctree_t<Node>::mctree_t( Node root ): root(root), nodes({root}){}
template<class Node> mctree_t<Node>::~mctree_t(){}


template<class Node>
int mctree_t<Node>::check_closure( int idNd )
{
  cut_t cutClose = nodes[idNd].net.check_closure( nodes[idNd].candidates );
  if( cutClose.size() > 0 )
  {
    Node nd( nodes[idNd].net );
    nd.setId(nodes.size());
    nd.isValid = false;
    nd.isExhausted = true;
    nd.net.add_cut( cutClose );
    nodes.push_back(nd);
    return nodes.size()-1;
  }
  else
    return -1;
}

template<class Node>
int mctree_t<Node>::check_closure( Node exNd )
{
  cut_t cutClose = exNd.net.check_closure( exNd.candidates );
  if( cutClose.size() > 0 )
  {
    Node nd( exNd.net );
    nd.setId(nodes.size());
    nd.isValid = false;
    nd.isExhausted = true;
    nd.net.add_cut( cutClose );
    nodes.push_back(nd);
    return nodes.size()-1;
  }
  else
    return -1;
}

template<class Node>
int mctree_t<Node>::select_random( )
{

  Node nd = nodes[0];

  while( true )
  {
    if( nd.isExhausted && nd.children.size() == 0 )
    {
      return -1;
    }
    std::uniform_int_distribution<> distrib(0, nd.children.size());
    int rnum = distrib(ccg_gen);
    if( !nd.isExhausted && rnum == nd.children.size() )
    {
      return nd.id;
    }
    else
    {
      uint32_t iNode = nd.children[rnum];
      nd = nodes[iNode];
      return iNode;
    }
  }
  return -1;
}

template<class Node>
int mctree_t<Node>::expand_random( int idPar )
{
  Node child( nodes[idPar].net );
  child.table.greedy_set_covering( -1 );  
  child.setId(nodes.size());
  child.isValid = false;
  std::uniform_int_distribution<> distrib(0, child.table.subsets.size()-1);
  for( int iTry{0}; iTry < 100; iTry++ )
  {
    int rnum = distrib(ccg_gen);
    std::vector<int> SelIds = child.table.subsets[rnum];
    std::sort(SelIds.begin(), SelIds.end()); 

    if( nodes[idPar].usedSets.find( SelIds ) == nodes[idPar].usedSets.end() )
    { 
      child.isValid = true; 
      nodes[idPar].usedSets.insert( SelIds );
      cut_t new_c;
      for( int i{0}; i < SelIds.size(); ++i )
      {
        new_c.add_node( child.candidates.nodes[SelIds[i]] );
      }
      child.net.add_cut( new_c );
      nodes.push_back(child);
      nodes[idPar].children.push_back( child.id );
      return nodes.size()-1;
    }
  }
  return -1;
}

template<class Node>
int mctree_t<Node>::simulate_random( int idPar )
{
  bool isSat{false};
  Node child( nodes[idPar].net );
  while( true )
  {
    Node nd( nodes[idPar].net );
    child = nd;
    child.setId(nodes.size());
    child.isValid = false;
    int idEnd = check_closure( child );
    if( idEnd >= 0 ) return idEnd;
    int isVanilla=0;

      child.table.greedy_set_covering( 10 );  
      std::uniform_int_distribution<> distrib(0, child.table.subsets.size()-1);
      for( int iTry{0}; iTry < 100; iTry++ )
      {
        int rnum = distrib(ccg_gen);
        std::vector<int> SelIds = child.table.subsets[rnum];
        std::sort(SelIds.begin(), SelIds.end()); 
        if( nodes[idPar].usedSets.find( SelIds ) == nodes[idPar].usedSets.end() )
        { 
          child.isValid = true; 
          nodes[idPar].children.push_back( child.id );
          nodes[idPar].usedSets.insert( SelIds );
          cut_t new_c;
          for( int i{0}; i < SelIds.size(); ++i )
            new_c.add_node( child.candidates.nodes[SelIds[i]] );
          child.net.add_cut( new_c );
          nodes.push_back(child);

          idPar = nodes.size()-1;
          break;
        }
      }
  }
  return -1;
  
}

} // namespace ccgame

} // namespace mockturtle