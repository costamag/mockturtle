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

#include "ccg_node.hpp"
#include <kitty/dynamic_truth_table.hpp>
#include <stdio.h>
#include <stack>

namespace mockturtle
{

namespace ccgame
{

using DTT = kitty::dynamic_truth_table;
class cut_t
{
public:
  /*! \brief cut identifier*/
  uint32_t id{0u};
  uint32_t shiftId{0u};
  /*! \brief number of nodes */
  uint32_t nNodes{0u};
  /*! \brief nodes stored in the cut */
  std::vector<node_t> nodes;
  /*! \brief cut functionality: !can differ from simulation pattern! */
  DTT tt;
  /*! \brief mask related to the cut functionality  */
  DTT mk;

  cut_t();
  ~cut_t();

  void set_id( uint32_t );
  void set_func( DTT );
  void set_mask( DTT );
  uint32_t get_id();
  uint32_t get_shifted_id();

  node_t add_node( TT, gate_t, uint32_t, uint32_t );
  node_t add_node( TT, gate_t );
  node_t add_node( node_t );
  int  size();
  void print();
};

#pragma region constructors
cut_t::cut_t(){}
cut_t::~cut_t(){}
#pragma endregion

void cut_t::set_id( uint32_t identifier )
{ 
  id = identifier;
  shiftId = id << 16u; 
}
void cut_t::set_func( DTT func ){ tt = func; }
void cut_t::set_mask( DTT mask ){ mk = mask; }
uint32_t cut_t::get_id(){ return id; }
uint32_t cut_t::get_shifted_id(){ return shiftId; }

#pragma region addition

/*! \brief Add leaf to the cut from complete specification for generic node. */
node_t cut_t::add_node( TT tt, gate_t gate, uint32_t idL, uint32_t idR )
{
  uint32_t nodeId = shiftId | nNodes;
  node_t node( tt, gate, nodeId, idL, idR );
  if( gate == PIS )
    node.idPi = nNodes;
  nNodes++;
  nodes.push_back( node );
  return node;
}
/*! \brief Add leaf to the cut after setting the identifier. */
node_t cut_t::add_node( node_t node )
{
  node.id = shiftId | nNodes++;
  nodes.push_back( node );
  return node;
}
#pragma endregion addition

#pragma region properties
int cut_t::size(){   return nodes.size();    }
#pragma endregion properties

#pragma region visualize
void cut_t::print()
{
    for( int j{0}; j < nodes.size(); ++j )
    {
        node_t node = nodes[j];
        uint32_t x = node.get_loc_id();
        uint32_t xL = node.get_loc_idL();
        uint32_t xR = node.get_loc_idR();
        uint32_t c = node.get_glb_id();
        uint32_t cL = node.get_glb_idL();
        uint32_t cR = node.get_glb_idR();
        switch ( node.gate )
        {
            case gate_t::PIS  : { printf("[ PI %d.%2d]", c, x); break; }
            case gate_t::CNTR : { printf("[00 %d]", x); break; }
            case gate_t::AI00 : { printf("[%d.%d=and( %d.%2d', %d.%2d' )]", c, x, cL, xL, cR, xR ); break; }
            case gate_t::AI01 : { printf("[%d.%d=and( %d.%2d', %d.%2d  )]", c, x, cL, xL, cR, xR ); break; }
            case gate_t::CMPL : { printf("[%d.%d=not(    %d.%2d     )]", c, x, xL ); break; }
            case gate_t::AI10 : { printf("[%d.%d=and( %d.%2d , %d.%2d' )]", c, x, cL, xL, cR, xR ); break; }
            case gate_t::CMPR : { printf("[%d.%d=not(    %d.%2d     )]", c, x, xR ); break; }
            case gate_t::EXOR : { printf("[%d.%d=xor( %d.%2d , %d.%2d  )]", c, x, cL, xL, cR, xR ); break; }
            case gate_t::OI11 : { printf("[%d.%d=and( %d.%2d', %d.%2d' )]", c, x, cL, xL, cR, xR ); break; }
            case gate_t::AI11 : { printf("[%d.%d=and( %d.%2d , %d.%2d  )]", c, x, cL, xL, cR, xR ); break; }
            case gate_t::XNOR : { printf("[%d.%d=xor( %d.%2d', %d.%2d' )]", c, x, cL, xL, cR, xR ); break; }
            case gate_t::PRJR : { printf("[%d.%d=buf(    %d.%2d     )]", c, x, xR ); break; }
            case gate_t::OI10 : { printf("[%d.%d=and( %d.%2d', %d.%2d  )]", c, x, cL, xL, cR, xR ); break; }
            case gate_t::PRJL : { printf("[%d.%d=buf(    %d.%2d     )]", c, x, xL ); break; }
            case gate_t::OI01 : { printf("[%d.%d=and( %d.%2d , %d.%2d' )]", c, x, cL, xL, cR, xR ); break; }
            case gate_t::OI00 : { printf("[%d.%d=and( %d.%2d , %d.%2d  )]", c, x, cL, xL, cR, xR ); break; }
            case gate_t::TAUT : { printf("[11 %d.%2d]", x); break; }
            case gate_t::POS  : { printf("[ PO %d.%2d]", x); break; }
            default:  break;
        }
    }
    printf("\n");
}
#pragma endregion visualize

} // namespace ccgame

} // namespace mockturtle