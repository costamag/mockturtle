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
  \file ccg_divisor.hpp
  \brief divisors for the ccgame

  \author Andrea Costamagna
*/
#pragma once

#include <kitty/partial_truth_table.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/bit_operations.hpp>
#include <kitty/operators.hpp>
#include <stdio.h>
#include <stack>

namespace mockturtle
{

namespace ccgame
{

using TT = kitty::dynamic_truth_table;

/*! \brief Gate type in the ccgame namespace. Convention Xl=1100, Xr=1010
 */
enum gate_t : uint8_t
{
    PIS  = 0xF0,
    // -------- // direct space 
    CNTR = 0x0, // 0000
    AI00 = 0x1, // 0001
    AI01 = 0x2, // 0010
    CMPL = 0x3, // 0011
    AI10 = 0x4, // 0100
    CMPR = 0X5, // 0101
    EXOR = 0X6, // 0110
    OI11 = 0X7, // 0111
    // -------- // dual space 
    AI11 = 0X8, // 1000
    XNOR = 0X9, // 1001
    PRJR = 0XA, // 1010
    OI10 = 0XB, // 1011
    PRJL = 0XC, // 1100
    OI01 = 0XD, // 1101
    OI00 = 0XE, // 1110
    TAUT = 0XF, // 1111
    POS  = 0xFF
};

class node_t
{

public:
  /*! \brief simulation pattern */
  TT tt;
  /*! \brief [8  bits: gate type] */
  gate_t gate;
  /*! \brief [16 bits external id][16 bits internal id] */
  uint32_t id;
  /*! \brief [32 bits: left-fanin identifier] */
  uint32_t idL;
  /*! \brief [32 bits: right-fanin identifier] */
  uint32_t idR;
  /*! \brief [1 bit NOT remapped 31 bits remapped pi] */
  uint32_t idPi;
  /*! \brief delay */
  int level;

public:
  node_t( TT, gate_t, uint32_t, uint32_t, uint32_t );
  node_t();
  ~node_t();

  /* properties */
  uint32_t get_loc_id();
  uint32_t get_glb_id();
  uint32_t get_loc_idL();
  uint32_t get_glb_idL();
  uint32_t get_loc_idR();
  uint32_t get_glb_idR();

  bool is_remapped();
  uint32_t remapped_pi();
  /* graph representation */
  TT graph();
};

#pragma region constructors
node_t::node_t(){}
node_t::node_t( TT tt, gate_t gate, uint32_t id, uint32_t idL, uint32_t idR ):
                                                                        tt(tt),
                                                                        gate(gate),
                                                                        id(id),
                                                                        idL(idL),
                                                                        idR(idR){}
node_t::~node_t(){}
#pragma endregion constructors

uint32_t node_t::get_loc_id(){ return id & 0x0000FFFF;  }
uint32_t node_t::get_glb_id(){ return 0x0000FFFF & (( id & 0xFFFF0000 ) >> 16u);  }
uint32_t node_t::get_loc_idL(){ return idL & 0x0000FFFF;  }
uint32_t node_t::get_glb_idL(){ return 0x0000FFFF & (( idL & 0xFFFF0000 ) >> 16u);  }
uint32_t node_t::get_loc_idR(){ return idR & 0x0000FFFF;  }
uint32_t node_t::get_glb_idR(){ return 0x0000FFFF & (( idR & 0xFFFF0000 ) >> 16u);  }

bool node_t::is_remapped(){ return ( idPi & 0x80000000 ) != 0x80000000; }
uint32_t node_t::remapped_pi(){ return idPi; }

#pragma region graph representation
/*! \brief { X : self } { Y : information graph as 2^(2n)-bits truth table } */
TT node_t::graph()
{
    int nBits = tt.num_bits();
    //int nBit2 = nBits*nBits;
    int nVars = tt.num_vars();
    TT  graph( 2*nVars );
    TT  tt2( 2*nVars );
    TT  mk2( 2*nVars );
    for( int iBit{0}; iBit < nBits; ++iBit )
    {
        kitty::set_bit( mk2, iBit );
        if( kitty::get_bit( tt, iBit ) == 1 )
            kitty::set_bit( tt2, iBit );
        else
            kitty::clear_bit( tt2, iBit );
    }

    for( int iBit{nBits-1}; iBit >= 0; --iBit )
    {
        if( kitty::get_bit( tt, iBit ) == 0 )
            graph |= ( tt2 << nBits*iBit );
        else
            graph |= ( ( tt2 ^ mk2 ) << nBits*iBit );
    }
    return graph;
}
#pragma endregion graph representation


} // namespace ccgame

} // namespace mockturtle