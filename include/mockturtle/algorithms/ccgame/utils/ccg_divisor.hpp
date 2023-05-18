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
  \brief data structure for storing the divisors for the ccgame

  \author Andrea Costamagna
*/
#pragma once

#include <kitty/partial_truth_table.hpp>
#include <stdio.h>
#include <stack>

namespace mockturtle
{

namespace ccgame
{

using TT = kitty::partial_truth_table;

/*! \brief Gate type in the ccgame namespace.
 *  [ Xl 1100 ] [ Xr 1010 ]
 */
enum gate_t : uint8_t
{
    PIS = 0xF0,
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
    POS = 0xFF
};

class divisor_t
{

public:
  TT tt;
  gate_t  gate;
  uint32_t inL;
  uint32_t inR;
  uint32_t id;
  uint32_t flags;

  divisor_t( TT, gate_t, uint32_t, uint32_t, uint32_t, uint32_t );
  divisor_t();
  ~divisor_t();

  TT graph();
  
};

#pragma region constructors
divisor_t::divisor_t(){}
divisor_t::divisor_t( TT tt, gate_t gate, uint32_t inL, uint32_t inR, uint32_t id, uint32_t flags ):
                                tt(tt),
                                gate(gate),
                                inL(inL),
                                inR(inR),
                                id(id),
                                flags(flags)
                                {}

divisor_t::~divisor_t(){}

/*! \brief represent a truth table as an informatyion graph */
TT divisor_t::graph()
{
    int nbits = tt.num_bits();
    TT  graph( nbits * nbits );
    TT  xlarge( nbits * nbits );
    TT  mlarge( nbits * nbits );
    assert( kitty::is_const0( graph ) );
    assert( kitty::is_const0( xlarge ) );
    for( int b{0}; b < nbits; ++b )
    {
        kitty::set_bit( mlarge, b );
        if( kitty::get_bit( tt, b ) == 1 )
            kitty::set_bit( xlarge, b );
        else
            kitty::clear_bit( xlarge, b );
    }

    for( int b{nbits-1}; b >= 0; --b )
    {
        if( kitty::get_bit( tt, b ) == 0 )
            graph |= ( xlarge << nbits*b );
        else
            graph |= ( ( xlarge ^ mlarge ) << nbits*b );
    }
    return graph;
}
#pragma endregion

} // namespace ccgame

} // namespace mockturtle