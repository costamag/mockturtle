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

#include "ccg_divisor.hpp"
#include <stdio.h>
#include <stack>

namespace mockturtle
{

namespace ccgame
{

class cut_t
{

public:
  std::vector<divisor_t> divisors;

  cut_t();
  ~cut_t();

  void add_divisor( TT, gate_t, uint32_t, uint32_t, uint32_t, uint32_t );
  void add_divisor( divisor_t );
  int  size();
  void print();
};

#pragma region constructors
cut_t::cut_t(){}
cut_t::~cut_t(){}
#pragma endregion

#pragma region addition
void cut_t::add_divisor( TT tt, gate_t gate, uint32_t inL, uint32_t inR, uint32_t id, uint32_t flags )
{
    divisor_t div( tt, gate, inL, inR, id, flags );
    divisors.push_back( div );
}

void cut_t::add_divisor( divisor_t div )
{
  divisors.push_back( div );
}
#pragma endregion addition

#pragma region properties
int cut_t::size(){   return divisors.size();    }
#pragma endregion properties

#pragma region visualize
void cut_t::print()
{
    for( int j{0}; j < divisors.size(); ++j )
    {
        divisor_t div = divisors[j];
        switch ( div.gate )
        {
            case gate_t::PIS  : { printf("[ PI %2d]", div.id ); break; }
            case gate_t::CNTR : { printf("[00 %d]", div.id ); break; }
            case gate_t::AI00 : { printf("[%d=and( %2d', %2d' )]", div.id, div.inL, div.inR ); break; }
            case gate_t::AI01 : { printf("[%d=and( %2d', %2d  )]", div.id, div.inL, div.inR ); break; }
            case gate_t::CMPL : { printf("[%d=not(    %2d     )]", div.id, div.inL ); break; }
            case gate_t::AI10 : { printf("[%d=and( %2d , %2d' )]", div.id, div.inL, div.inR ); break; }
            case gate_t::CMPR : { printf("[%d=not(    %2d     )]", div.id, div.inR ); break; }
            case gate_t::EXOR : { printf("[%d=xor( %2d , %2d  )]", div.id, div.inL, div.inR ); break; }
            case gate_t::OI11 : { printf("[%d=and( %2d', %2d' )]", div.id, div.inL, div.inR ); break; }
            case gate_t::AI11 : { printf("[%d=and( %2d , %2d  )]", div.id, div.inL, div.inR ); break; }
            case gate_t::XNOR : { printf("[%d=xor( %2d', %2d' )]", div.id, div.inL, div.inR ); break; }
            case gate_t::PRJR : { printf("[%d=buf(    %2d     )]", div.id, div.inR ); break; }
            case gate_t::OI10 : { printf("[%d=and( %2d', %2d  )]", div.id, div.inL, div.inR ); break; }
            case gate_t::PRJL : { printf("[%d=buf(    %2d     )]", div.id, div.inL ); break; }
            case gate_t::OI01 : { printf("[%d=and( %2d , %2d' )]", div.id, div.inL, div.inR ); break; }
            case gate_t::OI00 : { printf("[%d=and( %2d , %2d  )]", div.id, div.inL, div.inR ); break; }
            case gate_t::TAUT : { printf("[11 %2d]", div.id ); break; }
            case gate_t::POS  : { printf("[ PO %2d]", div.id ); break; }
            default:  break;
        }
    }
    printf("\n");
}
#pragma endregion visualize

} // namespace ccgame

} // namespace mockturtle