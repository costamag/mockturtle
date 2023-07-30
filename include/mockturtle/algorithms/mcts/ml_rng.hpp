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
  \brief data structure for storing the cuts for the mcts

  \author Andrea Costamagna
*/
#pragma once

#include <stdio.h>
#include <stack>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/partial_truth_table.hpp>
#include <kitty/operations.hpp>
#include <kitty/operators.hpp>
#include <random>

namespace mockturtle
{

namespace mcts
{


//detailed_gate_t not_( gate_t::CMPR, 1, 1*0.5, 1.0, &hpcompute_not );//0.5
//detailed_gate_t nor_( gate_t::AI00, 2, 1*1.0, 1.0, &hpcompute_ai00 );
//detailed_gate_t and_( gate_t::AI11, 2, 1*1.5, 1.0, &hpcompute_ai11 );
//detailed_gate_t xor_( gate_t::EXOR, 2, 1*2.0, 1.0, &hpcompute_exor );

std::mt19937 ml_gen(5);

} // namespace mcts

} // namespace mockturtle