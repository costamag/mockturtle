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

namespace mockturtle
{

namespace mcts
{

enum supp_selection_t
{
    SUP_RAND = 0,
    SUP_ENER = 1,
    SUP_ENUM = 2,
    SUP_GENE = 3
};

enum node_selection_t
{
    NODE_RAND = 0,
    NODE_UCT = 1
};

struct node_ps
{
    supp_selection_t sel_type;
    int nIters;
    double BETA0;
    double BETAZ;
    node_ps()
    {
        sel_type = supp_selection_t::SUP_ENER;
        nIters = 1;
        BETA0 = 1000;
        BETAZ = 0;    
    }
};

struct mct_method_ps
{
    bool verbose    = false;
    node_selection_t sel_type = node_selection_t::NODE_RAND;
    mct_method_ps(){};
};

enum class entropy_t
{
    MINF,
    GINI,
    SHAN,
    EN01
};

std::mt19937 ml_gen(5);

} // namespace mcts

} // namespace mockturtle