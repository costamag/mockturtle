/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2023  EPFL
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
  \file homo_xpx_analyzer.hpp
  \brief Inplace rewrite

  \author Andrea Costamagna
*/

#pragma once

#include "../database/mapped_list.hpp"

#include "../windowing/window_manager.hpp"
#include "../../algorithms/simulation.hpp"
#include "../list_simulator.hpp"
#include "../timing/sensing_analyzer.hpp"
#include "../timing/arrival_analyzer.hpp"
#include "../power/node_simulator.hpp"
#include "../../algorithms/detail/resub_utils.hpp"
#include "../../algorithms/reconv_cut.hpp"

#include <fmt/format.h>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>

namespace mockturtle
{

  template<class Ntk, uint32_t num_vars>
  class static_dependency_cut
  {
public:
    using node = typename Ntk::node;
    using truth_table_t = kitty::static_truth_table<num_vars>;

    static_dependency_cut(){}

    static_dependency_cut( node const& n, std::array<node, num_vars>const& leaves, truth_table_t const& func, truth_table_t const& care )
      : root( n ),
        leaves( leaves ),
        func( func ),
        care( care )
    {}

    truth_table_t const & get_function() const
    {
      return func;
    }

    truth_table_t const & get_careset() const
    {
      return care;
    } 

    std::array<node, num_vars> const & get_leaves() const
    {
      return leaves;
    } 

    node root;
    std::array<node, num_vars> leaves;
    truth_table_t func;
    truth_table_t care;
  };

} // namespace mockturtle