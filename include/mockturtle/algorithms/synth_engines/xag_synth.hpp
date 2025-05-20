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
  \file xag_synth.hpp
  \brief Synthesis through decomposition and database look-up for AIGs or XAGs.

  This header provides an efficient synthesis technique based on the following papers:
  - "An Enhanced Resubstitution Algorithm for Area-Oriented Logic Optimization", ISCAS 2024.
  - "Symmetry-Based Synthesis for Interpretable Boolean Evaluation", VLSID 2025.

  \author Andrea Costamagna
*/

#pragma once

#include "../../utils/index_list/index_list.hpp"
#include "../../utils/stopwatch.hpp"

#include <kitty/kitty.hpp>

namespace mockturtle
{

struct xag_synth_stats
{
  /*! \brief Time for finding 0-resub and collecting unate literals. */
  stopwatch<>::duration time_lookup{ 0 };

  /*! \brief Time for dividing the target and recursive call. */
  stopwatch<>::duration time_divide{ 0 };

  void report() const
  {
    fmt::print( "[i]         <xag_synth_decompose>\n" );
    fmt::print( "[i]             look-up      : {:>5.2f} secs\n", to_seconds( time_lookup ) );
    fmt::print( "[i]             dividing     : {:>5.2f} secs\n", to_seconds( time_divide ) );
  }
};

/*! \brief Logic synthesis engine for AIGs or XAGs.
 *
 * The algorithm combines symmetry-based synthesis, decomposition, and SPFD synthesis.
 *
 * The algorithm accepts an incompletely specified function, which must be provided 
 * by specifying the onset and the careset with two truth tables of the same type.
 * The algorithm then automatically synthesizes the function combining the following
 * techniques:
 * - Symmetry-based remapping.
 * - SPFD-based remapping.
 * - Database rewriting with don't cares.
 * - Top Disjoint support decomposition.
 * - Shannon Decomposition.
 *
   \verbatim embed:rst

   Example

   .. code-block:: c++

      constexpr uint32_t NumVars = ...;
      using TT = kitty::static_truth_table<NumVars>;
      const TT onset = ..., careset = ...;
      kitty::ternary_truth_table<TT> func( onset, careset );
      xag_synth_stats st;
      xag_synth_decompose resyn<NumVars>( st );
      uint32_t num_inputs = ...;
      auto result = resyn( func, num_inputs );
   \endverbatim
 */
template<uint32_t NumVars>
class xag_synth_decompose
{
public:
  using stats = xag_synth_stats;
  using index_list_t = large_xag_index_list;
  using list_element_t = index_list_t::element_type;

  using truth_table_t = kitty::static_truth_table<4u>;

  template<uint32_t N>
  using compl_specified_t = kitty::static_truth_table<N>;

  template<uint32_t N>
  using incompl_specified_t = kitty::ternary_truth_table<compl_specified_t<N>>;


public:
  explicit xag_synth_decompose( stats& st ) noexcept
      : st( st )
  {
    for ( auto i = 0u; i < 4u; ++i )
    {
      kitty::create_nth_var( proj_fns[i], i );
    }
    sims.resize( NumVars + 1 );
    sims.emplace_back();
    incompl_specified_t<NumVars> tmp;
    for ( auto i = 0u; i < NumVars; ++i )
    {
      kitty::create_nth_var( tmp, i );
      incompl_specified_t<NumVars> tt( tmp );
      sims[i+1] = tt;
    }
  }

  /*! \brief Perform XAIG synthesis from an i9ncompletely specified function.
   *
   * Resets the globally defined index list and calls the recursive synthesis engine.
   *
   * \param func Incompletely specified truth table ( `kitty::ternary_truth_table` ).
   * \param num_inputs Number of inputs of the index list to be synthesized.
   */
  void operator()( incompl_specified_t<NumVars> func )
  {
    index_list.clear();
    /* The number of inputs might be larger than the number of variables */
    index_list.add_inputs( NumVars );
    support.resize( NumVars );
    sims.resize( NumVars + 1 );

    std::iota( support.begin(), support.end(), 1u );
    list_element_t const lit = recursive_synthesis( support, func );
    index_list.add_output( lit );
  }

  list_element_t recursive_synthesis( std::vector<list_element_t> & support, incompl_specified_t<NumVars> & func )
  {
    func._bits &= func._care;
    auto supp = kitty::min_base_inplace( func );
    std::cout << supp.size() << std::endl;
    if ( supp.size() == 0 )
    {
      return index_list.get_constant( !kitty::is_const0( func._bits & func._care ) );
    }
    else if ( supp.size() == 1 )
    {
      auto const index = support[0];
      list_element_t const lit = index_list.get_literal( index );
      return kitty::equal( func, sims[index] ) ? lit : index_list.add_not( lit );
    }
    else
    {
      for ( auto it = supp.rbegin(); it != supp.rend(); ++it )
      {
        support.erase( support.begin() + *it );
      }
    }
  }

  index_list_t const& get_list() const
  {
    return index_list;
  }

private:
  /*! \brief Projection functions of the 4-dimensional Boolean space. */
  std::array<compl_specified_t<4u>, 4u> proj_fns;
  /*! \brief Projection functions of the 4-dimensional Boolean space. */
  std::vector<incompl_specified_t<NumVars>> sims;
  /*! \brief Global index list synthesized by a run of the engine. */
  index_list_t index_list;
  /*! \brief Index of a variable in the index list. */
  std::vector<list_element_t> support;

  stats& st;
};

} /* namespace mockturtle */