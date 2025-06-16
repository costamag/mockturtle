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
  \file bind_library.hpp
  \brief Implements methods for handling and evaluating a library of standard cells.

  This engine can be used for efficient Boolean evaluation of the gates in a
  standard cell library. Each gate is represented as an AIG index list for
  efficient evaluation. Additionally, in the presence of multiple-output cells,
  this engine identifies which gates belong to a multiple-output and allows handling
  this information.

  NOTE: The augmented library can be made arbitrarily complex adding technological
  information. This data structure can be modified to store detailed information from
  the liberty file.

  \author Andrea Costamagna
*/

#pragma once

#include <vector>

#include "../../algorithms/synth_engines/xag_synth.hpp"
#include "../../io/genlib_reader.hpp"
#include "../../utils/index_lists/lists/xag_index_list.hpp"

namespace mockturtle
{

template<typename Gate>
class augmented_library
{
public:
  /* truth table used to express the gate functionality */
  using func_t = kitty::dynamic_truth_table;
  using list_t = large_xag_index_list;
  /* gate type in the raw format */
  using raw_gate_t = Gate;

  /*! \brief Augmented gate.
   *
   * A raw gate is augmented by decomposing it into an index list for
   * efficient simulation.
   */
  struct aug_gate_t : raw_gate_t
  {
    list_t aig_list;

    aug_gate_t( const raw_gate_t& g, list_t const& list )
        : raw_gate_t( g ), aig_list( list )
    {}
  };

  /*! \brief Construction via specification of the simpler library.
   *
   * The gates should specify at least the gate's functionality, from
   * which this constructor can synthesize an index list for each gate.
   */
  augmented_library( std::vector<raw_gate_t> const& raw_gates )
  : synth( st )
  {
    aug_gates.reserve( raw_gates.size() );
    for ( gate const& g : raw_gates )
    {
      add_gate( g );
    }
  }

  /*! \brief Default construction expecting gate-by-gate library characterization.
   *
   * This constructor can be used for virtual libraries, useful when storing
   * look-up tables appearing in an LUT network.
   */
  augmented_library()
  : synth( st )
  {
    /* reserve initial space for the library */
    aug_gates.reserve( 128u );
  }

  /*! \brief Augment the gate and add it to the library.
   *
   * \param raw_gate Representation containing the gate's function
   */
  void add_gate( raw_gate_t const& g )
  {
    synth( g.function );
    list_t const list = synth.get_list();
    aug_gates.emplace_back( g, list );
  }

  /*! \brief Getter of gate containing detailed information. */
  aug_gate_t const& get_augmented_gate( uint32_t id ) const
  {
    return aug_gates[id];
  }

  /*! \brief Getter of the list synthesizing the gate's functionality. */
  list_t const& get_list( uint32_t id ) const
  {
    return aug_gates[id].aig_list;
  }

private:
  /*! \brief Augmented technology library */
  std::vector<aug_gate_t> aug_gates;
  /*! \brief Synthesis engine for AIG index lists */
  xag_synth_stats st;
  xag_synth_decompose<false, false> synth;
};

} // namespace mockturtle