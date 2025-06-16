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
  \file bound_node.hpp
  \brief Defines the core node structure for the bound network.
  \details
    This file introduces `bound::storage_node`, a data structure representing a logic node
    within the bound storage network. Each node maintains information about its fan-in,
    fan-out, user-defined metadata, and a list of output pins.

    The template parameter `MaxNumOutputBits` controls the maximum number of outputs
    a node can support, enabling compatibility with multi-output standard cells.

    Nodes can be structurally compared and marked as "dead" by setting the type of their output pins.

  \ingroup bound_storage
  \see mockturtle::bound::storage
  \see mockturtle::bound::storage_signal
  \see mockturtle::bound::storage_types

  \author Andrea Costamagna
*/

#pragma once

#include "bound_types.hpp"
#include "bound_signal.hpp"
#include <cstdint>
#include <vector>

namespace mockturtle
{

/*!
 * \namespace bound
 * \brief Types and utilities related to the bound network data structure.
 */
namespace bound
{

/*! \brief Node representation in the bound network.
 *
 * \tparam MaxNumOutputBits Maximum number of outputs supported per node.
 */
template<uint32_t MaxNumOutputBits>
struct storage_node
{
  using signal_t = storage_signal<MaxNumOutputBits>;

  /*! \brief Default constructor with a single (default) output. */
  storage_node()
  {
    outputs = decltype( outputs )( 1 );
  }

  /*! \brief Constructor that sets the pin type of the first output. */
  storage_node( pin_type_t const& type )
  {
    outputs = decltype( outputs )( 1 );
    outputs[0].type = type;
  }

  /*! \brief Marks the node as dead by setting all output pin types to DEAD. */
  void kill()
  {
    for ( auto& pin : outputs )
    {
      pin.type = pin_type_t::DEAD;
    }
  }

  /*! \brief Equality operator compares structural fan-in. */
  bool operator==( storage_node<MaxNumOutputBits> const& other ) const
  {
    return children == other.children;
  }

  /*! \brief Signals of the node's immediate fan-ins. */
  std::vector<signal_t> children;

  /*! \brief Custom user data for tagging or annotation. */
  uint32_t user_data { 0 };

  /*! \brief Traversal marker used in graph algorithms. */
  uint32_t traversal_id { 0 };

  /*! \brief Fan-out count; MSB may encode special flags (e.g., "dead"). */
  uint32_t fanout_count { 0 };

  /*! \brief Output pins associated with this node. */
  std::vector<output_pin_t> outputs;
};

} // namespace bound

} // namespace mockturtle