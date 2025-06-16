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
  \file bound_types.hpp
  \brief Basic types and enumerations used in the bound network data structure.

  This file defines types related to node indexing and output pin behavior in the
  bound storage network, including logic and mapping-related pin classifications.

  \ingroup bound_storage
  \see mockturtle::bound::storage
  \see mockturtle::bound::storage_node
  \see mockturtle::bound::storage_signal

  \author Andrea Costamagna
*/

#pragma once

#include <cstdint>
#include <limits>
#include <vector>

namespace mockturtle
{

/*!
 * \namespace bound
 * \brief Types and utilities related to the bound network data structure.
 */
namespace bound
{

/*!
 * \brief Describes the logical or structural role of a node’s output pin.
 *
 * These types are used to classify each output pin within the bound network.
 * Some types reflect logic roles (e.g., CONSTANT, PI), while others support
 * sequential mapping (e.g., CI/CO for flip-flop inputs/outputs).
 */
enum class pin_type_t : uint8_t
{
  CONSTANT,    //!< Constant node (logic 0 or 1)
  INTERNAL,    //!< Internal node within the network
  NONE,        //!< No type assigned or invalid
  DEAD,        //!< Node marked as dead (not used)
  PI,          //!< Primary input
  PO,          //!< Primary output
  CI,          //!< Combinational input (e.g., from flip-flop)
  CO           //!< Combinational output (e.g., to flip-flop)
};

/*! \brief Type used to identify a node within the bound network.
 *
 * Typically used as an index into node storage containers.
 */
using node_index_t = uint64_t;

/*! \brief Describes a specific output pin of a logic gate or node.
 *
 * Nodes can have multiple output pins to support multi-output gates.
 * Each output pin is identified by an `id` corresponding to its position in the
 * gate's output function list (as defined by the technology library).
 *
 * The `fanout` vector tracks which other nodes this output connects to.
 */
struct output_pin_t
{

  output_pin_t( uint32_t id, pin_type_t type, std::vector<node_index_t> const& fanout ) noexcept
   : id( id ), type( type ), fanout( fanout )
  {}

  output_pin_t( uint32_t id, pin_type_t type ) noexcept
   : output_pin_t( id, type, {} )
  {}

  output_pin_t() noexcept
   : output_pin_t( std::numeric_limits<uint32_t>::max(), pin_type_t::NONE, {} )
  {}

/*! \brief Identifier of the pin’s function in the gate (used for mapping) */
  uint32_t id;

  /*! \brief Logical type of the pin (PI, PO, constant, etc.) */
  pin_type_t type;

  /*! \brief List of nodes that receive this output as input */
  std::vector<node_index_t> fanout;
};

} // namespace bound

} // namespace mockturtle
