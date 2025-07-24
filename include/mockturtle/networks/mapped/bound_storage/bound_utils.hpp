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
 * \brief Design type adopted by the Bound network
 *
 */
enum class design_type_t : uint8_t
{
  ARRAY_BASED, //!< Gate array-based design type
  CELL_BASED,  //!< Standard cell-based design type
};

/*!
 * \brief Computes the number of bits required to represent a given number.
 *
 * This function calculates the minimum number of bits needed to represent
 * an integer `n` in binary. When wsriting this function, only two-output
 * gates are present in technology libraries, so the maximum number of outputs
 * is 4. The function uses static assertions to ensure that the maximum number
 * of outputs does not exceed 4 and is greater than or equal to 0. Modifying
 * the maximum number of outputs requires changing the static assertions and
 * the return values accordingly.
 *
 * \param n The number to compute bits for.
 * \return The number of bits required.
 */
template<uint32_t MaxNumOutputs>
constexpr uint32_t bits_required()
{
  static_assert( MaxNumOutputs <= 4u, "num_outputs must be less than or equal to 4" );
  static_assert( MaxNumOutputs > 0u, "num_outputs must be larger than or equal to 0" );
  if constexpr ( MaxNumOutputs <= 2u )
  {
    return 1u; // One bit is enough for up to two outputs
  }
  else
  {
    return 2u; // Two bits are enough for up to 4 outputs
  }
}

/*!
 * \brief Describes the logical or structural role of a node’s output pin.
 *
 * These types are used to classify each output pin within the bound network.
 * Some types reflect logic roles (e.g., CONSTANT, PI), while others support
 * sequential mapping (e.g., CI/CO for flip-flop inputs/outputs).
 */
enum class pin_type_t : uint8_t
{
  CONSTANT = 0b00000001, //!< Constant node (logic 0 or 1)
  INTERNAL = 0b00000010, //!< Internal node within the network
  NONE = 0b00000100,     //!< No type assigned or invalid
  DEAD = 0b00001000,     //!< Node marked as dead (not used)
  PI = 0b00010000,       //!< Primary input
  PO = 0b00100000,       //!< Primary output
  CI = 0b01000000,       //!< Combinational input (e.g., from flip-flop)
  CO = 0b10000000        //!< Combinational output (e.g., to flip-flop)
};

// Bitwise NOT
constexpr pin_type_t operator~( pin_type_t rhs )
{
  return static_cast<pin_type_t>(
      ~static_cast<uint8_t>( rhs ) );
}

// Bitwise OR operator for convenience
constexpr pin_type_t operator|( pin_type_t lhs, pin_type_t rhs )
{
  return static_cast<pin_type_t>(
      static_cast<uint8_t>( lhs ) | static_cast<uint8_t>( rhs ) );
}

// Bitwise AND operator
constexpr pin_type_t operator&( pin_type_t lhs, pin_type_t rhs )
{
  return static_cast<pin_type_t>(
      static_cast<uint8_t>( lhs ) & static_cast<uint8_t>( rhs ) );
}

inline pin_type_t& operator|=( pin_type_t& lhs, pin_type_t rhs )
{
  lhs = lhs | rhs; // Use previously defined operator|
  return lhs;
}

inline pin_type_t& operator&=( pin_type_t& lhs, pin_type_t rhs )
{
  lhs = lhs & rhs; // Use previously defined operator|
  return lhs;
}

constexpr bool has_intersection( pin_type_t target, pin_type_t query )
{
  return static_cast<uint8_t>( target & query ) > 0;
}

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
      : id( id ), fanout_count( 0 ), type( type ), fanout( fanout )
  {}

  output_pin_t( uint32_t id, pin_type_t type ) noexcept
      : output_pin_t( id, type, {} )
  {}

  output_pin_t() noexcept
      : output_pin_t( std::numeric_limits<uint32_t>::max(), pin_type_t::NONE, {} )
  {}

  /*! \brief Identifier of the pin’s function in the gate (used for mapping) */
  uint32_t id;

  /*! \brief Identifier of the pin’s function in the gate (used for mapping) */
  uint32_t fanout_count;

  /*! \brief Logical type of the pin (PI, PO, constant, etc.) */
  pin_type_t type;

  /*! \brief List of nodes that receive this output as input */
  std::vector<node_index_t> fanout;
};

} // namespace bound

} // namespace mockturtle
