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
  \file bound_storage.hpp
  \brief Storage for bound network specializing the operations on the nodes.
  \details
  This file defines the storage for the bound network, which is a specialized
  data structure designed to handle multiple-output gates and their bindings.
  It includes methods for creating primary inputs and outputs, managing nodes,
  and handling the functional properties of the network. The encapsulation of
  the storage allows for efficient manipulation of the network while maintaining
  the flexibility to support various gate functionalities and bindings.
  \note This storage is designed to work with the `bound_network` class, which
  provides a higher-level interface for interacting with the network.
  \note The `bound_network` class uses this storage to manage the nodes, inputs,
  outputs, and the library of gates. It provides methods for creating nodes,
  replacing nodes, and querying the network's structure and functionality.

  \ingroup bound_storage
  \see mockturtle::bound::storage_node
  \see mockturtle::bound::storage_signal
  \see mockturtle::bound::storage_types

  \author Andrea Costamagna
*/

#pragma once

#include "../../../io/genlib_reader.hpp"
#include "../../../utils/mapping/augmented_library.hpp"
#include "bound_node.hpp"
#include "bound_signal.hpp"
#include "bound_types.hpp"

#include <algorithm>
#include <cassert>
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

/*! \brief Compact storage for nodes in the bound network.
 *
 * This structure represents the storage in bound networks, enabling the
 * encapsulation of the detailed operations on nodes, inputs, and outputs.
 * It provides methods for creating primary inputs and outputs, managing nodes,
 * and handling the functional properties of the network. The storage is designed
 * to efficiently manage the nodes and their relationships, allowing for operations
 * such as creating nodes, replacing nodes, and querying the network's structure.
 *
 * \tparam NumBitsOutputs Number of bits used to represent the output pin specifier.
 */
template<uint32_t NumBitsOutputs>
class storage
{
public:
  using gate_t = mockturtle::gate;
  using node_t = storage_node<NumBitsOutputs>;
  using signal_t = storage_signal<NumBitsOutputs>;

  /*! \brief The storage constructor.
   *
   * This constructor initializes the storage with a given library of gates.
   * It reserves space for a maximum number of nodes and initializes the first
   * two nodes as constants (0 and 1).
   *
   * \param gates The gates to be included in the library.
   */
  storage( std::vector<gate> const& gates )
      : library( gates )
  {
    /* reserve space for nodes */
    nodes.reserve( 10000u ); // reserve to avoid frequent reallocations

    /* we reserve the first two nodes for constants */
    nodes.emplace_back( pin_type_t::CONSTANT ); // 0
    nodes.emplace_back( pin_type_t::CONSTANT ); // 1
  }

#pragma region Primary I / O and constants

  /*! \brief Creates a constant signal.
   *
   * This method creates a signal representing a constant value (0 or 1).
   * It returns a signal with the appropriate index and output pin.
   *
   * \param value The constant value to be represented (true for 1, false for 0).
   * \return A signal representing the constant value.
   */
  signal_t get_constant( bool value ) const
  {
    return value ? signal_t{ 1, 0 } : signal_t{ 0, 0 };
  }

  /*! \brief Creates a primary input signal.
   *
   * This method creates a primary input signal and adds it to the storage.
   * It returns a signal with the index of the newly created primary input.
   *
   * \return A signal representing the primary input.
   */
  signal_t create_pi()
  {
    const auto index = nodes.size();
    nodes.emplace_back( pin_type_t::PI );
    inputs.emplace_back( index );
    return signal_t{ index, 0 };
  }

  /*! \brief Creates a primary output signal.
   *
   * This method creates a primary output signal from a given signal.
   * It increases the reference count for the child node and sets the output pin type to PO.
   * It returns the index of the newly created primary output.
   *
   * \param f The signal representing the function to be used as primary output.
   * \return The index of the newly created primary output.
   */
  uint32_t create_po( signal_t const& f )
  {
    /* increase ref-count to children */
    nodes[f.index].fanout_count++;
    nodes[f.index].outputs[f.output].type = pin_type_t::PO;
    auto const po_index = static_cast<uint32_t>( outputs.size() );
    outputs.emplace_back( f.index, f.output );
    return po_index;
  }

  /*! \brief Check if the node is a multiple-output node
   */
  bool is_multioutput( node_index_t const& n ) const
  {
    return nodes[n].outputs.size() > 1;
  }

  /*! \brief Check if the node is a constant */
  bool is_constant( node_index_t const& n ) const
  {
    auto const& pins = nodes[n].outputs;
    return pins[0].type == pin_type_t::CONSTANT;
  }

  /*! \brief Check if the node is a combinational input (CI)
   *
   * \param n The node to check.
   * \return True if the node is a combinational input, false otherwise.
   */
  bool is_ci( node_index_t const& n ) const
  {
    auto const& pins = nodes[n].outputs;
    return ( pins[0].type == pin_type_t::PI ) ||
           ( pins[0].type == pin_type_t::CI );
  }

  /*! \brief Check if the node is a primary input (PI)
   *
   * \param n The node to check.
   * \return True if the node is a primary input, false otherwise.
   */
  bool is_pi( node_index_t const& n ) const
  {
    return is_ci( n );
  }

  /*! \brief Check if the node is a primary output (PO)
   *
   * \param n The node to check.
   * \param output The output pin index to check (default is 0).
   * \return True if the node is a primary output, false otherwise.
   */
  bool is_po( node_index_t const& n, uint32_t output = 0 ) const
  {
    auto const& pins = nodes[n].outputs;
    return ( pins[output].type == pin_type_t::PO ) ||
           ( pins[output].type == pin_type_t::CO );
  }

  /*! \brief Returns if node index is not 0 */
  bool constant_value( node_index_t const& n ) const
  {
    return n != 0;
  }
#pragma region

  /*! \brief Check if the node is dead
   *
   * \param n The node to check.
   * \return True if the node is dead, false otherwise.
   *
   * A dead node is one where all output pins are marked as DEAD.
   * This typically indicates that the node is no longer used in the network.
   */
  bool is_dead( node_index_t const& n ) const
  {
    bool all_dead{ true };
    bool one_dead{ false };
    for ( auto const& pin : nodes[n].outputs )
    {
      bool const dead = pin.type == pin_type_t::DEAD;
      all_dead &= dead;
      one_dead |= dead;
    }
    assert( !( all_dead ^ one_dead ) );
    /* A dead node is simply a dangling node */
    return all_dead;
  }

  /*! \brief Create a new node with multiple outputs.
   *
   * This method creates a new node with the specified children and output IDs.
   * It updates the fanout counts of the children and returns a signal representing
   * the new node.
   *
   * \param children The child signals that this node will depend on.
   * \param ids The IDs of the outputs for this node.
   * \return A signal representing the newly created node.
   */
  signal_t create_node( std::vector<signal_t> const& children, std::vector<uint32_t> const& ids )
  {
    node_t new_node;
    std::copy( children.begin(), children.end(), std::back_inserter( new_node.children ) );

    new_node.outputs = decltype( new_node.outputs )( ids.size() );

    for ( auto i = 0; i < ids.size(); ++i )
      new_node.outputs[i] = { ids[i], pin_type_t::INTERNAL };

    const auto index = nodes.size();
    nodes.push_back( new_node );

    /* increase ref-count to children */
    for ( auto c : children )
    {
      nodes[c.index].fanout_count++;
      nodes[c.index].outputs[c.output].fanout.push_back( index );
    }

    return { index, 0 };
  }

  /*! \brief Get the binding identifiers of the output pins in a node
   */
  std::vector<uint32_t> get_binding_ids( node_index_t const& n ) const
  {
    std::vector<uint32_t> ids;
    for ( auto const& pin : nodes[n].outputs )
      ids.push_back( pin.output );
    return ids;
  }

  /*! \brief Checks if a node is in the fanin of another one */
  bool in_fanin( node_index_t parent, node_index_t other ) const
  {
    bool in_fanin = false;
    auto& nobj = nodes[parent];
    for ( auto& f : nobj.children )
    {
      if ( f.index == other )
      {
        in_fanin = true;
        return in_fanin;
      }
    }

    return in_fanin;
  }

  /*! \brief Get the children of a node.
   *
   * This method retrieves the children of a specified node.
   * It returns a vector of signals representing the children.
   *
   * \param n The index of the node whose children are to be retrieved.
   * \return A vector of signals representing the children of the node.
   */
  std::vector<signal_t> const& get_children( node_index_t const& n ) const
  {
    return nodes[n].children;
  }

  /*! \brief Replace a node in the fanin of another node.
   *
   * This method replaces an old node with a new signal in the fanin of a specified node.
   * It updates the fanout count of the new signal and adjusts the outputs accordingly.
   *
   * \param n The index of the node where the replacement occurs.
   * \param old_node The index of the old node to be replaced.
   * \param new_signal The new signal to replace the old node.
   */
  void replace_in_node( node_index_t const& n, node_index_t const& old_node, signal_t new_signal )
  {
    auto& nobj = nodes[n];
    for ( auto& child : nobj.children )
    {
      if ( child.index == old_node )
      {
        child = signal_t{ new_signal.data };
        nodes[new_signal.index].fanout_count++;
        nodes[new_signal.index].outputs[new_signal.output].fanout.push_back( n );
      }
    }
  }

  /*! \brief Replace a node in the outputs of the storage.
   *
   * This method replaces an old node with a new signal in the outputs of the storage.
   * It increments the fanout count of the new signal and updates the output type.
   *
   * \param old_node The old node to be replaced.
   * \param new_signal The new signal to replace the old node.
   */
  void replace_in_outputs( node_index_t const& old_node, signal_t const& new_signal )
  {
    for ( auto& output : outputs )
    {
      if ( output.index == old_node )
      {
        if ( old_node != new_signal.index )
        {
          /* increment fan-in of new node */
          nodes[new_signal.index].fanout_count++;
          nodes[new_signal.index].outputs[new_signal.output].type = pin_type_t::PO;
        }
      }
    }
  }

  /*! \brief Traversal ID for graph algorithms.
   *
   * This ID is used to mark nodes during traversal operations.
   * It is initialized to zero and can be incremented for each traversal.
   */
  uint32_t trav_id = 0u;

  /*! \brief The nodes in the bound network.
   *
   * This vector stores all the nodes in the bound network, each represented
   * by a `storage_node` object. It includes primary inputs, outputs, and
   * internal nodes.
   */
  std::vector<node_t> nodes;

  /*! \brief The primary inputs of the bound network.
   *
   * This vector stores the indices of the primary input nodes in the network.
   * Each input corresponds to a node that can be used as a starting point for
   * logic operations.
   */
  std::vector<node_index_t> inputs;

  /*! \brief The primary outputs of the bound network.
   *
   * This vector stores the signals representing the primary outputs of the network.
   * Each output corresponds to a signal that can be used to observe the results
   * of logic operations in the network.
   */
  std::vector<signal_t> outputs;

  /*! \brief The library of gates used in the bound network.
   *
   * This library contains the gates that can be used to create nodes in the network.
   * It is initialized with a set of gates and provides methods for accessing and
   * manipulating the gates, as well as an AIG list representation for simulation.
   */
  augmented_library<gate> library;
};

} // namespace bound

} // namespace mockturtle
