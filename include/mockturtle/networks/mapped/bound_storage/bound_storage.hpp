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
#include "../../../utils/mapped/augmented_library.hpp"
#include "bound_node.hpp"
#include "bound_signal.hpp"
#include "bound_utils.hpp"

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
#pragma region Types and constructors

public:
  using gate_t = mockturtle::gate;
  using list_t = large_xag_index_list;
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
  storage( std::vector<gate_t> const& gates )
      : library( gates )
  {
    /* reserve space for nodes */
    nodes.reserve( 10000u ); // reserve to avoid frequent reallocations

    /* we reserve the first two nodes for constants */
    nodes.emplace_back( pin_type_t::CONSTANT ); // 0
    nodes.emplace_back( pin_type_t::CONSTANT ); // 1
  }

#pragma endregion

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
   * A PI stores its index in the only fanin it has.
   *
   * \return A signal representing the primary input.
   */
  signal_t create_pi()
  {
    const auto index = nodes.size();
    node_t input( pin_type_t::PI );
    input.children = { inputs.size() };
    nodes.emplace_back( input );
    inputs.emplace_back( index );
    return signal_t{ index, 0 };
  }

  /*! \brief Creates a primary output signal.
   *
   * This method creates a primary output signal from a given signal.
   * It increases the reference count for the node to avoid incorrect deletions
   * and updates the output pin type to indicate it is a primary output. A node
   * can be used as PO more than once, so the number of fanouts is equal to the
   * fanout size of the output pins, plus the number of times one of its pins
   * is used as a primary output.
   * It returns the index of the newly created primary output.
   *
   * \param f The signal representing the function to be used as primary output.
   * \return The index of the newly created primary output.
   */
  uint32_t create_po( signal_t const& f )
  {
    /* increase ref-count to children */
    nodes[f.index].fanout_count++;
    nodes[f.index].outputs[f.output].type |= pin_type_t::PO;
    auto const po_index = static_cast<uint32_t>( outputs.size() );
    outputs.emplace_back( f.index, f.output );
    return po_index;
  }

  /*! \brief Check if the node is a constant */
  bool is_constant( node_index_t const& n ) const
  {
    auto const& pins = nodes[n].outputs;
    return has_intersection( pins[0].type, pin_type_t::CONSTANT );
  }

  /*! \brief Check if the node is a combinational input (CI)
   *
   * \param n The node to check.
   * \return True if the node is a combinational input, false otherwise.
   */
  bool is_ci( node_index_t const& n ) const
  {
    auto const& pins = nodes[n].outputs;
    return has_intersection( pins[0].type, pin_type_t::PI ) ||
           has_intersection( pins[0].type, pin_type_t::CI );
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
    auto& pins = nodes[n].outputs;
    return has_intersection( pins[output].type, pin_type_t::PO );
  }

  bool is_po( signal_t f ) const
  {
    return is_po( f.index, f.output );
  }

  /*! \brief Returns if node index is not 0 */
  bool constant_value( node_index_t const& n ) const
  {
    return n != 0;
  }

#pragma endregion

#pragma region Create arbitrary functions

  /*! \brief Create a detailed node to be stored.
   *
   * \param children input signals.
   * \param ids binding identifier of the output pins.
   */
  node_t create_storage_node( std::vector<signal_t> const& children, std::vector<uint32_t> const& ids )
  {
    for ( auto i = 1; i < ids.size(); ++i )
    {
      assert( library.get_name( ids[i] ) == library.get_name( ids[0] ) &&
              "Multiple-output nodes are expected to have the same name" );
    }

    node_t new_node;
    std::copy( children.begin(), children.end(), std::back_inserter( new_node.children ) );

    new_node.outputs = decltype( new_node.outputs )( ids.size() );

    for ( auto i = 0; i < ids.size(); ++i )
      new_node.outputs[i] = { ids[i], pin_type_t::INTERNAL };

    return new_node;
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
  signal_t create_node( std::vector<signal_t> const& children, node_t const& n )
  {
    const auto index = nodes.size();
    nodes.push_back( n );

    /* increase ref-count to children */
    for ( auto c : children )
    {
      nodes[c.index].fanout_count++;
      nodes[c.index].outputs[c.output].fanout_count++;
      nodes[c.index].outputs[c.output].fanout.push_back( index );
    }

    if ( hash.find( n ) != hash.end() )
    {
      hash[n].push_back( index );
    }
    else
    {
      hash[n] = { index };
    }

    return { index, 0 };
  }

#pragma endregion

#pragma region Restructuring

  /*! \brief Update the list of POs when a signal in the old list is replaced.
   *
   * \param old_node The index of the old node to be replaced.
   * \param new_signals The new signals to replace the old node's outputs.
   */
  void replace_in_outputs( node_index_t const& old_node,
                           std::vector<signal_t> const& new_signals )
  {
    foreach_output_pin( old_node, [&]( auto const& pin, uint32_t i ) {
      signal_t const old_signal = signal_t{ old_node, i };
      if ( is_po( old_signal ) )
      {
        /* replace the output signals with the new signal */
        replace_output( old_signal, new_signals[i] );
      }
    } );
  }

  /*! \brief Replace an output signal in the outputs list.
   *
   * This method replaces an old signal with a new signal in the outputs list.
   * It updates the fanout count of the new signal and the old signal,
   *
   * \param old_signal The old signal to be replaced.
   * \param new_signal The new signal to replace the old signal.
   */
  void replace_output( signal_t const& old_signal,
                       signal_t const& new_signal )
  {
    bool found = false;
    for ( auto& output : outputs )
    {
      if ( output == old_signal )
      {
        found = true;
        /* replace the old signal with the new signal in the outputs */
        output = new_signal;
        /* update the fanout count of the new signal */
        nodes[new_signal.index].fanout_count++;
        nodes[old_signal.index].fanout_count--;
        nodes[old_signal.index].outputs[old_signal.output].type &= ~pin_type_t::PO;
        nodes[new_signal.index].outputs[new_signal.output].type |= pin_type_t::PO;
      }
    }
    assert( found && "Output signal not found in the outputs list" );
  }

  /*! \brief Insert a fanout node for a signal.
   *
   * This method inserts a fanout for a given signal and node index.
   * It updates the fanout count of the signal and the output pin accordingly.
   *
   * \param f The signal for which the fanout is to be inserted.
   * \param n The node index to be added as a fanout.
   */
  void insert_fanout( signal_t const& f, node_index_t const& n )
  {
    auto& fanout = nodes[f.index].outputs[f.output].fanout;
    uint32_t const occurrences = std::count( fanout.begin(),
                                             fanout.end(),
                                             n );
    if ( occurrences > 0 )
    {
      /* if the fanout already exists, we do not need to insert it again */
      return;
    }
    nodes[f.index].fanout_count++;
    nodes[f.index].outputs[f.output].fanout_count++;
    nodes[f.index].outputs[f.output].fanout.push_back( n );
  }

  /*! \brief Delete a fanout node for a signal.
   *
   * This method deletes a fanout for a given signal and node index.
   * It updates the fanout count of the signal and the output pin accordingly.
   *
   * \param f The signal for which the fanout is to be deleted.
   * \param n The node index to be removed from the fanout.
   */
  void delete_fanout( signal_t const& f, node_index_t const& n )
  {
    auto& fanout = nodes[f.index].outputs[f.output].fanout;
    uint32_t const occurrences = std::count( fanout.begin(),
                                             fanout.end(),
                                             n );
    nodes[f.index].fanout_count -= occurrences;
    nodes[f.index].outputs[f.output].fanout_count -= occurrences;
    fanout.erase( std::remove( fanout.begin(),
                               fanout.end(),
                               n ),
                  fanout.end() );
  }

  /*! \brief Update the interconnections of a node.
   *
   * This method updates the fanin-fanout information in the network by replacing an old signal
   * with a new signal in the fanin of a specified node. It updates the fanout count of the new signal
   * and applies events to notify about the modification.
   *
   * \param root The index of the node where the replacement occurs.
   * \param old_signal The old signal to be replaced.
   * \param new_signal The new signal to replace the old node.
   */
  void update_nets( node_index_t const& root,
                    signal_t const& old_signal,
                    signal_t new_signal )
  {
    auto& nd_root = nodes[root];
    for ( auto& child : nd_root.children )
    {
      if ( child == old_signal )
      {
        /* add the root node to the new signal's fanout */
        insert_fanout( new_signal, root );

        /* update the fanout of the old signal */
        delete_fanout( child, root );

        /* replace the old signal with the new signal in the children */
        child = new_signal;
      }
    }
  }

  /*! \brief Delete a node from the network.
   *
   * This method removes a node from the network, marking it as dead and
   * updating the fanout counts of its children. It also triggers events
   * for the deletion of the node.
   *
   * \param n The index of the node to be removed.
   */
  void delete_node( node_index_t const& n )
  {
    /* remove the node from the hash table if present */
    auto it = hash.find( nodes[n] );
    if ( it != hash.end() )
    {
      auto& list = it->second;
      std::remove( list.begin(), list.end(), n );
      if ( list.empty() )
      {
        hash.erase( it );
      }
    }
    /* mark the node as dead */
    for ( auto& pin : nodes[n].outputs )
    {
      pin.type |= pin_type_t::DEAD;
      pin.fanout.clear();
      pin.fanout_count = 0;
    }
    nodes[n].fanout_count = 0;
    nodes[n].children.clear();
  }

#pragma endregion

#pragma region Structural properties

  /*! \brief Returns true since the network is combinational.
   *
   * TODO: Add support for sequential elements in the future.
   * \return Always returns true.
   */
  bool is_combinational() const
  {
    return true;
  }

  bool is_multioutput( node_index_t const& n ) const
  {
    return num_outputs( n ) > 1;
  }

  bool is_multioutput( std::string const& name ) const
  {
    return library.is_multioutput( name );
  }

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
      bool const dead = has_intersection( pin.type, pin_type_t::DEAD );
      all_dead &= dead;
      one_dead |= dead;
    }
    assert( !( all_dead ^ one_dead ) );
    /* A dead node is simply a dangling node */
    return all_dead;
  }

  inline bool is_constant( signal_t const& f ) const
  {
    auto const& outputs = nodes[f.index].outputs;
    uint32_t i = f.output;
    return outputs.size() > 0 && has_intersection( outputs[i].type, pin_type_t::CONSTANT );
  }

  auto size() const
  {
    return static_cast<uint32_t>( nodes.size() );
  }

  auto num_cis() const
  {
    return static_cast<uint32_t>( inputs.size() );
  }

  auto num_cos() const
  {
    return static_cast<uint32_t>( outputs.size() );
  }

  auto num_pis() const
  {
    return static_cast<uint32_t>( inputs.size() );
  }

  auto num_pos() const
  {
    return static_cast<uint32_t>( outputs.size() );
  }

  auto num_gates() const
  {
    return static_cast<uint32_t>( nodes.size() - inputs.size() - 2 );
  }

  uint32_t num_outputs( node_index_t const& n ) const
  {
    return static_cast<uint32_t>( nodes[n].outputs.size() );
  }

  uint32_t fanin_size( node_index_t const& n ) const
  {
    return static_cast<uint32_t>( nodes[n].children.size() );
  }

  uint32_t fanout_size( node_index_t const& n ) const
  {
    return nodes[n].fanout_count;
  }

  uint32_t incr_fanout_size( node_index_t const& n ) const
  {
    return nodes[n].fanout_count++;
  }

  uint32_t decr_fanout_size( node_index_t const& n ) const
  {
    return --nodes[n].fanout_count;
  }

  uint32_t incr_fanout_size_pin( node_index_t const& n, uint32_t pin_index ) const
  {
    return ++nodes[n].outputs[pin_index].fanout_count;
  }

  uint32_t decr_fanout_size_pin( node_index_t const& n, uint32_t pin_index ) const
  {
    return --nodes[n].outputs[pin_index].fanout_count;
  }

  uint32_t fanout_size_pin( node_index_t const& n, uint32_t pin_index ) const
  {
    return nodes[n].outputs[pin_index].fanout_count;
  }

  bool is_function( node_index_t const& n ) const
  {
    auto const& outputs = nodes[n].outputs;
    return ( outputs.size() > 0 ) && ( has_intersection( outputs[0].type, bound::pin_type_t::INTERNAL ) ||
                                       has_intersection( outputs[0].type, bound::pin_type_t::PO ) );
  }

  /*! \brief Checks if a given node is already present in the storage.
   *
   * Uses structural hashing to check if a node is in the network.
   */
  std::optional<node_index_t> find( node_t const& n ) const
  {
    auto it = hash.find( n );
    if ( it != hash.end() )
    {
      assert( !is_dead( it->second[0] ) &&
              "The node should not be dead when looking for it" );
      return it->second[0];
    }
    return std::nullopt;
  }

  /*! \brief Checks if a node is in the fanin of another one.
   */
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

#pragma endregion

#pragma region Functional properties
  kitty::dynamic_truth_table signal_function( const signal_t& f ) const
  {
    auto const& outputs = nodes[f.index].outputs;
    auto const& id = outputs[f.output].id;
    return library[id].function;
  }
#pragma endregion

#pragma region Nodes and signals

  node_index_t ci_at( uint32_t index ) const
  {
    assert( index < inputs.size() );
    return *( inputs.begin() + index );
  }

  signal_t co_at( uint32_t index ) const
  {
    assert( index < outputs.size() );
    return *( outputs.begin() + index );
  }

  node_index_t pi_at( uint32_t index ) const
  {
    assert( index < inputs.size() );
    return *( inputs.begin() + index );
  }

  signal_t po_at( uint32_t index ) const
  {
    assert( index < outputs.size() );
    return *( outputs.begin() + index );
  }

  uint32_t pi_index( node_index_t const& n ) const
  {
    auto const& outputs = nodes[n].outputs;
    assert( has_intersection( outputs[0].type, bound::pin_type_t::PI ) );
    return static_cast<uint32_t>( nodes[n].children[0].data );
  }
#pragma endregion

#pragma region Node and signal iterators

  template<typename Fn>
  void foreach_node( Fn&& fn ) const
  {
    for ( node_index_t n = 2u; n < nodes.size(); ++n )
    {
      if ( !is_dead( n ) )
      {
        fn( n );
      }
    }
  }

  template<typename Fn>
  void foreach_ci( Fn&& fn ) const
  {
    detail::foreach_element( inputs.begin(), inputs.end(), fn );
  }

  template<typename Fn>
  void foreach_co( Fn&& fn ) const
  {
    using IteratorType = decltype( outputs.begin() );
    detail::foreach_element<IteratorType>(
        outputs.begin(), outputs.end(), []( auto f ) { return signal_t( f ); }, fn );
  }

  template<typename Fn>
  void foreach_pi( Fn&& fn ) const
  {
    detail::foreach_element( inputs.begin(), inputs.end(), fn );
  }

  template<typename Fn>
  void foreach_po( Fn&& fn ) const
  {
    for ( signal_t const& f : outputs )
    {
      fn( f );
    }
  }

  template<typename Fn>
  void foreach_gate( Fn&& fn ) const
  {
    auto r = range<uint64_t>( 2u, nodes.size() ); /* start from 2 to avoid constants */
    detail::foreach_element_if(
        r.begin(), r.end(),
        [this]( auto n ) { return !is_ci( n ) && !is_dead( n ); },
        fn );
  }

  template<typename Fn>
  void foreach_fanin( node_index_t const& n, Fn&& fn ) const
  {
    if ( is_constant( n ) || is_ci( n ) || is_pi( n ) )
      return;

    auto const& children = nodes[n].children;
    for ( auto i = 0u; i < children.size(); ++i )
    {
      fn( children[i], i );
    }
  }

  /*! \brief Iterate over the output pins of a node.
   *
   * This method iterates over the output pins of a specified node and applies
   * a function to each pin. The function receives the pin and its index as arguments.
   *
   * \param n The index of the node whose output pins are to be iterated.
   * \param fn The function to apply to each output pin.
   */
  template<typename Fn>
  void foreach_output_pin( node_index_t const& n, Fn&& fn ) const
  {
    auto& pins = nodes[n].outputs;
    for ( auto i = 0u; i < pins.size(); ++i )
    {
      fn( pins[i], i );
    }
  }

  /*! \brief Iterate over the outputs of a node.
   *
   * This method iterates over the output pins of a specified node represented
   * as signals, applying a function to each signal.
   *
   * \param n The index of the node whose output pins are to be iterated.
   * \param fn The function to apply to each output pin.
   */
  template<typename Fn>
  void foreach_output( node_index_t const& n, Fn&& fn ) const
  {
    auto& pins = nodes[n].outputs;
    for ( auto i = 0u; i < pins.size(); ++i )
    {
      auto const f = signal_t{ n, i };
      fn( f );
    }
  }

  template<typename Fn>
  void foreach_fanout( output_pin_t const& pin, Fn&& fn ) const
  {
    auto& fanout = pin.fanout;
    for ( uint32_t i = 0; i < fanout.size(); ++i )
    {
      node_index_t const& n = fanout[i];
      fn( n, i );
    }
  }

  template<typename Fn>
  void foreach_fanout( signal_t const& f, Fn&& fn ) const
  {
    auto& fanout = nodes[f.index].outputs[f.output].fanout;
    for ( uint32_t i = 0; i < fanout.size(); ++i )
    {
      node_index_t const& n = fanout[i];
      fn( n );
    }
  }

  template<typename Fn>
  void foreach_fanout( node_index_t const& n, Fn&& fn ) const
  {
    foreach_output_pin( n, [&]( auto const& pin, auto i ) {
      foreach_fanout( pin, [&]( auto const& fanout_node, auto j ) {
        fn( fanout_node );
      } );
    } );
  }

#pragma endregion

#pragma region Custom node values
  void clear_values()
  {
    std::for_each( nodes.begin(), nodes.end(), []( auto& n ) { n.user_data = 0; } );
  }

  uint32_t value( node_index_t const& n ) const
  {
    return nodes[n].user_data;
  }

  void set_value( node_index_t const& n, uint32_t v )
  {
    nodes[n].user_data = v;
  }

  uint32_t incr_value( node_index_t const& n )
  {
    return static_cast<uint32_t>( nodes[n].user_data++ );
  }

  uint32_t decr_value( node_index_t const& n )
  {
    return static_cast<uint32_t>( --nodes[n].user_data );
  }
#pragma endregion

#pragma region Visited flags
  void clear_visited()
  {
    std::for_each( nodes.begin(), nodes.end(), []( auto& n ) { n.traversal_id = 0; } );
  }

  auto visited( node_index_t const& n ) const
  {
    return nodes[n].traversal_id;
  }

  void set_visited( node_index_t const& n, uint32_t v )
  {
    nodes[n].traversal_id = v;
  }

  uint32_t get_trav_id() const
  {
    return trav_id;
  }

  void incr_trav_id()
  {
    if ( trav_id > ( std::numeric_limits<uint32_t>::max() - 10 ) )
    {
      std::cout << "[w] Traversal identifier exceeded safe treshold. Forced reset" << std::endl;
      clear_values();
      clear_visited();
      trav_id = 0;
    }
    ++trav_id;
  }
#pragma endregion

#pragma region Getters
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

  /*! \brief
   */
  list_t const& get_list( uint32_t id ) const
  {
    return library.get_list( id );
  }

  /*! \brief Get the binding identifiers of the output pins in a node
   */
  std::vector<uint32_t> get_binding_ids( node_index_t const& n ) const
  {
    std::vector<uint32_t> ids;
    for ( auto const& pin : nodes[n].outputs )
      ids.push_back( pin.id );
    return ids;
  }

  auto const& get_binding( signal_t const& f ) const
  {
    auto const& pin = nodes[f.index].outputs[f.output];
    return library.get_gate( pin.id );
  }

  bool has_binding( node_index_t const& n ) const
  {
    return std::is_same<gate_t, mockturtle::gate>::value;
  }
#pragma endregion

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
  augmented_library<gate_t> library;

  /*! \brief Hash map for fast node lookups.
   *
   * This hash map allows for quick access to nodes based on their indices.
   * It uses a custom hash function to ensure efficient storage and retrieval
   * of nodes in the network.
   */
  phmap::flat_hash_map<node_t, std::vector<node_index_t>, bound::node_hash<NumBitsOutputs>> hash;
};

} // namespace bound

} // namespace mockturtle
