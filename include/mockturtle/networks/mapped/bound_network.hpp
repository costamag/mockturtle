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
  \file bound.hpp
  \brief Bound network with multiple-output gates support.

  This data structure is a general logic representation which can be used to
  represent a network type with multiple-output cells. Natively, this network
  representation is designed for enabling efficient optimization after technology
  mapping. However, its generality allows us to use it to represent any network
  type in mockturtle, including AIGs, XAIGs, MIGs, XMGs, etc. To support these
  representations, signals can be complemented.

  \author Andrea Costamagna
*/

#pragma once

#include "../../io/genlib_reader.hpp"
#include "../../traits.hpp"
#include "../../utils/algorithm.hpp"
#include "../../utils/index_lists/list_simulator.hpp"
#include "../../utils/mapped/augmented_library.hpp"
#include "../../utils/truth_table_cache.hpp"
#include "../detail/foreach.hpp"
#include "../events.hpp"
#include "../storage.hpp"
#include "bound_storage/bound_node.hpp"
#include "bound_storage/bound_signal.hpp"
#include "bound_storage/bound_storage.hpp"

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>

#include <algorithm>
#include <memory>

namespace mockturtle
{

/*! \brief Network of gates from a technology library.
 *
 * \tparam MaxNumOutputs Maximum number of outputs of the cells in the library.
 */
template<uint32_t MaxNumOutputs = 2u>
class bound_network
{
public:
#pragma region Types and constructors

  static constexpr auto NumBitsOutputs = bound::bits_required<MaxNumOutputs>();

  /* aliases used in this class */
  using storage_t = std::shared_ptr<bound::storage<NumBitsOutputs>>;
  using list_t = large_xag_index_list;
  using node_t = bound::storage_node<NumBitsOutputs>;
  using signal_t = bound::storage_signal<NumBitsOutputs>;
  using node_index_t = bound::node_index_t;
  using hash_t = bound::signal_hash<NumBitsOutputs>;

  /* aliases for compatibility with the other network types */
  static constexpr auto min_fanin_size = 1;
  static constexpr auto max_fanin_size = 32;
  static constexpr auto max_num_outputs = 1u << NumBitsOutputs;
  using base_type = bound_network<MaxNumOutputs>;
  using storage = storage_t;
  using signal = signal_t;
  using node = node_index_t;

  /*! \brief Constructor from a technology library.
   *
   * \param gates The gates in the technology library.
   */
  bound_network( std::vector<gate> const& gates )
      : _storage( std::make_shared<bound::storage<NumBitsOutputs>>( gates ) ),
        _events( std::make_shared<typename decltype( _events )::element_type>() )
  {
  }

  /*! \brief Constructor from a storage object.
   *
   * This constructor is used to create a bound network from an existing storage
   * object, allowing for cloning and manipulation of the network without
   * needing to recreate the storage structure.
   *
   * \param storage The storage object containing the network data.
   */
  bound_network( std::shared_ptr<bound::storage<NumBitsOutputs>> storage )
      : _storage( storage ),
        _events( std::make_shared<typename decltype( _events )::element_type>() )
  {
  }

  /*! \brief Clone the current network.
   *
   * This method creates a new instance of the bound network with a copy of the
   * current storage. It is useful for creating a separate instance of the
   * network that can be modified independently of the original.
   *
   * \return A new bound_network instance with cloned storage.
   */
  bound_network<MaxNumOutputs> clone() const
  {
    return { std::make_shared<bound::storage<NumBitsOutputs>>( *_storage ) };
  }

#pragma endregion

#pragma region Primary I / O and constants
public:
  /*! \brief Returns a constant signal.
   *
   * This method returns a signal representing a constant value (0 or 1).
   * The value can be specified as an argument, with the default being false (0).
   *
   * \param value The constant value to be represented (true for 1, false for 0).
   * \return A signal representing the constant value.
   */
  signal_t get_constant( bool value = false ) const
  {
    return _storage->get_constant( value );
  }

  /*! \brief Creates a primary input signal.
   *
   * \return A signal representing the primary input.
   */
  signal_t create_pi()
  {
    return _storage->create_pi();
  }

  /*! \brief Label a signal as primary output.
   *
   * \param f The signal to be added to the primary outputs.
   * \return A unique identifier for the primary output.
   */
  uint32_t create_po( signal_t const& f )
  {
    return _storage->create_po( f );
  }

  /*! \brief Check if a node is a constant.
   *
   * \param n The node index to check.
   * \return True if the node is a constant, false otherwise.
   */
  bool is_constant( node_index_t const& n ) const
  {
    return _storage->is_constant( n );
  }

  /*! \brief Check if a node is a combinational input (CI).
   *
   * \param n The node index to check.
   * \return True if the node is a combinational input, false otherwise.
   */
  bool is_ci( node_index_t const& n ) const
  {
    return _storage->is_ci( n );
  }

  /*! \brief Check if a node is a primary input (PI).
   *
   * \param n The node index to check.
   * \return True if the node is a primary input, false otherwise.
   */
  bool is_pi( node_index_t const& n ) const
  {
    return _storage->is_pi( n );
  }

  /*! \brief Check if a node is a primary output (PO).
   *
   * \param n The node index to check.
   * \param output The output pin index to check (default is 0).
   * \return True if the node is a primary output, false otherwise.
   */
  bool is_po( node_index_t const& n, uint32_t output = 0 ) const
  {
    return _storage->is_po( n, output );
  }

  /*! \brief Check if a signal is a primary output (PO).
   *
   * \param f The signal to check.
   * \return True if the signal is a primary output, false otherwise.
   */
  bool is_po( signal_t const& f ) const
  {
    return is_po( f.index, f.output );
  }

  /*! \brief Check if a node is a constant 0 or not.
   *
   * \param n The node index to check.
   * \return False if the node represents the constant 0, true otherwise.
   */
  bool constant_value( node_index_t const& n ) const
  {
    return _storage->constant_value( n );
  }
#pragma endregion

#pragma region Create arbitrary functions

  /*! \brief Create a node from the fanin signals and binding IDs.
   *
   * This method creates a new node in the network with the specified children
   * and binding IDs. when more than one binding ID is provided, the node is a
   * multiple-output node, allowing for multiple outputs pins to be associated
   * with different functions.
   *
   * \param children The children signals of the new node.
   * \param ids The binding IDs for the outputs of the new node.
   * \tparam DoStrash If true, the node will be created with strashing enabled.
   * \return A signal representing the newly created node.
   */
  template<bool DoStrash = false>
  signal_t create_node( std::vector<signal_t> const& children,
                        std::vector<uint32_t> const& ids )
  {
    node_t const& n = _storage->create_storage_node( children, ids );

    /* structural hashing */
    if constexpr ( DoStrash )
    {
      const auto it = _storage->find( n );
      if ( it )
      {
        if ( !is_dead( *it ) )
          ;
        {
          return { *it, 0 };
        }
      }
    }

    signal_t const f = _storage->create_node( children, n );

    /* initialize the application specific value to 0 */
    set_value( f.index, 0 );

    for ( auto const& fn : _events->on_add )
    {
      ( *fn )( f.index );
    }
    return f;
  }

  /*! \brief Create a node with a single binding ID.
   *
   * This method is a convenience overload for creating a node with a single
   * binding ID. It calls the more general `create_node` method with a vector
   * containing the single ID.
   *
   * \param children The children signals of the new node.
   * \param id The binding ID for the output of the new node.
   * \tparam DoStrash If true, the node will be created with strashing enabled.
   * \return A signal representing the newly created node.
   */
  template<bool DoStrash = false>
  signal_t create_node( std::vector<signal_t> const& children,
                        uint32_t id )
  {
    return create_node<DoStrash>( children, std::vector<uint32_t>{ id } );
  }

  /*! \brief Clone a node from another bound network.
   *
   * This method creates a new node in the current network by cloning an existing
   * node from another bound network. It takes the source node and its children
   * signals, and creates a new node with the same binding IDs.
   *
   * \param other The other bound network from which to clone the node.
   * \param source The index of the source node to clone.
   * \param children The children signals of the new node.
   * \return A signal representing the newly cloned node.
   */

  template<bool DoStrash = false>
  signal_t clone_node( bound_network const& other,
                       node_index_t const& source,
                       std::vector<signal_t> const& children )
  {
    assert( !children.empty() );
    std::vector<uint32_t> const ids = other.get_binding_ids( source );
    return create_node<DoStrash>( children, ids );
  }

#pragma endregion

#pragma region Restructuring

  /*! \brief Substitute a node with signals equivalent to its output pins.
   *
   * This method replaces a node's output pins with functionally equivalent
   * signals. It updates the fanout of the old node's outputs to point to the
   * new signals, effectively substituting the old node in the network.
   *
   * [ pin 0 ] -> [ new_signals[0] ]
   * ...
   * [ pin j ] -> [ new_signals[j] ]
   *
   * \param old_node The index of the old node to be replaced.
   * \param new_signals The new signals to replace the old node's outputs.
   */
  void substitute_node( node_index_t const& old_node,
                        std::vector<signal_t> const& new_signals )
  {
    /* update the signals to be used as primary outputs.
     * Highest priority so that on_modified events operate on the correct POs.
     */
    _storage->replace_in_outputs( old_node, new_signals );

    /* update the fanins/fanout information and trigger modified events */
    replace_in_node( old_node, new_signals );

    /* remove the node and trigger on_delete events */
    take_out_node( old_node, new_signals );
  }

  /*! \brief Substitute a node with a new signal.
   *
   * This method replaces an old node with a new signal in the network, falling
   * back to the general case of vectorized substitute node.
   *
   * \param old_node The index of the old node to be replaced.
   * \param new_signal The new signal to replace the old node's outputs.
   */
  void substitute_node( node_index_t const& old_node, signal_t const& new_signal )
  {
    /* fall back to the general case */
    substitute_node( old_node, std::vector<signal_t>{ new_signal } );
  }

  /*! \brief Update the fanin-fanout information in the network.
   *
   * Iterate over all output pins of the node to be removed and replace
   * the old node's outputs with the new signals. This method updates the
   * fanout count of the new signals and adjusts its fanout accordingly.
   *
   * \param old_node The index of the old node to be replaced.
   * \param new_signals The new signals to replace the old node's outputs.
   */
  void replace_in_node( node_index_t const& old_node,
                        std::vector<signal_t> const& new_signals )
  {
    assert( num_outputs( old_node ) == new_signals.size() &&
            "Number of new signals must match the number of outputs" );

    /* iterate over all output pins of the node to be removed */
    _storage->foreach_output_pin( old_node, [&]( auto const& pin, auto i ) {
      signal_t const old_signal = signal_t{ old_node, i };
      /* replace the old signal in the fanout of the output pin */
      _storage->foreach_fanout( pin, [&]( auto const& fanout_node, auto j ) {
        (void)j; // unused variable
        /* replace the old signal with the new signal in the fanout */
        replace_in_node( fanout_node, old_signal, new_signals[i] );
      } );
    } );
  }

  /*! \brief Replace a node in the fanin of another node.
   *
   * This method replaces an old node with a new signal in the fanin of a specified node.
   * It updates the fanout count of the new signal and adjusts the outputs accordingly.
   *
   * \param root The index of the node where the replacement occurs.
   * \param old_signal The old signal to be replaced.
   * \param new_signal The new signal to replace the old node.
   */
  void replace_in_node( node_index_t const& root,
                        signal_t const& old_signal,
                        signal_t new_signal )
  {
    node_index_t const old_node = old_signal.index;

    if ( !_storage->in_fanin( root, old_node ) )
      return;

    auto const old_children = _storage->get_children( root );

    _storage->update_nets( root, old_signal, new_signal );

    /* provide the root and the old signals to the modified event.
     * This corresponds to all the information, since the new children can be
     * computed from the root.
     */
    for ( auto const& fn : _events->on_modified )
    {
      ( *fn )( root, old_children );
    }
  }

  /*! \brief Take out a node if it is not reused in the new nodes.
   *
   * This method checks if the old node is still used in the network.
   * If it is not, it removes the node and updates the fanout counts of its children.
   *
   * \param old_node The index of the old node to be removed.
   * \param new_signals The new signals to replace the old node's outputs.
   */
  void take_out_node( node_index_t const& old_node,
                      std::vector<signal_t> const& new_signals )
  {
    /* take out the node if it is not reused in the new nodes */
    for ( auto const& f : new_signals )
    {
      if ( f.index == old_node )
      {
        /* if the old node is still used, we cannot take it out */
        return;
      }
    }
    take_out_node( old_node );
  }

  /*! \brief Take out a node from the network.
   *
   * This method removes a node from the network, marking it as dead and
   * updating the fanout counts of its children. It also triggers events
   * for the deletion of the node.
   *
   * \param n The index of the node to be removed.
   */
  void take_out_node( node_index_t const& n )
  {
    /* we cannot delete CIs, constants, or already dead nodes */
    if ( is_constant( n ) || is_ci( n ) || is_dead( n ) )
      return;

    auto children = _storage->get_children( n );

    /* NOTE: the node's information is not cleared-up yet, so we can
     * access the node's outputs or the node's fanins. Not its old fanouts.
     */
    for ( auto const& fn : _events->on_delete )
    {
      ( *fn )( n );
    }

    /* mark the node as dead */
    _storage->delete_node( n );

    /* if the node has been deleted, then deref fanout_size of
       fanins and try to take them out if their fanout_size become 0 */
    for ( auto i = 0; i < children.size(); ++i )
    {
      auto& child = children[i];
      _storage->delete_fanout( child, n );
      if ( fanout_size( child.index ) == 0 )
      {
        take_out_node( child.index );
      }
    }
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
    return _storage->is_multioutput( n );
  }

  bool is_multioutput( std::string const& name ) const
  {
    return _storage->is_multioutput( name );
  }

  inline bool is_dead( node_index_t const& n ) const
  {
    return _storage->is_dead( n );
  }

  auto size() const
  {
    return _storage->size();
  }

  auto signal_size() const
  {
    return max_num_outputs * size();
  }

  auto num_cis() const
  {
    return _storage->num_cis();
  }

  auto num_cos() const
  {
    return _storage->num_cos();
  }

  auto num_pis() const
  {
    return _storage->num_pis();
  }

  auto num_pos() const
  {
    return _storage->num_pos();
  }

  auto num_gates() const
  {
    return _storage->num_gates();
  }

  uint32_t num_outputs( node_index_t const& n ) const
  {
    return _storage->num_outputs( n );
  }

  uint32_t fanin_size( node_index_t const& n ) const
  {
    return _storage->fanin_size( n );
  }

  uint32_t fanout_size( node_index_t const& n ) const
  {
    return _storage->fanout_size( n );
  }

  uint32_t incr_fanout_size( node_index_t const& n ) const
  {
    return _storage->incr_fanout_size( n );
  }

  uint32_t decr_fanout_size( node_index_t const& n ) const
  {
    return _storage->decr_fanout_size( n );
  }

  uint32_t incr_fanout_size_pin( node_index_t const& n, uint32_t pin_index ) const
  {
    return _storage->incr_fanout_size_pin( n, pin_index );
  }

  uint32_t decr_fanout_size_pin( node_index_t const& n, uint32_t pin_index ) const
  {
    return _storage->decr_fanout_size_pin( n, pin_index );
  }

  uint32_t fanout_size_pin( node_index_t const& n, uint32_t pin_index ) const
  {
    return _storage->fanout_size_pin( n, pin_index );
  }

  bool is_function( node_index_t const& n ) const
  {
    return _storage->is_function( n );
  }
#pragma endregion

#pragma region Functional properties
  kitty::dynamic_truth_table signal_function( const signal_t& f ) const
  {
    return _storage->signal_function( f );
  }

  kitty::dynamic_truth_table node_function( const node_index_t& n, uint32_t pin_index = 0 ) const
  {
    signal f = make_signal( n, pin_index );
    return signal_function( f );
  }
#pragma endregion

#pragma region Nodes and signals
  node_index_t get_node( signal_t const& f ) const
  {
    return f.index;
  }

  signal_t make_signal( node_index_t const& n, uint32_t output_pin ) const
  {
    return signal_t{ n, output_pin };
  }

  signal_t make_signal( node_index_t const& n ) const
  {
    return make_signal( n, 0 );
  }

  bool is_complemented( signal_t const& f ) const
  {
    (void)f;
    return false;
  }

  uint32_t get_output_pin( signal_t const& f ) const
  {
    return static_cast<uint32_t>( f.output );
  }

  signal_t next_output_pin( signal_t const& f ) const
  {
    return signal_t{ f.index, f.output + 1 };
  }

  uint32_t node_to_index( node_index_t const& n ) const
  {
    return static_cast<uint32_t>( n );
  }

  node_index_t index_to_node( uint32_t index ) const
  {
    return index;
  }

  uint64_t signal_to_index( signal_t const& f ) const
  {
    return static_cast<uint32_t>( f.data );
  }

  node_index_t ci_at( uint32_t index ) const
  {
    return _storage->ci_at( index );
  }

  signal_t co_at( uint32_t index ) const
  {
    return _storage->co_at( index );
  }

  node_index_t pi_at( uint32_t index ) const
  {
    return _storage->pi_at( index );
  }

  signal_t po_at( uint32_t index ) const
  {
    return _storage->po_at( index );
  }

  uint32_t pi_index( node_index_t const& n ) const
  {
    return _storage->pi_index( n );
  }

  uint32_t po_index( signal_t const& f ) const
  {
    return _storage->po_index( f );
  }
#pragma endregion

#pragma region Node and signal iterators
  template<typename Fn>
  void foreach_node( Fn&& fn ) const
  {
    _storage->foreach_node( fn );
  }

  template<typename Fn>
  void foreach_ci( Fn&& fn ) const
  {
    _storage->foreach_ci( fn );
  }

  template<typename Fn>
  void foreach_co( Fn&& fn ) const
  {
    _storage->foreach_co( fn );
  }

  template<typename Fn>
  void foreach_pi( Fn&& fn ) const
  {
    _storage->foreach_pi( fn );
  }

  template<typename Fn>
  void foreach_po( Fn&& fn ) const
  {
    _storage->foreach_po( fn );
  }

  template<typename Fn>
  void foreach_gate( Fn&& fn ) const
  {
    _storage->foreach_gate( fn );
  }

  template<typename Fn>
  void foreach_fanin( node_index_t const& n, Fn&& fn ) const
  {
    _storage->foreach_fanin( n, fn );
  }

  template<typename Fn>
  void foreach_fanout( node_index_t const& n, Fn&& fn ) const
  {
    _storage->foreach_fanout( n, fn );
  }

  template<typename Fn>
  void foreach_fanout( signal_t const& f, Fn&& fn ) const
  {
    _storage->foreach_fanout( f, fn );
  }

  template<typename Fn>
  void foreach_tfo_node( node_index_t const& n, Fn&& fn ) const
  {
    _storage->foreach_tfo_node( n, fn );
  }

  template<typename Fn>
  void foreach_output_pin( node_index_t const& n, Fn&& fn ) const
  {
    _storage->foreach_output_pin( n, fn );
  }

  template<typename Fn>
  void foreach_output( node_index_t const& n, Fn&& fn ) const
  {
    _storage->foreach_output( n, fn );
  }
#pragma endregion

#pragma region Simulate values

  /*! \brief Get the cached simulator for AIG index lists.
   *
   * Caching an unique simulator avoids reallocations of different simulation
   * engines, ensuring memory efficiency.
   */
  template<typename TT>
  std::shared_ptr<list_simulator<list_t, TT>> get_simulator() const
  {
    using simulator_t = list_simulator<list_t, TT>;
    static const std::shared_ptr<simulator_t> sim = std::make_shared<simulator_t>();
    return sim;
  }

  /*! \brief Simulation of the input patterns using the node's function.
   *
   * \param n index of the node to simulate
   * \param sim_ptrs vector of pointers to the simulation of the fanins.
   * \return A vector of truth-tables, one for each output pin of the node.
   */
  template<typename TT>
  std::vector<TT> compute( node_index_t const& n, std::vector<TT const*> sim_ptrs ) const
  {
    std::vector<TT> res;
    compute( res, n, sim_ptrs );
    return res;
  }

  /*! \brief Inline simulation of the input patterns using the node's function.
   *
   * \param n index of the node to simulate
   * \param sim_ptrs vector of pointers to the simulation of the fanins.
   */
  template<typename TT>
  void compute( std::vector<TT>& res, node_index_t const& n, std::vector<TT const*> sim_ptrs ) const
  {
    auto simulator_ptr = get_simulator<TT>();
    res.resize( num_outputs( n ) );
    const auto nfanin = fanin_size( n );
    assert( nfanin > 0 );
    assert( sim_ptrs.size() == nfanin );

    _storage->foreach_output_pin( n, [&]( auto const& pin, auto i ) {
      auto id = pin.id;
      auto const& list = _storage->get_list( id );
      ( *simulator_ptr )( list, sim_ptrs );
      simulator_ptr->get_simulation_inline( res[i], list, sim_ptrs, list.po_at( 0 ) );
    } );
  }

  /*! \brief Inline simulation of the input patterns using the node's function.
   *
   * \param n index of the node to simulate
   * \param sim_ptrs vector of pointers to the simulation of the fanins.
   */
  template<typename TT>
  void compute( TT& res, signal_t const& f, std::vector<TT const*> sim_ptrs ) const
  {
    auto simulator_ptr = get_simulator<TT>();
    const auto nfanin = fanin_size( get_node( f ) );
    assert( nfanin > 0 );
    assert( sim_ptrs.size() == nfanin );
    auto const& g = get_binding( f );
    auto const& list = _storage->get_list( g.id );
    ( *simulator_ptr )( list, sim_ptrs );
    simulator_ptr->get_simulation_inline( res, list, sim_ptrs, list.po_at( 0 ) );
  }
#pragma endregion

#pragma region Custom node values
  void clear_values() const
  {
    _storage->clear_values();
  }

  uint32_t value( node_index_t const& n ) const
  {
    return _storage->value( n );
  }

  void set_value( node_index_t const& n, uint32_t v ) const
  {
    _storage->set_value( n, v );
  }

  uint32_t incr_value( node_index_t const& n ) const
  {
    return _storage->incr_value( n );
  }

  uint32_t decr_value( node_index_t const& n ) const
  {
    return _storage->decr_value( n );
  }
#pragma endregion

#pragma region Visited flags
  void clear_visited() const
  {
    _storage->clear_visited();
  }

  auto visited( node_index_t const& n ) const
  {
    return _storage->visited( n );
  }

  void set_visited( node_index_t const& n, uint32_t v ) const
  {
    _storage->nodes[n].traversal_id = v;
  }

  uint32_t trav_id() const
  {
    return _storage->get_trav_id();
  }

  void incr_trav_id() const
  {
    _storage->incr_trav_id();
  }
#pragma endregion

#pragma region Getters
  std::vector<signal_t> const& get_children( node_index_t const& n ) const
  {
    return _storage->get_children( n );
  }
#pragma endregion

#pragma region General methods
  auto& events() const
  {
    return *_events;
  }
#pragma endregion

#pragma region Binding
  std::vector<uint32_t> get_binding_ids( node_index_t const& n ) const
  {
    return _storage->get_binding_ids( n );
  }

  auto const& get_binding( signal_t const& f ) const
  {
    return _storage->get_binding( f );
  }

  auto const& has_binding( node_index_t const& n ) const
  {
    return _storage->has_binding( n );
  }
#pragma endregion

public:
  std::shared_ptr<bound::storage<NumBitsOutputs>> _storage;
  std::shared_ptr<network_events<base_type>> _events;
};

} // namespace mockturtle
