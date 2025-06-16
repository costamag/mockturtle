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
  \brief Bound network for standard cell design with multiple-output support

  Similarly to the `block_network`, this data structure is designed to support
  mapping with multiple-output gates, but it introduces the following features:
  - Two nodes might have the same functionality, but different binding id. In
    traditional technology mappers group cells with the same functionality into
    equivalence classes. Supporting diversity across them allows us to consider
    load capacitance and sizing.
  - Each gate is combined with a Boolean chain for efficient Boolean evaluation.

  \author Andrea Costamagna
*/

#pragma once

#include "../io/genlib_reader.hpp"
#include "../traits.hpp"
#include "../utils/algorithm.hpp"
#include "../utils/index_lists/list_simulator.hpp"
#include "../utils/mapping/augmented_library.hpp"
#include "../utils/truth_table_cache.hpp"
#include "detail/foreach.hpp"
#include "events.hpp"
#include "storage.hpp"
#include "storage/bound_storage.hpp"

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

  /* aliases for compatibility with the other network types */
  using base_type = bound_network<MaxNumOutputs>;
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

  /*! \brief Returns true since the network is combinational.
   *
   * TODO: Add support for sequential elements in the future.
   * \return Always returns true.
   */
  bool is_combinational() const
  {
    return true;
  }

  /*! \brief Test if a node is a multiple-output node.
   *
   * \param n The node index to check.
   * \return True if the node has multiple outputs, false otherwise.
   */
  bool is_multioutput( node_index_t const& n ) const
  {
    return _storage->is_multioutput( n );
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
  signal_t create_node( std::vector<signal_t> const& children, uint32_t id )
  {
    return create_node( children, std::vector<uint32_t>{ id } );
  }

  signal_t create_node( std::vector<signal_t> const& children, std::vector<uint32_t> const& ids )
  {
    signal_t const f = _storage->create_node( children, ids );
    set_value( f.index, 0 );

    for ( auto const& fn : _events->on_add )
    {
      ( *fn )( f.index );
    }
    return f;
  }

  signal clone_node( bound_network const& other, node_index_t const& source, std::vector<signal_t> const& children )
  {
    assert( !children.empty() );
    std::vector<uint32_t> const ids = other.get_binding_ids( source );
    return create_node( children, ids );
  }
#pragma endregion

#pragma region Restructuring
  void replace_in_node( node_index_t const& n, node_index_t const& old_node, signal_t new_signal )
  {
    if ( !_storage->in_fanin( n, old_node ) )
      return;

    /* if here old_node is in the fanin of n. Store current children to apply events */
    std::vector<signal_t> const old_children = _storage->get_children( old_node );

    /* replace in n's fanin the new node to the old one */
    _storage->replace_in_node( n, old_node, new_signal );

    for ( auto const& fn : _events->on_modified )
    {
      ( *fn )( n, old_children );
    }
  }

  void replace_in_node_no_restrash( node_index_t const& n, node_index_t const& old_node, signal_t new_signal )
  {
    replace_in_node( n, old_node, new_signal );
  }

  void replace_in_outputs( node_index_t const& old_node, signal_t const& new_signal )
  {
    if ( is_dead( old_node ) || !is_po( old_node ) )
      return;

    _storage->replace_in_outputs( old_node, new_signal );
  }

  void take_out_node( node_index_t const& n )
  {
    /* we cannot delete CIs, constants, or already dead nodes */
    if ( n < 2 || is_ci( n ) )
      return;

    /* delete the node */
    auto& nobj = _storage->nodes[n];
    nobj.kill();

    for ( auto const& fn : _events->on_delete )
    {
      ( *fn )( n );
    }

    /* if the node has been deleted, then deref fanout_size of
       fanins and try to take them out if their fanout_size become 0 */
    for ( auto i = 0; i < nobj.children.size(); ++i )
    {
      auto& child = nobj.children[i];
      if ( fanout_size( nobj.children[i] ) == 0 )
      {
        continue;
      }

      decr_fanout_size_pin( child );
      if ( decr_fanout_size( child.index ) == 0 )
      {
        take_out_node( child.index );
      }
    }
  }

  void revive_node( node_index_t const& n )
  {
    assert( !is_dead( n ) );
    return;
  }

  void substitute_node( node_index_t const& old_node, signal_t const& new_signal )
  {
    /* find all parents from old_node */
    signal f = new_signal;
    auto const& outputs = _storage->nodes[old_node].outputs;
    for ( auto i = 0u; i < outputs.size(); ++i )
    {
      f.output = i;
      for ( uint64_t idx : outputs[i].fanout )
      {
        replace_in_node( idx, old_node, f );
      }
    }

    /* check outputs */
    replace_in_outputs( old_node, new_signal );

    /* recursively reset old node */
    if ( old_node != new_signal.index )
    {
      take_out_node( old_node );
    }
  }

  void substitute_node_no_restrash( node_index_t const& old_node, signal_t const& new_signal )
  {
    substitute_node( old_node, new_signal );
  }

  inline bool is_dead( node_index_t const& n ) const
  {
    return _storage->is_dead( n );
  }
#pragma endregion

#pragma region Structural properties
  auto size() const
  {
    return static_cast<uint32_t>( _storage->nodes.size() );
  }

  auto num_cis() const
  {
    return static_cast<uint32_t>( _storage->inputs.size() );
  }

  auto num_cos() const
  {
    return static_cast<uint32_t>( _storage->outputs.size() );
  }

  auto num_pis() const
  {
    return static_cast<uint32_t>( _storage->inputs.size() );
  }

  auto num_pos() const
  {
    return static_cast<uint32_t>( _storage->outputs.size() );
  }

  auto num_gates() const
  {
    return static_cast<uint32_t>( _storage->nodes.size() - _storage->inputs.size() - 2 );
  }

  uint32_t num_outputs( node_index_t const& n ) const
  {
    return static_cast<uint32_t>( _storage->nodes[n].fanout_count );
  }

  uint32_t fanin_size( node_index_t const& n ) const
  {
    return static_cast<uint32_t>( _storage->nodes[n].children.size() );
  }

  uint32_t fanout_size( node_index_t const& n ) const
  {
    return _storage->nodes[n].fanout_count;
  }

  uint32_t incr_fanout_size( node_index_t const& n ) const
  {
    return _storage->nodes[n].fanout_count++;
  }

  uint32_t decr_fanout_size( node_index_t const& n ) const
  {
    return --_storage->nodes[n].fanout_count;
  }

  uint32_t incr_fanout_size_pin( node_index_t const& n, uint32_t pin_index ) const
  {
    return ++_storage->nodes[n].outputs[pin_index].fanout_count;
  }

  uint32_t decr_fanout_size_pin( node_index_t const& n, uint32_t pin_index ) const
  {
    return --_storage->nodes[n].outputs[pin_index].fanout_count;
  }

  uint32_t fanout_size_pin( node_index_t const& n, uint32_t pin_index ) const
  {
    return _storage->nodes[n].outputs[pin_index].fanout_count;
  }

  bool is_function( node_index_t const& n ) const
  {
    auto const& outputs = _storage->nodes[n].outputs;
    return ( outputs.size() > 0 ) && ( outputs[0].status == bound::pin_type_t::INTERNAL );
  }
#pragma endregion

#pragma region Functional properties
  kitty::dynamic_truth_table signal_function( const signal_t& f ) const
  {
    auto const& outputs = _storage->nodes[f.index].outputs;
    auto const& id = outputs[f.output].id;
    return _storage->library[id].function;
  }

  kitty::dynamic_truth_table node_function_pin( const node_index_t& n, uint32_t pin_index ) const
  {
    signal f = make_signal( n, pin_index );
    return signal_function( f );
  }
#pragma endregion

#pragma region Nodes and signals
  node get_node( signal_t const& f ) const
  {
    return f.index;
  }

  signal make_signal( node_index_t const& n, uint32_t output_pin ) const
  {
    return { n, output_pin };
  }

  signal make_signal( node_index_t const& n ) const
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

  signal next_output_pin( signal_t const& f ) const
  {
    return { f.index, f.output + 1 };
  }

  uint32_t node_to_index( node_index_t const& n ) const
  {
    return static_cast<uint32_t>( n );
  }

  node index_to_node( uint32_t index ) const
  {
    return index;
  }

  node ci_at( uint32_t index ) const
  {
    assert( index < _storage->inputs.size() );
    return *( _storage->inputs.begin() + index );
  }

  signal co_at( uint32_t index ) const
  {
    assert( index < _storage->outputs.size() );
    return *( _storage->outputs.begin() + index );
  }

  node pi_at( uint32_t index ) const
  {
    assert( index < _storage->inputs.size() );
    return *( _storage->inputs.begin() + index );
  }

  signal po_at( uint32_t index ) const
  {
    assert( index < _storage->outputs.size() );
    return *( _storage->outputs.begin() + index );
  }
#pragma endregion

#pragma region Node and signal iterators
  template<typename Fn>
  void foreach_node( Fn&& fn ) const
  {
    auto r = range<uint64_t>( _storage->nodes.size() );
    detail::foreach_element_if(
        r.begin(), r.end(),
        [this]( auto n ) { return !is_dead( n ); },
        fn );
  }

  template<typename Fn>
  void foreach_ci( Fn&& fn ) const
  {
    detail::foreach_element( _storage->inputs.begin(), _storage->inputs.end(), fn );
  }

  template<typename Fn>
  void foreach_co( Fn&& fn ) const
  {
    using IteratorType = decltype( _storage->outputs.begin() );
    detail::foreach_element<IteratorType>(
        _storage->outputs.begin(), _storage->outputs.end(), []( auto f ) { return signal( f ); }, fn );
  }

  template<typename Fn>
  void foreach_pi( Fn&& fn ) const
  {
    detail::foreach_element( _storage->inputs.begin(), _storage->inputs.end(), fn );
  }

  template<typename Fn>
  void foreach_po( Fn&& fn ) const
  {
    using IteratorType = decltype( _storage->outputs.begin() );
    detail::foreach_element<IteratorType>(
        _storage->outputs.begin(), _storage->outputs.end(), []( auto f ) { return signal( f ); }, fn );
  }

  template<typename Fn>
  void foreach_gate( Fn&& fn ) const
  {
    auto r = range<uint64_t>( 2u, _storage->nodes.size() ); /* start from 2 to avoid constants */
    detail::foreach_element_if(
        r.begin(), r.end(),
        [this]( auto n ) { return !is_ci( n ) && !is_dead( n ); },
        fn );
  }

  template<typename Fn>
  void foreach_fanin( node_index_t const& n, Fn&& fn ) const
  {
    if ( n == 0 || is_ci( n ) )
      return;

    using IteratorType = decltype( _storage->outputs.begin() );
    detail::foreach_element<IteratorType>(
        _storage->nodes[n].children.begin(), _storage->nodes[n].children.end(), []( auto f ) { return signal( f ); }, fn );
  }
#pragma endregion

#pragma region Simulate values
  template<typename TT>
  std::shared_ptr<list_simulator<list_t, TT>> get_simulator()
  {
    using simulator_t = list_simulator<list_t, TT>;
    static std::shared_ptr<simulator_t> sim = std::make_shared<simulator_t>( simulator_t() );
    return sim;
  }

  template<typename TT>
  std::vector<TT> compute( node_index_t const& n, std::vector<TT const*> sim_ptrs ) const
  {
    std::vector<TT> res;
    compute( res, n, sim_ptrs );
    return res;
  }

  template<typename TT>
  void compute( std::vector<TT>& res, node_index_t const& n, std::vector<TT const*> sim_ptrs ) const
  {
    auto simulator_ptr = get_simulator<TT>();
    auto const& nd = _storage->nodes[n];
    res.resize( nd.outputs.size() );
    const auto nfanin = nd.children.size();
    assert( nfanin > 0 );
    assert( sim_ptrs.size() == nfanin );

    for ( auto i = 0u; i < nd.output.size(); ++i )
    {
      bound::output_pin_t const& pin = nd.outputs[i];
      auto id = pin.id;
      auto const& list = _storage->library[id].get_list();
      ( *simulator_ptr )( list, sim_ptrs );
      simulator_ptr->get_simulation_inline( res[i], list, sim_ptrs, list.po_at( i ) );
    }
  }
#pragma endregion

#pragma region Custom node values
  void clear_values() const
  {
    std::for_each( _storage->nodes.begin(), _storage->nodes.end(), []( auto& n ) { n.user_data = 0; } );
  }

  uint32_t value( node_index_t const& n ) const
  {
    return _storage->nodes[n].user_data;
  }

  void set_value( node_index_t const& n, uint32_t v ) const
  {
    _storage->nodes[n].user_data = v;
  }

  uint32_t incr_value( node_index_t const& n ) const
  {
    return static_cast<uint32_t>( _storage->nodes[n].user_data++ );
  }

  uint32_t decr_value( node_index_t const& n ) const
  {
    return static_cast<uint32_t>( --_storage->nodes[n].user_data );
  }
#pragma endregion

#pragma region Visited flags
  void clear_visited() const
  {
    std::for_each( _storage->nodes.begin(), _storage->nodes.end(), []( auto& n ) { n.traversal_id = 0; } );
  }

  auto visited( node_index_t const& n ) const
  {
    return _storage->nodes[n].traversal_id;
  }

  void set_visited( node_index_t const& n, uint32_t v ) const
  {
    _storage->nodes[n].traversal_id = v;
  }

  uint32_t trav_id() const
  {
    return _storage->trav_id;
  }

  void incr_trav_id() const
  {
    ++_storage->trav_id;
  }
#pragma endregion

#pragma region General methods
  auto& events() const
  {
    return *_events;
  }
#pragma endregion

#pragma region Binding
  std::vector<uint32_t> get_binding_ids( node_index_t const& n )
  {
    return _storage->get_binding_ids( n );
  }
#pragma endregion

public:
  std::shared_ptr<bound::storage<NumBitsOutputs>> _storage;
  std::shared_ptr<network_events<base_type>> _events;
};

} // namespace mockturtle
