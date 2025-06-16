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

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>

#include <algorithm>
#include <memory>

namespace mockturtle
{

/*! \brief Pointer to a node with output-pin specifier.
 *
 * This data structure contains the information to point to an output pin
 * of a node. The information is stored in an `uint64_t`, partitioned as
 * follows:
 * - `NumBitsOutputs` bits are used to indicate the output pin
 * - `64 - NumBitsOutputs` bits are used to specify the node index.
 *
 * \tparam NumBitsOutputs Number of bits for the output id 
 */
template<int NumBitsOutputs>
struct bound_storage_signal
{
public:
  bound_storage_signal() = default;
  bound_storage_signal( uint64_t index, uint64_t output ) : index( index ), output( output ) {}
  bound_storage_signal( uint64_t data ) : data( data ) {}

  union
  {
    struct
    {
      uint64_t output : NumBitsOutputs;
      uint64_t index : 64 - NumBitsOutputs;
    };
    uint64_t data;
  };

  bool operator==( bound_storage_signal<NumBitsOutputs> const& other ) const
  {
    return data == other.data;
  }

  bool operator!=( bound_storage_signal<NumBitsOutputs> const& other ) const
  {
    return data != other.data;
  }

  operator uint64_t() const
  {
    return data;
  }
};

/*! \brief Contains the information of an output pin of a gate.
 *
 * To support multiple-output gates, a node might have more than
 * one output pin. This structure contains the information which
 * specifies a specific output pin. This implementation assumes
 * that the technology library is implemented as a vector pin names
 * each implementing an unique functionality. The binding id is the
 * index of the pin's functionality.
 */
struct bound_output_pin
{
  using node_index_t = uint64_t;

  bound_output_pin( uint32_t id, pin_status status, std::vector<node_index_t> const& fanout )
   : id( id ), status( status ), fanout( fanout )
  {}

  bound_output_pin( uint32_t id, pin_status status )
   : bound_output_pin( id, status, {} )
  {}

  bound_output_pin()
   : bound_output_pin( std::numeric_limits<uint32_t>::max(), pin_status::NONE, {} )
  {}

  /*! \brief Binding id from the genlib-like library */
  uint32_t id;
  /*! \brief Specifies whether the node is dead, a constant, a PI, a PO, etc. */
  pin_status status;
  /*! \brief Contains the indeces of the nodes connected to the pin */
  std::vector<node_index_t> fanout;
};

/*! \brief Node information 
 *
 * \tparam N Maximum number of outputs of the cells in the library.
 */
template<uint32_t N>
struct bound_storage_node
{
  using pointer_type = bound_storage_signal<N>;

  bound_storage_node()
  {
    outputs = decltype( outputs )( 1 );
  }

  bound_storage_node( pin_status status )
  {
    outputs = decltype( outputs )( 1 );
    outputs[0].status = status;
  }
  
  /*! \brief Kill the node by setting the status of its output pins to dead.
  */
  void kill()
  {
    for ( auto const& pin : outputs )
    {
      pin.status = pin_status::DEAD;
    }
  }

  bool operator==( bound_storage_node<N> const& other ) const
  {
    return children == other.children;
  }

  /*! \brief Pointer to the signals in the immediate fanin. */
  std::vector<pointer_type> children;
  /*! \brief Application-specific value. */
  uint32_t value { 0 };
  /*! \brief Visited flag. */
  uint32_t visit { 0 };
  /*! \brief Total fan-out size ( we use MSB to indicate whether a node is dead ). */
  uint32_t num_fanouts { 0 };
  /*! \brief Vector of output pins. */
  std::vector<bound_output_pin> outputs;
};

template<uint32_t N>
struct bound_storage
{
  using storage_node = bound_storage_node<N>;
  using node = uint64_t;
  using signal = typename storage_node::pointer_type;
  using list_t = large_xag_index_list;

  bound_storage( std::vector<gate> const& gates )
  : library( gates )
  {
    nodes.reserve( 10000u );

    /* we generally reserve the first node for a constant */
    nodes.emplace_back();
    nodes.emplace_back();
  }

  signal get_constant( bool value )
  {
    return value ? signal{ 1, 0 } : signal{ 0, 0 };
  }

  signal create_pi()
  {
    const auto index = nodes.size();
    nodes.emplace_back( pin_status::PI );
    inputs.emplace_back( index );
    return signal{ index, 0 };
  }

  uint32_t create_po( signal const& f )
  {
    /* increase ref-count to children */
    nodes[f.index].num_fanouts++;
    nodes[f.index].outputs[f.output].status = pin_status::PO;
    auto const po_index = static_cast<uint32_t>( outputs.size() );
    outputs.emplace_back( f.index, f.output );
    return po_index;
  }

  bool is_multioutput( node const& n ) const
  {
    return nodes[n].outputs.size() > 1;
  }

  bool is_constant( node const& n )
  {
    auto const& pins = nodes[n].outputs;
    return pins[0].status == pin_status::CONSTANT;
  }

  bool is_ci( node const& n )
  {
    auto const& pins = nodes[n].outputs;
    return ( pins[0].status == pin_status::PI ) || 
    ( pins[0].status == pin_status::CI );
  }

  bool is_pi( node const& n ) const
  {
    auto const& pins = nodes[n].outputs;
    return ( pins[0].status == pin_status::PI );
  }

  bool is_po( node const& n, uint32_t output = 0 ) const
  {
    auto const& pins = nodes[n].outputs;
    return ( pins[output].status == pin_status::PO );
  }

  bool is_dead( node const& n )
  {
    bool all_dead { true };
    bool one_dead { false };
    for ( auto const& pin : nodes[n].outputs )
    {
      bool const dead = pin.status == pin_status::DEAD;
      all_dead &= dead;
      one_dead |= dead;
    }
    assert( !( all_dead ^ one_dead ) );
    /* A dead node is simply a dangling node */
    return all_dead;
  }

  bool constant_value( node const& n ) const
  {
    return n != 0;
  }

  signal create_node( std::vector<signal> const& children, std::vector<uint32_t> const& ids )
  {
    storage_node new_node;
    std::copy( children.begin(), children.end(), std::back_inserter( new_node.children ) );

    new_node = decltype( new_node.outputs )( ids.size() );

    for ( auto i = 0; i < ids.size(); ++i )
      new_node.outputs[i] = { ids[i], pin_status::INTERNAL };

    const auto index = nodes.size();
    nodes.push_back( new_node );

    /* increase ref-count to children */
    for ( auto c : children )
    {
      nodes[c.index].num_fanouts++;
      nodes[c.index].outputs[c.output].fanouts.push_back( index );
    }

    return { index, 0 };
  }

  std::vector<uint32_t> get_binding_ids( node const& n )
  {
    std::vector<uint32_t> ids;
    for ( auto const& pin : nodes[n].outputs )
      ids.push_back( pin.output );
    return ids;
  }

  bool in_fanin( node n, node old_node )
  {
    bool in_fanin = false;
    auto& nobj = nodes[n];
    for ( auto& child : nobj.children )
    {
      if ( child.index == old_node )
      {
        in_fanin = true;
        break;
      }
    }

    return in_fanin;
  }

  std::vector<signal> get_children( node const& n )
  {
    auto& nobj = nodes[n];
    std::vector<signal> children( nobj.children.size() );
    std::transform( nobj.children.begin(), nobj.children.end(), children.begin(), []( auto c ) { return signal{ c }; } );
    return children;
  }

  void replace_in_node( node const& n, node const& old_node, signal new_signal )
  {
    auto& nobj = nodes[n];
    for ( auto& child : nobj.children )
    {
      if ( child.index == old_node )
      {
        child = signal{ new_signal.data };
        nodes[new_signal.index].num_fanouts++;
        nodes[new_signal.index].outputs[new_signal.output].fanouts.push_back( n );
      }
    }
  }

  void replace_in_outputs( node const& old_node, signal const& new_signal )
  {
    for ( auto& output : outputs )
    {
      if ( output.index == old_node )
      {
        if ( old_node != new_signal.index )
        {
          /* increment fan-in of new node */
          nodes[new_signal.index].num_fanouts++;
          nodes[new_signal.index].outputs[new_signal.output].status = pin_status::PO;
        }
      }
    }
  }

  uint32_t trav_id = 0u;

  std::vector<storage_node> nodes;
  std::vector<node> inputs;
  std::vector<signal> outputs;
  augmented_library<gate> library;
};

/*! \brief Network of gates from a technology library.
 *
 * \tparam MaxNumOutputs Maximum number of outputs of the cells in the library.
 */
template<uint32_t MaxNumOutputs = 2u>
class bound_network
{
  public:
  #pragma region Types and constructors
  static constexpr auto min_fanin_size = 1;
  static constexpr auto max_fanin_size = 32;
  static constexpr auto min_gate_output_size = 1;
  static constexpr auto num_bits_output_code = static_cast<uint32_t>( ceil( log2( MaxNumOutputs ) ) );

  using base_type = bound_network<MaxNumOutputs>;
  using storage_node = bound_storage_node<num_bits_output_code>;
  using aig_list = typename bound_storage<num_bits_output_code>::list_t;
  using storage = std::shared_ptr<bound_storage<num_bits_output_code>>;
  using node = uint64_t;
  using signal = bound_storage_signal<num_bits_output_code>;

  bound_network( std::vector<gate> const& gates )
      : _storage( std::make_shared<bound_storage<num_bits_output_code>>( gates ) ),
        _events( std::make_shared<typename decltype( _events )::element_type>() )
  {
  }

  bound_network( std::shared_ptr<bound_storage<num_bits_output_code>> storage )
      : _storage( storage ),
        _events( std::make_shared<typename decltype( _events )::element_type>() )
  {
  }

  bound_network clone() const
  {
    return { std::make_shared<bound_storage<num_bits_output_code>>( *_storage ) };
  }

#pragma endregion

#pragma region Primary I / O and constants
public:
  /* */
  signal get_constant( bool value = false ) const
  {
    return _storage->get_constant( value );
  }

  signal create_pi()
  {
    return _storage->create_pi();
  }

  uint32_t create_po( signal const& f )
  {
    return _storage->create_po( f );
  }

  bool is_combinational() const
  {
    return true;
  }

  bool is_multioutput( node const& n ) const
  {
    return _storage->is_multioutput( n );
  }

  bool is_constant( node const& n ) const
  {
    return _storage->is_constant( n );
  }

  bool is_ci( node const& n ) const
  {
    return _storage->is_ci( n );
  }
  
  bool is_pi( node const& n ) const
  {
    return _storage->is_pi( n );
  }

  bool is_po( node const& n, uint32_t output = 0 ) const
  {
    return _storage->is_po( n, output );
  }

  bool constant_value( node const& n ) const
  {
    return _storage->constant_value( n );
  }
#pragma endregion

#pragma region Create arbitrary functions
  signal create_node( std::vector<signal> const& children, uint32_t id )
  {
    return create_node( children, { id } );
  }
  
  signal create_node( std::vector<signal> const& children, std::vector<uint32_t> const& ids )
  {
    signal const f = _storage->create_node( children, ids );
    set_value( f.index, 0 );
    
    for ( auto const& fn : _events->on_add )
    {
      ( *fn )( f.index );
    }
    return f;
  }

  signal clone_node( bound_network const& other, node const& source, std::vector<signal> const& children )
  {
    assert( !children.empty() );
    std::vector<uint32_t> const ids = other.get_binding_ids( source );
    return create_node( children, ids );
  }
#pragma endregion

#pragma region Restructuring
  void replace_in_node( node const& n, node const& old_node, signal new_signal )
  {
    if ( !_storage->in_fanin( n, old_node ) )
      return;
    
    /* if here old_node is in the fanin of n. Store current children to apply events */
    std::vector<signal> const old_children = _storage->get_children( old_node );

    /* replace in n's fanin the new node to the old one */
    _storage->replace_in_node( n, old_node, new_signal ); 

    for ( auto const& fn : _events->on_modified )
    {
      ( *fn )( n, old_children );
    }
  }

  void replace_in_node_no_restrash( node const& n, node const& old_node, signal new_signal )
  {
    replace_in_node( n, old_node, new_signal );
  }

  void replace_in_outputs( node const& old_node, signal const& new_signal )
  {
    if ( is_dead( old_node ) || !is_po( old_node ) )
      return;

    _storage->replace_in_outputs( old_node, new_signal );
  }

  void take_out_node( node const& n )
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

  void revive_node( node const& n )
  {
    assert( !is_dead( n ) );
    return;
  }

  void substitute_node( node const& old_node, signal const& new_signal )
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

  void substitute_node_no_restrash( node const& old_node, signal const& new_signal )
  {
    substitute_node( old_node, new_signal );
  }

  inline bool is_dead( node const& n ) const
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

  uint32_t num_outputs( node const& n ) const
  {
    return static_cast<uint32_t>( _storage->nodes[n].num_fanouts );
  }

  uint32_t fanin_size( node const& n ) const
  {
    return static_cast<uint32_t>( _storage->nodes[n].children.size() );
  }

  uint32_t fanout_size( node const& n ) const
  {
    return _storage->nodes[n].num_fanouts;
  }

  uint32_t incr_fanout_size( node const& n ) const
  {
    return _storage->nodes[n].num_fanouts++;
  }

  uint32_t decr_fanout_size( node const& n ) const
  {
    return --_storage->nodes[n].num_fanouts;
  }

  uint32_t incr_fanout_size_pin( node const& n, uint32_t pin_index ) const
  {
    return ++_storage->nodes[n].outputs[pin_index].num_fanouts;
  }

  uint32_t decr_fanout_size_pin( node const& n, uint32_t pin_index ) const
  {
    return --_storage->nodes[n].outputs[pin_index].num_fanouts;
  }

  uint32_t fanout_size_pin( node const& n, uint32_t pin_index ) const
  {
    return _storage->nodes[n].outputs[pin_index].num_fanouts;
  }

  bool is_function( node const& n ) const
  {
    auto const & outputs = _storage->nodes[n].outputs;
    return ( outputs.size() > 0) && ( outputs[0].status == pin_status::INTERNAL );
  }
#pragma endregion

#pragma region Functional properties
  kitty::dynamic_truth_table signal_function( const signal& f ) const
  {
    auto const& outputs = _storage->nodes[f.index].outputs;
    auto const& id = outputs[f.output].id;
    return _storage->library[id].function;
  }

  kitty::dynamic_truth_table node_function_pin( const node& n, uint32_t pin_index ) const
  {
    signal f = make_signal( n, pin_index );
    return signal_function( f );
  }
#pragma endregion

#pragma region Nodes and signals
  node get_node( signal const& f ) const
  {
    return f.index;
  }
  
  signal make_signal( node const& n, uint32_t output_pin ) const
  {
    return { n, output_pin };
  }

  signal make_signal( node const& n ) const
  {
    return make_signal( n, 0 );
  }

  bool is_complemented( signal const& f ) const
  {
    (void)f;
    return false;
  }

  uint32_t get_output_pin( signal const& f ) const
  {
    return static_cast<uint32_t>( f.output );
  }

  signal next_output_pin( signal const& f ) const
  {
    return { f.index, f.output + 1 };
  }

  uint32_t node_to_index( node const& n ) const
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
  void foreach_fanin( node const& n, Fn&& fn ) const
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
  std::shared_ptr<list_simulator<aig_list, TT>> get_simulator() {
    using simulator_t = list_simulator<aig_list, TT>;
    static std::shared_ptr<simulator_t> sim = std::make_shared<simulator_t>( simulator_t() );
    return sim;
  }

  template<typename TT>
  std::vector<TT> compute( node const& n, std::vector<TT const *> sim_ptrs ) const
  {
    std::vector<TT> res;
    compute( res, n, sim_ptrs );
    return res;
  }

  template<typename TT>
  void compute( std::vector<TT> & res, node const& n, std::vector<TT const *> sim_ptrs ) const
  {
    auto simulator_ptr = get_simulator<TT>();
    auto const& nd = _storage->nodes[n];
    res.resize( nd.outputs.size() );
    const auto nfanin = nd.children.size();
    assert( nfanin > 0 );
    assert( sim_ptrs.size() == nfanin );

    for ( auto i = 0u; i < nd.output.size(); ++i )
    {
      bound_output_pin const& pin = nd.outputs[i];
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
    std::for_each( _storage->nodes.begin(), _storage->nodes.end(), []( auto& n ) { n.value = 0; } );
  }

  uint32_t value( node const& n ) const
  {
    return _storage->nodes[n].value;
  }

  void set_value( node const& n, uint32_t v ) const
  {
    _storage->nodes[n].value = v;
  }

  uint32_t incr_value( node const& n ) const
  {
    return static_cast<uint32_t>( _storage->nodes[n].value++ );
  }

  uint32_t decr_value( node const& n ) const
  {
    return static_cast<uint32_t>( --_storage->nodes[n].value );
  }
#pragma endregion

#pragma region Visited flags
  void clear_visited() const
  {
    std::for_each( _storage->nodes.begin(), _storage->nodes.end(), []( auto& n ) { n.visit = 0; } );
  }

  auto visited( node const& n ) const
  {
    return _storage->nodes[n].visit;
  }

  void set_visited( node const& n, uint32_t v ) const
  {
    _storage->nodes[n].visit = v;
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
std::vector<uint32_t> get_binding_ids( node const& n )
{
  return _storage->get_binding_ids( n );
}
#pragma endregion

public:
  std::shared_ptr<bound_storage<num_bits_output_code>> _storage;
  std::shared_ptr<network_events<base_type>> _events;
};

} // namespace mockturtle
