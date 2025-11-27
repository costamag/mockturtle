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
  \brief Bound logic network implementation with multi-output support

  The first two nodes in the network are reserved for the constant values:
  - 0 : constant 0
  - 1 : constant 1
  Differently from other mapped network types this network allows for more than
  one gate implementing the same functionality, consequently allowing us to
  include sizing into the description.

  \author Andrea Costamagna
*/

#pragma once

#include "bindings.hpp"
#include "../../utils/algorithm.hpp"
#include "../../utils/index_list.hpp"
#include "../../utils/truth_table_cache.hpp"
#include "../../networks/detail/foreach.hpp"
#include "../../networks/events.hpp"
#include "../../networks/storage.hpp"
#include "../../traits.hpp"

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>

#include <algorithm>
#include <iostream>
#include <memory>

namespace mockturtle
{

struct bound_storage_data
{
  bound_storage_data( std::vector<gate> const gates )
    : library( gates )
  {}

  bindings library;
};

/*! \brief Bound node
 *
 * `data[0].h1`  : Application-specific value
 * `data[1].h1`  : Visited flags
 * `data[1].h2`  : Total fan-out size (we use the MSB to indicate whether a node is dead)
 * `data[2].h1`  : Binding identifier of the node
 * `data[2].h2`  : Multiple-output group (we use the MSB to indicate that a node is not a multiple-output gate)
 *
 */
struct bound_storage_node
{
  static int constexpr Size = 3;
  static int constexpr PointerFieldSize = 0;
  using pointer_type = node_pointer<PointerFieldSize>;

  std::vector<uint64_t> fanout;
  std::vector<pointer_type> children;
  std::array<cauint64_t, Size> data;

  bool operator==( bound_storage_node const& other ) const
  {
    return data[2].h1 == other.data[2].h1 && children == other.children;
  }
};

/*! \brief Block storage container

  ...
*/
//template<typename Node, typename T = empty_storage_data>
struct bound_storage
{
  using node_type = bound_storage_node;
  using data_type = bound_storage_data;

  bound_storage( std::vector<gate> const& gates )
    :  data( gates )
  {
    nodes.reserve( 30000u );

    /* reserve the first two nodes for the constant */
    nodes.emplace_back();
    nodes.emplace_back();
  }

  uint32_t trav_id = 0u;

  std::vector<node_type> nodes;
  std::vector<uint64_t> inputs;
  std::vector<typename node_type::pointer_type> outputs;

  data_type data;
};

class bound_network
{
public:
#pragma region Types and constructors
  static constexpr auto min_fanin_size = 1;
  static constexpr auto max_fanin_size = 32;
  static constexpr uint64_t null_node = -1;
  static constexpr uint32_t msb32 = ( 1 << 31 );

  using base_type = bound_network;
  using storage = std::shared_ptr<bound_storage>;
  using node = uint64_t;
  using index_list_t = large_xag_index_list;

  struct signal
  {
    signal() = default;


    signal( bound_storage::node_type::pointer_type const& p )
        : complement( 0 ), index( p.index )
    {
    }

    signal( uint64_t index, uint64_t complement )
        : complement( complement ), index( index )
    {
    }

    signal( uint32_t index )
        : complement( 0 ), index( index )
    {
    }

    explicit signal( uint64_t data )
        : data( data )
    {
    }

    union
    {
      struct
      {
        uint64_t complement : 1;
        uint64_t index : 63;
      };
      uint64_t data;
    };

    signal operator!() const
    {
      return signal( data ^ 1 );
    }

    signal operator+() const
    {
      return { index, 0 };
    }

    signal operator-() const
    {
      return { index, 1 };
    }

    signal operator^( bool complement ) const
    {
      return signal( data ^ ( complement ? 1 : 0 ) );
    }

    bool operator==( signal const& other ) const
    {
      return data == other.data;
    }

    bool operator==( bound_storage::node_type::pointer_type const& p ) const
    {
      return index == p.index;
    }

    bool operator!=( signal const& other ) const
    {
      return data != other.data;
    }

    bool operator<( signal const& other ) const
    {
      return data < other.data;
    }

    operator uint64_t() const
    {
      return data;
    }
  };

  bound_network( std::vector<gate> const& gates )
      : _storage( std::make_shared<bound_storage>( gates ) ),
        _events( std::make_shared<decltype( _events )::element_type>() )
  {
    _init();
  }

  bound_network( std::shared_ptr<bound_storage> storage )
      : _storage( storage ),
        _events( std::make_shared<decltype( _events )::element_type>() )
  {
    _init();
  }

  bound_network( bound_network const& other )
      : bound_network( other._storage )
  {
  }

  bound_network( bound_network const& other, std::vector<gate> const & gates )
    : _storage( std::make_shared<bound_storage>( gates ) ),
      _events( std::make_shared<decltype( _events )::element_type>() )
  {
    _init();
    unordered_node_map<node, bound_network> old_to_new( other );
    other.foreach_pi( [&]( auto ni ){
      old_to_new[ni] = get_node(create_pi());
    } );
    other.foreach_po( [&]( signal f ){
      recursive_build( old_to_new, other, f );
      create_po( make_signal( old_to_new[other.get_node( f )] ) );
    } );
  }

  void recursive_build( unordered_node_map<node, bound_network> & old_to_new, bound_network const& other, signal f )
  {
    node n = other.get_node( f );
    if ( old_to_new.has( n ) )
    {
      return;
    }
    std::vector<signal> children;
    other.foreach_fanin( n, [&]( signal fi, uint32_t ii ){
      recursive_build( old_to_new, other, fi );
      node ni = other.get_node( fi );
      children.push_back( make_signal( old_to_new[ni] ) );
    } );
    auto bind = other.get_binding( n );
    old_to_new[n] = get_node( create_node( children, bind.name )[0] );
  }

  bound_network clone() const
  {
    return { std::make_shared<bound_storage>( *_storage ) };
  }

protected:
  inline void _init()
  {
    /* TODO: identifiers for constants */
    auto lib = _storage->data.library.get_gates();
    int id = 0;
    for( auto g : lib )
    {
      if( kitty::is_const0( g.function ) )
      {
        _storage->nodes[0].data[2].h1 = id;
      }
      if( kitty::is_const0( ~g.function ) )
      {
        _storage->nodes[1].data[2].h1 = id;
      }
      id++;
    }
  }
#pragma endregion

#pragma region Storage
public:
  bindings get_bindings()
  {
    return _storage->data.library;
  }
#pragma endregion

#pragma region Primary I / O and constants
public:
  signal get_constant( bool value = false ) const
  {
    uint32_t const c = value ? 1 : 0;
    return signal( { c, 0 } );
  }

  signal create_pi()
  {
    const auto index = _storage->nodes.size();
    _storage->inputs.emplace_back( index );
    _storage->nodes.emplace_back();
    _storage->nodes[index].data[2].h1 = msb32;
    return { index, 0 };
  }

  node get_null() const
  {
    return null_node;
  }

  uint32_t create_po( signal const& f )
  {
    /* increase ref-count to children */
    _storage->nodes[f.index].data[1].h2++;
    auto const po_index = static_cast<uint32_t>( _storage->outputs.size() );
    _storage->outputs.emplace_back( f.index );
    return po_index;
  }

  uint32_t po_index( signal const& f )
  {
    /* increase ref-count to children */
    for ( uint32_t i = 0; i < _storage->outputs.size(); ++i )
    {
      if ( _storage->outputs[i] == f.index )
        return i;
    }
    return -1;
  }

  void add_fanout( node n, node fanout_node )
  {
    uint64_t s = size();
    assert( ( fanout_node < s ) );
    _storage->nodes[n].fanout.push_back( fanout_node );
  }

  bool is_combinational() const
  {
    return true;
  }

  bool is_multioutput( node const& n ) const
  {
    return ( _storage->nodes[n].data[2].h2 & msb32 ) != msb32;
  }

  bool is_constant( node const& n ) const
  {
    return n <= 1 ;
  }

  bool is_null( node const& n ) const
  {
    return n == null_node;
  }

  bool is_and( node const& n ) const
  {
    (void)n;
    return false;
  }
  bool is_nand( node const& n ) const
  {
    (void)n;
    return false;
  }

  bool is_or( node const& n ) const
  {
    (void)n;
    return false;
  }

  bool is_nor( node const& n ) const
  {
    (void)n;
    return false;
  }

  bool is_lt( node const& n ) const
  {
    (void)n;
    return false;
  }

  bool is_ge( node const& n ) const
  {
    (void)n;
    return false;
  }

  bool is_gt( node const& n ) const
  {
    (void)n;
    return false;
  }

  bool is_le( node const& n ) const
  {
    (void)n;
    return false;
  }

  bool is_xor( node const& n ) const
  {
    (void)n;
    return false;
  }

  bool is_xnor( node const& n ) const
  {
    (void)n;
    return false;
  }

  bool is_xor3( node const& n ) const
  {
    (void)n;
    return false;
  }

  bool is_maj( node const& n ) const
  {
    (void)n;
    return false;
  }

  bool is_ite( node const& n ) const
  {
    (void)n;
    return false;
  }

  bool is_ci( node const& n ) const
  {
    return n > 1 && ( _storage->nodes[n].data[2].h1 == msb32 );
  }

  bool is_pi( node const& n ) const
  {
    return n > 1 && ( _storage->nodes[n].data[2].h1 == msb32 );
  }

  bool is_po( node const& n ) const
  {
    const uint64_t num_stored = _storage->nodes[n].fanout.size();
    const uint64_t num_counted = fanout_size( n );
    bool found = false;
    {
      foreach_po( [&]( auto f ){
        node no = get_node( f );
        found |= ( no == n );
      } );
    }
    return found;
  }

  bool constant_value( node const& n ) const
  {
    return n <= 1 ? n : 2;
  }

#pragma endregion

#pragma region Create arbitrary functions
  std::vector<signal> create_node( std::vector<signal> const& children, std::string const& name )
  {

    auto tt = _storage->data.library.get_function( name );
    if( kitty::is_const0( tt ) )
    {
      return { get_constant( false ) };
    }
    else if( kitty::is_const0( ~tt ) )
    {
      return { get_constant( true ) };
    }


    for( auto k : children )
    {
      uint64_t const index = get_node( k );
      bool const cec = index < static_cast<uint64_t>( _storage->nodes.size() );
      assert( cec && "wrong" );
    }

    if( name == "_const1_" )
      assert(0);

    auto const ids = _storage->data.library.get_ids( name );
    std::vector<signal> signals;
    int i = 0;
    for( auto const& id : ids )
    {
      storage::element_type::node_type node;
      for( auto const& f : children )
      {
        node.children.emplace_back( f.index );
      }
      auto bind = get_binding_from_id( id );
      assert( ( children.size() == bind.num_vars ) );

      node.data[2].h1 = id;
      node.data[2].h2 = ( ids.size() > 1 ) ? _storage->data.library.get_group_id( name ) : msb32;

      const auto index = _storage->nodes.size();
      _storage->nodes.push_back( node );

      /* increase ref-count to children */
      for ( auto c : children )
      {
        _storage->nodes[c.index].data[1].h2++;
        add_fanout( c.index, index );
      }
      set_value( index, 0 );

      for ( auto const& fn : _events->on_add )
      {
        ( *fn )( index );
      }
      signals.emplace_back( index, 0 );
      assert( (index < static_cast<uint32_t>( _storage->nodes.size() )) );
      assert( (index!=get_node(children[0])) );
    }
    return signals;
  }

  std::vector<signal> clone_node( bound_network const& other, node const& source, std::vector<signal> const& children )
  {
    assert( !children.empty() );

    std::string const& name = other.get_name( source );
    return create_node( children, name );
  }
#pragma endregion

#pragma region Restructuring
  void replace_in_node( node const& n, node const& old_node, signal new_signal )
  {
    bool in_fanin = false;
    auto& nobj = _storage->nodes[n];
    for ( auto& child : nobj.children )
    {
      if ( child.index == old_node )
      {
        in_fanin = true;
        break;
      }
    }

    if ( !in_fanin )
      return;

    // remember before
    std::vector<signal> old_children( nobj.children.size() );
    std::transform( nobj.children.begin(), nobj.children.end(), old_children.begin(), []( auto c ) { return signal{ c }; } );

    /* replace in node */
    for ( auto& child : nobj.children )
    {
      if ( child.index == old_node )
      {
        child = { new_signal.index };
        // increment fan-out of new node
        _storage->nodes[new_signal.index].data[1].h2++;
      }
    }
  }

  void replace_in_node_no_restrash( node const& n, node const& old_node, signal new_signal )
  {
    replace_in_node( n, old_node, new_signal );
  }

  void replace_in_outputs( node const& old_node, signal const& new_signal )
  {
//    if ( is_dead( old_node ) )
//      return;

    for ( auto& output : _storage->outputs )
    {
      if ( output.index == old_node )
      {

        output = { new_signal.index };
        assert( (output.index < static_cast<uint32_t>( _storage->nodes.size() )) );

        if ( old_node != new_signal.index )
        {
          /* increment fan-in of new node */
          _storage->nodes[new_signal.index].data[1].h2++;
        }
      }
    }
  }

  void take_out_node( node const& n )
  {
    /* we cannot delete CIs, constants, or already dead nodes */
    if ( n < 2 || is_ci( n ) )
      return;

    /* delete the node */
    auto& nobj = _storage->nodes[n];
    nobj.data[1].h2 = UINT32_C( 0x80000000 ); /* fanout size 0, but dead */
    nobj.data[1].h1 = UINT32_C( 0x00000000 ); /* fanout size 0, but dead */
    _storage->nodes[n].fanout.clear();

    for ( auto const& fn : _events->on_delete )
    {
      ( *fn )( n );
    }

    /* if the node has been deleted, then deref fanout_size of
       fanins and try to take them out if their fanout_size become 0 */
    for ( auto i = 0; i < nobj.children.size(); ++i )
    {
      auto& child = nobj.children[i];

      if ( fanout_size( nobj.children[i].index ) == 0 )
      {
        continue;
      }

      uint32_t ie = -1;
      node nc = child.index;
      foreach_fanout( nc, [&]( auto fo, auto io ){
        auto no = get_node( fo );
        if( no == n )
        {
          ie = io;
          return;
        }
      } );
      assert( ie >= 0 );
      auto & fanout = _storage->nodes[nc].fanout;//.erase()      assert( ie < fanout.size() );
      if( ie < fanout.size() )
      {
        fanout.erase( fanout.begin() + ie );
        if ( decr_fanout_size( nc ) < 1 )
        {
          take_out_node( nc );
        }
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
    assert((old_node != new_signal.index));
    /* find all parents from old_node */
    std::vector<node> fanout = _storage->nodes[old_node].fanout;
    for ( auto idx : fanout )
    {
      if ( is_ci( idx ) || is_dead( idx ) )
        continue; /* ignore CIs and dead nodes */

      replace_in_node( idx, old_node, new_signal );
    }
    node nnew = get_node( new_signal );
    for( auto no : fanout )
    {
      add_fanout( nnew, no );
    }
    std::vector<signal> new_children;
    for( auto child : _storage->nodes[nnew].children )
    {
      new_children.push_back( make_signal(child.index) );
    }

    std::vector<signal> old_children;
    for( auto child : _storage->nodes[old_node].children )
    {
      old_children.push_back( make_signal(child.index) );
    }

    /* check outputs */
    replace_in_outputs( old_node, new_signal );

    for ( auto const& fn : _events->on_modified )
    {
      ( *fn )( old_node, old_children );
    }

    /* clear the fanout of this node */
    _storage->nodes[old_node].fanout.clear();
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
    /* A dead node is simply a dangling node */
    return ( _storage->nodes[n].data[1].h2 >> 31 ) & 1;
  }
#pragma endregion


#pragma region Structural properties
  uint32_t size() const
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

  uint32_t get_degree( node const& n ) const
  {
    if ( is_pi( n ) )
      return fanout_size( n );
    else
      return fanout_size( n ) + fanin_size( n );
  }

  double get_range( node const& n ) const
  {
    if ( is_constant( n ) || is_pi( n ) )
    {
      double max_range = 0;
      if ( is_pi( n ) )
      {
        return sqrt( _storage->data.library.get_max_area() / M_PI );
      }
      return max_range;
    }
    else
      return sqrt( get_area( n ) / M_PI );
  }

  double get_range_from_id( uint32_t const& id ) const
  {
    auto bind = get_binding_from_id( id );
    return sqrt( bind.area / M_PI );
  }

  double area()
  {
    /* set the visited flags to 0 */
    clear_visited();

    /* set the traversal id to 1 */
    _storage->trav_id = 1;
    double area { 0 };
    foreach_po( [&]( auto fo, auto io ){
      area += area_rec( fo );
    } );
    return area;
  }

  double area_rec( signal f )
  {
    node const n = get_node( f );
    uint32_t v = _storage->nodes[n].data[1].h1;
    if( is_pi( n ) || is_dead( n ) || is_constant( n ) || ( v == trav_id() ) )
    {
      return 0;
    }
    set_visited( n, trav_id() );
    /* iterate over the fanis of the node */
    double area{ 0 };
    foreach_fanin( n, [&]( auto fi, auto ii ){
      area += area_rec( fi );
    } );
    return area + get_area( n );
  }

  double depth()
  {
    /* set the visited flags to 0 */
    clear_visited();
    incomplete_node_map<double, bound_network> map( *this );

    /* set the traversal id to 1 */
    _storage->trav_id = 1;
    double levels { 0 };
    foreach_po( [&]( auto fo, auto io ){
      levels = std::max( levels, depth_rec( fo, map ) );
    } );
    return levels;
  }

  double depth_rec( signal f, incomplete_node_map<double, bound_network> & map )
  {
    node const n = get_node( f );
    uint32_t v = _storage->nodes[n].data[1].h1;
    if( is_pi( n ) || is_dead( n ) || is_constant( n ) )
    {
      return 0;
    }
    if ( ( v == trav_id() ) )
      return map[n];
    set_visited( n, trav_id() );
    /* iterate over the fanis of the node */
    double levels{ 0 };

    auto g = get_binding( n );
    foreach_fanin( n, [&]( auto fi, auto ii ){
      levels = std::max( levels, depth_rec( fi, map ) + 1 );
    } );
    map[n] = levels;
    return levels;
  }

  uint32_t num_outputs( node const& n ) const
  {
    return static_cast<uint32_t>( _storage->nodes[n].data.size() - 2 );
  }

  uint32_t fanin_size( node const& n ) const
  {
    return static_cast<uint32_t>( _storage->nodes[n].children.size() );
  }

  uint32_t fanout_size( node const& n ) const
  {
    return _storage->nodes[n].data[1].h2 & UINT32_C( 0x7FFFFFFF );
  }

  uint32_t incr_fanout_size( node const& n ) const
  {
    return _storage->nodes[n].data[1].h2++ & UINT32_C( 0x7FFFFFFF );
  }

  uint32_t decr_fanout_size( node const& n ) const
  {
    return --_storage->nodes[n].data[1].h2 & UINT32_C( 0x7FFFFFFF );
  }

  bool is_function( node const& n ) const
  {
    return n > 1 && !is_ci( n );
  }
  
  bool has_binding( node const& n ) const
  {
    return get_binding_id( n ) < _storage->data.library.size();
  }

  bool has_binding( signal const& f ) const
  {
    return has_binding( f.index );
  }
#pragma endregion

#pragma region Functional properties
  kitty::dynamic_truth_table node_function( const node& n ) const
  {
    unsigned int const id = _storage->nodes[n].data[2].h1;
    return _storage->data.library.get_function( id );
  }
#pragma endregion

#pragma region Nodes and signals
  node get_node( signal const& f ) const
  {
    return f.index;
  }

  std::string get_name( node const& n ) const
  {
    uint32_t id = get_binding_id( n );
    return _storage->data.library.get_name( id );
  }

  double get_area( node const& n ) const
  {
    if ( is_pi( n ) )
      return 1;

    uint32_t id = get_binding_id( n );
    return _storage->data.library.get_area( id );
  }

  double get_power( node const& n ) const
  {
    auto tt = node_function( n );
    double const n0 = kitty::count_zeros( tt );
    double const nb = tt.num_bits();
    double const n1 = nb - n0;
    return 2*n1*n0/( nb*nb );
  }

  unsigned int get_binding_id( node const& n ) const
  {
    return _storage->nodes[n].data[2].h1;
  }

  unsigned int get_binding_index( node const& n ) const
  {
    return get_binding_id( n );
  }

  gate_with_list get_binding( node const& n ) const
  {
    unsigned int id = get_binding_id( n );
    return _storage->data.library.get_binding( id );
  }

  auto get_library() const
  {
    return _storage->data.library.get_gates();
  }

  gate_with_list get_binding_from_id( uint32_t const& id ) const
  {
    return _storage->data.library.get_binding( id );
  }

  signal make_signal( node const& n ) const
  {
    return { n, 0 };
  }

  bool is_complemented( signal const& f ) const
  {
    return false;
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
    auto r = range<uint64_t>( 2, _storage->nodes.size() );
    detail::foreach_element_if(
        r.begin(), r.end(),
        [this]( auto n ) { return !is_dead( n ); },
        fn );
  }

  template<typename Fn>
  void foreach_node( node start, Fn&& fn ) const
  {
    auto r = range<uint64_t>( start, _storage->nodes.size() );
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
    detail::foreach_element_transform<IteratorType, signal>(
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
    detail::foreach_element_transform<IteratorType, signal>(
        _storage->outputs.begin(), _storage->outputs.end(), []( auto f ) { return signal( f ); }, fn );
  }

  template<typename Fn>
  void foreach_po_offset( uint32_t first, Fn&& fn ) const
  {
    using IteratorType = decltype( _storage->outputs.begin() );
    detail::foreach_element_transform<IteratorType, signal>(
        _storage->outputs.begin() + first, _storage->outputs.end(), []( auto f ) { return signal( f ); }, fn );
  }


  template<typename Fn>
  void foreach_gate( Fn&& fn ) const
  {
    auto r = range<uint64_t>( 2u, _storage->nodes.size() ); /* start from 2 to avoid constants */
    detail::foreach_element_if(
        r.begin(), r.end(),
        [this]( auto n ) { return !is_ci( n ); },
        fn );
  }

  template<typename Fn>
  void foreach_fanin( node const& n, Fn&& fn ) const
  {
    if ( n <= 0 || is_ci( n ) )
      return;

    using IteratorType = decltype( _storage->outputs.begin() );
    detail::foreach_element_transform<IteratorType, signal>(
        _storage->nodes[n].children.begin(), _storage->nodes[n].children.end(), []( auto f ) { return signal( f ); }, fn );
  }

  std::vector<signal> get_children( node const& n )
  {
    std::vector<signal> children;
    for ( auto c : _storage->nodes[n].children )
      children.push_back( signal( c ) );
    return children;
  }

  std::vector<node> get_nodes_fanin( node const& n )
  {
    std::vector<node> children;
    for ( auto c : _storage->nodes[n].children )
      children.push_back( get_node(signal( c )) );
    return children;
  }

  template<typename Fn>
  void foreach_fanout( node const& n, Fn&& fn ) const
  {
    using IteratorType = std::vector<node>::iterator;
    detail::foreach_element_transform<IteratorType, signal>(
      _storage->nodes[n].fanout.begin(), _storage->nodes[n].fanout.end(), []( auto no ) { return no; }, fn );
  }
#pragma endregion

#pragma region Simulate values
  template<typename Iterator>
  iterates_over_t<Iterator, bool>
  compute( node const& n, Iterator begin, Iterator end ) const
  {
    uint32_t index{ 0 };
    auto it = _storage->nodes[n].children.begin();
    while ( begin != end )
    {
      index <<= 1;
      index ^= *begin++;
      ++it;
    }
    auto const tt = node_function( n );
    return kitty::get_bit( tt, index );
  }

  template<class Tt_t>
  void simulate_list( Tt_t& result, index_list_t const& list, std::vector<Tt_t> & tts ) const
  {
    list.foreach_gate( [&]( uint32_t lit0, uint32_t lit1 ){
      uint32_t const index0 = ( lit0 >> 1u ) - 1;
      uint32_t const index1 = ( lit1 >> 1u ) - 1;
      bool const is_compl0 = ( lit0 & 0x1 ) == 0x1;
      bool const is_compl1 = ( lit1 & 0x1 ) == 0x1;
      Tt_t const sim0 = ( lit0 & 0x1 ) == 0x1 ? ~tts[index0] : tts[index0];
      Tt_t const sim1 = ( lit1 & 0x1 ) == 0x1 ? ~tts[index1] : tts[index1];
      tts.emplace_back( sim0 & sim1 );
    } );
    /* the last entry of the list is the output */
    uint32_t const lit = list.po_at(0);
    uint32_t const index = ( lit >> 1u ) - 1;
    bool const is_compl = ( lit & 0x1 ) == 0x1;
    result = is_compl ? ~tts[index] : tts[index];
  }

  template<typename Iterator>
  iterates_over_truth_table_t<Iterator>
  compute( node const& n, Iterator begin, Iterator end ) const
  {
    const auto nfanin = _storage->nodes[n].children.size();

    std::vector<typename std::iterator_traits<Iterator>::value_type> tts( begin, end );

    assert( nfanin != 0 );
    assert( tts.size() == nfanin );

    /* resulting truth table has the same size as any of the children */
    auto result = tts.front().construct();
    unsigned int id = get_binding_id( n );
    index_list_t list = _storage->data.library.get_list( id );
    simulate_list( result, list, tts );

    return result;
  }
#pragma endregion

#pragma region Custom node values
  void clear_values() const
  {
    std::for_each( _storage->nodes.begin(), _storage->nodes.end(), []( auto& n ) { n.data[0].h1 = 0; } );
  }

  uint32_t value( node const& n ) const
  {
    return _storage->nodes[n].data[0].h1;
  }

  void set_value( node const& n, uint32_t v ) const
  {
    _storage->nodes[n].data[0].h1 = v;
  }

  uint32_t incr_value( node const& n ) const
  {
    return static_cast<uint32_t>( _storage->nodes[n].data[0].h1++ );
  }

  uint32_t decr_value( node const& n ) const
  {
    return static_cast<uint32_t>( --_storage->nodes[n].data[0].h1 );
  }
#pragma endregion

#pragma region Visited flags
  void clear_visited() const
  {
    std::for_each( _storage->nodes.begin(), _storage->nodes.end(), []( auto& n ) { n.data[1].h1 = 0; } );
  }

  auto visited( node const& n ) const
  {
    return _storage->nodes[n].data[1].h1;
  }

  void set_visited( node const& n, uint32_t v ) const
  {
    assert( !is_null( n ) );
    _storage->nodes[n].data[1].h1 = v;
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

#pragma region Visualize
  void print() const
  {
    foreach_gate( [&]( auto n ){
      std::string name = get_name( n );
      std::cout << n << " " << name << " ";
      foreach_fanin( n, [&]( auto fi, auto ii ){
        std::cout << fi.index << " ";
      } );
      std::cout << " FO: ";
      foreach_fanout( n, [&]( auto fo ){
        std::cout << fo << " ";
      } );
      std::cout << std::endl;
      return true;
    } );

    foreach_po( [&]( auto fo ){
      std::cout << fo.index << " ";
    } );
    std::cout << " <- POs" << std::endl;
  }
#pragma endregion Visualize

public:
  std::shared_ptr<bound_storage> _storage;
  std::shared_ptr<network_events<base_type>> _events;
};

} // namespace mockturtle
