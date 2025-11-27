/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
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
  \file mapped.hpp
  \brief maped logic network implementation

  \author Andrea Costamagna
*/

#pragma once

#include "../algorithms/node_resynthesis/dsd.hpp"
#include "../algorithms/node_resynthesis/shannon.hpp"
#include "../algorithms/node_resynthesis/xag_npn.hpp"
#include "../algorithms/decomposition.hpp"
#include "../io/genlib_reader.hpp"
#include "../traits.hpp"
#include "../utils/algorithm.hpp"

#include "aig.hpp"
#include "detail/foreach.hpp"
#include "events.hpp"
#include "storage.hpp"

#include <kitty/properties.hpp>
#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>

#include <algorithm>
#include <bits/stdc++.h>
#include <memory>

namespace mockturtle
{

struct mapped_storage_data
{
  aig_network aig;
  std::vector<gate> library;
  /*! \brief mapped node
 *
 * `0`: Constant 0
 * `1`: Constant 1
 * `2`: Buffer
 * `3`: Inverter
 */
  std::array<unsigned int, 4u> special;
};

/*! \brief mapped node
 *
 * `data[0].h1`: Fan-out size
 * `data[0].h2`: Application-specific value
 * `data[1].h1`: Binding id equivalent to the library index
 * `data[1].h2`: Visited flags
 */
struct mapped_storage_node : mixed_fanin_node<2>
{
  bool operator==( mapped_storage_node const& other ) const
  {
    return data[1].h1 == other.data[1].h1 && children == other.children;
  }
};

/*! \brief k-LUT storage container

  ...
*/
using mapped_storage = storage<mapped_storage_node, mapped_storage_data>;

class mapped_network
{
public:
#pragma region Types and constructors
  static constexpr auto min_fanin_size = 1;
  static constexpr auto max_fanin_size = 32;

  using base_type = mapped_network;
  using storage = std::shared_ptr<mapped_storage>;
  using node = uint64_t;
  using signal = uint64_t;
  unsigned int id_pi = std::numeric_limits<unsigned int>::max(); 

  mapped_network()
      : _storage( std::make_shared<mapped_storage>() ),
        _events( std::make_shared<decltype( _events )::element_type>() )
  {
    _init();
  }

  mapped_network( std::vector<gate> const& library )
      : _storage( std::make_shared<mapped_storage>() ),
        _events( std::make_shared<decltype( _events )::element_type>() )
  {
    _init( library );
  }

  mapped_network( std::shared_ptr<mapped_storage> storage )
      : _storage( storage ),
        _events( std::make_shared<decltype( _events )::element_type>() )
  {
    _init();
  }

  template<class Ntk>
  void construct_rec( Ntk const& ntk, typename Ntk::signal f, unordered_node_map<node, Ntk> & old_to_new )
  {
    typename Ntk::node n = ntk.get_node( f );
    if( old_to_new.has( n ) )
    {
      return;
    }
    std::vector<signal> children;

    ntk.foreach_fanin( n, [&]( auto fi ){
      typename Ntk::node ni = ntk.get_node( fi );

      construct_rec( ntk, fi, old_to_new );
      children.push_back( make_signal(old_to_new[ni]) );
      return true;
    } );
    unsigned int id = ntk.get_binding( n ).id;
    auto fnew = create_node( children, id );
    old_to_new[n] = get_node( fnew );
  }

  template<class Ntk>
  mapped_network( std::vector<gate> const& library, Ntk const& other )
      : _storage( std::make_shared<mapped_storage>() ),
        _events( std::make_shared<decltype( _events )::element_type>() )
  {
    _init( library );

    unordered_node_map<node, Ntk> old_to_new( other ); 
    other.foreach_pi( [&]( auto f ){
      old_to_new[other.get_node( f )] = create_pi();
      return true;
    } );

    other.foreach_po( [&]( const auto& f, auto index ) {
      construct_rec( other, f, old_to_new );
      return true;
    } );

    other.foreach_po( [&]( const auto& f, auto index ) {
      create_po( old_to_new[f] );
      return true;
    } );
  }

  mapped_network clone() const
  {
    return { std::make_shared<mapped_storage>( *_storage ) };
  }

protected:
  inline void _init()
  {
    /* already initialized */
    if ( _storage->nodes.size() > 1 )
      return;
    _storage->nodes.clear();
  }

  inline void _init( std::vector<gate> const& library )
  {
    /* already initialized */
    if ( _storage->nodes.size() > 1 )
      return;

    //kitty::print_binary( new_tt ); printf(" %d\n", support.size() );
    xag_npn_resynthesis<aig_network, aig_network, xag_npn_db_kind::aig_complete> aig_resyn{};

    _storage->nodes.clear();
    /* reserve the second node for constant 0 and 1 */
    _storage->nodes.emplace_back();
    _storage->nodes.emplace_back();

    /* generate max_fanins pis in the AIG */
    std::vector<aig_network::signal> pis( max_fanin_size );
    std::generate( pis.begin(), pis.end(), [&]() { return _storage->data.aig.create_pi(); } );
    /* store the library */
    _storage->data.library = library;

    /* for each gate in the library generate the corresponding AIG */



    for( auto g : _storage->data.library )
    {
      /* save a reference to the constants */
      if( kitty::is_const0( g.function ) )
      {
        _storage->data.special[0] = g.id;
        /* truth tables for constant 0 */
        _storage->nodes[0].data[1].h1 = g.id;
      }
      else if( kitty::is_const0( ~g.function ) )
      {
        _storage->data.special[1] = g.id;
        /* truth tables for constant 0 */
        _storage->nodes[1].data[1].h1 = g.id;
      }
      else if( kitty::to_binary( g.function ) == "10" )
      {
        _storage->data.special[2] = g.id;
      }
      else if( kitty::to_binary( g.function ) == "01" )
      {
        _storage->data.special[3] = g.id;
      }

      /* save a reference to the constants */
      std::vector<aig_network::signal> vars( g.function.num_vars() );
      for( int i{0}; i<g.num_vars; ++i )
      {
        vars[i] = pis[i];
      }

      _storage->data.aig.create_po( synthesize_aig( vars, g.function, aig_resyn ) );

    }
  }

  template<class Resyn>
  aig_network::signal synthesize_aig( std::vector<aig_network::signal> aig_children, kitty::dynamic_truth_table const& tt, Resyn & fn )
  {
    if( kitty::is_const0( tt ) )
      return aig_network::signal{ 0, 0 };
    if( kitty::is_const0( ~tt ) )
      return aig_network::signal{ 0, 1 };
    if( aig_children.size() == 1u )
      return kitty::is_normal(tt) ? aig_children[0] : !aig_children[0];
    if( aig_children.size() <= 4u )
    {
      return match_twin( aig_children, tt, fn );
    }

    uint32_t idx = aig_children.size()-1;
    aig_network::signal x = aig_children[idx];
    aig_children.erase( aig_children.begin() + idx );
    aig_network::signal f1 = synthesize_aig( aig_children, kitty::cofactor1( tt, idx ), fn );
    aig_network::signal f0 = synthesize_aig( aig_children, kitty::cofactor0( tt, idx ), fn );
    
    if( f1.index == 0 )
    {
      auto res = f1.complement ? !_storage->data.aig.create_and( !x, !f0 ) : _storage->data.aig.create_and( !x, f0 );
      return res;
    }
    if( f0.index == 0 )
    {
      auto res = f0.complement ? !_storage->data.aig.create_and( x, !f1 ) : _storage->data.aig.create_and( x, f1 );
      return res;
    }

    return _storage->data.aig.create_ite( x, f1, f0 );
  }

  template<class Resyn>
  aig_network::signal match_twin( std::vector<aig_network::signal> aig_children, kitty::dynamic_truth_table tt, Resyn & fn )
  {
    const auto support = kitty::min_base_inplace( tt );
    kitty::dynamic_truth_table new_tt = kitty::shrink_to( tt, static_cast<unsigned int>( support.size() ) );

    for( int iVar{aig_children.size()-1}; iVar>=0; --iVar )
    {
      if( std::find( support.begin(), support.end(), iVar ) == support.end() )
        aig_children.erase( aig_children.begin()+iVar );//.push_back( inputs[support[iVar]] );
    }

    aig_network::signal out_sig;
    fn( _storage->data.aig, new_tt, aig_children.begin(), aig_children.end(), [&]( auto const& f ) {
            out_sig = f;
            return false;
          } );
    return out_sig;
  }

#pragma endregion

#pragma region Primary I / O and constants
public:
  signal get_constant( bool value = false ) const
  {
    return value ? 1 : 0;
  }

  signal create_pi()
  {
    const auto index = _storage->nodes.size();
    _storage->nodes.emplace_back();
    _storage->inputs.emplace_back( index );
    /* the identifier of a PI is a "virtual gate", the first available after storing the library */
    _storage->nodes[index].data[1].h1 = _storage->data.library.size();
    return index;
  }

  uint32_t create_po( signal const& f )
  {
    /* increase ref-count to children */
    _storage->nodes[f].data[0].h1++;
    auto const po_index = static_cast<uint32_t>( _storage->outputs.size() );
    _storage->outputs.emplace_back( f );
    return po_index;
  }

  bool is_combinational() const
  {
    return true;
  }

  bool is_constant( node const& n ) const
  {
    return n <= 1;
  }

  bool is_ci( node const& n ) const
  {
    return _storage->nodes[node_to_index(n)].data[1].h1 == _storage->data.library.size();
  }

  bool is_pi( node const& n ) const
  {
    return _storage->nodes[node_to_index(n)].data[1].h1 == _storage->data.library.size();
  }

  bool is_buffer( node const& n ) const
  {
    return kitty::to_binary( get_binding( n ).function ) == "10";
  }

  bool constant_value( node const& n ) const
  {
    assert( n <= 1 );
    return n == 1;
  }

  uint32_t po_index( signal const& s ) const
  {
    uint32_t i = -1;
    foreach_po( [&]( const auto& x, auto index ) {
      if ( x == s )
      {
        i = index;
        return false;
      }
      return true;
    } );
    return i;
  }
#pragma endregion

#pragma region Create unary functions
  signal create_buf( signal const& a )
  {
    return _create_node( { a }, _storage->data.special[2] );
  }

  signal create_not( signal const& a )
  {
    return _create_node( { a }, _storage->data.special[3] );
  }
#pragma endregion

#pragma region Create arbitrary functions
  signal _create_node( std::vector<signal> const& children, unsigned int id )
  {
    if ( children.size() == 0u )
    {
      assert( _storage->data.library[id].num_vars == 0u );
      return get_constant( !kitty::is_const0( _storage->data.library[id].function ) );
    }
    storage::element_type::node_type node;
    std::copy( children.begin(), children.end(), std::back_inserter( node.children ) );
    node.data[1].h1 = id;

    const auto it = _storage->hash.find( node );
    if ( it != _storage->hash.end() )
    {
      return it->second;
    }

    const auto index = _storage->nodes.size();
    _storage->nodes.push_back( node );
    _storage->hash[node] = index;

    /* increase ref-count to children */
    for ( auto c : children )
    {
      _storage->nodes[c].data[0].h1++;
    }

    set_value( index, 0 );

    for ( auto const& fn : _events->on_add )
    {
      ( *fn )( index );
    }

    return index;
  }

  signal create_node( std::vector<signal> const& children, unsigned int id )
  {
    return _create_node( children, id );
  }

  signal clone_node( mapped_network const& other, node const& source, std::vector<signal> const& children )
  {
    assert( !children.empty() );
    return _create_node( children, other._storage->nodes[source].data[1].h1 );
  }
#pragma endregion

#pragma region Restructuring
  void substitute_node( node const& old_node, signal const& new_signal )
  {
    /* find all parents from old_node */
    for ( auto i = 0u; i < _storage->nodes.size(); ++i )
    {
      auto& n = _storage->nodes[i];
      for ( auto& child : n.children )
      {
        if ( child == old_node )
        {
          std::vector<signal> old_children( n.children.size() );
          std::transform( n.children.begin(), n.children.end(), old_children.begin(), []( auto c ) { return c.index; } );
          child = new_signal;

          // increment fan-out of new node
          _storage->nodes[new_signal].data[0].h1++;

          for ( auto const& fn : _events->on_modified )
          {
            ( *fn )( i, old_children );
          }
        }
      }
    }

    /* check outputs */
    for ( auto& output : _storage->outputs )
    {
      if ( output == old_node )
      {
        output.index = new_signal;

        // increment fan-out of new node
        _storage->nodes[new_signal].data[0].h1++;
      }
    }

    // reset fan-out of old node
    _storage->nodes[old_node].data[0].h1 = 0;
    kill_node( old_node );
    // recursively delete mffc
    delete_redundant( old_node );
  }

  void delete_redundant( node const& n )
  {
    foreach_fanin( n, [&]( auto f ){
      auto ni = get_node( f );
      _storage->nodes[ni].data[0].h1--;
      if( _storage->nodes[ni].data[0].h1 == 0 )
      {
        kill_node( ni );
        delete_redundant( ni );
      }
    } );
  }

  inline bool is_dead( node const& n ) const
  {
    return ( _storage->nodes[n].data[0].h1 >> 31u & 0x1 ) == 0x1;
  }

  void kill_node( node const& n )
  {
    _storage->nodes[n].data[0].h1 = ( 0x1 << 31u );
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

  uint32_t fanin_size( node const& n ) const
  {
    return static_cast<uint32_t>( _storage->nodes[n].children.size() );
  }

  uint32_t fanout_size( node const& n ) const
  {
    return _storage->nodes[n].data[0].h1;
  }

  uint32_t incr_fanout_size( node const& n ) const
  {
    return _storage->nodes[n].data[0].h1++;
  }

  uint32_t decr_fanout_size( node const& n ) const
  {
    return --_storage->nodes[n].data[0].h1;
  }

  bool is_function( node const& n ) const
  {
    return n > 1 && !is_ci( n );
  }

  void compute_area_rec( signal f, double & area )
  {
    node n = get_node( f );
    if( visited( n ) == 1 )
    {
      return;
    }
    area += get_binding( n ).area;

    set_visited( n, 1 );
    foreach_fanin( n, [&]( auto fi ){
      compute_area_rec( fi, area );
    } );

  }

  double compute_area()
  {
    clear_visited();
    foreach_pi( [&]( auto n ){
      set_visited( n, 1u );
    } );

    double area{0};
    foreach_po( [&]( auto f ){
      compute_area_rec( f, area );
    } );
    clear_visited();
    return area;
  }
#pragma endregion

#pragma region Functional properties
  kitty::dynamic_truth_table node_function( const node& n ) const
  {
    return _storage->data.library[_storage->nodes[n].data[1].h1].function;
  }

  gate get_binding( const node& n ) const
  {
    return _storage->data.library[_storage->nodes[n].data[1].h1];
  }
#pragma endregion

#pragma region Nodes and signals
  node get_node( signal const& f ) const
  {
    return f;
  }

  signal make_signal( node const& n ) const
  {
    return n;
  }

  bool is_complemented( signal const& f ) const
  {
    (void)f;
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
    return ( _storage->outputs.begin() + index )->index;
  }

  node pi_at( uint32_t index ) const
  {
    assert( index < _storage->inputs.size() );
    return *( _storage->inputs.begin() + index );
  }

  signal po_at( uint32_t index ) const
  {
    assert( index < _storage->outputs.size() );
    return ( _storage->outputs.begin() + index )->index;
  }
#pragma endregion

#pragma region Node and signal iterators
  template<typename Fn>
  void foreach_node( Fn&& fn ) const
  {
    auto r = range<uint64_t>( _storage->nodes.size() );
    detail::foreach_element( r.begin(), r.end(), fn );
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
    detail::foreach_element_transform<IteratorType, uint32_t>(
        _storage->outputs.begin(), _storage->outputs.end(), []( auto o ) { return o.index; }, fn );
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
    detail::foreach_element_transform<IteratorType, uint32_t>(
        _storage->outputs.begin(), _storage->outputs.end(), []( auto o ) { return o.index; }, fn );
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
    if ( n == 0 || is_ci( n ) )
      return;

    using IteratorType = decltype( _storage->outputs.begin() );
    detail::foreach_element_transform<IteratorType, uint32_t>(
        _storage->nodes[n].children.begin(), _storage->nodes[n].children.end(), []( auto f ) { return f.index; }, fn );
  }
#pragma endregion

#pragma region Value simulation
  template<typename Iterator>
  iterates_over_t<Iterator, bool>
  compute( node const& n, Iterator begin, Iterator end ) const
  {
    uint32_t index{ 0 };
    while ( begin != end )
    {
      index <<= 1;
      index ^= *begin++ ? 1 : 0;
    }
    return kitty::get_bit( _storage->data.library[_storage->nodes[n].data[1].h1].function, index );
  }

  template<typename TT>
  TT compute_rec( aig_network::signal f, std::vector<aig_network::signal> const& children, std::vector<TT> const & tts, unordered_node_map<TT, aig_network>& node_to_tt ) const
  {
    aig_network const& aig = _storage->data.aig;
    aig_network::node n = aig.get_node( f );
    assert( !aig.is_constant(n) );

    if( node_to_tt.has( n ) )
    {
      return aig.is_complemented( f ) ? ~node_to_tt[n] : node_to_tt[n];
    }

    TT res = tts[0]^tts[0];

    if( aig.is_pi( n ) )
    {
      auto index = aig.pi_index( n );
      return aig.is_complemented( f ) ? ~tts[index] : tts[index];
    }
    else
    {
      aig_network::signal const & a = aig._storage->nodes[n].children[0];
      aig_network::signal const & b = aig._storage->nodes[n].children[1];
      TT sim_a = compute_rec( a, children, tts, node_to_tt );
      TT sim_b = compute_rec( b, children, tts, node_to_tt );
      auto res =  sim_a & sim_b;
      node_to_tt[n] = res;
      return aig.is_complemented(f) ? ~res : res;
    }
  }

  template<typename Iterator>
  iterates_over_truth_table_t<Iterator>
  compute( node const& n, Iterator begin, Iterator end ) const
  {
    (void)end;
    assert( n != 0 && !is_ci( n ) );

    const auto num_fanins = _storage->nodes[n].children.size();

    std::vector<typename std::iterator_traits<Iterator>::value_type> tts( begin, end );

    /* resulting truth table has the same size as any of the children */
    assert( tts.size() == num_fanins );

    std::vector<aig_network::signal> children( num_fanins );
    auto& aig_inputs = _storage->data.aig._storage->inputs;
    std::transform( children.begin(), children.end(), aig_inputs.begin(), [&]( auto c ) { return _storage->data.aig.get_node( c ); } );

    unordered_node_map<typename std::iterator_traits<Iterator>::value_type, aig_network> node_to_tt( _storage->data.aig );
    return compute_rec( _storage->data.aig.po_at( _storage->nodes[n].data[1].h1 ), children, tts, node_to_tt );
  }

  /*! \brief Re-compute the last block. */
  template<typename Iterator>
  void compute( node const& n, kitty::partial_truth_table& result, Iterator begin, Iterator end ) const
  {
    static_assert( iterates_over_v<Iterator, kitty::partial_truth_table>, "begin and end have to iterate over partial_truth_tables" );

    (void)end;
    assert( n != 0 && !is_ci( n ) );

    const auto num_fanins = _storage->nodes[n].children.size();
    std::vector<aig_network::signal> children( num_fanins );
    auto& aig_inputs = _storage->data.aig._storage->inputs;
    std::transform( aig_inputs.begin(), aig_inputs.begin() + num_fanins, children.begin(), [&]( auto c ) { return _storage->data.aig.make_signal( c ); } );

    std::vector<typename std::iterator_traits<Iterator>::value_type> tts( begin, end );

    assert( tts[0].num_bits() > 0 && "truth tables must not be empty" );
    assert( tts[0].num_bits() >= result.num_bits() );
    assert( result.num_blocks() == tts[0].num_blocks() || ( result.num_blocks() == tts[0].num_blocks() - 1 && result.num_bits() % 64 == 0 ) );

    result.resize( tts[0].num_bits() );

    unordered_node_map<uint64_t, aig_network> node_to_tt( _storage->data.aig );
    std::vector<uint64_t> bits( tts.size() );
    std::transform( tts.begin(), tts.end(), bits.begin(), [&]( auto tt ) { return tt._bits.back(); } );

    result._bits.back() = compute_rec( _storage->data.aig.po_at( _storage->nodes[n].data[1].h1 ), children, bits, node_to_tt );
    result.mask_bits();
  }
#pragma endregion

#pragma region delay


  double compute_arrival_rec( node n, unordered_node_map<double, mapped_network>& delays ) const
  {
    if( delays.has(n) )
    {
      return delays[n];
    }

    auto const& g = get_binding( n );
    double arrival = 0;

    foreach_fanin( n, [&]( auto const& f, auto i ) {
      double arr_fanin = compute_arrival_rec( get_node(f), delays ); 
      double pin_delay = std::max( g.pins[i].rise_block_delay, g.pins[i].fall_block_delay );
      arrival = std::max( arrival, (double)( arr_fanin + pin_delay ) );
    } );
    delays[n] = arrival;

    return delays[n];
  }



  double compute_worst_delay() const
  {
    unordered_node_map<double, mapped_network> delays( *this );

    foreach_pi( [&]( auto const& n, auto i ) {
      delays[n]=0.0;
    } );
    delays[0]=0.0;

    double max_delay=0;
    foreach_po( [&]( auto const& fo, auto i ) {
      node no = get_node(fo);
      double out_del = compute_arrival_rec( no, delays );
      if( out_del > max_delay )
      {
        max_delay = out_del;
      }
    } );
    
    return max_delay;
  }

#pragma endregion

#pragma region Custom node values
  void clear_values() const
  {
    std::for_each( _storage->nodes.begin(), _storage->nodes.end(), []( auto& n ) { n.data[0].h2 = 0; } );
  }

  uint32_t value( node const& n ) const
  {
    return _storage->nodes[n].data[0].h2;
  }

  void set_value( node const& n, uint32_t v ) const
  {
    _storage->nodes[n].data[0].h2 = v;
  }

  uint32_t incr_value( node const& n ) const
  {
    return static_cast<uint32_t>( _storage->nodes[n].data[0].h2++ );
  }

  uint32_t decr_value( node const& n ) const
  {
    return static_cast<uint32_t>( --_storage->nodes[n].data[0].h2 );
  }
#pragma endregion

#pragma region Visited flags
  void clear_visited() const
  {
    std::for_each( _storage->nodes.begin(), _storage->nodes.end(), []( auto& n ) { n.data[1].h2 = 0; } );
  }

  uint32_t visited( node const& n ) const
  {
    return _storage->nodes[n].data[1].h2;
  }

  void set_visited( node const& n, uint32_t v ) const
  {
    _storage->nodes[n].data[1].h2 = v;
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

public:
  std::shared_ptr<mapped_storage> _storage;
  std::shared_ptr<network_events<base_type>> _events;
};

} // namespace mockturtle