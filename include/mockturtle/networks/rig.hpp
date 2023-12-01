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
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYrigHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file rig.hpp
  \brief Representation Independent Graph

  \author Andrea Costamagna
*/

#pragma once

#include "../utils/truth_table_cache.hpp"
#include "../utils/algorithm.hpp"
#include "detail/foreach.hpp"
#include "../utils/tech_library.hpp"
#include "../algorithms/node_resynthesis/xag_npn.hpp"
#include "aig.hpp"
#include <kitty/constructors.hpp>
#include <kitty/operations.hpp>
#include "storage.hpp"
#include <kitty/kitty.hpp>
#include <algorithm>
#include <memory>
#include <list>

namespace mockturtle
{

namespace rils
{

/*! \brief identifier of the node category */
enum i_func_t : uint32_t
{
  i_CONST = 0u,
  i_PI = 1u,
  i_BUF = 2u
};

/*! \brief literals of the precomputed truth-tables in the truth-table cache.
  *
  * The default convention for literals is assumed.  That is an index \f$i\f$
  * (starting) from \f$0\f$ has positive literal \f$2i\f$ and negative literal
  * \f$2i + 1\f$.
  * 
  */
enum e_func_t : uint32_t
{
  e_CONST = 0u,
  e_PI = 1u,
  e_BUF = 2u,
  e_AND = 4u,
  e_OR = 6u,
  e_LT = 8u,
  e_GT = 10u,
  e_XOR = 12u
};

#pragma region utils
  template<typename PTR_T>
  struct signal_t
  {
    signal_t() = default;

    signal_t( uint64_t index, uint64_t complement )
        : complement( complement ), index( index )
    {
    }

    signal_t( uint32_t index )
        : complement( 0 ), index( index )
    {
    }

    signal_t( uint64_t index, uint64_t complement, uint64_t output )
        : complement( complement ), index( index )
    {
    }

    explicit signal_t( uint64_t data )
        : data( data )
    {
    }

    signal_t( PTR_T const& p )
        : complement( p.weight & 1 ), index( p.index )
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

    signal_t operator!() const
    {
      return signal_t( data ^ 1 );
    }

    signal_t operator+() const
    {
      return { index, 0 };
    }

    signal_t operator-() const
    {
      return { index, 1 };
    }

    signal_t operator^( bool complement ) const
    {
      return signal_t( data ^ ( complement ? 1 : 0 ) );
    }

    bool operator==( signal_t const& other ) const
    {
      return data == other.data;
    }

    bool operator!=( signal_t const& other ) const
    {
      return data != other.data;
    }

    bool operator<( signal_t const& other ) const
    {
      return data < other.data;
    }

    operator PTR_T() const
    {
      return { index, complement };
    }

    operator uint64_t() const
    {
      return data;
    }

  #if __cplusplus > 201703L
    bool operator==( PTR_T const& other ) const
    {
      return data == other.data;
    }
  #endif
  };

  /*! \brief rig luts storage
  */
  struct e_data_t
  {
    truth_table_cache<kitty::dynamic_truth_table> cache;
  };

  struct e_gate_t
  {
    using pointer_type = node_pointer<1>;
    using twin_pointer_type = node_pointer<1>;

    std::vector<pointer_type> children;

    /*! \brief number of fanouts */
    uint32_t nfos{0};
    /*! \brief id of the functionality stored in the tt-cache */
    uint32_t func{0};
    /*! \brief application specific bits */
    uint32_t bits{0};
    /*! \brief Twin internal signal */
    signal_t<twin_pointer_type> twin {0,0};

    bool operator==( e_gate_t const& other ) const
    {
      return func == other.func && children == other.children;
    }
  };

  struct i_gate_t
  {
    using pointer_type = node_pointer<1>;
    using twin_pointer_type = node_pointer<1>;

    std::array<pointer_type, 2> children;
    /*! \brief number of fanouts */
    uint32_t nfos{0};
    /*! \brief node type const:0 pi:1 */
    uint32_t func{0};
    /*! \brief application specific bits */
    uint32_t bits{0};
    /*! Twin external signal: W might not exit! */
    signal_t<twin_pointer_type> twin {0,0};

    bool operator==( i_gate_t const& other ) const
    {
      return children == other.children;
    }

  };

  /*! \brief Hash function for AIGs (from ABC) */
  struct i_hash_t
  {
    uint64_t operator()( i_gate_t const& n ) const
    {
      uint64_t seed = -2011;
      seed += n.children[0].index * 7937;
      seed += n.children[1].index * 2971;
      seed += n.children[0].weight * 911;
      seed += n.children[1].weight * 353;
      return seed;
    }
  };
  using e_storage_t = smart_storage<e_gate_t, e_data_t>;

  using i_data_t = empty_storage_data;
  struct i_hash_t;
  using i_storage_t = smart_storage<i_gate_t, i_data_t, i_hash_t>;
#pragma endregion utils

class rig_network
{

  #pragma region types
  public:
    using e_node_t = uint64_t;
    using e_signal_t = signal_t<e_gate_t::pointer_type>;

    using i_node_t = uint64_t;
    using i_signal_t = signal_t<i_gate_t::pointer_type>;

    static constexpr auto min_fanin_size = 1;
    static constexpr auto max_fanin_size = 32;
    using base_t = rig_network;
    using node = e_node_t;
    using signal = e_signal_t;
  #pragma endregion types

  #pragma region constructor
  public:
    rig_network();

  protected:
    inline void _init();
  #pragma endregion constructors

  #pragma region linking
  private:
    e_signal_t get_e_signal( i_signal_t const& );
    i_signal_t get_i_signal( e_signal_t const& );
  #pragma endregion linking

  #pragma region Primary I / O and constants
  public:
    signal get_constant( bool );
    signal create_pi();
    uint32_t create_po( signal const& );
    bool is_combinational() const;
    bool is_constant( node const& ) const;
    bool is_ci( node const& ) const;
    bool is_pi( node const& ) const;
  #pragma endregion Primary I / O and constants

  #pragma region nodes and signals
  public:
    node get_node( signal const& ) const;
    signal make_signal( node const& ) const;
    bool is_complemented( signal const& ) const;
    uint32_t node_to_index( node const& ) const;
    node index_to_node( uint32_t ) const;
    node pi_at( uint32_t ) const;
    node ci_at( uint32_t ) const;
    signal po_at( uint32_t ) const;
    signal co_at( uint32_t ) const;
    uint32_t ci_index( node const& ) const;
    uint32_t pi_index( node const& ) const;
    //uint32_t co_index( signal const& ) const;
    //uint32_t po_index( signal const& ) const;
  #pragma endregion nodes and signals

  #pragma region node and signal iterators
    template<typename Fn> void foreach_node( Fn&& fn ) const;
    template<typename Fn> void foreach_ci( Fn&& fn ) const;
    template<typename Fn> void foreach_co( Fn&& fn ) const;
    template<typename Fn> void foreach_pi( Fn&& fn ) const;
    template<typename Fn> void foreach_po( Fn&& fn ) const;
    template<typename Fn> void foreach_gate( Fn&& fn ) const;
    template<typename Fn> void foreach_fanin( node const& n, Fn&& fn ) const;
  #pragma endregion node and signal iterators

  #pragma region unary functions
    signal create_buf( signal const& );
    signal create_not( signal const& );
    bool is_buf( node const& );
    bool is_not( node const& );
  #pragma endregion unary functions

  #pragma region binary functions
    signal i_create_and( i_signal_t, i_signal_t );
    signal create_and( signal const&, signal const& );
    signal create_nand( signal const&, signal const& );
    signal create_or( signal const&, signal const& );
    signal create_nor( signal const&, signal const& );
    signal create_lt( signal const&, signal const& );
    signal create_ge( signal const&, signal const& );
    signal create_gt( signal const&, signal const& );
    signal create_le( signal const&, signal const& );
    signal create_xor( signal const&, signal const& );
    signal create_xnor( signal const&, signal const& );
    bool is_and( node const& );
    bool is_nand( node const& );
    bool is_or( node const& );
    bool is_nor( node const& );
    bool is_lt( node const& );
    bool is_ge( node const& );
    bool is_gt( node const& );
    bool is_le( node const& );
    bool is_xor( node const& );
    bool is_xnor( node const& );
  #pragma endregion binary functions

  #pragma region arbitrary function
    signal _create_node( std::vector<signal> const&, uint32_t );
    std::tuple<rig_network::signal, bool> _create_known_node( std::vector<signal> const&, uint32_t );
    bool is_function( node const& ) const;

    i_signal_t synthesize_twin( std::array<i_signal_t, 4u> &, uint32_t );
    i_signal_t synthesize_twin_rec( std::array<i_signal_t, 4u> &, kitty::dynamic_truth_table const& );
    i_signal_t boolean_matching( std::array<i_signal_t, 4u> &, kitty::dynamic_truth_table const& );
  #pragma endregion arbitrary function

  #pragma region structural properties
  public:
    size_t size() const;
    size_t num_cis() const;
    size_t num_cos() const;
    size_t num_pis() const;
    size_t num_pos() const;
    size_t num_gates() const;
    size_t fanin_size( node const& n ) const;
    size_t fanout_size( node const& n ) const;
    size_t incr_fanout_size( node const& n ) const;
    size_t decr_fanout_size( node const& n ) const;
    //bool is_and( node const& n ) const;
    //bool is_or( node const& n ) const;
    //bool is_xor( node const& n ) const;
    //bool is_maj( node const& n ) const;
    //bool is_ite( node const& n ) const;
    //bool is_xor3( node const& n ) const;
  #pragma endregion structural properties

  #pragma region functional properties
  public:
    kitty::dynamic_truth_table node_function( const node& ) const;
  #pragma endregion functional properties

public:
  std::shared_ptr<e_storage_t> _e_storage;
  std::shared_ptr<i_storage_t> _i_storage;

};


#pragma region constructors
  /*! \brief Network constructor

  * Construct the network using the init function
  */
  rig_network::rig_network()
    : _e_storage( std::make_shared<e_storage_t>() ), _i_storage( std::make_shared<i_storage_t>() )
  {
    _init();
  }
  /*! \brief Network initializer

  * At initialization, the network must have allocated only one node for constant 0.
  * This method stores the truth table of the constant function and connects the constant 0
  * node of the externale and the internal representations.
  */
  inline void rig_network::_init()
  {
    /* already initialized */
    if ( _e_storage->nodes.size() > 1 )
      return;

    /* constant node : #0 in the cache */
    kitty::dynamic_truth_table tt_zero( 0 );
    _e_storage->data.cache.insert( tt_zero );
    _i_storage->nodes[0].func = i_func_t::i_CONST;

    static uint64_t _not = 0x1;
    kitty::dynamic_truth_table tt_not( 1 );
    kitty::create_from_words( tt_not, &_not, &_not + 1 );
    _e_storage->data.cache.insert( tt_not );

    static uint64_t _and = 0x8;
    kitty::dynamic_truth_table tt_and( 2 );
    kitty::create_from_words( tt_and, &_and, &_and + 1 );
    _e_storage->data.cache.insert( tt_and );

    static uint64_t _or = 0xe;
    kitty::dynamic_truth_table tt_or( 2 );
    kitty::create_from_words( tt_or, &_or, &_or + 1 );
    _e_storage->data.cache.insert( tt_or );

    static uint64_t _lt = 0x4;
    kitty::dynamic_truth_table tt_lt( 2 );
    kitty::create_from_words( tt_lt, &_lt, &_lt + 1 );
    _e_storage->data.cache.insert( tt_lt );

    static uint64_t _le = 0xd;
    kitty::dynamic_truth_table tt_le( 2 );
    kitty::create_from_words( tt_le, &_le, &_le + 1 );
    _e_storage->data.cache.insert( tt_le );

    static uint64_t _xor = 0x6;
    kitty::dynamic_truth_table tt_xor( 2 );
    kitty::create_from_words( tt_xor, &_xor, &_xor + 1 );
    _e_storage->data.cache.insert( tt_xor );

    static uint64_t _maj = 0xe8;
    kitty::dynamic_truth_table tt_maj( 3 );
    kitty::create_from_words( tt_maj, &_maj, &_maj + 1 );
    _e_storage->data.cache.insert( tt_maj );

    static uint64_t _ite = 0xd8;
    kitty::dynamic_truth_table tt_ite( 3 );
    kitty::create_from_words( tt_ite, &_ite, &_ite + 1 );
    _e_storage->data.cache.insert( tt_ite );

    static uint64_t _xor3 = 0x96;
    kitty::dynamic_truth_table tt_xor3( 3 );
    kitty::create_from_words( tt_xor3, &_xor3, &_xor3 + 1 );
    _e_storage->data.cache.insert( tt_xor3 );

    /* truth tables for constants */
    _e_storage->nodes[0].func = 0;

  }
  #pragma endregion constructors

#pragma region linking
  rig_network::e_signal_t rig_network::get_e_signal( rig_network::i_signal_t const& f )
  {
    return _e_storage->nodes[f].twin;
  }

  rig_network::i_signal_t rig_network::get_i_signal( rig_network::e_signal_t const& f )
  {
    return _i_storage->nodes[f].twin;
  }
#pragma endregion linking

#pragma region Primary I / O and constants
  rig_network::signal rig_network::get_constant( bool value = false )
  {
    return { 0, static_cast<uint64_t>( value ? 1 : 0 ) };
  }

  rig_network::signal rig_network::create_pi()
  {
    const auto e_index = _e_storage->get_index();
    auto& e_node = _e_storage->nodes.emplace_back();
    e_node.children.emplace_back( static_cast<uint64_t>(_e_storage->inputs.size()) );
    _e_storage->inputs.push_back( e_index );
    _e_storage->nodes[e_index].func = e_func_t::e_PI;

    const auto i_index = _i_storage->get_index();
    auto& i_node = _i_storage->nodes.emplace_back();
    i_node.children[0].data = i_node.children[1].data = _i_storage->inputs.size();
    _i_storage->inputs.emplace_back( i_index );
    _i_storage->nodes[i_index].func = i_func_t::i_PI;

    _e_storage->nodes[e_index].twin = { i_index, 0 };
    _i_storage->nodes[i_index].twin = { e_index, 0 };
    
    return { e_index, 0 };
  }

  uint32_t rig_network::create_po( rig_network::signal const& e_signal )
  {
    /* increase ref-count to children */
    _e_storage->nodes[e_signal.index].nfos++;
    auto const e_po_index = _e_storage->outputs.size();
    _e_storage->outputs.emplace_back( e_signal );
    
    i_signal_t i_signal = get_i_signal( e_signal );
    _i_storage->nodes[i_signal.index].nfos++;
    i_signal_t i_out_signal { i_signal.index, i_signal.complement ^ e_signal.complement };
    _i_storage->outputs.emplace_back( i_out_signal );

    return e_po_index;
  }

  bool rig_network::is_combinational() const
  {
    return true;
  }

  bool rig_network::is_constant( node const& n ) const
  {
    return n == 0;
  }

  bool rig_network::is_ci( node const& n ) const
  {
    return (_e_storage->nodes[n].children.size() == 1) && (_e_storage->nodes[n].children[0].index < num_pis() );
  }

  bool rig_network::is_pi( node const& n ) const
  {
    return (_e_storage->nodes[n].children.size() == 1) && (_e_storage->nodes[n].children[0].index < num_pis() );
  }
#pragma endregion Primary I / O and constants

#pragma region nodes and signals
  rig_network::node rig_network::get_node( rig_network::signal const& f ) const
  {
    return f.index;
  }

  rig_network::signal rig_network::make_signal( rig_network::node const& n ) const
  {
    return signal( n, 0 );
  }

  bool rig_network::is_complemented( rig_network::signal const& f ) const
  {
    return f.complement;
  }

  uint32_t rig_network::node_to_index( rig_network::node const& n ) const
  {
    return static_cast<uint32_t>( n );
  }

  rig_network::node rig_network::index_to_node( uint32_t index ) const
  {
    return index;
  }

  rig_network::node rig_network::ci_at( uint32_t index ) const
  {
    assert( index < _e_storage->inputs.size() );
    return *( _e_storage->inputs.begin() + index );
  }

  rig_network::signal rig_network::co_at( uint32_t index ) const
  {
    assert( index < _e_storage->outputs.size() );
    return *( _e_storage->outputs.begin() + index );
  }

  rig_network::node rig_network::pi_at( uint32_t index ) const
  {
    assert( index < _e_storage->inputs.size() );
    return *( _e_storage->inputs.begin() + index );
  }

  rig_network::signal rig_network::po_at( uint32_t index ) const
  {
    assert( index < _e_storage->outputs.size() );
    return *( _e_storage->outputs.begin() + index );
  }

  uint32_t rig_network::ci_index( rig_network::node const& n ) const
  {
    assert( _e_storage->nodes[n].children[0].data == _e_storage->nodes[n].children[1].data );
    return static_cast<uint32_t>( _e_storage->nodes[n].children[0].data );
  }

  uint32_t rig_network::pi_index( node const& n ) const
  {
    return static_cast<uint32_t>( _e_storage->nodes[n].children[0].data );
  }
#pragma endregion nodes and signals

#pragma region node and signal iterators
  template<typename Fn>
  void rig_network::foreach_node( Fn&& fn ) const
  {
    auto r = range<uint64_t>( _e_storage->nodes.size() );
    detail::foreach_element( r.begin(), r.end(), fn );
  }

  template<typename Fn>
  void rig_network::foreach_ci( Fn&& fn ) const
  {
    detail::foreach_element( _e_storage->inputs.begin(), _e_storage->inputs.end(), fn );
  }

  template<typename Fn>
  void rig_network::foreach_co( Fn&& fn ) const
  {
    using IteratorType = decltype( _e_storage->outputs.begin() );
    detail::foreach_element_transform<IteratorType, uint32_t>(
        _e_storage->outputs.begin(), _e_storage->outputs.end(), []( auto o ) { return o.data; }, fn );
  }

  template<typename Fn>
  void rig_network::foreach_pi( Fn&& fn ) const
  {
    detail::foreach_element( _e_storage->inputs.begin(), _e_storage->inputs.end(), fn );
  }

  template<typename Fn>
  void rig_network::foreach_po( Fn&& fn ) const
  {
    using IteratorType = decltype( _e_storage->outputs.begin() );
    detail::foreach_element_transform<IteratorType, uint32_t>(
        _e_storage->outputs.begin(), _e_storage->outputs.end(), []( auto o ) { return o.data; }, fn );
  }

  template<typename Fn>
  void rig_network::foreach_gate( Fn&& fn ) const
  {
    auto r = range<uint64_t>( 2u, _e_storage->nodes.size() ); /* start from 2 to avoid constants */
    detail::foreach_element_if(
        r.begin(), r.end(),
        [this]( auto n ) { return !is_ci( n ); },
        fn );
  }

  template<typename Fn>
  void rig_network::foreach_fanin( node const& n, Fn&& fn ) const
  {
    if ( n == 0 || is_ci( n ) )
      return;

    using IteratorType = decltype( _e_storage->outputs.begin() );
    detail::foreach_element_transform<IteratorType, uint32_t>(
        _e_storage->nodes[n].children.begin(), _e_storage->nodes[n].children.end(), []( auto f ) { return f.index; }, fn );
  }
#pragma endregion node and signal iterators

#pragma region unary functions
  rig_network::signal rig_network::create_buf( signal const& f )
  {
  //  auto [e_signal, is_new] = _create_known_node( { f }, e_func_t::e_BUF );
  //  _e_storage->nodes[e_signal.index].twin = _e_storage->nodes[f.index].twin;
  // return e_signal;
    return f;
  }

  rig_network::signal rig_network::create_not( signal const& f )
  {
    //auto [e_signal, is_new] = _create_known_node( { f }, e_func_t::e_BUF ^ 0x1 );
    //_e_storage->nodes[e_signal.index].twin = _e_storage->nodes[f.index].twin;
    //return e_signal;
    return !f;
  }

  bool rig_network::is_buf( node const& n )
  {
    _e_storage->nodes[n].func == e_func_t::e_BUF;
  }

  bool rig_network::is_not( node const& n )
  {
    _e_storage->nodes[n].func == (e_func_t::e_BUF ^ 0x1);
  }
#pragma endregion unary functions

rig_network::i_signal_t rig_network::i_create_and( i_signal_t a, i_signal_t b )
{
  /* order inputs */
  if ( a.index > b.index )
  {
    std::swap( a, b );
  }

  /* trivial cases */
  if ( a.index == b.index )
  {
    return ( a.complement == b.complement ) ? a : get_constant( false );
  }
  else if ( a.index == 0 )
  {
    return a.complement ? b : get_constant( false );
  }

  std::shared_ptr<i_storage_t>::element_type::node_type node;
  node.children[0] = a;
  node.children[1] = b;

  /* structural hashing */
  const auto it = _i_storage->hash.find( node );
  if ( it != _i_storage->hash.end() )
  {
    assert( _i_storage->nodes[it->second].nfos > 0 );
    return { it->second, 0 };
  }

  const auto index = _i_storage->get_index();

  if ( index >= .9 * _i_storage->nodes.capacity() )
  {
    _i_storage->nodes.reserve( static_cast<uint64_t>( 3.1415f * index ) );
    _i_storage->hash.reserve( static_cast<uint64_t>( 3.1415f * index ) );
  }

  _i_storage->nodes.push_back( node );

  _i_storage->hash[node] = index;

  /* increase ref-count to children */
  _i_storage->nodes[a.index].nfos++;
  _i_storage->nodes[b.index].nfos++;

  return { index, 0 };
}

#pragma region binary functions
  rig_network::signal rig_network::create_and( signal const& a, signal const& b )
  {
    auto [e_signal, is_new] = _create_known_node( { a, b }, e_func_t::e_AND );
    i_signal_t twin_a = _e_storage->nodes[get_node( a )].twin;
    i_signal_t twin_b = _e_storage->nodes[get_node( b )].twin;
    i_signal_t i_signal;
    if( is_new )
    {
      i_signal = i_create_and( { twin_a.index, twin_a.complement ^ a.complement }, { twin_b.index, twin_b.complement ^ b.complement } );
    }
  }

  rig_network::signal rig_network::create_nand( signal const& a, signal const& b )
  {
    return _create_node( { a, b }, e_func_t::e_AND ^ 0x1 );
  }

  rig_network::signal rig_network::create_or( signal const& a, signal const& b )
  {
    return _create_node( { a, b }, e_func_t::e_OR );
  }

  rig_network::signal rig_network::create_nor( signal const& a, signal const& b )
  {
    return _create_node( { a, b }, e_func_t::e_OR ^ 0x1 );
  }

  rig_network::signal rig_network::create_lt( signal const& a, signal const& b )
  {
    return _create_node( { a, b }, e_func_t::e_LT );
  }

  rig_network::signal rig_network::create_ge( signal const& a, signal const& b )
  {
    return _create_node( { a, b }, e_func_t::e_LT ^ 0x1 );
  }

  rig_network::signal rig_network::create_gt( signal const& a, signal const& b )
  {
    return _create_node( { a, b }, e_func_t::e_GT );
  }

  rig_network::signal rig_network::create_le( signal const& a, signal const& b )
  {
    return _create_node( { a, b }, e_func_t::e_GT ^ 0x1 );
  }

  rig_network::signal rig_network::create_xor( signal const& a, signal const& b )
  {
    return _create_node( { a, b }, e_func_t::e_XOR );
  }

  rig_network::signal rig_network::create_xnor( signal const& a, signal const& b )
  {
    return _create_node( { a, b }, e_func_t::e_XOR ^ 0x1 );
  }

  bool rig_network::is_and( node const& n )
  {
    return _e_storage->nodes[n].func == e_AND;
  }
  bool rig_network::is_nand( node const& n )
  {
    return _e_storage->nodes[n].func == ( e_AND ^ 0x1 );
  }
  bool rig_network::is_or( node const& n )
  {
    return _e_storage->nodes[n].func == e_OR;
  }

  bool rig_network::is_nor( node const& n )
  {
    return _e_storage->nodes[n].func == ( e_OR ^ 0x1 );
  }

  bool rig_network::is_lt( node const& n )
  {
    return _e_storage->nodes[n].func == e_LT;
  }

  bool rig_network::is_ge( node const& n )
  {
    return _e_storage->nodes[n].func == ( e_LT ^ 0x1 );
  }

  bool rig_network::is_gt( node const& n )
  {
    return _e_storage->nodes[n].func == e_GT;
  }

  bool rig_network::is_le( node const& n )
  {
    return _e_storage->nodes[n].func == ( e_GT ^ 0x1 );
  }

  bool rig_network::is_xor( node const& n )
  {
    return _e_storage->nodes[n].func == e_XOR;
  }

  bool rig_network::is_xnor( node const& n )
  {
    return _e_storage->nodes[n].func == ( e_XOR ^ 0x1 );
  }
#pragma endregion binary functions

#pragma region arbitrary function
  std::tuple<rig_network::signal, bool> rig_network::_create_known_node( std::vector<signal> const& children, uint32_t literal )
  {
    std::shared_ptr<e_storage_t>::element_type::node_type node;
    std::copy( children.begin(), children.end(), std::back_inserter( node.children ) );
    
    node.func = literal;
    const auto it1 = _e_storage->hash.find( node );
    if ( it1 != _e_storage->hash.end() )
    {
      return std::make_pair( rig_network::signal{it1->second, 0} , false );
    }

    node.func = literal ^ 0x1;
    const auto it0 = _e_storage->hash.find( node );
    if ( it0 != _e_storage->hash.end() )
    {
      return std::make_pair( rig_network::signal{it0->second, 1} , false );
    }

    const auto e_index = _e_storage->get_index();
    _e_storage->nodes.push_back( node );
    _e_storage->hash[node] = e_index;

    /* increase ref-count to children */
    for ( auto c : children )
    {
      _e_storage->nodes[c].nfos++;
    }

    return std::make_pair<signal, bool>( rig_network::signal{ e_index, 0 }, true );
  }

  rig_network::signal rig_network::_create_node( std::vector<signal> const& children, uint32_t literal )
  {
    std::shared_ptr<e_storage_t>::element_type::node_type node;
    std::copy( children.begin(), children.end(), std::back_inserter( node.children ) );
    node.func = literal;

    const auto it = _e_storage->hash.find( node );
    if ( it != _e_storage->hash.end() )
    {
      return it->second;
    }

    const auto e_index = _e_storage->get_index();
    _e_storage->nodes.push_back( node );
    _e_storage->hash[node] = e_index;

    /* increase ref-count to children */
    for ( auto c : children )
    {
      _e_storage->nodes[c].nfos++;
    }

    // synthesize
    std::array<i_signal_t, 4u> i_children {};
    for( auto i{0}; children.size(); ++i )
      i_children[i] = _e_storage->nodes[children[i].index].twin;
    auto i_signal = synthesize_twin( i_children, literal );
    _e_storage->nodes[e_index].twin = i_signal;
    _i_storage->nodes[i_signal.index].twin = { e_index, 0 };

    return { e_index, 0 };
  }

  bool rig_network::is_function( node const& n ) const
  {
    return n > 0 && !is_ci( n );
  }

  rig_network::i_signal_t rig_network::synthesize_twin( std::array<i_signal_t, 4u> & i_children, uint32_t literal )
  {
    kitty::dynamic_truth_table const& tt = _e_storage->data.cache[literal];
    return synthesize_twin_rec( i_children, tt );
  }

  rig_network::i_signal_t rig_network::synthesize_twin_rec( std::array<i_signal_t, 4u> & i_children, kitty::dynamic_truth_table const& tt )
  {
    if( i_children.size() < 4u )
    {
      return boolean_matching( i_children, tt );
    }
    return synthesize_twin_rec( i_children, tt );
  }

  rig_network::i_signal_t rig_network::boolean_matching( std::array<i_signal_t, 4u> & i_children, kitty::dynamic_truth_table const& tt )
  {
    return {0,0};
//    kitty::static_truth_table<4u> tt_s = kitty::extend_to<4u>( tt );
//
//    auto [func_npn, neg, perm] = exact_npn_canonization( tt_s );
//    auto const structures = _database.get_supergates( func_npn );
//    bool phase = ( neg >> 4 == 1 ) ? true : false;
//
//    for( auto i{0}; i<i_children.size(); ++i )
//    {
//      if( ( neg >> i ) & 0x1 == 0x1 )
//        i_children[i] = !i_children[i];
//    }
//    std::array<i_signal_t, 4> leaves;
//    for( auto i{0}; i<4; ++i )
//    {
//      leaves[i] = i_children[perm[i]];
//    }
//
//    auto & db = _database.get_database();
//    i_signal_t i_signal = {0,0};//create_twin_network( db.get_node( structures->at(0).root ), leaves );
//    
//    return phase != db.is_complemented( structures->at(0).root ) ? !i_signal : i_signal;

  }
#pragma endregion arbitrary function

#pragma region structural properties
  size_t rig_network::size() const
  {
    return static_cast<uint32_t>( _e_storage->nodes.size() );
  }

  size_t rig_network::num_cis() const
  {
    return static_cast<uint32_t>( _e_storage->inputs.size() );
  }

  size_t rig_network::num_cos() const
  {
    return static_cast<uint32_t>( _e_storage->outputs.size() );
  }

  size_t rig_network::num_pis() const
  {
    return static_cast<uint32_t>( _e_storage->inputs.size() );
  }

  size_t rig_network::num_pos() const
  {
    return static_cast<uint32_t>( _e_storage->outputs.size() );
  }

  size_t rig_network::num_gates() const
  {
    return static_cast<uint32_t>( _e_storage->nodes.size() - _e_storage->inputs.size() - 1 );
  }

  size_t rig_network::fanin_size( node const& n ) const
  {
    return static_cast<uint32_t>( _e_storage->nodes[n].children.size() );
  }

  size_t rig_network::fanout_size( node const& n ) const
  {
    return _e_storage->nodes[n].nfos & UINT32_C( 0x7FFFFFFF );
  }

  size_t rig_network::incr_fanout_size( node const& n ) const
  {
    return _e_storage->nodes[n].nfos++ & UINT32_C( 0x7FFFFFFF );
  }

  size_t rig_network::decr_fanout_size( node const& n ) const
  {
    return --_e_storage->nodes[n].nfos & UINT32_C( 0x7FFFFFFF );
  }

#pragma endregion structural properties

#pragma region functional properties
  kitty::dynamic_truth_table rig_network::node_function( const node& n ) const
  {
    return _e_storage->data.cache[_e_storage->nodes[n].func];
  }
#pragma endregion functional properties

} // namespace rils

} // namespace mockturtle