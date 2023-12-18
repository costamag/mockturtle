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

  This network assumes that buffers, inverter and splitters are cost free.
  Everything you declare apart for these has a cost.
  The network is structurally hashed for gates of the same type.
  gates of different type are not hashed together even if related by negation.
  create_and(x1,x2) != !create_nand(x1,x2).
  but naturally reate_and(x1,x2) == !create_and(x1,x2).
  Any overwriting is a representation-dependent assumption.

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
#include <kitty/print.hpp>
#include "storage.hpp"
#include <kitty/kitty.hpp>

#include <algorithm>
#include <memory>
#include <list>

namespace mockturtle
{

namespace rils
{

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
  e_PI    = 1u,
  e_BUF   = 2u,
  e_AND   = 4u,
  e_OR    = 6u,
  e_LT    = 8u,
  e_GT    = 10u,
  e_XOR   = 12u,
  e_MAJ   = 14u,
  e_ITE   = 16u,
  e_XOR3  = 18u,
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

    std::vector<pointer_type> children;

    /*! \brief number of fanouts */
    uint32_t nfos{0};
    /*! \brief id of the functionality stored in the tt-cache */
    uint32_t func{0};
    /*! \brief application specific value */
    uint32_t value{0};
    /*! \brief visited flag 1 visited */
    uint32_t visited{0};
    /*! \brief aig signal */
    aig_network::signal twin;

    bool operator==( e_gate_t const& other ) const
    {
      return (func == other.func) && (children == other.children);
    }
  };

  using e_storage_t = smart_storage<e_gate_t, e_data_t>;

#pragma endregion utils

class rig_network
{

  #pragma region types
  public:
    using e_node_t = uint64_t;
    using e_signal_t = signal_t<e_gate_t::pointer_type>;

    static constexpr auto min_fanin_size = 1;
    static constexpr auto max_fanin_size = 32;
    using base_type = rig_network;
    using storage = std::shared_ptr<e_storage_t>;
    using node = e_node_t;
    using signal = e_signal_t;
  #pragma endregion types

  #pragma region constructor
  public:
    rig_network();
    rig_network( std::shared_ptr<e_storage_t> );

  protected:
    inline void _init();
  public:
    rig_network clone() const;
  #pragma endregion constructor

  #pragma region Primary I / O and constants
  public:
    signal get_constant( bool ) const;
    bool constant_value( node const&) const;
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
    signal create_and( signal, signal );
    signal create_nand( signal, signal );
    signal create_or( signal, signal );
    signal create_nor( signal, signal );
    signal create_lt( signal, signal );
    signal create_ge( signal, signal );
    signal create_gt( signal, signal );
    signal create_le( signal, signal );
    signal create_xor( signal, signal );
    signal create_xnor( signal, signal );
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

  #pragma region ternary functions
    signal create_ite( signal, signal, signal );
    signal create_xor3( signal, signal, signal );
    signal create_maj( signal, signal, signal );
  #pragma endregion ternary functions

  #pragma region arbitrary function
    signal clone_node( rig_network const&, node const&, std::vector<signal> const& );
    signal create_node( std::vector<signal> , kitty::dynamic_truth_table );
    signal create_node_in_cloning( std::vector<signal> , kitty::dynamic_truth_table const& );
    signal _create_node( std::vector<signal> const&, uint32_t );
    bool is_function( node const& ) const;

    aig_network::signal synthesize_twin( std::vector<signal> const&, uint32_t );
    aig_network::signal synthesize_twin_rec( std::vector<aig_network::signal>, kitty::dynamic_truth_table const& );

    void print_aig( signal ) const;
    void print_aig_rec( aig_network::node, std::vector<aig_network::signal> const& ) const;
  #pragma endregion arbitrary function

  #pragma region restructuring
    inline bool is_dead( node const& ) const;
    void take_out_node( node const& );
    void replace_in_outputs( node const&, signal const& );
    std::optional<std::pair<node, signal>> replace_in_node( node const&, node const&, signal );
    void replace_in_node_no_restrash( node const&, node const&, signal );
    void revive_node( node const& );
    void substitute_node( node const&, signal const& );
    void substitute_node_no_restrash( node const&, signal const& );
    void substitute_nodes( std::list<std::pair<node, signal>> );
    void normalize_node( e_gate_t & );
    void order_inputs( std::vector<signal> &, kitty::dynamic_truth_table & );
    void constants_propagation( std::vector<signal> &, kitty::dynamic_truth_table & );

  #pragma endregion restructuring

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
  #pragma endregion structural properties

  #pragma region functional properties
  public:
    kitty::dynamic_truth_table node_function( const node& ) const;
  #pragma endregion functional properties

  #pragma region simulation properties
    template<typename Iterator>
    iterates_over_t<Iterator, bool>
    compute( node const& n, Iterator, Iterator ) const;
    template<typename Iterator>
    iterates_over_truth_table_t<Iterator>
    compute( node const&, Iterator, Iterator ) const;
    template<typename Iterator>
    void compute( node const&, kitty::partial_truth_table&, Iterator, Iterator ) const;
    template<typename TT>
    TT compute( node , std::vector<TT> const & ) const;
    template<typename TT>
    TT compute_rec( aig_network::node, std::vector<aig_network::signal> const&, std::vector<TT> const & ) const;
  #pragma endregion simulation properties

  #pragma region application specific value
    void clear_values() const;
    uint32_t value( node const& ) const;
    void set_value( node const&, uint32_t ) const;
    uint32_t incr_value( node const& n ) const;
    uint32_t decr_value( node const& n ) const;
  #pragma endregion application specific value

  #pragma region visited flags
    void clear_visited() const;
    uint32_t visited( node const& n ) const;
    void set_visited( node const& n, uint32_t v ) const;
    uint32_t trav_id() const;
    void incr_trav_id() const;
  #pragma endregion visited flags
    void print();
  #pragma region general methods
    auto& events() const;
  #pragma endregion general methods

public:
  std::shared_ptr<e_storage_t> _e_storage;
  /* complete AIG database */
  std::shared_ptr<network_events<base_type>> _events;
  aig_network _aig;

};


#pragma region constructors
  /*! \brief Network constructor

  * Construct the network using the init function
  */
  rig_network::rig_network()
    : _e_storage( std::make_shared<e_storage_t>() ), _events( std::make_shared<decltype( _events )::element_type>() )
  {
    _init();
  }

  rig_network::rig_network( std::shared_ptr<e_storage_t> e_storage_ptr )
      : _e_storage( e_storage_ptr ), _events( std::make_shared<decltype( _events )::element_type>() )
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

    /* constant node : #1 in the cache => lit = 2 */
    static uint64_t _not = 0x1;
    kitty::dynamic_truth_table tt_not( 1 );
    kitty::create_from_words( tt_not, &_not, &_not + 1 );
    _e_storage->data.cache.insert( tt_not );

    /* constant node : #2 in the cache => lit = 4 */
    static uint64_t _and = 0x8;
    kitty::dynamic_truth_table tt_and( 2 );
    kitty::create_from_words( tt_and, &_and, &_and + 1 );
    _e_storage->data.cache.insert( tt_and );

    /* constant node : #3 in the cache => lit = 6 */
    static uint64_t _or = 0xe;
    kitty::dynamic_truth_table tt_or( 2 );
    kitty::create_from_words( tt_or, &_or, &_or + 1 );
    _e_storage->data.cache.insert( tt_or );

    /* constant node : #4 in the cache => lit = 8 */
    static uint64_t _lt = 0x2;
    kitty::dynamic_truth_table tt_lt( 2 );
    kitty::create_from_words( tt_lt, &_lt, &_lt + 1 );
    _e_storage->data.cache.insert( tt_lt );

    /* constant node : #5 in the cache => lit = 10 */
    static uint64_t _gt = 0x4;
    kitty::dynamic_truth_table tt_gt( 2 );
    kitty::create_from_words( tt_gt, &_gt, &_gt + 1 );
    _e_storage->data.cache.insert( tt_gt );

    /* constant node : #6 in the cache => lit = 12 */
    static uint64_t _xor = 0x6;
    kitty::dynamic_truth_table tt_xor( 2 );
    kitty::create_from_words( tt_xor, &_xor, &_xor + 1 );
    _e_storage->data.cache.insert( tt_xor );

    /* constant node : #7 in the cache => lit = 14 */
    static uint64_t _maj = 0xe8;
    kitty::dynamic_truth_table tt_maj( 3 );
    kitty::create_from_words( tt_maj, &_maj, &_maj + 1 );
    _e_storage->data.cache.insert( tt_maj );

    /* constant node : #8 in the cache => lit = 16 */
    static uint64_t _ite = 0xd8;
    kitty::dynamic_truth_table tt_ite( 3 );
    kitty::create_from_words( tt_ite, &_ite, &_ite + 1 );
    _e_storage->data.cache.insert( tt_ite );

    /* constant node : #9 in the cache => lit = 18 */
    static uint64_t _xor3 = 0x96;
    kitty::dynamic_truth_table tt_xor3( 3 );
    kitty::create_from_words( tt_xor3, &_xor3, &_xor3 + 1 );
    _e_storage->data.cache.insert( tt_xor3 );

    /* truth tables for constants */
    _e_storage->nodes[0].func = 0;

    for( uint32_t i{0}; i<32u; ++i )
      _aig.create_pi();

  }

  rig_network rig_network::clone() const
  {
    return { std::make_shared<e_storage_t>( *_e_storage ) };
  }

  #pragma endregion constructors

#pragma region Primary I / O and constants
  rig_network::signal rig_network::get_constant( bool value = false ) const
  {
    return { 0, static_cast<uint64_t>( value ? 1 : 0 ) };
  }

  bool rig_network::constant_value( node const& n ) const
  {
    return n == 0;
  }

  rig_network::signal rig_network::create_pi()
  {
    const auto e_index = _e_storage->get_index();
    auto& e_node = _e_storage->nodes.emplace_back();
    e_node.children.emplace_back( static_cast<uint64_t>(_e_storage->inputs.size()) );
    _e_storage->inputs.push_back( e_index );
    _e_storage->nodes[e_index].func = e_func_t::e_PI;
    
    return { e_index, 0 };
  }

  uint32_t rig_network::create_po( rig_network::signal const& e_signal )
  {
    /* increase ref-count to children */
    _e_storage->nodes[e_signal.index].nfos++;
    auto const e_po_index = _e_storage->outputs.size();
    _e_storage->outputs.emplace_back( e_signal.index, e_signal.complement );

    return static_cast<uint32_t>( e_po_index );
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
    return (_e_storage->nodes[n].func == 1) ;//&& (_e_storage->nodes[n].children[0].index < num_pis() );
  }

  bool rig_network::is_pi( node const& n ) const
  {
    return (_e_storage->nodes[n].func == 1) ;//&& (_e_storage->nodes[n].children[0].index < num_pis() );
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
    detail::foreach_element_if(
        r.begin(), r.end(),
        [this]( auto n ) { return !is_dead( n ); },
        fn );
  }

  template<typename Fn>
  void rig_network::foreach_ci( Fn&& fn ) const
  {
    detail::foreach_element( _e_storage->inputs.begin(), _e_storage->inputs.end(), fn );
  }

  template<typename Fn>
  void rig_network::foreach_co( Fn&& fn ) const
  {
    detail::foreach_element( _e_storage->outputs.begin(), _e_storage->outputs.end(), fn );
  }

  template<typename Fn>
  void rig_network::foreach_pi( Fn&& fn ) const
  {
    detail::foreach_element( _e_storage->inputs.begin(), _e_storage->inputs.end(), fn );
  }

  template<typename Fn>
  void rig_network::foreach_po( Fn&& fn ) const
  {
    detail::foreach_element( _e_storage->outputs.begin(), _e_storage->outputs.end(), fn );
  }

  template<typename Fn>
  void rig_network::foreach_gate( Fn&& fn ) const
  {
    auto r = range<uint64_t>( 1u, _e_storage->nodes.size() ); /* start from 1 to avoid constant */
    detail::foreach_element_if(
        r.begin(), r.end(),
        [this]( auto n ) { return !is_ci( n ) && !is_dead( n ); },
        fn );
  }

  template<typename Fn>
  void rig_network::foreach_fanin( node const& n, Fn&& fn ) const
  {
    if ( n == 0 || is_ci( n ) )
      return;

    detail::foreach_element( _e_storage->nodes[n].children.begin(), _e_storage->nodes[n].children.end(), fn );

  }
#pragma endregion node and signal iterators

#pragma region unary functions
  rig_network::signal rig_network::create_buf( signal const& f )
  {
    return f;
  }

  rig_network::signal rig_network::create_not( signal const& f )
  {
    return !f;
  }

  bool rig_network::is_buf( node const& n )
  {
    return false;//_e_storage->nodes[n].func == e_func_t::e_BUF;
  }

  bool rig_network::is_not( node const& n )
  {
    return false; //_e_storage->nodes[n].func == (e_func_t::e_BUF ^ 0x1);
  }
#pragma endregion unary functions

#pragma region binary functions
  rig_network::signal rig_network::create_and( signal a, signal b )
  {
    /* order inputs */
    if ( a.index > b.index )
    {
      std::swap( a, b );
    }

    if ( a.index == b.index )
    {
      return ( a.complement == b.complement ) ? a : get_constant( false );
    }
    else if ( a.index == 0 )
    {
      return a.complement ? b : get_constant( false );
    }

    return _create_node( { a, b }, e_func_t::e_AND );
  }

  rig_network::signal rig_network::create_nand( signal a, signal b )
  {
    /* order inputs */
    if ( a.index > b.index )
    {
      std::swap( a, b );
    }

    if ( a.index == b.index )
    {
      return ( a.complement == b.complement ) ? !a : get_constant( true );
    }
    else if ( a.index == 0 )
    {
      return a.complement ? !b : get_constant( true );
    }

    return _create_node( { a, b }, (e_func_t::e_AND ^ 1u ) );
  }

  rig_network::signal rig_network::create_or( signal a, signal b )
  {
    /* order inputs */
    if ( a.index > b.index )
    {
      std::swap( a, b );
    }

    if ( a.index == b.index )
    {
      return ( a.complement == b.complement ) ? a : get_constant( true );
    }
    else if ( a.index == 0 )
    {
      return a.complement ? get_constant( true ) : b;
    }

    return _create_node( { a, b }, e_func_t::e_OR );
  }

  rig_network::signal rig_network::create_nor( signal a, signal b )
  {
    /* order inputs */
    if ( a.index > b.index )
    {
      std::swap( a, b );
    } 

    if ( a.index == b.index )
    {
      return ( a.complement == b.complement ) ? !a : get_constant( false );
    }
    else if ( a.index == 0 )
    {
      return a.complement ? get_constant( false ) : !b;
    }

    return _create_node( { a, b }, e_func_t::e_OR ^ 1u );
  }

  rig_network::signal rig_network::create_lt( signal a, signal b )
  {
    if ( a.index == b.index )
    {
      return ( a.complement == b.complement ) ? get_constant( false ) : b;
    }
    else if ( a.index == 0 )
    {
      return a.complement ? get_constant( false ) : b;
    }
    else if( b.index == 0 )
    {
      return b.complement ? !a : get_constant( false );
    }

    return _create_node( { a, b }, e_func_t::e_LT );
  }

  rig_network::signal rig_network::create_ge( signal a, signal b )
  {
    if ( a.index == b.index )
    {
      return ( a.complement == b.complement ) ? get_constant( true ) : !b;
    }
    else if ( a.index == 0 )
    {
      return a.complement ? get_constant( true ) : !b;
    }
    else if( b.index == 0 )
    {
      return b.complement ? a : get_constant( true );
    }
    return _create_node( { a, b }, e_func_t::e_LT ^ 1u );
  }

  rig_network::signal rig_network::create_gt( signal a, signal b )
  {
    if ( a.index == b.index )
    {
      return ( a.complement == b.complement ) ? get_constant( false ) : a;
    }
    else if ( a.index == 0 )
    {
      return a.complement ? !b : get_constant( false );
    }
    else if( b.index == 0 )
    {
      return b.complement ? get_constant( false ) : a;
    }

    return _create_node( { a, b }, e_func_t::e_GT );
  }

  rig_network::signal rig_network::create_le( signal a, signal b )
  {
    if ( a.index == b.index )
    {
      return ( a.complement == b.complement ) ? get_constant( true ) : !a;
    }
    else if ( a.index == 0 )
    {
      return a.complement ? b : get_constant( true );
    }
    else if( b.index == 0 )
    {
      return b.complement ? get_constant( true ) : !a;
    }

    return _create_node( { a, b }, e_func_t::e_GT ^ 1u );
  }

  rig_network::signal rig_network::create_xor( signal a, signal b )
  {
    /* order inputs */
    if ( a.index < b.index )
    {
      std::swap( a, b );
    }

    bool f_compl = a.complement != b.complement;
    a.complement = b.complement = false;

    if ( a.index == b.index )
    {
      return get_constant( f_compl );
    }
    else if ( b.index == 0 )
    {
      return a ^ f_compl;
    }

    return _create_node( { a, b }, e_func_t::e_XOR ) ^ f_compl;
  }

  rig_network::signal rig_network::create_xnor( signal a, signal b )
  {
    /* order inputs */
    if ( a.index < b.index )
    {
      std::swap( a, b );
    }

    bool f_compl = a.complement != b.complement;
    a.complement = b.complement = false;

    if ( a.index == b.index )
    {
      return !get_constant( f_compl );
    }
    else if ( b.index == 0 )
    {
      return !( a ^ f_compl );
    }

    return _create_node( { a, b }, e_func_t::e_XOR ^ 1u );
  }

  bool rig_network::is_and( node const& n )
  {
    return _e_storage->nodes[n].func == e_AND;
  }
  bool rig_network::is_nand( node const& n )
  {
    return ( _e_storage->nodes[n].func == ( e_AND ^ 0x1 ) );
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

#pragma region ternary functions
  rig_network::signal rig_network::create_maj( signal a, signal b, signal c )
  {
    /* order inputs */
    if ( a.index > b.index )
    {
      std::swap( a, b );
      if ( b.index > c.index )
        std::swap( b, c );
      if ( a.index > b.index )
        std::swap( a, b );
    }
    else
    {
      if ( b.index > c.index )
        std::swap( b, c );
      if ( a.index > b.index )
        std::swap( a, b );
    }

    /* trivial cases */
    if ( a.index == b.index )
    {
      return ( a.complement == b.complement ) ? a : c;
    }
    else if ( b.index == c.index )
    {
      return ( b.complement == c.complement ) ? b : a;
    }

    /*  complemented edges minimization */
    auto node_complement = false;
    if ( static_cast<unsigned>( a.complement ) + static_cast<unsigned>( b.complement ) +
             static_cast<unsigned>( c.complement ) >=
         2u )
    {
      node_complement = true;
      a.complement = !a.complement;
      b.complement = !b.complement;
      c.complement = !c.complement;
    }
    return _create_node( { a, b, c }, e_func_t::e_MAJ ) ^ node_complement;
  }

  rig_network::signal rig_network::create_ite( signal x, signal cond1, signal cond0 )
  {
    bool complement{false};
    if( cond1.index > cond0.index )
    {
      std::swap( cond1, cond0 );
      complement = true;
    }
    return _create_node( { x, cond1, cond0 }, e_func_t::e_ITE );
  }

  rig_network::signal rig_network::create_xor3( signal a, signal b, signal c )
  {
    /* order inputs */
    if ( a.index > b.index )
    {
      std::swap( a, b );
      if ( b.index > c.index )
        std::swap( b, c );
      if ( a.index > b.index )
        std::swap( a, b );
    }
    else
    {
      if ( b.index > c.index )
        std::swap( b, c );
      if ( a.index > b.index )
        std::swap( a, b );
    }

    /* trivial cases */
    if ( a.index == b.index )
    {
      return ( a.complement == b.complement ) ? c : !c;
    }
    else if ( b.index == c.index )
    {
      return ( b.complement == c.complement ) ? a : !a;
    }
    else if ( a.index == c.index )
    {
      return ( a.complement == c.complement ) ? b : !b;
    }

    /*  complemented edges minimization */
    auto complement = false;
    if ( static_cast<unsigned>( a.complement ) + static_cast<unsigned>( b.complement ) +
             static_cast<unsigned>( c.complement ) >=
         2u )
    {
      complement = true;
      a.complement = !a.complement;
      b.complement = !b.complement;
      c.complement = !c.complement;
    }
    return _create_node( { a, b, c }, e_func_t::e_XOR3 );
  }

#pragma endregion ternary functions

#pragma region arbitrary function

  void rig_network::order_inputs( std::vector<signal> & inputs, kitty::dynamic_truth_table & function )
  {
    std::vector<std::pair<signal, uint32_t>> sorted;
    for( int i{0}; i<inputs.size(); ++i )
      sorted.push_back( std::make_pair( inputs[i], i ) );
    std::sort( sorted.begin(), sorted.end(), [](std::pair<signal,uint32_t> a, std::pair<signal,uint32_t> b)
                                  {
                                      return a.first < b.first;
                                  } );
    std::vector<uint32_t> perm;
    inputs.clear();
    for( int i{0}; i<sorted.size(); ++i )
    {
      perm.push_back(sorted[i].second);
      inputs.push_back( sorted[i].first );
    }

    auto tt_new = function.construct();
    for( int m{0}; m<function.num_bits(); ++m )
    {
      uint32_t p{0};
      for( int v{0}; v<function.num_vars(); ++v )
      {
        p |= ( ( m >> perm[v] ) & 1u ) << v;
      }
      kitty::get_bit( function, m ) ? kitty::set_bit( tt_new, p ) : kitty::clear_bit( tt_new, p );
    }
    function = tt_new;
  }

  void rig_network::constants_propagation( std::vector<signal> & inputs, kitty::dynamic_truth_table & function )
  {
    printf("|S|=%d\n", inputs.size());
    kitty::print_binary( function );
    printf("\n");

    for( int iVar{0}; iVar<inputs.size(); ++iVar )
    {
      if ( is_constant( inputs[iVar] ) )
      {
        if ( is_complemented( inputs[iVar] ) )
        {
          kitty::cofactor1_inplace( function, iVar );
        }
        else
        {
          kitty::cofactor0_inplace( function, iVar );
        }
      }
    }

    const auto support = kitty::min_base_inplace( function );

    printf("|S|=%d\n", support.size());

    auto new_func = kitty::shrink_to( function, static_cast<unsigned int>( support.size() ) );
    std::cout << 10 << std::endl;
    function = new_func;
    std::cout << 11 << std::endl;
    std::vector<signal> children;
    for( int iVar{inputs.size()-1}; iVar>=0; --iVar )
    {
      if( std::find( support.begin(), support.end(), iVar ) == support.end() )
        inputs.erase( inputs.begin()+iVar );//.push_back( inputs[support[iVar]] );
    }
  }

  rig_network::signal rig_network::create_node( std::vector<signal> children, kitty::dynamic_truth_table function )
  {
    if( children.size() > 0 )
    {
    /* order inputs */
      order_inputs( children, function );
    /* update in presence of trivial inputs constant nodes */
      constants_propagation( children, function );
    }

    // n-canonize

//    auto [n_repr, neg] = exact_n_canonization( function );
//    for( int iVar{0}; iVar<function.num_vars(); ++iVar )
//      children[iVar].complement ^= ( ( neg >> iVar ) & 1u );
    
    // add sorting of the variables
    if ( children.size() == 0u )
    {
    printf("children size is 0\n");
      assert( function.num_vars() == 0u );
      return get_constant( !kitty::is_const0( function ) );
    }
    printf("_create %d %d\n", children.size(), function.num_vars());
    return _create_node( children, _e_storage->data.cache.insert( function ) );
  }

  rig_network::signal rig_network::create_node_in_cloning( std::vector<signal> children, kitty::dynamic_truth_table const& function )
  {
    // add sorting of the variables
    if ( children.size() == 0u )
    {
      assert( function.num_vars() == 0u );
      return get_constant( !kitty::is_const0( function ) );
    }
    return _create_node( children, _e_storage->data.cache.insert( function ) );
  }

  rig_network::signal rig_network::clone_node( rig_network const& other, node const& source, std::vector<signal> const& children )
  {
    assert( children.size() == other._e_storage->nodes[source].children.size() );
    return create_node_in_cloning( children, other._e_storage->data.cache[other._e_storage->nodes[source].func] );
  }

  rig_network::signal rig_network::_create_node( std::vector<signal> const& children, uint32_t literal )
  {
    std::shared_ptr<e_storage_t>::element_type::node_type node;
    std::copy( children.begin(), children.end(), std::back_inserter( node.children ) );
    node.func = literal;

    const auto it = _e_storage->hash.find( node );
    if ( it != _e_storage->hash.end() )
    {
      return rig_network::signal{it->second, 0};
    }

    const auto e_index = _e_storage->get_index();

    _e_storage->nodes.push_back( node );
    _e_storage->hash[node] = e_index;

    /* increase ref-count to children */
    int i{0};
    for ( auto c : children )
    {
      _e_storage->nodes[c.index].nfos++;
    }

//    synthesize
    auto aig_signal = synthesize_twin( children, literal );
    _e_storage->nodes[e_index].twin = aig_signal;
 
    for ( auto const& fn : _events->on_add )
    {
      ( *fn )( e_index );
    }

    return rig_network::signal{ e_index, 0 };
  }

  bool rig_network::is_function( node const& n ) const
  {
    return n > 0 && !is_ci( n );
  }

  aig_network::signal rig_network::synthesize_twin( std::vector<signal> const& children, uint32_t literal )
  {
    kitty::dynamic_truth_table const& tt = _e_storage->data.cache[literal];
    size_t n_fanins = tt.num_vars();
    std::vector<aig_network::signal> aig_children;
    for( auto i{0}; i<n_fanins; ++i )
      aig_children.push_back( aig_network::signal{_aig.pi_at(i), 0} );
    return synthesize_twin_rec( aig_children, tt );
  }

  aig_network::signal rig_network::synthesize_twin_rec( std::vector<aig_network::signal> aig_children, kitty::dynamic_truth_table const& tt )
  {
    if( kitty::is_const0( tt ) )
      return aig_network::signal{ 0, 0 };
    if( kitty::is_const0( ~tt ) )
      return aig_network::signal{ 0, 1 };

    if( aig_children.size() == 1u )
      return kitty::is_normal(tt) ? aig_children[0] : !aig_children[0];

    uint32_t idx = aig_children.size()-1;
    aig_network::signal x = aig_children[idx];
    aig_children.erase( aig_children.begin() + idx );
    aig_network::signal f1 = synthesize_twin_rec( aig_children, kitty::cofactor1( tt, idx ) );
    aig_network::signal f0 = synthesize_twin_rec( aig_children, kitty::cofactor0( tt, idx ) );
    
    if( f1.index == 0 )
    {
      auto res = f1.complement ? !_aig.create_and( !x, !f0 ) : _aig.create_and( !x, f0 );
      return res;
    }
    if( f0.index == 0 )
    {
      auto res = f0.complement ? !_aig.create_and( x, !f1 ) : _aig.create_and( x, f1 );
      return res;
    }

    return _aig.create_ite( x, f1, f0 );
  }

#pragma endregion arbitrary function

#pragma region restructuring
  inline bool rig_network::is_dead( node const& n ) const
  {
    return ( _e_storage->nodes[n].nfos >> 31 ) & 1;
  }

  void rig_network::take_out_node( node const& n )
  {
    /* we cannot delete CIs, constants, or already dead nodes */
    if ( n == 0 || is_ci( n ) || is_dead( n ) )
    {
      return;
    }

    /* delete the node (ignoring its current fanout_size) */
    auto& nobj = _e_storage->nodes[n];
    nobj.nfos = UINT32_C( 0x80000000 ); /* fanout size 0, but dead */
    _e_storage->hash.erase( nobj );

    for ( auto const& fn : _events->on_delete )
    {
      ( *fn )( n );
    }

    /* if the node has been deleted, then deref fanout_size of
       fanins and try to take them out if their fanout_size become 0 */
    for ( auto i = 0u; i < nobj.children.size(); ++i )
    {
      if ( fanout_size( nobj.children[i].index ) == 0 )
      {
        continue;
      }
      if ( decr_fanout_size( nobj.children[i].index ) == 0 )
      {
        take_out_node( nobj.children[i].index );
      }
    }
  }

  void rig_network::replace_in_outputs( node const& old_node, signal const& new_signal )
  {
    if ( is_dead( old_node ) )
    {
      return;
    }

    for ( auto& output : _e_storage->outputs )
    {
      if ( output.index == old_node )
      {
        output.index = new_signal.index;
        output.weight ^= new_signal.complement;

        if ( old_node != new_signal.index )
        {
          /* increment fan-in of new node */
          _e_storage->nodes[new_signal.index].nfos++;
        }
      }
    }
  }

  std::optional<std::pair<rig_network::node, rig_network::signal>> rig_network::replace_in_node( node const& n, node const& old_node, signal new_signal )
  {  
    auto& node = _e_storage->nodes[n];
    auto oldnode = _e_storage->nodes[old_node];
    auto newsignode = _e_storage->nodes[get_node(new_signal)];

    // remember before
    const auto old_children = node.children;

    uint32_t fanin = 0u;
    while ( fanin < node.children.size() )
    {
      if ( node.children[fanin].index == old_node )
      {
        new_signal.complement ^= node.children[fanin].weight;
        break;
      }
      fanin++;
    }
    if( fanin == node.children.size() )
      return std::nullopt;

    std::vector<signal> children;
    for( uint32_t i{0}; i < _e_storage->nodes[n].children.size(); ++i )
    {
      if( i == fanin )
      {
        children.push_back( new_signal );
      }
      else
      {
        children.push_back( node.children[i] );
      }
    }
    // determine potential new children of node n
    // normalize (SORT THE VARIABLES)
    // check for trivial cases?
    kitty::dynamic_truth_table const& tt = _e_storage->data.cache[node.func];
    for( int i{0}; i < children.size(); ++i )
    {
      auto ttnew = tt;
      auto ttx = tt.construct();
      kitty::create_nth_var( ttx, i );

      if ( children[i].index == 0 )
      {
        ttnew = children[i].complement ? kitty::cofactor1( tt, i ) : kitty::cofactor0( tt, i ) ;
        children.erase( children.begin()+i );
      }
      else if ( ( i < children.size()-1 ) && ( children[i].index == children[i+1].index) )
      {
        if( children[i].complement == children[i+1].complement )
          ttnew = (ttx & kitty::cofactor1(kitty::cofactor1( tt, i ), i+1)) | ((~ttx) & kitty::cofactor0(kitty::cofactor0( tt, i ), i+1));
        else
          ttnew = (ttx & kitty::cofactor0(kitty::cofactor1( tt, i ), i+1)) | ((~ttx) & kitty::cofactor1(kitty::cofactor0( tt, i ), i+1));
        children.erase( children.begin()+i+1 );
      }
      if( kitty::is_const0( ttnew ) ) return std::make_pair( n, get_constant(false));
      if( kitty::is_const0( ~ttnew ) ) 
      {
        return std::make_pair( n, get_constant(true));
      }
      if( children.size() == 1 )
      {
        //_e_storage->nodes[children[0].index].nfos++;
        if( kitty::is_normal( ttnew ) ) // BEWARE
        {
          return std::make_pair( n, children[0] );
        }
        else
        {
          return std::make_pair( n, !children[0] );
        }
      }
    }
    // node already in hash table
    storage::element_type::node_type _hash_obj;
    for( uint32_t i{0}; i < children.size(); ++i )
    {
      _hash_obj.children.push_back( children[i] );
    }
    _hash_obj.func = node.func;
    if ( const auto it = _e_storage->hash.find( _hash_obj ); it != _e_storage->hash.end() && it->second != old_node )
    {
      return std::make_pair( n, signal( it->second, 0 ) );
    }

    // erase old node in hash table
    _e_storage->hash.erase( node );

    // insert updated node into hash table
    node.children = _hash_obj.children;
    _e_storage->hash[node] = n;

    // update the reference counter of the new signal
    _e_storage->nodes[new_signal.index].nfos++;

  //  for ( auto const& fn : _events->on_modified )
  //  {
  //    ( *fn )( n, _hash_obj.children );
  //  }

    return std::nullopt;
  }

  void rig_network::normalize_node( e_gate_t & n )
  {
    auto children = n.children;
    std::sort( children.begin(), children.end() );
    n.children = children;
  }

  void rig_network::replace_in_node_no_restrash( node const& n, node const& old_node, signal new_signal )
  {
    auto& node = _e_storage->nodes[n];

    uint32_t fanin = 0u;
    while ( fanin < node.children.size() )
    {
      if ( node.children[fanin].index == old_node )
      {
        new_signal.complement ^= node.children[fanin].weight;
        break;
      }
      fanin++;
    }
    if( fanin == node.children.size() )
      return ;

    // determine potential new children of node n
    std::vector<e_gate_t::pointer_type> children;
    for( uint32_t i{0}; i < _e_storage->nodes[n].children.size(); ++i )
    {
      if( i == fanin )
      {
        children.push_back( new_signal );
      }
      else
      {
        children.push_back( node.children[i] );
      }
    }

    std::sort( children.begin(), children.end() );//WARNING

    // don't check for trivial cases

    // remember before
    //std::vector<signal> old_children;


    // erase old node in hash table
    _e_storage->hash.erase( node );

    // insert updated node into the hash table
    node.children = children;
    if ( _e_storage->hash.find( node ) == _e_storage->hash.end() )
    {
      _e_storage->hash[node] = n;
    }

    // update the reference counter of the new signal
    _e_storage->nodes[new_signal.index].nfos++;

//    for ( auto const& fn : _events->on_modified )
//    {
//      ( *fn )( n, { old_child0, old_child1 } );
//    }
  }

  void rig_network::revive_node( node const& n )
  {
    if ( !is_dead( n ) )
      return;
    
    assert( n < _e_storage->nodes.size() );
    auto& nobj = _e_storage->nodes[n];
    nobj.nfos = UINT32_C( 0 ); /* fanout size 0, but not dead (like just created) */
    _e_storage->hash[nobj] = n;

    for ( auto const& fn : _events->on_add )
    {
      ( *fn )( n );
    }

    /* revive its children if dead, and increment their fanout_size */
    for ( auto i = 0u; i < nobj.children.size(); ++i )
    {
      if ( is_dead( nobj.children[i].index ) )
      {
        revive_node( nobj.children[i].index );
      }
      incr_fanout_size( nobj.children[i].index );
    }
  }

  void rig_network::substitute_node( node const& old_node, signal const& new_signal )
  {
    std::unordered_map<node, signal> old_to_new;
    std::stack<std::pair<node, signal>> to_substitute;
    to_substitute.push( { old_node, new_signal } );

    while ( !to_substitute.empty() )
    {
      const auto [_old, _curr] = to_substitute.top();
      to_substitute.pop();

      signal _new = _curr;
      /* find the real new node */
      if ( is_dead( get_node( _new ) ) )
      {
        auto it = old_to_new.find( get_node( _new ) );
        while ( it != old_to_new.end() )
        {
          _new = is_complemented( _new ) ? create_not( it->second ) : it->second;
          it = old_to_new.find( get_node( _new ) );
        }
      }
      /* revive */
      if ( is_dead( get_node( _new ) ) )
      {
        revive_node( get_node( _new ) );
      }

      for ( auto idx = 1u; idx < _e_storage->nodes.size(); ++idx )
      {
        if ( is_ci( idx ) || is_dead( idx ) )
        {
          continue; /* ignore CIs */
        }

        if ( const auto repl = replace_in_node( idx, _old, _new ); repl )
        {
          to_substitute.push( *repl );
        }
      }

      /* check outputs */
      replace_in_outputs( _old, _new );
      /* recursively reset old node */
      if ( _old != _new.index )
      {
        old_to_new.insert( { _old, _new } );
        take_out_node( _old );
      }
    }
  }

  void rig_network::substitute_node_no_restrash( node const& old_node, signal const& new_signal )
  {
    if ( is_dead( get_node( new_signal ) ) )
    {
      revive_node( get_node( new_signal ) );
    }

    for ( auto idx = 1u; idx < _e_storage->nodes.size(); ++idx )
    {
      if ( is_ci( idx ) || is_dead( idx ) )
        continue; /* ignore CIs and dead nodes */

      replace_in_node_no_restrash( idx, old_node, new_signal );
    }

    /* check outputs */
    replace_in_outputs( old_node, new_signal );

    /* recursively reset old node */
    if ( old_node != new_signal.index )
    {
      take_out_node( old_node );
    }
  }

  void rig_network::substitute_nodes( std::list<std::pair<node, signal>> substitutions )
  {
    auto clean_substitutions = [&]( node const& n ) {
      substitutions.erase( std::remove_if( std::begin( substitutions ), std::end( substitutions ),
                                           [&]( auto const& s ) {
                                             if ( s.first == n )
                                             {
                                               node const nn = get_node( s.second );
                                               if ( is_dead( nn ) )
                                                 return true;

                                               /* deref fanout_size of the node */
                                               if ( fanout_size( nn ) > 0 )
                                               {
                                                 decr_fanout_size( nn );
                                               }
                                               /* remove the node if it's fanout_size becomes 0 */
                                               if ( fanout_size( nn ) == 0 )
                                               {
                                                 take_out_node( nn );
                                               }
                                               /* remove substitution from list */
                                               return true;
                                             }
                                             return false; /* keep */
                                           } ),
                           std::end( substitutions ) );
    };

    /* register event to delete substitutions if their right-hand side
       nodes get deleted */
    auto clean_sub_event = _events->register_delete_event( clean_substitutions );

    /* increment fanout_size of all signals to be used in
       substitutions to ensure that they will not be deleted */
    for ( const auto& s : substitutions )
    {
      incr_fanout_size( get_node( s.second ) );
    }

    while ( !substitutions.empty() )
    {
      auto const [old_node, new_signal] = substitutions.front();
      substitutions.pop_front();

      for ( auto index = 1u; index < _e_storage->nodes.size(); ++index )
      {
        /* skip CIs and dead nodes */
        if ( is_ci( index ) || is_dead( index ) )
          continue;

        /* skip nodes that will be deleted */
        if ( std::find_if( std::begin( substitutions ), std::end( substitutions ),
                           [&index]( auto s ) { return s.first == index; } ) != std::end( substitutions ) )
          continue;

        /* replace in node */
        if ( const auto repl = replace_in_node( index, old_node, new_signal ); repl )
        {
          incr_fanout_size( get_node( repl->second ) );
          substitutions.emplace_back( *repl );
        }
      }

      /* replace in outputs */
      replace_in_outputs( old_node, new_signal );

      /* replace in substitutions */
      for ( auto& s : substitutions )
      {
        if ( get_node( s.second ) == old_node )
        {
          s.second = is_complemented( s.second ) ? !new_signal : new_signal;
          incr_fanout_size( get_node( new_signal ) );
        }
      }

      /* finally remove the node: note that we never decrement the
         fanout_size of the old_node. instead, we remove the node and
         reset its fanout_size to 0 knowing that it must be 0 after
         substituting all references. */
      assert( !is_dead( old_node ) );
      take_out_node( old_node );

      /* decrement fanout_size when released from substitution list */
      decr_fanout_size( get_node( new_signal ) );
    }

    _events->release_delete_event( clean_sub_event );
  }
#pragma endregion restructuring

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
    return static_cast<uint32_t>( _e_storage->hash.size() );
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

#pragma region simulation properties

  template<typename Iterator>
  iterates_over_t<Iterator, bool>
  rig_network::compute( node const& n, Iterator begin, Iterator end ) const
  {
    (void)end;
    uint32_t child{ 0 };
    uint32_t index{ 0 };
    while ( begin != end )
    {
      index <<= 1;
      index ^= *begin++ ? 1 : 0;
      if( is_complemented(_e_storage->nodes[n].children[child]) )
        index ^= 1;
      child++;
    }
    return kitty::get_bit( _e_storage->data.cache[_e_storage->nodes[n].func], index );
  }

    template<typename Iterator>
    iterates_over_truth_table_t<Iterator>
    rig_network::compute( node const& n, Iterator begin, Iterator end ) const
    {
      (void)end;

      assert( n != 0 && !is_ci( n ) );

      const auto nfanin = _e_storage->nodes[n].children.size();

      std::vector<typename std::iterator_traits<Iterator>::value_type> tts( begin, end );

      //assert( nfanin != 0 );
      assert( tts.size() == nfanin );

      /* resulting truth table has the same size as any of the children */

      std::vector<aig_network::signal> children;
      
      int i{0};
      foreach_fanin( n, [&]( auto const& fi ) {
        aig_network::signal i_signal = {_aig.pi_at(i), is_complemented(fi) };//is_complemented(fi) ? !_aig.pi_at() : _e_storage->nodes[get_node(fi)].twin;
        children.push_back( i_signal );
        i++;
      } );

      auto res = compute_rec( _aig.get_node(_e_storage->nodes[n].twin), children, tts );
      if( _aig.is_complemented(_e_storage->nodes[n].twin) )
        res = ~res;

      return res;
    }

    template<typename Iterator>
    void rig_network::compute( node const& n, kitty::partial_truth_table& result, Iterator begin, Iterator end ) const
    {
      static_assert( iterates_over_v<Iterator, kitty::partial_truth_table>, "begin and end have to iterate over partial_truth_tables" );

      (void)end;

      const auto nfanin = _e_storage->nodes[n].children.size();

      std::vector<typename std::iterator_traits<Iterator>::value_type> tts( begin, end );

      assert( nfanin != 0 );
      assert( tts.size() == nfanin );

      /* resulting truth table has the same size as any of the children */

      std::vector<aig_network::signal> children;
      
      int i{0};
      foreach_fanin( n, [&]( auto const& fi ) {
        aig_network::signal i_signal = {_aig.pi_at(i), is_complemented(fi) };
        children.push_back( i_signal );
        i++;
      } );

      result = compute_rec( _e_storage->nodes[n].twin.index, children, tts );
      if( _e_storage->nodes[n].twin.complement )
        result = ~result;

    }

    template<typename TT>
    TT rig_network::compute( node n, std::vector<TT> const & tts ) const
    {
      const auto nfanin = _e_storage->nodes[n].children.size();
      assert( nfanin == tts.size() );

      std::vector<aig_network::signal> children;
      int i{0};
      foreach_fanin( n, [&]( auto const& fi ) {
        aig_network::signal i_signal = {_aig.pi_at(i), is_complemented(fi) };//is_complemented(fi) ? !_aig.pi_at() : _e_storage->nodes[get_node(fi)].twin;
        children.push_back( i_signal );
        i++;
      } );

      TT res = compute_rec( _e_storage->nodes[n].twin.index, children, tts );

      if( _e_storage->nodes[n].twin.complement )
        res = ~res;

      return res;
    }

    template<typename TT>
    TT rig_network::compute_rec( aig_network::node i_node, std::vector<aig_network::signal> const& children, std::vector<TT> const & tts ) const
    {
      auto const& i_gate = _aig._storage->nodes[i_node];
      TT res;
      if( _aig.is_constant(i_node) )
        return res;

      if( _aig.is_pi(i_node) )
      {
        return _aig.is_complemented(children[_aig.pi_index(i_node)]) ? ~tts[_aig.pi_index(i_node)] : tts[_aig.pi_index(i_node)];
      }
      else
      {
        aig_network::signal const & a = i_gate.children[0];
        aig_network::signal const & b = i_gate.children[1];
        TT sim_a = a.complement ? ~compute_rec( a.index, children, tts ) : compute_rec( a.index, children, tts );
        TT sim_b = b.complement ? ~compute_rec( b.index, children, tts ) : compute_rec( b.index, children, tts );
        res =  sim_a & sim_b;
        return res;
      }
    }

    void rig_network::print_aig( signal f ) const
    {
      node n = get_node(f);
      const auto nfanin = _e_storage->nodes[n].children.size();

      std::vector<aig_network::signal> children;
      int i{0};
      foreach_fanin( n, [&]( auto const& fi ) {
        aig_network::signal i_signal = {_aig.pi_at(i), is_complemented(fi) };//is_complemented(fi) ? !_aig.pi_at() : _e_storage->nodes[get_node(fi)].twin;
        children.push_back( i_signal );
        i++;
      } );

      print_aig_rec( _aig.get_node(_e_storage->nodes[n].twin), children );

      if( _e_storage->nodes[n].twin.complement )
        printf(" invert\n");
      else
        printf(" don't invert\n");
    }

    void rig_network::print_aig_rec( aig_network::node i_node, std::vector<aig_network::signal> const& children ) const
    {
      auto const& i_gate = _aig._storage->nodes[i_node];
      if( _aig.is_constant(i_node) )
      {
        printf("[%d=%d]", i_node, 0 );
        return;
      }

      if( _aig.is_pi(i_node) )
      {
        printf( "[%d = %c%d]", i_node, children[_aig.pi_index(i_node)].complement ? '!' : ' ', children[_aig.pi_index(i_node)].index );
      }
      else
      {
        aig_network::signal const & a = i_gate.children[0];
        aig_network::signal const & b = i_gate.children[1];
        print_aig_rec( _aig.get_node(a), children );
        print_aig_rec( _aig.get_node(b), children );
        printf("[%d=(%c%d, %c%d)]", i_node, _aig.is_complemented(a) ? '!' : ' ', a.index, _aig.is_complemented(b) ? '!' : ' ', b.index );
        return;
      }
    }
  #pragma endregion simulation properties

  #pragma region application specific value
    void rig_network::clear_values() const
    {
      std::for_each( _e_storage->nodes.begin(), _e_storage->nodes.end(), []( auto& n ) { n.value = 0; } );
    }

    uint32_t rig_network::value( node const& n ) const
    {
      return _e_storage->nodes[n].value;
    }

    void rig_network::set_value( node const& n, uint32_t v ) const
    {
      _e_storage->nodes[n].value = v;
    }

    uint32_t rig_network::incr_value( node const& n ) const
    {
      return _e_storage->nodes[n].value++;
    }

    uint32_t rig_network::decr_value( node const& n ) const
    {
      return --_e_storage->nodes[n].value;
    }
  #pragma endregion application specific value

  #pragma region visited flags
    void rig_network::clear_visited() const
    {
      std::for_each( _e_storage->nodes.begin(), _e_storage->nodes.end(), []( auto& n ) { n.visited = 0; } );
    }

    uint32_t rig_network::visited( node const& n ) const
    {
      return _e_storage->nodes[n].visited;
    }

    void rig_network::set_visited( node const& n, uint32_t v ) const
    {
      _e_storage->nodes[n].visited = v;
    }

    uint32_t rig_network::trav_id() const
    {
      return _e_storage->trav_id;
    }

    void rig_network::incr_trav_id() const
    {
      ++_e_storage->trav_id;
    }
  #pragma endregion visited flags

  #pragma region general methods
    auto& rig_network::events() const
    {
      return *_events;
    }
  #pragma endregion general methods

    void rig_network::print()
    {
      printf("POs: ");
      foreach_po( [&]( auto s, auto i ) {
        printf("%c%d ", is_complemented(s) ? '!' : ' ', s.index );
      } );
      foreach_gate( [&]( auto n ) 
      { 
        printf("[%d=");
        foreach_fanin( n, [&]( auto const& fi ) {
          printf("%c%d ", is_complemented(fi) ? '!' : ' ' , fi.index );
        } );
        printf("]");
      } );
      printf("\n");
    } 

} // namespace rils

} // namespace mockturtle