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
#include <kitty/constructors.hpp>
#include "storage.hpp"

namespace mockturtle
{

namespace rils
{

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
  using e_storage_t = storage<e_gate_t, e_data_t>;

  using i_data_t = empty_storage_data;
  struct i_hash_t;
  using i_storage_t = storage<i_gate_t, i_data_t, i_hash_t>;
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

  #pragma region constructors
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

#pragma region structural properties
public:
  auto size() const;
  auto num_cis() const;
  auto num_cos() const;
  auto num_pis() const;
  auto num_pos() const;
  auto num_gates() const;
  uint32_t fanin_size( node const& n ) const;
  uint32_t fanout_size( node const& n ) const;
  uint32_t incr_fanout_size( node const& n ) const;
  uint32_t decr_fanout_size( node const& n ) const;
  bool is_function( node const& ) const;
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

    /* reserve some truth tables for nodes */
    kitty::dynamic_truth_table tt_zero( 0 );
    _e_storage->data.cache.insert( tt_zero );
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
    const auto e_index = _e_storage->nodes.size();
    _e_storage->nodes.emplace_back();
    _e_storage->inputs.emplace_back( e_index );
    _e_storage->nodes[e_index].func = 1;
    
    const auto i_index = _i_storage->nodes.size();
    _i_storage->nodes.emplace_back();
    _i_storage->inputs.emplace_back( i_index );
    _i_storage->nodes[i_index].func = 1;

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
    return _e_storage->nodes[n].func == 1;
  }

  bool rig_network::is_pi( node const& n ) const
  {
    return _e_storage->nodes[n].func == 1;
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

  uint32_t rig_network::pi_index( rig_network::node const& n ) const
  {
    assert( _e_storage->nodes[n].children[0].data == _e_storage->nodes[n].children[1].data );
    return static_cast<uint32_t>( _e_storage->nodes[n].children[0].data );
  }
#pragma endregion nodes and signals

#pragma region structural properties
  auto rig_network::size() const
  {
    return static_cast<uint32_t>( _e_storage->nodes.size() );
  }

  auto rig_network::num_cis() const
  {
    return static_cast<uint32_t>( _e_storage->inputs.size() );
  }

  auto rig_network::num_cos() const
  {
    return static_cast<uint32_t>( _e_storage->outputs.size() );
  }

  auto rig_network::num_pis() const
  {
    return static_cast<uint32_t>( _e_storage->inputs.size() );
  }

  auto rig_network::num_pos() const
  {
    return static_cast<uint32_t>( _e_storage->outputs.size() );
  }

  auto rig_network::num_gates() const
  {
    return static_cast<uint32_t>( _e_storage->nodes.size() - _e_storage->inputs.size() - 1 );
  }

  uint32_t rig_network::fanin_size( node const& n ) const
  {
    return static_cast<uint32_t>( _e_storage->nodes[n].children.size() );
  }

  uint32_t rig_network::fanout_size( node const& n ) const
  {
    return _e_storage->nodes[n].nfos & UINT32_C( 0x7FFFFFFF );
  }

  uint32_t rig_network::incr_fanout_size( node const& n ) const
  {
    return _e_storage->nodes[n].nfos++ & UINT32_C( 0x7FFFFFFF );
  }

  uint32_t rig_network::decr_fanout_size( node const& n ) const
  {
    return --_e_storage->nodes[n].nfos & UINT32_C( 0x7FFFFFFF );
  }

  bool rig_network::is_function( node const& n ) const
  {
    return n > 0 && !is_ci( n );
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