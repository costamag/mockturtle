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


/*! \brief rig luts storage
*/
struct e_data_t
{
  truth_table_cache<kitty::dynamic_truth_table> cache;
};

/*! \brief rig node
 *
 * `data[0].h1`: Fan-out size
 * `data[0].h2`: To be defined
 * `data[1].h1`: Function literal in truth table cache
 * `data[1].h2`: To be defined
 * `data[2].h1`: Twin internal node: ( i_index << 1 ) | ( i-to-e negation )
 * `data[2].h2`: To be defined
 */
struct e_gate_t : mixed_fanin_node<3, 1>
{
  bool operator==( e_gate_t const& other ) const
  {
    return data[1].h1 == other.data[1].h1 && children == other.children;
  }
};

/*! \brief Hash function for AIGs (from ABC) */
struct i_hash_t
{
  uint64_t operator()( regular_node<2, 3, 1> const& n ) const
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

/*! \brief AIG storage container

 * AIGs have nodes with fan-in 2.  We split of one bit of the index pointer to
 * store a complemented attribute.  Every node has 64-bit of additional data
 * used for the following purposes:
 *
 * `data[0].h1`: Fan-out size (we use MSB to indicate whether a node is dead)
 * `data[0].h2`: To be defined
 * `data[1].h1`: Node type ( CONST:0 PI:1 )
 * `data[1].h2`: To be defined
 * `data[2].h1`: Twin external node: ( e_index << 1 ) | ( e-to-i negation )
 * `data[2].h2`: To be defined
*/
using i_gate_t = regular_node<2, 3, 1>;
using i_data_t = empty_storage_data;
struct i_hash_t;
using i_storage_t = storage<i_gate_t, i_data_t, i_hash_t>;

class rig_network
{
  public:

  #pragma region types
    using e_node_t = uint64_t;
    using i_node_t = uint64_t;

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

    using e_signal_t = signal_t<e_gate_t::pointer_type>;
    using i_signal_t = signal_t<i_gate_t::pointer_type>;

    static constexpr auto min_fanin_size = 1;
    static constexpr auto max_fanin_size = 32;
    using base_t = rig_network;
    using node = e_node_t;
    using signal = e_signal_t;
  #pragma endregion types

  #pragma region constructors
    rig_network()
      : _e_storage( std::make_shared<e_storage_t>() ), _i_storage( std::make_shared<i_storage_t>() )
    {
      _init();
    }

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

public:
  std::shared_ptr<e_storage_t> _e_storage;
  std::shared_ptr<i_storage_t> _i_storage;

};


#pragma region constructors
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

    /* truth tables for constants */
    _e_storage->nodes[1].data[0].h1 = 0;
    _i_storage->nodes[1].data[0].h1 = 0;
    _e_storage->nodes[2].data[0].h1 = 0;
    _i_storage->nodes[2].data[0].h1 = 0;
  }
  #pragma endregion constructors

#pragma region linking
  rig_network::e_signal_t rig_network::get_e_signal( rig_network::i_signal_t const& f )
  {
    return { _e_storage->nodes[f].data[2].h1 >> 1, _e_storage->nodes[f].data[2].h1 & 1u };
  }

  rig_network::i_signal_t rig_network::get_i_signal( rig_network::e_signal_t const& f )
  {
    return { _i_storage->nodes[f].data[2].h1 >> 1, _i_storage->nodes[f].data[2].h1 & 1u };
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
    _e_storage->nodes[e_index].data[1].h1 = 1;
    
    const auto i_index = _i_storage->nodes.size();
    _i_storage->nodes.emplace_back();
    _i_storage->inputs.emplace_back( i_index );
    _e_storage->nodes[e_index].data[1].h1 = 1;

    _e_storage->nodes[e_index].data[2].h1 = i_index << 1;
    _i_storage->nodes[i_index].data[2].h1 = e_index << 1;

    return { e_index, 0 };
  }

  uint32_t rig_network::create_po( rig_network::signal const& e_signal )
  {
    /* increase ref-count to children */
    _e_storage->nodes[e_signal.index].data[0].h1++;
    auto const e_po_index = _e_storage->outputs.size();
    _e_storage->outputs.emplace_back( e_signal );
    
    i_signal_t i_signal = get_i_signal( e_signal );
    _i_storage->nodes[i_signal.index].data[0].h1++;
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
    return _e_storage->nodes[n].data[1].h1++ == 1;
  }

  bool rig_network::is_pi( node const& n ) const
  {
    return _e_storage->nodes[n].data[1].h1++ == 1;
  }
#pragma endregion Primary I / O and constants

} // namespace mockturtle