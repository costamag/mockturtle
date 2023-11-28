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
  \file cig.hpp
  \brief Configurable Inverter Graph

  \author Andrea Costamagna
*/

#pragma once

#include "../utils/truth_table_cache.hpp"
#include "storage.hpp"

namespace mockturtle
{

struct cig_storage_data;
struct cig_storage_node;
using  cig_storage = storage<cig_storage_node, cig_storage_data>;

class cig_network
{
  public:
  #pragma region Types and constructors
  static constexpr auto min_fanin_size = 1;
  static constexpr auto max_fanin_size = 32;

  using base_type = cig_network;
  using storage = std::shared_ptr<cig_storage>;
  using node_t = uint64_t;
  struct signal_t;


};

#pragma region storage

struct cig_storage_data
{
  truth_table_cache<kitty::dynamic_truth_table> cache;
};

/*! \brief cig node
 *
 * `data[0].h1`: Fan-out size
 * `data[0].h2`: Application-specific value
 * `data[1].h1`: Function literal in truth table cache
 * `data[1].h2`: Visited flags
 */
struct cig_storage_node : mixed_fanin_node<2, 1>
{
  bool operator==( cig_storage_node const& other ) const
  {
    return data[1].h1 == other.data[1].h1 && children == other.children;
  }
};

/*! \brief cig storage container

  ...
*/
using cig_storage = storage<cig_storage_node, cig_storage_data>;

struct cig_network::signal_t
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

  signal_t( cig_storage::node_type::pointer_type const& p )
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

  operator cig_storage::node_type::pointer_type() const
  {
    return { index, complement };
  }

  operator uint64_t() const
  {
    return data;
  }

#if __cplusplus > 201703L
  bool operator==( block_storage::node_type::pointer_type const& other ) const
  {
    return data == other.data;
  }
#endif
};
#pragma endregion storage

} // namespace mockturtle