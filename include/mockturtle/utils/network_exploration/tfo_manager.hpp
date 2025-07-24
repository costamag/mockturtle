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
 * \file tfo_manager.hpp
 * \brief Collects the TFO of a node and maintains flags on their exploration
 *
 * \author Andrea Costamagna
 */

#pragma once

#include "../node_map.hpp"
#include <limits>

namespace mockturtle
{

/*! \brief Manager for the transitive fanout
 *
 * Data structure to extract and manipulate the TFO of a network.
 */
template<typename Ntk>
class tfo_manager
{
public:
  using node_index_t = typename Ntk::node;
  using signal_t = typename Ntk::signal;
  /*! \brief Node information
   *
   * Contains flags to label a node based on:
   * - The belonging to the TFO of another node
   * - The readiness of its information ( arrival time )
   * - Whether it was seen or not in the process
   */
  struct node_info_t
  {
    node_info_t( uint64_t index, uint64_t ready, uint64_t seen )
        : index( index ), ready( ready ), seen( seen )
    {}

    node_info_t( uint64_t index )
        : node_info_t( index, 0u, 0u )
    {}

    node_info_t() = default;

    union
    {
      struct
      {
        /* index of the root node for the TFO */
        uint64_t index : 62u;
        /* set to 1 when the arrival time is up-to-date */
        uint64_t ready : 1u;
        /* set to 1 when the node has already been explored in the TFO */
        uint64_t seen : 1u;
      };
      uint64_t data;
    };
  };

public:
  tfo_manager( Ntk& ntk )
      : ntk_( ntk ),
        map_( ntk, node_info_t( 0 ) )
  {}

  /*! \brief Mark the nodes in the root's TFO */
  void init( node_index_t const& root )
  {
    root_ = root;
    mark_tfo( root );
  }

  void resize()
  {
    map_.resize();
  }

  /*! \brief Returns true if the node is in the root's TFO */
  bool belongs_to_tfo( node_index_t const& n )
  {
    return map_[n].index == root_;
  }

  bool is_marked_ready( node_index_t const& n )
  {
    return map_[n].ready > 0u;
  }

  void mark_ready( node_index_t const& n )
  {
    map_[n].ready = 1;
  }

  bool is_marked_seen( node_index_t const& n )
  {
    return ntk_.is_pi( n ) || ( map_[n].seen > 0u );
  }

  void mark_seen( node_index_t const& n )
  {
    map_[n].seen = 1;
  }

private:
  void make_tfo( node_index_t const& n )
  {
    map_[n] = node_info_t( root_ );
  }

  void mark_tfo( node_index_t const& n )
  {
    if ( belongs_to_tfo( n ) || ntk_.is_pi( n ) )
    {
      return;
    }

    make_tfo( n );

    ntk_.foreach_fanout( n, [&]( auto u ) {
      mark_tfo( u );
    } );
  }

private:
  /*! \brief Root node defining the TFO */
  uint64_t root_;
  /*! \brief Network where the TFO is analyzed */
  Ntk& ntk_;
  /*! \brief Container of the information of each node */
  incomplete_node_map<node_info_t, Ntk> map_;
};
} // namespace mockturtle