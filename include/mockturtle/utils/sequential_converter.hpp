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
  \file network_converters.hpp
  \brief Convert network types

  \author Andrea Costamagna
*/

#pragma once

#include "node_map.hpp"
#include "../views/topo_view.hpp"

#include <iostream>
#include <type_traits>
#include <vector>

namespace mockturtle
{

struct network_converters_stats
{
  uint32_t num_pos{0};
  uint32_t num_pis{0};
};

namespace detail
{

template<typename Ntk>
void generate_combinational_inputs( sequential<Ntk> const& sntk, Ntk& cntk, std::vector<signal<Ntk>>& cis)
{
  /* network name */
  if constexpr ( has_get_network_name_v<sequential<Ntk>> && has_set_network_name_v<Ntk> )
  {
    cntk.set_network_name( sntk.get_network_name() );
  }

  /* PIs & PI names */
  sntk.foreach_pi( [&]( auto n ) {
    cis.push_back( cntk.create_pi() );
    if constexpr ( has_has_name_v<sequential<Ntk>> && has_get_name_v<sequential<Ntk>> && has_set_name_v<Ntk> )
    {
      auto const s = sntk.make_signal( n );
      if ( sntk.has_name( s ) )
      {
        cntk.set_name( cis.back(), sntk.get_name( s ) );
      }
      if ( sntk.has_name( !s ) )
      {
        cntk.set_name( !cis.back(), sntk.get_name( !s ) );
      }
    }
  } );

  /* ROs & RO names & register information */

  sntk.foreach_ro( [&]( auto const& n, auto i ) {
    cis.push_back( cntk.create_pi() );
    if constexpr ( has_has_name_v<sequential<Ntk>> && has_get_name_v<sequential<Ntk>> && has_set_name_v<Ntk> )
    {
      auto const s = sntk.make_signal( n );
      if ( sntk.has_name( s ) )
      {
        cntk.set_name( cis.back(), sntk.get_name( s ) );
      }
      if ( sntk.has_name( !s ) )
      {
        cntk.set_name( !cis.back(), sntk.get_name( !s ) );
      }
    }
  } );
}

template<typename Ntk, typename LeavesIterator>
void generate_combinational_nodes( sequential<Ntk> const& sntk, Ntk& cntk, LeavesIterator begin, LeavesIterator end, unordered_node_map<signal<Ntk>, sequential<Ntk>>& old_to_new )
{
  /* constants */
  old_to_new[sntk.get_constant( false )] = cntk.get_constant( false );
  if ( sntk.get_node( sntk.get_constant( true ) ) != sntk.get_node( sntk.get_constant( false ) ) )
  {
    old_to_new[sntk.get_constant( true )] = cntk.get_constant( true );
  }

  /* create inputs in the same order */
  auto it = begin;
  sntk.foreach_pi( [&]( auto node ) {
    old_to_new[node] = *it++;
  } );
  assert( it == end );
  (void)end;

  /* foreach node in topological order */
  topo_view topo{ sntk };
  topo.foreach_node( [&]( auto node ) {
    if ( sntk.is_constant( node ) || sntk.is_ci( node ) )
      return;

    /* collect children */
    std::vector<signal<Ntk>> children;
    sntk.foreach_fanin( node, [&]( auto child, auto ) {
      const auto f = old_to_new[child];
      if ( sntk.is_complemented( child ) )
      {
        children.push_back( cntk.create_not( f ) );
      }
      else
      {
        children.push_back( f );
      }
    } );

    /* clone node */
    old_to_new[node] = cntk.clone_node( sntk, node, children );

    /* copy name */
    if constexpr ( has_has_name_v<sequential<Ntk>> && has_get_name_v<sequential<Ntk>> && has_set_name_v<Ntk> )
    {
      auto const s = sntk.make_signal( node );
      if ( sntk.has_name( s ) )
      {
        cntk.set_name( old_to_new[node], sntk.get_name( s ) );
      }
      if ( sntk.has_name( !s ) )
      {
        cntk.set_name( !old_to_new[node], sntk.get_name( !s ) );
      }
    }
  } );
}

template<typename Ntk>
void generate_combinational_outputs( sequential<Ntk> const& sntk, Ntk& cntk, unordered_node_map<signal<Ntk>, sequential<Ntk>> const& old_to_new )
{
  /* POs */
  sntk.foreach_po( [&]( auto const& po ) {
    auto const f = old_to_new[po];
    auto const n = cntk.get_node( f );
    cntk.create_po( sntk.is_complemented( po ) ? cntk.create_not( f ) : f );
  } );

  sntk.foreach_ri( [&]( auto const& ri ) {
    auto const f = old_to_new[ri];
    auto const n = cntk.get_node( f );
    cntk.create_po( sntk.is_complemented( ri ) ? cntk.create_not( f ) : f );
  } );

  /* CO names */
  if constexpr ( has_has_output_name_v<sequential<Ntk>> && has_get_output_name_v<sequential<Ntk>> && has_set_output_name_v<Ntk> )
  {
    sntk.foreach_co( [&]( auto co, auto index ) {
      (void)co;
      if ( sntk.has_output_name( index ) )
      {
        cntk.set_output_name( index, sntk.get_output_name( index ) );
      }
    } );
  }
}

template<typename Ntk>
void generate_sequential_inputs( Ntk const& cntk, sequential<Ntk>& sntk, uint32_t num_pis, std::vector<signal<sequential<Ntk>>>& cis, std::vector<signal<sequential<Ntk>>>& ros)
{
  /* network name */
  if constexpr ( has_get_network_name_v<sequential<Ntk>> && has_set_network_name_v<Ntk> )
  {
    sntk.set_network_name( cntk.get_network_name() );
  }

  /* PIs & PI names */
  cntk.foreach_pi( [&]( auto n, uint32_t i ) {
    auto const s = cntk.make_signal( n );

    if( i < num_pis )
    {
      cis.push_back( sntk.create_pi() );
      if constexpr ( has_has_name_v<sequential<Ntk>> && has_get_name_v<sequential<Ntk>> && has_set_name_v<Ntk> )
      {
        if ( cntk.has_name( s ) )
        {
          sntk.set_name( cis.back(), cntk.get_name( s ) );
        }
        if ( cntk.has_name( !s ) )
        {
          sntk.set_name( !cis.back(), cntk.get_name( !s ) );
        }
      }
    }
    else
    {
      ros.push_back( sntk.create_ro() );
      if constexpr ( has_has_name_v<sequential<Ntk>> && has_get_name_v<sequential<Ntk>> && has_set_name_v<Ntk> )
      {
        if ( cntk.has_name( s ) )
        {
          sntk.set_name( ros.back(), cntk.get_name( s ) );
        }
        if ( cntk.has_name( !s ) )
        {
          sntk.set_name( !ros.back(), cntk.get_name( !s ) );
        }
      }
    }
  } );
}

template<typename Ntk>
signal<sequential<Ntk>> generate_sequential_rec( Ntk const& cntk, sequential<Ntk>& sntk, typename Ntk::node const& n, uint32_t num_pis, unordered_node_map<signal<sequential<Ntk>>, Ntk>& old_to_new )
{
  if( cntk.is_constant(n) )
  {
    return old_to_new[n];
  }
  if( cntk.is_pi(n) )
  {
    return old_to_new[n];
  }

  if( old_to_new.has(n) ) return old_to_new[n];
  return sntk.get_constant(false);

  std::vector<signal<sequential<Ntk>>> children;
  cntk.foreach_fanin( n, [&]( auto child, auto ) {
    const auto nd_child = cntk.get_node( child );
    if( old_to_new.has( nd_child ) )
    {
      const auto f = old_to_new[child];
      if ( sntk.is_complemented( child ) )
      {
        children.push_back( sntk.create_not( f ) );
      }
      else
      {
        children.push_back( f );
      }
    }
    else
    {
      auto f = generate_sequential_rec( cntk, sntk, nd_child, num_pis, old_to_new );
      if ( sntk.is_complemented( child ) )
      {
        children.push_back( sntk.create_not( f ) );
      }
      else
      {
        children.push_back( f );
      }
    }
  } );
//
//    /* clone node */
  old_to_new[n] = sntk.clone_node( cntk, n, children );
//
  /* copy name */
  if constexpr ( has_has_name_v<sequential<Ntk>> && has_get_name_v<sequential<Ntk>> && has_set_name_v<Ntk> )
  {
    auto const s = cntk.make_signal( n );
    if ( cntk.has_name( s ) )
    {
      sntk.set_name( old_to_new[n], cntk.get_name( s ) );
    }
    if ( cntk.has_name( !s ) )
    {
      sntk.set_name( !old_to_new[n], cntk.get_name( !s ) );
    }
  }

  return old_to_new[n];
//
//  return generate_sequential_rec( cntk, sntk, n, unordered_node_map<signal<sequential<Ntk>>, Ntk>& old_to_new )
}

template<typename Ntk, typename LeavesIterator>
void generate_sequential_from_outputs( Ntk const& cntk, sequential<Ntk>& sntk, uint32_t num_pis, uint32_t num_pos, LeavesIterator begin, LeavesIterator end, LeavesIterator rbegin, LeavesIterator rend, unordered_node_map<signal<sequential<Ntk>>, Ntk>& old_to_new )
{
  old_to_new[cntk.get_constant( false )] = sntk.get_constant( false );
  if ( cntk.get_node( cntk.get_constant( true ) ) != cntk.get_node( cntk.get_constant( false ) ) )
  {
    old_to_new[cntk.get_constant( true )] = sntk.get_constant( true );
  }

  /* create inputs in the same order */
  auto it = begin;
  sntk.foreach_pi( [&]( auto node ) {
    old_to_new[node] = *it++;
  } );
  assert( it == end );
  (void)end;

  auto rit = rbegin;
  sntk.foreach_ro( [&]( auto node ) {
    old_to_new[node] = *rit++;
  } );
  assert( rit == rend );
  (void)rend;

  cntk.foreach_po( [&]( auto const& po, uint32_t index ) {
    auto const n = cntk.get_node( po );
    auto const f = detail::generate_sequential_rec( cntk, sntk, n, num_pis, old_to_new );

    if( index < num_pos )
    {
      sntk.create_po( cntk.is_complemented( po ) ? sntk.create_not( f ) : f );
    }
    else // this is a register input
    {
      sntk.create_ri( cntk.is_complemented( po ) ? sntk.create_not( f ) : f );
    }
  } );
}

} // namespace detail


/*! \brief Converts a sequential network to a combinational one.
 *
 * This method transforms a sequential network into a combinational one. 
 * The network types of the source and destination network are the same.
 *
   \verbatim embed:rst

   .. note::

      This method returns the combinational version of a sequential network as a return value.  It does
      *not* modify the input network. The user must provide the converter stats, needed to restore the 
      sequential network after combinational optimization.
   \endverbatim
 *
 * **Required network functions:**
 * - `get_node`
 * - `node_to_index`
 * - `get_constant`
 * - `create_pi`
 * - `create_po`
 * - `create_not`
 * - `is_complemented`
 * - `foreach_node`
 * - `foreach_pi`
 * - `foreach_po`
 * - `clone_node`
 * - `is_pi`
 * - `is_constant`
 */
template<class Ntk>
[[nodiscard]] Ntk sequential_to_combinatorial( sequential<Ntk> const& sntk, network_converters_stats& st )
{
  static_assert( is_network_type_v<sequential<Ntk>>, "Ntk is not a network type" );
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_get_node_v<sequential<Ntk>>, "Ntk does not implement the get_node method" );
  static_assert( has_node_to_index_v<sequential<Ntk>>, "Ntk does not implement the node_to_index method" );
  static_assert( has_get_constant_v<sequential<Ntk>>, "Ntk does not implement the get_constant method" );
  static_assert( has_foreach_node_v<sequential<Ntk>>, "Ntk does not implement the foreach_node method" );
  static_assert( has_foreach_pi_v<sequential<Ntk>>, "Ntk does not implement the foreach_pi method" );
  static_assert( has_foreach_po_v<sequential<Ntk>>, "Ntk does not implement the foreach_po method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_is_constant_v<sequential<Ntk>>, "Ntk does not implement the is_constant method" );
  static_assert( has_clone_node_v<Ntk>, "Ntk does not implement the clone_node method" );
  static_assert( has_create_pi_v<Ntk>, "Ntk does not implement the create_pi method" );
  static_assert( has_create_po_v<Ntk>, "Ntk does not implement the create_po method" );
  static_assert( has_create_not_v<Ntk>, "Ntk does not implement the create_not method" );
  static_assert( has_is_complemented_v<sequential<Ntk>>, "Ntk does not implement the is_complemented method" );

  Ntk cntk;

  if constexpr ( is_crossed_network_type_v<Ntk> )
  {
    printf("[w]: Crossed networks not supported in sequential_to_combinatorial\n");
    return cntk;
  }

  st.num_pis = sntk._sequential_storage->num_pis;
  st.num_pos = sntk._sequential_storage->num_pos;

  std::vector<signal<Ntk>> cis; 

  detail::generate_combinational_inputs( sntk, cntk, cis );

  unordered_node_map<signal<Ntk>, sequential<Ntk>> old_to_new( sntk );

  detail::generate_combinational_nodes( sntk, cntk, cis.begin(), cis.end(), old_to_new );

  detail::generate_combinational_outputs( sntk, cntk, old_to_new );

  return cntk;
}

/*! \brief Converts a sequential network to a combinational one.
 *
 * This method transforms a sequential network into a combinational one. 
 * The network types of the source and destination network are the same.
 *
   \verbatim embed:rst

   .. note::

      This method returns the combinational version of a sequential network as a return value.  It does
      *not* modify the input network. The user must provide the converter stats, needed to restore the 
      sequential network after combinational optimization.
   \endverbatim
 *
 * **Required network functions:**
 * - `get_node`
 * - `node_to_index`
 * - `get_constant`
 * - `create_pi`
 * - `create_po`
 * - `create_not`
 * - `is_complemented`
 * - `foreach_node`
 * - `foreach_pi`
 * - `foreach_po`
 * - `clone_node`
 * - `is_pi`
 * - `is_constant`
 */
template<class Ntk>
[[nodiscard]] sequential<Ntk> combinatorial_to_sequential( Ntk const& cntk, network_converters_stats const& st )
{
  static_assert( is_network_type_v<sequential<Ntk>>, "Ntk is not a network type" );
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_get_node_v<sequential<Ntk>>, "Ntk does not implement the get_node method" );
  static_assert( has_node_to_index_v<sequential<Ntk>>, "Ntk does not implement the node_to_index method" );
  static_assert( has_get_constant_v<sequential<Ntk>>, "Ntk does not implement the get_constant method" );
  static_assert( has_foreach_node_v<sequential<Ntk>>, "Ntk does not implement the foreach_node method" );
  static_assert( has_foreach_pi_v<sequential<Ntk>>, "Ntk does not implement the foreach_pi method" );
  static_assert( has_foreach_po_v<sequential<Ntk>>, "Ntk does not implement the foreach_po method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_is_constant_v<sequential<Ntk>>, "Ntk does not implement the is_constant method" );
  static_assert( has_clone_node_v<Ntk>, "Ntk does not implement the clone_node method" );
  static_assert( has_create_pi_v<Ntk>, "Ntk does not implement the create_pi method" );
  static_assert( has_create_po_v<Ntk>, "Ntk does not implement the create_po method" );
  static_assert( has_create_not_v<Ntk>, "Ntk does not implement the create_not method" );
  static_assert( has_is_complemented_v<sequential<Ntk>>, "Ntk does not implement the is_complemented method" );

  sequential<Ntk> sntk;

  if constexpr ( is_crossed_network_type_v<Ntk> )
  {
    printf("[w]: Crossed networks not supported in sequential_to_combinatorial\n");
    return sntk;
  }

  std::vector<signal<sequential<Ntk>>> cis; 
  std::vector<signal<sequential<Ntk>>> ros; 

  detail::generate_sequential_inputs( cntk, sntk, st.num_pis, cis, ros );
  unordered_node_map<signal<sequential<Ntk>>, Ntk> old_to_new( cntk );
  //detail::generate_sequential_nodes( cntk, sntk, cis.begin(), cis.end(), old_to_new );
  detail::generate_sequential_from_outputs( cntk, sntk, st.num_pis, st.num_pos, cis.begin(), cis.end(), ros.begin(), ros.end(), old_to_new );
  return sntk;
}

} // namespace mockturtle