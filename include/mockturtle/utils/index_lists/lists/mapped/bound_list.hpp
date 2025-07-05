#pragma once

#include "../../../../networks/mapped/bound_storage/bound_utils.hpp"
#include "../../../symm_utils.hpp"
#include <cassert>
#include <vector>

namespace mockturtle
{

/*! \brief Boolean chain of gates from a technology library.
 *
 * The inputs are associated with the literals 0, ..., num_inputs - 1.
 * The subsequent literals identify the nodes in the chain.
 */
template<bound::design_type_t DesignType>
class bound_list
{
public:
  using element_type = uint32_t;

private:
  /*! \brief Node of a mapped index list */
  struct node_t
  {
    node_t() = default;
    node_t( const node_t& ) = default;
    node_t( node_t&& ) noexcept = default;
    node_t& operator=( const node_t& ) = default;
    node_t& operator=( node_t&& ) noexcept = default;

    node_t( std::vector<element_type> const& fanins, uint32_t id )
        : fanins( fanins ), id( id )
    {}

    bool operator==( node_t const& other ) const
    {
      return fanins == other.fanins && id == other.id;
    }

    bool operator!=( node_t const& other ) const
    {
      return fanins != other.fanins || id == other.id;
    }

    std::vector<element_type> fanins; /**< literals of fanins */
    uint32_t id;                      /**< binding id */
  };

public:
#pragma region Types and Constructors

  bound_list( uint32_t num_inputs, size_t reserve_size )
      : num_inputs( num_inputs )
  {
    nodes.reserve( reserve_size );
  }

  explicit bound_list( uint32_t num_inputs )
      : bound_list( num_inputs, 10u )
  {}

  bound_list()
      : bound_list( 0 )
  {}

  ~bound_list() = default;
  bound_list( bound_list const& ) = default;
  bound_list( bound_list&& ) noexcept = default;
  bound_list& operator=( bound_list const& ) = default;
  bound_list& operator=( bound_list&& ) noexcept = default;

#pragma endregion

#pragma region Equality

  bool operator==( bound_list const& other ) const
  {
    if ( num_inputs != other.num_inputs || outputs != other.outputs || nodes.size() != other.nodes.size() )
      return false;

    for ( size_t i = 0; i < nodes.size(); ++i )
    {
      if ( nodes[i] != other.nodes[i] )
        return false;
    }
    return true;
  }

#pragma endregion

#pragma region Primary I/O and node creation

  void add_inputs( uint32_t n = 1 )
  {
    num_inputs += n;
  }

  void add_output( element_type const v )
  {
    outputs.push_back( v );
  }

  element_type pi_at( uint32_t index ) const
  {
    assert( index < num_inputs );
    return index;
  }

  element_type po_at( uint32_t index ) const
  {
    assert( index < outputs.size() );
    return outputs[index];
  }

  bool is_pi( element_type f ) const
  {
    return f < num_inputs;
  }

  /*! \brief Create a node in the list. The returned literal uniquely identifies it. */
  element_type add_gate( std::vector<element_type> const& fanins, uint32_t id )
  {
    const element_type f = static_cast<element_type>( nodes.size() + num_inputs );
    nodes.emplace_back( fanins, id );
    return f;
  }

  /*! \brief Replace in node */
  void replace_in_node( uint32_t node, uint32_t fanin, element_type other )
  {
    auto& n = nodes[node];
    n.fanins[fanin] = other;
  }

  /*! \brief Replace in node */
  void replace_output( uint32_t index, element_type other )
  {
    outputs[index] = other;
  }
#pragma endregion

#pragma region Iterators

  template<typename Fn>
  void foreach_pi( Fn&& fn ) const
  {
    for ( uint32_t i = 0; i < num_inputs; ++i )
    {
      fn( static_cast<element_type>( i ) );
    }
  }

  template<typename Fn>
  void foreach_gate( Fn&& fn ) const
  {
    for ( size_t i = 0; i < nodes.size(); ++i )
    {
      auto const& n = nodes[i];
      fn( n.fanins, n.id, i );
    }
  }

  template<typename Fn>
  void foreach_gate_rev( Fn&& fn ) const
  {
    for ( int i = static_cast<int>( nodes.size() ) - 1; i >= 0; --i )
    {
      auto const& n = nodes[i];
      fn( n.fanins, n.id, i );
    }
  }

  template<typename Fn>
  void foreach_po( Fn&& fn ) const
  {
    for ( size_t i = 0; i < outputs.size(); ++i )
    {
      fn( outputs[i], i );
    }
  }

#pragma endregion

#pragma region Structural properties

  [[nodiscard]] uint32_t num_gates() const
  {
    return static_cast<uint32_t>( nodes.size() );
  }

  [[nodiscard]] uint32_t num_pis() const
  {
    return num_inputs;
  }

  [[nodiscard]] uint32_t num_pos() const
  {
    return static_cast<uint32_t>( outputs.size() );
  }

  [[nodiscard]] uint32_t size() const
  {
    return num_inputs + static_cast<uint32_t>( nodes.size() );
  }

#pragma endregion

#pragma region Getters

  template<typename Lib>
  [[nodiscard]] double get_area( Lib const& lib ) const
  {
    double area = 0;
    foreach_gate( [&]( auto const& fanin, auto const& id, auto i ) {
      area += lib.get_area( id );
    } );
    return area;
  }

  [[nodiscard]] std::vector<node_t> get_nodes() const
  {
    return nodes;
  }

  [[nodiscard]] std::vector<element_type> get_outputs() const
  {
    return outputs;
  }

  [[nodiscard]] element_type get_num_inputs() const
  {
    return num_inputs;
  }

  [[nodiscard]] auto get_pi_index( element_type lit ) const
  {
    return static_cast<uint32_t>( lit );
  }

  [[nodiscard]] auto get_node_index( element_type lit ) const
  {
    return static_cast<uint32_t>( lit ) - num_inputs;
  }

#pragma endregion

private:
  std::vector<node_t> nodes;
  std::vector<element_type> outputs;
  element_type num_inputs = 0;
};

/*! \brief Returns the longest path from the inputs to any output */
template<bound::design_type_t DesignType, typename Lib>
std::vector<double> get_longest_paths( bound_list<DesignType>& list, Lib& library )
{
  std::vector<double> input_delays( list.num_pis(), 0 );
  std::vector<double> nodes_delays( list.num_pis() + list.num_gates(), 0 );
  list.foreach_gate_rev( [&]( auto const& fanin, auto id, auto i ) {
    for ( auto j = 0u; j < fanin.size(); ++j )
    {
      nodes_delays[fanin[j]] = std::min( nodes_delays[fanin[j]], nodes_delays[list.num_pis() + i] - library.get_max_pin_delay( id, j ) );
      if ( list.is_pi( fanin[j] ) )
        input_delays[fanin[j]] = nodes_delays[fanin[j]];
    }
  } );
  for ( auto& d : input_delays )
    d = -d;
  return input_delays;
}

/*! \brief Permutes the input variables to have the ones closest to the output first
 *
 *    *
 *   ***
 *   \**
 *    \**
 *   0 \*
 *    1 \*
 *     2
 *      3
 *
 * */
template<bound::design_type_t DesignType, typename Lib>
void time_canonize( bound_list<DesignType>& list, Lib& library, symmetries_t const& symm )
{
  std::vector<uint8_t> inputs( list.num_pis() );
  std::iota( inputs.begin(), inputs.end(), 0 );

  auto delays = get_longest_paths( list, library );

  for ( uint8_t i = 0; i < list.num_pis(); ++i )
  {
    if ( symm.has_symmetries( i ) )
    {
      uint8_t k = i;
      int j = i - 1;
      double delay = delays[inputs[i]];
      bool swapped = true;
      while ( swapped && ( j >= 0 ) )
      {
        if ( symm.symmetric( inputs[k], inputs[j] ) )
        {
          if ( delay > delays[j] )
          {
            std::swap( inputs[k], inputs[j] );
            std::swap( delays[k], delays[j] );
            k = j;
            swapped = true;
          }
          else
          {
            swapped = false;
          }
        }
        j--;
      }
    }
  }

  permutation_t perm( inputs );

  list.foreach_gate( [&]( auto& fanin, auto id, auto i ) {
    for ( auto j = 0u; j < fanin.size(); ++j )
    {
      auto const lit = fanin[j];
      if ( list.is_pi( lit ) )
      {
        list.replace_in_node( i, j, perm.inverse( lit ) );
      }
    }
  } );

  list.foreach_po( [&]( auto& lit, auto i ) {
    if ( list.is_pi( lit ) )
    {
      list.replace_output( i, perm.inverse( lit ) );
    }
  } );
}

template<bound::design_type_t DesignType>
void perm_canonize( bound_list<DesignType>& list, permutation_t const& perm )
{
  list.foreach_gate( [&]( auto& fanins, auto const& id, auto i ) {
    for ( int j = 0; j < fanins.size(); ++j )
    {
      if ( list.is_pi( fanins[j] ) )
      {
        auto k = perm.inverse( fanins[j] );
        list.replace_in_node( i, j, k );
      }
    }
  } );

  list.foreach_po( [&]( auto& output, auto i ) {
    if ( list.is_pi( output ) )
    {
      auto k = perm.inverse( output );
      list.replace_output( i, k );
    }
  } );
}

template<typename Ntk, bound::design_type_t DesignType, bool DoStrash = true>
typename Ntk::signal_t insert( Ntk& ntk,
                               std::vector<typename Ntk::signal_t> const& inputs,
                               bound_list<DesignType> const& list )
{
  std::vector<typename Ntk::signal_t> fs;

  for ( auto i = 0; i < list.num_pis(); ++i )
    fs.emplace_back( inputs[i] );

  list.foreach_gate( [&]( auto const& fanin, auto id, auto i ) {
    std::vector<typename Ntk::signal_t> children( fanin.size() );
    for ( auto i = 0u; i < fanin.size(); ++i )
    {
      children[i] = fs[fanin[i]];
    }
    fs.push_back( ntk.template create_node<DoStrash>( children, id ) );
  } );
  return fs[list.po_at( 0 )];
}

} // namespace mockturtle
