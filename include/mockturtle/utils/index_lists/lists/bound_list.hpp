#pragma once

#include <vector>
#include <cassert>

namespace mockturtle
{

/*! \brief Boolean chain of gates from a technology library.
 *
 * The inputs are associated with the literals 0, ..., num_inputs - 1.
 * The subsequent literals identify the nodes in the chain.
 */
template<typename Gate>
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
    uint32_t id;                   /**< binding id */
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

} // namespace mockturtle
