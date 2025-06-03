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
  \file lib_index_list.hpp
  \brief List of indices to represent small networks in a certain library.

  \author Andrea Costamagna
*/

#pragma once

#include <fmt/format.h>

#include <vector>

namespace mockturtle
{

/*! \brief Index list for graphs of nodes from a technology library.
 *
 * Small network represented as a list of literals. Supports standard cells in
 * an arbitrary technology library. The list has the following 32-bit unsigned
 * integer elements. The first two entries are the number of inputs `num_pis`
 * and the number of outputs `num_pos`. Afterwards, the gates are characterized
 * by specifying the fanin size, the list of the fanins literals, and the id of
 * the gate in the technology library. Contrarily to other lists, no literal is
 * reserved for the constants, as a standard cell technology library should
 * provide gates for representing the constant 0 function and the constant 1
 * function. The last `num_pos` literals contain the literals corresponding to
 * the outputs of the index list.
 */
template<typename Gate>
struct lib_index_list
{
public:
  using element_type = uint32_t;

public:
  explicit lib_index_list( uint32_t num_pis = 0 )
      : values( { num_pis } )
  {
    /* add the entry to record the number of outputs */
    values.emplace_back( 0 );
    /* add the entry to record the number of gates */
    values.emplace_back( 0 );
  }

  explicit lib_index_list( std::vector<element_type> const& values )
      : values( std::begin( values ), std::end( values ) )
  {}

  /*! \brief Move constructor.
   */
  lib_index_list( lib_index_list&& ) noexcept = default;

  /*! \brief Assignment operator.
   */
  lib_index_list& operator=( lib_index_list&& ) noexcept = default;

  /*! \brief Checker for equivalence with another lis.
   *
   * \param other List to be compared with `this` for equivalence
   * \return `true` if the lists are equivalent
   */
  bool operator==( const lib_index_list& other ) const
  {
    return values == other.values;
  }

  /*! \brief Getter for the raw information of the list.
   *
   * \return The list of literals stored in the list
   */
  std::vector<element_type> raw() const
  {
    return values;
  }

  /*! \brief Getter for the size of the raw information in the list.
   *
   * \return The number of literals defining the index list
   */
  uint64_t size() const
  {
    return values.size();
  }

  /*! \brief Getter for the number of gates in the list.
   *
   * \return Returns the number of gates in the index list
   */
  uint64_t num_gates() const
  {
    return values.at( 2 );
  }

  /*! \brief Getter for the number of inputs in the list.
   *
   * \return the number of inputs in the index list
   */
  uint64_t num_pis() const
  {
    return values.at( 0 );
  }

  /*! \brief Getter for the number of outputs in the list.
   *
   * \return The number of outputs in the index list
   */
  uint64_t num_pos() const
  {
    return values.at( 1 );
  }

  /*! \brief Iterator over the gates in the index list.
   *
   * This function allows applying a user defined function to the gates in the
   * list. Each such function is expected to take three arguments:
   * - An iterator to the first children index
   * - An iterator to the last children index
   * - A binding identifier of the gate corresponding to a standard cell output
   *
   * \param fn Function to be applied to the gate
   */
  template<typename Fn>
  void foreach_gate( Fn&& fn ) const
  {
    using iterate_type = typename std::vector<element_type>::const_iterator;
    auto i = 3u;
    while ( i < values.size() - num_pos() )
    {
      element_type const num_fanins = values[i];
      iterate_type start_iter = values.begin() + i + 1;
      iterate_type final_iter = start_iter + num_fanins;
      element_type const identifier = values[i + num_fanins + 1];
      fn( start_iter, final_iter, identifier );
      i += num_fanins + 2;
    }
  }

  /*! \brief Iterator over the output literals of the index list.
   *
   * This function allows applying a user defined function to the outputs of
   * the list. Each such function is expected to take the output literal as the
   * only argument.
   *
   * \param fn Function to be applied to the output literals
   */
  template<typename Fn>
  void foreach_po( Fn&& fn ) const
  {
    for ( auto i = values.size() - num_pos(); i < values.size(); ++i )
    {
      fn( values.at( i ) );
    }
  }

  /*! \brief Reset the index list by deleting all the information it contains.
   */
  void clear()
  {
    values.clear();
    values = { 0, 0, 0 };
  }

  /*! \brief Add new inputs to the index list.
   *
   * \param n (default 1) Number of inputs to add
   */
  void add_inputs( uint32_t n = 1u )
  {
    values.at( 0u ) += n;
  }

  /*! \brief Add a new gate to the index list.
   *
   * \param children Vector of literals of the inputs
   * \param id Binding identifier from the technology library
   * \return The literal identifying the new gate added in the list
   */
  element_type add_gate( std::vector<element_type> const& children,
                         element_type id )
  {
    /* in lib-lists constants are gates --> no offset */
    element_type const lit = num_pis() + num_gates();
    /* increment the number of gates */
    values.at( 2u ) += 1;
    /* store the fanin size for the gate */
    values.push_back( children.size() );
    /* insert the children literals */
    values.insert( values.end(), children.begin(), children.end() );
    /* insert the binding id of the gate */
    values.push_back( id );
    return lit;
  }

  /*! \brief Add a literal to the outputs of the list.
   *
   * \param lit Literal to be considered as an output
   */
  void add_output( element_type lit )
  {
    values.at( 1u ) += 1;
    values.push_back( lit );
  }

  /*! \brief Check if a literal is a PI
   *
   * The literals identifying an input are [0, ..., num_pis() - 1]
   *
   * \param lit Literal to be analyzed to see if is an input
   * \return `true` if `lit` is an input
   */
  inline bool is_pi( element_type const& lit ) const
  {
    return lit < num_pis();
  }

  /*! \brief Returns the output at a given index
   *
   * \param index Index of the desired output
   * \return The output literal of the `index`-th output
   */
  inline element_type po_at( uint32_t index ) const
  {
    if ( index >= num_pos() )
      throw std::out_of_range( "Output index out of bounds" );
    return *( values.end() - num_pos() + index );
  }

  /*! \brief Returns the input literal at a given index
   *
   * The literals identifying an input are [0, ..., num_pis() - 1]
   *
   * \param index Index of the desired input
   * \return The input literal of the `index`-th input
   */
  inline element_type pi_at( uint32_t index ) const
  {
    if ( index < num_pos() )
      throw std::out_of_range( "Input index out of bounds" );
    return index;
  }

  /*! \brief Returns the node index excluding the constants and the inputs */
  inline uint32_t get_node_index( element_type const& lit ) const
  {
    return lit - num_pis();
  }

  /*! \brief Returns index of an input literal */
  inline uint32_t get_pi_index( element_type const& lit ) const
  {
    return lit;
  }

private:
  /* raw information stored in the list  */
  std::vector<element_type> values;
};

/*! \brief Converts an lib_index_list to a string
 *
 * \param list An lib index list
 * \return A string representation of the index list
 */
template<typename Gate>
inline std::string to_index_list_string( lib_index_list<Gate> const& list )
{

  auto s = fmt::format( "{{{}, {}, {}", list.num_pis(), list.num_pos(), list.num_gates() );

  list.foreach_gate( [&]( auto start, auto end, auto id ) {
    /* fanin size */
    s += fmt::format( ", {}", end - start );
    while ( start != end )
    {
      s += fmt::format( ", {}", *start++ );
    }
  } );

  list.foreach_po( [&]( uint32_t lit ) {
    s += fmt::format( ", {}", lit );
  } );

  s += "}";

  return s;
}

} // namespace mockturtle
