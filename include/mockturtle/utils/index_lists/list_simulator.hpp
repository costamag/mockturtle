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
  \file list_simulator.hpp
  \brief Simulator engine for index lists.

  \author Andrea Costamagna
*/

#pragma once

#include "../mapping/augmented_library.hpp"
#include "index_list.hpp"

#include <algorithm>
#include <vector>

namespace mockturtle
{

/*! \brief Simulator engine for XAG, AIG, and MIG-index lists.
 *
 * This engine can be used to efficiently simulate many index lists.
 * In this context, a simulation pattern is a truth table corresponding
 * to a node’s Boolean vector under the given input assignments.
 * The simulator pre-allocates the memory necessary to store the simulation patterns
 * of an index list, and extends this memory when the evaluator is required
 * to perform Boolean evaluation of a list that is larger than the current capacity
 * of the simulator. To avoid unnecessary copies, the input simulation patterns
 * must be passed as a vector of raw pointers.
 *
 * \tparam List Type of the index list to be simulated.
 * \tparam TT Truth table type.
 *
 * \verbatim embed:rst

   Example

   .. code-block:: c++

      using list_t = large_xag_index_list;
      using truth_table_t = kitty::static_truth_table<4u>;
      list_simulator<list_t, truth_table_t> sim;
      sim( list1, inputs1 );
      sim( list2, inputs2 );
   \endverbatim
 */
template<typename List, typename TT>
class list_simulator
{
public:
  /* type of the literals of the index list */
  using element_type = typename List::element_type;

  list_simulator()
      : const0( TT().construct() )
  {
    /* The value 20 is an upper-bound of the size of most practical lists */
    sims.resize( 20u );
  }

  /*! \brief Simulate the list in topological order.
   *
   * This method updates the internal state of the simulator by
   * storing in `sims` the simulation patterns of the nodes in the list.
   *
   * \param list Index list to be simulated.
   * \param inputs Vector of TT raw pointers to the input simulation patterns.
   */
  void operator()( List const& list, std::vector<TT const*> const& inputs )
  {
    /* update the allocated memory */
    if ( sims.size() < list.num_gates() )
      sims.resize( std::max<size_t>( sims.size(), list.num_gates() ) ); /* ensure that the constant 0 simulation is correct ( for dynamic truth tables ) */
    if ( ( inputs.size() ) > 0 && ( const0.num_vars() != inputs[0]->num_vars() ) )
      const0 = inputs[0]->construct();
    /* traverse the list in topological order and simulate each node */
    size_t i = 0;
    if constexpr ( std::is_same<List, xag_index_list<true>>::value || std::is_same<List, xag_index_list<false>>::value )
    {
      list.foreach_gate( [&]( element_type const& lit_lhs, element_type const& lit_rhs ) {
        auto const [tt_lhs_ptr, is_lhs_compl] = get_simulation( list, inputs, lit_lhs );
        auto const [tt_rhs_ptr, is_rhs_compl] = get_simulation( list, inputs, lit_rhs );
        sims[i++] = list.is_and( lit_lhs, lit_rhs )
                        ? complement( *tt_lhs_ptr, is_lhs_compl ) & complement( *tt_rhs_ptr, is_rhs_compl )
                        : complement( *tt_lhs_ptr, is_lhs_compl ) ^ complement( *tt_rhs_ptr, is_rhs_compl );
      } );
    }
    else if constexpr ( std::is_same<List, mig_index_list>::value )
    {
      list.foreach_gate( [&]( element_type const& lit0, element_type const& lit1, element_type const& lit2 ) {
        auto const [tt_0_ptr, is_0_compl] = get_simulation( list, inputs, lit0 );
        auto const [tt_1_ptr, is_1_compl] = get_simulation( list, inputs, lit1 );
        auto const [tt_2_ptr, is_2_compl] = get_simulation( list, inputs, lit2 );
        sims[i++] = maj( complement( *tt_0_ptr, is_0_compl ), complement( *tt_1_ptr, is_1_compl ), complement( *tt_2_ptr, is_2_compl ) );
      } );
    }
  }

  /*! \brief Invert a truth table if needed.
   *
   *  Defined for readability.
   *
   * \param tt Truth table.
   * \param is_compl True when the truth table should be complemented.
   * \return The input truth table complemented when `is_compl=true`.
   */
  inline TT complement( TT const& tt, bool is_compl )
  {
    return is_compl ? ~tt : tt;
  }

  /*! \brief Compute the majority function of three truth tables.
   *
   * \param tt0 First truth table.
   * \param tt1 Second truth table.
   * \param tt2 Third truth table.
   * \return Majority-of-3 of `tt0`, `tt1`, and `tt2`.
   */
  inline TT maj( TT const& tt0, TT const& tt1, TT const& tt2 )
  {
    return ( tt0 & tt1 ) | ( tt0 & tt2 ) | ( tt1 & tt2 );
  }

  /*! \brief Return the simulation associated to the literal
   *
   * \param list An XAIG index list, with or without separated header.
   * \param inputs A vector of pointers to the input simulation patterns.
   * \param lit The literal whose simulation we want to extract.
   *
   * The size of `inputs` should be equal to the number of inputs of the list.
   * Keep private to avoid giving external access to memory that could be later corrupted.
   *
   * \return A tuple containing a pointer to the simulation pattern and a flag for complementation.
   */
  [[nodiscard]] std::tuple<TT const*, bool> get_simulation( List const& list, std::vector<TT const*> const& inputs, element_type const& lit )
  {
    if ( list.is_constant( lit ) )
    {
      return { &const0, list.is_complemented( lit ) };
    }
    if ( list.num_pis() != inputs.size() )
      throw std::invalid_argument( "Mismatch between number of PIs and input simulations." );
    if ( list.is_pi( lit ) )
    {
      uint32_t index = list.get_pi_index( lit );
      TT const& sim = *inputs[index];
      return { &sim, list.is_complemented( lit ) };
    }
    uint32_t index = list.get_node_index( lit );
    TT const& sim = sims[index];
    return { &sim, list.is_complemented( lit ) };
  }

  /*! \brief Extract the simulation of a literal
   *
   * \param res Truth table where to store the result.
   * \param list Index list to be simulated.
   * \param inputs Vector of pointers to the input truth tables.
   * \param lit Literal whose simulation we want to extract.
   */
  inline void get_simulation_inline( TT& res, List const& list, std::vector<TT const*> const& inputs, element_type const& lit )
  {
    auto const [tt, is_compl] = get_simulation( list, inputs, lit );
    res = is_compl ? ~( *tt ) : *tt;
  }

private:
  /*! \brief Simulation of the internal nodes ( no inputs and constants ) */
  std::vector<TT> sims;
  /*! \brief Constant 0 simulation */
  TT const0;

}; /* list simulator */

/*! \brief Specialized simulator engine for index lists using a gate library.
 *
 * This engine can be used to efficiently simulate index lists representing
 * small netlists where each gate is taken from a technology library.
 * A simulation pattern is a truth table corresponding to a node’s Boolean
 * behavior under the given input assignments. It pre-allocates the memory
 * necessary to store the simulation patterns of an index list, and extends
 * this memory when the evaluator is required to perform Boolean evaluation of
 * a list that is larger than the current capacity of the simulator. To avoid
 * unnecessary copies, the input simulation patterns must be passed as a vector
 * of raw pointers. The netlist is simulated in topological order by simulating
 * each node using an AIG index list corresponding to a decomposition of its
 * functionality.
 *
 * \tparam Gate Gate-type. Must specify at least the gate's functionality.
 * \tparam TT Truth table type.
 *
 * \verbatim embed:rst

   Example

   .. code-block:: c++

      using list_t = lib_index_list;
      using truth_table_t = kitty::static_truth_table<4u>;
      list_simulator<list_t, truth_table_t> sim;
      sim( list1, inputs1 );
      sim( list2, inputs2 );
   \endverbatim
 */
template<typename Gate, typename TT>
class list_simulator<lib_index_list<Gate>, TT>
{
public:
  using outer_list_t = lib_index_list<Gate>;
  using inner_list_t = large_xag_index_list;
  using element_type = typename outer_list_t::element_type;

  /*! \brief Construction requires the specification of the gate-library.
   *
   * \param library Vector of gates where at least the functionality is specified.
   *
   * The size of `inputs` should be equal to the number of inputs of the list.
   * Keep private to avoid giving external access to memory that could be later corrupted.
   *
   * \return A tuple containing a pointer to the simulation pattern and a flag for complementation.
   */
  list_simulator( std::vector<Gate> const& library )
      : library( library ),
        inner_simulator()
  {
    /* The value 20 allows us to store practical lists */
    sims.resize( 20u );
  }

  /*! \brief Simulate the list in topological order.
   *
   * This method updates the internal state of the simulator by
   * storing in `sims` the simulation patterns of the nodes in the list.
   *
   * \param outer_list Index list of gates to be simulated.
   * \param inputs Vector of TT raw pointers to the input patterns.
   */
  void operator()( outer_list_t const& outer_list,
                   std::vector<TT const*> const& inputs )
  {
    /* update the allocated memory */
    if ( sims.size() < outer_list.num_gates() )
      sims.resize( std::max<size_t>( sims.size(), outer_list.num_gates() ) );

    /* traverse the list in topological order and simulate each node */
    size_t i = 0;
    using iterate_type = typename std::vector<element_type>::iterator;
    std::vector<TT const*> sims_ptrs;
    outer_list.foreach_gate( [&]( auto const& start, auto const& end, element_type const& id ) {
      sims_ptrs.clear();

      for ( auto it = start; it != end; it++ )
      {
        element_type lit = *it;
        if ( outer_list.is_pi( lit ) )
        {
          auto const index = outer_list.get_pi_index( lit );
          sims_ptrs.push_back( inputs[index] );
        }
        else
        {
          auto const index = outer_list.get_node_index( lit );
          sims_ptrs.push_back( &sims[index] );
        }
      }
      inner_list_t const& inner_list = library.get_list( id );
      inner_simulator( inner_list, sims_ptrs );
      /* each gate has a single outout, multiple-output are represented as two gates. */
      auto const lit = inner_list.po_at( 0 );
      inner_simulator.get_simulation_inline( sims[i++],
                                             inner_list,
                                             sims_ptrs,
                                             lit );
    } );
  }

  /*! \brief Return the simulation associated to the literal
   *
   * \param list An XAIG index list, with or without separated header.
   * \param inputs A vector of pointers to the input simulation patterns.
   * \param lit The literal whose simulation we want to extract.
   * \return The simulation pattern.
   */
  [[nodiscard]] TT const& get_simulation( outer_list_t const& list,
                                          std::vector<TT const*> const& inputs,
                                          element_type const& lit )
  {
    if ( list.num_pis() != inputs.size() )
      throw std::invalid_argument( "Mismatch between number of PIs and input simulations." );
    if ( list.is_pi( lit ) )
    {
      uint32_t index = list.get_pi_index( lit );
      return *inputs[index];
    }
    uint32_t index = list.get_node_index( lit );
    return sims[index];
  }

  /*! \brief Extract the simulation of a literal
   *
   * Inline specifier used to ensure manual inlining.
   *
   * \param res Truth table where to store the result.
   * \param list Index list to be simulated.
   * \param inputs Vector of pointers to the input truth tables.
   * \param lit Literal whose simulation we want to extract.
   */
  inline void get_simulation_inline( TT& res,
                                     outer_list_t const& list,
                                     std::vector<TT const*> const& inputs,
                                     element_type const& lit )
  {
    res = get_simulation( list, inputs, lit );
  }

private:
  /*! \brief Simulation patterns of the list's nodes */
  std::vector<TT> sims;
  /*! \brief Augmented library */
  augmented_library<Gate> library;
  /*! \brief Simulator engine for the individual nodes */
  list_simulator<inner_list_t, TT> inner_simulator;
};

} /* namespace mockturtle */