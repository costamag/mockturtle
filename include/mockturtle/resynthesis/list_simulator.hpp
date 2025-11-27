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

  This engine can be use in contexts in which many index lists
  need to be evaluated. It pre-allocates a memory necessary to
  store the simulation patterns of an index list, and this memory
  is extended when the evaluator is required to perform Boolean
  evaluation of a list that is larger than the current capacity
  of the simulator.

  \author Andrea Costamagna
*/

#pragma once

#include "../utils/index_list.hpp"

namespace mockturtle
{

template<class Tt_t, uint32_t Capacity = 20>
class aig_list_simulator
{
public:
  using index_list_t = large_xag_index_list;

  aig_list_simulator()
  {
    static_assert( std::is_same_v<typename index_list_t::element_type, uint32_t>, "Incompatible element type with the synthesis implementation" );
    sims.reserve( Capacity );
  }

  /*! \brief Prepare the input simulation patterns and call the simulation method.
   *
   * `tts[*begin]` must be of type `Tt_t`.
   *
   * \param res Truth table where to store the result.
   * \param list Index list to be simulated.
   * \param begin Begin iterator to the list's inputs.
   * \param end End iterator to the list's outputs.
   * \param tts A data structure (e.g. std::vector<Tt_t>) 
   *            storing the simulation patterns and accessible
   *            with the iterator.
   */
  template<class Iterator_t, class Storage_t>
  void operator()( Tt_t & res, index_list_t const& list, Iterator_t begin, Iterator_t end, Storage_t const& tts )
  {
    /* update the capacity if necessary */
    sims.reserve( list.num_gates() + list.num_pis() + 2 );
    /* reserve one dummy simulation pattern for the constant 0 */
    sims.resize( 1 );
    sims[0] = res.construct();
    while( begin != end )
    {
      sims.emplace_back( tts[*begin] );
      ++begin;
    }
    simulate( res, list );
  }

  void operator()( Tt_t & res, index_list_t const& list, std::vector<Tt_t*> const& sim_ptrs )
  {
    /* update the capacity if necessary */
    sims.reserve( list.num_gates() + list.num_pis() + 2 );
    /* reserve one dummy simulation pattern for the constant 0 */
    sims.resize( 1 );
    sims[0] = res.construct();
    for( auto ptr : sim_ptrs )
    {
      sims.emplace_back( *ptr );
    }
    simulate( res, list );
  }

  void operator()( Tt_t & res, index_list_t const& list, std::vector<Tt_t> const& tts )
  {
    /* update the capacity if necessary */
    sims.reserve( list.num_gates() + list.num_pis() + 2 );
    /* reserve one dummy simulation pattern for the constant 0 */
    sims.resize( 1 );
    sims[0] = res.construct();
    for( auto tt : tts )
    {
      sims.emplace_back( tt );
    }
    simulate( res, list );
  }

  /*! \brief Simulate the index list.
   *
   * \param res Truth table where to store the result.
   * \param list Index list to be simulated.
   */
  void simulate( Tt_t & res, index_list_t const& list )
  {
    list.foreach_gate( [&]( uint32_t lit0, uint32_t lit1 ){
      uint32_t const index0 = lit0 >> 1u;
      uint32_t const index1 = lit1 >> 1u;
      bool const is_compl0 = ( lit0 & 0x1 ) == 0x1;
      bool const is_compl1 = ( lit1 & 0x1 ) == 0x1;
      Tt_t const sim0 = ( lit0 & 0x1 ) == 0x1 ? ~sims[index0] : sims[index0];
      Tt_t const sim1 = ( lit1 & 0x1 ) == 0x1 ? ~sims[index1] : sims[index1];
      sims.emplace_back( sim0 & sim1 );
    } );
    /* the last entry of the list is the output */
    uint32_t const lit = list.po_at( 0 );
    uint32_t const index = lit >> 1u;
    bool const is_compl = ( lit & 0x1 ) == 0x1;
    res = is_compl ? ~sims[index] : sims[index];
  }

private:
  std::vector<Tt_t> sims;
}; /* list simulator */

} /* namespace mockturtle */