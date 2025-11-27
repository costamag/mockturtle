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
  \file mapped_list_simulator.hpp
  \brief Simulator engine for mapped index lists.

  This engine can be use in contexts in which many index lists
  need to be evaluated. It pre-allocates a memory necessary to
  store the simulation patterns of an index list, and this memory
  is extended when the evaluator is required to perform Boolean
  evaluation of a list that is larger than the current capacity
  of the simulator.

  \author Andrea Costamagna
*/

#pragma once

#include "list_simulator.hpp"
#include "database/mapped_list.hpp"
#include "../utils/index_list.hpp"

namespace mockturtle
{

template<class Bindings_t, class Tt_t, uint32_t Capacity = 20>
class mapped_list_simulator
{
public:
  explicit mapped_list_simulator( Bindings_t const& bindings ) noexcept
    : bindings( bindings ),
    simulator()
  {
    sims.reserve( Capacity );
  }

  void operator()( Tt_t& res, mapped_list const& mlist, std::vector<Tt_t*> sim_ptrs )
  {
    /* update the capacity if necessary */
    sims.reserve( mlist.num_gates() + mlist.num_pis() + 2 );
    /* reserve one dummy simulation pattern for the constant 0 */
    sims.resize( 2 );
    sims[0] = res.construct();
    sims[1] = ~res.construct();
    for ( auto ptr : sim_ptrs )
    {
      sims.emplace_back( *ptr );
    }
    for( int i = sim_ptrs.size(); i < mlist.num_pis(); ++i )
    {
      sims.emplace_back( sims[0] );
    }

    mlist.foreach_gate( [&]( auto children, auto id ) {
      sim_ptrs.clear();
      for ( auto f : children )
      {
        sim_ptrs.emplace_back( &sims[f] );
      }
      auto const list = bindings.get_list( id );
      uint32_t const index = sims.size();
      sims.emplace_back();
      simulator( sims[index], list, sim_ptrs );
    } );

    res = sims[mlist.po_at(0)];
  }

  double operator()( Tt_t& res, mapped_list const& mlist, std::vector<Tt_t> const& input_sims )
  {
    double switching = 0;
    /* update the capacity if necessary */
    sims.reserve( mlist.num_gates() + mlist.num_pis() + 2 );
    /* reserve one dummy simulation pattern for the constant 0 */
    sims.resize( 2 );
    sims[0] = res.construct();
    sims[1] = ~res.construct();
    for ( int i = 0; i < input_sims.size(); ++i )
    {
      sims.emplace_back( input_sims[i] );
    }
    for( int i = input_sims.size(); i < mlist.num_pis(); ++i )
    {
      sims.emplace_back( sims[0] );
    }

    double num_bits = res.num_bits();
    mlist.foreach_gate( [&]( auto children, auto id ) {
      auto const list = bindings.get_list( id );
      uint32_t const index = sims.size();
      sims.emplace_back();
      simulator( sims[index], list, children.begin(), children.end(), sims );
      double num_zeros = kitty::count_zeros( sims[index] );
      switching += 2.0*( num_zeros*( num_bits - num_zeros ) )/( num_bits * num_bits );
    } );
    res = sims[mlist.po_at(0)];

    return switching;
  }


private:
  std::vector<Tt_t> sims;
  Bindings_t bindings;
  aig_list_simulator<Tt_t, Capacity> simulator;
}; /* list simulator */

} /* namespace mockturtle */