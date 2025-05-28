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
  \file index_list.hpp
  \brief List of indices to represent small networks.

  \author Andrea Costamagna
  \author Heinz Riener
  \author Siang-Yun (Sonia) Lee
*/

#pragma once

#include "../../traits.hpp"

#include "lists/mig_index_list.hpp"
#include "lists/muxig_index_list.hpp"
#include "lists/xag_index_list.hpp"

#include <vector>

namespace mockturtle
{

/*! \brief Generates a network from an index_list
 *
 * **Required network functions:**
 * - `create_pi`
 * - `create_po`
 *
 * \param ntk A logic network
 * \param indices An index list
 */
template<typename Ntk, typename IndexList>
void decode( Ntk& ntk, IndexList const& indices )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_create_pi_v<Ntk>, "Ntk does not implement the create_pi method" );
  static_assert( has_create_po_v<Ntk>, "Ntk does not implement the create_po method" );

  using signal = typename Ntk::signal;

  std::vector<signal> signals( indices.num_pis() );
  std::generate( std::begin( signals ), std::end( signals ),
                 [&]() { return ntk.create_pi(); } );

  insert( ntk, std::begin( signals ), std::end( signals ), indices,
          [&]( signal const& s ) { ntk.create_po( s ); } );
}

template<class T>
struct is_index_list : std::false_type
{
};

template<>
struct is_index_list<xag_index_list<true>> : std::true_type
{
};

template<>
struct is_index_list<xag_index_list<false>> : std::true_type
{
};

template<>
struct is_index_list<mig_index_list> : std::true_type
{
};

template<class T>
inline constexpr bool is_index_list_v = is_index_list<T>::value;

} // namespace mockturtle
