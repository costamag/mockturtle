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
  \file timing.hpp
  \brief Compute timing information of a network.

  This engine can be used for power analysis of mapped network.
  For each node, the following information is stored:
  - The sensing time : the first time at which a transition can happen
  - The arrival time : the first time at which the output is stable
  - A vector of simulation patterns identifying quantized timesteps in this interval
  The class has the following template parameters:
  \param N : number of timesteps used for simulating each signal
  \param I : number of inputs to be used for simulating the static truth table ( 2^I input pairs )

  \author Andrea Costamagna
*/

#pragma once

#include <array>

namespace mockturtle
{


struct bin_t
{
  int x;
  int y;
};

template<class T>
struct vector_t
{
  vector_t( T const& a, T const& b )
   : data( { a, b } )
  {
  }

  vector_t( std::array<T, 2> const& v )
   : data( v )
  {
  }

  vector_t( vector_t const& other )
   : vector_t( other.data )
  {
  }

  double squared_norm2()
  {
    return pow( data[0], 2 ) + pow( data[1], 2 );
  }

  double norm2()
  {
    return sqrt( squared_norm2() );
  }

  T & operator[]( int i )
  {
    assert( ( i == 0 ) || ( i == 1 ) );
    return data[i];
  }

  vector_t<T> operator-( vector_t<T> const & rhs ) const
  {
    return vector_t<T>( data[0] - rhs.data[0], data[1] - rhs.data[1] );
  }

  vector_t<T> operator+( vector_t<T> const & rhs ) const
  {
    return vector_t<T>( data[0] + rhs.data[0], data[1] + rhs.data[1] );
  }

  vector_t<T> times( T const & factor ) const
  {
    return vector_t<T>( data[0] * factor, data[1] * factor );
  }

  T dot ( vector_t<T> & v ) const
  {
    return ( data[0] * v[0] ) + ( data[1] * v[1] );
  }


  std::array<T, 2> data;
};

} // namespace mockturtle

