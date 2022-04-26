/* kitty: C++ truth table library
 * Copyright (C) 2017-2021  EPFL
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
  \file statistics.hpp
  \brief Computes statistical properties of truth tables

  \author Andrea Costamagna
*/

#pragma once

#include "partial_truth_table.hpp"
#include "dynamic_truth_table.hpp"
#include "static_truth_table.hpp"
#include "bit_operations.hpp"
#include "constructors.hpp"
#include <algorithm>
#include <iostream>
#include <vector>

#include <cmath>

namespace kitty
{

#pragma region probability
/*! \brief Compute the probability for the binary value of a single truth table

  \param tt Truth table
*/
template<typename TT>
std::vector<double> probability( TT& tt )
{
  double p = (double)count_ones(tt)/tt.num_bits();
  std::vector<double> probs = { 1-p, p };
  return probs;
}

/*! \brief Compute the probability for the binary value of a vector of truth tables

  \param tts Vector of Truth tables
*/
template<typename TT>
std::vector<double> probability( std::vector<TT>& tts )
{
  auto ntts = tts.size();
  size_t nvars;
  nvars = tts[0].num_bits();
  std::vector<double> probs;
  TT const0( nvars );
  TT ttc( nvars );

  for( uint32_t p = 0; p < pow( 2, ntts ); ++p )
  {
    ttc = ~const0;
    for( uint32_t i = 0; i < ntts; ++i )
    {
      bool bit = ( p >> i ) & 1u;

      if( bit == 1 )
        ttc &= tts[i];
      else  
        ttc &= ~tts[i];
    }
    probs.push_back( (double)count_ones(ttc)/ttc.num_bits() );
  }

  return probs;
}
/*
template <typename TT>
std::vector<double> probability( std::vector<TT>& tt )
{
  if constexpr( std::is_same_v< std::vector<kitty::partial_truth_table>, TT> )
    return probability( tt );
  else
    return probability_v( tt );
}
*/
/*! \brief Compute the probability for the binary value of a single truth table

  \param tt Truth table
  \param index Bit index
*/
template<typename TT>
std::vector<double> probability( std::vector<TT>& tt1, std::vector<TT>& tt2 )
{
  auto n1 = tt1[0].num_vars();
  auto n2 = tt2[0].num_vars();
  assert( ( n1 == n2 ) );

  std::vector<kitty::dynamic_truth_table> tts;

  for( auto tt : tt2 )
    tt1.push_back( tt );

  return probability( tt1 );
}


/*! \brief Compute the probability for the binary value of a single truth table

  \param tt Truth table
  \param index Bit index
*/
template<typename TT>
std::vector<double> probability( TT& tt1, TT& tt2 )
{
  std::vector<TT> dbv = { tt1, tt2 };

  return probability( dbv );
}

/*! \brief Compute the probability for the binary value of a single truth table

  \param tt Truth table
  \param index Bit index
*/
template<typename TT>
std::vector<double> probability( std::vector<TT>& tts, TT& tt )
{
  auto n1 = tts[0].num_bits();
  auto n2 = tt.num_bits();
  assert( ( n1 == n2 ) );

  tts.push_back( tt );

  return probability( tts );
}

#pragma endregion probability

#pragma region entropy
/*! \brief Compute the entropy of a single truth table

  \param tt Truth table
*/
template<typename X>
double entropy( X& x )
{
  auto probs = probability( x );
  double H = 0;
  for( auto i = 0; i < probs.size(); ++i )
    H += ( ( probs[i] > 0 ) ? -probs[i]*log2( probs[i] ) : 0 );
  return H;  
}

/*! \brief Compute the entropy of a single truth table

  \param tt Truth table
*/
template<typename X, typename Y>
double entropy( X& x, Y& y )
{
  auto probs = probability( x, y );
  double H = 0;
  for( auto i = 0; i < probs.size(); ++i )
    H += ( ( probs[i] > 0 ) ? -probs[i]*log2( probs[i] ) : 0 );
  return H;  
}


#pragma endregion entropy

#pragma region mutual information
/*! \brief Compute the mutual information of two 

  \param tt Truth table
*/
template<typename X, typename Y>
double mutual_information( X& x, Y& y )
{
  return entropy( x ) + entropy( y ) - entropy( x, y );  
}
#pragma endregion mutual information

} // namespace kitty