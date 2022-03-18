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
  \file statistical_bit_operations.hpp
  \brief Implements statistical computations on truth tables

  \author Andrea Costamagna
*/

#pragma once

#include "bit_operations.hpp"
#include <algorithm>
#include <iostream>
#include <vector>
#include <bitset>

#include <cmath>

namespace kitty
{
using DBS = boost::dynamic_bitset<>;
#pragma region probability


/*! \brief Compute the probability for the binary value of a vector of truth tables

  \param tts Vector of Truth tables
*/
std::vector<double> probability( std::vector<DBS> const& X )
{
  auto ntts = X.size();
  size_t nvars;
  nvars = X[0].size();
  std::vector<double> probs;
  DBS const0( nvars, 0u );
  DBS ttc( nvars, 0u );

  for( uint32_t p = 0; p < pow( 2, ntts ); ++p )
  {
    ttc = ~const0;
    for( uint32_t i = 0; i < ntts; ++i )
    {
      bool bit = ( p >> i ) & 1u;

      if( bit == 1 )
        ttc &= X[i];
      else  
        ttc &= ~X[i];
    }
    probs.push_back( (double)ttc.count()/ttc.size() );
  }

  return probs;
}
#pragma endregion probability

#pragma region entropy
/*! \brief Compute the entropy of a single truth table

  \param x first signal 
  \param y second signal 
*/
double entropy( DBS const& x, DBS const& y )
{
  std::vector<DBS> X = { x, y };
  auto probs = probability( X );
  double H = 0;
  for( auto i = 0; i < probs.size(); ++i )
    H += ( ( probs[i] > 0 ) ? -probs[i]*log2( probs[i] ) : 0 );
  return H;  
}

/*! \brief Compute the entropy of a single truth table

  \param x signal 
*/
double entropy( DBS const& x )
{
  std::vector<DBS> X = {x};
  auto probs = probability( X );
  double H = 0;
  for( auto i = 0; i < probs.size(); ++i )
    H += ( ( probs[i] > 0 ) ? -probs[i]*log2( probs[i] ) : 0 );
  return H;  
}

/*! \brief Compute the entropy of a single truth table

  \param x first signal 
  \param y second signal 
*/
double entropy( std::vector<DBS> const& X )
{
  auto probs = probability( X );
  double H = 0;
  for( auto i = 0; i < probs.size(); ++i )
    H += ( ( probs[i] > 0 ) ? -probs[i]*log2( probs[i] ) : 0 );
  return H;  
}

/*! \brief Compute the entropy of a single truth table

  \param x first signal 
  \param y second signal 
*/
double entropy( std::vector<DBS> const& X, DBS const& y )
{
  std::vector<DBS> Xn = X;
  Xn.push_back( y );
  auto probs = probability( Xn );
  double H = 0;
  for( auto i = 0; i < probs.size(); ++i )
    H += ( ( probs[i] > 0 ) ? -probs[i]*log2( probs[i] ) : 0 );
  return H;  
}
#pragma endregion entropy

#pragma region mutual information
/*! \brief Compute the mutual information of two 

  \param x first signal
  \param y second signal
*/
double mutual_information( DBS const& x, DBS const& y )
{
  return entropy( x ) + entropy( y ) - entropy( x, y );  
}

/*! \brief Compute the mutual information of two 

  \param x first signal
  \param y second signal
*/
double mutual_information( std::vector<DBS> const& X, DBS const& y )
{
  return entropy( X ) + entropy( y ) - entropy( X, y );  
}
#pragma endregion mutual information

} // namespace kitty