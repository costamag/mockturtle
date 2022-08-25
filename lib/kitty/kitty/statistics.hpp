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
  \brief Computes statistical properties of Boolean functions
  \author Andrea Costamagna
*/

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>
#include <cstdlib>

#include "bit_operations.hpp"
#include "constructors.hpp"
#include "traits.hpp"

namespace kitty
{
/*! \brief Computes the probabilities of the random variable X

  The random variable \f$X\in\mathbb{B}^{n}\f$ is represented using \f$n\f$ Truth tables. 
  The collection of the \f$i\f$-th bit of all truth tables identifies a pattern \f$\pi_i\in\mathbb{B}^{n}\f$, i.e., a sampling of \f$X\f$.
  This function computes the probability \f$\mathbb{P}(\pi)\f$ of each pattern \f$\pi\in\mathbb{B}^{n}\f$.

  \param X \f$n\f$-dimensional vector of Truth tables or of pointers to Truth tables
*/
template<typename TT>
std::vector<double> probabilities( std::vector<TT> const& X )
{

  size_t ntts = X.size();
  size_t nconstr = std::is_same<TT, kitty::dynamic_truth_table>::value ? (int)log2( X[0].num_bits() ) : X[0].num_bits();

  std::vector<double> probs;
  TT const0( nconstr );
  TT mask( nconstr );

  for ( uint32_t p = 0; p < pow( 2, ntts ); ++p )
  {
    mask = ~const0;
    for ( uint32_t i = 0; i < ntts; ++i )
    {
      bool bit = ( p >> i ) & 1u;

      if ( bit == 1 )
        mask &= X[i];
      else
        mask &= ~X[i];
    }
    probs.push_back( (double)count_ones( mask ) / mask.num_bits() );
  }

  return probs;
}

template<typename TT>
std::vector<double> probabilities( std::vector<TT*> const& X )
{

  size_t ntts = X.size();
  size_t nconstr = std::is_same<TT, kitty::dynamic_truth_table>::value ? (int)log2( (*X[0]).num_bits() ) : (*X[0]).num_bits();

  std::vector<double> probs;
  TT const0( nconstr );
  TT mask( nconstr );

  for ( uint32_t p = 0; p < pow( 2, ntts ); ++p )
  {
    mask = ~const0;
    for ( uint32_t i = 0; i < ntts; ++i )
    {
      bool bit = ( p >> i ) & 1u;

      if ( bit == 1 )
        mask &= *X[i];
      else
        mask &= ~*X[i];
    }
    probs.push_back( (double)count_ones( mask ) / mask.num_bits() );
  }

  return probs;
}


/*! \brief Computes the joint probabilities of random variables X and Y

  The random variables \f$X\in\mathbb{B}^{n_1}\f$ and \f$Y\in\mathbb{B}^{n_2}\f$ are represented using \f$n_1\f$ and \f$n_2\f$ Truth tables, respectively. 
  The collection of the \f$i\f$-th bit of all truth tables identifies a pattern \f$\pi_i=\pi^{(1)}_i\pi^{(2)}_i\in\mathbb{B}^{n}\f$, i.e., a sampling of \f$(X,Y)\f$.
  This function computes the probability \f$\mathbb{P}(\pi)\f$ of each pattern \f$\pi\in\mathbb{B}^{n}\f$.
  
  \param X \f$n_1\f$-dimensional vector of Truth tables or of pointers to Truth tables
  \param Y \f$n_2\f$-dimensional vector of Truth tables or of pointers to Truth tables
*/
template<typename TT>
std::vector<double> probabilities( std::vector<TT> const& X, std::vector<TT> const& Y )
{

  std::vector<TT> U;

  U = X;
  U.insert( U.end(), Y.begin(), Y.end() );

  return probabilities( U );
}

/*! \brief Computes the entropy of the random variable X.

  The random variable \f$X\in\mathbb{B}^{n}\f$ is represented using \f$n\f$ Truth tables. 
  The collection of the \f$i\f$-th bit of all truth tables identifies a pattern \f$\pi_i\in\mathbb{B}^{n}\f$, i.e., a sampling of \f$X\f$.
  The entropy quantifies the uncertainty on the value of \f$X\f$:

  \f$H(X)=-\sum_{\pi\in\mathbb{B}^n}\mathbb{P}(\pi)\cdot\log_2\mathbb{P}(\pi)\f$

  \param X \f$n\f$-dimensional vector of Truth tables or of pointers to Truth tables
*/
template<typename TT>
double entropy( std::vector<TT> const& X )
{
  std::vector<double> probs = probabilities( X );
  double H = 0;
  for ( size_t i = 0; i < probs.size(); ++i )
    H += ( ( probs[i] > (double)0.0 ) ? -probs[i] * log2( probs[i] ) : (double)0.0 );
  return H;
}

/*!  \brief Computes the entropy of the patterns in X and Y.

  The random variables \f$X\in\mathbb{B}^{n_1}\f$ and \f$Y\in\mathbb{B}^{n_2}\f$ are represented using \f$n_1\f$ and \f$n_2\f$ Truth tables, respectively. 
  The collection of the \f$i\f$-th bit of all truth tables identifies a pattern \f$\pi_i=\pi^{(1)}_i\pi^{(2)}_i\in\mathbb{B}^{n}\f$, i.e., a sampling of \f$(X,Y)\f$.
  The entropy quantifies the uncertainty on the value of \f$(X,Y)\f$:
  
  \f$H(X,Y)=-\sum_{\pi\in\mathbb{B}^n}\mathbb{P}(\pi)\cdot\log_2\mathbb{P}(\pi)\f$

  \param X \f$n_1\f$-dimensional vector of Truth tables or of pointers to Truth tables
  \param Y \f$n_2\f$-dimensional vector of Truth tables or of pointers to Truth tables
*/
template<typename TT>
double entropy( std::vector<TT> const& X, std::vector<TT> const& Y )
{
  std::vector<double> probs = probabilities( X, Y );
  double H = 0;
  for ( size_t i = 0; i < probs.size(); ++i )
    H += ( ( probs[i] > 0 ) ? -probs[i] * log2( probs[i] ) : 0 );
  return H;
}

/*! \brief Computes the mutual information of random variables X and Y

  The random variables \f$X\in\mathbb{B}^{n_1}\f$ and \f$Y\in\mathbb{B}^{n_2}\f$ are represented using \f$n_1\f$ and \f$n_2\f$ Truth tables, respectively. 
  The mutual information quantifies the reduction in uncertainty on \f$Y\f$, given that \f$X\f$ is known:
  
  \f$I(X;Y)=H(Y)-H(Y|X)=H(X)+H(Y)-H(X,Y)\f$

  \param X \f$n_1\f$-dimensional Vector of Truth tables or of pointers to Truth tables
  \param Y \f$n_2\f$-dimensional Vector of Truth tables or of pointers to Truth tables
*/
template<typename TT>
double mutual_information( std::vector<TT> const& X, std::vector<TT> const& Y )
{
  double mi = entropy( X ) + entropy( Y ) - entropy( X, Y );
  return ( mi > 1e-14 ) ? mi : 0;
}

/*! \brief Computes the mutual information of random variables X and y

  The random variables \f$X\in\mathbb{B}^{n_1}\f$ and \f$y\in\mathbb{B}\f$ are represented using \f$n_1\f$ and \f$1\f$ Truth tables, respectively. 
  The mutual information quantifies the reduction in uncertainty on \f$y\f$, given that \f$X\f$ is known:
  
  \f$I(X;y)=H(y)-H(y|X)=H(X)+H(y)-H(X,y)\f$

  \param X \f$n_1\f$-dimensional Vector of Truth tables or of pointers to Truth tables
  \param y Truth table or of pointer to Truth tables
*/
template<typename TT>
double mutual_information( std::vector<TT> const& X, TT const& y )
{
  std::vector<TT> Y = { y };

  return mutual_information( X, Y );
}

/*! \brief Computes the normalized mutual information of random variables X and Y

  The random variables \f$X\in\mathbb{B}^{n_1}\f$ and \f$Y\in\mathbb{B}^{n_2}\f$ are represented using \f$n_1\f$ and \f$n_2\f$ Truth tables, respectively. 
  The normalized mutual information quantifies the reduction in uncertainty on \f$Y\f$, given that \f$X\f$ is known:
  
  \f$NI(X;Y)=\frac{H(X)+H(Y)}{H(X,Y)}\f$

  \param X \f$n_1\f$-dimensional Vector of Truth tables or of pointers to Truth tables
  \param Y \f$n_2\f$-dimensional Vector of Truth tables or of pointers to Truth tables
*/
template<typename TT>
double normalized_mutual_information( std::vector<TT> const& X, std::vector<TT> const& Y )
{
  double mi = ( entropy( X ) + entropy( Y ) ) / entropy( X, Y );
  return ( mi > 1e-14 ) ? mi : 0;
}

/*! \brief Computes the normalized mutual information of random variables X and y

  The random variables \f$X\in\mathbb{B}^{n_1}\f$ and \f$y\in\mathbb{B}\f$ are represented using \f$n_1\f$ and \f$1\f$ Truth tables, respectively. 
  The normalized mutual information quantifies the reduction in uncertainty on \f$y\f$, given that \f$X\f$ is known:
  
  \f$NI(X;y)=\frac{H(X)+H(y)}{H(X,y)}\f$

  \param X \f$n_1\f$-dimensional Vector of Truth tables or of pointers to Truth tables
  \param y Truth table or pointer to Truth table
*/
template<typename TT>
double normalized_mutual_information( std::vector<TT> const& X, TT const& y )
{
  std::vector<TT> Y = { y };
  return normalized_mutual_information( X, Y );
}


/*! \brief Computes the correlation of random variables X and y

  The random variables \f$X\in\mathbb{B}^{n_1}\f$ and \f$y\in\mathbb{B}\f$ are represented using \f$n_1\f$ and \f$1\f$ Truth tables, respectively. 
  The correlation quantifies the agrement/disagrement relationship of \f$y\f$ and \f$X\f$:
  
  \f$C(X;y)=agreements-disagreements\f$

  \param X Truth table or of pointer to Truth tables
  \param y Truth table or of pointer to Truth tables
*/
//template<typename TT>
int correlation( kitty::partial_truth_table* const& X, kitty::partial_truth_table* const& Y )
{
  return abs((int)((*X).num_bits()-2*( count_ones( (*X)^(*Y) ) )));
}

//template<typename TT>
int correlation( kitty::partial_truth_table const& X, kitty::partial_truth_table const& Y )
{
  return abs((int)(X.num_bits()-2*( count_ones( X^Y ) )));
}

int norm_correlation( kitty::partial_truth_table* const& X, kitty::partial_truth_table* const& Y )
{
  return std::abs(kitty::count_ones( (*X)&(*Y) )-(double)kitty::count_ones( *X )*kitty::count_ones( *Y )/((*Y).num_bits())); 

}

//template<typename TT>
int norm_correlation( kitty::partial_truth_table const& X, kitty::partial_truth_table const& Y )
{
  return std::abs(kitty::count_ones( X&Y )-(double)kitty::count_ones( X )*kitty::count_ones( Y )/(Y.num_bits()));
}



} // namespace kitty