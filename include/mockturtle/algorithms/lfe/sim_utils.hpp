/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2021  EPFL
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
  \file utils.hpp
  \brief Helper functions for the LFE algorithms

  \author Andrea Costamagna
*/

#pragma once

#include <vector>

namespace mockturtle
{
  template<typename T>
  void swap_T( T& a, T& b)
  {
    T t = a;
    a = b;
    b = t;
  }
  template<typename TA, typename TB>
  int partition ( std::vector<TA>& support, std::vector<TB>& attribute , uint32_t low, uint32_t high )
  {
  TB pivot = attribute[high];    // pivot
  int first_small = ( low - 1 );  // Index of smaller element

  for ( uint32_t j = low; j < high; j++ )
    {
      if ( attribute[j] >= pivot )
      {
          first_small++;    // increment index of larger element
          swap_T(attribute[first_small], attribute[j]);
          swap_T(support[first_small], support[j]);
      }
    }
    swap_T(attribute[first_small + 1], attribute[high]);
    swap_T(support[first_small + 1], support[high]);
    return (first_small + 1);
  }

  template<typename TA, typename TB>
  void rquicksort_by_attribute( std::vector<TA>& support, std::vector<TB>& attribute,  int low, int high )
  {
    if( low >= high )
      return;
    else
    {
      int pi = partition( support, attribute, low, high);
      rquicksort_by_attribute( support, attribute, low, pi - 1);
      rquicksort_by_attribute( support, attribute, pi+1, high);
    }
  }

  template<typename TA, typename TB>
  void quicksort_by_attribute( std::vector<TA>& support, std::vector<TB>& attribute )
  {
    assert( support.size() == attribute.size() );
    int low = 0;
    int high = support.size() - 1;
    rquicksort_by_attribute( support, attribute, low, high );
  }

  #pragma region cum sum
  double Pk_f( uint32_t const& k, uint32_t const& N0, uint32_t const& N1, uint const& n )
  {
    
    uint32_t Nh = std::max(N0,N1);
    uint32_t Nl = std::min(N0,N1);
    uint32_t n_inf = 10;
    if( n>n_inf || Nl==0 || Nh == 0)
      return (k==0 ? (double)1:(double)0);
    double Pk = 1;
    if(k>Nl)
    {
      return 0;
    }
    if( pow(2,n-1)+k < Nh+Nl) 
    {
      return 0;
    }
    if(Nh==pow(2,n-1))
    {
      if( k==Nl )
      {
        return 1;
      } 
    }

    for(uint32_t j=0;j<Nl-k;++j)
    {
      Pk*=(1-(double)Nh/(pow(2,n-1)-j));
    }
    
    if(k!=0)
    {
      for( uint32_t j = 0; j<k ; ++j )
      {
        double Ak = (double)(Nl-j)/(j+1);
        double Bk = (double)(Nh-j)/(pow(2,n-1)-Nl+j+1);
        Pk *= Ak*Bk;
      }
    }

    return Pk;
  }

  std::pair<double,double> M1M2k(uint32_t const& N0, uint32_t const& N1, uint32_t const& n )
  {
    uint32_t Nh = std::max(N0,N1);
    uint32_t Nl = std::min(N0,N1);
    uint32_t n_inf = 32;
    if( n>n_inf )
      return std::make_pair(0,0);
    uint32_t kmin = std::max(1, (int)((Nh+Nl)-pow(2,n-1)));
    double Pk = Pk_f( kmin, N0, N1, n );
    double M1 = kmin*Pk;
    double M2 = kmin*kmin*Pk;
    for( uint32_t k=kmin+1; k< Nl+1 ;++k)
    {
      double Ak = k*Pk_f( k, N0, N1, n );
      M1 += Ak;
      M2 += Ak*k;
    }
    return std::make_pair(M1,sqrt(M2-M1*M1));
  }

  double CumSum(uint32_t const& kmax, uint32_t const& N0, uint32_t const& N1, uint32_t const& n )
  {
    double CS =0;
    for(uint32_t k=0; k<kmax+1; ++k)
    {
      double dP=Pk_f( k, N0, N1, n );
      CS+=dP;
    }
    return CS;
  }
#pragma endregion cum sum
}