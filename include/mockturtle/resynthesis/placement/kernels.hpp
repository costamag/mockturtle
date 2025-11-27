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

#include "placement_utils.hpp"
#include <array>
#include <limits>
#include <math.h>


namespace mockturtle
{
template<int N>
class gaussian_kernel_t
{
public:
  gaussian_kernel_t():
    sigma( 1 ), delta( 1 ), border( 0 )
  {
  }

  void set_sigma( double sigma_new )
  {
    sigma = sigma_new;
  }

  void set_delta( double delta_new )
  {
    delta = delta_new;
  }

  void set_border( double border_new )
  {
    border = border_new;
  }


  vector_t<double> slot_to_coord( int i, int j )
  {
    double const dx = delta * ( i - ( N / 2 ) );
    double const dy = delta * ( j - ( N / 2 ) );
    vector_t<double> const coord { dx, dy };
    return coord;
  }

  vector_t<double> get_delta_min( vector_t<double> & q )
  {
    double dx = -delta * ( N / 2 );
    double dy = -delta * ( N / 2 );
    dx = dx+q[0] < 0 ? -q[0] + delta : dx;
    dy = dy+q[1] < 0 ? -q[1] + delta : dy;
    vector_t<double> const coord { dx, dy };
    return coord;
  }

  vector_t<double> get_delta_max( vector_t<double> & q )
  {
    double delta = 0.5 * ( delta / ( (double)N ) );
    double dx = +delta * ( N / 2 );
    double dy = +delta * ( N / 2 );
    dx = dx+q[0] > border ? border - q[0] : dx;
    dy = dy+q[1] > border ? border - q[1] : dy;
    vector_t<double> const coord { dx, dy };
    return coord;
  }


  void set( double delta_new, double sigma_new, double border_new )
  {
    set_delta( delta_new );
    set_sigma( sigma_new );
    set_border( border_new );
    vector_t<double> o = { 0, 0 };
    for ( int i = 0; i < N; i++ )
    {
      for ( int j = 0; j < N; ++j )
      {
        vector_t<double> q = slot_to_coord( i, j );
        K_field[i][j] = field( q, o );
        K_gradx[i][j] = gradx( q, o );
        K_grady[i][j] = grady( q, o );
      }
    }
  }

  double field( vector_t<double> & q, vector_t<double> & m )
  {
    vector_t<double> r = q - m;
    double const den = 2 * M_PI * pow( sigma, 2 );
    double const exponent = -r.squared_norm2() / ( 2 * pow( sigma, 2 ) );
    double const measure = pow( delta, 2 );
    double const res = exp( exponent ) * measure / den;
    return res;
  }

  double gradx( vector_t<double> & q, vector_t<double> & m )
  {
    vector_t<double> r = q - m;
    return field( q, m )*( r[0]/pow( sigma, 2 ) );
  }

  double grady( vector_t<double> & q, vector_t<double> & m )
  {
    vector_t<double> r = q - m;
    return field( q, m )*( r[1]/pow( sigma, 2 ) );
  }

  std::tuple<int, int> get_slot( vector_t<double> & q, vector_t<double> & m )
  {
    vector_t<double> dq = q - m;

    int i = int( round( dq[0]/delta) ) + N / 2;
    int j = int( round( dq[1]/delta) ) + N / 2;
    return std::make_tuple( i, j );
  }

  bool is_out( int i, int j )
  {
    return ( i < 0 ) || ( i >= N ) || ( j < 0 ) || ( j >= N );
  }

  double kernel_field( vector_t<double> & q, vector_t<double> & m )
  {
    auto [i, j] = get_slot( q, m );
    if ( is_out( i, j )  )
      return 0;
    return K_field[i][j];
  }

  double kernel_gradx( vector_t<double> & q, vector_t<double> & m )
  {
    auto [i, j] = get_slot( q, m );
    if ( is_out( i, j )  )
      return 0;
    return K_gradx[i][j];
  }

  double kernel_grady( vector_t<double> & q, vector_t<double> & m )
  {
    auto [i, j] = get_slot( q, m );
    if ( is_out( i, j )  )
      return 0;
    return K_grady[i][j];
  }

  vector_t<double> kernel_grads( vector_t<double> & q, vector_t<double> & m )
  {
    auto [i, j] = get_slot( q, m );
    return vector_t<double>{ K_gradx[i][j], K_grady[i][j] };
  }

  void print_field()
  {
    std::cout << "K_field" << std::endl;

    for ( int y = N-1; y >= 0; --y )
    {
      for ( int x = 0; x < N; ++x )
      {
        std::cout << K_field[x][y] << " ";
      }      
      std::cout << std::endl;
    }
  }

  void print_gradx()
  {
    std::cout << "K_gradx" << std::endl;
    for ( int y = N-1; y >= 0; --y )
    {
      for ( int x = 0; x < N; ++x )
      {
        std::cout << K_gradx[x][y] << " ";
      }      
      std::cout << std::endl;
    }
  }

  void print_grady()
  {
    std::cout << "K_grady" << std::endl;

    for ( int y = N-1; y >= 0; --y )
    {
      for ( int x = 0; x < N; ++x )
      {
        std::cout << K_grady[x][y] << " ";
      }      
      std::cout << std::endl;
    }
  }


private:
  /*! \brief position vector */
  std::array<std::array<double, N>, N> K_field;
  std::array<std::array<double, N>, N> K_gradx;
  std::array<std::array<double, N>, N> K_grady;
  double sigma;
  double delta;
  double border;
};

} // namespace mockturtle

