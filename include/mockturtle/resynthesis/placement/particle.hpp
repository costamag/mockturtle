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

class particle_t
{
  static constexpr double DOUBLE_MAX = std::numeric_limits<double>::max();
  vector_t<double> _0{ 0, 0 };

public:
  particle_t( vector_t<double> & position, vector_t<double> & moment, vector_t<double> & force, double degree, double range, int anchor )
    : q( position ), p( moment ), F( force ), Q( degree ), R( range ), anchored( anchor )
  {
    A = M_PI * pow( R, 2 );
  }

  particle_t( particle_t & other )
    : particle_t( other.get_position(), other.get_moment(), other.get_force(), other.get_charge(), other.get_range(), other.get_anchored() )
  {
  }

  particle_t( particle_t && other )
    : particle_t( other.get_position(), other.get_moment(), other.get_force(), other.get_charge(), other.get_range(), other.get_anchored() )
  {
  }

  particle_t()
    : particle_t ( _0, _0, _0, -1, -1, 0 )
  {
  }
  
  particle_t( vector_t<double> & q )
    : particle_t ( q, _0, _0, -1, -1, 0 )
  {
  }

  bool is_valid()
  {
    return Q >= 0;
  }

  void make_invalid()
  {
    Q = -1;
  }

  particle_t operator=( particle_t & rhs)
  {
    if(this == &rhs)
       return *this; // calls copy constructor SimpleCircle(*this)
    q = rhs.get_position();
    p = rhs.get_moment();
    F = rhs.get_force();
    Q = rhs.get_charge();
    A = rhs.get_area();
    R = rhs.get_range();
    C = rhs.get_cost();
    anchored = rhs.get_anchored();
    return *this; // calls copy constructor
  }
  
  particle_t operator=( particle_t && rhs)
  {
    if(this == &rhs)
       return *this; // calls copy constructor SimpleCircle(*this)
    q = rhs.get_position();
    p = rhs.get_moment();
    F = rhs.get_force();
    Q = rhs.get_charge();
    A = rhs.get_area();
    R = rhs.get_range();
    C = rhs.get_cost();
    anchored = rhs.get_anchored();
    return *this; // calls copy constructor
  }

  /*! \brief Returns the position vector */
  vector_t<double> & get_position()
  {
    return q;
  }

  /*! \brief Returns the position vector */
  vector_t<double> & get_moment()
  {
    return p;
  }

  /*! \brief Returns the position vector */
  vector_t<double> & get_force()
  {
    return F;
  }

  /*! \brief Returns the position vector */
  double & get_charge()
  {
    return Q;
  }

  /*! \brief Returns the position vector */
  double & get_range()
  {
    return R;
  }

  int & get_anchored()
  {
    return anchored;
  }

  /*! \brief Returns the position vector */
  void set_range( double const& Rnew )
  {
    R = Rnew;
    A = M_PI * pow( R, 2 );
  }

  void set_charge( double const& Qnew )
  {
    Q = Qnew;
  }

  void set_position( vector_t<double> const& qnew )
  {
    q = qnew;
  }

  /*! \brief Returns the position vector */
  double & get_area()
  {
    return A;
  }

  /*! \brief Returns the position vector */
  double & get_cost()
  {
    return C;
  }

  void reset_moment()
  {
    p[0] = 0;
    p[1] = 0;
  }

  void reset_force()
  {
    F[0] = 0;
    F[1] = 0;
  }

  void add_charge( double dQ )
  {
    Q += dQ;
  }

  bool is_anchored()
  {
    return anchored > 0;
  }

  /* pi 1 po 2*/
  void make_pi()
  {
    anchored = 1;
  }

  void make_po()
  {
    anchored = 2;
  }


  bool is_pi()
  {
    return ( anchored == 1 );
  }

  bool is_po()
  {
    return ( anchored == 2 );
  }

private:
  /*! \brief position vector */
  vector_t<double> q;
  /*! \brief position vector */
  vector_t<double> p;
  /*! \brief grads : -grads of the cost */
  vector_t<double> F;
  /*! \brief charge */
  double Q;
  /*! \brief charge */
  double A;
  /*! \brief range */
  double R;
  /*! \brief contribution to the cost */
  double C;
  int anchored = 0;
};

} // namespace mockturtle

