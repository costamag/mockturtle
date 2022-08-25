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
  \file partial_truth_table.hpp
  \brief Implements partial_truth_table

  \author Andrea Costamagna
*/

#pragma once

#include <kitty/partial_truth_table.hpp>

namespace mockturtle
{

/*! Truth table with resizable, arbitrary number of bits
*/
template<typename Ntk>
struct sim_pattern
{
  /*! \brief Standard constructor.
    \param num_bits Number of bits in use initially
  */
  sim_pattern( kitty::partial_truth_table pat, signal<Ntk> sig = 0, bool simulated = false ) : 
               pat(pat), sig(sig), simulated(simulated) 
  {
  }

  sim_pattern()  
  {
    simulated = false;
    sig = 0;
    pat = kitty::partial_truth_table(1u);
  }

  /*! \brief Constructs a new partial truth table instance with the same number of bits and blocks. */
  inline sim_pattern construct() const
  {
    return sim_pattern( pat, sig, simulated );
  }


  /*! \cond PRIVATE */
public: /* fields */
  kitty::partial_truth_table pat;
  signal<Ntk> sig;
  bool simulated;
  bool flag{false};
  bool flag_sized{false};
  uint32_t layer{0};
  double weight{-1.};
  std::vector<uint32_t> oclass;
  /*! \endcond */
};


} // namespace mockturtle