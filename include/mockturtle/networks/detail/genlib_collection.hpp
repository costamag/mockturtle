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
  \file genlib_collection.hpp
  \brief For each functor utilities

  \author Andrea Costamagna
*/

#pragma once

#include <string>


namespace mockturtle::rils::detail
{
  std::string const mcnc_library =  "GATE   inv1    1  O=!a;             PIN * INV 1 999 0.9 0.3 0.9 0.3\n"
                                    "GATE   inv2    2  O=!a;             PIN * INV 2 999 1.0 0.1 1.0 0.1\n"
                                    "GATE   inv3    3  O=!a;             PIN * INV 3 999 1.1 0.09 1.1 0.09\n"
                                    "GATE   inv4    4  O=!a;             PIN * INV 4 999 1.2 0.07 1.2 0.07\n"
                                    "GATE   nand2   2  O=!(a*b);         PIN * INV 1 999 1.0 0.2 1.0 0.2\n"
                                    "GATE   nand3   3  O=!(a*b*c);       PIN * INV 1 999 1.1 0.3 1.1 0.3\n"
                                    "GATE   nand4   4  O=!(a*b*c*d);     PIN * INV 1 999 1.4 0.4 1.4 0.4\n"
                                    "GATE   nor2    2  O=!(a+b);         PIN * INV 1 999 1.4 0.5 1.4 0.5\n"
                                    "GATE   nor3    3  O=!(a+b+c);       PIN * INV 1 999 2.4 0.7 2.4 0.7\n"
                                    "GATE   nor4    4  O=!(a+b+c+d);     PIN * INV 1 999 3.8 1.0 3.8 1.0\n"
                                    "GATE   and2    3  O=a*b;            PIN * NONINV 1 999 1.9 0.3 1.9 0.3\n"
                                    "GATE   or2     3  O=a+b;            PIN * NONINV 1 999 2.4 0.3 2.4 0.3\n"
                                    "GATE   xor2a   5  O=a*!b+!a*b;      PIN * UNKNOWN 2 999 1.9 0.5 1.9 0.5\n"
                                    "#GATE  xor2b   5  O=!(a*b+!a*!b);   PIN * UNKNOWN 2 999 1.9 0.5 1.9 0.5\n"
                                    "GATE   xnor2a  5  O=a*b+!a*!b;      PIN * UNKNOWN 2 999 2.1 0.5 2.1 0.5\n"
                                    "#GATE  xnor2b  5  O=!(a*!b+!a*b);   PIN * UNKNOWN 2 999 2.1 0.5 2.1 0.5\n"
                                    "GATE   aoi21   3  O=!(a*b+c);       PIN * INV 1 999 1.6 0.4 1.6 0.4\n"
                                    "GATE   aoi22   4  O=!(a*b+c*d);     PIN * INV 1 999 2.0 0.4 2.0 0.4\n"
                                    "GATE   oai21   3  O=!((a+b)*c);     PIN * INV 1 999 1.6 0.4 1.6 0.4\n"
                                    "GATE   oai22   4  O=!((a+b)*(c+d)); PIN * INV 1 999 2.0 0.4 2.0 0.4\n"
                                    "GATE   buf     2  O=a;              PIN * NONINV 1 999 1.0 0.0 1.0 0.0\n"
                                    "GATE   zero    0  O=CONST0;\n"
                                    "GATE   one     0  O=CONST1;";

  std::string const aig_library = "GATE   and2    1  O=a*b;            PIN * NONINV 1 999 1.0 0.0 1.0 0.0\n"
                                  "GATE   inv     0  O=!a;             PIN * INV 1 999 0.0 0.0 0.0 0.0\n"
                                  "GATE   buf     0  O=a;              PIN * NONINV 1 999 0.0 0.0 0.0 0.0\n"
                                  "GATE   zero    0  O=CONST0;\n"
                                  "GATE   one     0  O=CONST1;";

  std::string const xaig_library =  "GATE   and2    1  O=a*b;            PIN * NONINV 1 999 1.0 0.0 1.0 0.0\n"
                                    "GATE   inv     0  O=!a;             PIN * INV 1 999 0.0 0.0 0.0 0.0\n"
                                    "GATE   buf     0  O=a;              PIN * NONINV 1 999 0.0 0.0 0.0 0.0\n"
                                    "GATE   xor2    1  O=a*!b+!a*b;      PIN * NONINV 1 999 1.0 0.0 1.0 0.0\n"
                                    "GATE   zero    0  O=CONST0;\n"
                                    "GATE   one     0  O=CONST1;";

  std::string const mig_library =   "GATE   mig     1  O=(a*b)+(a*c)+(b*c);    PIN * NONINV 1 999 1.0 0.0 1.0 0.0\n"
                                    "GATE   inv     0  O=!a;             PIN * INV 1 999 0.0 0.0 0.0 0.0\n"
                                    "GATE   buf     0  O=a;              PIN * NONINV 1 999 0.0 0.0 0.0 0.0\n"
                                    "GATE   zero    0  O=CONST0;\n"
                                    "GATE   one     0  O=CONST1;";

} /* namespace mockturtle::rils::detail */