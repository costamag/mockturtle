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

#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/sim_resub.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/xag.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, float, bool> exp( "spfd_resubstitution_xag", "benchmark", "size", "gain", "runtime", "equivalent" );

  for ( auto const& benchmark : iscas_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    xag_network xag;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xag ) ) != lorina::return_code::success )
    {
      continue;
    }

    resubstitution_params ps;
    resubstitution_stats st;

    // ps.pattern_filename = "1024sa1/" + benchmark + ".pat";
    ps.max_inserts = 20;
    ps.max_pis = 8;
    ps.use_dont_cares = true;
    ps.max_divisors = std::numeric_limits<uint32_t>::max();

    const uint32_t size_before = xag.num_gates();
    sim_resubstitution( xag, ps, &st );
    xag = cleanup_dangling( xag );

    const auto cec = benchmark == "hyp" ? true : abc_cec( xag, benchmark );

    exp( benchmark, size_before, size_before - xag.num_gates(), to_seconds( st.time_total ), cec );
  }

  exp.save();
  exp.table();

  return 0;
}

//4
//|       c17 |    6 |    0 |    0.00 |       true |
//|      c432 |  208 |   40 |    0.04 |       true |
//|      c499 |  398 |   50 |    0.74 |       true |
//|      c880 |  325 |   51 |    0.21 |       true |
//|     c1355 |  502 |  111 |    0.71 |       true |
//|     c1908 |  341 |  118 |    0.34 |       true |
//|     c2670 |  716 |  144 |    0.83 |       true |
//|     c3540 | 1024 |  145 |    1.53 |       true |
//|     c5315 | 1776 |  273 |    2.24 |       true |
//|     c6288 | 2337 |   46 |   36.72 |       true |
//|     c7552 | 1469 |  208 |    1.51 |       true |
//5
//| benchmark | size | gain | runtime | equivalent |
//|       c17 |    6 |    0 |    0.00 |       true |
//|      c432 |  208 |   40 |    0.07 |       true |
//|      c499 |  398 |   62 |    0.83 |       true |
//|      c880 |  325 |   60 |    0.33 |       true |
//|     c1355 |  502 |   91 |    0.98 |       true |
//|     c1908 |  341 |  122 |    0.39 |       true |
//|     c2670 |  716 |  164 |    0.95 |       true |
//|     c3540 | 1024 |  165 |    1.65 |       true |
//|     c5315 | 1776 |  293 |    2.52 |       true |
//|     c6288 | 2337 |   46 |   40.83 |       true |
//|     c7552 | 1469 |  237 |    2.60 |       true |
//6
//|       c17 |    6 |    0 |    0.00 |       true |
//|      c432 |  208 |   40 |    0.09 |       true |
//|      c499 |  398 |   80 |    0.92 |       true |
//|      c880 |  325 |   58 |    0.35 |       true |
//|     c1355 |  502 |  122 |    0.98 |       true |
//|     c1908 |  341 |  136 |    0.56 |       true |
//|     c2670 |  716 |  173 |    1.75 |       true |
//|     c3540 | 1024 |  170 |    2.54 |       true |
//|     c5315 | 1776 |  344 |    4.14 |       true |
//|     c6288 | 2337 |   46 |   47.03 |       true |
//|     c7552 | 1469 |  273 |    2.83 |       true |


