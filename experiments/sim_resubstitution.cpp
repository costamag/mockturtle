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
#include <mockturtle/networks/aig.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, float, bool> exp( "sim_resubstitution", "benchmark", "size", "gain", "runtime", "equivalent" );

  for ( auto const& benchmark : iscas_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    resubstitution_params ps;
    resubstitution_stats st;

    // ps.pattern_filename = "1024sa1/" + benchmark + ".pat";
    ps.max_inserts = 20;
    ps.max_pis = 8;
    ps.max_divisors = std::numeric_limits<uint32_t>::max();

    const uint32_t size_before = aig.num_gates();
    sim_resubstitution( aig, ps, &st );
    aig = cleanup_dangling( aig );

    const auto cec = benchmark == "hyp" ? true : abc_cec( aig, benchmark );

    exp( benchmark, size_before, size_before - aig.num_gates(), to_seconds( st.time_total ), cec );
  }

  exp.save();
  exp.table();

  return 0;
}
//                   |   sim-resub     | spfd-k4-s1-i1   | spfd-k4-s10-i1  | spfd-k5-s1-i1   | spfd-k5-s1-i10 | spfd-k5-s10-i10 |
//| benchmark | size | gain | runtime  | gain | runtime  | gain | runtime  | gain | runtime  | gain | runtime  |
//|       c17 |    6 |    0 |    0.00  |    0 |    0.00  |    0 |    0.03  |    0 |    0.00  |    0 |    0.00 |    0 |    0.00  |
//|      c432 |  208 |   41 |    0.00  |   40 |    0.01  |   40 |    1.22  |   40 |    0.02  |   40 |    0.02 |   40 |    0.15  |
//|      c499 |  398 |   10 |    0.01  |   43 |    0.73  |   66 |    3.35  |   59 |    0.76  |   65 |    0.88 |   54 |    1.16  |
//|      c880 |  325 |   29 |    0.01  |   45 |    0.14  |   52 |    2.27  |   55 |    0.16  |   56 |    0.15 |   60 |    0.49  |
//|     c1355 |  502 |   82 |    0.01  |   95 |    0.77  |   93 |    4.17  |  106 |    0.88  |  100 |    0.99 |   99 |    1.57  |
//|     c1908 |  341 |   58 |    0.01  |  109 |    0.29  |  105 |    1.89  |  122 |    0.37  |  125 |    0.47 |  120 |    0.54  |
//|     c2670 |  716 |  174 |    0.02  |  122 |    0.57  |  138 |    5.68  |  137 |    0.65  |  158 |    0.88 |  147 |    1.40  |
//|     c3540 | 1024 |  214 |    0.09  |  141 |    0.94  |  151 |   11.44  |  172 |    1.34  |  172 |    1.16 |  172 |    3.11  |
//|     c5315 | 1776 |  467 |    0.07  |  261 |    1.37  |  271 |   20.65  |  278 |    1.69  |  285 |    2.07 |  300 |    5.56  |
//|     c6288 | 2337 |  451 |    0.12  |   20 |   38.37  |   51 |   67.85  |   22 |   46.98  |   39 |   38.85 |   52 |   39.23  |
//|     c7552 | 1469 |  147 |    0.04  |  194 |    1.20  |  188 |   15.51  |  208 |    1.18  |  239 |    1.81 |  242 |    3.88  |

//spfd-k5-s100-i100
//|       c17 |    6 |    0 |    0.16 |       true |
//|      c432 |  208 |   40 |    4.79 |       true |
//|      c499 |  398 |   53 |    4.55 |       true |
//|      c880 |  325 |   60 |    8.63 |       true |
//|     c1355 |  502 |  115 |    5.31 |       true |
//|     c1908 |  341 |  109 |    4.16 |       true |
//|     c2670 |  716 |  171 |   20.36 |       true |
//|     c3540 | 1024 |  161 |   41.12 |       true |
//|     c5315 | 1776 |  325 |   99.46 |       true |
//|     c6288 | 2337 |   55 |   67.68 |       true |
//|     c7552 | 1469 |  235 |   41.91 |       true |


// STATISTICAL SUPPORT SELECTION
//                   |   sim-resub     | spfd-k4-s1-i1   | spfd-k4-s10-i1  |
//| benchmark | size | gain | runtime  | gain | runtime  | gain | runtime  |
//|       c17 |    6 |    0 |    0.00  |    0 |    0.00  |    0 |    0.03  |
//|      c432 |  208 |   41 |    0.00  |   40 |    0.01  |   40 |    1.22  |
//|      c499 |  398 |   10 |    0.01  |   43 |    0.73  |   66 |    3.35  |
//|      c880 |  325 |   29 |    0.01  |   45 |    0.14  |   52 |    2.27  |
//|     c1355 |  502 |   82 |    0.01  |   95 |    0.77  |   93 |    4.17  |
//|     c1908 |  341 |   58 |    0.01  |  109 |    0.29  |  105 |    1.89  |
//|     c2670 |  716 |  174 |    0.02  |  122 |    0.57  |  138 |    5.68  |
//|     c3540 | 1024 |  214 |    0.09  |  141 |    0.94  |  151 |   11.44  |
//|     c5315 | 1776 |  467 |    0.07  |  261 |    1.37  |  271 |   20.65  |
//|     c6288 | 2337 |  451 |    0.12  |   20 |   38.37  |   51 |   67.85  |
//|     c7552 | 1469 |  147 |    0.04  |  194 |    1.20  |  188 |   15.51  |


