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
#include <ctime>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, float, uint32_t, double, float, bool, bool> exp( "spfd_resubstitution_xag_infinite_ISCAS", "benchmark", "size", "u-size", "u-runtime", "i-size", "i-gain",  "i-runtime", "u-equivalent" , "i-equivalent" );

  for ( auto const& benchmark : iscas_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    xag_network xag;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xag ) ) != lorina::return_code::success )
    {
      continue;
    }

    resubstitution_params ps;
    resubstitution_stats ust;
    resubstitution_stats ist;

    // ps.pattern_filename = "1024sa1/" + benchmark + ".pat";
    ps.max_inserts = 20;
    ps.max_pis = 8;
    ps.progress = true;
    ps.max_divisors = std::numeric_limits<uint32_t>::max();

    const uint32_t size_before = xag.num_gates();

    /* urs x infinite */
    uint32_t size_new = xag.num_gates();
    uint32_t size_old = std::numeric_limits<uint32_t>::max();

    std::clock_t start_simresub;
    double duration_simresub;
    start_simresub = std::clock();

    sim_resubstitution( xag, ps, &ust );
    xag = cleanup_dangling( xag );

    duration_simresub = ( std::clock() - start_simresub ) / (double) CLOCKS_PER_SEC;

    double size_urs = xag.num_gates();
    printf("urs=%d\n", xag.num_gates());
    const auto cecu = benchmark == "hyp" ? true : abc_cec( xag, benchmark );
    
    /* irs x 1 */

    std::clock_t start_spfdresub;
    double duration_spfdresub;
    start_spfdresub = std::clock();


    xag_network xagA;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xagA ) ) != lorina::return_code::success )
    {
      continue;
    }

    gen_resubstitution( xagA, ps, &ist );
    xagA = cleanup_dangling( xagA );

    duration_spfdresub = ( std::clock() - start_spfdresub ) / (double) CLOCKS_PER_SEC;

    double size_irs = xagA.num_gates();

    const auto ceci = benchmark == "hyp" ? true : abc_cec( xagA, benchmark );
    double gain = 100*(size_irs-size_urs)/size_urs;
    printf("irs=%d --> %f%\n", xagA.num_gates(), gain );

    exp( benchmark, size_before, (uint32_t)size_urs, duration_simresub, (uint32_t)size_irs, gain, duration_spfdresub, cecu, ceci );
  }

  exp.save();
  exp.table();

  return 0;
}


//| benchmark | size | u-size | u-runtime | i-size | i-gain | i-runtime | u-equivalent | i-equivalent |
//|       c17 |    6 |      6 |      0.00 |      6 |   0.00 |      0.00 |         true |         true |
//|      c432 |  208 |    167 |      0.01 |    168 |   0.60 |      0.09 |         true |         true |
//|      c499 |  398 |    388 |      0.01 |    353 |  -9.02 |      6.44 |         true |         true |
//|      c880 |  325 |    296 |      0.01 |    279 |  -5.74 |      1.18 |         true |         true |
//|     c1355 |  502 |    420 |      0.02 |    408 |  -2.86 |      6.14 |         true |         true |
//|     c1908 |  341 |    283 |      0.01 |    221 | -21.91 |      2.67 |         true |         true |
//|     c2670 |  716 |    542 |      0.02 |    586 |   8.12 |      3.64 |         true |         true |
//|     c3540 | 1024 |    810 |      0.10 |    876 |   8.15 |      7.78 |         true |         true |
//|     c5315 | 1776 |   1309 |      0.11 |   1515 |  15.74 |      8.96 |         true |         true |
//|     c6288 | 2337 |   1886 |      0.15 |   2294 |  21.63 |    120.59 |         true |         true |
//|     c7552 | 1469 |   1322 |      0.07 |   1256 |  -4.99 |      7.10 |         true |         true |