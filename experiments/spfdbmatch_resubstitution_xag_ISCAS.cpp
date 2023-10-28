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
    xag_network xagA;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xagA ) ) != lorina::return_code::success )
    {
      continue;
    }

    xag_network xagB;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xagB ) ) != lorina::return_code::success )
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

    const uint32_t size_before = xagA.num_gates();

    /* urs x infinite */
    uint32_t size_new = xagA.num_gates();
    uint32_t size_old = std::numeric_limits<uint32_t>::max();

    std::clock_t start_simresub;
    double duration_simresub;
    start_simresub = std::clock();

    sim_resubstitution( xagA, ps, &ust );
    xagA = cleanup_dangling( xagA );

    duration_simresub = ( std::clock() - start_simresub ) / (double) CLOCKS_PER_SEC;

    double size_urs = xagA.num_gates();
    printf("urs=%d\n", xagA.num_gates());
    const auto cecu = benchmark == "hyp" ? true : abc_cec( xagA, benchmark );
    
    /* irs x 1 */

    std::clock_t start_spfdresub;
    double duration_spfdresub;
    start_spfdresub = std::clock();
    //int it{10};
    //while(it-->0)
    //{
    bmatch_resubstitution( xagB, ps, &ist );
    xagB = cleanup_dangling( xagB );
    //}
    duration_spfdresub = ( std::clock() - start_spfdresub ) / (double) CLOCKS_PER_SEC;

    double size_irs = xagB.num_gates();

    const auto ceci = benchmark == "hyp" ? true : abc_cec( xagB, benchmark );
    double gain = 100*(size_irs-size_urs)/size_urs;
    printf("irs=%d --> %f%\n", xagB.num_gates(), gain );

    exp( benchmark, size_before, (uint32_t)size_urs, duration_simresub, (uint32_t)size_irs, gain, duration_spfdresub, cecu, ceci );
  }

  exp.save();
  exp.table();

  return 0;
}
