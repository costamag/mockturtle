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
#include <mockturtle/algorithms/network_analyzer.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/xag.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, float, uint32_t, float, float> exp( "sim_resub_he_xag", "benchmark", "size", "#LMFFC", "rs", "t(A)", "hers", "t(B)", "d(gates)"  );

  for ( auto const& benchmark : all_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    
    xag_network xag;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xag ) ) != lorina::return_code::success )
    {
      continue;
    }
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

    analyzer_params anPs;
    analyzer_stats anSt;
    // ps.pattern_filename = "1024sa1/" + benchmark + ".pat";
    anPs.max_inserts = 20;
    anPs.max_pis = 8;
    anPs.max_divisors = std::numeric_limits<uint32_t>::max();
    default_analyzer( xag, anPs, &anSt );

    uint32_t nLargeMffc = anSt.nXXLMFFC;
    const uint32_t size_before = xag.num_gates();

    resubstitution_params psA;
    resubstitution_stats stA;
    psA.max_inserts = 20;
    psA.max_pis = 8;
    psA.max_divisors = std::numeric_limits<uint32_t>::max();
    sim_resubstitution( xagA, psA, &stA );
    xagA = cleanup_dangling( xagA );
    float nA = xagA.num_gates();
    float timeA = to_seconds( stA.time_total );
    const auto cecA = benchmark == "hyp" ? true : abc_cec( xagA, benchmark );
    assert(cecA);

    resubstitution_params psB;
    resubstitution_stats stB;
    psB.max_inserts = 20;
    psB.max_pis = 8;
    psB.max_divisors = std::numeric_limits<uint32_t>::max();
    psB.useInfo = true;
    sim_resubstitution( xagB, psB, &stB );
    xagB = cleanup_dangling( xagB );
    float nB = xagB.num_gates();
    float timeB = to_seconds( stB.time_total );

    const auto cecB = benchmark == "hyp" ? true : abc_cec( xagB, benchmark );
    assert(cecB);
    float deltaG = 100*(nB-nA)/nA;

    float deltaT = (timeB-timeA)/(timeA);

    exp( benchmark, size_before, nLargeMffc, xagA.num_gates(), timeA , xagB.num_gates(), timeB, deltaG );
  }

  exp.save();
  exp.table();

  return 0;
}
