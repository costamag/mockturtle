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

  experiment<std::string, uint32_t, uint32_t, uint32_t, double, float, float, bool, bool> exp( "spfd_resubstitution_aig_infinite", "benchmark", "size", "u-size", "i-size", "i-gain", "u-runtime", "i-runtime", "u-equivalent" , "i-equivalent" );

  for ( auto const& benchmark : iscas_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
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

    const uint32_t size_before = aig.num_gates();

    /* urs x infinite */
    uint32_t size_new = aig.num_gates();
    uint32_t size_old = std::numeric_limits<uint32_t>::max();
    while( size_new < size_old )
    {
      sim_resubstitution( aig, ps, &ust );
      aig = cleanup_dangling( aig );
      size_old = size_new;
      size_new = aig.num_gates();
    }
    double size_urs = aig.num_gates();
    const auto cecu = benchmark == "hyp" ? true : abc_cec( aig, benchmark );
    
    /* irs x 1 */
    spfd_resubstitution( aig, ps, &ist );
    aig = cleanup_dangling( aig );
    double size_irs = aig.num_gates();

    const auto ceci = benchmark == "hyp" ? true : abc_cec( aig, benchmark );
    double gain = (size_irs-size_urs)/size_urs;

    exp( benchmark, size_before, (uint32_t)size_urs, (uint32_t)size_irs, gain, to_seconds( ust.time_total ), to_seconds( ist.time_total ), cecu, ceci );
  }

  exp.save();
  exp.table();

  return 0;
}