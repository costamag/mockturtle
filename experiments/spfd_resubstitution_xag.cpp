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
    ps.use_dont_cares = true;
    ps.max_pis = 8;
    ps.max_divisors = std::numeric_limits<uint32_t>::max();

    const uint32_t size_before = xag.num_gates();
    spfd_resubstitution( xag, ps, &st );
    xag = cleanup_dangling( xag );

    const auto cec = benchmark == "hyp" ? true : abc_cec( xag, benchmark );

    exp( benchmark, size_before, size_before - xag.num_gates(), to_seconds( st.time_total ), cec );
  }

  exp.save();
  exp.table();

  return 0;
}