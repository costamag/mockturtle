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


  for ( auto const& benchmark : resub_benchmarks( iscas ))
  {
    fmt::print( "[i] processing {}\n", benchmark );

    #pragma region SOA
    
    aig_network aig_soa;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig_soa ) ) != lorina::return_code::success )
    {
      continue;
    }

    resubstitution_params ps_soa;
    resubstitution_stats st_soa;

    // ps.pattern_filename = "1024sa1/" + benchmark + ".pat";
    ps_soa.max_inserts = 20;
    ps_soa.max_pis = 10;
    ps_soa.max_trials = 1;
    ps_soa.max_divisors = std::numeric_limits<uint32_t>::max();

    const uint32_t size_before = aig_soa.num_gates();
    sim_resubstitution( aig_soa, ps_soa, &st_soa );
    aig_soa = cleanup_dangling( aig_soa );

    const auto cec_soa = benchmark == "hyp" ? true : abc_cec( aig_soa, benchmark );
    if( !cec_soa )
    {
      printf("[e] not equivalent\n");
    }
    
    aig_network aig_spfd;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig_spfd ) ) != lorina::return_code::success )
    {
      continue;
    }

    //printf( ".g %d\n", size_before - aig_soa.num_gates() );

  }

  return 0;
}

