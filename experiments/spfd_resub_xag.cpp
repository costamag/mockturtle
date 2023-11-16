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

  double gain_soa{0};
  double gain_spfd{0};

  experiment<std::string, uint32_t, uint32_t, uint32_t, float, float, bool, bool> exp( "spfd_aig", "benchmark", "size", "gates(SOA)", "gates(SPFD)", "time(SOA)", "time(SPFD)", "eq(SOA)", "eq(SPFD)" );

  double cnt{0};

  for ( auto const& benchmark : resub_benchmarks( ))//experiments::c499 ))
  {
    fmt::print( "[i] processing {}\n", benchmark );

    #pragma region SOA
    
    xag_network xag_soa;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xag_soa ) ) != lorina::return_code::success )
    {
      continue;
    }

    resubstitution_params ps_soa;
    resubstitution_stats st_soa;

    // ps.pattern_filename = "1024sa1/" + benchmark + ".pat";
    ps_soa.max_inserts = 20;
    ps_soa.max_pis = 10;
    ps_soa.max_trials = 100;
    ps_soa.max_divisors = std::numeric_limits<uint32_t>::max();

    const uint32_t size_before = xag_soa.num_gates();
    sim_resubstitution( xag_soa, ps_soa, &st_soa );
    xag_soa = cleanup_dangling( xag_soa );

    const auto cec_soa = benchmark == "hyp" ? true : abc_cec( xag_soa, benchmark );
    
    #pragma endregion SOA
    printf("=================\n");
    #pragma region SPFD
    
    xag_network xag_spfd;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xag_spfd ) ) != lorina::return_code::success )
    {
      continue;
    }

    resubstitution_params ps_spfd;
    resubstitution_stats st_spfd;

    // ps.pattern_filename = "1024sa1/" + benchmark + ".pat";
    ps_spfd.max_inserts = 20;
    ps_spfd.max_pis = 10;
    ps_spfd.max_trials = 100;
    ps_spfd.max_divisors = std::numeric_limits<uint32_t>::max();

    static constexpr uint32_t K = 8u;
    static constexpr uint32_t S = 1u;
    static constexpr uint32_t I = 10u;
    static constexpr bool use_bmatch = false;
    static constexpr bool use_greedy = true;
    static constexpr bool use_lsearch = true;

    sim_resubstitution_spfd<K, S, I, use_bmatch, use_greedy, use_lsearch>( xag_spfd, ps_spfd, &st_spfd );
    xag_spfd = cleanup_dangling( xag_spfd );

    const auto cec_spfd = benchmark == "hyp" ? true : abc_cec( xag_spfd, benchmark );
    
    #pragma endregion SPFD

    cnt++;
    gain_soa += (double)(size_before - xag_soa.num_gates())/((double)size_before);
    gain_spfd += (double)(size_before - xag_spfd.num_gates())/((double)size_before);
    printf( "gain(SOA)=%d gain(SPFD)=%d\n", size_before - xag_soa.num_gates(), size_before - xag_spfd.num_gates() );

    exp( benchmark, size_before, xag_soa.num_gates(), xag_spfd.num_gates(), to_seconds( st_soa.time_total ), to_seconds( st_spfd.time_total ), cec_soa, cec_spfd );
    cnt+=1;
  }
  printf("<gain(SOA)>=%.2f <gain(SPFD)>=%.2f\n", 100*gain_soa/cnt, 100*gain_spfd/cnt );

  exp.save();
  exp.table();

  return 0;
}

