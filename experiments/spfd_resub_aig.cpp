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

  uint32_t gain_soa{0};
  uint32_t gain_spfd{0};
  uint32_t cnt{0};

  experiment<std::string, uint32_t, uint32_t, uint32_t, float, float, bool, bool> exp( "spfd_aig", "benchmark", "size", "gates(SOA)", "gates(SPFD)", "time(SOA)", "time(SPFD)", "eq(SOA)", "eq(SPFD)" );

  for ( auto const& benchmark : resub_benchmarks( iscas | epfl ) )//experiments::c499
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
    ps_soa.max_pis = 8;
    ps_soa.max_trials = 100;
    ps_soa.max_divisors = std::numeric_limits<uint32_t>::max();

    const uint32_t size_before = aig_soa.num_gates();
    sim_resubstitution( aig_soa, ps_soa, &st_soa );
    aig_soa = cleanup_dangling( aig_soa );

    const auto cec_soa = benchmark == "hyp" ? true : abc_cec( aig_soa, benchmark );
    
    #pragma endregion SOA

    #pragma region SPFD
    
    aig_network aig_spfd;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig_spfd ) ) != lorina::return_code::success )
    {
      continue;
    }

    resubstitution_params ps_spfd;
    resubstitution_stats st_spfd;

    // ps.pattern_filename = "1024sa1/" + benchmark + ".pat";
    ps_spfd.max_inserts = 20;
    ps_spfd.max_pis = 8;
    ps_spfd.max_trials = 2;
    ps_spfd.max_divisors = std::numeric_limits<uint32_t>::max();

    static constexpr uint32_t K = 8u;
    static constexpr uint32_t S = 10u;
    static constexpr uint32_t I = 10u;
    static constexpr bool use_bmatch = false;
    static constexpr bool use_greedy = false;
    static constexpr bool use_lsearch = true;

    sim_resubstitution_spfd<K, S, I, use_bmatch, use_greedy, use_lsearch>( aig_spfd, ps_spfd, &st_spfd );
    aig_spfd = cleanup_dangling( aig_spfd );

    const auto cec_spfd = benchmark == "hyp" ? true : abc_cec( aig_spfd, benchmark );
    
    #pragma endregion SPFD

    cnt++;
    gain_soa += size_before - aig_soa.num_gates();
    gain_spfd += size_before - aig_spfd.num_gates();
    printf( "gain(SOA)=%d gain(SPFD)=%d\n", size_before - aig_soa.num_gates(), size_before - aig_spfd.num_gates() );

    exp( benchmark, size_before, aig_soa.num_gates(), aig_spfd.num_gates(), to_seconds( st_soa.time_total ), to_seconds( st_spfd.time_total ), cec_soa, cec_spfd );
  }


  exp.save();
  exp.table();

  return 0;
}
