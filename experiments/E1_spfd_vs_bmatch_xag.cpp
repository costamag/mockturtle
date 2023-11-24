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

  double gain_bmatch{0};
  double gain_spfd{0};
  double cum_gain_bmatch{0};
  double cum_gain_spfd{0};

  experiment<std::string, uint32_t, float, float, float, float, bool, bool> exp( "spfd_xag", "benchmark", "size", "gain(BMATCH)", "gain(SPFD)", "time(BMATCH)", "time(SPFD)", "eq(BMATCH)", "eq(SPFD)" );

  double cnt{0};

  for ( auto const& benchmark : resub_benchmarks( iscas ))
  {
    fmt::print( "[i] processing {}\n", benchmark );
    printf("BMATCH\n");

    #pragma region BMATCH
    
    xag_network xag_bmatch;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xag_bmatch ) ) != lorina::return_code::success )
    {
      continue;
    }

    double const size_before = xag_bmatch.num_gates();

    resubstitution_params ps_bmatch;
    resubstitution_stats st_bmatch;

    ps_bmatch.max_inserts = 20;
    ps_bmatch.max_pis = 10;
    ps_bmatch.max_trials = 1;
    ps_bmatch.progress = false;
    ps_bmatch.max_divisors = std::numeric_limits<uint32_t>::max();

    static constexpr uint32_t K2 = 10u;
    static constexpr uint32_t K2 = 10u;
    static constexpr uint32_t S2 = 1u;
    static constexpr uint32_t I2 = 1u;

    sim_resubstitution_spfd<K2, S2, I2, true>( xag_bmatch, ps_bmatch, &st_bmatch );
    xag_bmatch = cleanup_dangling( xag_bmatch );

    const auto cec_bmatch = benchmark == "hyp" ? true : abc_cec( xag_bmatch, benchmark );
    
    #pragma endregion BMATCH
    printf("=================\n");
    printf("SPFD\n");
    #pragma region SPFD
    
    xag_network xag_spfd;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xag_spfd ) ) != lorina::return_code::success )
    {
      continue;
    }

    resubstitution_params ps_spfd;
    resubstitution_stats st_spfd;


    ps_spfd.max_inserts = 20;
    ps_spfd.max_pis = 10;
    ps_spfd.max_trials = 1;
    ps_spfd.progress = false;
    ps_spfd.max_divisors = std::numeric_limits<uint32_t>::max();

    static constexpr uint32_t K = 5u;
    static constexpr uint32_t S = 1u;
    static constexpr uint32_t I = 1u;

    sim_resubstitution_spfd<K, S, I, false>( xag_spfd, ps_spfd, &st_spfd );
    xag_spfd = cleanup_dangling( xag_spfd );

    const auto cec_spfd = benchmark == "hyp" ? true : abc_cec( xag_spfd, benchmark );
    
    #pragma endregion SPFD

    cnt++;
    gain_bmatch = (double)(size_before - xag_bmatch.num_gates())/((double)size_before);
    gain_spfd = (double)(size_before - xag_spfd.num_gates())/((double)size_before);

    cum_gain_bmatch += gain_bmatch;
    cum_gain_spfd += gain_spfd;

    printf( "gain(BMATCH)=%f gain(SPFD)=%f\n", size_before - xag_bmatch.num_gates(), size_before - xag_spfd.num_gates() );

    exp( benchmark, size_before, 100*gain_bmatch, 100*gain_spfd, to_seconds( st_bmatch.time_total ), to_seconds( st_spfd.time_total ), cec_bmatch, cec_spfd );
    cnt+=1;
  }
  printf("<gain(BMATCH)>=%.2f <gain(SPFD)>=%.2f\n", 100*cum_gain_bmatch/cnt, 100*cum_gain_spfd/cnt );

  exp.save();
  exp.table();

  return 0;
}
