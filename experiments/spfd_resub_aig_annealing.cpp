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

  double gain_spfd{0};
  double cnt{0};
  std::vector<double> gains;

  for ( auto const& benchmark : resub_benchmarks( experiments::c2670 ))
  {
    fmt::print( "[i] processing {}\n", benchmark );

    #pragma region SPFD
    
    aig_network aig_spfd;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig_spfd ) ) != lorina::return_code::success )
    {
      continue;
    }

    const double size_before = aig_spfd.num_gates();

    resubstitution_params ps_spfd;
    resubstitution_stats st_spfd;

    // ps.pattern_filename = "1024sa1/" + benchmark + ".pat";
    ps_spfd.max_inserts = 20;
    ps_spfd.max_pis = 10;
    ps_spfd.max_trials = 1;
    ps_spfd.progress = true;
    ps_spfd.max_divisors = std::numeric_limits<uint32_t>::max();

    static constexpr uint32_t K = 7u;
    static constexpr uint32_t S = 100u;
    static constexpr uint32_t I = 10u;
    static constexpr bool use_bmatch = false;
    static constexpr bool use_greedy = true;
    static constexpr bool use_lsearch = true;

    sim_resubstitution_spfd<K, S, I, use_bmatch, use_greedy, use_lsearch>( aig_spfd, ps_spfd, &st_spfd );
    aig_spfd = cleanup_dangling( aig_spfd );

    const auto cec_spfd = benchmark == "hyp" ? true : abc_cec( aig_spfd, benchmark );
    
    #pragma endregion SPFD

    cnt++;
    gain_spfd = 100*(size_before - aig_spfd.num_gates())/(size_before);

    gains.push_back( gain_spfd );
    printf("%f\n", gain_spfd );
    cnt+=1;
  }
  printf("[");
  for( int i{0}; i<gains.size()-1; ++i )
  {
    printf("%.2f, ", gains[i] );
  }
  printf("%.2f ]", gains.back() );


  return 0;
}

