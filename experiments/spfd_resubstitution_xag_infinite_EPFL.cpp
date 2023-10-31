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

struct experiments_stats_t
{
  uint32_t num_gates{0};
  double time{0};
  double gain{0};
  bool cec{false};
};

template<uint32_t K, uint32_t S, uint32_t I, class Ntk> void spfd_resub( std::string const& benchmark, Ntk& ntk, experiments_stats_t& stats )
{
  using namespace mockturtle;
  using namespace experiments;

  double size0 = ntk.num_gates();

  resubstitution_params ps;
  resubstitution_stats st;

  ps.max_inserts = 20;
  ps.max_pis = 8;
  ps.progress = true;
  ps.max_divisors = std::numeric_limits<uint32_t>::max();

  std::clock_t start = std::clock();

  spfd_resubstitution<K,S,I>( ntk, ps, &st );
  ntk = cleanup_dangling( ntk );

  stats.num_gates = ntk.num_gates();
  stats.time = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  stats.gain = 100*( (double)ntk.num_gates()-size0 )/size0;
  stats.cec = benchmark == "hyp" ? true : abc_cec( ntk, benchmark );
}

template<class Ntk>
void infinite_sim_resub( std::string const& benchmark, Ntk& ntk, experiments_stats_t & stats )
{
  using namespace mockturtle;
  using namespace experiments;

  resubstitution_params ps;
  resubstitution_stats st;

  ps.max_inserts = 20;
  ps.max_pis = 8;
  ps.progress = true;
  ps.max_divisors = std::numeric_limits<uint32_t>::max();

  uint32_t size_new = ntk.num_gates();
  uint32_t size_old = std::numeric_limits<uint32_t>::max();

  std::clock_t start = std::clock();

  while( size_new < size_old )
  {
    sim_resubstitution( ntk, ps, &st );
    ntk = cleanup_dangling( ntk );
    size_old = size_new;
    size_new = ntk.num_gates();
  }
  stats.num_gates = ntk.num_gates();
  stats.time = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  stats.gain = 0;
  stats.cec = benchmark == "hyp" ? true : abc_cec( ntk, benchmark );
}


int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, float, uint32_t, double, float, uint32_t, double, float, bool, bool, bool> exp( "spfd_resubstitution_xag_infinite_EPFL", "benchmark", "size", "size(u)", "time(u)", "size(4,1,1)", "gain(4,1,1)", "time(4,1,1)", "size(7,10,10)", "gain(7,10,10)", "time(7,10,10)", "cec(u)", "cec(4)", "cec(7)" );

  double gain1{0};
  double gain2{0};
  uint32_t cnt{0};

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    /* low effort K=4 S=1 I=1 */

    xag_network xag1;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xag1 ) ) != lorina::return_code::success )
    {
      continue;
    }
    auto size0 = xag1.num_gates();

    experiments_stats_t stU;
    infinite_sim_resub( benchmark, xag1, stU );

    experiments_stats_t st1;
    spfd_resub<4u,1u,1u>( benchmark, xag1, st1 );
    
    /* high effort K=7 S=10 I=10 */

    xag_network xag2;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xag2 ) ) != lorina::return_code::success )
    {
      continue;
    }

    infinite_sim_resub( benchmark, xag2, stU );

    experiments_stats_t st2;
    spfd_resub<7u,10u,100u>( benchmark, xag2, st2 );

    printf("[4,1,1]=%f [7,10,100]=%f\n", st1.gain, st2.gain );

    exp( benchmark, size0, stU.num_gates, stU.time, st1.num_gates, st1.gain, st1.time, st2.num_gates, st2.gain, st2.time, stU.cec, st1.cec, st2.cec );

    gain1 += st1.gain;
    gain2 += st2.gain;
    cnt++;
  }

  exp.save();
  exp.table();

  printf("[4,1,1]=%f [7,10,100]=%f\n", gain1/cnt, gain2/cnt );


  return 0;
}
