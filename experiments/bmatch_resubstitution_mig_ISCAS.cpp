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
#include <mockturtle/algorithms/mig_resub.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/views/fanout_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <ctime>

#include <experiments.hpp>

struct experiments_stats_t
{
  uint32_t num_gates{0};
  double time{0};
  double gain{0};
  bool cec{false};
};

template<uint32_t K, uint32_t S, uint32_t I, class Ntk> void bmatch_resub( std::string const& benchmark, Ntk& ntk, experiments_stats_t& stats )
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

  bmatch_resubstitution<K,S,I>( ntk, ps, &st );
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

    using view_t = depth_view<fanout_view<Ntk>>;


  resubstitution_params ps;
  resubstitution_stats st;

  ps.max_inserts = 20;
  ps.max_pis = 8;
  ps.progress = true;
  ps.max_divisors = std::numeric_limits<uint32_t>::max();

  uint32_t size_old = ntk.num_gates();
  uint32_t size_new;

  std::clock_t start = std::clock();

  fanout_view<Ntk> fview{ ntk };
  view_t vntk{ fview };

  //while( size_new < size_old )
  {
    mig_resubstitution( vntk );
    ntk = cleanup_dangling( ntk );
    size_new = ntk.num_gates();
  }
  stats.num_gates = ntk.num_gates();
  stats.time = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  stats.gain = 100*(((double)size_new)/((double)size_old) - 1);
  stats.cec = benchmark == "hyp" ? true : abc_cec( ntk, benchmark );
}


int main()
{
  using namespace experiments;
  using namespace mockturtle;

  static constexpr uint32_t K1 = 6u;
  static constexpr uint32_t S1 = 1u;
  static constexpr uint32_t I1 = 1u;

  static constexpr uint32_t K2 = 4u;
  static constexpr uint32_t S2 = 10u;
  static constexpr uint32_t I2 = 1u;

  std::string e1 = "(SOA)";
  std::string e2 = "(" + std::to_string(K2) + "," + std::to_string(S2) + "," + std::to_string(I2) + ")";

  experiment<std::string, uint32_t, uint32_t, float, uint32_t, double, float, uint32_t, double, float, bool, bool> exp( "bmatch_resubstitution_mig_infinite_ISCAS", "benchmark", "size", "size(u)", "time(u)", "i-size"+e1, "gain"+e1, "time"+e1, "size"+e2, "gain"+e2, "time"+e2, "cec(u)", "cec"+e2 );

  double gain1{0};
  double gain2{0};
  uint32_t cnt{0};

  for ( auto const& benchmark : resub_benchmarks( iscas ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    /* low effort K=4 S=1 I=1 */

    mig_network mig1;


    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( mig1 ) ) != lorina::return_code::success )
    {
      continue;
    }

    auto size0 = mig1.num_gates();

    experiments_stats_t stU;
    infinite_sim_resub( benchmark, mig1, stU );

    experiments_stats_t st1;


    //bmatch_resub<K1,S1,I1>( benchmark, mig1, st1 );
    
    /* high effort K=7 S=10 I=10 */

    mig_network mig2;
    
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( mig2 ) ) != lorina::return_code::success )
    {
      continue;
    }

    //infinite_sim_resub( benchmark, mig2, stU );

    experiments_stats_t st2;

    bmatch_resub<K2,S2,I2>( benchmark, mig2, st2 );

    exp( benchmark, size0, stU.num_gates, stU.time, stU.num_gates, stU.gain, stU.time, st2.num_gates, st2.gain, st2.time, stU.cec, st2.cec );
    printf("[SOA]=%f [%d,%d,%d]=%f\n",  stU.gain, K2, S2, I2, st2.gain );

    gain1 += stU.gain;
    gain2 += st2.gain;
    cnt++;
  }

  exp.save();
  exp.table();

  printf("[SOA]=%f [%d,%d,%d]=%f\n",  gain1/cnt, K2, S2, I2, gain2/cnt );


  return 0;
}


