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
#include <mockturtle/algorithms/rewrite.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/topo_view.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  double gain_rs{0};
  double gain_rw{0};
  double gain_spfd{0};
  double gain_bmatch{0};
  double cum_gain_rs{0};
  double cum_gain_rw{0};
  double cum_gain_spfd{0};
  double cum_gain_bmatch{0};

  std::vector<uint32_t> GATES_RS;
  std::vector<uint32_t> GATES_BM;
  std::vector<uint32_t> GATES_SP;
  std::vector<uint32_t> GATES_N0;

  experiment<std::string, uint32_t, uint32_t, float, uint32_t, float, uint32_t, float, uint32_t, float, bool, bool, bool, bool> exp( "spfd_xag", "benchmark", "size", "gates(RS)", "time(RS)", "gates(BMATCH)", "time(BMATCH)", "gates(SPFD)", "time(SPFD)", "gates(RW)", "time(RW)",  "eq(RS)", "eq(RW)", "eq(BMATCH)", "eq(SPFD)" );

  double cnt{0};

  static constexpr uint32_t S = 1u;
  static constexpr uint32_t I = 1u;
  static constexpr uint32_t N = 1u;
  static constexpr uint32_t Ks = 10u;
  static constexpr uint32_t Kb = 6u;

  xag_npn_resynthesis<xag_network, xag_network, xag_npn_db_kind::xag_complete> resyn;
  exact_library_params eps;
  eps.np_classification = false;
  eps.compute_dc_classes = true;
  exact_library<xag_network, decltype( resyn )> exact_lib( resyn, eps );

  for ( auto const& benchmark : all_benchmarks( iscas | epfl | iwls))
  {
    fmt::print( "[i] processing {}\n", benchmark );

    #pragma region RS
    
    xag_network xag_rs;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xag_rs ) ) != lorina::return_code::success )
    {
      continue;
    }
    if( xag_rs.num_gates() > 300000 ) continue;

    resubstitution_params ps_rs;
    resubstitution_stats st_rs;

    // ps.pattern_filename = "1024sa1/" + benchmark + ".pat";
    ps_rs.max_inserts = 20;
    ps_rs.max_pis = Ks;
    ps_rs.max_trials = N;
    ps_rs.use_dont_cares = true;
    ps_rs.max_divisors = std::numeric_limits<uint32_t>::max();

    const uint32_t size_before = xag_rs.num_gates();

    GATES_N0.push_back( size_before );

    sim_resubstitution( xag_rs, ps_rs, &st_rs );
    xag_rs = cleanup_dangling( xag_rs );

    const auto cec_rs = benchmark == "hyp" ? true : abc_cec( xag_rs, benchmark );
    
    #pragma endregion RS
    
    #pragma region RW
    
    xag_network xag_rw;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xag_rw ) ) != lorina::return_code::success )
    {
      continue;
    }

    rewrite_params ps_rw;
    rewrite_stats st_rw;
    ps_rw.use_dont_cares = true;

    rewrite( xag_rw, exact_lib, ps_rw, &st_rw );
    xag_rw = cleanup_dangling( xag_rw );

    const auto cec_rw = benchmark == "hyp" ? true : abc_cec( xag_rw, benchmark );
    
    #pragma endregion RW
    
    printf("=================\n");

    #pragma region BMATCH
    xag_network xag_bmatch;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xag_bmatch ) ) != lorina::return_code::success )
    {
      continue;
    }

    resubstitution_params ps_bmatch;
    resubstitution_stats st_bmatch;


    ps_bmatch.max_inserts = 20;
    ps_bmatch.max_pis = Ks;
    ps_bmatch.max_trials = N;
    ps_bmatch.progress = true;
    ps_bmatch.use_dont_cares = true;
    ps_bmatch.max_divisors = std::numeric_limits<uint32_t>::max();


    sim_resubstitution_spfd<Kb, S, I, true>( xag_bmatch, ps_bmatch, &st_bmatch );
    xag_bmatch = cleanup_dangling( xag_bmatch );

    const auto cec_bmatch = benchmark == "hyp" ? true : abc_cec( xag_bmatch, benchmark );
    #pragma endregion BMATCH

    #pragma region SPFD
    xag_network xag_spfd;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xag_spfd ) ) != lorina::return_code::success )
    {
      continue;
    }

    resubstitution_params ps_spfd;
    resubstitution_stats st_spfd;


    ps_spfd.max_inserts = 20;
    ps_spfd.max_pis = Ks;
    ps_spfd.max_trials = N;
    ps_spfd.progress = true;
    ps_spfd.use_dont_cares = true;
    ps_spfd.max_divisors = std::numeric_limits<uint32_t>::max();

    sim_resubstitution_spfd<Kb, S, I, false>( xag_spfd, ps_spfd, &st_spfd );
    xag_spfd = cleanup_dangling( xag_spfd );

    const auto cec_spfd = benchmark == "hyp" ? true : abc_cec( xag_spfd, benchmark );
    #pragma endregion SPFD

    cnt++;
    gain_rs = (double)(size_before - xag_rs.num_gates())/((double)size_before);
    gain_rw = (double)(size_before - xag_rw.num_gates())/((double)size_before);
    gain_spfd = (double)(size_before - xag_spfd.num_gates())/((double)size_before);
    gain_bmatch = (double)(size_before - xag_bmatch.num_gates())/((double)size_before);

    GATES_RS.push_back(xag_rs.num_gates());
    GATES_BM.push_back(xag_bmatch.num_gates());
    GATES_SP.push_back(xag_spfd.num_gates());

    cum_gain_rs += gain_rs;
    cum_gain_rw += gain_rw;
    cum_gain_spfd += gain_spfd;
    cum_gain_bmatch += gain_bmatch;

    printf( "gates(RS)=%d gates(RW)=%d gates(BMATCH) = %d gates(SPFD)=%d\n", xag_rs.num_gates(), xag_rw.num_gates(), xag_bmatch.num_gates(), xag_spfd.num_gates() );
    //experiment<std::string, uint32_t, float, float, float, float, float, float, float, float, bool, bool, bool, bool> exp( "spfd_xag", "benchmark", "size", "gain(RS)", "time(SOA)", "gain(BMATCH)", "time(BMATCH)", "gain(SPFD)", "time(SPFD)", "gain(RW)", "time(RW)",  "eq(RS)", "eq(RW)", "eq(BMATCH)", "eq(SPFD)" );
    exp( benchmark, size_before, xag_rs.num_gates(), to_seconds( st_rs.time_total ),  xag_bmatch.num_gates(), to_seconds( st_bmatch.time_total ), xag_spfd.num_gates(),   to_seconds( st_spfd.time_total ), xag_rw.num_gates(), to_seconds( st_rw.time_total ),cec_rs, cec_rw, cec_bmatch, cec_spfd );
    cnt+=1;
  }

  exp.save();
  exp.table();

  printf("gates_rs=np.array([");
  for( auto x : GATES_RS )
  {
    printf("%d, ", x );
  }
  printf("])\n");

  printf("gates_bmatch=np.array([");
  for( auto x : GATES_BM )
  {
    printf("%d, ", x );
  }
  printf("])\n");

  printf("gates_spfd=np.array([");
  for( auto x : GATES_SP )
  {
    printf("%d, ", x );
  }
  printf("])\n");

  printf("gates_0=np.array([");
  for( auto x : GATES_N0 )
  {
    printf("%d, ", x );
  }
  printf("])\n");


  return 0;
}
//| benchmark | size | gain(RS) | gain(RW) | gain(BMATCH) | gain(SPFD) | time(RS) | time(RW) | time(BMATCH) | time(SPFD) | eq(RS) | eq(RW) | eq(BMATCH) | eq(SPFD) |
//|       c17 |    6 |     0.00 |     0.00 |         0.00 |       0.00 |     0.00 |     0.00 |         0.00 |       0.00 |   true |   true |       true |     true |
//|      c432 |  208 |    19.71 |    20.19 |        19.23 |      19.23 |     0.00 |     0.00 |         0.04 |       0.04 |   true |   true |       true |     true |
//|      c499 |  398 |     1.51 |     1.51 |         1.01 |       0.75 |     0.01 |     0.02 |         0.14 |       0.13 |   true |   true |       true |     true |
//|      c880 |  325 |     5.85 |     4.00 |         4.00 |       3.69 |     0.00 |     0.01 |         0.09 |       0.08 |   true |   true |       true |     true |
//|     c1355 |  502 |     9.16 |    21.91 |        19.92 |       8.17 |     0.01 |     0.01 |         0.11 |       0.15 |   true |   true |       true |     true |
//|     c1908 |  341 |    15.84 |     6.74 |        15.54 |      14.37 |     0.01 |     0.01 |         0.10 |       0.16 |   true |   true |       true |     true |
//|     c2670 |  716 |    21.23 |    20.53 |        19.27 |      16.34 |     0.02 |     0.01 |         0.38 |       0.33 |   true |   true |       true |     true |
//|     c3540 | 1024 |    18.46 |    11.33 |        18.85 |      15.72 |     0.04 |     0.02 |         0.25 |       0.28 |   true |   true |       true |     true |
//|     c5315 | 1776 |    21.06 |    20.61 |        22.18 |      12.95 |     0.03 |     0.04 |         0.74 |       0.82 |   true |  false |       true |     true |
//|     c6288 | 2337 |    19.21 |    19.34 |        18.14 |       0.47 |     0.06 |     0.06 |         0.11 |       0.08 |   true |   true |       true |     true |
//|     c7552 | 1469 |     6.74 |     3.47 |         6.19 |       4.08 |     0.02 |     0.04 |         0.81 |       0.79 |   true |   true |       true |     true |


