/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2019  EPFL
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
#include <lorina/blif.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/rig.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/mapping_view.hpp>
#include <mockturtle/algorithms/sim_resub.hpp>


#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;
  using namespace rils;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, double, bool> exp( "rig_exp_2", "benchmark", "rigs0", "rigs0_depth", "rigs1", "rigs1_depth", "t(RS)", "eq(RS)" );

  for ( auto const& benchmark : epfl_benchmarks( ~experiments::sin ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    rig_network rig0;
    auto path = benchmark + "_size.blif";

    if ( lorina::read_blif( path, blif_reader( rig0 ) ) != lorina::return_code::success )
    {
      printf("rig0 unsuccessful\n");
      continue;
    }
    rig0 = cleanup_dangling( rig0 );
    depth_view<rig_network> rig0_d{ rig0 };

    rig_network rig1;
    if ( lorina::read_blif( path, blif_reader( rig1 ) ) != lorina::return_code::success )
    {
      printf("rig1 unsuccessful\n");
      continue;
    }

    resubstitution_params rps;
    resubstitution_stats rst;
    // ps.pattern_filename = "1024sa1/" + benchmark + ".pat";
    rps.max_inserts = 20;
    rps.max_trials = 100;
    rps.max_pis = 10;
    rps.max_divisors = std::numeric_limits<uint32_t>::max();

    sim_resubstitution( rig1, rps, &rst );
    rig1 = cleanup_dangling( rig1 );
    depth_view<rig_network> rig1_d{ rig1 };

    const auto cec1 = benchmark == "hyp" ? true : abc_cec( rig1, benchmark );

    exp( benchmark, rig0.num_gates(), rig0_d.depth(), rig1.num_gates(), rig1_d.depth(), to_seconds(rst.time_total), cec1 );
  }

  exp.save();
  exp.table();

  return 0;
}


//| benchmark | luts | lut_depth | rigs | rigs_depth | rs rigs | rs rigs_depth | eq(LUT) | eq(RIG) | eq(RS) |
//|       c17 |    6 |         3 |    6 |          3 |       6 |             3 |    true |    true |   true |
//|      c432 |  172 |        24 |  172 |         24 |     171 |            25 |    true |    true |   true |
//|      c499 |  190 |        13 |  190 |         13 |     190 |            13 |    true |    true |   true |
//|      c880 |  271 |        24 |  271 |         24 |     271 |            24 |    true |    true |   true |
//|     c1355 |  190 |        13 |  190 |         13 |     190 |            13 |    true |    true |   true |
//|     c1908 |  184 |        19 |  184 |         19 |     166 |            16 |    true |    true |   true |
//|     c2670 |  582 |        19 |  566 |         19 |     561 |            19 |    true |    true |   true |
//|     c3540 |  917 |        37 |  914 |         37 |     890 |            37 |    true |    true |   true |
//|     c5315 | 1508 |        27 | 1489 |         27 |    1427 |            27 |    true |    true |   true |
//|     c6288 | 1408 |        73 | 1408 |         73 |    1407 |            73 |    true |    true |   true |
//|     c7552 | 1093 |        24 | 1092 |         24 |    1085 |            24 |    true |    true |   true |