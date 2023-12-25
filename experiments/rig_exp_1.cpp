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

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool> exp( "rig_exp_1", "benchmark", "luts", "lut_depth", "rigs", "rigs_depth", "rs rigs", "rs rigs_depth", "eq(LUT)", "eq(RIG)", "eq(RS)" );

  for ( auto const& benchmark : epfl_benchmarks( ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    lut_map_params lps;
    lps.cut_enumeration_ps.cut_size = 4u;
    lps.cut_enumeration_ps.cut_limit = 8u;
    lps.recompute_cuts = true;
    lps.area_oriented_mapping = true;
    lps.cut_expansion = true;
    lut_map_stats st;
    auto klut = lut_map( aig, lps, &st );

    klut = cleanup_luts(klut);
    depth_view<klut_network> klut_d{ klut };

    auto const cec = benchmark == "hyp" ? true : abc_cec( klut, benchmark );

    std::string tmp = benchmark + "tmp.blif";
    write_blif( klut, tmp );
    
    rig_network rig;
    if ( lorina::read_blif( tmp, blif_reader( rig ) ) != lorina::return_code::success )
    {
      printf("rig unsuccessful\n");
      continue;
    }
    depth_view<rig_network> rig_d{ rig };
    auto const rig_cec = benchmark == "hyp" ? true : abc_cec( rig, benchmark );

    uint32_t rig_num_gates = rig.num_gates();
    uint32_t rig_depth = rig_d.depth();

    resubstitution_params rps;
    resubstitution_stats rst;
    // ps.pattern_filename = "1024sa1/" + benchmark + ".pat";
    rps.max_inserts = 20;
    rps.max_trials = 1;
    rps.max_pis = 10;
    rps.max_divisors = std::numeric_limits<uint32_t>::max();

    sim_resubstitution( rig, rps, &rst );
    rig = cleanup_dangling( rig );
    depth_view<rig_network> rs_rig_d{ rig };

    uint32_t rs_rig_num_gates = rig.num_gates();
    uint32_t rs_rig_depth = rs_rig_d.depth();

    const auto cec_rs = benchmark == "hyp" ? true : abc_cec( rig, benchmark );

    exp( benchmark, klut.num_gates(), klut_d.depth(), rig_num_gates, rig_depth, rs_rig_num_gates, rs_rig_depth, cec, rig_cec, cec_rs );
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