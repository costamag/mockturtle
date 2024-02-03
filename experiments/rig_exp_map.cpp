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
#include <lorina/genlib.hpp>
#include <mockturtle/algorithms/mapper_rig.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/rig.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/algorithms/sim_resub.hpp>
#include <mockturtle/algorithms/representation_switching.hpp>
#include <mockturtle/networks/detail/genlib_collection.hpp>


#include <experiments.hpp>

//std::string const mcnc_library = "GATE   inv1    1  O=!a;             PIN * INV 1 999 0.9 0.3 0.9 0.3\n"
//                                 "GATE   inv2    2  O=!a;             PIN * INV 2 999 1.0 0.1 1.0 0.1\n"
//                                 "GATE   inv3    3  O=!a;             PIN * INV 3 999 1.1 0.09 1.1 0.09\n"
//                                 "GATE   inv4    4  O=!a;             PIN * INV 4 999 1.2 0.07 1.2 0.07\n"
//                                 "GATE   nand2   2  O=!(a*b);         PIN * INV 1 999 1.0 0.2 1.0 0.2\n"
//                                 "GATE   nand3   3  O=!(a*b*c);       PIN * INV 1 999 1.1 0.3 1.1 0.3\n"
//                                 "GATE   nand4   4  O=!(a*b*c*d);     PIN * INV 1 999 1.4 0.4 1.4 0.4\n"
//                                 "GATE   nor2    2  O=!(a+b);         PIN * INV 1 999 1.4 0.5 1.4 0.5\n"
//                                 "GATE   nor3    3  O=!(a+b+c);       PIN * INV 1 999 2.4 0.7 2.4 0.7\n"
//                                 "GATE   nor4    4  O=!(a+b+c+d);     PIN * INV 1 999 3.8 1.0 3.8 1.0\n"
//                                 "GATE   and2    3  O=a*b;            PIN * NONINV 1 999 1.9 0.3 1.9 0.3\n"
//                                 "GATE   or2     3  O=a+b;            PIN * NONINV 1 999 2.4 0.3 2.4 0.3\n"
//                                 "GATE   xor2a   5  O=a*!b+!a*b;      PIN * UNKNOWN 2 999 1.9 0.5 1.9 0.5\n"
//                                 "#GATE  xor2b   5  O=!(a*b+!a*!b);   PIN * UNKNOWN 2 999 1.9 0.5 1.9 0.5\n"
//                                 "GATE   xnor2a  5  O=a*b+!a*!b;      PIN * UNKNOWN 2 999 2.1 0.5 2.1 0.5\n"
//                                 "#GATE  xnor2b  5  O=!(a*!b+!a*b);   PIN * UNKNOWN 2 999 2.1 0.5 2.1 0.5\n"
//                                 "GATE   aoi21   3  O=!(a*b+c);       PIN * INV 1 999 1.6 0.4 1.6 0.4\n"
//                                 "GATE   aoi22   4  O=!(a*b+c*d);     PIN * INV 1 999 2.0 0.4 2.0 0.4\n"
//                                 "GATE   oai21   3  O=!((a+b)*c);     PIN * INV 1 999 1.6 0.4 1.6 0.4\n"
//                                 "GATE   oai22   4  O=!((a+b)*(c+d)); PIN * INV 1 999 2.0 0.4 2.0 0.4\n"
//                                 "GATE   buf     2  O=a;              PIN * NONINV 1 999 1.0 0.0 1.0 0.0\n"
//                                 "GATE   zero    0  O=CONST0;\n"
//                                 "GATE   one     0  O=CONST1;";

int main()
{
  using namespace experiments;
  using namespace mockturtle;
  using namespace rils;

  experiment<std::string, double, double, double, double, float, float, bool> exp(
      "rig-mapper", "benchmark", "a(map)", "a(opt)", "d(map)", "d(opt)", "t(map)", "t(opt)", "cec" );

  fmt::print( "[i] processing technology library\n" );

  /* library to map to MIGs */
  mig_npn_resynthesis resyn{ true };
  exact_library_params eps;
  exact_library<mig_network, mig_npn_resynthesis> exact_lib( resyn, eps );

  /* library to map to technology */
  std::vector<gate> gates;
  std::istringstream in( rils::detail::mcnc_library );

  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  tech_library_params tps;
  tech_library<5, classification_type::np_configurations> tech_lib( gates, tps );

  for ( auto const& benchmark : all_benchmarks( iscas ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    const uint32_t size_before = aig.num_gates();
    const uint32_t depth_before = depth_view( aig ).depth();

    map_params ps2;
    ps2.cut_enumeration_ps.minimize_truth_table = true;
    ps2.cut_enumeration_ps.cut_limit = 24;
    map_stats st2;

    rig_network res2 = rils::map( aig, tech_lib, ps2, &st2 );

    if( res2.num_gates() > 50000 ) continue;

    double area_before = res2.compute_area();
    double delay_before = res2.compute_worst_delay();

    res2.report_binding_stats();
    res2.report_gates_usage();

    const auto cec2 = benchmark == "hyp" ? true : abc_cec( res2, benchmark );

    resubstitution_params rps;
    resubstitution_stats rst;
    // ps.pattern_filename = "1024sa1/" + benchmark + ".pat";
    rps.progress =true;
    rps.max_inserts = 20;
    rps.max_trials = 1;
    rps.max_pis = 10;
    rps.max_divisors = 20;//std::numeric_limits<uint32_t>::max();

    printf("spf %d\n", res2.num_gates() );  
    rig_resubstitution<rils::network_t::MAPPED, rils::support_selection_t::PIVOT, 4>( res2, rps, &rst );
    res2 = cleanup_dangling( res2 );

    printf("spf %d\n", res2.num_gates() );  
    res2.report_binding_stats();
    res2.report_gates_usage();
    printf("\n\n");

    exp( benchmark, area_before, res2.compute_area(), delay_before, res2.compute_worst_delay(), to_seconds( st2.time_total ), to_seconds( rst.time_total ), cec2 );


  }

  exp.save();
  exp.table();

  return 0;
}

//| benchmark |  a(map) |  a(opt) | d(map) | d(opt) | t(map) | t(opt) |  cec |
//|       c17 |   12.00 |   12.00 |   3.00 |   3.00 |   0.00 |   0.00 | true |
//|      c432 |  359.00 |  357.00 |  20.00 |  22.60 |   0.00 |   0.08 | true |
//|      c499 |  808.00 |  776.00 |  15.60 |  16.50 |   0.00 |   0.17 | true |
//|      c880 |  590.00 |  590.00 |  18.20 |  18.20 |   0.00 |   0.14 | true |
//|     c1355 |  776.00 |  744.00 |  15.60 |  16.50 |   0.01 |   0.15 | true |
//|     c1908 |  670.00 |  616.00 |  22.00 |  23.60 |   0.00 |   0.10 | true |
//|     c2670 | 1192.00 | 1191.00 |  16.10 |  16.10 |   0.01 |   0.26 | true |
//|     c3540 | 1802.00 | 1790.00 |  30.40 |  31.90 |   0.01 |   0.86 | true |
//|     c5315 | 3062.00 | 3027.00 |  31.20 |  32.50 |   0.02 |   1.08 | true |
//|     c6288 | 4562.00 | 4311.00 |  80.00 |  82.80 |   0.04 |   0.27 | true |
//|     c7552 | 2894.00 | 2892.00 |  22.10 |  23.40 |   0.02 |   1.36 | true |