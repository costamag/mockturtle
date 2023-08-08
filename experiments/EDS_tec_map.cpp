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

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <mockturtle/algorithms/balancing.hpp>
#include <mockturtle/algorithms/balancing/mct1_balancing.hpp>
#include <mockturtle/algorithms/balancing/mct1_balancing.hpp>
#include <mockturtle/algorithms/sim_resub.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/algorithms/collapse_mapped.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/mapping_view.hpp>

#include <mockturtle/algorithms/mapper.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/views/depth_view.hpp>


#include <string>

#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>

#include <experiments.hpp>

std::string const mcnc_library = "GATE   inv1    1  O=!a;             PIN * INV 1 999 0.9 0.3 0.9 0.3\n"
                                 "GATE   inv2    2  O=!a;             PIN * INV 2 999 1.0 0.1 1.0 0.1\n"
                                 "GATE   inv3    3  O=!a;             PIN * INV 3 999 1.1 0.09 1.1 0.09\n"
                                 "GATE   inv4    4  O=!a;             PIN * INV 4 999 1.2 0.07 1.2 0.07\n"
                                 "GATE   nand2   2  O=!(a*b);         PIN * INV 1 999 1.0 0.2 1.0 0.2\n"
                                 "GATE   nand3   3  O=!(a*b*c);       PIN * INV 1 999 1.1 0.3 1.1 0.3\n"
                                 "GATE   nand4   4  O=!(a*b*c*d);     PIN * INV 1 999 1.4 0.4 1.4 0.4\n"
                                 "GATE   nor2    2  O=!(a+b);         PIN * INV 1 999 1.4 0.5 1.4 0.5\n"
                                 "GATE   nor3    3  O=!(a+b+c);       PIN * INV 1 999 2.4 0.7 2.4 0.7\n"
                                 "GATE   nor4    4  O=!(a+b+c+d);     PIN * INV 1 999 3.8 1.0 3.8 1.0\n"
                                 "GATE   and2    3  O=a*b;            PIN * NONINV 1 999 1.9 0.3 1.9 0.3\n"
                                 "GATE   or2     3  O=a+b;            PIN * NONINV 1 999 2.4 0.3 2.4 0.3\n"
                                 "GATE   xor2a   5  O=a*!b+!a*b;      PIN * UNKNOWN 2 999 1.9 0.5 1.9 0.5\n"
                                 "#GATE  xor2b   5  O=!(a*b+!a*!b);   PIN * UNKNOWN 2 999 1.9 0.5 1.9 0.5\n"
                                 "GATE   xnor2a  5  O=a*b+!a*!b;      PIN * UNKNOWN 2 999 2.1 0.5 2.1 0.5\n"
                                 "#GATE  xnor2b  5  O=!(a*!b+!a*b);   PIN * UNKNOWN 2 999 2.1 0.5 2.1 0.5\n"
                                 "GATE   aoi21   3  O=!(a*b+c);       PIN * INV 1 999 1.6 0.4 1.6 0.4\n"
                                 "GATE   aoi22   4  O=!(a*b+c*d);     PIN * INV 1 999 2.0 0.4 2.0 0.4\n"
                                 "GATE   oai21   3  O=!((a+b)*c);     PIN * INV 1 999 1.6 0.4 1.6 0.4\n"
                                 "GATE   oai22   4  O=!((a+b)*(c+d)); PIN * INV 1 999 2.0 0.4 2.0 0.4\n"
                                 "GATE   buf     2  O=a;              PIN * NONINV 1 999 1.0 0.0 1.0 0.0\n"
                                 "GATE   zero    0  O=CONST0;\n"
                                 "GATE   one     0  O=CONST1;";

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, double, double, uint32_t, uint32_t, double, double, float, float, bool, bool> exp(
      "mcts", "benchmark", "size", "size_after", "area", "area_after", "depth", "depth_after", "delay", "delay_after", "runtime1", "runtime2", "equivalent1", "equivalent2" );
  mcts_rebalancing<xag_network> mct_balancing;

  /* library to map to technology */
  std::vector<gate> gates;
  std::istringstream in( mcnc_library );

  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  tech_library_params tps;
  tech_library<5, classification_type::np_configurations> tech_lib( gates, tps );


  for ( auto const& benchmark : epfl_benchmarks( ~experiments::hyp  ) )//bms__ bms__
  {
    fmt::print( "[i] processing {}\n", benchmark );
    xag_network xag;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xag ) ) != lorina::return_code::success )
    {
      continue;
    }


    const uint32_t size_before = xag.num_gates();
    const uint32_t depth_before = depth_view( xag ).depth();

    map_params ps0;
    ps0.skip_delay_round = false;
    ps0.required_time = std::numeric_limits<double>::max();
    ps0.cut_enumeration_ps.minimize_truth_table = true;
    ps0.cut_enumeration_ps.cut_limit = 24;

    map_stats st0;
    binding_view<klut_network> res0 = map( xag, tech_lib, ps0, &st0 );
    const auto cec0 = benchmark == "hyp" ? true : abc_cec( res0, benchmark );

    using namespace std::chrono;

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    uint32_t DEPTH = std::numeric_limits<uint32_t>::max();
    uint32_t SIZE = std::numeric_limits<uint32_t>::max();

    depth_view dxag{ xag };

    balancing_params ps;
    balancing_stats st;
    double time_sop=0;
    ps.progress = true;
    ps.only_on_critical_path = true;
    ps.cut_enumeration_ps.cut_size = 4u;
    auto xag_opt = balancing( xag, { mct_balancing }, ps, &st );
    auto xag_bst = xag;
    depth_view dxag_opt{ xag_opt };

    if( (dxag_opt.depth() < DEPTH) || ((dxag_opt.depth() == DEPTH) && (dxag_opt.num_gates() < SIZE)) )
    {
      DEPTH = dxag_opt.depth();
      SIZE = dxag_opt.num_gates();
      xag_bst = xag_opt;
    }

    auto depth_old = dxag_opt.depth()+1;
    auto depth_new = dxag_opt.depth();
    
    uint32_t K{0u};
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

    bool isUpdated0{false};
    bool isUpdated1{false};
    bool isUpdated2{false};
    bool isUpdated3{false};
    bool isUpdated4{false};

    int IT = 0;

    while( (time_span.count() < 600)&&( time_span.count() < 120 || IT < 4 || (isUpdated0||isUpdated1||isUpdated2||isUpdated3||isUpdated4)) )
    {
      IT++;
      if(depth_old == depth_new && K < 3u )
        ps.cut_enumeration_ps.cut_size = 4u+(K++);
      else if( depth_old == depth_new && K >=3u )
      {
        ps.only_on_critical_path = false;
        K=0;
        ps.cut_enumeration_ps.cut_size = 4u;
      }
      else
      {
        ps.only_on_critical_path = true;
        K=0;   
      }

      xag_opt = balancing( xag_opt, { mct_balancing }, ps, &st );

      depth_view dloc{ xag_opt };
      //daig_mct = d_mct;
      printf("SOPi: d=%d/%d g=%d/%d\n",dloc.depth(), dxag.depth(), dloc.num_gates(), dxag.num_gates() );
      depth_old = depth_new;
      depth_new = dloc.depth();
      dxag_opt = dloc;

      if( (dxag_opt.depth() < DEPTH) || ((dxag_opt.depth() == DEPTH) && (dxag_opt.num_gates() < SIZE)) )
      {
        DEPTH = dxag_opt.depth();
        SIZE = dxag_opt.num_gates();
        xag_bst = xag_opt;
      }

      t2 = high_resolution_clock::now();
      time_span = duration_cast<duration<double>>(t2 - t1);
      ps.only_on_critical_path = true;

      isUpdated0 = isUpdated1;
      isUpdated1 = isUpdated2;
      isUpdated2 = isUpdated3;
      isUpdated3 = isUpdated4;
      isUpdated4 = (depth_old > depth_new);
    }

    resubstitution_params res_ps;
    resubstitution_stats res_st;

    // ps.pattern_filename = "1024sa1/" + benchmark + ".pat";

    depth_view dxag_bst{ xag_bst };

    map_params ps1;
    ps1.skip_delay_round = false;
    ps1.required_time = std::numeric_limits<double>::max();
    ps1.cut_enumeration_ps.minimize_truth_table = true;
    ps1.cut_enumeration_ps.cut_limit = 24;
    map_stats st1;

    binding_view<klut_network> res1 = map( xag_bst, tech_lib, ps1, &st1 );
    const auto cec1 = benchmark == "hyp" ? true : abc_cec( res1, benchmark );


    printf("-->: d=%d/%d g=%d/%d\n", dxag_bst.depth(), dxag.depth(), dxag_bst.num_gates(), dxag.num_gates() );
    printf("%d : [NONE: d=%f g=%f] [MCTS: d=%f g=%f]\n", cec1, st0.delay, st0.area, st1.delay, st1.area );


    t2 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t2 - t1);

    const auto cec = abc_cec( xag, benchmark );
    const auto cec_opt = abc_cec( xag_opt, benchmark );

    exp( benchmark, dxag.num_gates(), dxag_bst.num_gates(), st0.area, st1.area, dxag.depth(), dxag_bst.depth(), st0.delay, st1.delay, to_seconds( st0.time_total ), to_seconds( st1.time_total ), cec0, cec1 );
  
  
  }

  exp.save();
  exp.table();

  return 0;
}