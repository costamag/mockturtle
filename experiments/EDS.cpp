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
#include <string>

#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, bool, uint32_t, uint32_t, double, bool> exp( "eds", "benchmark", "s(ORI)", "d(ORI)", "c(ORI)", "s(MCT)", "d(MCT)", "t(MCT)", "c(MCT)" );

  mcts_rebalancing<xag_network> mct_balancing;


  for ( auto const& benchmark : iscas_benchmarks( c6288 ) )//bms__ 
  {
    fmt::print( "[i] processing {}\n", benchmark );
    xag_network xag;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xag ) ) != lorina::return_code::success )
    {
      continue;
    }

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
    depth_view dxag_opt{ xag_opt };

    if( (dxag_opt.depth() < DEPTH) || ((dxag_opt.depth() == DEPTH) && (dxag_opt.num_gates() < SIZE)) )
    {
      DEPTH = dxag_opt.depth();
      SIZE = dxag_opt.num_gates();
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

    while( time_span.count() < 180 || (isUpdated0||isUpdated1||isUpdated2||isUpdated3||isUpdated4) )
    {
      if(depth_old == depth_new && K < 3u)
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


    printf("-->: d=%d/%d g=%d/%d\n",DEPTH, dxag.depth(), SIZE, dxag.num_gates() );


    t2 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t2 - t1);

    const auto cec = abc_cec( xag, benchmark );
    const auto cec_opt = abc_cec( xag_opt, benchmark );

    exp( benchmark, xag.num_gates(), dxag.depth(), cec, SIZE, DEPTH, time_span.count(), cec_opt  );
  }

  exp.save();
  exp.table();

  return 0;
}