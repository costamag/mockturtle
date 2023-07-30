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
#include <mockturtle/algorithms/balancing/sop_balancing.hpp>
#include <mockturtle/algorithms/balancing/mct_balancing.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <string>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, double, bool, uint32_t, uint32_t, double, bool> exp( "mct_balancing", "benchmark", "size", "depth", "size 4", "depth 4", "RT 4", "cec 4", "size 6", "depth 6", "RT 6", "cec 6" );

  mct_rebalancing<xag_network> mct_balancing;
  sop_rebalancing<xag_network> sop_balancing;

  for ( auto const& benchmark : iscas_benchmarks())//epfl_benchmarks( ~experiments::hyp ) )//
  {
    fmt::print( "[i] processing {}\n", benchmark );
    xag_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    depth_view daig{ aig };

    balancing_params ps;
    balancing_stats st4_mct, st4_sop;
    double time_sop=0;
    ps.progress = true;
    ps.only_on_critical_path = true;
    ps.cut_enumeration_ps.cut_size = 4u;
    const auto aig4_sop = balancing( aig, { sop_balancing }, ps, &st4_sop );
    depth_view daig_sop{ aig4_sop };

    ps.cut_enumeration_ps.cut_size = 4u;
    const auto aig4_mct = balancing( aig, { mct_balancing }, ps, &st4_mct );
    depth_view daig_mct{ aig4_mct };

    printf("SOP: d=%d g=%d\n",daig_sop.depth(), aig4_sop.num_gates() );
    printf("SYM: d=%d g=%d\n",daig_mct.depth(), aig4_mct.num_gates() );

    const auto cec4_mct = abc_cec( aig4_mct, benchmark );
    const auto cec4_sop = abc_cec( aig4_sop, benchmark );

    exp( benchmark,
         aig.num_gates(), daig.depth(),
         aig4_sop.num_gates(), daig_sop.depth(),
         to_seconds( st4_sop.time_total ), cec4_sop,
         aig4_mct.num_gates(), daig_mct.depth(),
         to_seconds( st4_mct.time_total ), cec4_mct );
  }

  exp.save();
  exp.table();

  return 0;
}