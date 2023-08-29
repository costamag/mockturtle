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
#include <mockturtle/algorithms/balancing/sym_balancing.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <string>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, double, bool, uint32_t, uint32_t, double, bool> exp( "sym_balancing", "benchmark", "size", "depth", "size 4", "depth 4", "RT 4", "cec 4", "size 6", "depth 6", "RT 6", "cec 6" );

  sym_rebalancing<xag_network> sym_balancing;

  for ( auto const& benchmark : iscas_benchmarks( ) )//~experiments::hyp
  {
    fmt::print( "[i] processing {}\n", benchmark );
    xag_network xaig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xaig ) ) != lorina::return_code::success )
    {
      continue;
    }

    balancing_params ps;
    balancing_stats st4, st6;

    ps.progress = true;
    ps.cut_enumeration_ps.cut_size = 4u;
    const auto xaig4 = balancing( xaig, { sym_balancing }, ps, &st4 );

    ps.cut_enumeration_ps.cut_size = 6u;
    const auto xaig6 = balancing( xaig, { sym_balancing }, ps, &st6 );

    depth_view dxaig{ xaig };
    depth_view dxaig4{ xaig4 };
    depth_view dxaig6{ xaig6 };

    const auto cec4 = abc_cec( xaig4, benchmark );
    const auto cec6 = abc_cec( xaig6, benchmark );

    exp( benchmark,
         xaig.num_gates(), dxaig.depth(),
         xaig4.num_gates(), dxaig4.depth(),
         to_seconds( st4.time_total ), cec4,
         xaig6.num_gates(), dxaig6.depth(),
         to_seconds( st6.time_total ), cec6 );
  }

  exp.save();
  exp.table();

  return 0;
}