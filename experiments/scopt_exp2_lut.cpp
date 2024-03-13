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

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <mockturtle/algorithms/lut_mapper2.hpp>
#include <mockturtle/algorithms/post_mapping/lut_resub.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/lig.hpp>
#include <mockturtle/views/mapping_view.hpp>
#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, float, bool> exp( "lutopt", "benchmark", "a(map)", "a(opt)", "d(map)", "d(opt)", "runtime", "equivalent" );

  for ( auto const& benchmark : epfl_benchmarks( experiments::priority ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    scopt::lut_map2_params ps2;
    ps2.cut_enumeration_ps.cut_size = 4u;
    ps2.cut_enumeration_ps.cut_limit = 8u;
    ps2.recompute_cuts = true;
    ps2.area_oriented_mapping = true;
    ps2.cut_expansion = true;
    scopt::lut_map2_stats st2;
    lig_network lig = scopt::lut_map2( aig, ps2, &st2 );
    depth_view<lig_network> lig_d{ lig };

    uint32_t const initial_depth = lig_d.depth();
    uint32_t const initial_size = lig.num_gates();

    lut_resub_params ps;
    //ps.cut_enumeration_ps.cut_size = 6;
    //ps.cut_enumeration_ps.cut_limit = 16;
    //ps.conflict_limit = 100;
    //ps.progress = true;
    lut_resub_stats st;

    lut_resub( lig, ps, &st );

    auto cec = benchmark == "hyp" ? true : abc_cec( lig, benchmark );

    exp( benchmark, initial_size, lig.num_gates(), initial_depth, lig_d.depth(), to_seconds( st.time_total ), cec );
  }

  exp.save();
  exp.table();

  return 0;
}