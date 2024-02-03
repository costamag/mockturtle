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
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/rig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/algorithms/collapse_mapped.hpp>
#include <mockturtle/algorithms/lut_mapping.hpp>
#include <mockturtle/views/mapping_view.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;
  using namespace rils;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, double, double, bool, bool> exp( "rig_exp_0", "benchmark", "g(LUT)", "g(RIG)", "d(LUT)", "d(RIG)", "t(LUT)", "t(RIG)", "eq(LUT)", "eq(RIG)" );

  for ( auto const& benchmark : iscas_benchmarks( ) ) //
  {
    fmt::print( "[i] processing {}\n", benchmark );
    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      printf("aig unsuccessful\n");
      continue;
    }

    std::string tmp = benchmark + "tmp.blif";
    write_blif( aig, tmp );

    klut_network klut;
    if ( lorina::read_blif( tmp, blif_reader( klut ) ) != lorina::return_code::success )
    {
      printf("klut unsuccessful\n");
      continue;
    }
    depth_view<klut_network> klut_d{ klut };
    auto const lut_cec = benchmark == "hyp" ? true : abc_cec( klut, benchmark );

    rig_network rig;
    rig._is_smart = true;
    if ( lorina::read_blif( tmp, blif_reader( rig ) ) != lorina::return_code::success )
    {
      printf("rig unsuccessful\n");
      continue;
    }

    depth_view<rig_network> rig_d{ rig };
    auto const rig_cec = benchmark == "hyp" ? true : abc_cec( rig, benchmark );

    exp( benchmark, klut.num_gates(), rig.num_gates(), klut_d.depth(), rig_d.depth(), 0, 0, lut_cec, rig_cec );
  }

  exp.save();
  exp.table();

  return 0;
}
