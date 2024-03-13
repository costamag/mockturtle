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
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/algorithms/lut_mapper2.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/lig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/mapping_view.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;


  for ( auto const& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    lut_map_params ps;
    ps.cut_enumeration_ps.cut_size = 6u;
    ps.cut_enumeration_ps.cut_limit = 8u;
    ps.recompute_cuts = true;
    ps.area_oriented_mapping = false;
    ps.cut_expansion = true;
    lut_map_stats st;
    const auto klut = lut_map( aig, ps, &st );

    depth_view<klut_network> klut_d{ klut };


    scopt::lut_map2_params ps2;
    ps2.cut_enumeration_ps.cut_size = 6u;
    ps2.cut_enumeration_ps.cut_limit = 8u;
    ps2.recompute_cuts = true;
    ps2.area_oriented_mapping = false;
    ps2.cut_expansion = true;
    scopt::lut_map2_stats st2;
    const auto lig = scopt::lut_map2( aig, ps2, &st2 );

    depth_view<lig_network> lig_d{ lig };

    printf("%d %d\n", klut.num_gates(), lig.num_gates());
    printf("%d %d\n", klut_d.depth(), lig_d.depth());


    auto const cec = benchmark == "hyp" ? true : abc_cec( klut, benchmark );
    auto const cec2 = benchmark == "hyp" ? true : abc_cec( lig, benchmark );
    if( !cec )  printf("[e] klut not equivalent\n");
    if( !cec2 )  printf("[e] lig not equivalent\n");

  }

  return 0;
}
