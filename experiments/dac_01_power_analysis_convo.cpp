/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2024  EPFL
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

#include <iostream>
#include <string>
#include <vector>
#include <chrono>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/genlib.hpp>
#include <mockturtle/algorithms/aig_balancing.hpp>
#include <mockturtle/algorithms/emap.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/block.hpp>
#include <mockturtle/utils/name_utils.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/cell_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/names_view.hpp>

#include <mockturtle/resynthesis/networks/bound.hpp>
#include <mockturtle/resynthesis/networks/mapped_to_bound.hpp>
#include <mockturtle/resynthesis/power/analyzers/convo_analyzer.hpp>
#include <mockturtle/resynthesis/power/node_simulator.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, double, double, double, double, double, float, bool> exp(
      "dac_001", "benchmark", "size", "area", "delay", "zero-delay", "glitching", "switching", "runtime", "cec" );

  /* library to map to technology */
  //fmt::print( "[i] processing technology library\n" );
  std::string library = "multioutput";
  std::vector<gate> gates;
  std::ifstream in( cell_libraries_path( library ) );

  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  tech_library_params tps;
  tps.ignore_symmetries = false; // set to true to drastically speed-up mapping with minor delay increase
  tps.verbose = true;
  tech_library<9> tech_lib( gates, tps );

  for ( auto const& benchmark : all_benchmarks( iscas ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    names_view<aig_network> aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    /* remove structural redundancies */
    aig_balancing_params bps;
    bps.minimize_levels = false;
    bps.fast_mode = true;
    aig_balance( aig, bps );

    const uint32_t size_before = aig.num_gates();
//    const uint32_t depth_before = depth_view( aig ).depth();

    emap_params ps;
    ps.matching_mode = emap_params::hybrid;
    ps.area_oriented_mapping = false;
    ps.map_multioutput = false;
    ps.relax_required = 0;
    emap_stats st;
    cell_view<block_network> res = emap<9>( aig, tech_lib, ps, &st );

    auto bound = convert_mapped_to_bound( res, gates );

    uint32_t constexpr N = 10;
    uint32_t constexpr I = 10;
    using Ntk = bound_network;
    using node_simulator_t = mockturtle::node_simulator<Ntk, N, I>;
    
    using dNtk = depth_view<bound_network>;
    dNtk dbound{bound};
    mockturtle::convo_analyzer<dNtk, N, I> power( dbound );

    auto beg = std::chrono::high_resolution_clock::now();
    power( 1 );
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - beg);

    double const zero_delay = power.get_zero_delay();
    double const dyn_switch = power.get_dyn_switch();
    double const glitching = dyn_switch - zero_delay;

    power.print_ntk();

    std::cout << "ZD: train " << zero_delay << " " << zero_delay << std::endl;
    std::cout << "DS: train " << dyn_switch << " " << dyn_switch << std::endl;

    const auto cec = false;//benchmark == "hyp" ? true : abc_cec_mapped_cell( res_names, benchmark, library );

    /* write verilog netlist */
    // write_verilog_with_cell( res_names, benchmark + "_mapped.v" );
    exp( benchmark, size_before, bound.area(), res.compute_worst_delay(), zero_delay, glitching, dyn_switch, duration.count(), cec );
  }

  exp.save();
  exp.table();

  return 0;
}