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
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/rewrite.hpp>
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

#include <mockturtle/resynthesis/aig_syn.hpp>
#include <mockturtle/resynthesis/algorithms/resynthesize.hpp>
//#include <mockturtle/resynthesis/analyzers/homo_xpx_analyzer.hpp>
#include <mockturtle/resynthesis/power/node_simulator.hpp>
#include <mockturtle/resynthesis/power/analyzers/power_analyzer.hpp>
#include <mockturtle/resynthesis/database/homo_database.hpp>
#include <mockturtle/resynthesis/database/mapped_list.hpp>
#include <mockturtle/resynthesis/list_simulator.hpp>
#include <mockturtle/resynthesis/networks/bindings.hpp>
#include <mockturtle/resynthesis/networks/bound.hpp>
#include <mockturtle/resynthesis/database/homo_storage.hpp>
#include <mockturtle/resynthesis/networks/mapped_to_bound.hpp>
#include <mockturtle/resynthesis/networks/cleanup_bound.hpp>
#include <mockturtle/resynthesis/power/analyzers/convo_analyzer.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, double, double, double, double, double, double, double, double, float, bool> exp(
      "dac_001", "benchmark", "A(bef)", "A(aft)", "D(bef)", "D(aft)", "G(bef)", "G(aft)", "P(bef)", "P(aft)", "runtime", "cec" );

  /* library to map to technology */
  //fmt::print( "[i] processing technology library\n" );
  std::string library = "asap7";
  std::vector<gate> gates;
  std::ifstream in( cell_libraries_path( library ) );

  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  bindings binds( gates );
  /* get the database */
  static constexpr uint32_t num_vars = 4;
  static constexpr uint32_t capacity = 10;
  static constexpr mockturtle::choice_order order = power_area_delay;
  using database_t = mockturtle::homo_database<num_vars, capacity, order>;
  /*
  _base_
  _emap_area_
  _lazy_dap_
  */
  std::string path = "/home/costamag/POST_MAPPING/mockturtle/build/" + library + "_lazy_pad_" + std::to_string( num_vars ) + ".db";

  database_t database( binds, path );// = create_resyn_database<num_vars, capacity, order>();
//database.print();
  tech_library_params tps;
  tps.ignore_symmetries = false; // set to true to drastically speed-up mapping with minor delay increase
  tps.verbose = true;
  tech_library<9> tech_lib( gates, tps );

  /* run resynthesis */
  for ( auto const& benchmark : all_benchmarks( iscas | epfl ) )
  {
    if( benchmark == "hyp" || benchmark == "div" || benchmark == "mem_ctrl" )
      continue;
    fmt::print( "[i] processing {}\n", benchmark );

    names_view<aig_network> aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    mockturtle::xag_npn_resynthesis<aig_network, aig_network, mockturtle::xag_npn_db_kind::aig_complete> resyn;
    mockturtle::exact_library_params eps;
    eps.np_classification = false;
    mockturtle::exact_library<aig_network> exact_lib( resyn, eps );

//    mockturtle::rewrite_params rwps;
//    rwps.preserve_depth = true;
//    std::cout << aig.num_gates() << " " << aig.num_pos() << std::endl;
//    mockturtle::rewrite( aig, exact_lib, rwps );
//    mockturtle::rewrite( aig, exact_lib, rwps );
//    std::cout << aig.num_gates() << " " << aig.num_pos() << std::endl;
//    /* remove structural redundancies */
//    aig_balancing_params bps;
//    bps.minimize_levels = false;
//    bps.fast_mode = true;
//    aig_balance( aig, bps );
//    const uint32_t size_before = aig.num_gates();
//    const uint32_t depth_before = depth_view( aig ).depth();

    emap_params ps;
    ps.matching_mode = emap_params::hybrid;
    ps.area_oriented_mapping = true;
    ps.eswp_rounds = 4;
    ps.map_multioutput = false;
    emap_stats st;
    cell_view<block_network> res = emap<9>( aig, tech_lib, ps, &st );

    auto bound = convert_mapped_to_bound( res, gates );

    uint32_t constexpr num_steps = 10;
    uint32_t constexpr cube_size = 10;
    uint32_t constexpr max_leaves = 10;

    using Ntk = bound_network;
    
    using dNtk = depth_view<bound_network>;
    dNtk dbound{bound};

    using node_simulator_t = mockturtle::node_simulator<Ntk, num_steps, cube_size>;

    std::vector<node_simulator_t> pis_sims;
    std::default_random_engine::result_type seed = 1;
    std::array<kitty::static_truth_table<cube_size>, 2> workload;
    int i = 0;
    bound.foreach_pi([&]( auto n ){
      kitty::create_random( workload[0], seed + i );
      kitty::create_random( workload[1], seed + i + 1 );
      i += 2;
      node_simulator_t sim( n, workload , -0.01, -0.02, 0 );
      pis_sims.push_back( sim );
    });

    double power_simple_1_zero_delay;
    double power_simple_1_dyn_before;
    double power_simple_1_glitching ;
    double delay_bef ;
    double area_bef       ;
    double zero_delay_bef ;
    double dyn_switch_bef ;
    double glitching_bef  ;

    {
      mockturtle::power_analyzer<Ntk, num_steps, cube_size> power_simple_1( bound );
      power_simple_1( pis_sims );
      power_simple_1_zero_delay = power_simple_1.get_zero_delay();
      power_simple_1_dyn_before = power_simple_1.get_dyn_switch();
      power_simple_1_glitching = power_simple_1_dyn_before - power_simple_1_zero_delay;
      arrival_analyzer<Ntk> arrival1( bound );
      arrival1.compute_arrival_times();
      delay_bef = arrival1.worst_delay();
      area_bef = bound.area();
    }
    std::cout << "ZD(power_simple_1_bef) " << power_simple_1_zero_delay << std::endl;
    std::cout << "DS(power_simple_1_bef) " << power_simple_1_dyn_before << std::endl;

    auto beg = std::chrono::high_resolution_clock::now();
    mockturtle::resynthesis_params rps;
    rps.preserve_depth = true;
    rps.dynamic_database = false;
    rps.try_struct = true;
    rps.try_window = true;
    rps.try_simula = true;

    mockturtle::power_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves>( bound, database, rps );

    auto end = std::chrono::high_resolution_clock::now();

    dNtk dbound2{bound};

    double delay_aft ;
    double area_aft ;
    double power_simple_2_zero_delay ;
    double power_simple_2_dyn_before ;
    double power_simple_2_glitching ;

    auto result = mockturtle::cleanup_bound( bound, gates );
    {
      arrival_analyzer<Ntk> arrival2( result );
      arrival2.compute_arrival_times();
      delay_aft = arrival2.worst_delay();
      area_aft = result.area();

      mockturtle::power_analyzer<Ntk, num_steps, cube_size> power_simple_2( result );
      power_simple_2( pis_sims );
      power_simple_2_zero_delay = power_simple_2.get_zero_delay();
      power_simple_2_dyn_before = power_simple_2.get_dyn_switch();
      power_simple_2_glitching = power_simple_2_dyn_before - power_simple_2_zero_delay;
    }
    std::cout << "ZD(power_simple_2_bef) " << power_simple_2_zero_delay << std::endl;
    std::cout << "DS(power_simple_2_bef) " << power_simple_2_dyn_before << std::endl;

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - beg);

    const auto cec = benchmark == "hyp" ? true : abc_cec( bound, benchmark );

    /* write verilog netlist */
    // write_verilog_with_cell( res_names, benchmark + "_mapped.v" );
    exp( benchmark, area_bef, area_aft, delay_bef, delay_aft, power_simple_1_glitching, power_simple_2_glitching, power_simple_1_dyn_before, power_simple_2_dyn_before, duration.count(), cec );
  }

  exp.save();
  exp.table();

  return 0;
}