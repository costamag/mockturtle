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
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/rewrite.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/block.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/resynthesis/aig_syn.hpp>
#include <mockturtle/resynthesis/database/homo_database.hpp>
#include <mockturtle/resynthesis/database/mapped_list.hpp>
#include <mockturtle/resynthesis/list_simulator.hpp>
#include <mockturtle/resynthesis/networks/bindings.hpp>
#include <mockturtle/resynthesis/networks/bound.hpp>
#include <mockturtle/resynthesis/networks/cleanup_bound.hpp>
#include <mockturtle/resynthesis/networks/mapped_to_bound.hpp>
#include <mockturtle/resynthesis/power/optimization/power_rewrite.hpp>
#include <mockturtle/resynthesis/power/power_analyzer.hpp>
#include <mockturtle/resynthesis/power/node_simulator.hpp>
#include <mockturtle/resynthesis/database/homo_database.hpp>
#include <mockturtle/utils/name_utils.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/cell_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/names_view.hpp>

#include <experiments.hpp>

/*! \brief construc a database from a standard cell library */
template<uint32_t num_vars>
mockturtle::homo_database<num_vars> create_database( std::string const& );

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, double, double, double, double, float, bool> exp(
      "dac_001", "benchmark", "size", "area", "delay", "SW(bef)", "SW(after)", "runtime", "cec" );

  /* library to map to technology */
  fmt::print( "[i] processing technology library\n" );
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

  for ( auto const& benchmark : all_benchmarks( experiments::c17 ) )
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

    uint32_t constexpr num_vars = 4;
    uint32_t constexpr num_steps = 30;
    uint32_t constexpr cube_size = 10;
    using Ntk = bound_network;
    using node_simulator_t = mockturtle::node_simulator<Ntk, num_steps, cube_size>;
    
    mockturtle::power_analyzer<Ntk, num_steps, cube_size> power( bound );

    std::vector<node_simulator_t> pis_sims;
    std::default_random_engine::result_type seed = 1;
    std::array<kitty::static_truth_table<cube_size>, 2> workload;
    int i = 0;
    bound.foreach_pi([&]( auto n ){
      kitty::create_random( workload[0], seed + i );
      kitty::create_random( workload[1], seed + i + 1 );
      i += 2;
      node_simulator_t sim( n, workload , -0.5, -1, 0 );
      pis_sims.push_back( sim );
    });

    auto beg = std::chrono::high_resolution_clock::now();
    power( pis_sims );
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - beg);

    double const zero_delay = power.get_zero_delay();
    double const dyn_before = power.get_dyn_switch();
    double const glitching = dyn_before - zero_delay;
    std::cout << "ZD(bef) " << zero_delay << std::endl;
    std::cout << "DS(bef) " << dyn_before << std::endl;

    bound.print();
    /* write verilog netlist */
    // write_verilog_with_cell( res_names, benchmark + "_mapped.v" );
    std::string const library{ "asap7" };
    auto database = create_database<num_vars>( library );
    mockturtle::power_rewrite<Ntk, num_vars, num_steps, cube_size>( bound, pis_sims, database );

    mockturtle::bound_network bound2 = mockturtle::cleanup_bound( bound, gates );

    bound.print();

    mockturtle::power_analyzer<Ntk, num_steps, cube_size> power2( bound );
    power2( pis_sims );
    double const dyn_after = power2.get_dyn_switch();

    double const zero_delay2 = power2.get_zero_delay();
    double const dyn_before2 = power2.get_dyn_switch();
    double const glitching2 = dyn_before2 - zero_delay2;
    std::cout << "ZD(aft) " << zero_delay2 << std::endl;
    std::cout << "DS(aft) " << dyn_before2 << std::endl;
    const auto cec = benchmark == "hyp" ? true : abc_cec( bound, benchmark );

    exp( benchmark, size_before, bound.area(), res.compute_worst_delay(), dyn_before, dyn_after, duration.count(), cec );
  }

  exp.save();
  exp.table();

  return 0;
}

template<uint32_t num_vars = 4>
mockturtle::homo_database<num_vars> create_database( std::string const& library )
{
  using namespace experiments;
  using namespace mockturtle;

  std::vector<mockturtle::gate> gates;

  std::ifstream in( "/home/costamag/POST_MAPPING/mockturtle/experiments/cell_libraries/" + library + ".genlib" );//cell_libraries_path( library ));
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  bindings binds( gates );
  mockturtle::homo_database<num_vars> database( binds );

  if( result != lorina::return_code::success )
  {
    std::cerr << "[e] error in the library parsing" << std::endl;
    return database;
  }

  using Tt_t = kitty::dynamic_truth_table;
  using tt_hash = kitty::hash<Tt_t>;
  std::unordered_set<Tt_t, tt_hash> classes;
  Tt_t tt(4);
  do {
    const auto res = kitty::exact_p_canonization( tt );
    classes.insert( std::get<0>( res ) );
    kitty::next_inplace( tt );
  } while( !kitty::is_const0( tt ) );

  klut_network klut;
  std::vector<typename klut_network::signal> pis;
  for( int i = 0; i < 4; ++i )
  {
    pis.push_back( klut.create_pi() );
  }

  for( auto const& entry : classes )
  {
    auto const f = klut.create_node( pis, entry );
    klut.create_po( f );
  }
  aig_network aig = convert_klut_to_graph<aig_network, klut_network>( klut );

  xag_npn_resynthesis<aig_network> resyn;
  exact_library_params eps;
  eps.np_classification = false;
  exact_library<aig_network> exact_lib( resyn, eps );

  rewrite_params ps;
  ps.preserve_depth = false;
  rewrite( aig, exact_lib, ps );
  rewrite( aig, exact_lib, ps );
  rewrite( aig, exact_lib, ps );

  tech_library_params tps;
  tps.ignore_symmetries = false; // set to true to drastically speed-up mapping with minor delay increase
  tps.verbose = true;
  tech_library<9> tech_lib( gates, tps );
  emap_params mps;
  mps.matching_mode = emap_params::hybrid;
  mps.area_oriented_mapping = true;
  mps.map_multioutput = false;
  mps.relax_required = 0;
  emap_stats mst;
  cell_view<block_network> mapped = emap<9>( aig, tech_lib, mps, &mst );
  auto bound = convert_mapped_to_bound( mapped, gates );
  std::cout << bound.num_gates() << std::endl;

  bound.foreach_po( [&]( auto fo ){
    database.add_entry( bound, pis, fo );
  } );
  return database;
}