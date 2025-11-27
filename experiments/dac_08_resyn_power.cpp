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
#include <mockturtle/resynthesis/analyzers/homo_xpx_analyzer.hpp>
#include <mockturtle/resynthesis/power/node_simulator.hpp>
#include <mockturtle/resynthesis/power/analyzers/power_analyzer.hpp>
#include <mockturtle/resynthesis/database/homo_database.hpp>
#include <mockturtle/resynthesis/database/mapped_list.hpp>
#include <mockturtle/resynthesis/list_simulator.hpp>
#include <mockturtle/resynthesis/networks/bindings.hpp>
#include <mockturtle/resynthesis/networks/bound.hpp>
#include <mockturtle/resynthesis/networks/mapped_to_bound.hpp>
#include <mockturtle/resynthesis/power/analyzers/convo_analyzer.hpp>

#include <experiments.hpp>

std::string const test_library = "GATE   inv1    1 O=!a;            PIN * INV 1 999 0.9 0.3 0.9 0.3\n"
                                 "GATE   inv2    2 O=!a;            PIN * INV 2 999 1.0 0.1 1.0 0.1\n"
                                 "GATE   nand2   2 O=!(a*b);        PIN * INV 1 999 1.0 0.2 1.0 0.2\n"
                                 "GATE   and2    3 O=a*b;           PIN * INV 1 999 1.7 0.2 1.7 0.2\n"
                                 "GATE   xor2    4 O=a^b;           PIN * UNKNOWN 2 999 1.9 0.5 1.9 0.5\n"
                                 "GATE   mig3    3 O=a*b+a*c+b*c;   PIN * INV 1 999 2.0 0.2 2.0 0.2\n"
                                 "GATE   xor3    5 O=a^b^c;         PIN * UNKNOWN 2 999 3.0 0.5 3.0 0.5\n"
                                 "GATE   buf     2 O=a;             PIN * NONINV 1 999 1.0 0.0 1.0 0.0\n"
                                 "GATE   zero    0 O=CONST0;\n"
                                 "GATE   one     0 O=CONST1;\n"
                                 "GATE   ha      5 C=a*b;           PIN * INV 1 999 1.7 0.4 1.7 0.4\n"
                                 "GATE   ha      5 S=!a*b+a*!b;     PIN * INV 1 999 2.1 0.4 2.1 0.4\n"
                                 "GATE   fa      6 C=a*b+a*c+b*c;   PIN * INV 1 999 2.1 0.4 2.1 0.4\n"
                                 "GATE   fa      6 S=a^b^c;         PIN * INV 1 999 3.0 0.4 3.0 0.4";


auto create_resyn_database()
{

  std::string library = "multioutput";
  std::vector<mockturtle::gate> gates;
  std::ifstream in( experiments::cell_libraries_path( library ) );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );

  mockturtle::bound_network bound( gates );
  /* network's definition */

  auto const bs = bound.get_bindings();
  mockturtle::homo_database<4> db( bs );
  using aig_network = mockturtle::aig_network;
  mockturtle::xag_npn_resynthesis<aig_network, aig_network, mockturtle::xag_npn_db_kind::aig_complete> resyn;
  mockturtle::exact_library_params eps;
  eps.np_classification = false;
  mockturtle::exact_library<aig_network> exact_lib( resyn, eps );

  using Tt_t = kitty::dynamic_truth_table;
  using tt_hash = kitty::hash<Tt_t>;
  std::unordered_set<Tt_t, tt_hash> classes;
  Tt_t tt( 4u );
  do
  {
    const auto res = kitty::exact_p_canonization( tt );
    classes.insert( std::get<0>( res ) );
    kitty::next_inplace( tt );
  } while ( !kitty::is_const0( tt ) );

  using klut_network = mockturtle::klut_network;
  using aig_network = mockturtle::aig_network;
  klut_network klut;
  std::vector<typename mockturtle::klut_network::signal> pis;
  for ( int i = 0; i < 4; ++i )
  {
    pis.push_back( klut.create_pi() );
  }

  std::vector<Tt_t> tts;
  for ( auto const& entry : classes )
  {
    auto const f = klut.create_node( pis, entry );
    klut.create_po( f );
    tts.push_back( entry );
  }
  aig_network aig = mockturtle::convert_klut_to_graph<aig_network, klut_network>( klut );

  mockturtle::tech_library_params tps;
  tps.ignore_symmetries = false; // set to true to drastically speed-up mapping with minor delay increase
  tps.verbose = true;
  mockturtle::tech_library<9> tech_lib( gates, tps );
  mockturtle::emap_params mps;
  mps.matching_mode = mockturtle::emap_params::hybrid;
  mps.area_oriented_mapping = true;
  mps.map_multioutput = false;
  mps.relax_required = 0;
  mockturtle::emap_stats mst;
  mockturtle::cell_view<mockturtle::block_network> mapped = mockturtle::emap<9>( aig, tech_lib, mps, &mst );

  int i = 0;
  using node_t = mockturtle::cell_view<mockturtle::block_network>::node;
  std::array<node_t, 4> vars;
  vars[0] = mapped.pi_at( 0 );
  vars[1] = mapped.pi_at( 1 );
  vars[2] = mapped.pi_at( 2 );
  vars[3] = mapped.pi_at( 3 );
  mapped.foreach_po( [&]( auto fo ) {
    kitty::static_truth_table<4> func;
    kitty::create_from_binary_string( func, kitty::to_binary( tts[i++] ) );
    db.add_entry( mapped, vars, fo );
  } );

  return db;
}

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

  tech_library_params tps;
  tps.ignore_symmetries = false; // set to true to drastically speed-up mapping with minor delay increase
  tps.verbose = true;
  tech_library<9> tech_lib( gates, tps );

  bindings binds( gates );
  /* get the database */
  static constexpr uint32_t num_vars = 4;
  static constexpr uint32_t capacity = 10;
  static constexpr mockturtle::choice_order order = power_area_delay;
  using database_t = mockturtle::homo_database<num_vars, capacity, order>;
  database_t database( binds );

  /* run resynthesis */
  for ( auto const& benchmark : all_benchmarks( epfl | iscas ) )
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
    ps.eswp_rounds = 4;
    emap_stats st;
    cell_view<block_network> res = emap<9>( aig, tech_lib, ps, &st );

    auto bound = convert_mapped_to_bound( res, gates );

    uint32_t constexpr N = 30;
    uint32_t constexpr I = 10;
    using Ntk = bound_network;
    using node_simulator_t = mockturtle::node_simulator<Ntk, N, I>;
    
    using dNtk = depth_view<bound_network>;
    dNtk dbound{bound};
    mockturtle::convo_analyzer<dNtk, N, I> power( dbound );
    power(1);

    /* build the analyzer */

    arrival_analyzer<Ntk> arrival1( bound );
    arrival1.compute_arrival_times();
    double const delay_bef = arrival1.worst_delay();

    double const area_bef = bound.area();
    double const zero_delay_bef = power.get_zero_delay();
    double const dyn_switch_bef = power.get_dyn_switch();
    double const glitching_bef = dyn_switch_bef - zero_delay_bef;

    resynthesis_params rps;
    rps.dynamic_database = true;
    rps.try_simula = true;
    rps.try_window = true;
    rps.try_struct = true;
    auto beg = std::chrono::high_resolution_clock::now();
    mockturtle::power_resynthesize<Ntk, database_t, num_vars>( bound, database, rps );
    auto end = std::chrono::high_resolution_clock::now();

    dNtk dbound2{bound};
    mockturtle::convo_analyzer<dNtk, N, I> power2( dbound2 );
    power2(1);
    arrival_analyzer<Ntk> arrival2( bound );
    arrival2.compute_arrival_times();
    double const delay_aft = arrival2.worst_delay();

    double const area_aft = bound.area();
    double const zero_delay_aft = power2.get_zero_delay();
    double const dyn_switch_aft = power2.get_dyn_switch();
    double const glitching_aft = dyn_switch_aft - zero_delay_aft;

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - beg);


    std::cout << "ZD: train " << zero_delay_bef << " " << zero_delay_aft << std::endl;
    std::cout << "DS: train " << dyn_switch_bef << " " << dyn_switch_aft << std::endl;

    const auto cec = benchmark == "hyp" ? true : abc_cec( bound, benchmark );

    /* write verilog netlist */
    // write_verilog_with_cell( res_names, benchmark + "_mapped.v" );
    exp( benchmark, area_bef, area_aft, delay_bef, delay_aft, glitching_bef, glitching_aft, dyn_switch_bef, dyn_switch_aft, duration.count(), cec );
  }

  exp.save();
  exp.table();

  return 0;
}