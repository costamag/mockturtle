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
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/block.hpp>
#include <mockturtle/networks/klut.hpp>
 
#include <mockturtle/utils/name_utils.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/cell_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/names_view.hpp>
 
#include <mockturtle/resynthesis/aig_syn.hpp>
#include <mockturtle/resynthesis/algorithms/resynthesize.hpp>
#include <mockturtle/resynthesis/algorithms/glitch_substitution.hpp>
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
 
  experiment<std::string, double, double, double, double, double, double, double, double,double, double, double,double, double, double, double, float, float, bool, bool> exp(
      "DDECS25-PANDA", "benchmark", "A(bef)", "A(sota)","A(opt)", "D(bef)", "D(sota)","D(opt)", "G(bef)", "G(sota)", "G(opt)", "S(bef)", "S(sota)", "S(opt)", "P(bef)", "P(sota)", "P(opt)", "t(sota)", "t(opt)", "cec1", "cec2" );
 
  /* Library to map to technology */
  std::string library = "asap7";
  std::vector<gate> gates;
  std::ifstream in( cell_libraries_path( library ) );
 
  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }
 
  bindings binds( gates );
 
  /* Load the database */
  static constexpr uint32_t num_vars = 4;
  static constexpr uint32_t capacity = 10;
  static constexpr mockturtle::choice_order order = power_area_delay;
  using database_t = mockturtle::homo_database<num_vars, capacity, order>;
  std::string path = "/home/costamag/POST_MAPPING/mockturtle/build/" + library + "_lazy_pad_" + std::to_string( num_vars ) + ".db";
  database_t database( binds, path );
 
  /* Resynthesis on the benchmarks */
  for ( auto const& benchmark : all_benchmarks( iscas | epfl ) )//experiments::c499 | experiments::c1355 | experiments::arbiter | experiments::router ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    //if ( (benchmark == "div")||(benchmark == "hyp")||(benchmark == "log2") )
    //  continue;
 
    names_view<aig_network> aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }
 
    /* Technology-independent optimization */
    mockturtle::xag_npn_resynthesis<aig_network, aig_network, mockturtle::xag_npn_db_kind::aig_complete> resyn;
    mockturtle::exact_library_params eps;
    eps.np_classification = false;
    mockturtle::exact_library<aig_network> exact_lib( resyn, eps );
    /* Nodes-count minimization to reduce switching elements */
    mockturtle::rewrite_params rwps;
    rwps.preserve_depth = true;
    mockturtle::rewrite( aig, exact_lib, rwps );
    /* Paths-balancing to reduce glitches before mapping */
    aig_balancing_params bps;
    bps.minimize_levels = true;
    bps.fast_mode = true;
    aig_balance( aig, bps );
 
    /* Technology mapping */
    tech_library_params tps;
    tps.ignore_symmetries = false;
    tps.verbose = false;
    tech_library<9> tech_lib( gates, tps );
    emap_params mps;
    mps.matching_mode = emap_params::hybrid;
    mps.area_oriented_mapping = true;
    mps.eswp_rounds = 3;
    mps.map_multioutput = false;
    emap_stats mst;
    cell_view<block_network> res = emap<9>( aig, tech_lib, mps, &mst );
     
    /* Network to be optimized with Panda's method */
    bound_network ntk1( gates );
    bound_network ntk2( gates );
    /* Network to be optimized with our glitching-aware resynthesis */
    convert_mapped_to_bound( ntk1, res, gates );
    convert_mapped_to_bound( ntk2, res, gates );
 
    uint32_t constexpr num_steps_eval = 20;// good one is 20; 10
    uint32_t constexpr cube_size_eval = 12;// good one is 12 10
    uint32_t constexpr num_steps = 12;//good one is 12 10
    uint32_t constexpr cube_size = 12;//good one is 12 10
 
    using Ntk = bound_network;
     
    using dNtk = depth_view<bound_network>;
    dNtk dbound{ntk1};
 
    using node_simulator_t = mockturtle::node_simulator<Ntk, num_steps_eval, cube_size_eval>;
 
    /* Create the input simulations */
    std::vector<node_simulator_t> pis_sims;
    std::default_random_engine::result_type seed = 1;
    std::array<kitty::static_truth_table<cube_size_eval>, 2> workload;
    int i = 1;
    ntk1.foreach_pi([&]( auto n ){
      kitty::create_random( workload[0], seed*i );
      kitty::create_random( workload[1], seed*( i + 1 ) );
      i += 2;
      node_simulator_t sim( n, workload, -1, 0 );
      pis_sims.push_back( sim );
    });
 
    double zero_delay_init;
    double dyn_switch_init;
    double dyn_power_init;
    double glitching_init;
    double delay_init;
    double area_init;
     
     write_verilog_with_binding( ntk1, "glitch/" + benchmark + "_o0_2.v" );
     {
        dNtk dbound{ntk1};
        mockturtle::power_analyzer<Ntk, num_steps_eval, cube_size_eval> power_analyzer_init( ntk1 );
        power_analyzer_init( pis_sims );
 
        zero_delay_init = power_analyzer_init.get_zero_delay();
        dyn_switch_init = power_analyzer_init.get_dyn_switch();
        dyn_power_init = power_analyzer_init.get_dyn_power();
        glitching_init = dyn_switch_init - zero_delay_init;
        arrival_analyzer<Ntk> arrival_analyzer_init( ntk1 );
        arrival_analyzer_init.compute_arrival_times();
        delay_init = arrival_analyzer_init.worst_delay();
        area_init = ntk1.area();
      }
      std::cout << "ZD(0) " << zero_delay_init << std::endl;
      std::cout << "GL(0) " << glitching_init << std::endl;
      std::cout << "DS(0) " << dyn_switch_init << std::endl;
      std::cout << "PW(0) " << dyn_power_init << std::endl;
     
      /* State-of-the-art glitching minimization */
      mockturtle::resynthesis_params rps;
      mockturtle::resub_params sps;
      rps.preserve_depth = false;
      rps.dynamic_database = false;
      
      uint32_t constexpr max_leaves_sota = 10;
      rps.try_struct = false;
      rps.try_window = false;
      rps.try_simula = false;
      rps.try_rewire = true;
      rps.preserve_depth = true;
      auto start_time_sota = std::chrono::high_resolution_clock::now();
      mockturtle::glitch_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves_sota>( ntk1, database, rps );
      mockturtle::glitch_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves_sota>( ntk1, database, rps );
      //mockturtle::glitch_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves_sota>( ntk1, database, rps );
      //auto ntk1_2 = mockturtle::cleanup_bound( ntk1, gates );
      //mockturtle::glitch_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves_sota>( ntk1, database, rps );
      //auto ntk1_3 = mockturtle::cleanup_bound( ntk1_2, gates );
      auto end_time_sota = std::chrono::high_resolution_clock::now();
      write_verilog_with_binding( ntk1, "glitch/" + benchmark + "_panda_sota_2.v" );
 
      /* Proposed glitching minimization */
      rps.try_struct = true;
      rps.try_window = true;
      rps.try_simula = true;
      rps.try_rewire = true;
      rps.preserve_depth = true;

      uint32_t constexpr max_leaves_ours = 10;
      auto start_time_ours = std::chrono::high_resolution_clock::now();
      //ntk2.print();
      mockturtle::glitch_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves_ours>( ntk2, database, rps );
      mockturtle::glitch_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves_ours>( ntk2, database, rps );
      //ntk2.print();
      //mockturtle::glitch_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves_ours>( ntk2, database, rps );
      //auto ntk2_2 = mockturtle::cleanup_bound( ntk2, gates );
      //mockturtle::glitch_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves_ours>( ntk2, database, rps );
      //auto ntk2_3 = mockturtle::cleanup_bound( ntk2_2, gates );
      auto end_time_ours = std::chrono::high_resolution_clock::now();
      write_verilog_with_binding( ntk2, "glitch/" + benchmark + "_panda_ours_2.v" ); 
 
      double delay_sota;
      double area_sota;
      double zero_delay_sota;
      double dyn_switch_sota;
      double dyn_power_sota;
      double glitching_sota;
      double delay_ours;
      double area_ours;
      double zero_delay_ours;
      double dyn_switch_ours;
      double dyn_power_ours;
      double glitching_ours;

      {
        dNtk dresult1{ntk1};
        arrival_analyzer<Ntk> arrival_analyzer_sota( ntk1 );
        arrival_analyzer_sota.compute_arrival_times();
        delay_sota = arrival_analyzer_sota.worst_delay();
        area_sota = ntk1.area();
        mockturtle::power_analyzer<Ntk, num_steps_eval, cube_size_eval> power_analyzer_sota( ntk1 );
        power_analyzer_sota( pis_sims );
        zero_delay_sota = power_analyzer_sota.get_zero_delay();
        dyn_switch_sota = power_analyzer_sota.get_dyn_switch();
        dyn_power_sota = power_analyzer_sota.get_dyn_power();
        glitching_sota = dyn_switch_sota - zero_delay_sota;
      }
      {
        dNtk dresult2{ntk2};
        arrival_analyzer<Ntk> arrival_analyzer_ours( ntk2 );
        arrival_analyzer_ours.compute_arrival_times();
        delay_ours = arrival_analyzer_ours.worst_delay();
        area_ours = ntk2.area();
        mockturtle::power_analyzer<Ntk, num_steps_eval, cube_size_eval> power_analyzer_ours( ntk2 );
        power_analyzer_ours( pis_sims );
        zero_delay_ours = power_analyzer_ours.get_zero_delay();
        dyn_switch_ours = power_analyzer_ours.get_dyn_switch();
        dyn_power_ours = power_analyzer_ours.get_dyn_power();
        glitching_ours = dyn_switch_ours - zero_delay_ours;
      }

      printf( "AR | init : %7.2f", area_init );
      printf( "   | sota : %7.2f", area_sota );
      printf( "   | ours : %7.2f\n", area_ours );
      printf( "DE | init : %7.2f", delay_init );
      printf( "   | sota : %7.2f", delay_sota );
      printf( "   | ours : %7.2f\n", delay_ours );
      printf( "ZD | init : %7.2f", zero_delay_init );
      printf( "   | sota : %7.2f", zero_delay_sota );
      printf( "   | ours : %7.2f\n", zero_delay_ours );
      printf( "GL | init : %7.2f", glitching_init );
      printf( "   | sota : %7.2f", glitching_sota );
      printf( "   | ours : %7.2f\n", glitching_ours );
      printf( "DS | init : %7.2f", dyn_switch_init );
      printf( "   | sota : %7.2f", dyn_switch_sota );
      printf( "   | ours : %7.2f\n", dyn_switch_ours );
      printf( "DP | init : %7.2f", dyn_power_init );
      printf( "   | sota : %7.2f", dyn_power_sota );
      printf( "   | ours : %7.2f\n", dyn_power_ours );

      auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_sota - start_time_sota);
      auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_ours - start_time_ours);
  
      const bool cec1 = benchmark == "hyp" ? true : abc_cec( ntk1, benchmark );
      const bool cec2 = benchmark == "hyp" ? true : abc_cec( ntk2, benchmark );
 
      exp( benchmark, area_init, area_sota, area_ours, delay_init, delay_sota, delay_ours, glitching_init, glitching_sota, glitching_ours, dyn_switch_init, dyn_switch_sota, dyn_switch_ours, dyn_power_init, dyn_power_sota, dyn_power_ours, duration1.count()*0.001, duration2.count()*0.001, cec1, cec2 );
   }
 
    exp.save();
    exp.table();
 
    return 0;
 }