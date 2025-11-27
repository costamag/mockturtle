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
#include <mockturtle/io/write_aiger.hpp>
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
//#include <mockturtle/resynthesis/analyzers/homo_xpx_analyzer.hpp>
#include <mockturtle/resynthesis/power/node_simulator.hpp>
#include <mockturtle/resynthesis/power/analyzers/power_analyzer.hpp>
#include <mockturtle/resynthesis/database/easy_database.hpp>
#include <mockturtle/resynthesis/database/mapped_list.hpp>
#include <mockturtle/resynthesis/list_simulator.hpp>
#include <mockturtle/resynthesis/networks/bindings.hpp>
#include <mockturtle/resynthesis/networks/bound.hpp>
#include <mockturtle/resynthesis/database/homo_storage.hpp>
#include <mockturtle/resynthesis/networks/mapped_to_bound.hpp>
#include <mockturtle/resynthesis/networks/cleanup_bound.hpp>
#include <mockturtle/resynthesis/power/analyzers/convo_analyzer.hpp>

#include <experiments.hpp>

using namespace mockturtle;

template<class Ntk>
Ntk abc_opto( Ntk const& ntk, std::string str_code, std::string abc_script = "resyn2rs" )
{
  write_aiger( ntk, "/tmp/"+str_code+".aig" );
  std::string command = "abc -q \"r /tmp/"+str_code+".aig; " + abc_script + "; write_aiger /tmp/" + str_code + ".aig\"";

  std::array<char, 128> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }

  Ntk res;
  std::string string_path = ("/tmp/"+str_code+".aig");
  if( lorina::read_aiger( string_path, aiger_reader( res ) ) != lorina::return_code::success )
    std::cerr << "read_aiger failed" << std::endl;

  return res;
}

int main()
{
  using namespace experiments;
  using namespace mockturtle;//klut_to_bound

  experiment<std::string, double, double, double, double, double, double, float, float, bool, bool> exp(
      "delay", "benchmark", "A(bef)", "A(sota)", "A(ours)", "D(bef)", "D(sota)", "D(ours)", "T(sota)", "T(ours)", "cec1", "cec2" );

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
  using database_t = mockturtle::easy_database<num_vars, capacity>;
  std::string file1 = "/home/costamag/POST_MAPPING/mockturtle/build/" + library + "_enum3_" + std::to_string( num_vars ) + ".db";
  std::string file2 = "/home/costamag/POST_MAPPING/mockturtle/build/" + library + "_enum6_" + std::to_string( num_vars ) + ".db";

  database_t database1( binds, file1 );// = create_resyn_database<num_vars, capacity, order>();
  database_t database2( binds, file2 );// = create_resyn_database<num_vars, capacity, order>();
  tech_library_params tps;
  tps.ignore_symmetries = false; // set to true to drastically speed-up mapping with minor delay increase
  tps.verbose = true;
  tech_library<9> tech_lib( gates, tps );
  bool do_cont = false;
  double avg_area_sota=0;
  double avg_delay_sota=0;
  double avg_area_ours=0;
  double avg_delay_ours=0;
  double avg_time_ours=0;
  double avg_time_sota=0;
  double number = 0;
  /* run resynthesis */
  for ( auto const& benchmark : all_benchmarks( epfl ) ) // experiments::voter ))//
  {
    //if( benchmark == "square" )
    //  do_cont = true;
    //if( benchmark == "hyp"||benchmark == "sqrt" )
    //  continue;
    fmt::print( "[i] processing {}\n", benchmark );
    //if ( do_cont )
    //  continue;
    names_view<aig_network> aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    emap_params ps;
    ps.matching_mode = emap_params::hybrid;
    ps.area_oriented_mapping = true;
    ps.map_multioutput = false;
    emap_stats st;

    aig_balancing_params bps;
    bps.minimize_levels = true;
    bps.fast_mode = true;
    //aig_balance( aig, bps );
    aig = abc_opto( aig, "op1"+benchmark, "resyn2rs; resyn2rs; if -g;" );
    if( aig.num_gates()>300000 )
      continue;
    std::cout << "done opto" << std::endl;
    cell_view<block_network> res = emap<9>( aig, tech_lib, ps, &st );
    auto ntk = convert_mapped_to_bound( res, gates );
    auto ntk2 = convert_mapped_to_bound( res, gates );

    double area_bef = res.compute_area();
    double delay_bef = res.compute_worst_delay();

    uint32_t constexpr num_steps = 10;
    uint32_t constexpr cube_size = 10;

    using Ntk = bound_network;
    
    using dNtk = depth_view<bound_network>;
    dNtk dbound{ntk};

    mockturtle::resynthesis_params rps;
    rps.preserve_depth = true;
    rps.dynamic_database = false;

    auto beg_sota = std::chrono::high_resolution_clock::now();
    rps.try_struct = true;
    rps.try_window = false;
    rps.try_simula = false;
    rps.decompose = 0;
    uint32_t constexpr max_leaves1 = 8;
    uint32_t constexpr max_leaves2 = 8;
    uint32_t constexpr max_leaves3 = 8;
    uint32_t constexpr max_leaves4 = 8;
    uint32_t constexpr max_leaves5 = 8;
    uint32_t constexpr max_leaves6 = 8;
    uint32_t constexpr max_leaves7 = 8;
    mockturtle::delay_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves1, 0>( ntk, database1, rps );
    mockturtle::delay_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves2, 1>( ntk, database1, rps );
    mockturtle::delay_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves3, 2>( ntk, database1, rps );
    mockturtle::delay_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves4, 2>( ntk, database1, rps );
    //mockturtle::delay_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves5, 4>( ntk, database1, rps );
    //mockturtle::delay_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves6, 5>( ntk, database1, rps );
    //mockturtle::delay_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves7, 6>( ntk, database1, rps );

    //mockturtle::delay_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves, 2>( ntk, database1, rps );
    
    arrival_analyzer<bound_network> arrival( ntk );
    arrival.compute_arrival_times();
    double area_sota = ntk.area();
    double delay_sota = arrival.worst_delay();
    auto end_sota = std::chrono::high_resolution_clock::now();
    auto duration_sota = std::chrono::duration_cast<std::chrono::milliseconds>(end_sota - beg_sota);
    float time_sota = duration_sota.count()*0.001;

    auto beg_ours = std::chrono::high_resolution_clock::now();
    rps.try_struct = true;
    rps.try_window = true;
    rps.try_simula = false;
    rps.try_rewire = false;
    rps.decompose = 1;

    mockturtle::delay_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves1, 0>( ntk2, database2, rps );
    mockturtle::delay_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves2, 1>( ntk2, database2, rps );
    mockturtle::delay_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves3, 2>( ntk2, database2, rps );
    mockturtle::delay_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves4, 2>( ntk2, database2, rps );
    //mockturtle::delay_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves5, 2>( ntk2, database2, rps );
    //mockturtle::delay_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves6, 5>( ntk2, database2, rps );
    //mockturtle::delay_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves7, 6>( ntk2, database2, rps );
    
    dNtk dntk2{ntk2};
    rps.try_struct = false;
    
    arrival_analyzer<bound_network> arrival2( ntk2 );
    arrival2.compute_arrival_times();
    double area_ours = ntk2.area();
    double delay_ours = arrival2.worst_delay();
    auto end_ours = std::chrono::high_resolution_clock::now();
    auto duration_ours = std::chrono::duration_cast<std::chrono::milliseconds>(end_ours - beg_ours);
    float time_ours = duration_ours.count()*0.001;

    const auto cec = benchmark == "hyp" ? true : abc_cec( ntk, benchmark );
    const auto cec2 = benchmark == "hyp" ? true : abc_cec( ntk2, benchmark );

    std::cout << "A: " << area_bef << " " << area_sota << " " << area_ours << std::endl;
    std::cout << "D: " << delay_bef << " " << delay_sota << " " << delay_ours << std::endl;
    std::cout << "E: " << cec << " " << cec2 << std::endl;
    exp( benchmark, area_bef, 100*(area_sota-area_bef)/area_bef, 100*(area_ours-area_bef)/area_bef, delay_bef, 100*(delay_sota-delay_bef)/delay_bef, 100*(delay_ours-delay_bef)/delay_bef, time_sota, time_ours, cec, cec2 );
  
    avg_area_sota+=(area_sota-area_bef)/area_bef;
    avg_delay_sota+=(delay_sota-delay_bef)/delay_bef;
    avg_area_ours+=(area_ours-area_bef)/area_bef;
    avg_delay_ours+=(delay_ours-delay_bef)/delay_bef;
    avg_time_sota += time_sota;
    avg_time_ours += time_ours;
    number += 1;
  
  }

  std::cout << " avg_area_sota  " << avg_area_sota/number  << std::endl;
  std::cout << " avg_area_ours  " << avg_area_ours/number  << std::endl;
  std::cout << " avg_delay_sota " << avg_delay_sota/number << std::endl;
  std::cout << " avg_delay_ours " << avg_delay_ours/number << std::endl;
  std::cout << " avg_time_sota " << avg_time_sota/number << std::endl;
  std::cout << " avg_time_ours " << avg_time_ours/number << std::endl;

  exp.save();
  exp.table();

  return 0;
}