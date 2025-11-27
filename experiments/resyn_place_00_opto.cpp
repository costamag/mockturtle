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
#include <mockturtle/resynthesis/placement/placement_analyzer.hpp>
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
  std::string command = "abc -q \"read_aiger /tmp/"+str_code+".aig; " + abc_script + "; write_aiger /tmp/" + str_code + ".aig\"";

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

  experiment<std::string, double, double, double, double, double, double, double, double, float, bool> exp(
      "HWPL", "benchmark", "A(place)", "A(resyn)", "D(place)", "D(resyn)", "U(place)", "U(resyn)", "W(place)", "W(resyn)", "t(det)", "cec" );



  /* library to map to technology */
  //fmt::print( "[i] processing technology library\n" );
  std::string library = "asap7";
  std::vector<gate> gates;
  std::ifstream in( cell_libraries_path( "asap7" ) );

  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  bindings binds( gates );
  /* get the database */
  static constexpr uint32_t cube_size = 8;
  static constexpr uint32_t num_steps = 5;
  static constexpr uint32_t num_vars = 4;
  static constexpr uint32_t capacity = 10;
  static constexpr uint32_t Size1 = 200;
  static constexpr uint32_t Size2 = 200;
  using database_t = mockturtle::easy_database<num_vars, capacity>;
  std::string file = "/home/costamag/POST_MAPPING/mockturtle/build/" + library + "_lazy_dap_" + std::to_string( num_vars ) + ".db";
  database_t database( binds, file );// = create_resyn_database<num_vars, capacity, order>();

  mockturtle::placement_params pps;
  /* get the database */
  tech_library_params tps;
  tps.ignore_symmetries = false; // set to true to drastically speed-up mapping with minor delay increase
  tps.verbose = true;
  tech_library<9> tech_lib( gates, tps );
  bool do_cont = false;
  /* run resynthesis */
  for ( auto const& benchmark : all_benchmarks( iscas ) ) // experiments::voter ))//
  {
    fmt::print( "[i] processing {}\n", benchmark );
    if ( benchmark == "div" || benchmark == "hyp" )
      continue;
    names_view<aig_network> aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }
  //  if ( aig.num_gates() > 2000 )
  //    continue;

    emap_params ps;
    ps.matching_mode = emap_params::hybrid;
    ps.area_oriented_mapping = false;
    ps.map_multioutput = false;
    emap_stats st;

    aig_balancing_params bps;
    bps.minimize_levels = true;
    bps.fast_mode = true;
    aig = abc_opto( aig, "op1"+benchmark, "resyn2rs; resyn2rs; if -g; if -g; balance;" );
    cell_view<block_network> res = emap<9>( aig, tech_lib, ps, &st );
    auto ntk = convert_mapped_to_bound( res, gates );
    auto libr = ntk.get_bindings();

    binding_view bnd{res, gates};
    write_verilog_with_binding( ntk, benchmark + "_p_i.v" );
    pps.area_norm = ntk.area();
    pps.gates_norm = ntk.size();

    mockturtle::placement_analyzer<mockturtle::bound_network, Size1, Size2> analyzer( ntk, pps );
    double hpwl_bef = analyzer.compute_cost();
    double hpwl_aft1 = analyzer.compute_cost();
    double hpwl_bef1 = hpwl_aft1 + 10;
    int s = 0;

  double NUM=1;

  double U1=0;
  double W1=0;
  double R1=0;
  double P1=0;
  //  for( int i = 0; i < NUM; i++ )
    while( ( hpwl_aft1 < (hpwl_bef1-1) ) )
    {
      hpwl_bef1 = hpwl_aft1;
      analyzer.fire();
      hpwl_aft1 = analyzer.compute_cost();
      std::cout << hpwl_bef1 << "-->" << hpwl_aft1 << std::endl;
    }
    W1 = analyzer.compute_hpwl();
    U1 = analyzer.compute_cost();
    analyzer.save_result( "/home/costamag/POST_MAPPING/mockturtle/build/graphs/"+benchmark+"_opto_0" );

    //analyzer.print();

    std::cout << "OPTO" << std::endl;

    mockturtle::resynthesis_params rps;
    rps.preserve_depth = true;
    rps.dynamic_database = false;
    rps.try_struct = true;
    rps.try_window = true;
    rps.try_simula = true;
    rps.try_rewire = true;
    rps.analyzer_ps.placement_ps = pps;
    uint32_t constexpr max_leaves = 10;
    using Ntk = mockturtle::bound_network;
    double Abef = ntk.area();
    arrival_analyzer<Ntk> arrival_bef( ntk );
    arrival_bef.compute_arrival_times();;
    double Dbef = arrival_bef.worst_delay();

    auto beg = std::chrono::high_resolution_clock::now();
    mockturtle::place_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves, Size1, Size2>( ntk, database, rps );
    auto end = std::chrono::high_resolution_clock::now();

    double Aaft = ntk.area();
    arrival_analyzer<Ntk> arrival_aft( ntk );
    arrival_aft.compute_arrival_times();
    double Daft = arrival_aft.worst_delay();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - beg);
    mockturtle::placement_analyzer<mockturtle::bound_network, Size1, Size2> analyzer2( ntk, pps );

    double hpwl_aft2 = analyzer2.compute_cost();
    double hpwl_bef2 = hpwl_aft2 + 10;

    s=0;
    double U2 = 0;
    double W2 = 0;
    double P2 = 0;
    double R2 = 0;
    //for( int i = 0; i < NUM; ++i )
    while( hpwl_aft2 < (hpwl_bef2-1) )
    {
      hpwl_bef2 = hpwl_aft2;
      analyzer2.fire();
      hpwl_aft2 = analyzer2.compute_cost();
      std::cout << hpwl_bef2 << "-->" << hpwl_aft2 << std::endl;
    }
    W2 = analyzer2.compute_hpwl();
    U2 = analyzer2.compute_cost();
    analyzer2.print_params_placer();
    analyzer2.save_result( "/home/costamag/POST_MAPPING/mockturtle/build/graphs/"+benchmark+"_opto_1" );

    std::cout << "A " << Abef << "-->" << Aaft << std::endl;
    std::cout << "D " << Dbef << "-->" << Daft << std::endl;
    std::cout << "U " << U1 << "-->" << U2 << std::endl;
    std::cout << "W " << W1 << "-->" << W2 << std::endl;

    std::cout << "***************************************************" << std::endl;
    std::cout << "***************************************************" << std::endl;
    std::cout << "***************************************************" << std::endl;
    write_verilog_with_binding( ntk, benchmark + "_p_e.v" );

    const auto cec = benchmark == "hyp" ? true : abc_cec( ntk, benchmark );

    exp( benchmark, Abef, Aaft, Dbef, Daft, U1, U2, W1, W2, 0.001*duration.count(), cec );
  
  }

  exp.save();
  exp.table();

  return 0;
}