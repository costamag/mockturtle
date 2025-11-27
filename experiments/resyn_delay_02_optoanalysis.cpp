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

  tech_library_params tps;
  tps.ignore_symmetries = false; // set to true to drastically speed-up mapping with minor delay increase
  tps.verbose = true;
  tech_library<9> tech_lib( gates, tps );
  bool do_cont = false;
  /* run resynthesis */
  for ( auto const& benchmark : all_benchmarks( iscas ) ) // experiments::voter ))//
  {
    //if( benchmark == "square" )
    //  do_cont = true;
    if( benchmark == "hyp" )
      continue;
    fmt::print( "[i] processing {}\n", benchmark );
    if ( do_cont )
      continue;
    names_view<aig_network> aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }
    if( aig.num_gates()>100000 )
    continue;

    emap_params ps;
    ps.matching_mode = emap_params::hybrid;
    ps.area_oriented_mapping = false;
    ps.map_multioutput = false;
    emap_stats st;

    aig_balancing_params bps;
    bps.minimize_levels = true;
    bps.fast_mode = true;
    aig_balance( aig, bps );
    aig = abc_opto( aig, "op1"+benchmark, "resyn2rs; if -g; balance;" );

    cell_view<block_network> res = emap<9>( aig, tech_lib, ps, &st );
    auto ntk = convert_mapped_to_bound( res, gates );

    double area_bef = res.compute_area();
    double delay_bef = res.compute_worst_delay();

    uint32_t constexpr num_steps = 12;
    uint32_t constexpr cube_size = 12;
    uint32_t constexpr max_leaves = 10;

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
    std::cout << "a1" << std::endl;
    
    arrival_analyzer<bound_network> arrival( ntk );
    arrival.compute_arrival_times();
    uint32_t num_crit = 0;
    uint32_t num_tot = 0;
    ntk.foreach_gate( [&]( auto n ){
      num_tot += 1;
      if ( arrival.get_arrival( n ) >= ( arrival.get_required( n ) - 0.1 ) )
        num_crit += 1;
    } );
    std::cout << num_crit << "/" << num_tot << std::endl;
  }

  return 0;
}