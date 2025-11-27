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

#include <SFML/Graphics.hpp>
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


//class renderer_t
//{
//public:
//  renderer_t( sf::RenderTarget& target_ )
//    : target{ target_ }
//    {}
//
//  void render()
//  {
//    target.clear( sf::Color::Black );
//    sf::CircleShape circle{ 1.0f };
//    circle.setPointCount( 32 );
//    circle.setPosition( { 420.0f, 420.0f } );
//    circle.setScale( 10.0f, 10.0f );
//    circle.setFillColor( sf::Color::Blue );
//    target.draw( circle );
//  }
//private:
//  sf::RenderTarget& target;
//};

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

  mockturtle::placement_params pps;
  pps.area_density = 0.04;
  pps.num_placers = 1;

  experiment<std::string, double, double> exp(
      "HWPL", "benchmark", "cost(bef)", "cost(aft)" );

  /* run resynthesis */
  for ( auto const& benchmark : all_benchmarks( experiments::c880 ) ) // experiments::voter ))//
  {
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
  static constexpr uint32_t Nx = 100;
  static constexpr uint32_t Ny = 100;

  /* get the database */
  tech_library_params tps;
  tps.ignore_symmetries = false; // set to true to drastically speed-up mapping with minor delay increase
  tps.verbose = true;
  tech_library<9> tech_lib( gates, tps );
  bool do_cont = false;
    fmt::print( "[i] processing {}\n", benchmark );
    names_view<aig_network> aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    emap_params ps;
    ps.matching_mode = emap_params::hybrid;
    ps.area_oriented_mapping = false;
    ps.map_multioutput = false;
    emap_stats st;

    aig_balancing_params bps;
    bps.minimize_levels = true;
    bps.fast_mode = true;
    aig_balance( aig, bps );
    aig = abc_opto( aig, "op1"+benchmark, "resyn2rs; resyn2rs; balance; if -g; if -g;" );
    cell_view<block_network> res = emap<9>( aig, tech_lib, ps, &st );
    auto ntk = convert_mapped_to_bound( res, gates );
    write_verilog_with_binding( ntk, benchmark + "_exp0.v" );

    mockturtle::placement_analyzer<mockturtle::bound_network, Nx, Ny> analyzer( ntk, pps );
    double cost_i = analyzer.compute_cost();
    
    double cost_aft = analyzer.compute_cost();
    double cost_bef = cost_aft + 10;
    //analyzer.print_params_placer();

    analyzer.save_result( "/home/costamag/POST_MAPPING/mockturtle/build/graphs/"+benchmark+"_0" );
    while( cost_aft < (cost_bef - 1) )
    {
      cost_bef = cost_aft;
      analyzer.fire();
      cost_aft = analyzer.compute_cost();
      std::cout << cost_bef << "-fire->" << cost_aft << std::endl;
    }
    analyzer.save_result( "/home/costamag/POST_MAPPING/mockturtle/build/graphs/"+benchmark+"_1" );
    analyzer.save_nets_stats( "/home/costamag/POST_MAPPING/mockturtle/build/graphs/"+benchmark+"_nets" );
    double cost_e = analyzer.compute_cost();

    exp( benchmark, cost_i, cost_e );
  }


  exp.save();
  exp.table();

  return 0;
}