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

#include <mockturtle/generators/random_network.hpp>

#include <experiments.hpp>

using namespace mockturtle;

mockturtle::bound_network random_network( std::vector<gate> gates, int num_pis, int num_pos, int ntk_size, int ntk_index )
{
  std::cout << "ntk_size " << ntk_size << " ntk_index " << ntk_index << std::endl;
  bound_network ntk( gates );
  std::vector<typename bound_network::signal> xs;
  for ( int i = 0; i < num_pis; ++i )
    xs.push_back( ntk.create_pi() );
  std::mt19937 gen( num_pis + num_pos );
  std::uniform_int_distribution<> Ug ( 0, gates.size()-1 );
  while ( ntk.num_gates() < ntk_size )
  {
    std::uniform_int_distribution<> U ( 0, xs.size()-1 );
    int id = Ug( gen );
    auto g = gates[id];
    
    if ( (g.num_vars < 1) || (g.num_vars > xs.size()) )
      continue;
    std::vector<typename bound_network::signal> children;
    std::set<int> idxs;
    while ( children.size() < g.num_vars )
    {
      int idx = U ( gen );
      if ( idxs.find( idx ) == idxs.end() )
      {
        idxs.insert( idx );
        children.push_back( xs[idx] );
      }
    }
    auto fnew = ntk.create_node( children, gates[id].name )[0];
    xs.push_back( fnew );
  }
  std::uniform_int_distribution<> Uint ( 0, xs.size()-1 );
  for ( int i = 0; i < num_pos; ++i )
    ntk.create_po( xs[Uint(gen)] );

  bound_network ntk1( ntk, gates );
  return ntk1;
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
  std::mt19937 gen( 28 );
  bindings binds( gates );
  /* get the database */
  static constexpr uint32_t cube_size = 8;
  static constexpr uint32_t num_steps = 5;
  static constexpr uint32_t num_vars = 4;
  static constexpr uint32_t capacity = 10;
  static constexpr uint32_t Size1 = 200;
  static constexpr uint32_t Size2 = 200;
  using database_t = mockturtle::homo_database<num_vars, capacity>;
  std::string file = "/home/costamag/POST_MAPPING/mockturtle/build/" + library + "_lazy_dap_" + std::to_string( num_vars ) + ".db";
  database_t database( binds, file );// = create_resyn_database<num_vars, capacity, order>();

  mockturtle::placement_params pps;
  /* run resynthesis */
  int ntk_pis = 4;
  int ntk_pos = 3;
  for ( int ntk_size = 5; ntk_size < 21; ntk_size++ ) // experiments::voter ))//
  {
    for ( int ntk_idx = 0; ntk_idx < 10; ntk_idx++ )
    {
      fmt::print( "[i] size {} network {}\n", ntk_size, ntk_idx );
      mockturtle::bound_network ntk = random_network( gates, ntk_pis, ntk_pos, ntk_size, ntk_idx );

      pps.area_norm = ntk.area();
      pps.gates_norm = ntk.size();
      mockturtle::placement_analyzer<mockturtle::bound_network, Size1, Size2> analyzer( ntk, pps );
      double U0_aft = analyzer.compute_cost();
      double U0_bef = U0_aft + 10;
      while ( ( U0_aft < (U0_bef-0.01) ) )
      {
        U0_bef = U0_aft;
        analyzer.fire();
        U0_aft = analyzer.compute_cost();
        std::cout << U0_bef << "-->" << U0_aft << std::endl;
      }
      analyzer.print_costs();
      analyzer.save_result( "/home/costamag/POST_MAPPING/mockturtle/build/graphs/fuzz/tmp0" );

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
      
      ntk.print();
  
      auto beg = std::chrono::high_resolution_clock::now();
      mockturtle::place_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves, Size1, Size2>( ntk, database, rps );
      auto end = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - beg);
      mockturtle::placement_analyzer<mockturtle::bound_network, Size1, Size2> analyzer2( ntk, pps );
  
      double hpwl_aft2 = analyzer2.compute_cost();
      double hpwl_bef2 = hpwl_aft2 + 10;
  
      double U1_aft = analyzer2.compute_cost();
      double U1_bef = U1_aft + 10;
      while ( ( U1_aft < (U1_bef-0.01) ) )
      {
        U1_bef = U1_aft;
        analyzer2.fire();
        U1_aft = analyzer2.compute_cost();
        std::cout << U1_bef << "-->" << U1_aft << std::endl;
      }
      analyzer2.save_result( "/home/costamag/POST_MAPPING/mockturtle/build/graphs/fuzz/tmp1" );
      analyzer2.print_costs();
      if ( U1_aft > U0_aft )
      {
        std::cout << "[d] fuzzed" << std::endl;
        std::cout << "new" << std::endl;
        ntk.print();
        assert(0);
      }
    }
  }

  return 0;
}