/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2023  EPFL
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

#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>

#include <mockturtle/algorithms/emap.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/rewrite.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/resynthesis/aig_syn.hpp>
#include <mockturtle/resynthesis/database/homo_database.hpp>
#include <mockturtle/resynthesis/database/mapped_list.hpp>
#include <mockturtle/resynthesis/networks/bindings.hpp>
#include <mockturtle/resynthesis/networks/bound.hpp>
#include <mockturtle/resynthesis/list_simulator.hpp>
#include <mockturtle/utils/tech_library.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  std::vector<mockturtle::gate> gates;

  std::string library = "asap7";
  std::ifstream in( cell_libraries_path( library ));
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  if( result != lorina::return_code::success )
  {
    std::cerr << "[e] error in the library parsing" << std::endl;
    return 1;
  }

  bindings binds( gates );
  mockturtle::homo_database<6> database( binds );

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
  ps.preserve_depth = true;
  std::cout << aig.num_gates() << " " << aig.num_pos() << std::endl;
  rewrite( aig, exact_lib, ps );
  rewrite( aig, exact_lib, ps );
  rewrite( aig, exact_lib, ps );
  std::cout << aig.num_gates() << " " << aig.num_pos() << std::endl;
  

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
  std::cout << mapped.compute_area() << " " << mapped.num_pos() << std::endl;

  mapped.foreach_po( [&]( auto fo ){
    database.add_entry( mapped, pis, fo );
  } );
  write_verilog_with_cell( mapped, library + "_database.v" );

  return 0;
}