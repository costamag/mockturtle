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
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/resynthesis/aig_syn.hpp>
#include <mockturtle/resynthesis/algorithms/resynthesize.hpp>
#include <mockturtle/resynthesis/database/homo_database.hpp>
#include <mockturtle/resynthesis/networks/bindings.hpp>
#include <mockturtle/resynthesis/networks/bound.hpp>
#include <mockturtle/resynthesis/database/mapped_list.hpp>
#include <mockturtle/resynthesis/networks/mapped_to_bound.hpp>
#include <mockturtle/resynthesis/list_simulator.hpp>
#include <mockturtle/utils/tech_library.hpp>

#include <experiments.hpp>

enum database_construction_t : uint32_t
{
  db_base,
  db_emap,
  db_lazy
};

enum database_metrics_t : uint32_t
{
  db_adp,
  db_apd,
  db_dap,
  db_dpa,
  db_pad,
  db_pda
};

template<uint32_t num_vars, database_metrics_t metric, class Database_t, class Gates_t>
void fill_with_emap( Database_t&, Gates_t const& gates );

template<uint32_t num_vars, database_metrics_t metric, class Database_t, class Gates_t>
void fill_with_lazy( Database_t&, Gates_t const& gates );

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  /* setup the parameters */
  std::string library = "asap7";
  uint32_t constexpr num_vars = 6;
  uint32_t constexpr capacity = 10;
  database_metrics_t constexpr metric = database_metrics_t::db_adp;
  database_construction_t constr = database_construction_t::db_lazy;
  std::string file = library + "_lazy_adp_" + std::to_string( num_vars ) + ".db";

  /* parse the gates */
  std::vector<mockturtle::gate> gates;
  std::ifstream in( cell_libraries_path( library ));
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  if( result != lorina::return_code::success )
  {
    std::cerr << "[e] error in the library parsing" << std::endl;
    return 1;
  }

  /* build the basic database */
  bindings binds( gates );
  mockturtle::homo_database<num_vars, capacity> database( binds );

  if( constr == database_construction_t::db_emap )
  {
    fill_with_emap<num_vars, metric>( database, gates );
  }
  else if( constr == database_construction_t::db_lazy )
  {
    fill_with_lazy<num_vars, metric>( database, gates );
  }
  database.dump( file );
  return 0;

}

template<uint32_t num_vars, database_metrics_t metric, class Database_t, class Gates_t>
void fill_with_emap( Database_t& database, Gates_t const& gates )
{
  using namespace mockturtle;
  emap_params mps;
  if( ( metric == database_metrics_t::db_adp ) || ( metric == database_metrics_t::db_apd ) )
  {
    mps.area_oriented_mapping = true;
  }
  else if( ( metric == database_metrics_t::db_dap ) || ( metric == database_metrics_t::db_dpa ) )
  {
    mps.area_oriented_mapping = false;
  }
  else
  {
    mps.eswp_rounds = 4u;
  }
  
  using Tt_t = kitty::dynamic_truth_table;
  using tt_hash = kitty::hash<Tt_t>;
  std::unordered_set<Tt_t, tt_hash> classes;
  Tt_t tt( num_vars );
  do {
    const auto res = kitty::exact_p_canonization( tt );
    classes.insert( std::get<0>( res ) );
    kitty::next_inplace( tt );
  } while( !kitty::is_const0( tt ) );

  klut_network klut;
  std::array<typename klut_network::signal, num_vars> pis;
  std::vector<typename klut_network::signal> children;
  for( int i = 0; i < num_vars; ++i )
  {
    pis[i] = klut.create_pi();
    children.push_back( pis[i] );
  }

  for( auto const& entry : classes )
  {
    auto const f = klut.create_node( children, entry );
    klut.create_po( f );
  }
  aig_network aig = convert_klut_to_graph<aig_network, klut_network>( klut );

  xag_npn_resynthesis<aig_network> resyn;
  exact_library_params eps;
  eps.np_classification = false;
  exact_library<aig_network> exact_lib( resyn, eps );

  rewrite_params ps;
  ps.preserve_depth = true;
  rewrite( aig, exact_lib, ps );
  rewrite( aig, exact_lib, ps );
  rewrite( aig, exact_lib, ps );
  
  tech_library_params tps;
  tps.ignore_symmetries = false; // set to true to drastically speed-up mapping with minor delay increase
  tps.verbose = true;
  tech_library<9> tech_lib( gates, tps );
  mps.matching_mode = emap_params::hybrid;
  mps.map_multioutput = false;
  emap_stats mst;
  cell_view<block_network> mapped = emap<9>( aig, tech_lib, mps, &mst );
  std::cout << mapped.compute_area() << " " << mapped.num_pos() << std::endl;

  mapped.foreach_po( [&]( auto fo ){
    database.add_entry( mapped, pis, fo );
  } );
}

template<database_metrics_t metric>
const auto find_order( mockturtle::emap_params & ps )
{
  using namespace mockturtle;
  if constexpr ( metric == db_adp )
  {
    ps.area_oriented_mapping = true;
    return area_delay_power;
  }
  else if constexpr ( metric == db_apd )
  {
    ps.area_oriented_mapping = true;
    return area_power_delay;
  }
  else if constexpr ( metric == db_dap )
  {
    ps.area_oriented_mapping = false;
    return delay_area_power;
  }
  else if constexpr ( metric == db_dpa )
  {
    ps.area_oriented_mapping = false;
    return delay_power_area;
  }
  else if constexpr ( metric == db_pda )
  {
    ps.eswp_rounds = 4u;
    return power_delay_area;
  }
  else
  {
    ps.eswp_rounds = 4u;
    return power_area_delay;
  }
}

template<uint32_t num_vars, database_metrics_t metric, class Database_t, class Gates_t>
void fill_with_lazy( Database_t& database, Gates_t const& gates )
{
  using namespace mockturtle;
  using namespace experiments;

  emap_params ps;
  ps.matching_mode = emap_params::hybrid;
  ps.map_multioutput = false;
  ps.relax_required = 0;

  auto order = find_order<metric>( ps );

  tech_library_params tps;
  tps.ignore_symmetries = false; // set to true to drastically speed-up mapping with minor delay increase
  tps.verbose = true;
  tech_library<9> tech_lib( gates, tps );

  /* run resynthesis */
  for ( auto const& benchmark : all_benchmarks( iscas ) )
  {
    if( benchmark == "hyp" || benchmark == "div" )
      continue;
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    mockturtle::xag_npn_resynthesis<aig_network, aig_network, mockturtle::xag_npn_db_kind::aig_complete> resyn;
    mockturtle::exact_library_params eps;
    eps.np_classification = false;
    mockturtle::exact_library<aig_network> exact_lib( resyn, eps );

    emap_stats st;
    cell_view<block_network> res = emap<9>( aig, tech_lib, ps, &st );

    auto bound = convert_mapped_to_bound( res, gates );

    using Ntk = bound_network;
    
    using dNtk = depth_view<bound_network>;
    dNtk dbound{bound};

    mockturtle::resynthesis_params rps;
    rps.collect_mode = true;
    mockturtle::area_resynthesize<Ntk, Database_t, num_vars>( bound, database, rps );
  }
}