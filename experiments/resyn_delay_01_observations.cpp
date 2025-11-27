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

#include <mockturtle/algorithms/aig_balancing.hpp>
#include <mockturtle/algorithms/cut_enumeration.hpp>
#include <mockturtle/algorithms/cut_enumeration/rewrite_cut.hpp>
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


template<class Gate_t>
double derive_mux_delay( std::vector<Gate_t> const& gates );

template<class Gate_t>
uint32_t derive_num_equioptimizable( std::vector<Gate_t> const& gates );

template<class Ntk>
void extract_arrival_patterns( Ntk & ntk, double mux_delay, std::ofstream & fout );

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  /* setup the parameters */
  std::string library = "asap7";
  std::ofstream fout( "patterns.eatp" );

  /* parse the gates */
  std::vector<mockturtle::gate> gates;
  std::ifstream in( cell_libraries_path( library ));
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  if( result != lorina::return_code::success )
  {
    std::cerr << "[e] error in the library parsing" << std::endl;
    return 1;
  }

  auto const mux_delay = derive_mux_delay( gates );
  std::cout << "The worst delay is " << mux_delay << std::endl;

  uint32_t num_equi = derive_num_equioptimizable( gates );

  /* Structural cut enumeration and statistics collection */
  tech_library_params tps;
  tps.ignore_symmetries = false; // set to true to drastically speed-up mapping with minor delay increase
  tech_library<9> tech_lib( gates, tps );

  /* run resynthesis */
  for ( auto const& benchmark : all_benchmarks( epfl ) ) // experiments::voter ))//
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    aig_balancing_params bps;
    bps.minimize_levels = true;
    bps.fast_mode = true;
    aig_balance( aig, bps );

    emap_params ps;
    ps.matching_mode = emap_params::hybrid;
    ps.area_oriented_mapping = false;
    ps.map_multioutput = false;
    emap_stats st;
    cell_view<block_network> res = emap<9>( aig, tech_lib, ps, &st );

    auto ntk = convert_mapped_to_bound( res, gates );

    extract_arrival_patterns( ntk, mux_delay, fout );
  }

  fout.close();

  return 0;
}

template<class Gate_t>
double derive_mux_delay( std::vector<Gate_t> const& gates )
{
  using namespace mockturtle;
  emap_params mps;
  mps.area_oriented_mapping = false;

  using Tt_t = kitty::dynamic_truth_table;
  Tt_t tt( 3 );
  kitty::create_from_binary_string( tt, "10101100" );

  klut_network klut;
  std::vector<typename klut_network::signal> xs;
  for( int i = 0; i < 3; ++i )
  {
    xs.push_back( klut.create_pi() );
  }

  auto const f = klut.create_node( xs, tt );
  klut.create_po( f );
  aig_network aig = convert_klut_to_graph<aig_network, klut_network>( klut );

  tech_library_params tps;
  tps.ignore_symmetries = false; // set to true to drastically speed-up mapping with minor delay increase
  tps.verbose = true;
  tech_library<9> tech_lib( gates, tps );
  mps.matching_mode = emap_params::hybrid;
  mps.map_multioutput = false;
  emap_stats mst;
  cell_view<block_network> mapped = emap<9>( aig, tech_lib, mps, &mst );

  return mapped.compute_worst_delay();
}

template<class Gate_t>
uint32_t derive_num_equioptimizable( std::vector<Gate_t> const& gates )
{
  double mux_delay = derive_mux_delay( gates );
  double min_delay = std::numeric_limits<double>::max();
  for( auto g : gates )
  {
    for( int i = 0; i < g.num_vars; ++i )
    {
      auto const pin_delay = std::min( g.pins[i].rise_block_delay, g.pins[i].fall_block_delay );
      if( pin_delay < min_delay )
      {
        min_delay = pin_delay;
      }
    }
  }
  uint32_t delta = std::ceil( mux_delay / min_delay );
  std::cout << "ratio is " << mux_delay / min_delay ;
  std::cout << " -> delta = " << delta << std::endl;
  uint32_t num_equioptimizable = std::ceil(pow(( 1 + 4*mux_delay / min_delay ), 4));
  std::cout << "num_equioptimizable = " << num_equioptimizable << std::endl;
  return num_equioptimizable;
}

template<class Ntk>
void extract_arrival_patterns( Ntk & ntk, double mux_delay, std::ofstream & fout )
{
  using namespace mockturtle;
  arrival_analyzer<bound_network> arrival( ntk );
  arrival.compute_arrival_times();

  using network_cuts_t = dynamic_network_cuts<Ntk, 4, true, cut_enumeration_rewrite_cut>;
  using cut_manager_t = detail::dynamic_cut_enumeration_impl<Ntk, 4, true, cut_enumeration_rewrite_cut>;
  using cut_t = typename network_cuts_t::cut_t;
  cut_enumeration_params ps;
  ps.cut_limit = 13;
  ps.minimize_truth_table = false;
  cut_enumeration_stats st;
  network_cuts_t cuts( ntk.size() + ( ntk.size() >> 1 ));
  cut_manager_t cut_manager( ntk, ps, st, cuts );

  ntk.foreach_gate( [&]( auto n, auto i ){
    cut_manager.clear_cuts( n );
    cut_manager.compute_cuts( n );

    for ( auto& cut : cuts.cuts( ntk.node_to_index( n ) ) )
    {
      if( cut->size() == 4 )
      {
        std::vector<double> times;
        std::vector<typename Ntk::node> leaves;
        for ( auto const leaf : *cut )
        {
          times.push_back( arrival.get_time( leaf ) );
          leaves.push_back( leaf );
        }

        /* STEP 1 : find the maximum arrival time in the cut */
        double max_T = 0;
        for ( int i = 0; i < 4; ++i )
        {
          if( times[i] > max_T )
          {
            max_T = times[i];
          }
        }
        /* STEP 2 : bring forward the patterns that are too low */
        double delta = 4*mux_delay;
        double min_delay = std::numeric_limits<double>::max();
        for ( int i = 0; i < 4; ++i )
        {
          if( times[i] < ( max_T - delta ) )
          {
            times[i] = max_T - delta;
          }
          min_delay = std::min( min_delay, times[i] );
        }

        /* STEP 3 : shift everything to have the lowest value to zero */
        for ( int i = 0; i < 4; ++i )
        {
          times[i] -= min_delay;
        }

       for ( int i = 0; i < 4; ++i )
        {
          fout << times[i] << " ";
        }
        fout << std::endl;
      }
    }
  } );
  

}