/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
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
#include <fstream>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/genlib.hpp>
#include <mockturtle/algorithms/emap2.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/scg.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/algorithms/cleanup.hpp>


#include <experiments.hpp>

using namespace experiments;
using namespace mockturtle;


int main()
{
  using namespace experiments;
  using namespace mockturtle;

  /* library to map to technology */
  std::vector<gate> gates;
  std::ifstream in( cell_libraries_path( "asap7" ) );

  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  tech_library_params tps;
  tech_library<5, classification_type::np_configurations> tech_lib( gates, tps );

  xag_npn_resynthesis<aig_network, aig_network, xag_npn_db_kind::aig_complete> resyn;
  exact_library_params eps;
  eps.np_classification = false;
  exact_library<aig_network> exact_lib( resyn, eps );

  //rewrite( aig, exact_lib );

  using TT = kitty::dynamic_truth_table;
  using tt_hash = kitty::hash<TT>;

  std::unordered_set<TT, tt_hash> classes;
  TT tt(4u);
  do
  {
    const auto res = kitty::exact_p_canonization( tt );
    classes.insert( std::get<0>( res ) );
    kitty::next_inplace( tt );
  } while ( !kitty::is_const0( tt ) );

    scopt::emap2_params ps2;
    scopt::emap2_stats st2;
    ps2.required_time = std::numeric_limits<float>::max();
    ps2.area_oriented_mapping= true;
    //ps2.cut_enumeration_ps.cut_limit = 24;

    std::vector<TT> tts;
    std::vector<double> areas;
    std::vector<std::vector<uint32_t>> id_lists;

    int i{0};
    for ( auto const& entry : classes )
    {
        printf("%d ", i++);

        kitty::print_binary( entry );
        printf(" ");


        aig_network aig;
        std::vector<aig_network::signal> pis;
        pis.push_back(aig.create_pi());
        pis.push_back(aig.create_pi());
        pis.push_back(aig.create_pi());
        pis.push_back(aig.create_pi());

        resyn( aig, entry, pis.begin(), pis.end(), [&]( auto const& f_new ) {
        return aig.create_po(f_new);
        });

        scopt::scg_network scg = scopt::emap2_klut( aig, tech_lib, ps2, &st2 );

        default_simulator<kitty::dynamic_truth_table> sim( 4 );
        const auto tt = simulate<kitty::dynamic_truth_table>( scg, sim )[0];
        if( tt._bits[0] != entry._bits[0] )
        {
          printf("ERROR\n");
        }

        std::vector<uint32_t> idlist;
        scg.foreach_gate( [&]( auto n ) 
          { 
            idlist.push_back( scg.fanin_size(n) );

            scg.foreach_fanin( n, [&]( auto const& fi ) {
              idlist.push_back( fi.index );
            } );
            idlist.push_back( scg.get_binding(n).id ); 
          } );

        tts.push_back( entry );
        id_lists.push_back( idlist );
        areas.push_back( scg.compute_area() );

        for( auto lit : id_lists.back() )
        {
            printf("%d ", lit );
        }
        printf("- %.2f\n", areas.back() );

    }

    std::ofstream fTt;
    fTt.open ("asap7_2.tts");
    for( auto tt : tts )
    {
        fTt << kitty::to_binary(tt)+ "\n";
    }
    fTt.close();

    std::ofstream fList;
    fList.open ("asap7_2.list");
    for( auto list : id_lists )
    {
        for( auto lit : list )
        {
            fList << lit << " ";
        }
        fList << "\n";
    }
    fList.close();

    std::ofstream fArea;
    fArea.open ("asap7_2.area");
    for( auto area : areas )
    {
        fArea << area << "\n";
    }
    fArea.close();


  return 0;
}


