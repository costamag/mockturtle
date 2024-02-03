/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2019  EPFL
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

#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/blif.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/rig.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/bench_reader.hpp>
#include <mockturtle/algorithms/rewrite.hpp>

#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/io/write_bench.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/mapping_view.hpp>
#include <mockturtle/algorithms/sim_resub.hpp>


#include <experiments.hpp>

using namespace mockturtle;
using namespace rils;

int main()
{
  using namespace experiments;
  using namespace mockturtle;
  using namespace rils;


  xag_npn_resynthesis<aig_network, aig_network, xag_npn_db_kind::aig_complete> resyn;
  exact_library_params eps;
  eps.np_classification = false;
  exact_library<aig_network, decltype( resyn )> exact_lib( resyn, eps );


  std::vector<kitty::dynamic_truth_table> targets;
  for( uint32_t i{0}; i<8u; ++i )
    targets.emplace_back(8u);

  kitty::create_from_binary_string(targets[0], "1010101010101010000000000000000010101010101010100000000000000000101010101010101000000000000000001010101010101010000000000000000010101010101010100000000000000000101010101010101000000000000000001010101010101010000000000000000010101010101010100000000000000000" );
  kitty::create_from_binary_string(targets[1], "0110011001100110101010101010101011001100110011000000000000000000011001100110011010101010101010101100110011001100000000000000000001100110011001101010101010101010110011001100110000000000000000000110011001100110101010101010101011001100110011000000000000000000" );
  kitty::create_from_binary_string(targets[2], "0001111000011110011001100110011001011010010110101010101010101010101101001011010011001100110011001111000011110000000000000000000000011110000111100110011001100110010110100101101010101010101010101011010010110100110011001100110011110000111100000000000000000000" );
  kitty::create_from_binary_string(targets[3], "0000000111111110000111100001111000111001110001100110011001100110011011011001001001011010010110100101010110101010101010101010101010101011010101001011010010110100100100110110110011001100110011001100011100111000111100001111000011111111000000000000000000000000" );
  kitty::create_from_binary_string(targets[4], "0101010101010100101010110101010001010010100101001011010010110100010010010010010010010011011011000110011011001100110011001100110000110011100110001100011100111000000111000111000011110000111100000000011111000000111111110000000000000000000000000000000000000000" );
  kitty::create_from_binary_string(targets[5], "1001100110011000001100111001100001100011000110001100011100111000100011100011100000011100011100000111100011110000111100001111000011000011111000000000011111000000000111111000000011111111000000001111100000000000000000000000000000000000000000000000000000000000" );
  kitty::create_from_binary_string(targets[6], "1110000111100000110000111110000010000011111000000000011111000000000011111100000000011111100000000111111100000000111111110000000011111100000000001111100000000000111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000" );
  kitty::create_from_binary_string(targets[7], "1111111000000000111111000000000011111100000000001111100000000000111100000000000011100000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000" );


  klut_network klut;
  std::vector<klut_network::signal> inputs;
  for( uint32_t i{0}; i<8; ++i )
  {
    inputs.push_back( klut.create_pi() );
  }
  for( uint32_t i{0}; i<targets.size(); ++i )
  {
    klut.create_po( klut.create_node( inputs, targets[i] ) );
  }

  aig_network aig = convert_klut_to_graph<aig_network>(klut);

  printf("%d\n", aig.num_gates());

  rewrite_params ps;
  rewrite_stats st;

  rewrite( aig, exact_lib, ps, &st );

  printf("%d\n", aig.num_gates());


  static constexpr uint32_t K = 4;
  lut_map_params lps;
  lps.cut_enumeration_ps.cut_size = K;
  lps.cut_enumeration_ps.cut_limit = 8;
  lps.recompute_cuts = true;
  lps.area_oriented_mapping = true;
  lps.cut_expansion = true;
  lut_map_stats lst;
  auto klut1 = lut_map( aig, lps, &lst );
  klut1 = cleanup_luts(klut1);

  std::string tmp = "m8_1.blif";
  write_blif( klut1, tmp );
    
  rig_network rig;
  if ( lorina::read_blif( tmp, blif_reader( rig ) ) != lorina::return_code::success )
  {
    printf("rig unsuccessful\n");
  }

  uint32_t rig_num_gates = rig.num_gates();
  printf("%d\n", rig.num_gates());

  resubstitution_params rps;
  resubstitution_stats rst;
  // ps.pattern_filename = "1024sa1/" + benchmark + ".pat";
  rps.max_inserts = 20;
  rps.max_trials = 100;
  rps.max_pis = 10;
  rps.max_divisors = std::numeric_limits<uint32_t>::max();

  rig_resubstitution<rils::support_selection_t::PIVOT, 4>( rig, rps, &rst );
  rig = cleanup_dangling( rig );
  printf("%d\n", rig.num_gates());

  rig_resubstitution<rils::support_selection_t::PIVOT, 4>( rig, rps, &rst );
  rig = cleanup_dangling( rig );
  printf("%d\n", rig.num_gates());

  rig_resubstitution<rils::support_selection_t::PIVOT, 4>( rig, rps, &rst );
  rig = cleanup_dangling( rig );
  printf("%d\n", rig.num_gates());

  rig_resubstitution<rils::support_selection_t::PIVOT, 4>( rig, rps, &rst );
  rig = cleanup_dangling( rig );
  printf("%d\n", rig.num_gates());

  rig_resubstitution<rils::support_selection_t::PIVOT, 4>( rig, rps, &rst );
  rig = cleanup_dangling( rig );
  printf("%d\n", rig.num_gates());

  tmp = "m8_2.blif";
  write_blif( rig, tmp );
    
  klut_network klut2;
  if ( lorina::read_blif( tmp, blif_reader( klut2 ) ) != lorina::return_code::success )
  {
    printf("klut unsuccessful\n");
  }

  aig_network aig2 = convert_klut_to_graph<aig_network>(klut2);

  printf("%d\n", aig2.num_gates());

  sim_resubstitution(aig2);
  printf("1 %d\n", aig2.num_gates());




  tmp = "m8_2.blif";
  write_blif( aig2, tmp );
  printf("written\n");
  
  rig_network rig2;
  if ( lorina::read_blif( tmp, blif_reader( rig2 ) ) != lorina::return_code::success )
  {
    printf("rig unsuccessful\n");
  }

  printf("%d\n", rig2.num_gates());

  rig_resubstitution<rils::support_selection_t::PIVOT, 3>( rig2, rps, &rst );
  rig2 = cleanup_dangling( rig2 );
  printf("%d\n", rig2.num_gates());

  rig_resubstitution<rils::support_selection_t::PIVOT, 3>( rig2, rps, &rst );
  rig2 = cleanup_dangling( rig2 );
  printf("%d\n", rig2.num_gates());

  rig_resubstitution<rils::support_selection_t::PIVOT, 3>( rig2, rps, &rst );
  rig2 = cleanup_dangling( rig2 );
  printf("%d\n", rig2.num_gates());

  rig_resubstitution<rils::support_selection_t::PIVOT, 3>( rig2, rps, &rst );
  rig2 = cleanup_dangling( rig2 );
  printf("%d\n", rig2.num_gates());

  rig_resubstitution<rils::support_selection_t::PIVOT, 3>( rig2, rps, &rst );
  rig2 = cleanup_dangling( rig2 );
  printf("%d\n", rig2.num_gates());

  tmp = "m8_3.blif";
  write_blif( rig2, tmp );
    
  klut_network klut4;
  if ( lorina::read_blif( tmp, blif_reader( klut4 ) ) != lorina::return_code::success )
  {
    printf("klut unsuccessful\n");
  }

  aig_network aig5 = convert_klut_to_graph<aig_network>(klut4);

  printf("%d\n", aig5.num_gates());

  sim_resubstitution(aig5);
  printf("%d\n", aig5.num_gates());

  tmp = "m8_4.blif";
  write_blif( aig5, tmp );
  printf("written\n");
  
  rig_network rig3;
  if ( lorina::read_blif( tmp, blif_reader( rig3 ) ) != lorina::return_code::success )
  {
    printf("rig unsuccessful\n");
  }

  printf("%d\n", rig3.num_gates());

  rig_resubstitution<rils::support_selection_t::PIVOT, 2>( rig3, rps, &rst );
  rig3 = cleanup_dangling( rig3 );
  printf("%d\n", rig3.num_gates());

  rig_resubstitution<rils::support_selection_t::PIVOT, 2>( rig3, rps, &rst );
  rig3 = cleanup_dangling( rig3 );
  printf("%d\n", rig3.num_gates());

  rig_resubstitution<rils::support_selection_t::PIVOT, 2>( rig3, rps, &rst );
  rig3 = cleanup_dangling( rig3 );
  printf("%d\n", rig3.num_gates());

  tmp = "m8_4.blif";
  write_blif( rig3, tmp );
    
  klut_network klut5;
  if ( lorina::read_blif( tmp, blif_reader( klut5 ) ) != lorina::return_code::success )
  {
    printf("klut unsuccessful\n");
  }

  aig_network aig6 = convert_klut_to_graph<aig_network>(klut5);

  printf("%d\n", aig6.num_gates());

  sim_resubstitution(aig6);
  printf("%d\n", aig6.num_gates());
  sim_resubstitution(aig6);
  printf("%d\n", aig6.num_gates());
  sim_resubstitution(aig6);
  printf("%d\n", aig6.num_gates());

  rewrite( aig6, exact_lib, ps, &st );

  //default_simulator<kitty::dynamic_truth_table> sim( targets[0].num_vars() );
  //const auto tt = simulate<kitty::dynamic_truth_table>( rep, sim )[0];
  //assert( kitty::equal( tt, targets[0] ) );

  return 0;
}
