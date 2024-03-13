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
#include <vector>

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
#include <mockturtle/networks/scg.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/boptimizer.hpp>

#include <experiments.hpp>

using namespace experiments;
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
  using namespace mockturtle;

  experiment<std::string, double, double, double, double, double> exp( "SCOPT", "benchmark", "a(map)", "a(opt)", "d(map)", "d(opt)", "t(opt)");


  fmt::print( "[i] processing technology library\n" );


  /* library to map to technology */
  std::vector<gate> gates;
  std::ifstream in( cell_libraries_path( "sky130" ) );

  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  tech_library_params tps;
  tech_library<5, classification_type::np_configurations> tech_lib( gates, tps );

  for ( auto const& benchmark : epfl_benchmarks( experiments::multiplier ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    scopt::emap2_params ps2;
    ps2.cut_enumeration_ps.minimize_truth_table = true;
    ps2.cut_enumeration_ps.cut_limit = 24;
    ps2.area_flow_rounds=2;
    ps2.area_oriented_mapping = true;
    scopt::emap2_stats st2;


    scopt::scg_network scg_p1 = scopt::emap2_klut( aig, tech_lib, ps2, &st2 );

    double const a_map = scg_p1.compute_area();
    double const d_map = scg_p1.compute_worst_delay();

    printf("a(start) -> %f\n", scg_p1.compute_area());

    std::cout << std::endl;

    /* set up the parameters */
    boptimizer_params rps;
    rps.progress =true;
    rps.max_inserts = 100;
    rps.max_trials = 1;
    rps.max_pis = 10;
    rps.verbose = false;
    rps.max_divisors = 32;

    boptimizer_stats rst_p1;
    double aOld = scg_p1.compute_area()+1;
    while( scg_p1.compute_area() < aOld )
    {
      aOld = scg_p1.compute_area();

      boptimize_sc<scopt::support_selection_t::GREEDY, 4u, 4u>( scg_p1, rps, &rst_p1 );
      scg_p1 = cleanup_dangling( scg_p1 );
      printf("GRE[4,4]: %6f \n", scg_p1.compute_area() );

//      if( aOld == scg_p1.compute_area() )
//      {
//        boptimize_sc<scopt::support_selection_t::GREEDY, 7u, 4u>( scg_p1, rps, &rst_p1 );
//        scg_p1 = cleanup_dangling( scg_p1 );
//        printf("GRE[7,4]: %6d [%6d]\n", scg_p1.num_gates(), scg_p1.max_num_fanins);
//      }
//      if( aOld == scg_p1.compute_area() )
//      {
//        boptimize_sc<scopt::support_selection_t::PIVOT, 4u, 4u>( scg_p1, rps, &rst_p1 );
//        scg_p1 = cleanup_dangling( scg_p1 );
//        printf("PIV[4,4]: %6d [%6d]\n", scg_p1.num_gates(), scg_p1.max_num_fanins);
//      }
//      if( aOld == scg_p1.compute_area() )
//      {
//        boptimize_sc<scopt::support_selection_t::PIVOT, 7u, 4u>( scg_p1, rps, &rst_p1 );
//        scg_p1 = cleanup_dangling( scg_p1 );
//        printf("PIV[7,4]: %6d [%6d]\n", scg_p1.num_gates(), scg_p1.max_num_fanins);
//      }
      printf("%f %f\n", aOld, scg_p1.compute_area());
    }

    double const d_map = scg_p1.compute_area();
    double const d_opt = scg_p1.compute_worst_delay();

    printf("a( end ) -> %f\n", scg_p1.compute_area());
    std::cout << std::endl;

    auto const cec = benchmark == "hyp" ? true : abc_cec( res, benchmark );
    auto const cec2 = benchmark == "hyp" ? true : abc_cec( scg_p1, benchmark );
    if( !cec )  printf("[e] klut not equivalent\n");
    if( !cec2 )  printf("[e] scg not equivalent\n");

  }

  return 0;
}