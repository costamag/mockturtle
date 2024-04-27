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
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/scg.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/boptimizer.hpp>

#include <experiments.hpp>

using namespace experiments;
using namespace mockturtle;

aig_network abc_if( aig_network const& ntk, std::string str_code, uint32_t K=4u )
{
  write_aiger( ntk, "/tmp/"+str_code+".aig" );
  std::string command = "abc -q \"r /tmp/"+str_code+".aig; dch -f; if -g; strash; fraig; write_aiger /tmp/" + str_code + ".aig\"";

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

  aig_network res;

  std::string string_path = ("/tmp/"+str_code+".aig");
  if( lorina::read_aiger( string_path, aiger_reader( res ) ) != lorina::return_code::success )
    std::cerr << "read_blif failed" << std::endl;

  return res;
}

template<class Ntk>
Ntk abc_opto( Ntk const& ntk, std::string str_code, std::string abc_script = "resyn2rs" )
{
  write_aiger( ntk, "/tmp/"+str_code+".aig" );
  std::string command = "abc -q \"r /tmp/"+str_code+".aig; fraig;" + abc_script + "; write_aiger /tmp/" + str_code + ".aig\"";

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

  fmt::print( "[i] processing technology library\n" );


  /* library to map to technology */
  std::vector<gate> gates;
  std::ifstream in( cell_libraries_path( "asap7" ) );

  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  tech_library_params tps;
  tech_library<5, classification_type::np_configurations> tech_lib( gates, tps );

  for ( auto const& benchmark : epfl_benchmarks( ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    std::vector<uint32_t> aig_size;
    std::vector<uint32_t> aig_depth;
    std::vector<double> map_size;
    std::vector<double> map_delay;
    std::vector<double> opt_size;
    std::vector<double> opt1_size;
    std::vector<double> optN_size;
    std::vector<double> opt_delay;
    std::vector<uint32_t> vheuristic;
    uint32_t heuristic=0;
    bool start=true;
    bool close=false;

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    

    #pragma region initial
    scopt::emap2_params ps;
    ps.cut_enumeration_ps.minimize_truth_table = true;
    ps.cut_enumeration_ps.cut_limit = 24;
    ps.area_flow_rounds=2;
    ps.area_oriented_mapping = false;
    scopt::emap2_stats st;

    scopt::scg_network scg = emap2_klut( aig, tech_lib, ps, &st );

    printf("a(map) -> %f\n", scg.compute_area());
    printf("d(map) -> %f\n", scg.compute_worst_delay());

    depth_view<aig_network> daig {aig};
    aig_size.push_back( aig.num_gates() );
    aig_depth.push_back( daig.depth() );
    map_size.push_back( scg.compute_area() );
    map_delay.push_back( scg.compute_worst_delay() );
    vheuristic.push_back(10);

    boptimizer_params rps;
    rps.progress =true;
    rps.max_inserts = 300;
    rps.max_trials = 10;
    rps.max_pis = 16;
    rps.verbose = false;
    rps.max_divisors = 128;
    boptimizer_stats rst_p1;

    double aOld = scg.compute_area() + 1;
    while( aOld > scg.compute_area() )
    {
      aOld = scg.compute_area();

      //boptimize_sc<scopt::support_selection_t::NGREEDY, 4u, 4u>( scg, rps, &rst_p1 );
      boptimize_sc<PV2, 4u, 4u>( scg, rps, &rst_p1 );
      scg = cleanup_scg( scg );
      printf("GRE[4,4]: %6f %6f", scg.compute_area(), scg.compute_worst_delay() );

  //    boptimize_sc<scopt::support_selection_t::PIVOT, 6u, 4u>( scg, rps, &rst_p1 );
  //    scg = cleanup_dangling( scg );
  //    printf("PIV[6,4]: %6f ", scg.compute_area() );
  //    std::cout << std::endl;
    }

    printf("a(start) -> %f\n", scg.compute_area());
    printf("d(start) -> %f\n", scg.compute_worst_delay());

    opt_size.push_back(scg.compute_area());
    opt_delay.push_back(scg.compute_worst_delay());

    std::cout << std::endl;

    auto const cecMP = benchmark == "hyp" ? true : abc_cec( scg, benchmark );
    if(!cecMP)
      printf("ERROR\n");
    std::cout << std::endl;

    #pragma endregion initial


    depth_view<aig_network> daig1{aig};
    uint32_t depth_old = daig1.depth()+1;
    uint32_t depth_new = daig1.depth();

    int it=0;
    while( depth_old > depth_new )
    {
        depth_old = depth_new;
        it++;
        aig_network res2;

        aig = abc_if(aig, benchmark);
        aig=cleanup_dangling(aig);
        aig = abc_opto(aig, benchmark, "resyn2rs");

        scopt::scg_network scg = emap2_klut( aig, tech_lib, ps, &st );

        printf( "%d -> %f %f", aig.num_gates(), st.area, st.delay );
        std::cout << std::endl;
        

        aig_size.push_back( aig.num_gates() );
        aig_depth.push_back( daig.depth() );
        map_size.push_back( scg.compute_area() );
        map_delay.push_back( scg.compute_worst_delay() );

     //   double aOld = scg.compute_area()+1;
     //   while( aOld > scg.compute_area() )
     //   {
      //    aOld = scg.compute_area();
        depth_view<aig_network> daig2{aig};
        depth_new = daig2.depth();
        if( depth_old >= depth_new )
        {
          
          scopt::scg_network scg = emap2_klut( aig, tech_lib, ps, &st );

          boptimizer_params rps;
          rps.progress =true;
          rps.max_inserts = 300;
          rps.max_trials = 1;
          rps.max_pis = 16;
          rps.verbose = false;
          rps.max_divisors = 300; 

          double aOld = scg.compute_area() + 1;
          int it=0;
          while( aOld > scg.compute_area() && it++ < 3 )
          {

            aOld = scg.compute_area();
            boptimize_sc<PV2, 4u, 4u>( scg, rps, &rst_p1 );
            scg = cleanup_scg( scg );
            printf("Ex3[4,4]: %6f %6f", scg.compute_area(), scg.compute_worst_delay() );

//            boptimize_sc<scopt::support_selection_t::PIVOT, 6u, 4u>( scg, rps, &rst_p1 );
//            scg = cleanup_dangling( scg );
//            printf("PIV[6,4]: %6f ", scg.compute_area() );
//            std::cout << std::endl;
          }

          std::cout << std::endl;
          start=false;

          printf("a( end ) -> %f\n", scg.compute_area());
          printf("d( end ) -> %f\n", scg.compute_worst_delay());
          std::cout << std::endl;
          opt_size.push_back(scg.compute_area());
          opt_delay.push_back(scg.compute_worst_delay());
          auto const cecMP = benchmark == "hyp" ? true : abc_cec( scg, benchmark );
          if(!cecMP)
            printf("ERROR\n");
          std::cout << std::endl;
        }
    //    }
//
    //    aOld = scg.compute_area()+1;
    //    while( aOld > scg.compute_area() )
    //    {
    //      aOld = scg.compute_area();
    //      boptimize_sc<scopt::support_selection_t::PIVOT, 4u, 4u>( scg, rps, &rst_p1 );
    //      scg = cleanup_dangling( scg );
    //      printf("PIV[4,4]: %6f ", scg.compute_area() );
    //      std::cout << std::endl;
    //    }

    //    double aOld = scg.compute_area()+1;
    //    while( aOld > scg.compute_area() )
    //    {
    //      aOld = scg.compute_area();
    //      boptimize_sc<scopt::support_selection_t::GREEDY, 7u, 4u>( scg, rps, &rst_p1 );
    //      scg = cleanup_dangling( scg );
    //      printf("GRE[7,4]: %6f ", scg.compute_area() );
    //      std::cout << std::endl;
    //    }

      //  aOld = scg.compute_area()+1;
      //  while( aOld > scg.compute_area() )
      //  {
      //    aOld = scg.compute_area();
      //    boptimize_sc<scopt::support_selection_t::PIVOT, 7u, 4u>( scg, rps, &rst_p1 );
      //    scg = cleanup_dangling( scg );
      //    printf("PIV[7,4]: %6f ", scg.compute_area() );
      //    std::cout << std::endl;
      //  }


    }

    write_aiger( aig, benchmark+"_optmap.aig" );

    printf("aaig=np.array([");
    for( int i{0}; i< aig_size.size()-1; ++i )
    {
        if( i == aig_size.size()-2 )
            printf("%d])\n", aig_size[i]);
        else
            printf("%d, ", aig_size[i]);
    }

    printf("amap=np.array([");
    for( int i{0}; i< map_size.size()-1; ++i )
    {
        if( i == map_size.size()-2 )
            printf("%f])\n", map_size[i]);
        else
            printf("%f, ", map_size[i]);
    }

    printf("color=np.array([");
    for( int i{0}; i< vheuristic.size()-1; ++i )
    {
        if( i == vheuristic.size()-2 )
            printf("%d])\n", vheuristic[i]);
        else
            printf("%d, ", vheuristic[i]);
    }

    printf("aopt=np.array([");
    for( int i{0}; i< opt_size.size()-1; ++i )
    {
        if( i == opt_size.size()-2 )
            printf("%f])\n", opt_size[i]);
        else
            printf("%f, ", opt_size[i]);
    }

    printf("d(aig)=[");
    for( int i{0}; i< aig_depth.size()-1; ++i )
    {
        if( i == aig_depth.size()-2 )
            printf("%d])\n", aig_depth[i]);
        else
            printf("%d, ", aig_depth[i]);
    }

    printf("d(map)=[");
    for( int i{0}; i< map_delay.size()-1; ++i )
    {
        if( i == map_delay.size()-2 )
            printf("%f])\n", map_delay[i]);
        else
            printf("%f, ", map_delay[i]);
    }

    printf("d(opt)=[");
    for( int i{0}; i< opt_delay.size()-1; ++i )
    {
        if( i == opt_delay.size()-2 )
            printf("%f])\n", opt_delay[i]);
        else
            printf("%f, ", opt_delay[i]);
    }

    const auto cec = benchmark == "hyp" ? true : abc_cec( aig, benchmark );
    assert( cec && "[e] not equivalent" );




  }

  return 0;
}