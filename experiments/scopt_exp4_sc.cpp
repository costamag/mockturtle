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

template<class Ntk>
Ntk abc_opto( Ntk const& ntk, std::string str_code, std::string abc_script = "resyn2rs" )
{
  write_aiger( ntk, "/tmp/"+str_code+".aig" );
  std::string command = "abc -q \"r /tmp/"+str_code+".aig; fraig; " + abc_script + "; write_aiger /tmp/" + str_code + ".aig\"";

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

  experiment<std::string, double, double, double, double, double, bool> exp( "SCOPT", "benchmark", "a(map)", "a(opt)", "d(map)", "d(opt)", "t(opt)", "cec");

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

  for ( auto const& benchmark : all_benchmarks( ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    if( aig.num_gates() > 300000 ) continue;
    std::vector<uint32_t> aig_size;
    std::vector<double> map_size;
    std::vector<uint32_t> aig_depth;
    std::vector<double> map_delay;
    std::vector<uint32_t> vheuristic;
    uint32_t heuristic=0;


    scopt::scg_network scg;

    double best_area = std::numeric_limits<double>::max();
    double best_delay = std::numeric_limits<double>::max();

    uint32_t num_old = aig.num_gates()+1;
    while( num_old > aig.num_gates() )
    {
      num_old = aig.num_gates();
      aig = abc_opto( aig, benchmark, "compress2rs" );
      aig = cleanup_dangling( aig );
      printf( "aig>%d\n", aig.num_gates() );
    }
        scopt::emap2_params ps;
        ps.cut_enumeration_ps.minimize_truth_table = true;
        ps.cut_enumeration_ps.cut_limit = 24;
        ps.area_flow_rounds=2;
        ps.area_oriented_mapping = true;
        scopt::emap2_stats st;
//
//        aig_network res2;
//
//        aig = abc_opto( aig, benchmark, "rw" );
//        aig = cleanup_dangling( aig );
//        if( aig.num_gates() != num_old ){ vheuristic.push_back(0); }
//
//        if( aig.num_gates() == num_old )
//        {
//            aig = abc_opto( aig, benchmark, "rs" );
//            aig = cleanup_dangling( aig );
//            if( aig.num_gates() != num_old ){ vheuristic.push_back(1); }
//
//        }
//        if( aig.num_gates() == num_old )
//        {
//            aig = abc_opto( aig, benchmark, "rf" );
//            aig = cleanup_dangling( aig );
//            if( aig.num_gates() != num_old ){ vheuristic.push_back(2); }
//
//        }
//        if( aig.num_gates() == num_old )
//        {
//            aig = abc_opto( aig, benchmark, "resyn2rs" );
//            aig = cleanup_dangling( aig );
//            if( aig.num_gates() != num_old ){ vheuristic.push_back(3); }
//
//        }
//        if( aig.num_gates() == num_old )
//        {


//            if( aig.num_gates() != num_old ){ vheuristic.push_back(4); }
//        }

        scopt::scg_network res = emap2_klut( aig, tech_lib, ps, &st );
        double area = res.compute_area();
//        if( area < best_area )
//        {
            vheuristic.push_back(5);
            best_area = area;
            best_delay = res.compute_worst_delay();
            scg = res;
      printf( "map>%f\n", scg.compute_area() );

//        }
//
//
//        printf( "%d -> %f %f", aig.num_gates(), st.area, st.delay );
//        std::cout << std::endl;
//        depth_view<aig_network> daig {aig};

//        aig_size.push_back( aig.num_gates() );
//        aig_depth.push_back( daig.depth() );
//        map_size.push_back( res.compute_area() );
//        map_delay.push_back( res.compute_worst_delay() );
 //   }

    write_aiger( aig, benchmark+"_optmap.aig" );

//    printf("aaig=np.array([");
//    for( int i{0}; i< aig_size.size()-1; ++i )
//    {
//        if( i == aig_size.size()-2 )
//            printf("%d])\n", aig_size[i]);
//        else
//            printf("%d, ", aig_size[i]);
//    }
//
//    printf("amap=np.array([");
//    for( int i{0}; i< map_size.size()-1; ++i )
//    {
//        if( i == map_size.size()-2 )
//            printf("%f])\n", map_size[i]);
//        else
//            printf("%f, ", map_size[i]);
//    }
//
//    printf("color=np.array([");
//    for( int i{0}; i< vheuristic.size()-1; ++i )
//    {
//        if( i == vheuristic.size()-2 )
//            printf("%d])\n", vheuristic[i]);
//        else
//            printf("%d, ", vheuristic[i]);
//    }
//
//    printf("d(aig)=[");
//    for( int i{0}; i< aig_depth.size()-1; ++i )
//    {
//        if( i == aig_depth.size()-2 )
//            printf("%d])\n", aig_depth[i]);
//        else
//            printf("%d, ", aig_depth[i]);
//    }
//
//    printf("d(map)=[");
//    for( int i{0}; i< map_delay.size()-1; ++i )
//    {
//        if( i == map_delay.size()-2 )
//            printf("%f])\n", map_delay[i]);
//        else
//            printf("%f, ", map_delay[i]);
//    }
//

    const auto cec = benchmark == "hyp" ? true : abc_cec( aig, benchmark );
    assert( cec && "[e] not equivalent" );


    boptimizer_params rps;
    rps.progress =true;
    rps.max_inserts = 300;
    rps.max_trials = 1;
    rps.max_pis = 16;
    rps.verbose = false;
    rps.max_divisors = 128;

    boptimizer_stats rst_p1;
    int effort{0};
    double aOld = scg.compute_area()+1;
    while( scg.compute_area() < aOld )
    {
      aOld = scg.compute_area();

      //if( aOld == scg.compute_area() && effort <= 0 )
      //{
        boptimize_sc<scopt::support_selection_t::GREEDY, 4u, 4u>( scg, rps, &rst_p1 );
        printf("GRE[4,4]: %6f ", scg.compute_area() );
//        scg = cleanup_dangling( scg );
//        printf("GRE[4,4]: %6f ", scg.compute_area() );
        std::cout << std::endl;
        effort = 0;
      //}
//      if( aOld == scg.compute_area() && effort <= 1 )
//      {
//        boptimize_sc<scopt::support_selection_t::GREEDY, 7u, 4u>( scg, rps, &rst_p1 );
//        scg = cleanup_dangling( scg );
//        printf("GRE[7,4]: %6f \n", scg.compute_area() );
//        std::cout << std::endl;
//        effort = 1;
//      }
      //if( aOld == scg.compute_area() && effort <= 2 )
      //{
//        boptimize_sc<scopt::support_selection_t::PIVOT, 4u, 4u>( scg, rps, &rst_p1 );
//        scg = cleanup_dangling( scg );
//        printf("PIV[4,4]: %6f \n", scg.compute_area() );
//        std::cout << std::endl;
//        effort = 2;
//      //}
//      if( aOld == scg.compute_area() && effort <= 3 )
//      {
//        boptimize_sc<scopt::support_selection_t::PIVOT, 7u, 4u>( scg, rps, &rst_p1 );
//        scg = cleanup_dangling( scg );
//        printf("PIV[7,4]: %6f \n", scg.compute_area() );
//        std::cout << std::endl;
//        effort = 3;
//      }
      printf("%f %f\n", aOld, scg.compute_area());
    }

    double const d_map = scg.compute_area();
    double const d_opt = scg.compute_worst_delay();

    printf("a( end ) -> %f\n", scg.compute_area());
    std::cout << std::endl;

    auto const cecMP = benchmark == "hyp" ? true : abc_cec( scg, benchmark );
    if(!cecMP)
      printf("ERROR\n");
    std::cout << std::endl;

    exp( benchmark, best_area, scg.compute_area(), best_delay, scg.compute_worst_delay(), 0.0, cecMP );

  }

  exp.save();
  exp.table();
  return 0;
}
//ISCAS compress2rs
// benchmark   &     a(map)   &     a(opt)   &    d(map)   &    d(opt)   &    t(opt) \\
//       c17   & $    25.03 $ & $    25.03 $ & $  209.99 $ & $  209.99 $ & $   0.00 $\\
//      c432   & $   833.41 $ & $\color{green}   800.87 $ & $ 2372.10 $ & $ 2478.70 $ & $   0.00 $\\
//      c499   & $  2220.34 $ & $  2220.34 $ & $ 1734.40 $ & $ 1734.40 $ & $   0.00 $\\
//      c880   & $  1669.46 $ & $  1669.46 $ & $ 2267.41 $ & $ 2267.41 $ & $   0.00 $\\
//     c1355   & $  2259.14 $ & $  2259.14 $ & $ 1621.88 $ & $ 1621.88 $ & $   0.00 $\\
//     c1908   & $  1823.47 $ & $\color{green}  1775.91 $ & $ 1983.17 $ & $ 1983.17 $ & $   0.00 $\\
//     c2670   & $  3193.57 $ & $\color{green}  3083.49 $ & $ 1954.60 $ & $ 1954.60 $ & $   0.00 $\\
//     c3540   & $  4455.14 $ & $\color{green}  4324.99 $ & $ 2793.62 $ & $ 2793.62 $ & $   0.00 $\\
//     c5315   & $  7408.77 $ & $\color{green}  7293.63 $ & $ 2757.67 $ & $ 2602.58 $ & $   0.00 $\\
//     c6288   & $ 11673.98 $ & $ 11673.98 $ & $ 7377.53 $ & $ 7377.53 $ & $   0.00 $\\
//     c7552   & $  7906.75 $ & $\color{green}  7655.18 $ & $ 4097.25 $ & $ 4417.80 $ & $   0.00 $\\

//opt=np.array([25.03,800.87,2220.34,1669.46,2259.14,1775.91,3083.49,4324.99,7293.63,11673.98,7655.18])
//map=np.array([25.03,833.41,2220.34,1669.46,2259.14,1823.47,3193.57,4455.14,7408.77,11673.98,7906.75])

//ISCAS compress2rs X INF

//| benchmark |   a(map) |   a(opt) |  d(map) |  d(opt) | t(opt) |
//|       c17 & $    25.03 $ & $\color{colk}    25.03 $ & $  209.99 $ & $  209.99 $ & $  0.00 $\\
//|      c432 & $   833.41 $ & $\color{colg}   800.87 $ & $ 2372.10 $ & $ 2478.70 $ & $  0.00 $\\
//|      c499 & $  2289.17 $ & $\color{colk}  2289.17 $ & $ 1690.19 $ & $ 1690.19 $ & $  0.00 $\\
//|      c880 & $  1664.45 $ & $\color{colg}  1635.66 $ & $ 2267.41 $ & $ 2267.41 $ & $  0.00 $\\
//|     c1355 & $  2240.38 $ & $\color{colk}  2240.38 $ & $ 1708.91 $ & $ 1708.91 $ & $  0.00 $\\
//|     c1908 & $  1637.05 $ & $\color{colk}  1637.05 $ & $ 2026.36 $ & $ 2026.36 $ & $  0.00 $\\
//|     c2670 & $  3122.24 $ & $\color{colg}  3080.96 $ & $ 1841.29 $ & $ 1841.29 $ & $  0.00 $\\
//|     c3540 & $  4481.43 $ & $\color{colg}  4466.41 $ & $ 2932.07 $ & $ 2932.07 $ & $  0.00 $\\
//|     c5315 & $  7396.22 $ & $\color{colg}  7346.15 $ & $ 3000.88 $ & $ 3000.88 $ & $  0.00 $\\
//|     c6288 & $ 11650.20 $ & $\color{colk} 11650.20 $ & $ 7377.53 $ & $ 7377.53 $ & $  0.00 $\\
//|     c7552 & $  7814.19 $ & $\color{colg}  7558.87 $ & $ 4736.52 $ & $ 4836.40 $ & $  0.00 $\\

//opt=np.array([  25.03,   800.87,  2289.17,  1635.66,  2240.38,  1637.05,  3080.96,  4466.41,  7346.15, 11650.20,  7558.87])
//map=np.array([    25.03,   833.41,  2289.17,  1664.45,  2240.38,  1637.05,  3122.24,  4481.43,  7396.22, 11650.20,  7814.19])