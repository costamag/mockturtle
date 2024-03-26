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
#include <ctime>
#include <experiments.hpp>

using namespace experiments;
using namespace mockturtle;
using namespace std::chrono;

aig_network abc_if( aig_network const& ntk, std::string str_code, uint32_t K=4u )
{
  write_aiger( ntk, "/tmp/"+str_code+".aig" );
  std::string command = "abc -q \"r /tmp/"+str_code+".aig; dch -f; if -g; strash; dfraig; write_aiger /tmp/" + str_code + ".aig\"";

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
  std::string command = "abc -q \"r /tmp/"+str_code+".aig;" + abc_script + "; write_aiger /tmp/" + str_code + ".aig\"";

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

  experiment<std::string, double, double, double, double, double, double, double, double, bool> exp( "SCOPTA", "benchmark", "a(map)", "a(opt1)", "a(optN)", "d(map)", "d(opt1)", "d(optN)", "t(opt1)", "t(optN)", "cec");

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
    if( benchmark == "hyp" ) continue;


    fmt::print( "[i] processing {}\n", benchmark );

    bool start=true;
    bool close=false;

    aig_network aig;
    if ( lorina::read_aiger( experiments::benchmark_path(benchmark), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }
    auto aaig_old = aig.num_gates()+1;
    auto aaig_new = aig.num_gates();
    depth_view<aig_network> daig{aig};
    uint32_t daig_old = daig.depth()+1;
    uint32_t daig_new = daig.depth();

    if( aig.num_gates()>300000 ) continue;

    //if( benchmark != "hyp" )
    {
      aig = abc_opto( aig, benchmark, "rw; rs;" );
      aig = abc_opto( aig, benchmark, "rw; rs;" );
      aig = abc_opto( aig, benchmark, "rw; rs;" );
      aig = abc_opto( aig, benchmark, "rw; rs;" );
      aig = abc_opto( aig, benchmark, "resyn2rs" );
      aig = abc_opto( aig, benchmark, "compress2rs" );
      while( aaig_new < aaig_old )
      {
        aig = abc_opto( aig, benchmark, "resyn2rs" );
        aig = cleanup_dangling( aig );
        aig = abc_opto( aig, benchmark, "compress2rs" );
        aig = cleanup_dangling( aig );

        aaig_old = aaig_new;
        aaig_new = aig.num_gates();
        printf("%d\n", aaig_new);

      }
    }

    const auto cec = benchmark == "hyp" ? true : abc_cec( aig, benchmark );
    assert( cec && "[e] not equivalent" );

    scopt::emap2_params ps;
    ps.cut_enumeration_ps.minimize_truth_table = true;
    ps.cut_enumeration_ps.cut_limit = 24;
    ps.area_flow_rounds=2;
    ps.area_oriented_mapping = true;
    scopt::emap2_stats st;

    printf("map..\n");

    scopt::scg_network scg = emap2_klut( aig, tech_lib, ps, &st );
    scg = cleanup_scg( scg );

    double const aold = scg.compute_area();
    double const dold = scg.compute_worst_delay();

    printf("a0)%6f ", aold );
    std::cout << std::endl;
    printf("d0)%6f ", dold );
    std::cout << std::endl;

    boptimizer_params rps;
    rps.progress =true;
    rps.max_inserts = 300;
    rps.max_trials = 10;
    rps.max_pis = 16;
    rps.verbose = false;
    rps.use_delay_constraints = false;
    rps.max_divisors = 128;

    boptimizer_stats rst_p1;

    double aold1 = scg.compute_area();
    double dold1 = scg.compute_worst_delay();


    std::clock_t begin1 = clock();
    std::clock_t beginN = clock();

    boptimize_sc<scopt::support_selection_t::NGREEDY, 4u, 4u>( scg, rps, &rst_p1 );
    scg = cleanup_scg( scg );
    std::clock_t end1 = clock();

    double aopt1 = scg.compute_area();
    double dopt1 = scg.compute_worst_delay();
    printf("[a]%6f ", aold );
    printf("-> %6f ", scg.compute_area() );
    printf("[d]%6f ", dold );
    printf("-> %6f ", scg.compute_worst_delay() );

    std::cout << std::endl;

    scopt::scg_network scg_copy=scg;

    int it = 0;
    while( scg.compute_area() < aold1 )
    {
      it++;
      aold1 = scg.compute_area();
      boptimize_sc<scopt::support_selection_t::NGREEDY, 4u, 4u>( scg, rps, &rst_p1 );
      scg = cleanup_scg( scg );
//      if( aold1 == scg.compute_area() )
//      {
//        boptimize_sc<scopt::support_selection_t::PIVOT, 4u, 4u>( scg, rps, &rst_p1 );
//        scg = cleanup_scg( scg );
//      }
      printf("[a]%6f ", aold );
      printf("-> %6f ", scg.compute_area() );
      printf("[d]%6f ", dold );
      printf("-> %6f ", scg.compute_worst_delay() );
      std::cout << std::endl;
    }

//    printf("{a}%6f ", aold );
//    printf("-> %6f ", scg.compute_area() );
//    printf("{d}%6f ", dold );
//    printf("-> %6f ", scg.compute_worst_delay() );

    std::clock_t endN = clock();
    double timeN = double(endN - beginN) / CLOCKS_PER_SEC;
    double time1 = double(end1 - begin1) / CLOCKS_PER_SEC;

    start=false;       

    double const d_opt = scg.compute_worst_delay();

    printf("[d]%6f ", dold );
    printf("-> %6f ", d_opt );
    std::cout << std::endl;

    auto const cecMP = benchmark == "hyp" ? true : abc_cec( scg, benchmark );
    if(!cecMP)
      printf("ERROR\n");
    std::cout << std::endl;

    exp( benchmark, aold, aopt1, scg.compute_area(), dold, dopt1, scg.compute_worst_delay(), time1, timeN, cecMP );


  }

  exp.save();
  exp.table();
  return 0;
}

// only asap
//|  benchmark |     a(map) |    a(opt1) |    a(optN) |    d(map) |   d(opt1) |   d(optN) | t(opt1) | t(optN) |  cec |
//|      adder |   16849.03 |   16835.27 |   16816.50 |    750.11 |    750.11 |    750.11 |    0.59 |    1.75 | true |
//|        bar |   23180.62 |   23180.62 |   23180.62 |    637.58 |    637.58 |    637.58 |    1.05 |    1.05 | true |
//|        div |  383396.66 |  382931.19 |  382585.94 |  32671.86 |  32671.86 |  32671.86 |   21.22 |  146.00 | true |
//|        hyp | 1540569.12 | 1538194.25 | 1536930.62 | 913722.88 | 913722.88 | 913722.88 |   31.31 |  526.34 | true |
//|       log2 |  277906.72 |  277388.78 |  276824.53 |   8922.27 |   8922.27 |   9081.71 |   14.89 |  161.80 | true |
//|        max |   30661.42 |   30478.74 |   30469.99 |   1535.95 |   1535.95 |   1535.95 |    0.71 |    3.65 | true |
//| multiplier |  252636.78 |  252449.12 |  252213.89 |   4512.04 |   4632.32 |   4558.32 |   10.05 |   70.82 | true |
//|        sin |   56218.87 |   56135.05 |   56082.51 |   4689.99 |   4689.99 |   4689.99 |    2.65 |   12.65 | true |
//|       sqrt |  212138.59 |  211858.48 |  211583.41 | 153025.02 | 153025.02 | 153024.05 |    5.05 |   55.09 | true |
//|     square |  170233.20 |  170039.30 |  169954.25 |   1828.73 |   1828.73 |   1828.73 |    4.52 |   18.34 | true |
//|    arbiter |   43046.72 |   43017.96 |   42999.20 |    655.88 |    655.88 |    655.88 |    0.99 |    2.95 | true |
//|      cavlc |    3779.19 |    3774.19 |    3774.19 |    578.83 |    578.83 |    578.83 |    0.48 |    0.96 | true |
//|       ctrl |     634.50 |     628.25 |     628.25 |    359.95 |    359.95 |    359.95 |    0.35 |    0.70 | true |
//|        dec |    2125.58 |    2125.58 |    2125.58 |    235.66 |    235.66 |    235.66 |    0.41 |    0.42 | true |
//|        i2c |    7154.19 |    7154.19 |    7154.19 |    474.63 |    474.63 |    474.63 |    0.41 |    0.41 | true |
//|  int2float |    1341.53 |    1341.53 |    1341.53 |    479.04 |    479.04 |    479.04 |    0.36 |    0.37 | true |
//|   mem_ctrl |  194325.36 |  191044.81 |  185621.89 |   1821.07 |   1983.36 |   1881.12 |    6.04 |  173.33 | true |
//|   priority |    8546.42 |    8542.67 |    8542.67 |   5488.67 |   5488.67 |   5488.67 |    0.45 |    0.90 | true |
//|     router |    1312.71 |    1311.46 |    1311.46 |    530.48 |    530.48 |    530.48 |    0.34 |    0.69 | true |
//|      voter |   92939.71 |   92929.71 |   92925.96 |   2500.19 |   2500.19 |   2500.19 |    2.60 |   11.12 | true |

//map=np.array([16849.03,23180.62,383396.66,1540569.12,277906.72,30661.42,252636.78,56218.87,212138.59,170233.20,43046.72,3779.19,634.50,2125.58,7154.19,1341.53,194325.36,8546.42,1312.71,92939.71])
//opt1=np.array([16835.27,23180.62,382931.19,1538194.25,277388.78,30478.74,252449.12,56135.05,211858.48,170039.30,43017.96,3774.19,628.25,2125.58,7154.19,1341.53,191044.81,8542.67,1311.46,92929.71])
//optN=np.array([16816.50,23180.62,382585.94,1536930.62,276824.53,30469.99,252213.89,56082.51,211583.41,169954.25,42999.20,3774.19,628.25,2125.58,7154.19,1341.53,185621.89,8542.67,1311.46,92925.96])


//|       c17 |    40.05 |    40.05 |    40.05 |  162.13 |  162.13 |  162.13 |    0.28 |    0.29 | true |
//|      c432 |  2178.60 |  2178.60 |  2178.60 |  995.92 |  995.92 |  995.92 |    0.31 |    0.31 | true |
//|      c499 |  4197.39 |  4182.39 |  4156.12 |  831.76 |  831.76 |  831.76 |    0.32 |    1.29 | true |
//|      c880 |  2650.53 |  2650.53 |  2650.53 |  735.43 |  735.43 |  735.43 |    0.30 |    0.31 | true |
//|     c1355 |  5039.42 |  4949.40 |  4916.90 |  865.26 |  865.26 |  865.26 |    0.32 |    1.28 | true |
//|     c1908 |  3067.22 |  3059.71 |  3059.71 |  859.19 |  859.19 |  859.19 |    0.31 |    0.62 | true |
//|     c2670 |  4667.89 |  4647.87 |  4630.35 |  786.08 |  786.08 |  786.08 |    0.32 |    0.96 | true |
//|     c3540 |  6140.56 |  6119.28 |  6119.28 | 1277.77 | 1277.77 | 1277.77 |    0.36 |    0.74 | true |
//|     c5315 |  9923.81 |  9902.55 |  9887.53 |  914.05 |  914.05 |  914.05 |    0.39 |    1.19 | true |
//|     c6288 | 26171.34 | 26137.55 | 26076.28 | 2426.04 | 2426.04 | 2426.04 |    0.76 |    2.96 | true |
//|     c7552 | 10607.22 | 10570.96 | 10559.70 |  970.33 |  970.33 |  970.33 |    0.42 |    1.93 | true |

//iscas_map=np.array([40.05,2178.60,4197.39,2650.53,5039.42,3067.22,4667.89,6140.56,9923.81,26171.34,10607.22])
//iscas_opt1=np.array([40.05,2178.60,4182.39,2650.53,4949.40,3059.71,4647.87,6119.28,9902.55,26137.55,10570.96])
//iscas_optN=np.array([40.05,2178.60,4156.12,2650.53,4916.90,3059.71,4630.35,6119.28,9887.53,26076.28,10559.70])


//|  benchmark &      a(map) &    a(opt1) &    a(optN) &    d(map) &   d(opt1) |   d(optN) | t(opt1) | t(optN) |  cec |
//|      adder &    16849.03 &   16835.27 &   16816.50 &    750.11 &    750.11 |    750.11 |    0.59 |    1.75 | true |
//|        bar &    23180.62 &   23180.62 &   23180.62 &    637.58 &    637.58 |    637.58 |    1.05 |    1.05 | true |
//|        div &   383396.66 &  382931.19 &  382585.94 &  32671.86 &  32671.86 |  32671.86 |   21.22 |  146.00 | true |
//|        hyp &  1540569.12 & 1538194.25 & 1536930.62 & 913722.88 & 913722.88 | 913722.88 |   31.31 |  526.34 | true |
//|       log2 &   277906.72 &  277388.78 &  276824.53 &   8922.27 &   8922.27 |   9081.71 |   14.89 |  161.80 | true |
//|        max &    30661.42 &   30478.74 &   30469.99 &   1535.95 &   1535.95 |   1535.95 |    0.71 |    3.65 | true |
//| multiplier &   252636.78 &  252449.12 &  252213.89 &   4512.04 &   4632.32 |   4558.32 |   10.05 |   70.82 | true |
//|        sin &    56218.87 &   56135.05 &   56082.51 &   4689.99 &   4689.99 |   4689.99 |    2.65 |   12.65 | true |
//|       sqrt &   212138.59 &  211858.48 &  211583.41 & 153025.02 & 153025.02 | 153024.05 |    5.05 |   55.09 | true |
//|     square &   170233.20 &  170039.30 &  169954.25 &   1828.73 &   1828.73 |   1828.73 |    4.52 |   18.34 | true |
//|    arbiter &    43046.72 &   43017.96 &   42999.20 &    655.88 &    655.88 |    655.88 |    0.99 |    2.95 | true |
//|      cavlc &     3779.19 &    3774.19 &    3774.19 &    578.83 &    578.83 |    578.83 |    0.48 |    0.96 | true |
//|       ctrl &      634.50 &     628.25 &     628.25 &    359.95 &    359.95 |    359.95 |    0.35 |    0.70 | true |
//|        dec &     2125.58 &    2125.58 &    2125.58 &    235.66 &    235.66 |    235.66 |    0.41 |    0.42 | true |
//|        i2c &     7154.19 &    7154.19 &    7154.19 &    474.63 &    474.63 |    474.63 |    0.41 |    0.41 | true |
//|  int2float &     1341.53 &    1341.53 &    1341.53 &    479.04 &    479.04 |    479.04 |    0.36 |    0.37 | true |
//|   mem_ctrl &   194325.36 &  191044.81 &  185621.89 &   1821.07 &   1983.36 |   1881.12 |    6.04 |  173.33 | true |
//|   priority &     8546.42 &    8542.67 &    8542.67 &   5488.67 &   5488.67 |   5488.67 |    0.45 |    0.90 | true |
//|     router &     1312.71 &    1311.46 &    1311.46 &    530.48 &    530.48 |    530.48 |    0.34 |    0.69 | true |
//|      voter &    92939.71 &   92929.71 &   92925.96 &   2500.19 &   2500.19 |   2500.19 |    2.60 |   11.12 | true |






//| benchmark |   a(map) |  a(opt1) |  a(optN) |  d(map) | d(opt1) | d(optN) | t(opt1) | t(optN) |  cec |
//|       c17 |    40.05 |    40.05 |    40.05 |  162.13 |  162.13 |  162.13 |    0.35 |    0.35 | true |
//|      c432 |  2184.85 |  2177.35 |  2177.35 | 1061.17 | 1061.17 | 1061.17 |    0.57 |    0.89 | true |
//|      c499 |  4144.75 |  4119.75 |  4083.46 |  841.26 |  841.26 |  841.26 |    1.08 |    1.95 | true |
//|      c880 |  2704.26 |  2701.76 |  2701.76 |  739.10 |  739.10 |  739.10 |    0.44 |    0.81 | true |
//|     c1355 |  5396.19 |  5351.15 |  5323.63 |  855.77 |  855.77 |  855.77 |    0.54 |    1.42 | true |
//|     c1908 |  3064.67 |  3064.67 |  3064.67 |  859.19 |  859.19 |  859.19 |    0.38 |    0.38 | true |
//|     c2670 |  4760.56 |  4760.56 |  4760.56 |  756.89 |  756.89 |  756.89 |    0.45 |    0.46 | true |
//|     c3540 |  5897.87 |  5875.34 |  5857.82 | 1303.61 | 1303.61 | 1303.61 |    0.70 |    1.53 | true |
//|     c5315 |  9269.28 |  9248.00 |  9248.00 |  939.53 |  939.53 |  939.53 |    0.62 |    1.13 | true |
//|     c6288 | 27492.70 | 27428.91 | 27355.09 | 2447.45 | 2447.45 | 2447.45 |    1.76 |    7.31 | true |
//|     c7552 | 12416.79 | 12346.71 | 12317.93 |  896.68 |  896.68 |  896.68 |    1.16 |    2.35 | true |
//1 best
//|       c17 |    40.05 |    40.05 |    40.05 |  162.13 |  162.13 |  162.13 |    0.35 |    0.35 | true |
//|      c432 |  2184.85 |  2177.35 |  2159.83 | 1061.17 | 1061.17 | 1061.17 |    0.60 |    1.33 | true |
//|      c499 |  4144.75 |  4121.00 |  4097.23 |  841.26 |  841.26 |  841.26 |    1.10 |    2.47 | true |
//|      c880 |  2704.26 |  2704.26 |  2704.26 |  739.10 |  739.10 |  739.10 |    0.39 |    0.39 | true |
//|     c1355 |  5396.19 |  5321.13 |  5318.63 |  855.77 |  855.77 |  855.77 |    0.85 |    1.68 | true |
//|     c1908 |  3064.67 |  3033.39 |  3033.39 |  859.19 |  859.19 |  859.19 |    0.53 |    0.87 | true |
//|     c2670 |  4760.56 |  4743.04 |  4743.04 |  756.89 |  756.89 |  756.89 |    0.51 |    0.85 | true |
//|     c3540 |  5897.87 |  5861.59 |  5849.07 | 1303.61 | 1303.61 | 1303.61 |    0.72 |    1.58 | true |
//|     c5315 |  9269.28 |  9232.99 |  9229.24 |  939.53 |  939.53 |  939.53 |    0.82 |    1.92 | true |
//|     c6288 | 27492.70 | 27448.91 | 27368.85 | 2447.45 | 2447.45 | 2447.45 |    1.65 |    7.74 | true |
//|     c7552 | 12416.79 | 12341.71 | 12325.44 |  896.68 |  896.68 |  896.68 |    1.25 |    2.50 | true |

// n=2
//|       c17 |    40.05 |    40.05 |    40.05 |  162.13 |  162.13 |  162.13 |    0.35 |    0.35 | true |
//|      c432 |  2184.85 |  2121.04 |  2113.53 | 1061.17 | 1061.17 | 1061.17 |    0.84 |    2.02 | true |
//|      c499 |  4144.75 |  4134.75 |  4112.24 |  841.26 |  841.26 |  841.26 |    0.69 |    1.81 | true |
//|      c880 |  2704.26 |  2698.00 |  2679.22 |  739.10 |  739.10 |  739.10 |    0.52 |    1.31 | true |
//|     c1355 |  5396.19 |  5223.53 |  5099.68 |  855.77 |  855.77 |  855.77 |    1.62 |    4.50 | true |
//|     c1908 |  3064.67 |  3044.65 |  3033.39 |  859.19 |  859.19 |  859.19 |    0.53 |    1.28 | true |
//|     c2670 |  4760.56 |  4720.52 |  4720.52 |  756.89 |  756.89 |  756.89 |    0.60 |    0.95 | true |
//|     c3540 |  5897.87 |  5851.57 |  5842.81 | 1303.61 | 1303.61 | 1303.61 |    1.07 |    2.05 | true |
//|     c5315 |  9269.28 |  9179.18 |  9129.15 |  939.53 |  939.53 |  939.53 |    1.44 |    5.45 | true |
//|     c6288 | 27492.70 | 26957.16 | 26314.09 | 2447.45 | 2447.45 | 2447.45 |    5.44 |   28.88 | true |
//|     c7552 | 12416.79 | 12321.69 | 12105.23 |  896.68 |  896.68 |  896.68 |    1.29 |    5.67 | true |

//  benchmark &     a(map) &   a(opt1) &   a(optN) &    d(map) &   d(opt1) &   d(optN) & t(opt1) & t(optN) \\
//      adder &   16849.03 &  16698.87 &  16638.81 &    750.11 &    750.11 &    750.11 &    0.77 &    3.01 \\
//        bar &   23180.62 &  23179.37 &  23179.37 &    637.58 &    637.58 &    637.58 &    1.36 &    2.77 \\
//        div &  376611.84 & 373741.41 & 368049.56 &  32615.78 &  32635.42 &  32665.47 &   24.17 &  100.06 \\
//       log2 &  274561.25 & 272531.75 & 269261.12 &   8940.15 &   8940.15 &   8940.15 &   18.14 &   74.02 \\
//        max &   30661.42 &  30542.56 &  30424.92 &   1535.95 &   1535.95 &   1535.95 &    1.12 &    4.08 \\
// multiplier &  252636.78 & 250519.61 & 247225.06 &   4512.04 &   4512.04 &   4512.04 &   15.78 &   64.30 \\
//        sin &   55348.14 &  55049.10 &  54487.31 &   4675.70 &   4675.70 &   4675.70 &    3.85 &   16.24 \\
//       sqrt &  212138.59 & 209938.80 & 207103.38 & 153025.02 & 153057.30 & 153057.16 &   10.85 &   45.78 \\
//     square &  173812.38 & 173000.41 & 171561.62 &   1822.10 &   1822.10 &   1822.10 &    7.56 &   29.70 \\
//    arbiter &   43046.72 &  42922.87 &  42623.86 &    655.88 &    655.88 &    655.88 &    1.35 &    5.94 \\
//      cavlc &    3779.19 &   3732.90 &   3729.15 &    578.83 &    578.83 &    578.83 &    0.50 &    1.55 \\
//       ctrl &     634.50 &    634.50 &    634.50 &    359.95 &    359.95 &    359.95 &    0.34 &    0.41 \\
//        dec &    2125.58 &   2125.58 &   2125.58 &    235.66 &    235.66 &    235.66 &    0.47 &    0.53 \\
//        i2c &    7154.19 &   7145.43 &   7135.42 &    474.63 &    474.63 &    474.63 &    0.43 &    1.40 \\
//  int2float &    1341.53 &   1341.53 &   1341.53 &    479.04 &    479.04 &    479.04 &    0.37 &    0.43 \\
//   mem_ctrl &  195729.22 & 192972.67 & 188172.91 &   1847.50 &   1847.50 &   1847.50 &    9.50 &   37.22 \\
//   priority &    8546.42 &   8535.16 &   8522.65 &   5488.67 &   5488.67 &   5488.67 &    0.54 &    2.14 \\
//     router &    1312.71 &   1303.95 &   1293.94 &    530.48 &    530.48 &    530.48 &    0.34 &    1.08 \\
//      voter &   92939.71 &  92563.13 &  91941.31 &   2500.19 &   2500.19 &   2500.19 &    3.93 &   15.12 \\
//
//amap=np.array([  16849.03,  23180.62, 376611.84, 274561.25,  30661.42, 252636.78,  55348.14, 212138.59, 173812.38,  43046.72,   3779.19,    634.50,   2125.58,   7154.19,   1341.53, 195729.22,   8546.42,   1312.71,  92939.71])
//aopt1=np.array([16698.87,23179.37,373741.41,272531.75,30542.56,250519.61,55049.10,209938.80,173000.41,42922.87,3732.90,634.50,2125.58,7145.43,1341.53,192972.67,8535.16,1303.95,92563.13])
//aopt4=([16638.81,23179.37,368049.56,269261.12,30424.92,247225.06,54487.31,207103.38,171561.62,42623.86,3729.15,634.50,2125.58,7135.42,1341.53,188172.91,8522.65,1293.94,91941.31])


//|  benchmark |    a(map) |   a(opt1) |   a(optN) |    d(map) |   d(opt1) |   d(optN) | t(opt1) | t(optN) |  cec |
//      adder & $  12624.30 $ & $  12417.85 $ & $  12335.26 $ & $    775.40 $ & $    775.40 $ & $    775.40 $ & $    0.67 $ & $    3.38 $\\
//        bar & $  23026.63 $ & $  23015.37 $ & $  22995.35 $ & $    637.58 $ & $    637.58 $ & $    637.58 $ & $    1.47 $ & $    4.60 $\\
//        div & $ 374409.62 $ & $ 371596.69 $ & $ 362587.78 $ & $  32674.05 $ & $  32735.75 $ & $  32765.45 $ & $   26.72 $ & $  178.02 $\\
//       log2 & $ 276389.19 $ & $ 274345.94 $ & $ 268238.75 $ & $   9028.77 $ & $   9056.90 $ & $   9105.03 $ & $   40.84 $ & $  234.71 $\\
//        max & $  35189.79 $ & $  34955.82 $ & $  34669.24 $ & $   1566.76 $ & $   1566.76 $ & $   1566.76 $ & $    1.18 $ & $    6.94 $\\
// multiplier & $ 242398.56 $ & $ 240300.23 $ & $ 231669.69 $ & $   4640.00 $ & $   4649.30 $ & $   4718.32 $ & $   17.64 $ & $  113.60 $\\
//        sin & $  53765.88 $ & $  53289.17 $ & $  52038.04 $ & $   4700.92 $ & $   4700.92 $ & $   4700.92 $ & $    9.54 $ & $   53.87 $\\
//       sqrt & $ 206618.39 $ & $ 203979.34 $ & $ 199147.98 $ & $ 146805.59 $ & $ 146805.59 $ & $ 146805.59 $ & $   12.76 $ & $   73.24 $\\
//     square & $ 163537.81 $ & $ 162411.86 $ & $ 160061.06 $ & $   1977.66 $ & $   1977.66 $ & $   1977.66 $ & $    6.94 $ & $   39.23 $\\
//    arbiter & $  44021.13 $ & $  43944.82 $ & $  43700.88 $ & $    653.73 $ & $    653.73 $ & $    653.73 $ & $    1.73 $ & $   10.58 $\\
//      cavlc & $   3767.86 $ & $   3725.32 $ & $   3725.32 $ & $    582.35 $ & $    582.35 $ & $    582.35 $ & $    0.52 $ & $    1.10 $\\
//       ctrl & $    673.26 $ & $    673.26 $ & $    673.26 $ & $    360.43 $ & $    360.43 $ & $    360.43 $ & $    0.38 $ & $    0.44 $\\
//        dec & $   2125.58 $ & $   2125.58 $ & $   2125.58 $ & $    235.66 $ & $    235.66 $ & $    235.66 $ & $    0.50 $ & $    0.57 $\\
//        i2c & $   7168.03 $ & $   7145.51 $ & $   7117.99 $ & $    455.45 $ & $    455.45 $ & $    455.45 $ & $    0.57 $ & $    1.65 $\\
//  int2float & $   1377.85 $ & $   1355.33 $ & $   1355.33 $ & $    476.85 $ & $    476.85 $ & $    476.85 $ & $    0.37 $ & $    0.82 $\\
//  mem\_ctrl & $ 188258.36 $ & $ 185730.91 $ & $ 178542.45 $ & $   1868.68 $ & $   1893.36 $ & $   1915.00 $ & $    9.82 $ & $   66.41 $\\
//   priority & $   3963.30 $ & $   3937.02 $ & $   3915.75 $ & $    603.74 $ & $    603.74 $ & $    603.74 $ & $    0.40 $ & $    1.62 $\\
//     router & $   1165.01 $ & $   1165.01 $ & $   1165.01 $ & $    552.01 $ & $    552.01 $ & $    552.01 $ & $    0.35 $ & $    0.42 $\\
//      voter & $  92372.17 $ & $  91896.77 $ & $  90852.19 $ & $   2491.19 $ & $   2491.19 $ & $   2491.19 $ & $    5.01 $ & $   30.22 $\\
//
//amap=np.array([12624.30,23026.63,374409.62,276389.19,35189.79,242398.56,53765.88,206618.39,163537.81,44021.13,3767.86,673.26,2125.58,7168.03,1377.85,188258.36,3963.30,1165.01,92372.17])
//aopt1=np.array([12417.85,23015.37,371596.69,274345.94,34955.82,240300.23,53289.17,203979.34,162411.86,43944.82,3725.32,673.26,2125.58,7145.51,1355.33,185730.91,3937.02,1165.01,91896.77])
//aopt5=np.array([12335.26,22995.35,362587.78,268238.75,34669.24,231669.69,52038.04,199147.98,160061.06,43700.88,3725.32,673.26,2125.58,7117.99,1355.33,178542.45,3915.75,1165.01,90852.19])

//dmap=np.array([775.40,637.58,32674.05,9028.77,1566.76,4640.00,4700.92,146805.59,1977.66,653.73,582.35,360.43,235.66,455.45,476.85,1868.68,603.74,552.01,2491.19])
//dopt1=np.array([775.40,637.58,32735.75,9056.90,1566.76,4649.30,4700.92,146805.59,1977.66,653.73,582.35,360.43,235.66,455.45,476.85,1893.36,603.74,552.01,2491.19])
//dopt5=np.array([775.40,637.58,32765.45,9105.03,1566.76,4718.32,4700.92,146805.59,1977.66,653.73,582.35,360.43,235.66,455.45,476.85,1915.00,603.74,552.01,2491.19])


//       c17 & $     40.05 $ & $    40.05 $ & $    40.05 $ & $  162.13 $ & $  162.13 $ & $  162.13 $ & $  0.34 $ & $  0.40 $\\
//      c432 & $   2471.42 $ & $  2462.66 $ & $  2395.09 $ & $ 1053.92 $ & $ 1053.92 $ & $ 1053.92 $ & $  0.46 $ & $  3.24 $\\
//      c499 & $   4586.66 $ & $  4561.63 $ & $  4541.61 $ & $  813.68 $ & $  813.68 $ & $  813.68 $ & $  0.74 $ & $  2.37 $\\
//      c880 & $   2890.75 $ & $  2871.97 $ & $  2871.97 $ & $  709.17 $ & $  709.17 $ & $  709.17 $ & $  0.57 $ & $  1.03 $\\
//     c1355 & $   4514.17 $ & $  4479.13 $ & $  4479.13 $ & $  837.04 $ & $  837.04 $ & $  837.04 $ & $  0.75 $ & $  1.24 $\\
//     c1908 & $   3079.76 $ & $  3071.00 $ & $  3071.00 $ & $  871.81 $ & $  871.81 $ & $  871.81 $ & $  0.52 $ & $  0.98 $\\
//     c2670 & $   4461.38 $ & $  4446.37 $ & $  4431.35 $ & $  769.98 $ & $  769.98 $ & $  769.98 $ & $  0.59 $ & $  2.29 $\\
//     c3540 & $   6966.47 $ & $  6925.17 $ & $  6915.16 $ & $ 1266.27 $ & $ 1266.27 $ & $ 1266.27 $ & $  0.89 $ & $  2.80 $\\
//     c5315 & $   9739.83 $ & $  9671.01 $ & $  9592.20 $ & $  948.36 $ & $  948.36 $ & $  948.36 $ & $  1.31 $ & $  4.57 $\\
//     c6288 & $  26484.01 $ & $ 25915.92 $ & $ 25354.16 $ & $ 2390.41 $ & $ 2390.41 $ & $ 2390.41 $ & $  5.89 $ & $ 19.17 $\\
//     c7552 & $  11274.19 $ & $ 11136.55 $ & $ 11060.23 $ & $  991.88 $ & $  991.88 $ & $  991.88 $ & $  1.85 $ & $  4.32 $\\