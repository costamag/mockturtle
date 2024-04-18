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

  experiment<std::string, double, double, double, double, double, double, double, double, double, double, double, double, uint32_t, bool> exp( "SCOPT", "benchmark", "a(map)", "a(opt1)", "a(optN)", "da(opt1)", "da(optN)", "d(map)", "d(opt1)", "d(optN)", "dd(opt1)", "dd(optN)", "t(opt1)", "t(optN)", "n(iters)", "cec");

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

  double N{0};
  double rarea1{0};
  double rareaN{0};
  double rdept1{0};
  double rdeptN{0};

  for ( auto const& benchmark : all_benchmarks( iwls ) )
  {
    if( benchmark == "hyp" ) continue;
    N+=1;

    fmt::print( "[i] processing {}\n", benchmark );

    bool start=true;
    bool close=false;

    double darea1;
    double dareaN;
    double ddept1;
    double ddeptN;

    aig_network aig;//experiments::benchmark_path(benchmark)
    if ( lorina::read_aiger( experiments::benchmark_path(benchmark), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    if( aig.num_gates() > 300000 )  continue;

    auto aaig_old = aig.num_gates()+1;
    auto aaig_new = aig.num_gates();
    depth_view<aig_network> daig{aig};
    uint32_t daig_old = daig.depth()+1;
    uint32_t daig_new = daig.depth();

    //if( benchmark != "hyp" )
    //{
    //aig = abc_opto( aig, benchmark, "rw; rs;" );
    //aig = abc_opto( aig, benchmark, "rw; rs;" );
    //aig = abc_opto( aig, benchmark, "rw; rs;" );
    //aig = abc_opto( aig, benchmark, "rw; rs;" );

    int it2{0};
    while( it2++<3 )
    {
      aig = abc_opto( aig, benchmark, "resyn2rs" );
      aig = cleanup_dangling( aig );
      aig = abc_if( aig, benchmark );
      aig = cleanup_dangling( aig );
      depth_view<aig_network> daig{aig};
      daig_old = daig_new;
      daig_new = daig.depth();
      printf("%d\n", daig_new);
    }
//
  //  write_aiger(aig, "ERR.aig");
   // const auto cec = benchmark == "hyp" ? true : abc_cec( aig, benchmark );
   // assert( cec && "[e] not equivalent" );

    scopt::emap2_params ps;
    //ps.cut_enumeration_ps.minimize_truth_table = true;
    //ps.cut_enumeration_ps.cut_limit = 24;
    ps.area_oriented_mapping = false;
    scopt::emap2_stats st;

    printf("map..\n");

    scopt::scg_network scg = emap2_klut( aig, tech_lib, ps, &st );
    scg = cleanup_scg( scg );

    if( scg.compute_area() > 400000 )  continue;


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
    rps.max_pis = 12;
    rps.verbose = false;
    rps.max_divisors = 128;

    boptimizer_stats rst_p1;
    rps.use_delay_constraints = true;

    double aold1 = scg.compute_area();
    double dold1 = scg.compute_worst_delay();


    std::clock_t begin1 = clock();
    std::clock_t beginN = clock();

    scopt::scg_network scg_copy1=scg;


    boptimize_sc<support_selection_t::EX1, 4u, 4u>( scg, rps, &rst_p1 );
    scg = cleanup_scg( scg );
    if( scg.compute_worst_delay() > dold1 )
    {
      scg=scg_copy1;
    }

    std::clock_t end1 = clock();

    double aopt1 = scg.compute_area();
    double dopt1 = scg.compute_worst_delay();
    printf("[a]%6f ", aold );
    printf("-> %6f ", scg.compute_area() );
    printf("[d]%6f ", dold );
    printf("-> %6f ", scg.compute_worst_delay() );


    darea1=(scg.compute_area()-aold)/(aold);
    ddept1=(scg.compute_worst_delay()-dold)/(dold);

    rarea1 = rarea1*(N-1)/N+darea1/N;
    rdept1 = rdept1*(N-1)/N+ddept1/N;

    std::cout << std::endl;

    scopt::scg_network scg_copy=scg;

    double aoldN = scg.compute_area()+1;
    double doldN = scg.compute_worst_delay()+1;

    int it = 1;
    //while( scg.compute_area() < aoldN && scg.compute_worst_delay() <= dold && it < 2u )
    {
      it++;
      aoldN = scg.compute_area();
      boptimize_sc<support_selection_t::EX1, 4u, 4u>( scg, rps, &rst_p1 );
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
      if( scg.compute_worst_delay() <= dold )
        scg_copy=scg;
    }
    scg=scg_copy;

    dareaN=(scg.compute_area()-aold)/(aold);
    ddeptN=(scg.compute_worst_delay()-dold)/(dold);

    rareaN = rareaN*(N-1)/N+dareaN/N;
    rdeptN = rdeptN*(N-1)/N+ddeptN/N;

    printf(" a1 %6f  d1 %6f\n", 100*darea1, 100*ddept1 );
    printf(" aN %6f  dN %6f\n", 100*dareaN, 100*ddeptN );
    printf("<a1>%6f <d1>%6f\n", 100*rarea1, 100*rdept1 );
    printf("<aN>%6f <dN>%6f\n", 100*rareaN, 100*rdeptN );
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

    exp( benchmark, aold, aopt1, scg.compute_area(), 100*darea1, 100*dareaN, dold, dopt1, scg.compute_worst_delay(), 100*ddept1, 100*ddeptN, time1, timeN, it, cecMP );


  }

  exp.save();
  exp.table();
  return 0;
}


//|      adder |  12624.30 |  12563.01 |  12539.24 |    775.40 |    775.40 |    775.40 |    0.67 |    2.73 | true |
//|        bar |  23026.63 |  23026.63 |  23026.63 |    637.58 |    637.58 |    637.58 |    1.35 |    1.42 | true |
//|        div | 386844.47 | 385968.72 | 384846.47 |  32282.48 |  32282.48 |  32282.46 |   43.19 |  253.33 | true |
//|       log2 | 279122.59 | 278169.16 | 277190.72 |   8958.88 |   8958.88 |   8958.88 |   63.23 |  387.85 | true |
//|        max |  35189.79 |  34834.45 |  34710.59 |   1566.76 |   1566.76 |   1566.76 |    1.08 |    5.69 | true |
//| multiplier | 240776.30 | 239845.44 | 238917.05 |   4651.26 |   4651.26 |   4651.23 |   18.26 |  110.94 | true |
//|        sin |  55466.51 |  55271.30 |  55044.82 |   4755.72 |   4755.72 |   4755.72 |   16.19 |   85.28 | true |
//|       sqrt | 206527.38 | 206372.27 | 206150.81 | 147180.48 | 147180.48 | 147180.48 |   12.99 |   66.99 | true |
//|     square | 161457.44 | 160964.58 | 160240.30 |   1959.73 |   1959.73 |   1959.73 |    7.93 |   44.95 | true |
//|    arbiter |  44021.13 |  43969.84 |  43866.03 |    653.73 |    653.73 |    653.73 |    1.70 |   10.01 | true |
//|      cavlc |   3853.12 |   3853.12 |   3853.12 |    558.22 |    558.22 |    558.22 |    0.50 |    0.57 | true |
//|       ctrl |    710.83 |    710.83 |    710.83 |    364.44 |    364.44 |    364.44 |    0.35 |    0.42 | true |
//|        dec |   2125.58 |   2125.58 |   2125.58 |    235.66 |    235.66 |    235.66 |    0.45 |    0.51 | true |
//|        i2c |   7004.07 |   7002.82 |   6999.07 |    456.90 |    456.90 |    456.90 |    0.48 |    1.57 | true |
//|  int2float |   1336.51 |   1332.76 |   1332.76 |    477.25 |    477.25 |    477.25 |    0.42 |    0.87 | true |
//|   mem_ctrl | 188560.59 | 185935.48 | 179511.31 |   1883.51 |   1883.51 |   1883.51 |   11.68 |   67.02 | true |
//|   priority |   3963.30 |   3952.04 |   3952.04 |    603.74 |    603.74 |    603.74 |    0.39 |    0.87 | true |
//|     router |   1165.01 |   1165.01 |   1165.01 |    552.01 |    552.01 |    552.01 |    0.35 |    0.42 | true |
//|      voter |  90069.11 |  89812.70 |  89336.12 |   2586.57 |   2586.57 |   2586.57 |    7.12 |   39.86 | true |



//|  benchmark |    a(map) |   a(opt1) |   a(optN) |    d(map) |   d(opt1) |   d(optN) | t(opt1) | t(optN) |  cec |
//      adder & $  12624.30 $ & $  12542.97 $ & $  12432.85 $ & $    775.40 $ & $    775.40 $ & $    775.40 $ & $    0.86 $ & $    4.49 $\\
//        bar & $  23026.63 $ & $  23012.86 $ & $  23012.86 $ & $    637.58 $ & $    637.58 $ & $    637.58 $ & $    1.42 $ & $    2.89 $\\
//        div & $ 380682.81 $ & $ 377828.59 $ & $ 371708.44 $ & $  32348.71 $ & $  32355.10 $ & $  32354.92 $ & $   39.03 $ & $  209.91 $\\
//       log2 & $ 274103.12 $ & $ 272010.88 $ & $ 267996.66 $ & $   8956.95 $ & $   8956.95 $ & $   8956.95 $ & $   41.53 $ & $  232.83 $\\
//        max & $  35189.79 $ & $  34925.75 $ & $  34730.55 $ & $   1566.76 $ & $   1566.76 $ & $   1566.76 $ & $    1.13 $ & $    6.25 $\\
// multiplier & $ 242398.56 $ & $ 240603.00 $ & $ 236612.64 $ & $   4640.00 $ & $   4640.00 $ & $   4640.00 $ & $   17.56 $ & $   94.81 $\\
//        sin & $  53272.17 $ & $  52781.70 $ & $  51714.34 $ & $   4743.26 $ & $   4743.26 $ & $   4742.63 $ & $    7.90 $ & $   48.11 $\\
//       sqrt & $ 206527.38 $ & $ 204353.75 $ & $ 200419.62 $ & $ 147180.48 $ & $ 147180.48 $ & $ 147180.48 $ & $   15.54 $ & $   76.53 $\\
//     square & $ 163537.81 $ & $ 162794.62 $ & $ 161238.19 $ & $   1977.66 $ & $   1977.66 $ & $   1977.66 $ & $    7.67 $ & $   40.94 $\\
//    arbiter & $  44021.13 $ & $  43976.09 $ & $  43717.11 $ & $    653.73 $ & $    653.73 $ & $    653.73 $ & $    1.46 $ & $    9.26 $\\
//      cavlc & $   3853.12 $ & $   3820.59 $ & $   3790.57 $ & $    558.22 $ & $    558.22 $ & $    558.22 $ & $    0.58 $ & $    1.70 $\\
//       ctrl & $    710.83 $ & $    699.57 $ & $    699.57 $ & $    364.44 $ & $    364.44 $ & $    364.44 $ & $    0.35 $ & $    0.78 $\\
//        dec & $   2125.58 $ & $   2125.58 $ & $   2125.58 $ & $    235.66 $ & $    235.66 $ & $    235.66 $ & $    0.45 $ & $    0.51 $\\
//        i2c & $   7004.07 $ & $   7001.57 $ & $   7001.57 $ & $    456.90 $ & $    456.90 $ & $    456.90 $ & $    0.46 $ & $    0.99 $\\
//  int2float & $   1336.51 $ & $   1306.48 $ & $   1300.23 $ & $    477.25 $ & $    477.25 $ & $    477.25 $ & $    0.44 $ & $    1.35 $\\
//  mem\_ctrl & $ 191124.98 $ & $ 188287.05 $ & $ 179819.34 $ & $   1914.17 $ & $   1914.17 $ & $   1914.17 $ & $   15.46 $ & $   75.59 $\\
//   priority & $   3963.30 $ & $   3929.52 $ & $   3920.76 $ & $    603.74 $ & $    603.74 $ & $    603.74 $ & $    0.43 $ & $    1.30 $\\
//     router & $   1165.01 $ & $   1165.01 $ & $   1165.01 $ & $    552.01 $ & $    552.01 $ & $    552.01 $ & $    0.35 $ & $    0.42 $\\
//      voter & $  96837.03 $ & $  96457.90 $ & $  95976.22 $ & $   2534.20 $ & $   2534.20 $ & $   2534.20 $ & $    6.07 $ & $   35.48 $\\

//       ac97_ctrl & $  68035.37 $ & $  67490.99 $ & $  65726.26 $ & $  449.04 $ & $  449.04 $ & $  449.04 $ & $    1.54 $ & $    9.74 \\
//        aes_core & $ 138505.59 $ & $ 137873.64 $ & $ 136816.06 $ & $  964.26 $ & $  964.26 $ & $  964.26 $ & $    8.90 $ & $   43.08 \\
//        des_area & $  31261.40 $ & $  31122.49 $ & $  30931.04 $ & $ 1261.47 $ & $ 1261.47 $ & $ 1261.47 $ & $    2.25 $ & $   10.95 \\
//             DMA & $ 136259.86 $ & $ 135959.55 $ & $ 130898.91 $ & $  841.34 $ & $  841.34 $ & $ 1144.94 $ & $    5.80 $ & $   35.32 \\
//             DSP & $ 222689.50 $ & $ 221390.69 $ & $ 218960.80 $ & $ 1848.97 $ & $ 1848.97 $ & $ 1848.97 $ & $    9.76 $ & $   46.11 \\
//        ethernet & $ 342269.50 $ & $ 341303.66 $ & $ 339139.09 $ & $  767.53 $ & $  767.53 $ & $  767.53 $ & $   38.70 $ & $  208.20 \\
//      iwls05_i2c & $   6464.83 $ & $   6408.51 $ & $   6399.74 $ & $  453.34 $ & $  453.34 $ & $  453.34 $ & $    0.45 $ & $    1.47 \\
// iwls05_mem_ctrl & $  44006.05 $ & $  43799.61 $ & $  43533.09 $ & $  965.52 $ & $  965.52 $ & $  965.52 $ & $    1.57 $ & $    7.18 \\
//    pci_bridge32 & $ 108588.46 $ & $ 107849.12 $ & $ 106132.57 $ & $  854.12 $ & $  854.12 $ & $  854.12 $ & $    4.17 $ & $   22.09 \\
//            RISC & $ 345506.94 $ & $ 344341.97 $ & $ 341546.38 $ & $ 1858.51 $ & $ 1858.51 $ & $ 1858.51 $ & $   13.91 $ & $   67.32 \\
//            sasc & $   4496.33 $ & $   4463.79 $ & $   4435.00 $ & $  375.74 $ & $  375.74 $ & $  375.74 $ & $    0.40 $ & $    1.40 \\
//      simple_spi & $   4955.49 $ & $   4936.71 $ & $   4936.71 $ & $  475.81 $ & $  475.81 $ & $  475.81 $ & $    0.49 $ & $    0.98 \\
//             spi & $  24766.01 $ & $  24553.30 $ & $  24150.46 $ & $  931.73 $ & $  931.73 $ & $  931.73 $ & $    1.48 $ & $    6.44 \\
//          ss_pcm & $   3208.64 $ & $   3196.13 $ & $   3156.09 $ & $  315.55 $ & $  315.55 $ & $  315.55 $ & $    0.37 $ & $    1.61 \\
//      systemcaes & $  57760.39 $ & $  57433.80 $ & $  56774.44 $ & $ 1450.12 $ & $ 1450.12 $ & $ 1450.07 $ & $    3.04 $ & $   13.53 \\
//      systemcdes & $  20267.15 $ & $  20069.48 $ & $  19946.87 $ & $ 1086.82 $ & $ 1086.82 $ & $ 1086.82 $ & $    1.10 $ & $    4.17 \\
//            tv80 & $  42871.13 $ & $  42466.99 $ & $  41927.76 $ & $ 1607.90 $ & $ 1607.90 $ & $ 1607.90 $ & $    2.48 $ & $   11.89 \\
//       usb_funct & $  79712.68 $ & $  79468.62 $ & $  78919.25 $ & $  873.86 $ & $  873.86 $ & $  873.86 $ & $    2.48 $ & $   12.68 \\
//         usb_phy & $   2749.13 $ & $   2737.87 $ & $   2737.87 $ & $  370.37 $ & $  370.37 $ & $  370.37 $ & $    0.37 $ & $    0.80 \\
//       wb_conmax & $ 180592.08 $ & $ 180092.83 $ & $ 179144.38 $ & $  891.71 $ & $  891.71 $ & $  891.71 $ & $    7.33 $ & $   36.94 \\
//
//
//amap=np.array([12624.30,23026.63,380682.81,274103.12,35189.79,242398.56,53272.17,206527.38,163537.81,44021.13,3853.12,710.83,2125.58,7004.07,1336.51,191124.98,3963.30,1165.01,96837.03,  68035.37, 138505.59,  31261.40, 136259.86, 222689.50, 342269.50,   6464.83,  44006.05, 108588.46, 345506.94,   4496.33,   4955.49,  24766.01,   3208.64,  57760.39,  20267.15,  42871.13,  79712.68,   2749.13, 180592.08])
//aopt1=np.array([12542.97,23012.86,377828.59,272010.88,34925.75,240603.00,52781.70,204353.75,162794.62,43976.09,3820.59,699.57,2125.58,7001.57,1306.48,188287.05,3929.52,1165.01,96457.90,  67490.99, 137873.64,  31122.49, 135959.55, 221390.69, 341303.66,   6408.51,  43799.61, 107849.12, 344341.97,   4463.79,   4936.71,  24553.30,   3196.13,  57433.80,  20069.48,  42466.99,  79468.62,   2737.87, 180092.83])
//aoptN=np.array([12432.85,23012.86,371708.44,267996.66,34730.55,236612.64,51714.34,200419.62,161238.19,43717.11,3790.57,699.57,2125.58,7001.57,1300.23,179819.34,3920.76,1165.01,95976.22,  65726.26, 136816.06,  30931.04, 130898.91, 218960.80, 339139.09,   6399.74,  43533.09, 106132.57, 341546.38,   4435.00,   4936.71,  24150.46,   3156.09,  56774.44,  19946.87,  41927.76,  78919.25,   2737.87, 179144.38])

