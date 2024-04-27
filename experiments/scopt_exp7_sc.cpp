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

  experiment<std::string, double, double, double, double, double, double, double, double, double, double, double, double, bool> exp( "SCOPTA", "benchmark", "a(map)", "a(opt1)", "a(optN)", "da(opt1)", "da(optN)", "d(map)", "d(opt1)", "d(optN)", "dd(opt1)", "dd(optN)", "t(opt1)", "t(optN)", "cec");

  fmt::print( "[i] processing technology library\n" );


  /* library to map to technology */
  std::vector<gate> gates;
  std::ifstream in( cell_libraries_path( "asap7" ) );//sky130

  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  tech_library_params tps;
  tech_library<5, classification_type::np_configurations> tech_lib( gates, tps );


  double N{1};
  double rarea1{0};
  double rareaN{0};
  double rdept1{0};
  double rdeptN{0};
  for ( auto const& benchmark : all_benchmarks( epfl | iwls ) )
  {


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

    if( aig.num_gates()>300000 || benchmark == "hyp"  ) continue;

    if( benchmark != "hyp" )
    {

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
    ps.required_time = std::numeric_limits<float>::max();
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
    rps.max_trials = 1;
    rps.max_pis = 16;
    rps.verbose = false;
    rps.use_delay_constraints = true;
    rps.max_divisors = 128u;

    boptimizer_stats rst_p1;

    double aold1 = scg.compute_area();
    double dold1 = scg.compute_worst_delay();


    std::clock_t begin1 = clock();
    std::clock_t beginN = clock();

    boptimize_sc<EX2, 4u, 4u>( scg, rps, &rst_p1 );
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

  //  rps.max_trials = 10;
  //  rps.max_pis = 8;
  //  rps.max_divisors = 32;256u
    double aold_N = scg.compute_area()+1;

    double time_now=0;
    while( time_now < 300 && aold_N > scg.compute_area() )
    {
      aold_N = scg.compute_area();
      it++;
      aold1 = scg.compute_area();
      boptimize_sc<EX2, 4u, 4u>( scg, rps, &rst_p1 );
      scg = cleanup_scg( scg );

      
      printf("[a]%6f ", aold );
      printf("-> %6f ", scg.compute_area() );
      printf("[d]%6f ", dold );
      printf("-> %6f ", scg.compute_worst_delay() );

      std::clock_t end_now = clock();
      time_now = double(end_now - beginN) / CLOCKS_PER_SEC;
      if( (aold_N-scg.compute_area() )/aold_N < 0.01 )
      {
        rps.max_pis = std::min( 16u, rps.max_pis + 1 );
        rps.max_divisors = std::min( 300u, 2*rps.max_divisors );
      }

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


    rarea1 = rarea1*(N-1)/N+(aopt1-aold)/(N*aold);
    rdept1 = rdept1*(N-1)/N+(dopt1-dold)/(N*dold);

    rareaN = rareaN*(N-1)/N+(scg.compute_area()-aold)/(N*aold);
    rdeptN = rdeptN*(N-1)/N+(scg.compute_worst_delay()-dold)/(N*dold);

    double dA1 = (aopt1-aold)/(aold);
    double dAN = (scg.compute_area()-aold)/(aold);

    double dD1 = (dopt1-dold)/(dold);
    double dDN = (scg.compute_worst_delay()-dold)/(dold);

    printf(" n1 =%f  aN =%f\n", dA1, dAN );
    printf("<a1>=%f <d1>=%f\n", rarea1, rdept1 );
    printf("<aN>=%f <dN>=%f\n", rareaN, rdeptN );
    std::cout << std::endl;
    N+=1;

    exp( benchmark, aold, aopt1, scg.compute_area(), 100*dA1, 100*dAN, dold, dopt1, scg.compute_worst_delay(), 100*dD1, 100*dDN, time1, timeN, cecMP );

  }

  exp.save();
  exp.table();
  return 0;
}

//|       benchmark |    a(map) |   a(opt1) |   a(optN) |    d(map) |   d(opt1) |   d(optN) | t(opt1) | t(optN) |  cec |
//|           adder |   4942.62,    4942.62,    4942.62, |  17141.27 |  17141.27 |  17141.27 |    0.44 |    0.51 | true |
//|             bar |  13358.60,   13326.06,   13293.52, |   1563.40 |   1563.40 |   1563.40 |    0.64 |    2.50 | true |
//|             div | 217558.84,  217128.33,  212440.59, | 312464.03 | 313519.12 | 326271.31 |   11.47 |  314.48 | true |
//|            log2 | 150191.22,  148210.02,  147154.97, |  30475.64 |  30129.36 |  30402.94 |   34.31 |  374.31 | true |
//|             max |  16234.48,   16234.48,   16234.48, |  19885.02 |  19885.02 |  19885.02 |    0.54 |    0.62 | true |
//|      multiplier | 126844.97,  126710.99,  126539.47, |  18577.35 |  18577.35 |  18734.22 |    7.04 |   42.47 | true |
//|             sin |  26194.90,   26077.25,   25841.99, |  15116.20 |  15116.20 |  15282.53 |    5.00 |   34.14 | true |
//|            sqrt |  99702.40,   99646.08,   99502.21, | 374553.91 | 374848.12 | 374364.56 |    4.38 |   58.40 | true |
//|          square |  89569.38,   89379.17,   88700.87, |  17783.07 |  18091.07 |  18329.58 |    3.54 |   43.76 | true |
//|         arbiter |  63432.36,   63427.36,   63411.11, |   6484.11 |   6511.00 |   6554.34 |    1.43 |    5.84 | true |
//|           cavlc |   2949.64,    2945.89,    2945.89, |   1537.65 |   1607.04 |   1607.04 |    0.46 |    0.95 | true |
//|            ctrl |    569.42,     566.92,     566.92, |    830.05 |    830.05 |    830.05 |    0.29 |    0.64 | true |
//|             dec |   2023.04,    2023.04,    2023.04, |    484.51 |    484.51 |    484.51 |    0.38 |    0.44 | true |
//|             i2c |   6270.94,    6265.94,    6265.94, |   1653.97 |   1653.97 |   1653.97 |    0.42 |    0.87 | true |
//|       int2float |   1101.22,    1101.22,    1101.22, |   1262.34 |   1262.34 |   1262.34 |    0.32 |    0.39 | true |
//|        mem_ctrl | 191598.19,  184993.12,  151548.97, |   9352.22 |   9445.59 |   9177.44 |    9.65 |  410.90 | true |
//|        priority |   2713.05,    2701.79,    2701.79, |   4819.59 |   4819.59 |   4819.59 |    0.32 |    0.69 | true |
//|          router |    868.49,     868.49,     868.49, |   1665.39 |   1665.39 |   1665.39 |    0.28 |    0.35 | true |
//|           voter |  51364.55,   51294.45,   51045.39, |   4769.73 |   4769.73 |   4924.25 |    2.30 |   22.86 | true |
//|       ac97_ctrl |  53463.37,   53382.06,   53216.85, |   1043.31 |   1043.31 |   1153.01 |    1.07 |    5.70 | true |
//|        aes_core |  95644.55,   95401.81,   94262.04, |   2499.02 |   2499.02 |   2543.41 |    7.29 |  122.39 | true |
//|        des_area |  21333.63,   21283.58,   21164.73, |   3293.60 |   3485.19 |   3323.94 |    1.15 |    6.63 | true |
//|        des_perf | 397213.78,  394451.72,  383134.06, |   2877.39 |   2877.39 |   3275.57 |   22.96 |  496.79 | true |
//|             DMA | 107967.54,  107759.85,  107080.48, |   2256.69 |   2256.69 |   2372.95 |    3.56 |   42.61 | true |
//|             DSP | 199650.59,  198797.17,  194705.56, |   7267.58 |   7589.77 |   7472.18 |    7.92 |  332.46 | true |
//|        ethernet | 228245.88,  227007.39,  221485.17, |   3322.95 |   3322.95 |   3322.95 |   56.84 | 1088.39 | true |
//|      iwls05_i2c |   5241.05,    5241.04,    5241.04, |   1441.10 |   1441.10 |   1441.10 |    0.37 |    0.81 | true |
//| iwls05_mem_ctrl |  42827.79,   42210.92,   41148.54, |   4633.83 |   4633.83 |   4587.35 |    1.50 |   19.82 | true |
//|    pci_bridge32 |  99602.79,   98266.77,   92842.50, |   3359.16 |   3359.16 |   3359.16 |    3.60 |   60.93 | true |
//|            RISC | 343479.25,  336523.94,  320604.56, |   7697.43 |   7697.43 |   8153.10 |   12.24 |  307.01 | true |
//|            sasc |   2929.34,    2914.32,    2914.32, |    973.09 |    973.09 |    973.09 |    0.31 |    0.69 | true |
//|      simple_spi |   3877.83,    3856.55,    3850.29, |   1876.28 |   1876.28 |   1876.28 |    0.34 |    1.43 | true |
//|             spi |  15882.86,   15819.05,   15737.73, |   2956.47 |   2956.47 |   2956.47 |    0.84 |    4.31 | true |
//|          ss_pcm |   2242.48,    2242.48,    2242.48, |    670.15 |    670.15 |    670.15 |    0.29 |    0.36 | true |
//|      systemcaes |  52909.36,   52041.01,   49859.82, |   3260.15 |   3245.72 |   3220.24 |    2.02 |   30.52 | true |
//|      systemcdes |  11760.04,   11699.98,   11556.04, |   3252.88 |   3464.80 |   3604.19 |    0.84 |    6.05 | true |
//|            tv80 |  34781.51,   34433.63,   34026.93, |   5541.92 |   5728.06 |   5766.63 |    2.10 |   33.62 | true |
//|       usb_funct |  69748.36,   69486.82,   69009.86, |   3073.09 |   3073.09 |   3073.09 |    2.14 |   29.87 | true |
//|         usb_phy |   2350.01,    2350.01,    2350.01, |    947.78 |    947.78 |    947.78 |    0.30 |    0.36 | true |
//|         vga_lcd | 539188.81,  537435.06,  530003.12, |   2590.79 |   2590.79 |   2590.79 |  106.91 | 3171.05 | true |
//|       wb_conmax | 158355.09,  158211.19,  158143.62, |   2630.61 |   2630.61 |   2630.61 |    5.35 |   32.69 | true |
//|             c17 |     25.03,      25.03,      25.03, |    209.99 |    209.99 |    209.99 |    0.28 |    0.34 | true |
//|            c432 |    788.36,     788.36,     788.36, |   2584.10 |   2584.10 |   2584.10 |    0.29 |    0.38 | true |
//|            c499 |   2416.83,    2403.05,    2394.28, |   1791.48 |   1791.48 |   1791.48 |    0.34 |    1.42 | true |
//|            c880 |   1762.07,    1762.07,    1762.07, |   2113.71 |   2113.71 |   2113.71 |    0.31 |    0.37 | true |
//|           c1355 |   2357.99,    2344.24,    2306.68, |   1848.83 |   1848.83 |   1933.16 |    0.34 |    2.44 | true |
//|           c1908 |   1763.43,    1707.12,    1699.61, |   1803.16 |   1803.16 |   1803.16 |    0.32 |    1.00 | true |
//|           c2670 |   3088.50,    3074.72,    3074.72, |   1759.63 |   1759.63 |   1759.63 |    0.32 |    0.70 | true |
//|           c3540 |   4522.71,    4472.63,    4413.81, |   2758.93 |   2758.93 |   2758.93 |    0.42 |    1.64 | true |
//|           c5315 |   7781.85,    7760.57,    7723.04, |   2594.00 |   2594.00 |   2594.00 |    0.42 |    1.78 | true |
//|           c6288 |  11637.67,   11616.42,   11593.89, |   7516.67 |   7523.28 |   7516.67 |    0.78 |    4.10 | true |
//|           c7552 |   7908.05,    7842.97,    7716.49, |   3886.02 |   4098.74 |   5197.30 |    0.46 |    3.28 | true |
//
////amap=np.array([   4942.62,  13358.60, 217558.84, 150191.22,  16234.48, 126844.97,  26194.90,  99702.40,  89569.38,  63432.36,   2949.64,    569.42,   2023.04,   6270.94,   1101.22, 191598.19,   2713.05,    868.49,  51364.55,  53463.37,  95644.55,  21333.63, 397213.78, 107967.54, 199650.59, 228245.88,   5241.05,  42827.79,  99602.79, 343479.25,   2929.34,   3877.83,  15882.86,   2242.48,  52909.36,  11760.04,  34781.51,  69748.36,   2350.01, 539188.81, 158355.09,     25.03,    788.36,   2416.83,   1762.07,   2357.99,   1763.43,   3088.50,   4522.71,   7781.85,  11637.67,   7908.05,])
//aop1=np.array([   4942.62,  13326.06, 217128.33, 148210.02,  16234.48, 126710.99,  26077.25,  99646.08,  89379.17,  63427.36,   2945.89,    566.92,   2023.04,   6265.94,   1101.22, 184993.12,   2701.79,    868.49,  51294.45,  53382.06,  95401.81,  21283.58, 394451.72, 107759.85, 198797.17, 227007.39,   5241.04,  42210.92,  98266.77, 336523.94,   2914.32,   3856.55,  15819.05,   2242.48,  52041.01,  11699.98,  34433.63,  69486.82,   2350.01, 537435.06, 158211.19,     25.03,    788.36,   2403.05,   1762.07,   2344.24,   1707.12,   3074.72,   4472.63,   7760.57,  11616.42,   7842.97,])
//aopN=np.array([   4942.62,   13293.52,  212440.59,  147154.97,   16234.48,  126539.47,   25841.99,   99502.21,   88700.87,   63411.11,    2945.89,     566.92,    2023.04,    6265.94,    1101.22,  151548.97,    2701.79,     868.49,   51045.39,   53216.85,   94262.04,   21164.73,  383134.06,  107080.48,  194705.56,  221485.17,    5241.04,   41148.54,   92842.50,  320604.56,    2914.32,    3850.29,   15737.73,    2242.48,   49859.82,   11556.04,   34026.93,   69009.86,    2350.01,  530003.12,  158143.62,      25.03,     788.36,    2394.28,    1762.07,    2306.68,    1699.61,    3074.72,    4413.81,    7723.04,   11593.89,    7716.49,])