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
#include <chrono>
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

  experiment<std::string, double, double, double, double, double, double, double, double, bool> exp( "SCOPT", "benchmark", "a(map)", "a(opt1)", "a(optN)", "d(map)", "d(opt1)", "d(optN)", "t(opt1)", "t(optN)", "cec");

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

  double dArea1;
  double dAreaN;
  std::vector<double> dAreas1;
  std::vector<double> dAreasN;

  for ( auto const& benchmark : all_benchmarks( iscas ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    dAreas1.clear();
    dAreasN.clear();

    if( aig.num_gates() > 300000 ) continue;

    std::vector<uint32_t> aig_size;
    std::vector<uint32_t> aig_depth;

    std::vector<double> map_size;
    std::vector<double> map_delay;

    std::vector<double> opt1_size;
    std::vector<double> opt1_delay;

    std::vector<double> optN_size;
    std::vector<double> optN_delay;

    std::vector<uint32_t> vheuristic;
    uint32_t heuristic=0;

    double amap;
    double dmap;
    double aopt1;
    double dopt1;
    double aoptN;
    double doptN;
    double time1;
    double timeN;

    bool cecMP;

    scopt::scg_network scg;

    uint32_t num_old = aig.num_gates()+1;
    uint32_t iter=0;
    while( num_old > aig.num_gates() )
    {
      if( iter == 0 )
      {
        num_old = aig.num_gates()+1;
        iter++;
        vheuristic.push_back(5);
      }
      else
      {
        num_old = aig.num_gates();

        aig = abc_opto( aig, benchmark, "resyn2rs" );
        aig = cleanup_dangling( aig );
      //  aig = abc_opto( aig, benchmark, "compress2rs" );
      //  aig = cleanup_ligs( aig );
        if( aig.num_gates() != num_old ){ vheuristic.push_back(1); }

        if( aig.num_gates() == num_old )
        {
          vheuristic.push_back(2);
        }
      }
      
        printf( "aig>>>%d\n", aig.num_gates() );

        depth_view<aig_network> daig{aig};
        aig_size.push_back( aig.num_gates() );
        aig_depth.push_back( daig.depth() );


        /* map the result */
        scopt::emap2_params ps;
        ps.cut_enumeration_ps.minimize_truth_table = true;
        ps.cut_enumeration_ps.cut_limit = 24;
        ps.area_flow_rounds=2;
        ps.area_oriented_mapping = true;
        scopt::emap2_stats st;

        scopt::scg_network scg = emap2_klut( aig, tech_lib, ps, &st );

        printf( "map>>>%f\n", scg.compute_area() );
        scg=cleanup_scg( scg );
        printf( "map*>>%f\n", scg.compute_area() );


        printf("%d %d\n", scg.num_pis(), scg.num_pos());

        cecMP = benchmark == "hyp" ? true : abc_cec( scg, benchmark );
        if(!cecMP)
          printf("ERROR\n");

        amap = scg.compute_area();
        dmap = scg.compute_worst_delay();
        map_size.push_back( amap );
        map_delay.push_back( dmap );

        /* optimize the design */
        boptimizer_params rps;
        rps.progress =false;
        rps.max_inserts = 300;
        rps.max_trials = 1;
        rps.max_pis = 16;
        rps.verbose = false;
        rps.max_divisors = 128;
        boptimizer_stats rst_p1;

        uint32_t step{1};
        auto beg1 = std::chrono::high_resolution_clock::now();
        auto begN = std::chrono::high_resolution_clock::now();

        boptimize_sc<scopt::support_selection_t::GREEDY, 4u, 4u>( scg, rps, &rst_p1 );
        printf("opt%2d>: %6f \n", step, scg.compute_area() );
        scg=cleanup_scg( scg );
        printf("opt1*>: %6f \n", step++, scg.compute_area() );
        printf("%d %d\n", scg.num_pis(), scg.num_pos());
        std::cout << std::endl;

        auto end1 = std::chrono::high_resolution_clock::now();
        auto duration1 = std::chrono::duration<double>(end1 - beg1).count();

        aopt1 = scg.compute_area();
        dopt1 = scg.compute_worst_delay();
        opt1_size.push_back( aopt1 );
        opt1_delay.push_back( dopt1 );

        double aOld = scg.compute_area()+1;
        for( int iA{0}; iA<4; ++iA )
        {
          aOld = scg.compute_area();

          boptimize_sc<scopt::support_selection_t::GREEDY, 4u, 4u>( scg, rps, &rst_p1 );

          printf("opt%2d>: %6f\n", step, scg.compute_area() );
          scg=cleanup_scg( scg );
          printf("opt1%2d>: %6f \n", step++, scg.compute_area() );

          std::cout << std::endl;


        }
        auto endN = std::chrono::high_resolution_clock::now();
        auto durationN = std::chrono::duration<double>(endN - begN).count();

        time1 = (double)duration1;
        timeN = (double)durationN;

        aoptN = scg.compute_area();
        doptN = scg.compute_worst_delay();
        optN_size.push_back( aoptN );
        optN_delay.push_back( doptN );

        cecMP = benchmark == "hyp" ? true : abc_cec( scg, benchmark );
        if(!cecMP)
          printf("ERROR\n");
        std::cout << std::endl;

    }


    printf("aaig=np.array([");
    for( int i{0}; i< aig_size.size(); ++i )
    {
        if( i == aig_size.size()-1 )
            printf("%d])\n", aig_size[i]);
        else
            printf("%d, ", aig_size[i]);
    }

    printf("amap=np.array([");
    for( int i{0}; i< map_size.size(); ++i )
    {
        if( i == map_size.size()-1 )
            printf("%f])\n", map_size[i]);
        else
            printf("%f, ", map_size[i]);
    }


    printf("aopt1=np.array([");
    for( int i{0}; i< opt1_size.size(); ++i )
    {
        if( i == opt1_size.size()-1 )
            printf("%f])\n", opt1_size[i]);
        else
            printf("%f, ", opt1_size[i]);
    }

    printf("aoptN=np.array([");
    for( int i{0}; i< optN_size.size(); ++i )
    {
        if( i == optN_size.size()-1 )
            printf("%f])\n", optN_size[i]);
        else
            printf("%f, ", optN_size[i]);
    }

    printf("color=np.array([");
    for( int i{0}; i< vheuristic.size(); ++i )
    {
        if( i == vheuristic.size()-1 )
            printf("%d])\n", vheuristic[i]);
        else
            printf("%d, ", vheuristic[i]);
    }

    printf("d(aig)=[");
    for( int i{0}; i< aig_depth.size(); ++i )
    {
        if( i == aig_depth.size()-1 )
            printf("%d])\n", aig_depth[i]);
        else
            printf("%d, ", aig_depth[i]);
    }

    printf("d(map)=[");
    for( int i{0}; i< map_delay.size(); ++i )
    {
        if( i == map_delay.size()-1 )
            printf("%f])\n", map_delay[i]);
        else
            printf("%f, ", map_delay[i]);
    }

    printf("d(opt1)=[");
    for( int i{0}; i< opt1_delay.size(); ++i )
    {
        if( i == opt1_delay.size()-1 )
            printf("%f])\n", opt1_delay[i]);
        else
            printf("%f, ", opt1_delay[i]);
    }

    printf("d(optN)=[");
    for( int i{0}; i< optN_delay.size(); ++i )
    {
        if( i == optN_delay.size()-1 )
            printf("%f])\n", optN_delay[i]);
        else
            printf("%f, ", optN_delay[i]);
    }

    dArea1 = ( 100*opt1_size[opt1_size.size()-1]-opt1_size[0] )/opt1_size[0];
    printf("  1)%.2f\n", dArea1 );

    dAreaN = ( 100*optN_size[optN_size.size()-1]-optN_size[0] )/optN_size[0];
    printf("INF)%.2f\n", dAreaN );

    dAreas1.push_back( dArea1 );
    dAreasN.push_back( dAreaN );

    const auto cec = benchmark == "hyp" ? true : abc_cec( aig, benchmark );
    assert( cec && "[e] not equivalent" );

    std::cout << std::endl;

    if(!cecMP)
      printf("ERROR\n");
    std::cout << std::endl;

    write_aiger( aig, benchmark+"_optmap.aig" );
    exp( benchmark, map_size.back(), opt1_size.back(), optN_size.back(), map_delay.back(), opt1_delay.back(), optN_delay.back(), time1, timeN, cecMP );

  }

  exp.save();
  exp.table();

  double avg1{0};
  for( int i{0}; i< dAreas1.size(); ++i )
  {
    avg1+=dAreas1[i]/((double)dAreas1.size());
  }

  double avgN{0};
  for( int i{0}; i< dAreasN.size(); ++i )
  {
    avgN+=dAreasN[i]/((double)dAreasN.size());
  }

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

//           c17   & $     25.03 $ & $     25.03 $ & $     25.03 $ & $    209.99 $ & $     209.99 $ & $     209.99 $ & $   2.71 $ & $  13.62 $\\
//          c432   & $    838.41 $ & $    820.89 $ & $    820.89 $ & $   2524.15 $ & $    2524.15 $ & $    2524.15 $ & $   2.73 $ & $  13.69 $\\
//          c499   & $   2443.07 $ & $   2438.07 $ & $   2401.79 $ & $   1841.29 $ & $    1841.29 $ & $    1841.29 $ & $   2.74 $ & $  13.80 $\\
//          c880   & $   1762.07 $ & $   1762.07 $ & $   1762.07 $ & $   2113.71 $ & $    2113.71 $ & $    2113.71 $ & $   2.77 $ & $  13.74 $\\
//         c1355   & $   2275.43 $ & $   2275.43 $ & $   2275.43 $ & $   1745.31 $ & $    1745.31 $ & $    1745.31 $ & $   2.76 $ & $  13.85 $\\
//         c1908   & $   1712.14 $ & $   1689.62 $ & $   1687.12 $ & $   1826.03 $ & $    1826.03 $ & $    1826.03 $ & $   2.75 $ & $  13.67 $\\
//         c2670   & $   3077.19 $ & $   3032.15 $ & $   3032.15 $ & $   1746.49 $ & $    1746.49 $ & $    1746.49 $ & $   2.72 $ & $  13.71 $\\
//         c3540   & $   4451.42 $ & $   4426.40 $ & $   4387.60 $ & $   2887.67 $ & $    2930.85 $ & $    2944.98 $ & $   2.79 $ & $  13.87 $\\
//         c5315   & $   7233.58 $ & $   7179.75 $ & $   7150.97 $ & $   2724.68 $ & $    3215.42 $ & $    3215.42 $ & $   2.79 $ & $  13.94 $\\
//         c6288   & $  11637.67 $ & $  11617.67 $ & $  11582.65 $ & $   7516.67 $ & $    7516.67 $ & $    7585.48 $ & $   2.92 $ & $  14.62 $\\
//         c7552   & $   7951.84 $ & $   7880.51 $ & $   7784.13 $ & $   4164.32 $ & $    4841.63 $ & $    4812.35 $ & $   2.81 $ & $  14.10 $\\
//          adder  & $   4942.62 $ & $   4942.62 $ & $   4942.62 $ & $  17141.27 $ & $   17141.27 $ & $   17141.27 $ & $   2.78 $ & $  14.00 $\\
//            bar  & $  13358.60 $ & $  13294.76 $ & $  13289.75 $ & $   1563.40 $ & $    1563.40 $ & $    1563.40 $ & $   2.89 $ & $  14.38 $\\
//            div  & $ 112126.52 $ & $ 111801.96 $ & $ 111769.43 $ & $ 307669.62 $ & $  308302.66 $ & $  309277.50 $ & $   5.58 $ & $  27.69 $\\
//            hyp  & $1108772.38 $ & $1104516.38 $ & $1101271.12 $ & $1545227.12 $ & $ 1543656.50 $ & $ 1540590.00 $ & $  20.52 $ & $  93.36 $\\
//           log2  & $ 147611.86 $ & $ 147108.53 $ & $ 146846.98 $ & $  28357.35 $ & $   28417.24 $ & $   28583.88 $ & $  11.05 $ & $  53.48 $\\
//            max  & $  16234.48 $ & $  16178.21 $ & $  16130.71 $ & $  19885.02 $ & $   19863.55 $ & $   19863.55 $ & $   2.86 $ & $  14.27 $\\
//     multiplier  & $ 127161.61 $ & $ 126727.03 $ & $ 126589.35 $ & $  18574.54 $ & $   18574.54 $ & $   18544.39 $ & $   5.64 $ & $  27.76 $\\
//            sin  & $  26179.92 $ & $  26035.99 $ & $  25974.66 $ & $  15367.77 $ & $   15348.76 $ & $   15355.94 $ & $   4.34 $ & $  21.65 $\\
//           sqrt  & $  99993.92 $ & $  99846.38 $ & $  99746.34 $ & $ 381985.88 $ & $  381836.34 $ & $  381692.25 $ & $   4.63 $ & $  21.96 $\\
//         square  & $  89616.99 $ & $  88929.99 $ & $  88436.99 $ & $  18173.36 $ & $   18353.45 $ & $   18276.24 $ & $   4.01 $ & $  19.61 $\\
//        arbiter  & $  63432.36 $ & $  63393.61 $ & $  63356.11 $ & $   6484.11 $ & $    6511.00 $ & $    6511.00 $ & $   3.26 $ & $  16.23 $\\
//          cavlc  & $   3018.45 $ & $   2999.68 $ & $   2990.93 $ & $   1461.20 $ & $    1552.85 $ & $    1575.84 $ & $   2.80 $ & $  14.02 $\\
//           ctrl  & $    576.95 $ & $    574.45 $ & $    574.45 $ & $    670.87 $ & $     814.56 $ & $     814.56 $ & $   2.73 $ & $  13.72 $\\
//            dec  & $   2023.04 $ & $   2023.04 $ & $   2023.04 $ & $    484.51 $ & $     484.51 $ & $     484.51 $ & $   2.78 $ & $  13.88 $\\
//            i2c  & $   6028.19 $ & $   6018.18 $ & $   5989.40 $ & $   1756.36 $ & $    1756.36 $ & $    1756.36 $ & $   2.77 $ & $  13.93 $\\
//      int2float  & $   1076.22 $ & $   1074.97 $ & $   1074.97 $ & $   1103.79 $ & $    1103.79 $ & $    1103.79 $ & $   2.76 $ & $  13.75 $\\
//       priority  & $   2701.78 $ & $   2678.00 $ & $   2652.97 $ & $   5253.41 $ & $    5253.41 $ & $    5506.55 $ & $   2.73 $ & $  13.74 $\\
//         router  & $    893.55 $ & $    893.55 $ & $    893.55 $ & $   1863.52 $ & $    1863.52 $ & $    1863.52 $ & $   2.74 $ & $  13.76 $\\
//          voter  & $  50786.24 $ & $  50508.43 $ & $  49972.60 $ & $   4951.83 $ & $    5015.06 $ & $    5085.40 $ & $   3.84 $ & $  18.69 $\\
//       ac97_ctrl & $  52416.48 $ & $  49904.27 $ & $  49765.43 $ & $   1043.31 $ & $    1043.31 $ & $    1043.31 $ & $   3.10 $ & $  15.28 $\\
//        aes_core & $  95731.90 $ & $  94970.10 $ & $  94438.53 $ & $   2427.92 $ & $    2427.92 $ & $    2536.11 $ & $   5.70 $ & $  27.98 $\\
//        des_area & $  21899.37 $ & $  21823.04 $ & $  21731.69 $ & $   3115.60 $ & $    3261.92 $ & $    3261.92 $ & $   3.18 $ & $  15.89 $\\
//        des_perf & $ 399137.53 $ & $ 393128.50 $ & $ 388168.06 $ & $   2988.11 $ & $    2980.23 $ & $    3176.13 $ & $  13.38 $ & $  65.36 $\\
//             DMA & $ 108449.49 $ & $ 106867.04 $ & $ 106208.92 $ & $   2500.26 $ & $    2402.90 $ & $    2651.39 $ & $   4.05 $ & $  20.22 $\\
//             DSP & $ 199671.31 $ & $ 196423.53 $ & $ 194724.53 $ & $   6208.11 $ & $    6109.41 $ & $    6214.93 $ & $   5.90 $ & $  28.86 $\\
//        ethernet & $ 223277.80 $ & $ 220501.53 $ & $ 219803.27 $ & $   3384.71 $ & $    3384.71 $ & $    3451.63 $ & $  11.70 $ & $  60.46 $\\
//      iwls05_i2c & $   5239.77 $ & $   5213.50 $ & $   5182.22 $ & $   1617.45 $ & $    1583.87 $ & $    1583.87 $ & $   2.78 $ & $  13.93 $\\
// iwls05_mem_ctrl & $  38946.99 $ & $  38357.81 $ & $  38171.43 $ & $   5116.92 $ & $    5116.92 $ & $    5116.92 $ & $   3.07 $ & $  15.43 $\\
//    pci_bridge32 & $ 101899.92 $ & $  94492.83 $ & $  92512.58 $ & $   3375.91 $ & $    3371.61 $ & $    3371.61 $ & $   4.11 $ & $  20.28 $\\
//            sasc & $   2909.32 $ & $   2876.78 $ & $   2863.02 $ & $    836.50 $ & $     836.50 $ & $     836.50 $ & $   2.78 $ & $  13.78 $\\
//      simple_spi & $   3934.15 $ & $   3926.65 $ & $   3920.39 $ & $   1491.10 $ & $    1590.00 $ & $    1590.00 $ & $   2.78 $ & $  13.88 $\\
//             spi & $  16162.87 $ & $  16035.26 $ & $  15921.44 $ & $   2832.56 $ & $    2876.25 $ & $    2876.25 $ & $   2.98 $ & $  15.07 $\\
//          ss_pcm & $   2242.48 $ & $   2152.40 $ & $   2152.40 $ & $    670.15 $ & $     670.15 $ & $     670.15 $ & $   2.76 $ & $  13.77 $\\
//      systemcaes & $  52571.62 $ & $  50299.06 $ & $  49162.77 $ & $   3417.80 $ & $    3641.01 $ & $    3493.54 $ & $   3.50 $ & $  17.35 $\\
//      systemcdes & $  14030.64 $ & $  13805.45 $ & $  13705.38 $ & $   2928.76 $ & $    3264.41 $ & $    3264.41 $ & $   2.97 $ & $  14.85 $\\
//            tv80 & $  33896.60 $ & $  33746.45 $ & $  33577.53 $ & $   5774.83 $ & $    5774.83 $ & $    5763.81 $ & $   3.42 $ & $  16.87 $\\
//       usb_funct & $  69071.21 $ & $  68896.04 $ & $  68749.62 $ & $   3649.50 $ & $    3342.98 $ & $    3342.98 $ & $   3.47 $ & $  17.26 $\\
//         usb_phy & $   2332.48 $ & $   2316.20 $ & $   2306.19 $ & $    875.32 $ & $     873.69 $ & $     873.69 $ & $   2.77 $ & $  13.86 $\\
//         vga_lcd & $ 546428.50 $ & $ 524047.78 $ & $ 508026.56 $ & $   2714.57 $ & $    2714.57 $ & $    2725.27 $ & $  40.38 $ & $ 193.98 $\\
//       wb_conmax & $ 164650.84 $ & $ 164237.98 $ & $ 163716.23 $ & $   2399.34 $ & $    2404.18 $ & $    2404.18 $ & $   5.29 $ & $  26.63 $\\
//
//amap=np.array([25.03,838.41,2443.07,1762.07,2275.43,1712.14,3077.19,4451.42,7233.58,11637.67,7951.84,4942.62,13358.60,112126.52,1108772.38,147611.86,16234.48,127161.61,26179.92,99993.92,89616.99,63432.36,3018.45,576.95,2023.04,6028.19,1076.22,2701.78,893.55,50786.24,52416.48,95731.90,21899.37,399137.53,108449.49,199671.31,223277.80,5239.77,38946.99,101899.92,2909.32,3934.15,16162.87,2242.48,52571.62,14030.64,33896.60,69071.21,2332.48,546428.50,164650.84])
//aopt1=np.array([     25.03,    820.89,   2438.07,   1762.07,   2275.43,   1689.62,   3032.15,   4426.40,   7179.75,  11617.67,   7880.51,   4942.62,  13294.76, 111801.96,1104516.38, 147108.53,  16178.21, 126727.03,  26035.99,  99846.38,  88929.99,  63393.61,   2999.68,    574.45,   2023.04,   6018.18,   1074.97,   2678.00,    893.55,  50508.43,  49904.27,  94970.10,  21823.04, 393128.50, 106867.04, 196423.53, 220501.53,   5213.50,  38357.81,  94492.83,   2876.78,   3926.65,  16035.26,   2152.40,  50299.06,  13805.45,  33746.45,  68896.04,   2316.20, 524047.78, 164237.98])
//aoptN=np.array([ 25.03,    820.89,   2401.79,   1762.07,   2275.43,   1687.12,   3032.15,   4387.60,   7150.97,  11582.65,   7784.13,   4942.62,  13289.75, 111769.43,1101271.12, 146846.98,  16130.71, 126589.35,  25974.66,  99746.34,  88436.99,  63356.11,   2990.93,    574.45,   2023.04,   5989.40,   1074.97,   2652.97,    893.55,  49972.60,  49765.43,  94438.53,  21731.69, 388168.06, 106208.92, 194724.53, 219803.27,   5182.22,  38171.43,  92512.58,   2863.02,   3920.39,  15921.44,   2152.40,  49162.77,  13705.38,  33577.53,  68749.62,   2306.19, 508026.56, 163716.23])
//

//dmap=np.array([209.99,2524.15,1841.29,2113.71,1745.31,1826.03,1746.49,2887.67,2724.68,7516.67,4164.32,17141.27,1563.40,307669.62,1545227.12,28357.35,19885.02,18574.54,15367.77,381985.88,18173.36,6484.11,1461.20,670.87,484.51,1756.36,1103.79,5253.41,1863.52,4951.83,1043.31,2427.92,3115.60,2988.11,2500.26,6208.11,3384.71,1617.45,5116.92,3375.91,836.50,1491.10,2832.56,670.15,3417.80,2928.76,5774.83,3649.50,875.32,2714.57,2399.34])
//dopt1=np.array([209.99 ,2524.15 ,1841.29 ,2113.71 ,1745.31 ,1826.03 ,1746.49 ,2930.85 ,3215.42 ,7516.67 ,4841.63 ,17141.27 ,1563.40 ,308302.66 ,1543656.50 ,28417.24 ,19863.55 ,18574.54 ,15348.76 ,381836.34 ,18353.45 ,6511.00 ,1552.85 ,814.56 ,484.51 ,1756.36 ,1103.79 ,5253.41 ,1863.52 ,5015.06 ,1043.31 ,2427.92 ,3261.92 ,2980.23 ,2402.90 ,6109.41 ,3384.71 ,1583.87 ,5116.92 ,3371.61 ,836.50 ,1590.00 ,2876.25 ,670.15 ,3641.01 ,3264.41 ,5774.83 ,3342.98 ,873.69 ,2714.57 ,2404.18 ])
//doptN=np.array([209.99,2524.15,1841.29,2113.71,1745.31,1826.03,1746.49,2944.98,3215.42,7585.48,4812.35,17141.27,1563.40,309277.50,1540590.00,28583.88,19863.55,18544.39,15355.94,381692.25,18276.24,6511.00,1575.84,814.56,484.51,1756.36,1103.79,5506.55,1863.52,5085.40,1043.31,2536.11,3261.92,3176.13,2651.39,6214.93,3451.63,1583.87,5116.92,3371.61,836.50,1590.00,2876.25,670.15,3493.54,3264.41,5763.81,3342.98,873.69,2725.27,2404.18])