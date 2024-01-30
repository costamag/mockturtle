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
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/rig.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/bench_reader.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/io/write_bench.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/mapping_view.hpp>
#include <mockturtle/algorithms/sim_resub.hpp>


#include <experiments.hpp>

using namespace mockturtle;
using namespace rils;

template<class Ntk>
std::tuple<uint32_t,uint32_t,float> abc_mfs( Ntk const& ntk, std::string benchmark )
{
  mockturtle::write_bench( ntk, "/tmp/mfsin_"+benchmark+".bench" );
  //std::string command = fmt::format( "abc -q \"read_bench /tmp/mfsin_{}.bench; mfs -D 0 -dea; write_bench /tmp/mfsout_{}.bench;\"", benchmark, benchmark );
  std::string command = fmt::format( "abc -q \"read_bench /tmp/mfsin_{}.bench; mfs2 -L 5 -ea; time; &get -mn; &ps;\"", benchmark );

  std::array<char, 128> buffer;
  std::string result;
#if WIN32
  std::unique_ptr<FILE, decltype( &_pclose )> pipe( _popen( command.c_str(), "r" ), _pclose );
#else
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
#endif
  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }

  /* search for one line which says "Networks are equivalent" and ignore all other debug output from ABC */
  std::stringstream ss( result );
  std::string line;
  std::tuple<uint32_t,uint32_t,float> res;
  while ( std::getline( ss, line, '\n' ) )
  {

    std::vector<std::string> words;

    for (int i = 0; i < line.size(); i++)
    {
        std::string tmp = "";
        while (line[i] != ' ' && i < line.size())
        {
            if ((isalpha(line.at(i)) || isdigit(line.at(i)) || line.at(i)=='.' )){
                tmp += line[i];
            }
            else if (line.at(i) == ' ' )
            {
                i++;
                break;
            }
            i++;
        }
        words.push_back(tmp);
    }
    if( words[0]=="elapse" )
    {
      //std::cout << words[1] << std::endl;
      std::get<2>(res) = std::stof(words[1]);
    }
    if ( line.size() >= 23u && line.substr( 25u, 3u ) == "lut" )
    {
      std::get<0>(res) = std::stoi(line.substr( 30u, 9u ));
      std::get<1>(res) = std::stoi(line.substr( 82u, 15u ));
      return res;
    }
  }

  return res;

}

int main()
{
  using namespace experiments;
  using namespace mockturtle;
  using namespace rils;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, double, uint32_t, uint32_t, double, bool, bool, bool> exp( "rig_exp_3", "benchmark", "luts", "lut_depth", "rigs", "rigs_depth", "rs rigs", "rs rigs_depth", "t(spf)", "rs-mfs rigs", "rs-mfs rigs_depth", "t(mfs)", "eq(LUT)", "eq(RIG)", "eq(RS)" );

  static constexpr uint32_t K = 4;

  for ( auto const& benchmark : all_benchmarks( iscas ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }
    if( aig.size() > 300000 ) continue;


    lut_map_params lps;
    lps.cut_enumeration_ps.cut_size = K;
    lps.cut_enumeration_ps.cut_limit = 8u;
    lps.recompute_cuts = true;
    lps.area_oriented_mapping = true;
    lps.cut_expansion = true;
    lut_map_stats st;
    auto klut = lut_map( aig, lps, &st );

    klut = cleanup_luts(klut);
    depth_view<klut_network> klut_d{ klut };
    printf("#LUTS[0]=%d\n", klut.num_gates());

    auto const cec = benchmark == "hyp" ? true : abc_cec( klut, benchmark );

    std::string tmp = benchmark + "tmp.blif";
    write_blif( klut, tmp );
    
    rig_network rig;
    if ( lorina::read_blif( tmp, blif_reader( rig ) ) != lorina::return_code::success )
    {
      printf("rig unsuccessful\n");
      continue;
    }
    depth_view<rig_network> rig_d{ rig };
    auto const rig_cec = benchmark == "hyp" ? true : abc_cec( rig, benchmark );

    uint32_t rig_num_gates = rig.num_gates();
    uint32_t rig_depth = rig_d.depth();

    resubstitution_params rps;
    resubstitution_stats rst;
    // ps.pattern_filename = "1024sa1/" + benchmark + ".pat";
    rps.progress =true;
    rps.max_inserts = 20;
    rps.max_trials = 100;
    rps.max_pis = 10;
    rps.max_divisors = std::numeric_limits<uint32_t>::max();

    rig_resubstitution<rils::network_t::kLUT, rils::support_selection_t::PIVOT, K>( rig, rps, &rst );
    rig = cleanup_dangling( rig );

    printf("spf %d\n", rig.num_gates() );

    depth_view<rig_network> rs_rig_d{ rig };

    uint32_t rs_rig_num_gates = rig.num_gates();
    uint32_t rs_rig_depth = rs_rig_d.depth();

    const auto cec_rs = klut.num_gates() > 25000 ? true : abc_cec( rig, benchmark );

    auto mfs_res = abc_mfs( klut, benchmark );

    printf("mfs %d\n", std::get<0>(mfs_res));

    exp( benchmark, klut.num_gates(), klut_d.depth(), rig_num_gates, rig_depth, rs_rig_num_gates, rs_rig_depth, to_seconds(rst.time_total), std::get<0>(mfs_res), std::get<1>(mfs_res), std::get<2>(mfs_res), cec, rig_cec, cec_rs );
  }

  exp.save();
  exp.table();

  return 0;
}


//
//|  benchmark & $  rigs | rigs_depth | rs rigs | rs rigs_depth | t(spf) | rs-mfs rigs | rs-mfs rigs_depth | t(mfs) |
//|       c432 & $    66 |         14 |      66 |            14 |   0.03 |          64 |                14 |   0.09 |
//|       c880 & $   112 |         15 |     107 |            14 |   0.08 |         112 |                15 |   0.11 |
//|      c1908 & $    90 |         12 |      79 |            10 |   0.03 |          87 |                12 |   0.10 |
//|      c2670 & $   174 |          9 |     168 |            10 |   0.08 |         174 |                 9 |   0.17 |
//|      c3540 & $   343 |         18 |     332 |            18 |   0.24 |         322 |                18 |   0.75 |
//|      c5315 & $   471 |         14 |     428 |            14 |   0.29 |         455 |                14 |   0.47 |
//|      c6288 & $   495 |         30 |     494 |            30 |   0.01 |         495 |                30 |   0.21 |
//|      c7552 & $   444 |         15 |     423 |            19 |   0.28 |         439 |                15 |   0.87 |
//|        div & $ 24482 |       2185 |   24242 |          2173 |   2.15 |       24411 |              2185 |   6.00 |
//|        hyp & $ 63017 |       8529 |   62837 |          8521 |   6.81 |       62936 |              8529 |  64.74 |
//|       log2 & $  9698 |        166 |    9589 |           154 |  15.16 |        9693 |               166 |   7.04 |
//| multiplier & $  7216 |        129 |    7215 |           129 |  19.26 |        7216 |               129 |  11.33 |
//|        sin & $  1798 |         94 |    1791 |            94 |   2.05 |        1796 |                94 |   2.91 |
//|       sqrt & $  6626 |       2121 |    6527 |          2121 |   0.54 |        6598 |              2121 |   1.57 |
//|     square & $  6159 |        125 |    6006 |           125 |   1.44 |        6014 |               125 |   1.64 |
//|      cavlc & $   286 |          9 |     269 |            10 |   7.30 |         286 |                 9 |   0.20 |
//|       ctrl & $    53 |          3 |      41 |             8 |   0.04 |          51 |                 3 |   0.03 |
//|        i2c & $   527 |          9 |     491 |             9 |   0.35 |         506 |                 8 |   0.23 |
//|  int2float & $    95 |          8 |      88 |             8 |   0.27 |          89 |                 8 |   0.06 |
//|   mem_ctrl & $ 17667 |         70 |   16051 |            69 |  17.90 |       16245 |                70 |  47.12 |
//|   priority & $   327 |        122 |     294 |           122 |   0.03 |         297 |               122 |   0.26 |
//|     router & $    88 |         32 |      85 |            32 |   0.05 |          78 |                30 |   0.06 |
//
//[66,107,79,168,332,428,494,423,24242,62837,9589,7215,1791,6527,6006,269,41,491,88,16051,294,85]
//[64,112,87,174,322,455,495,439,24411,62936,9693,7216,1796,6598,6014,286,51,506,89,16245,297,78]
//2.13%

//| benchmark | luts | lut_depth | rigs | rigs_depth | rs rigs | rs rigs_depth | t(spf) | rs-mfs rigs | rs-mfs rigs_depth | t(mfs) | eq(LUT) | eq(RIG) | eq(RS) |
//|       c17 |    2 |         1 |    2 |          1 |       2 |             1 |   0.00 |           2 |                 1 |   0.03 |    true |    true |   true |
//|      c432 |   66 |        14 |   66 |         14 |      66 |            14 |   0.03 |          64 |                14 |   0.03 |    true |    true |   true |
//|      c499 |   78 |         7 |   78 |          7 |      78 |             7 |   0.02 |          78 |                 7 |   0.03 |    true |    true |   true |
//|      c880 |  112 |        15 |  112 |         15 |     107 |            14 |   0.08 |         112 |                15 |   0.03 |    true |    true |   true |
//|     c1355 |   78 |         7 |   78 |          7 |      78 |             7 |   0.02 |          78 |                 7 |   0.03 |    true |    true |   true |
//|     c1908 |   90 |        12 |   90 |         12 |      79 |            10 |   0.03 |          86 |                12 |   0.03 |    true |    true |   true |
//|     c2670 |  185 |         9 |  174 |          9 |     168 |            10 |   0.08 |         170 |                 9 |   0.04 |    true |    true |   true |
//|     c3540 |  344 |        18 |  343 |         18 |     332 |            18 |   0.26 |         323 |                18 |   0.06 |    true |    true |   true |
//|     c5315 |  487 |        14 |  471 |         14 |     428 |            14 |   0.32 |         455 |                14 |   0.06 |    true |    true |   true |
//|     c6288 |  495 |        30 |  495 |         30 |     494 |            30 |   0.04 |         495 |                30 |   0.05 |    true |    true |   true |
//|     c7552 |  445 |        15 |  444 |         15 |     423 |            19 |   0.31 |         432 |                17 |   0.06 |    true |    true |   true |