/* mockturtle: C++ logic network library
 * Copylight (C) 2018-2019  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the lights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copylight notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYligHT
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
#include <mockturtle/networks/lig.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/bench_reader.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/write_bench.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/mapping_view.hpp>
#include <mockturtle/algorithms/boptimizer.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/algorithms/rewrite.hpp>

#include <experiments.hpp>

using namespace mockturtle;
using namespace rils;

template<class Ntk>
std::tuple<uint32_t,uint32_t,float> abc_eval( Ntk const& ntk, std::string benchmark )
{
  mockturtle::write_bench( ntk, "/tmp/mfsin2_"+benchmark+".bench" );
  //std::string command = fmt::format( "abc -q \"read_bench /tmp/mfsin_{}.bench; mfs -D 0 -dea; write_bench /tmp/mfsout_{}.bench;\"", benchmark, benchmark );
  std::string command = fmt::format( "abc -q \"read_bench /tmp/mfsin2_{}.bench; &get -mn; &ps;\"", benchmark );

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

  static constexpr uint32_t K = 4;

  for ( auto const& benchmark : all_benchmarks(  ) )
  {
    if( benchmark == "hyp" ) continue;
    //fmt::print( "[i] processing {}\n", benchmark );

    std::string tmp = benchmark + "_exp1.blif";

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    if( aig.num_gates() > 100000 ) continue;

    /* set up the parameters */
    boptimizer_stats rst;
    boptimizer_params rps;
    rps.progress = false;
    rps.max_inserts = 20;
    rps.max_trials = 100;
    rps.max_pis = 16;
    rps.verbose = false;
    rps.max_divisors = 64;

    /* CASE 0: parameters */
    //klut_network klut=abc_if( aig, benchmark, 4u );

    scopt::lut_map2_params ps;
    ps.cut_enumeration_ps.cut_size = K;
    ps.cut_enumeration_ps.cut_limit = 8u;
    ps.recompute_cuts = true;
    ps.area_oriented_mapping = true;
    ps.cut_expansion = true;
    scopt::lut_map2_stats st;

    auto lig = scopt::lut_map2( aig, ps, &st );
    if( lig.num_gates() > 2000 ) continue;

    depth_view<lig_network> lig_d{ lig };
    auto const lig_cec = benchmark == "hyp" ? true : abc_cec( lig, benchmark );

    uint32_t map_num_gates = lig.num_gates();
    printf("MP : %6d\n", map_num_gates );
    uint32_t map_depth = lig_d.depth();

    boptimize_klut<ENU, K, K>( lig, rps, &rst );
    lig = cleanup_luts( lig );
    printf("ENU[4,4]: %6d [%6d]\n", lig.num_gates(), lig.max_num_fanins);

    depth_view<lig_network> lig1_d{ lig };

    auto piv_res = abc_eval( lig, benchmark );
    uint32_t p1_num_gates = std::get<0>(piv_res);
    uint32_t p1_depth = std::get<1>(piv_res);
    double p1_time = to_seconds( rst.time_total );

    bool const lig_p1_cec = benchmark == "hyp" ? true : abc_cec( lig, benchmark );
    if( !lig_p1_cec ) 
    {
      printf("NEQ\n");
      continue;
    }
  }

  return 0;
}

//| benchmark | a(map) | a(mfs) | a(mfs2) | a(g,1) | a(p,1) | d(map) | d(mfs) | d(mfs2) | d(g,1) | d(p,1) | t(mfs) | t(mfs2) | t(g,1) | t(p,1) |
//|       c17 |      2 |      2 |       2 |      2 |      2 |      1 |      1 |       1 |      1 |      1 |   0.00 |    0.00 |   0.00 |   2.66 |
//|      c432 |     66 |     64 |      64 |     66 |     64 |     14 |     14 |      14 |     14 |     15 |   0.00 |    0.00 |   0.00 |   2.65 |
//|      c499 |     78 |     78 |      78 |     78 |     78 |      7 |      7 |       7 |      7 |      8 |   0.00 |    0.00 |   0.00 |   2.74 |
//|      c880 |    112 |    112 |     112 |    111 |    107 |     15 |     15 |      15 |     15 |     16 |   0.00 |    0.00 |   0.00 |   2.73 |
//|     c1355 |     78 |     78 |      78 |     78 |     78 |      7 |      7 |       7 |      7 |      8 |   0.00 |    0.00 |   0.00 |   3.33 |
//|     c1908 |     90 |     84 |      84 |     90 |     85 |     12 |     10 |      10 |     12 |     10 |   0.00 |    0.00 |   0.00 |   2.96 |
//|     c2670 |    174 |    173 |     172 |    155 |    170 |      9 |      9 |       9 |     16 |     10 |   0.00 |    0.00 |   0.00 |   2.70 |
//|     c3540 |    343 |    342 |     322 |    329 |    340 |     18 |     18 |      18 |     18 |     19 |   0.00 |    0.00 |   0.00 |   2.79 |
//|     c5315 |    471 |    455 |     447 |    452 |    433 |     14 |     14 |      14 |     13 |     15 |   0.00 |    0.00 |   0.00 |   2.87 |
//|     c6288 |    495 |    495 |     495 |    495 |    495 |     30 |     30 |      30 |     30 |     57 |   0.00 |    0.00 |   0.00 |   3.16 |
//|     c7552 |    444 |    431 |     431 |    410 |    420 |     15 |     17 |      17 |     19 |     34 |   0.00 |    0.00 |   0.00 |   3.00 |