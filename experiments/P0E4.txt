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
#include <mockturtle/networks/lig.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/bench_reader.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/io/write_bench.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/mapping_view.hpp>
#include <mockturtle/algorithms/boptimizer.hpp>
#include <chrono>

#include <experiments.hpp>

using namespace std::chrono;

using namespace mockturtle;

std::tuple<uint32_t,uint32_t,float> abc_mfs( lig_network& ntk, std::string benchmark, uint32_t M = 5000 )
{//mfs2 -L 20 -ea
  mockturtle::write_bench( ntk, "/tmp/mfsin_"+benchmark+".bench" );
  //std::string command = fmt::format( "abc -q \"read_bench /tmp/mfsin_{}.bench; mfs -D 0 -dea; write_bench /tmp/mfsout_{}.bench;\"", benchmark, benchmark );
  std::string command = fmt::format( "abc -q \"read_bench /tmp/mfsin_{}.bench; mfs -W 4 -M " + std::to_string(M) + " -L 200; time; write_blif /tmp/mfsin_{}.blif; &get -mn; &ps;\"", benchmark, benchmark );

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

  lig_network ntk_res;

  std::string string_path = ("/tmp/mfsin_"+benchmark+".blif");
  if( lorina::read_blif( string_path, blif_reader( ntk_res ) ) != lorina::return_code::success )
    std::cerr << "read_blif failed" << std::endl;
  ntk = ntk_res;

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

std::tuple<uint32_t,uint32_t,float> abc_mfs2( lig_network& ntk, std::string benchmark, uint32_t M = 5000 )
{
  mockturtle::write_bench( ntk, "/tmp/mfsin2_"+benchmark+".bench" ); 

  //std::string command = fmt::format( "abc -q \"read_bench /tmp/mfsin_{}.bench; mfs -D 0 -dea; write_bench /tmp/mfsout_{}.bench;\"", benchmark, benchmark );
  std::string command = fmt::format( "abc -q \"read_bench /tmp/mfsin2_{}.bench; mfs2 -e -W 4 -M " + std::to_string(M) + " -L 200; time; write_blif /tmp/mfsin2_{}.blif; &get -mn; &ps;\"", benchmark, benchmark );

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

  lig_network ntk_res;

  std::string string_path = ("/tmp/mfsin2_"+benchmark+".blif");
  if( lorina::read_blif( string_path, blif_reader( ntk_res ) ) != lorina::return_code::success )
    std::cerr << "read_blif failed" << std::endl;
  ntk = ntk_res;

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

std::tuple<uint32_t,uint32_t,float> abc_lutpack( lig_network& ntk, std::string benchmark )
{
  mockturtle::write_bench( ntk, "/tmp/mfsin2_"+benchmark+".bench" );
  //std::string command = fmt::format( "abc -q \"read_bench /tmp/mfsin_{}.bench; mfs -D 0 -dea; write_bench /tmp/mfsout_{}.bench;\"", benchmark, benchmark );
  std::string command = fmt::format( "abc -q \"read_bench /tmp/mfsin2_{}.bench; lutpack -N 3 -S 3 -L 200; time; write_blif /tmp/mfsin2_{}.blif; &get -mn; &ps;\"", benchmark, benchmark );

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

  lig_network ntk_res;

  std::string string_path = ("/tmp/mfsin2_"+benchmark+".blif");
  if( lorina::read_blif( string_path, blif_reader( ntk_res ) ) != lorina::return_code::success )
    std::cerr << "read_blif failed" << std::endl;
  ntk = ntk_res;

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
  static constexpr uint32_t K{6};
  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, double, bool> exp( "lig_exp_2", "benchmark", "a(init)", "d(init)", "a(new)", "d(new)", "t(new)", "eq(RS)" );

  for ( auto const& benchmark : epfl_benchmarks( ) )
  {
    //if( benchmark == "log2" || benchmark == "hyp" ) continue;
    fmt::print( "[i] processing {}\n", benchmark );
    auto path = "best_results/size/" + benchmark + "_sizen.blif";

    
    klut_network klut_olig;

    if ( lorina::read_blif( benchmark+"_dse8.blif", blif_reader( klut_olig ) ) != lorina::return_code::success )
    {
      printf("unsuccessful read 1\n");
      continue;
    }
    printf("|klut_olig|=%d\n", klut_olig.num_gates());

    lig_network lig0( klut_olig );
    depth_view<lig_network> lig0_d{ lig0 };

    printf("|lig0|=%d ", lig0.num_gates());

    lig_network lig1(klut_olig );

    std::string tmp0 = benchmark + "tmp0.bench";
    write_bench( lig0, tmp0 );
    klut_network klut0;
    if ( lorina::read_bench( tmp0,   bench_reader( klut0 ) ) != lorina::return_code::success )
    {
      printf("unsuccessful read 2\n");
      continue;
    }


    boptimizer_params rps;
    boptimizer_stats rst;
    rps.progress = true;
    // ps.pattern_filename = "1024sa1/" + benchmark + ".pat";
    rps.max_inserts = 30;
    rps.max_pis = 8;
    rps.max_divisors = 256;
    rps.verbose = false;
    uint32_t nGatesOld = lig1.num_gates()+1;
    printf("%d\n", lig1.num_gates());
    rps.max_trials = 1;

    auto beg = high_resolution_clock::now();

    uint32_t runtime_limit = 600;
    while( lig1.num_gates() < nGatesOld )
    {
      //abc_lutpack( lig1, benchmark );
      //lig1 = cleanup_dangling( lig1 );    
      //printf("LPK  %d\n", lig1.num_gates());
      //abc_mfs2( lig1, benchmark );
      //lig1 = cleanup_dangling( lig1 );    
      //printf("MFS2 %d\n", lig1.num_gates());
      //abc_mfs( lig1, benchmark );
      //lig1 = cleanup_dangling( lig1 );    
      //printf("MFS  %d\n", lig1.num_gates());

      rps.max_trials = 1;
      auto sto = high_resolution_clock::now();
      if( (duration_cast<seconds>(sto - beg)).count() > runtime_limit ) break;
      
      nGatesOld=lig1.num_gates();
      boptimize_klut<EX2, 6u, 6u>( lig1, rps, &rst );
      lig1 = cleanup_dangling( lig1 );    
      printf("P66L %d\n", lig1.num_gates());
      sto = high_resolution_clock::now();
      if( (duration_cast<seconds>(sto - beg)).count() > runtime_limit ) break;

      if( lig1.num_gates() == nGatesOld )
      {
        boptimize_klut<EX2, 8u, 6u>( lig1, rps, &rst );
        lig1 = cleanup_dangling( lig1 );    
        printf("P76L %d\n", lig1.num_gates());
        sto = high_resolution_clock::now();
        if( (duration_cast<seconds>(sto - beg)).count() > runtime_limit ) break;
      }

      if( lig1.num_gates() == nGatesOld )
      {
        rps.max_trials = 100;
        boptimize_klut<EX2, 8u, 6u>( lig1, rps, &rst );
        lig1 = cleanup_dangling( lig1 );    
        printf("P86H %d\n", lig1.num_gates());
        sto = high_resolution_clock::now();
        if( (duration_cast<seconds>(sto - beg)).count() > runtime_limit ) break;
      }
    }

    auto sto = high_resolution_clock::now();
    auto duration = (duration_cast<milliseconds>(sto - beg)).count()*0.001;;

//    }
    depth_view<lig_network> lig1_d{ lig1 };


    printf("|lig1|=%d  ", lig1.num_gates());

    std::string tmp1 = benchmark + experiments/P0E5.cpp"_dse9.blif";
    write_blif( lig1, tmp1 );

    const auto cec1 = true;//(benchmark == "hyp" || benchmark == "log2") ? true : abc_cec( lig1, benchmark );
    if( !cec1 )
    {
        printf("NOT EQUIVALENT\n");
    }

    exp( benchmark, lig0.num_gates(), lig0_d.depth(), lig1.num_gates(), lig1_d.depth(), duration, cec1 );
  }

  exp.save();
  exp.table();

  return 0;
}


//RUNTIME LIMIT AT 10 minutes
//|  benchmark | a(init) | d(init) | a(new) | d(new) | t(new) | eq(RS) |
//|      adder |     129 |     126 |    129 |    128 |   1.88 |   true |
//|        bar |     512 |       4 |    512 |      5 |   1.88 |   true |
//|        div |    3090 |    1100 |   3090 |   1103 |  30.68 |   true |
//|        hyp |   36836 |    4384 |  36814 |   4535 | 919.81 |   true |
//|       log2 |    6076 |     278 |   6043 |    252 | 627.60 |   true |
//|        max |     511 |     135 |    511 |    135 |   3.97 |   true |
//| multiplier |    4330 |     195 |   4316 |    211 | 595.46 |   true |
//|        sin |    1053 |      92 |   1025 |    112 | 244.22 |   true |
//|       sqrt |    2983 |    1526 |   2978 |   1180 | 159.69 |   true |
//|     square |    2959 |     172 |   2939 |    206 | 345.30 |   true |
//|    arbiter |     268 |      70 |    268 |     80 |   1.64 |   true |
//|      cavlc |      50 |       6 |     50 |      7 |   1.84 |   true |
//|       ctrl |      25 |       2 |     25 |      7 |   1.00 |   true |
//|        dec |     264 |       2 |    264 |      7 |   4.93 |   true |
//|        i2c |     177 |       9 |    176 |     10 |   2.98 |   true |
//|  int2float |      19 |       5 |     19 |      6 |   1.02 |   true |
//|   mem_ctrl |    1708 |      14 |   1696 |     14 |  39.42 |   true |
//|   priority |      93 |      30 |     92 |     30 |   2.99 |   true |
//|     router |      19 |       5 |     19 |      5 |   0.93 |   true |
//|      voter |    1180 |      30 |   1178 |     30 |  24.88 |   true |

// one round of aggressive 100 100 opt
//|  benchmark | a(init) | d(init) | a(new) | d(new) |  t(new) | eq(RS) |
//|      adder |     129 |     128 |    129 |    128 |    1.23 |   true |
//|        bar |     512 |       5 |    512 |      5 |    1.50 |   true |
//|        div |    3090 |    1103 |   3089 |   1105 |   46.17 |   true |
//|        hyp |   36814 |    4535 |  36776 |   4589 | 1283.50 |   true |
//|       log2 |    6043 |     252 |   6041 |    249 |  578.04 |   true |
//|        max |     511 |     135 |    511 |    135 |    9.36 |   true |
//| multiplier |    4316 |     211 |   4316 |    211 |  219.88 |   true |
//|        sin |    1025 |     112 |   1024 |    106 |   40.62 |   true |
//|       sqrt |    2978 |    1180 |   2978 |   1182 |  112.58 |   true |
//|     square |    2939 |     206 |   2939 |    204 |   57.50 |   true |
//|    arbiter |     268 |      80 |    268 |     74 |    1.64 |   true |
//|      cavlc |      50 |       7 |     50 |      6 |    3.79 |   true |
//|       ctrl |      25 |       7 |     25 |      9 |    0.68 |   true |
//|        dec |     264 |       7 |    264 |      8 |   19.37 |   true |
//|        i2c |     176 |      10 |    176 |      8 |    2.86 |   true |
//|  int2float |      19 |       6 |     19 |      6 |    1.38 |   true |
//|   mem_ctrl |    1696 |      14 |   1695 |     15 |   13.72 |   true |
//|   priority |      92 |      30 |     92 |     30 |    2.38 |   true |
//|     router |      19 |       5 |     19 |      5 |    1.12 |   true |
//|      voter |    1178 |      30 |   1177 |     30 |   52.35 |   true |

//un'altro round di 10 minuti. Questa volta 10 come base max_trials
//|  benchmark | a(init) | d(init) | a(new) | d(new) | t(new) | eq(RS) |
//|      adder |     129 |     128 |    129 |    128 |   1.83 |   true |
//|        bar |     512 |       5 |    512 |      4 |   1.88 |   true |
//|        div |    3089 |    1105 |   3089 |   1102 |  45.55 |   true |
//|        hyp |   36776 |    4589 |  36753 |   4592 | 940.64 |   true |
//|       log2 |    6041 |     249 |   6030 |    248 | 616.57 |   true |
//|        max |     511 |     135 |    511 |    135 |   3.76 |   true |
//| multiplier |    4316 |     211 |   4316 |    210 | 102.56 |   true |
//|        sin |    1024 |     106 |   1024 |    101 |  40.87 |   true |
//|       sqrt |    2978 |    1182 |   2976 |   1184 | 179.96 |   true |
//|     square |    2939 |     204 |   2939 |    205 |  49.75 |   true |
//|    arbiter |     268 |      74 |    268 |     78 |   1.94 |   true |
//|      cavlc |      50 |       6 |     50 |      6 |   1.84 |   true |
//|       ctrl |      25 |       9 |     25 |      8 |   1.01 |   true |
//|        dec |     264 |       8 |    264 |      7 |   4.83 |   true |
//|        i2c |     176 |       8 |    176 |     11 |   1.87 |   true |
//|  int2float |      19 |       6 |     19 |      6 |   1.14 |   true |
//|   mem_ctrl |    1695 |      15 |   1695 |     16 |  11.35 |   true |
//|   priority |      92 |      30 |     92 |     30 |   1.45 |   true |
//|     router |      19 |       5 |     19 |      5 |   1.06 |   true |
//|      voter |    1177 |      30 |   1177 |     29 |  21.97 |   true |