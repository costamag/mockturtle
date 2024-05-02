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
#include <chrono>
#include <experiments.hpp>

using namespace mockturtle;
using namespace rils;

template<class Ntk>
Ntk abc_opto( Ntk const& ntk, std::string str_code, std::string abc_script = "resyn2rs" )
{
  write_aiger( ntk, "/tmp/"+str_code+".aig" );
  std::string command = "abc -q \"r /tmp/"+str_code+".aig; " + abc_script + "; fraig; write_aiger /tmp/" + str_code + ".aig\"";

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
//dch; if -a –C 12 -K 
void abc_satlut( lig_network& ntk, std::string benchmark )
{
  mockturtle::write_bench( ntk, "/tmp/satlut_"+benchmark+".bench" );
  //std::string command = fmt::format( "abc -q \"read_bench /tmp/mfsin_{}.bench; mfs -D 0 -dea; write_bench /tmp/mfsout_{}.bench;\"", benchmark, benchmark );
  std::string command = fmt::format( "abc -q \"read_bench /tmp/satlut_{}.bench; &get; &if -a - K 4; &satlut; &put; write_blif /tmp/satlut_{}.blif;\"", benchmark, benchmark );

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

  std::string string_path = ("/tmp/satlut_"+benchmark+".blif");
  if( lorina::read_blif( string_path, blif_reader( ntk_res ) ) != lorina::return_code::success )
    std::cerr << "read_blif failed" << std::endl;
  ntk = ntk_res;
}

std::tuple<uint32_t,uint32_t,float> abc_lutpack( lig_network& ntk, std::string benchmark )
{
  mockturtle::write_bench( ntk, "/tmp/mfsin2_"+benchmark+".bench" );
  //std::string command = fmt::format( "abc -q \"read_bench /tmp/mfsin_{}.bench; mfs -D 0 -dea; write_bench /tmp/mfsout_{}.bench;\"", benchmark, benchmark );
  std::string command = fmt::format( "abc -q \"read_bench /tmp/mfsin2_{}.bench; lutpack -L 200 -S 3; time; write_blif /tmp/mfsin2_{}.blif; &get -mn; &ps;\"", benchmark, benchmark );

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

template<class Ntk>
lig_network abc_if( Ntk const& ntk, std::string str_code, uint32_t K=4u )
{
  write_aiger( ntk, "/tmp/"+str_code+".aig" );
  std::string command = "abc -q \"r /tmp/"+str_code+".aig; fraig; st; dch; if -a –C 12 -K "+std::to_string(K)+"; write_blif /tmp/" + str_code + ".blif\"";

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

  lig_network res;

  std::string string_path = ("/tmp/"+str_code+".blif");
  if( lorina::read_blif( string_path, blif_reader( res ) ) != lorina::return_code::success )
    std::cerr << "read_blif failed" << std::endl;

  return res;
}

int main()
{
  using namespace experiments;
  using namespace mockturtle;
  using namespace rils;
  using namespace std::chrono;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, double, double> exp( "lig_exp_3", "benchmark", "a(map)", "a(flow1)", "a(flow2)", "d(flow1)", "d(flow2)", "t(flow1)", "t(flow2)" );

  static constexpr uint32_t K = 4;


  double  avg_flow1_g(0);
  double  avg_flow2_g(0);

  double N{0};
  double Tmax=600;

  for ( auto const& benchmark : all_benchmarks( ) )
  {
    N++;
    if( benchmark == "hyp" ) continue;
    fmt::print( "[i] processing {}\n", benchmark );

    std::string tmp = benchmark + "_exp1.blif";

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    if( aig.num_gates() > 300000 ) continue;
    aig = abc_opto( aig, benchmark );

    /* set up the parameters */
    boptimizer_params rps;
    rps.progress =false;
    rps.max_inserts = 20;
    rps.max_trials = 1;
    rps.max_pis = 16;
    rps.verbose = false;
    rps.max_divisors = 300;

    scopt::lut_map2_params ps;
    ps.cut_enumeration_ps.cut_size = K;
    ps.cut_enumeration_ps.cut_limit = 8u;
    ps.recompute_cuts = true;
    ps.area_oriented_mapping = true;
    ps.cut_expansion = true;
    scopt::lut_map2_stats st;

    //const auto lig_map = scopt::lut_map2( aig, ps, &st );
    const auto lig_map0 = abc_if( aig, benchmark, K );
    const auto lig_map = cleanup_ligs( lig_map0 );

    write_blif( lig_map, tmp );

    lig_network lig;
    if ( lorina::read_blif( tmp, blif_reader( lig ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    depth_view<lig_network> lig_d{ lig };
    auto const lig_cec = benchmark == "hyp" ? true : abc_cec( lig, benchmark );

    uint32_t map_num_gates = lig.num_gates();
    uint32_t map_depth = lig_d.depth();
    printf("MAP: %d %d\n", map_num_gates, map_depth );

    /* MFS2 */
    lig_network lig1;
    if ( lorina::read_blif( tmp, blif_reader( lig1 ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    if( lig1.num_gates() != map_num_gates ) continue;
    
    uint32_t num_gates_1 = lig1.num_gates()+1;

    auto beg1 = high_resolution_clock::now();
    double duration1{0};
    printf("ABC-FLOW\n");
    while( lig1.num_gates() < num_gates_1 && duration1 < Tmax )
    {
      printf(">>>>: %d \n", lig1.num_gates() );
      num_gates_1 = lig1.num_gates();
      abc_mfs2( lig1, benchmark, num_gates_1 );
      printf(".mfs: %d \n", lig1.num_gates() );
      abc_lutpack( lig1, benchmark );
      printf(".lpk: %d \n", lig1.num_gates() );
      //abc_satlut( lig1, benchmark );
      //printf(".slt: %d \n", lig1.num_gates() );
      lig1 = cleanup_dangling( lig1 );
      depth_view<lig_network> lig1_d{ lig1 };

      printf("FLOW1: %d %d\n", lig1.num_gates(), lig1_d.depth() );

      auto sto1 = high_resolution_clock::now();
      duration1 = (double)(duration_cast<seconds>(sto1 - beg1)).count()*1.0;

    }
    auto sto1 = high_resolution_clock::now();
    duration1 = (double)(duration_cast<seconds>(sto1 - beg1)).count()*1.0;
    depth_view<lig_network> lig1_d{ lig1 };

    /* FLOW 2 */    
    lig_network lig2;
    if ( lorina::read_blif( tmp, blif_reader( lig2 ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    if( lig2.num_gates() != map_num_gates ) continue;

    uint32_t num_gates_2 = lig2.num_gates()+1;
    printf("IRS-FLOW\n");
    auto beg2 = high_resolution_clock::now();
    double duration2{0};

    while( lig2.num_gates() < num_gates_2 && duration2 < Tmax )
    {
      num_gates_2 = lig2.num_gates();
      boptimizer_stats rst1;

      printf(">>>>: %d \n", lig2.num_gates() );
      boptimize_klut<EX2, K, K>( lig2, rps, &rst1 );
      lig2 = cleanup_dangling( lig2 );
      printf(".r44: %d \n", lig2.num_gates() );
      if( lig2.num_gates() == num_gates_2 )
      {
        boptimize_klut<EX2, K+2, K>( lig2, rps, &rst1 );
        lig2 = cleanup_dangling( lig2 );
        printf(".r64: %d \n", lig2.num_gates() );
      }

      abc_lutpack( lig2, benchmark );
      printf(".lpk: %d \n", lig2.num_gates() );
      //abc_satlut( lig2, benchmark );
      //printf(".slt: %d \n", lig2.num_gates() );
      lig2 = cleanup_dangling( lig2 );
      depth_view<lig_network> lig2_d{ lig2 };

      printf("FLOW2: %d %d\n", lig2.num_gates(), lig2_d.depth() );
      auto sto2 = high_resolution_clock::now();
      duration2 = (double)(duration_cast<seconds>(sto2 - beg2)).count()*1.0;
    }
    auto sto2 = high_resolution_clock::now();
    duration2 = (double)(duration_cast<seconds>(sto2 - beg2)).count()*1.0;

    depth_view<lig_network> lig2_d{ lig2 };

    bool cec = lig2.num_gates() < 9000 ? abc_cec( lig2, benchmark ) : true;
    if( !abc_cec( lig2, benchmark ) ) 
    {
      printf("NEQ\n");
      continue;
    }

    exp( benchmark, map_num_gates, lig1.num_gates(), lig2.num_gates(), lig1_d.depth(), lig2_d.depth(), duration1, duration2 );


    avg_flow1_g = avg_flow1_g*(N-1)/N+((double)lig1.num_gates()-(double)map_num_gates)/((double)map_num_gates)/N;
    avg_flow2_g = avg_flow2_g*(N-1)/N+((double)lig2.num_gates()-(double)map_num_gates)/((double)map_num_gates)/N;


    printf( "<g(flow1)> : %f\n", avg_flow1_g );
    printf( "<g(flow2)> : %f\n", avg_flow2_g );

    //printf("\n");
  }






  exp.save();
  exp.table();

  return 0;
}


//             div |   4335 |     4299 |     4262 |     2126 |     2126 |     3.00 |    50.00 |
//            log2 |   9777 |     9519 |     9474 |      149 |      283 |   142.00 |   689.00 |
//             max |   1006 |      772 |      780 |      250 |      242 |     3.00 |    16.00 |
//      multiplier |   7231 |     7147 |     7110 |      130 |      244 |     4.00 |   332.00 |
//             sin |   1853 |     1755 |     1706 |       78 |      119 |    28.00 |   167.00 |
//            sqrt |   6601 |     5207 |     4225 |     2028 |     2018 |    18.00 |   240.00 |
//          square |   5312 |     5233 |     4813 |      123 |      231 |     7.00 |   251.00 |
//         arbiter |   4139 |     4139 |     4011 |       30 |       30 |     6.00 |    42.00 |
//           cavlc |    288 |      250 |      223 |        8 |       11 |     3.00 |    40.00 |
//            ctrl |     48 |       43 |       38 |        4 |        8 |     1.00 |     3.00 |
//             i2c |    400 |      359 |      356 |        9 |        9 |     2.00 |     6.00 |
//       int2float |     88 |       68 |       69 |        7 |        7 |     1.00 |     5.00 |
//        mem_ctrl |  15940 |    11835 |     8617 |       59 |       49 |   176.00 |   609.00 |
//          router |     62 |       49 |       54 |        9 |       11 |     1.00 |     4.00 |
//           voter |   2475 |     2339 |     1981 |       19 |       29 |     9.00 |   107.00 |
//       ac97_ctrl |   3799 |     3403 |     3378 |        7 |        7 |     3.00 |    31.00 |
//        aes_core |   8350 |     7537 |     6918 |       16 |       28 |   132.00 |   614.00 |
//        des_area |   1831 |     1642 |     1662 |       20 |       24 |    16.00 |    40.00 |
//        des_perf |  25167 |    17501 |    19506 |       11 |       20 |   358.00 |   618.00 |
//             DMA |   7310 |     6954 |     6914 |       13 |       15 |    23.00 |    97.00 |
//             DSP |  14129 |    13020 |    12720 |       43 |       49 |    73.00 |   392.00 |
//        ethernet |  19943 |    19523 |    19551 |       18 |       19 |    39.00 |   162.00 |
//      iwls05_i2c |    350 |      329 |      319 |        9 |       11 |     2.00 |     7.00 |
// iwls05_mem_ctrl |   2763 |     2364 |     2430 |       20 |       19 |     7.00 |   121.00 |
//    pci_bridge32 |   5626 |     5512 |     5482 |       18 |       18 |    11.00 |   123.00 |
//            RISC |  20724 |    19176 |    18783 |       36 |       37 |    79.00 |   612.00 |
//            sasc |    190 |      182 |      181 |        5 |        5 |     0.00 |     2.00 |
//      simple_spi |    270 |      254 |      253 |        6 |       10 |     1.00 |     3.00 |
//             spi |   1175 |     1006 |      973 |       15 |       22 |     7.00 |    55.00 |
//          ss_pcm |    117 |      116 |      115 |        4 |        4 |     0.00 |     2.00 |
//      systemcaes |   2712 |     2521 |     2534 |       16 |       17 |     6.00 |    17.00 |
//      systemcdes |    899 |      624 |      634 |       13 |       24 |    11.00 |    87.00 |
//            tv80 |   2617 |     2286 |     2313 |       29 |       36 |    18.00 |    56.00 |
//       usb_funct |   4311 |     4049 |     3957 |       20 |       30 |     8.00 |   145.00 |
//         usb_phy |    139 |      128 |      123 |        5 |        7 |     1.00 |     4.00 |
//         vga_lcd |  31184 |    30944 |    30925 |       17 |       17 |    99.00 |   327.00 |
//       wb_conmax |  14719 |    13223 |    13440 |       13 |       13 |    70.00 |   602.00 |