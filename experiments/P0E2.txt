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

using namespace std::chrono;

using namespace mockturtle;
using namespace rils;

template<class Ntk>
Ntk abc_opto( Ntk const& ntk, std::string str_code, std::string abc_script = "resyn2rs" )
{
  write_aiger( ntk, "/tmp/"+str_code+".aig" );
  std::string command = "abc -q \"r /tmp/"+str_code+".aig; fraig; write_aiger /tmp/" + str_code + ".aig\"";
  //std::string command = "abc -q \"r /tmp/"+str_code+".aig; " + abc_script + "; fraig; write_aiger /tmp/" + str_code + ".aig\"";

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

std::tuple<uint32_t,uint32_t,float> abc_lutpack( lig_network& ntk, std::string benchmark )
{
  mockturtle::write_bench( ntk, "/tmp/mfsin2_"+benchmark+".bench" );
  //std::string command = fmt::format( "abc -q \"read_bench /tmp/mfsin_{}.bench; mfs -D 0 -dea; write_bench /tmp/mfsout_{}.bench;\"", benchmark, benchmark );
  std::string command = fmt::format( "abc -q \"read_bench /tmp/mfsin2_{}.bench; lutpack -L 200; time; write_blif /tmp/mfsin2_{}.blif; &get -mn; &ps;\"", benchmark, benchmark );

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

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t,double, double, double,  double, double, double, double, double, double, double> exp( "lig_exp_3", "benchmark", "a(map)", "a(mfs)", "a(mfs2)", "a(rnd)", "a(gre)", "a(piv1)", "a(piv2)", "a(piv3)", "a(exp1)", "a(exp2)", "a(exp3)", "t(mfs)", "t(mfs2)", "t(rnd)", "t(gre)", "t(piv1)", "t(piv2)", "t(piv3)", "t(exp1)", "t(exp2)", "t(exp3)" );

  static constexpr uint32_t K = 4;

  std::vector<double> aMAP;
  std::vector<double> aMFS;
  std::vector<double> aMF2;
  std::vector<double> aRND;
  std::vector<double> aGRE;
  std::vector<double> aEX1;
  std::vector<double> aEX2;
  std::vector<double> aEX3;
  std::vector<double> aPV1;
  std::vector<double> aPV2;
  std::vector<double> aPV3;

  std::vector<double> tMFS;
  std::vector<double> tMF2;
  std::vector<double> tRND;
  std::vector<double> tGRE;
  std::vector<double> tEX1;
  std::vector<double> tEX2;
  std::vector<double> tEX3;
  std::vector<double> tPV1;
  std::vector<double> tPV2;
  std::vector<double> tPV3;

  for ( auto const& benchmark : all_benchmarks() )
  {
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
    rps.progress = true;
    rps.max_inserts = 20;
    rps.max_trials = 1;
    rps.max_pis = 16;
    rps.verbose = false;
    rps.max_divisors = 300;

    scopt::lut_map2_params ps;
    ps.cut_enumeration_ps.cut_size = K;
    ps.cut_enumeration_ps.cut_limit = 8u;
    ps.required_delay = std::numeric_limits<uint32_t>::max();
    ps.area_oriented_mapping = true;
    ps.cut_expansion = true;
    scopt::lut_map2_stats st;

    //const auto lig_map0 = scopt::lut_map2( aig, ps, &st );
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
    printf("MAP: %6d\n", map_num_gates );
    uint32_t map_depth = lig_d.depth();

    /* MFS */
    lig_network lig_mfs;
    if ( lorina::read_blif( tmp, blif_reader( lig_mfs ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    if( lig_mfs.num_gates() != map_num_gates ) continue;

    auto beg_mfs = high_resolution_clock::now();
    //auto mfs_res = abc_mfs( lig_mfs, benchmark, map_num_gates );
    auto end_mfs = high_resolution_clock::now();
    auto duration_mfs = (duration_cast<milliseconds>(end_mfs - beg_mfs)).count()*0.001;

    uint32_t mfs_num_gates = 0;//std::get<0>(mfs_res);
    printf("MFS: %6d\n", mfs_num_gates );

    /* MFS2 */
    lig_network lig_mfs2;
    if ( lorina::read_blif( tmp, blif_reader( lig_mfs2 ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    if( lig_mfs2.num_gates() != map_num_gates ) continue;
    auto beg_mfs2 = high_resolution_clock::now();
    //auto mfs2_res = abc_mfs2( lig_mfs2, benchmark, map_num_gates );
    auto end_mfs2 = high_resolution_clock::now();
    auto duration_mfs2 = (duration_cast<milliseconds>(end_mfs2 - beg_mfs2)).count()*0.001;

    uint32_t mfs2_num_gates = 0;//std::get<0>(mfs2_res);
    printf("MF2: %6d\n", mfs2_num_gates );

    /* RANDOM */    
    lig_network lig_rand;
    if ( lorina::read_blif( tmp, blif_reader( lig_rand ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    if( lig_rand.num_gates() != map_num_gates ) continue;
    boptimizer_stats rst_rand;
    auto beg_rand = high_resolution_clock::now();
    boptimize_klut<RND, K, K>( lig_rand, rps, &rst_rand );
    auto end_rand = high_resolution_clock::now();
    auto duration_rand = (duration_cast<milliseconds>(end_rand - beg_rand)).count()*0.001;

    lig_rand = cleanup_ligs( lig_rand );
    auto rand_res = abc_eval( lig_rand, benchmark );
    uint32_t rand_num_gates = std::get<0>(rand_res);
    if( false && !abc_cec( lig_rand, benchmark ) ) 
    {
      printf("NEQ\n");
      continue;
    }
    printf("RND: %6d\n", lig_rand.num_gates() );

    /* GREEDY */    
    lig_network lig_gre;
    if ( lorina::read_blif( tmp, blif_reader( lig_gre ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    if( lig_gre.num_gates() != map_num_gates ) continue;
    boptimizer_stats rst_gre;
    auto beg_gre = high_resolution_clock::now();
    boptimize_klut<GRE, K, K>( lig_gre, rps, &rst_gre );
    auto end_gre = high_resolution_clock::now();
    auto duration_gre = (duration_cast<milliseconds>(end_gre - beg_gre)).count()*0.001;


    lig_gre = cleanup_ligs( lig_gre );
    auto gre_res = abc_eval( lig_gre, benchmark );
    uint32_t gre_num_gates = std::get<0>(gre_res);
    if( false && !abc_cec( lig_gre, benchmark ) ) 
    {
      printf("NEQ\n");
      continue;
    }
    printf("GRE: %6d\n", lig_gre.num_gates() );

    /* PIV1 */    
    lig_network lig_piv1;
    if ( lorina::read_blif( tmp, blif_reader( lig_piv1 ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    if( lig_piv1.num_gates() != map_num_gates ) continue;
    boptimizer_stats rst_piv1;
    auto beg_piv1 = high_resolution_clock::now();
    boptimize_klut<PV1, K, K>( lig_piv1, rps, &rst_piv1 );
    auto end_piv1 = high_resolution_clock::now();
    auto duration_piv1 = (duration_cast<milliseconds>(end_piv1 - beg_piv1)).count()*0.001;
    lig_piv1 = cleanup_ligs( lig_piv1 );
    auto piv1_res = abc_eval( lig_piv1, benchmark );
    uint32_t piv1_num_gates = std::get<0>(piv1_res);
    if( false && !abc_cec( lig_piv1, benchmark ) ) 
    {
      printf("NEQ\n");
      continue;
    }
    printf("PV1: %6d\n", lig_piv1.num_gates() );

    /* PIV2 */    
    lig_network lig_piv2;
    if ( lorina::read_blif( tmp, blif_reader( lig_piv2 ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    if( lig_piv2.num_gates() != map_num_gates ) continue;
    boptimizer_stats rst_piv2;

    auto beg_piv2 = high_resolution_clock::now();
    boptimize_klut<PV2, K, K>( lig_piv2, rps, &rst_piv2 );
    auto end_piv2 = high_resolution_clock::now();
    auto duration_piv2 = (duration_cast<milliseconds>(end_piv2 - beg_piv2)).count()*0.001;

    lig_piv2 = cleanup_ligs( lig_piv2 );
    auto piv2_res = abc_eval( lig_piv2, benchmark );
    uint32_t piv2_num_gates = std::get<0>(piv2_res);
    if( false && !abc_cec( lig_piv2, benchmark ) ) 
    {
      printf("NEQ\n");
      continue;
    }
    printf("PV2: %6d\n", lig_piv2.num_gates() );

    /* PIV3 */    
    lig_network lig_piv3;
    if ( lorina::read_blif( tmp, blif_reader( lig_piv3 ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    if( lig_piv3.num_gates() != map_num_gates ) continue;
    boptimizer_stats rst_piv3;
    auto beg_piv3 = high_resolution_clock::now();
    boptimize_klut<PV3, K, K>( lig_piv3, rps, &rst_piv3 );
    auto end_piv3 = high_resolution_clock::now();
    auto duration_piv3 = (duration_cast<milliseconds>(end_piv3 - beg_piv3)).count()*0.001;
    lig_piv3 = cleanup_ligs( lig_piv3 );
    auto piv3_res = abc_eval( lig_piv3, benchmark );
    uint32_t piv3_num_gates = std::get<0>(piv3_res);
    if( false && !abc_cec( lig_piv3, benchmark ) ) 
    {
      printf("NEQ\n");
      continue;
    }
    printf("PV3: %6d\n", lig_piv3.num_gates() );

    /* EXP1 */    
    lig_network lig_exp1;
    if ( lorina::read_blif( tmp, blif_reader( lig_exp1 ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    if( lig_exp1.num_gates() != map_num_gates ) continue;
    boptimizer_stats rst_exp1;
    
    auto beg_exp1 = high_resolution_clock::now();
    boptimize_klut<EX1, K, K>( lig_exp1, rps, &rst_exp1 );
    auto end_exp1 = high_resolution_clock::now();
    auto duration_exp1 = (duration_cast<milliseconds>(end_exp1 - beg_exp1)).count()*0.001;

    lig_exp1 = cleanup_ligs( lig_exp1 );
    auto exp1_res = abc_eval( lig_exp1, benchmark );
    uint32_t exp1_num_gates = std::get<0>(exp1_res);
    if( false && !abc_cec( lig_exp1, benchmark ) ) 
    {
      printf("NEQ\n");
      continue;
    }
    printf("EX1: %6d\n", lig_exp1.num_gates() );

    /* EXP2 */    
    lig_network lig_exp2;
    if ( lorina::read_blif( tmp, blif_reader( lig_exp2 ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    if( lig_exp2.num_gates() != map_num_gates ) continue;
    boptimizer_stats rst_exp2;
    auto beg_exp2 = high_resolution_clock::now();
    boptimize_klut<EX2, K, K>( lig_exp2, rps, &rst_exp2 );
    auto end_exp2 = high_resolution_clock::now();
    auto duration_exp2 = (duration_cast<milliseconds>(end_exp2 - beg_exp2)).count()*0.001;

    lig_exp2 = cleanup_ligs( lig_exp2 );
    auto exp2_res = abc_eval( lig_exp2, benchmark );
    uint32_t exp2_num_gates = std::get<0>(exp2_res);
    if( false && !abc_cec( lig_exp2, benchmark ) ) 
    {
      printf("NEQ\n");
      continue;
    }
    printf("EX2: %6d\n", lig_exp2.num_gates() );

    /* EXP3 */    
    lig_network lig_exp3;
    if ( lorina::read_blif( tmp, blif_reader( lig_exp3 ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    if( lig_exp3.num_gates() != map_num_gates ) continue;
    boptimizer_stats rst_exp3;
    auto beg_exp3 = high_resolution_clock::now();
    boptimize_klut<EX3, K, K>( lig_exp3, rps, &rst_exp3 );
    auto end_exp3 = high_resolution_clock::now();
    auto duration_exp3 = (duration_cast<milliseconds>(end_exp3 - beg_exp3)).count()*0.001;

    lig_exp3 = cleanup_ligs( lig_exp3 );
    auto exp3_res = abc_eval( lig_exp3, benchmark );
    uint32_t exp3_num_gates = std::get<0>(exp3_res);
    if( false && !abc_cec( lig_exp3, benchmark ) ) 
    {
      printf("NEQ\n");
      continue;
    }
    printf("EX3: %6d\n", lig_exp3.num_gates() );


    exp( benchmark, map_num_gates, mfs_num_gates, mfs2_num_gates, rand_num_gates, gre_num_gates, piv1_num_gates, piv2_num_gates, piv3_num_gates, exp1_num_gates, exp2_num_gates, exp3_num_gates,
    duration_mfs, duration_mfs2, duration_rand, duration_gre, duration_piv1, duration_piv2, duration_piv3, duration_exp1, duration_exp2, duration_exp3 );

//    if( (map_num_gates != mfs_num_gates)||(mfs_num_gates != mfs2_num_gates)||(mfs2_num_gates != rand_num_gates)||(rand_num_gates != gre_num_gates)||
//        (gre_num_gates != piv1_num_gates)||(piv1_num_gates != piv2_num_gates)||(piv1_num_gates != piv2_num_gates) )
    {
      aMAP.push_back( (double)map_num_gates );
      aMFS.push_back( (double)mfs_num_gates );
      aMF2.push_back( (double)mfs2_num_gates );
      aRND.push_back( (double)rand_num_gates );
      aGRE.push_back( (double)gre_num_gates );
      aEX1.push_back( (double)exp1_num_gates );
      aEX2.push_back( (double)exp2_num_gates );
      aEX3.push_back( (double)exp3_num_gates );
      aPV1.push_back( (double)piv1_num_gates );
      aPV2.push_back( (double)piv2_num_gates );
      aPV3.push_back( (double)piv3_num_gates );

      tMFS.push_back( duration_mfs );
      tMF2.push_back( duration_mfs2);
      tRND.push_back( duration_rand);
      tGRE.push_back( duration_gre );
      tEX1.push_back( duration_exp1);
      tEX2.push_back( duration_exp2);
      tEX3.push_back( duration_exp3);
      tPV1.push_back( duration_piv1);
      tPV2.push_back( duration_piv2);
      tPV3.push_back( duration_piv3);
    }

    printf("\n");
  }

  double  avg_mfs_g(0);
  double  avg_mf2_g(0);
  double  avg_gre_g(0);
  double  avg_rnd_g(0);
  double  avg_ex1_g(0);
  double  avg_ex2_g(0);
  double  avg_ex3_g(0);
  double  avg_pv1_g(0);
  double  avg_pv2_g(0);
  double  avg_pv3_g(0);

  double  avg_mfs_t(0);
  double  avg_mf2_t(0);
  double  avg_gre_t(0);
  double  avg_rnd_t(0);
  double  avg_ex1_t(0);
  double  avg_ex2_t(0);
  double  avg_ex3_t(0);
  double  avg_pv1_t(0);
  double  avg_pv2_t(0);
  double  avg_pv3_t(0);


  for( int i{0}; i<aMFS.size(); ++i )
  {
    avg_mfs_g += ( - aMAP[i] + aMFS[i] ) /aMAP[i]/( (double)aMAP.size() );
    avg_mf2_g += ( - aMAP[i] + aMF2[i] ) /aMAP[i]/( (double)aMAP.size() );
    avg_gre_g += ( - aMAP[i] + aGRE[i] ) /aMAP[i]/( (double)aMAP.size() );
    avg_rnd_g += ( - aMAP[i] + aRND[i] ) /aMAP[i]/( (double)aMAP.size() );
    avg_ex1_g += ( - aMAP[i] + aEX1[i] ) /aMAP[i]/( (double)aMAP.size() );
    avg_ex2_g += ( - aMAP[i] + aEX2[i] ) /aMAP[i]/( (double)aMAP.size() );
    avg_ex3_g += ( - aMAP[i] + aEX3[i] ) /aMAP[i]/( (double)aMAP.size() );
    avg_pv1_g += ( - aMAP[i] + aPV1[i] ) /aMAP[i]/( (double)aMAP.size() );
    avg_pv2_g += ( - aMAP[i] + aPV2[i] ) /aMAP[i]/( (double)aMAP.size() );
    avg_pv3_g += ( - aMAP[i] + aPV3[i] ) /aMAP[i]/( (double)aMAP.size() );

    avg_mfs_t += ( tMFS[i] ) / ( (double)aMAP.size() );
    avg_mf2_t += ( tMF2[i] ) / ( (double)aMAP.size() );
    avg_gre_t += ( tGRE[i] ) / ( (double)aMAP.size() );
    avg_rnd_t += ( tRND[i] ) / ( (double)aMAP.size() );
    avg_ex1_t += ( tEX1[i] ) / ( (double)aMAP.size() );
    avg_ex2_t += ( tEX2[i] ) / ( (double)aMAP.size() );
    avg_ex3_t += ( tEX3[i] ) / ( (double)aMAP.size() );
    avg_pv1_t += ( tPV1[i] ) / ( (double)aMAP.size() );
    avg_pv2_t += ( tPV2[i] ) / ( (double)aMAP.size() );
    avg_pv3_t += ( tPV3[i] ) / ( (double)aMAP.size() );

  }

  printf( "<g(mfs)> : %f\n", avg_mfs_g );
  printf( "<g(mf2)> : %f\n", avg_mf2_g );
  printf( "<g(rnd)> : %f\n", avg_rnd_g );
  printf( "<g(gre)> : %f\n", avg_gre_g );
  printf( "<g(pv1)> : %f\n", avg_pv1_g );
  printf( "<g(pv2)> : %f\n", avg_pv2_g );
  printf( "<g(pv3)> : %f\n", avg_pv3_g );
  printf( "<g(ex1)> : %f\n", avg_ex1_g );
  printf( "<g(ex2)> : %f\n", avg_ex2_g );
  printf( "<g(ex3)> : %f\n", avg_ex3_g );

  printf( "<t(mfs)> : %f\n", avg_mfs_t );
  printf( "<t(mf2)> : %f\n", avg_mf2_t );
  printf( "<t(rnd)> : %f\n", avg_rnd_t );
  printf( "<t(gre)> : %f\n", avg_gre_t );
  printf( "<t(pv1)> : %f\n", avg_pv1_t );
  printf( "<t(pv2)> : %f\n", avg_pv2_t );
  printf( "<t(pv3)> : %f\n", avg_pv3_t );
  printf( "<t(ex1)> : %f\n", avg_ex1_t );
  printf( "<t(ex2)> : %f\n", avg_ex2_t );
  printf( "<t(ex3)> : %f\n", avg_ex3_t );
  

  exp.save();
  exp.table();

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