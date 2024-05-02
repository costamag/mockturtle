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
Ntk abc_opto( Ntk const& ntk, std::string str_code, std::string abc_script = "resyn2rs" )
{
  write_aiger( ntk, "/tmp/"+str_code+".aig" );
  std::string command = "abc -q \"r /tmp/"+str_code+".aig; " + abc_script + "; write_aiger /tmp/" + str_code + ".aig\"";

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

std::tuple<uint32_t,uint32_t,float> abc_mfs( lig_network& ntk, std::string benchmark )
{//mfs2 -L 20 -ea
  mockturtle::write_bench( ntk, "/tmp/mfsin_"+benchmark+".bench" );
  //std::string command = fmt::format( "abc -q \"read_bench /tmp/mfsin_{}.bench; mfs -D 0 -dea; write_bench /tmp/mfsout_{}.bench;\"", benchmark, benchmark );
  std::string command = fmt::format( "abc -q \"read_bench /tmp/mfsin_{}.bench; mfs -e -M 0 -L 200; time; write_blif /tmp/mfsin_{}.blif; &get -mn; &ps;\"", benchmark, benchmark );

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

std::tuple<uint32_t,uint32_t,float> abc_mfs2( lig_network& ntk, std::string benchmark )
{
  mockturtle::write_bench( ntk, "/tmp/mfsin2_"+benchmark+".bench" );
  //std::string command = fmt::format( "abc -q \"read_bench /tmp/mfsin_{}.bench; mfs -D 0 -dea; write_bench /tmp/mfsout_{}.bench;\"", benchmark, benchmark );
  std::string command = fmt::format( "abc -q \"read_bench /tmp/mfsin2_{}.bench; mfs2 -e -M 0 -L 200; time; write_blif /tmp/mfsin2_{}.blif; &get -mn; &ps;\"", benchmark, benchmark );

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
klut_network abc_if( Ntk const& ntk, std::string str_code, uint32_t K=4u )
{
  write_aiger( ntk, "/tmp/"+str_code+".aig" );
  std::string command = "abc -q \"r /tmp/"+str_code+".aig; dch -f; if -a -K "+std::to_string(K)+"; write_blif /tmp/" + str_code + ".blif\"";

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

  klut_network res;

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

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, double, double, double, double, double> exp( "lig_exp_3", "benchmark", "a(map)", "a(mfs)", "a(mfs2)", "a(gre)", "a(exp1)", "a(exp3)", "d(map)", "d(mfs)", "d(mfs2)", "d(gre)", "d(exp1)", "d(exp3)", "t(mfs)", "t(mfs2)", "t(gre)", "t(exp1)", "t(exp3)" );

  static constexpr uint32_t K = 4;

  std::vector<double> aMAP;
  std::vector<double> aMFS;
  std::vector<double> aMFS2;
  std::vector<double> aGRE;
  std::vector<double> aEXP1;
  std::vector<double> aEXP3;

  std::vector<double> dMAP;
  std::vector<double> dMFS;
  std::vector<double> dMFS2;
  std::vector<double> dGRE;
  std::vector<double> dEXP1;
  std::vector<double> dEXP3;

  std::vector<double> tMFS;
  std::vector<double> tMFS2;
  std::vector<double> tGRE;
  std::vector<double> tEXP1;
  std::vector<double> tEXP3;

  for ( auto const& benchmark : iscas_benchmarks( ) )
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

    /* set up the parameters */
    boptimizer_params rps;
    rps.progress =false;
    rps.max_inserts = 20;
    rps.max_trials = 100;
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

    const auto lig_map = scopt::lut_map2( aig, ps, &st );

    if(lig_map.num_gates()>1000) continue;
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
    printf("aMAP : %6d\n", map_num_gates );
    uint32_t map_depth = lig_d.depth();

    /* MFS */
    lig_network lig_mfs;
    if ( lorina::read_blif( tmp, blif_reader( lig_mfs ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    if( lig_mfs.num_gates() != map_num_gates ) continue;
    auto mfs_res = abc_mfs( lig_mfs, benchmark );
    printf("MFS: %6d\n", lig_mfs.num_gates() );
    uint32_t mfs_num_gates = std::get<0>(mfs_res);
    uint32_t mfs_depth = std::get<1>(mfs_res);
    double mfs_time = 0;

    /* MFS2 */
    lig_network lig_mfs2;
    if ( lorina::read_blif( tmp, blif_reader( lig_mfs2 ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    if( lig_mfs2.num_gates() != map_num_gates ) continue;
    auto mfs2_res = abc_mfs2( lig_mfs2, benchmark );
    printf("MF2: %6d\n", lig_mfs2.num_gates() );
    uint32_t mfs2_num_gates = std::get<0>(mfs2_res);
    uint32_t mfs2_depth = std::get<1>(mfs2_res);
    double mfs2_time = 0;

    /* CASE 1 : GREEDY */    
    lig_network lig_gre;
    if ( lorina::read_blif( tmp, blif_reader( lig_gre ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    if( lig_gre.num_gates() != map_num_gates ) continue;

    boptimizer_stats rst_gre;
    boptimize_klut<GREEDY, K, K>( lig_gre, rps, &rst_gre );
    lig_gre = cleanup_dangling( lig_gre );
    printf("GRE: %6d\n", lig_gre.num_gates() );
    auto gre_res = abc_eval( lig_gre, benchmark );
    uint32_t gre_num_gates = std::get<0>(gre_res);
    uint32_t gre_depth = std::get<1>(gre_res);
    double gre_time = 0;

    bool const gre_cec = benchmark == "hyp" ? true : abc_cec( lig_gre, benchmark );
    if( !gre_cec ) 
    {
      printf("NEQ\n");
      continue;
    }

    /* CASE 3 : EXP1 1 */    
    lig_network lig_exp1;
    if ( lorina::read_blif( tmp, blif_reader( lig_exp1 ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    if( lig_exp1.num_gates() != map_num_gates ) continue;

    boptimizer_stats rst_exp1;

    boptimize_klut<EXP1, K, K>( lig_exp1, rps, &rst_exp1 );
    lig_exp1 = cleanup_dangling( lig_exp1 );
    printf("EXP1: %6d [%6d]\n", lig_exp1.num_gates(), lig_exp1.max_num_fanins);
    auto exp1_res = abc_eval( lig_exp1, benchmark );
    uint32_t exp1_num_gates = std::get<0>(exp1_res);
    uint32_t exp1_depth = std::get<1>(exp1_res);
    double exp1_time = 0;

    bool const exp1_cec = benchmark == "hyp" ? true : abc_cec( lig_exp1, benchmark );
    if( !exp1_cec ) 
    {
      printf("NEQ\n");
      continue;
    }

    lig_network lig_exp3;
    if ( lorina::read_blif( tmp, blif_reader( lig_exp3 ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    if( lig_exp3.num_gates() != map_num_gates ) continue;

    boptimizer_stats rst_exp3;

    boptimize_klut<GREEDY, K, K>( lig_exp3, rps, &rst_exp3 );
    lig_exp3 = cleanup_dangling( lig_exp3 );
    printf("EXP3: %6d [%6d]\n", lig_exp3.num_gates(), lig_exp3.max_num_fanins);
    auto exp3_res = abc_eval( lig_exp3, benchmark );
    uint32_t exp3_num_gates = std::get<0>(exp3_res);
    uint32_t exp3_depth = std::get<1>(exp3_res);
    double exp3_time = 0;

    bool const exp3_cec = benchmark == "hyp" ? true : abc_cec( lig_exp3, benchmark );
    if( !exp3_cec ) 
    {
      printf("NEQ\n");
      continue;
    }

    exp( benchmark, map_num_gates, mfs_num_gates, mfs2_num_gates, gre_num_gates, exp1_num_gates, exp3_num_gates, map_depth, mfs_depth, mfs2_depth, gre_depth, exp1_depth, exp3_depth, mfs_time, mfs2_time, gre_time, exp1_time, exp3_time );

    {
      aMAP.push_back( (double)map_num_gates );
      dMAP.push_back( (double)map_depth );

      aMFS.push_back( (double)mfs_num_gates );
      dMFS.push_back( (double)mfs_depth );
      tMFS.push_back( mfs_time );

      aMFS2.push_back( (double)mfs2_num_gates );
      dMFS2.push_back( (double)mfs2_depth );
      tMFS2.push_back( mfs2_time );

      aGRE.push_back( (double)gre_num_gates );
      dGRE.push_back( (double)gre_depth );
      tGRE.push_back( gre_time );

      aEXP1.push_back( (double)exp1_num_gates );
      dEXP1.push_back( (double)exp1_depth );
      tEXP1.push_back( exp1_time );

      aEXP3.push_back( (double)exp3_num_gates );
      dEXP3.push_back( (double)exp3_depth );
      tEXP3.push_back( exp3_time );

    }

    printf("\n");
  }

  double   avg_mf_g(0);
  double  avg_mf2_g(0);
  double  avg_gre_g(0);
  double avg_exp1_g(0);
  double avg_exp3_g(0);

  double   avg_mf_d(0);
  double  avg_mf2_d(0);
  double  avg_gre_d(0);
  double avg_exp1_d(0);
  double avg_exp3_d(0);

  double   avg_mf_t(0);
  double  avg_mf2_t(0);
  double  avg_gre_t(0);
  double avg_exp1_t(0);
  double avg_exp3_t(0);

  for( int i{0}; i<aMFS.size(); ++i )
  {
    avg_mf_g += ( - aMAP[i] + aMFS[i] ) /aMAP[i]/( (double)aMAP.size() );
   avg_mf2_g += ( - aMAP[i] + aMFS2[i] )/aMAP[i]/( (double)aMAP.size() );
   avg_gre_g += ( - aMAP[i] + aGRE[i] ) /aMAP[i]/( (double)aMAP.size() );
  avg_exp1_g += ( - aMAP[i] + aEXP1[i] )/aMAP[i]/( (double)aMAP.size() );
  avg_exp3_g += ( - aMAP[i] + aEXP3[i] )/aMAP[i]/( (double)aMAP.size() );
    
    avg_mf_d += ( - dMAP[i] + dMFS[i] ) /dMAP[i]/( (double)aMAP.size() );
   avg_mf2_d += ( - dMAP[i] + dMFS2[i] )/dMAP[i]/( (double)aMAP.size() );
   avg_gre_d += ( - dMAP[i] + dGRE[i] ) /dMAP[i]/( (double)aMAP.size() );
  avg_exp1_d += ( - dMAP[i] + dEXP1[i] )/dMAP[i]/( (double)aMAP.size() );
  avg_exp3_d += ( - dMAP[i] + dEXP3[i] )/dMAP[i]/( (double)aMAP.size() );

    avg_mf_t += ( tMFS[i] ) /( (double)aMAP.size() );
   avg_mf2_t += ( tMFS2[i] )/( (double)aMAP.size() );
   avg_gre_t += ( tGRE[i] ) /( (double)aMAP.size() );
  avg_exp1_t += ( tEXP1[i] )/( (double)aMAP.size() );
  avg_exp3_t += ( tEXP3[i] )/( (double)aMAP.size() );
  }

  printf( "<g(mfs) > : %f\n", avg_mf_g );
  printf( "<g(mfs2)> : %f\n", avg_mf2_g );
  printf( "<g(gre) > : %f\n", avg_gre_g );
  printf( "<g(exp1)> : %f\n", avg_exp1_g );
  printf( "<g(exp3)> : %f\n", avg_exp3_g );
  
  printf( "<d(mfs) > : %f\n",   avg_mf_d );
  printf( "<d(mfs2)> : %f\n",  avg_mf2_d );
  printf( "<d(gre) > : %f\n",  avg_gre_d );
  printf( "<d(exp1)> : %f\n", avg_exp1_d );
  printf( "<d(exp3)> : %f\n", avg_exp3_d );

  printf( "<t(mfs) > : %f\n",   avg_mf_t );
  printf( "<t(mfs2)> : %f\n",  avg_mf2_t );
  printf( "<t(gre) > : %f\n",  avg_gre_t );
  printf( "<t(exp1)> : %f\n", avg_exp1_t );
  printf( "<t(exp3)> : %f\n", avg_exp3_t );

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