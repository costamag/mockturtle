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
  std::string command = fmt::format( "abc -q \"read_bench /tmp/mfsin_{}.bench; mfs -e -W 20 -L 200; time; write_blif /tmp/mfsin_{}.blif; &get -mn; &ps;\"", benchmark, benchmark );

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
  std::string command = fmt::format( "abc -q \"read_bench /tmp/mfsin2_{}.bench; mfs2 -e -W 20 -L 200; time; write_blif /tmp/mfsin2_{}.blif; &get -mn; &ps;\"", benchmark, benchmark );

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

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, double, double, double, double> exp( "lig_exp_3", "benchmark", "a(map)", "a(mfs)", "a(mfs2)", "a(g,1)", "a(p,1)", "d(map)", "d(mfs)", "d(mfs2)", "d(g,1)", "d(p,1)", "t(mfs)", "t(mfs2)", "t(g,1)", "t(p,1)" );

  static constexpr uint32_t K = 4;

  std::vector<double> mp_areas;
  std::vector<double> mf_areas;
  std::vector<double> mf2_areas;
  std::vector<double> lp_areas;
  std::vector<double> p1_areas;

  std::vector<double> mp_depths;
  std::vector<double> mf_depths;
  std::vector<double> mf2_depths;
  std::vector<double> lp_depths;
  std::vector<double> p1_depths;

  std::vector<double> mf_times;
  std::vector<double> mf2_times;
  std::vector<double> lp_times;
  std::vector<double> p1_times;

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

    int nOld = aig.num_gates()+1;
//    while( aig.num_gates() < nOld )
//    {
//      nOld = aig.num_gates();
//      aig = abc_opto( aig, benchmark, "resyn2rs" );
//      printf("..%d\n", aig.num_gates());
//    }

    if( aig.num_gates() > 300000 ) continue;

    /* aig optimization */
    uint32_t nGates_old = aig.num_gates()+1;
  //  while( nGates_old > aig.num_gates() )
  //  {
  //    nGates_old = aig.num_gates();
  //    //aig = abc_opto( aig, benchmark );
  //    //printf("resyn2rs      : %d\n", aig.num_gates() );
  //    aig = abc_opto( aig, benchmark, "compress2rs" );
  //    printf("compress2rs   : %d\n", aig.num_gates() );
  //  }

    /* set up the parameters */
    boptimizer_params rps;
    rps.progress =true;
    rps.max_inserts = 20;
    rps.max_trials = 1;
    rps.max_pis = 16;
    rps.verbose = false;
    rps.max_divisors = 64;

    /* CASE 0: parameters */
    //klut_network klut=abc_if( aig, benchmark, 4u );

    scopt::lut_map2_params ps;
    ps.cut_enumeration_ps.cut_size = 4u;
    ps.cut_enumeration_ps.cut_limit = 8u;
    ps.recompute_cuts = true;
    ps.area_oriented_mapping = true;
    ps.cut_expansion = true;
    scopt::lut_map2_stats st;
    const auto lig0 = scopt::lut_map2( aig, ps, &st );


    write_blif( lig0, tmp );

    lig_network lig;
    if ( lorina::read_blif( tmp, blif_reader( lig ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    depth_view<lig_network> lig_d{ lig };
    auto const lig_cec = benchmark == "hyp" ? true : abc_cec( lig, benchmark );

    uint32_t map_num_gates = lig.num_gates();
    printf("MP : %6d\n", map_num_gates );
    uint32_t map_depth = lig_d.depth();

    /* MFS */
    lig_network lig_mfs;
    if ( lorina::read_blif( tmp, blif_reader( lig_mfs ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }

    if( lig_mfs.num_gates() != map_num_gates ) continue;

    int nOldMfs = lig_mfs.num_gates()+1;
    while ( lig_mfs.num_gates() < nOldMfs )
    {
      nOldMfs = lig_mfs.num_gates();
      auto mfs2_res = abc_mfs( lig_mfs, benchmark );
      printf("MFS: %6d\n", lig_mfs.num_gates() );
    }
    
    depth_view<lig_network> dlig_mfs{lig_mfs};

    uint32_t mfs_num_gates = lig_mfs.num_gates();//std::get<0>(mfs2_res);
    uint32_t mfs_depth = dlig_mfs.depth();
    double mfs_time = 0;//(double)std::get<2>(mfs2_res);


    /* MFS2 */
    lig_network lig_mfs2;
    if ( lorina::read_blif( tmp, blif_reader( lig_mfs2 ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    if( lig_mfs2.num_gates() != map_num_gates ) continue;
    int nOldMfs2 = lig_mfs2.num_gates()+1;
    while ( lig_mfs2.num_gates() < nOldMfs2 )
    {
      nOldMfs2 = lig_mfs2.num_gates();
      auto mfs2_res = abc_mfs2( lig_mfs2, benchmark );
      printf("MF2: %6d\n", lig_mfs2.num_gates() );
    }
    
    depth_view<lig_network> dlig_mfs2{lig_mfs2};

    uint32_t mfs2_num_gates = lig_mfs2.num_gates();//std::get<0>(mfs2_res);
    uint32_t mfs2_depth = dlig_mfs2.depth();
    double mfs2_time = 0;//(double)std::get<2>(mfs2_res);

    /* CASE 1 : GREEDY 1 */    
//    int nOld = aig.num_gates()+1;
//    while( aig.num_gates() < nOld )
//    {
//      nOld = aig.num_gates();
//      aig = abc_opto( aig, benchmark, "resyn2rs" );
//      printf("..%d\n", aig.num_gates());
//    }

    lig_network lig_lp;
    if ( lorina::read_blif( tmp, blif_reader( lig_lp ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    if( lig_lp.num_gates() != map_num_gates ) continue;
    int nOldLP = lig_lp.num_gates()+1;
    while ( lig_lp.num_gates() < nOldLP )
    {
      nOldLP = lig_lp.num_gates();
      auto lp_res = abc_lutpack( lig_lp, benchmark );
      printf("LPK: %6d\n", lig_lp.num_gates() );
    }
    
    depth_view<lig_network> dlig_lp{lig_lp};

    uint32_t lp_num_gates = lig_lp.num_gates();//std::get<0>(mfs2_res);
    uint32_t lp_depth = dlig_lp.depth();
    double lp_time = 0;//(double)std::get<2>(mfs2_res);

    /* CASE 3 : PIVOT 1 */    
    lig_network lig_p1;
    if ( lorina::read_blif( tmp, blif_reader( lig_p1 ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    if( lig_p1.num_gates() != map_num_gates ) continue;

    boptimizer_stats rst_p1;
    nOld = lig_p1.num_gates()+1;
    while( lig_p1.num_gates() < nOld )
    {
      nOld = lig_p1.num_gates();
      boptimize_klut<GREEDY, 4u, 4u>( lig_p1, rps, &rst_p1 );
      lig_p1 = cleanup_dangling( lig_p1 );
      printf("GRE[4,4]: %6d [%6d]\n", lig_p1.num_gates(), lig_p1.max_num_fanins);
      if( nOld == lig_p1.num_gates() )
      {
        boptimize_klut<GREEDY, 7u, 4u>( lig_p1, rps, &rst_p1 );
        lig_p1 = cleanup_dangling( lig_p1 );
        printf("GRE[7,4]: %6d [%6d]\n", lig_p1.num_gates(), lig_p1.max_num_fanins);
      }
      if( nOld == lig_p1.num_gates() )
      {
        boptimize_klut<PIVOT, 4u, 4u>( lig_p1, rps, &rst_p1 );
        lig_p1 = cleanup_dangling( lig_p1 );
        printf("PIV[4,4]: %6d [%6d]\n", lig_p1.num_gates(), lig_p1.max_num_fanins);
      }
      if( nOld == lig_p1.num_gates() )
      {
        boptimize_klut<PIVOT, 7u, 4u>( lig_p1, rps, &rst_p1 );
        lig_p1 = cleanup_dangling( lig_p1 );
        printf("PIV[7,4]: %6d [%6d]\n", lig_p1.num_gates(), lig_p1.max_num_fanins);
      }
    }

    depth_view<lig_network> lig_p1_d{ lig_p1 };

    //uint32_t lp_num_gates = lig_lp.num_gates();
    //uint32_t lp_depth = lig_lp_d.depth();
    //double p1_time = to_seconds( rst_p1.time_total );

    auto piv_res = abc_eval( lig_p1, benchmark );
    uint32_t p1_num_gates = std::get<0>(piv_res);
    uint32_t p1_depth = std::get<1>(piv_res);
    //double mfs_time = (double)std::get<2>(mfs_res);
    double p1_time = to_seconds( rst_p1.time_total );

    bool const lig_p1_cec = benchmark == "hyp" ? true : abc_cec( lig_p1, benchmark );
    if( !lig_p1_cec ) 
    {
      printf("NEQ\n");
      continue;
    }


    exp( benchmark, map_num_gates, mfs_num_gates, mfs2_num_gates, lp_num_gates, p1_num_gates, map_depth, mfs_depth, mfs2_depth, lp_depth, p1_depth, mfs_time, mfs2_time, lp_time, p1_time );

    if( !((map_num_gates == mfs2_num_gates)&&(map_num_gates == mfs_num_gates)&&( map_num_gates==lp_num_gates )&&( map_num_gates==p1_num_gates ) ) )
    {
      mp_areas.push_back( (double)map_num_gates );
      mp_depths.push_back( (double)map_depth );

      mf_areas.push_back( (double)mfs_num_gates );
      mf_depths.push_back( (double)mfs_depth );
      mf_times.push_back( mfs_time );

      mf2_areas.push_back( (double)mfs2_num_gates );
      mf2_depths.push_back( (double)mfs2_depth );
      mf2_times.push_back( mfs2_time );

      lp_areas.push_back( (double)lp_num_gates );
      lp_depths.push_back( (double)lp_depth );
      lp_times.push_back( lp_time );

      p1_areas.push_back( (double)p1_num_gates );
      p1_depths.push_back( (double)p1_depth );
      p1_times.push_back( p1_time );
    }

    printf("\n");
  }

  double avg_mf_g(0);
  double avg_mf2_g(0);
  double avg_lp_g(0);
  double avg_p1_g(0);

  double avg_mf_d(0);
  double avg_mf2_d(0);
  double avg_lp_d(0);
  double avg_p1_d(0);

  double avg_mf_t(0);
  double avg_mf2_t(0);
  double avg_lp_t(0);
  double avg_p1_t(0);

  for( int i{0}; i<mf_areas.size(); ++i )
  {
    avg_mf_g += ( - mp_areas[i] + mf_areas[i] )/mp_areas[i]/( (double)mf_areas.size() );
    avg_mf2_g += ( - mp_areas[i] + mf2_areas[i] )/mp_areas[i]/( (double)mf2_areas.size() );
    avg_lp_g += ( - mp_areas[i] + lp_areas[i] )/mp_areas[i]/( (double)mf_areas.size() );
    avg_p1_g += ( - mp_areas[i] + p1_areas[i] )/mp_areas[i]/( (double)mf_areas.size() );
    
    avg_mf_d += ( - mp_depths[i] + mf_depths[i] )/mp_depths[i]/( (double)mf_areas.size() );
    avg_mf2_d += ( - mp_depths[i] + mf2_depths[i] )/mp_depths[i]/( (double)mf2_areas.size() );
    avg_lp_d += ( - mp_depths[i] + lp_depths[i] )/mp_depths[i]/( (double)mf_areas.size() );
    avg_p1_d += ( - mp_depths[i] + p1_depths[i] )/mp_depths[i]/( (double)mf_areas.size() );

    avg_mf_t += mf_times[i]/( (double)mf_areas.size() );
    avg_mf2_t += mf2_times[i]/( (double)mf2_areas.size() );
    avg_lp_t += lp_times[i]/( (double)mf_areas.size() );
    avg_p1_t += p1_times[i]/( (double)mf_areas.size() );
  }

  printf( "<g(mfs)> : %f\n", avg_mf_g );
  printf( "<g(mfs2)>: %f\n", avg_mf2_g );
  printf( "<g(lp)>  : %f\n", avg_lp_g );
  printf( "<g(p1)>  : %f\n", avg_p1_g );
  
  printf( "<d(mfs)> : %f\n", avg_mf_d );
  printf( "<d(mfs2)>: %f\n", avg_mf2_d );
  printf( "<d(lp)>  : %f\n", avg_lp_d );
  printf( "<d(p1)>  : %f\n", avg_p1_d );

  printf( "<t(mfs)> : %f\n", avg_mf_t );
  printf( "<t(mfs2)>: %f\n", avg_mf2_t );
  printf( "<t(lp)>  : %f\n", avg_lp_t );
  printf( "<t(p1)>  : %f\n", avg_p1_t );

  exp.save();
  exp.table();

  return 0;
}