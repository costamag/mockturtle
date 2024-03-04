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
#include <mockturtle/algorithms/lig_optimization.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/algorithms/rewrite.hpp>
#include <time.h>
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
  ntk_res._is_smart = ntk._is_smart;

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
  ntk_res._is_smart = ntk._is_smart;

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
  ntk_res._is_smart = ntk._is_smart;

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

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, double, double> exp( "lig_exp_2", "benchmark", "a(map)", "a(abc)", "a(new)", "d(map)", "d(abc)", "d(new)", "t(abc)", "t(new)" );

  static constexpr uint32_t K = 4;

  std::vector<double> mp_areas;
  std::vector<double> abc_areas;
  std::vector<double> dse_areas;

  std::vector<double> mp_depths;
  std::vector<double> abc_depths;
  std::vector<double> dse_depths;

  std::vector<double> abc_times;
  std::vector<double> dse_times;

  for ( auto const& benchmark : all_benchmarks( ) )
  {
    //if( benchmark == "hyp" ) continue;
    fmt::print( "[i] processing {}\n", benchmark );

    std::string tmp = benchmark + "_exp2.blif";

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    int nOld = aig.num_gates()+1;
    while( aig.num_gates() < nOld )
    {
      nOld = aig.num_gates();
      aig = abc_opto( aig, benchmark, "resyn2rs" );
      printf("[aig]%6d\n", aig.num_gates());
    }

    if( aig.num_gates() > 300000 ) continue;

    uint32_t nGates_old = aig.num_gates()+1;

    /* set up the parameters */
    lig_optimizer_params rps;
    rps.progress =true;
    rps.max_inserts = 20;
    rps.max_trials = 100;
    rps.max_pis = 16;
    rps.verbose = false;
    rps.max_divisors = 64;

    /* CASE 0: parameters */
    klut_network klut=abc_if( aig, benchmark, 4u );
    write_blif( klut, tmp );

    lig_network lig;
    lig._is_smart = true;
    if ( lorina::read_blif( tmp, blif_reader( lig ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    depth_view<lig_network> lig_d{ lig };
    auto const lig_cec = benchmark == "hyp" ? true : abc_cec( lig, benchmark );

    const uint32_t map_num_gates = lig.num_gates();
    printf("MP : %6d\n", map_num_gates );
    uint32_t map_depth = lig_d.depth();

    /* MFS */
    lig_network lig_abc;
    lig_abc._is_smart = true;
    if ( lorina::read_blif( tmp, blif_reader( lig_abc ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }

    if( lig_abc.num_gates() != map_num_gates ) 
    {
      printf("ERR ABC\n");
      continue;
    }
    int nOldMfs = lig_abc.num_gates()+1;

    clock_t t_abc;
    t_abc = clock();
    while ( lig_abc.num_gates() < nOldMfs )
    {
      nOldMfs = lig_abc.num_gates();
      abc_mfs( lig_abc, benchmark + "_mfs" );
      lig_abc = cleanup_dangling( lig_abc );
      printf("MFS: %6d\n", lig_abc.num_gates() );
      abc_mfs2( lig_abc, benchmark + "_mfs2" );
      lig_abc = cleanup_dangling( lig_abc );
      printf("MF2: %6d\n", lig_abc.num_gates() );
      abc_lutpack( lig_abc, benchmark + "_lpack" );
      lig_abc = cleanup_dangling( lig_abc );
      printf("LPK: %6d\n", lig_abc.num_gates() );
    }
    t_abc = clock() - t_abc;

    depth_view<lig_network> dlig_abc{lig_abc};

    auto abc_res = abc_eval( lig_abc, benchmark );
    uint32_t abc_num_gates = std::get<0>(abc_res);
    uint32_t abc_depth = std::get<1>(abc_res);
    double abc_time = ((double)t_abc)/CLOCKS_PER_SEC;


    /* CASE 3 : PIVOT 1 */    
    lig_network lig_dse;
    lig_dse._is_smart = true;
    if ( lorina::read_blif( tmp, blif_reader( lig_dse ) ) != lorina::return_code::success )
    {
      printf("lig unsuccessful\n");
      continue;
    }
    if( lig_dse.num_gates() != map_num_gates )
    {
      printf("ERR DSE\n");
      continue;
    }

    lig_optimizer_stats rst_dse;
    nOld = lig_dse.num_gates()+1;

    uint32_t OPTO{0};
    clock_t t_new;
    t_new = clock();
    while( lig_dse.num_gates() < nOld )
    {
      nOld = lig_dse.num_gates();
      if( OPTO == 0 || ( OPTO == 3 && lig_dse.num_gates() == nOld ) )
      {
        optimize_lig<GREEDY, 4u, 4u>( lig_dse, rps, &rst_dse );
        lig_dse = cleanup_dangling( lig_dse );
        OPTO = lig_dse.num_gates() == nOld ? 1 : 0;
        printf("GRE[4,4]: %6d [%6d]\n", lig_dse.num_gates(), lig_dse.max_num_fanins);
      }
      if( OPTO == 1 || ( OPTO == 0 && lig_dse.num_gates() == nOld ) )
      {
        optimize_lig<GREEDY, 7u, 4u>( lig_dse, rps, &rst_dse );
        lig_dse = cleanup_dangling( lig_dse );
        OPTO = lig_dse.num_gates() == nOld ? 2 : 0;
        printf("GRE[7,4]: %6d [%6d]\n", lig_dse.num_gates(), lig_dse.max_num_fanins);
      }
      if( OPTO == 2 || ( OPTO == 1 && lig_dse.num_gates() == nOld ) )
      {
        optimize_lig<PIVOT, 4u, 4u>( lig_dse, rps, &rst_dse );
        lig_dse = cleanup_dangling( lig_dse );
        OPTO = lig_dse.num_gates() == nOld ? 3 : 1;
        printf("PIV[4,4]: %6d [%6d]\n", lig_dse.num_gates(), lig_dse.max_num_fanins);
      }
      if( OPTO == 3 || ( OPTO == 2 && lig_dse.num_gates() == nOld ) )
      {
        optimize_lig<PIVOT, 7u, 4u>( lig_dse, rps, &rst_dse );
        lig_dse = cleanup_dangling( lig_dse );
        OPTO = lig_dse.num_gates() == nOld ? 3 : 2;
        printf("PIV[7,4]: %6d [%6d]\n", lig_dse.num_gates(), lig_dse.max_num_fanins);
      }
    }
    t_new = clock() - t_new;

    depth_view<lig_network> lig_dse_d{ lig_dse };

    auto piv_res = abc_eval( lig_dse, benchmark );
    uint32_t dse_num_gates = std::get<0>(piv_res);
    uint32_t dse_depth = std::get<1>(piv_res);
    //double mfs_time = (double)std::get<2>(mfs_res);
    double dse_time = ((double)t_new)/CLOCKS_PER_SEC;


    bool const lig_dse_cec = benchmark == "hyp" ? true : abc_cec( lig_dse, benchmark );
    if( !lig_dse_cec ) 
    {
      printf("NEQ\n");
      continue;
    }


    exp( benchmark, map_num_gates, abc_num_gates, dse_num_gates, map_depth, abc_depth, dse_depth, abc_time, dse_time );

    if( !((map_num_gates == abc_num_gates)&&( map_num_gates==dse_num_gates ) ) )
    {
      mp_areas.push_back( (double)map_num_gates );
      mp_depths.push_back( (double)map_depth );

      abc_areas.push_back( (double)abc_num_gates );
      abc_depths.push_back( (double)abc_depth );
      abc_times.push_back( abc_time );

      dse_areas.push_back( (double)dse_num_gates );
      dse_depths.push_back( (double)dse_depth );
      dse_times.push_back( dse_time );
    }

    printf("\n");
  }

  double avg_abc_g(0);
  double avg_dse_g(0);

  double avg_abc_d(0);
  double avg_dse_d(0);

  double avg_abc_t(0);
  double avg_dse_t(0);

  for( int i{0}; i<mp_areas.size(); ++i )
  {
    avg_abc_g += ( - mp_areas[i] + abc_areas[i] )/mp_areas[i]/( (double)abc_areas.size() );
    avg_dse_g += ( - mp_areas[i] + dse_areas[i] )/mp_areas[i]/( (double)mp_areas.size() );
    
    avg_abc_d += ( - mp_depths[i] + abc_depths[i] )/mp_depths[i]/( (double)abc_areas.size() );
    avg_dse_d += ( - mp_depths[i] + dse_depths[i] )/mp_depths[i]/( (double)mp_areas.size() );

    avg_abc_t += abc_times[i]/( (double)abc_areas.size() );
    avg_dse_t += dse_times[i]/( (double)mp_areas.size() );
  }

  printf( "<g(abc)> : %f\n", avg_abc_g );
  printf( "<g(dse)>  : %f\n", avg_dse_g );
  
  printf( "<d(abc)> : %f\n", avg_abc_d );
  printf( "<d(dse)>  : %f\n", avg_dse_d );

  printf( "<t(abc)> : %f\n", avg_abc_t );
  printf( "<t(dse)>  : %f\n", avg_dse_t );

  exp.save();
  exp.table();

  return 0;
}