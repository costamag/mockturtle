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
#include <mockturtle/io/write_bench.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/mapping_view.hpp>
#include <mockturtle/algorithms/lig_optimization.hpp>
#include <time.h>


#include <experiments.hpp>


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
  static constexpr uint32_t K{6};
  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, double, uint32_t, uint32_t, double, bool> exp( "lig_exp_2", "benchmark", "a(init)", "d(init)", "a(mfs)", "d(mfs)", "t(mfd)", "a(new)", "d(new)", "t(new)", "eq(RS)" );

  for ( auto const& benchmark : epfl_benchmarks( ) )
  {
    //if( benchmark == "log2" ) continue;
    fmt::print( "[i] processing {}\n", benchmark );
    auto path = "benchmarks/best_results/size/" + benchmark + "_sizen.blif";

    klut_network klut_olig;

    if ( lorina::read_blif( path, blif_reader( klut_olig ) ) != lorina::return_code::success )
    {
      printf("unsuccessful read \n");
      continue;
    }
    printf("|klut_olig|=%d\n", klut_olig.num_gates());

    lig_network lig0( klut_olig );
    depth_view<lig_network> lig0_d{ lig0 };

    printf("|lig0|=%d ", lig0.num_gates());

    lig_network lig1(klut_olig );
    lig1._is_smart=true;

    std::string tmp0 = benchmark + "tmp0.bench";
    write_bench( lig0, tmp0 );
    klut_network klut0;
    if ( lorina::read_bench( tmp0, bench_reader( klut0 ) ) != lorina::return_code::success )
    {
      printf("unsuccessful read\n");
      continue;
    }

    lig_optimizer_params rps;
    lig_optimizer_stats rst;
    rps.progress = true;
    // ps.pattern_filename = "1024sa1/" + benchmark + ".pat";
    rps.max_inserts = 50;
    rps.max_trials = 100;
    rps.max_pis = 20;
    rps.max_divisors = 256;
    rps.verbose = false;
    uint32_t nGatesOld = lig1.num_gates()+1;
    printf("%d\n", lig1.num_gates());

    clock_t t_new;
    t_new = clock();

    while( lig1.num_gates() < nGatesOld )
    {
      nGatesOld=lig1.num_gates();
      optimize_lig<GREEDY, 6u, 6u>( lig1, rps, &rst );
      lig1 = cleanup_dangling( lig1 );    
      printf("%d\n", lig1.num_gates());
//      if( lig1.num_gates() == nGatesOld )
//      {
//        nGatesOld=lig1.num_gates();
//        optimize_lig<GREEDY, 10u, 6u>( lig1, rps, &rst );
//        lig1 = cleanup_dangling( lig1 );    
//        printf("%d\n", lig1.num_gates());
//      }
//      if( lig1.num_gates() == nGatesOld )
//      {
//        nGatesOld=lig1.num_gates();
//        optimize_lig<PIVOT, 6u, 6u>( lig1, rps, &rst );
//        lig1 = cleanup_dangling( lig1 );    
//        printf("%d\n", lig1.num_gates());
//      }
//      if( lig1.num_gates() == nGatesOld )
//      {
//        nGatesOld=lig1.num_gates();
//        optimize_lig<PIVOT, 11u, 6u>( lig1, rps, &rst );
//        lig1 = cleanup_dangling( lig1 );    
//        printf("%d\n", lig1.num_gates());
//      }
    }
    t_new = clock() - t_new;

    depth_view<lig_network> lig1_d{ lig1 };


    printf("|lig1|=%d  ", lig1.num_gates());

    std::string tmp1 = benchmark + "_lig.blif";
    write_blif( lig1, tmp1 );
    klut_network klut1;
    if ( lorina::read_blif( tmp1, blif_reader( klut1 ) ) != lorina::return_code::success )
    {
      continue;
    }
    printf("|klut1|=%d\n", klut1.num_gates()-klut1.num_pos());
    printf("\n");
    const auto cec1 = true;//benchmark == "hyp" ? true : abc_cec( klut1, benchmark );

    auto mfs_res = abc_mfs( klut0, benchmark );

    exp( benchmark, lig0.num_gates(), lig0_d.depth(), std::get<0>(mfs_res), std::get<1>(mfs_res), std::get<2>(mfs_res), lig1.num_gates(), lig1_d.depth(), ((double)t_new)/CLOCKS_PER_SEC, cec1 );
  }

  exp.save();
  exp.table();

  return 0;
}


