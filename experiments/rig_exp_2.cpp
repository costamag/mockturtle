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


template<class Ntk>
inline bool abc_mfs( Ntk const& ntk, std::string const& benchmark_fullpath )
{
  std::string command = fmt::format( "abc -q \"read_blif {}; mfs -ea {} /tmp/test.bench\"", benchmark_fullpath );

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
  while ( std::getline( ss, line, '\n' ) )
  {
    if ( line.size() >= 23u && line.substr( 0u, 23u ) == "Networks are equivalent" )
    {
      return true;
    }
  }

  return false;
}

int main()
{
  using namespace experiments;
  using namespace mockturtle;
  using namespace rils;
  static constexpr uint32_t K{6};
  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, double, bool> exp( "rig_exp_2", "benchmark", "rigs0", "rigs0_depth", "rigs1", "rigs1_depth", "t(RS)", "eq(RS)" );

  for ( auto const& benchmark : epfl_benchmarks( epfl & ~(experiments::div) ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    auto path = "benchmarks/best_results/size/" + benchmark + "_sizen.blif";

    klut_network klut_orig;

    if ( lorina::read_blif( path, blif_reader( klut_orig ) ) != lorina::return_code::success )
    {
      continue;
    }
    printf("|klut_orig|=%d\n", klut_orig.num_gates());

    rig_network rig0;

    if ( lorina::read_blif( path, blif_reader( rig0 ) ) != lorina::return_code::success )
    {
      printf("rig0 unsuccessful\n");
      continue;
    }
    rig0 = cleanup_dangling( rig0 );
    depth_view<rig_network> rig0_d{ rig0 };

    printf("|rig0|=%d ", rig0.num_gates());

    rig_network rig1;
    if ( lorina::read_blif( path, blif_reader( rig1 ) ) != lorina::return_code::success )
    {
      printf("rig1 unsuccessful\n");
      continue;
    }
    std::string tmp0 = benchmark + "tmp0.bench";
    write_bench( rig0, tmp0 );
    klut_network klut0;
    if ( lorina::read_bench( tmp0, bench_reader( klut0 ) ) != lorina::return_code::success )
    {
      continue;
    }
    printf("|klut0*|=%d\n", klut0.num_gates()-klut0.num_pos());

    resubstitution_params rps;
    resubstitution_stats rst;
    rps.progress = true;
    // ps.pattern_filename = "1024sa1/" + benchmark + ".pat";
    rps.max_inserts = 30;
    rps.max_trials = 100;
    rps.max_pis = 10;
    rps.max_divisors = std::numeric_limits<uint32_t>::max();
printf("a\n");
    rig_resubstitution<rils::network_t::kLUT, rils::support_selection_t::STRUCT_PIVOT, K>( rig1, rps, &rst );
printf("b\n");

    rig1 = cleanup_dangling( rig1 );
    depth_view<rig_network> rig1_d{ rig1 };


    printf("|rig1|=%d  ", rig1.num_gates());

    std::string tmp1 = benchmark + "_rig.blif";
    write_blif( rig1, tmp1 );
    klut_network klut1;
    if ( lorina::read_blif( tmp1, blif_reader( klut1 ) ) != lorina::return_code::success )
    {
      continue;
    }
    printf("|klut1|=%d\n", klut1.num_gates()-klut1.num_pos());
    printf("\n");
    const auto cec1 = true;//benchmark == "hyp" ? true : abc_cec( klut1, benchmark );

    exp( benchmark, rig0.num_gates(), rig0_d.depth(), rig1.num_gates(), rig1_d.depth(), to_seconds(rst.time_total), cec1 );
  }

  exp.save();
  exp.table();

  return 0;
}