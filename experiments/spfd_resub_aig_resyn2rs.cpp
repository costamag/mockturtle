/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
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
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/sim_resub.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/rewrite.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/topo_view.hpp>

#include <experiments.hpp>

template<class Ntk>
mockturtle::aig_network abc_opt( Ntk const& ntk, std::string str_code, std::string abc_script = "share; resyn2rs" )
{
  mockturtle::aig_network res;
  write_blif( ntk, "/tmp/pre" + str_code + ".blif" );

  std::string command = "abc -q \"r /tmp/pre" + str_code + ".blif; " + abc_script + "; write_aiger /tmp/pre" + str_code + ".aig\"";

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

  std::string string_path = ( "/tmp/pre" + str_code + ".aig" );
  if( lorina::read_aiger( string_path, mockturtle::aiger_reader( res ) ) != lorina::return_code::success )
    std::cerr << "read_aiger failed" << std::endl;

  return res;
}

std::string SCRIPT = "compress2rs; compress2rs";

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  double gain_rs{0};
  double gain_rw{0};
  double gain_spfd{0};
  double gain_bmatch{0};
  double cum_gain_rs{0};
  double cum_gain_rw{0};
  double cum_gain_spfd{0};
  double cum_gain_bmatch{0};

  experiment<std::string, uint32_t, float, float, float, float, float, float, float, float, bool, bool, bool, bool> exp( "spfd_aig", "benchmark", "size", "gain(RS)", "gain(RW)", "gain(BMATCH)", "gain(SPFD)", "time(RS)", "time(RW)", "time(BMATCH)", "time(SPFD)", "eq(RS)", "eq(RW)", "eq(BMATCH)", "eq(SPFD)" );

  double cnt{0};

  static constexpr uint32_t S = 10u;
  static constexpr uint32_t I = 10u;
  static constexpr uint32_t N = 100u;
  static constexpr uint32_t Ks = 10u;
  static constexpr uint32_t Kb = 10u;

  xag_npn_resynthesis<aig_network, aig_network, xag_npn_db_kind::aig_complete> resyn;
  exact_library_params eps;
  eps.np_classification = false;
  eps.compute_dc_classes = true;
  exact_library<aig_network, decltype( resyn )> exact_lib( resyn, eps );

  for ( auto const& benchmark : all_benchmarks( iscas | epfl | iwls ))
  {
    fmt::print( "[i] processing {}\n", benchmark );

    #pragma region RS
    
    aig_network aig_rs;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig_rs ) ) != lorina::return_code::success )
    {
      continue;
    }
    if( aig_rs.num_gates() > 300000 ) continue;

    resubstitution_params ps_rs;
    resubstitution_stats st_rs;

    // ps.pattern_filename = "1024sa1/" + benchmark + ".pat";
    ps_rs.max_inserts = 20;
    ps_rs.max_pis = Ks;
    ps_rs.max_trials = N;
    ps_rs.use_dont_cares = true;

    ps_rs.odc_levels = -1;
    ps_rs.max_divisors = std::numeric_limits<uint32_t>::max();
    //aig_rs = abc_opt( aig_rs, benchmark, SCRIPT );
    const uint32_t size_before = aig_rs.num_gates();
    sim_resubstitution( aig_rs, ps_rs, &st_rs );
    aig_rs = cleanup_dangling( aig_rs );

    const auto cec_rs = benchmark == "hyp" ? true : abc_cec( aig_rs, benchmark );
    
    #pragma endregion RS
    
    #pragma region RW
    
    aig_network aig_rw;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig_rw ) ) != lorina::return_code::success )
    {
      continue;
    }
    //aig_rw = abc_opt( aig_rw, benchmark, SCRIPT );

    rewrite_params ps_rw;
    rewrite_stats st_rw;
    ps_rw.use_dont_cares = true;

    //aig_rw = abc_opt( aig_rw, benchmark, SCRIPT );
    rewrite( aig_rw, exact_lib, ps_rw, &st_rw );
    aig_rw = cleanup_dangling( aig_rw );

    const auto cec_rw = benchmark == "hyp" ? true : abc_cec( aig_rw, benchmark );
    
    #pragma endregion RW
    
    printf("=================\n");

    #pragma region BMATCH
    aig_network aig_bmatch;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig_bmatch ) ) != lorina::return_code::success )
    {
      continue;
    }

    resubstitution_params ps_bmatch;
    resubstitution_stats st_bmatch;


    ps_bmatch.max_inserts = 20;
    ps_bmatch.max_pis = Ks;
    ps_bmatch.max_trials = N;
    ps_bmatch.progress = true;
    ps_bmatch.use_dont_cares = true;
    ps_bmatch.odc_levels = -1;

    ps_bmatch.max_divisors = std::numeric_limits<uint32_t>::max();

    //aig_bmatch = abc_opt( aig_bmatch, benchmark, SCRIPT );
    sim_resubstitution_spfd<Kb, S, I, true>( aig_bmatch, ps_bmatch, &st_bmatch );
    aig_bmatch = cleanup_dangling( aig_bmatch );

    const auto cec_bmatch = benchmark == "hyp" ? true : abc_cec( aig_bmatch, benchmark );
    #pragma endregion BMATCH

    #pragma region SPFD
    aig_network aig_spfd;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig_spfd ) ) != lorina::return_code::success )
    {
      continue;
    }

    resubstitution_params ps_spfd;
    resubstitution_stats st_spfd;


    ps_spfd.max_inserts = 20;
    ps_spfd.max_pis = Ks;
    ps_spfd.max_trials = N;
    ps_spfd.progress = true;
    ps_spfd.use_dont_cares = true;
    ps_spfd.max_divisors = std::numeric_limits<uint32_t>::max();

    //aig_spfd = abc_opt( aig_spfd, benchmark, SCRIPT );
    sim_resubstitution_spfd<Kb, S, I, false>( aig_spfd, ps_spfd, &st_spfd );
    aig_spfd = cleanup_dangling( aig_spfd );

    const auto cec_spfd = benchmark == "hyp" ? true : abc_cec( aig_spfd, benchmark );
    #pragma endregion SPFD

    cnt++;
    gain_rs = (double)(size_before - aig_rs.num_gates())/((double)size_before);
    gain_rw = (double)(size_before - aig_rw.num_gates())/((double)size_before);
    gain_spfd = (double)(size_before - aig_spfd.num_gates())/((double)size_before);
    gain_bmatch = (double)(size_before - aig_bmatch.num_gates())/((double)size_before);

    cum_gain_rs += gain_rs;
    cum_gain_rw += gain_rw;
    cum_gain_spfd += gain_spfd;
    cum_gain_bmatch += gain_bmatch;

    printf( "gain(RS)=%d gain(RW)=%d gain(BMATCH) = %d gain(SPFD)=%d\n", size_before - aig_rs.num_gates(), size_before - aig_rw.num_gates(), size_before - aig_bmatch.num_gates(), size_before - aig_spfd.num_gates() );

    exp( benchmark, size_before, 100*gain_rs, 100*gain_rw, 100*gain_bmatch, 100*gain_spfd, to_seconds( st_rs.time_total ), to_seconds( st_rw.time_total ), to_seconds( st_bmatch.time_total ), to_seconds( st_spfd.time_total ), cec_rs, cec_rw, cec_bmatch, cec_spfd );
    cnt+=1;
  }
  printf("<gain(RS)>=%.2f <gain>(RW)=%.2f <gain(BMATCH)>=%.2f <gain(SPFD)>=%.2f\n", 100*cum_gain_rs/cnt, 100*cum_gain_rw/cnt, 100*cum_gain_bmatch/cnt, 100*cum_gain_spfd/cnt );

  exp.save();
  exp.table();

  return 0;
}


