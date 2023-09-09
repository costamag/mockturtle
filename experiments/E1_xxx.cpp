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

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/algorithms/balancing.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <string>
#include <ctime>
#include <experiments.hpp>

template<class Ntk> void abc_map( Ntk const& );

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool> exp( "xxx_vs_sop_epfl", "benchmark", "s(sop)", "s(xxx)", "d(sop)", "d(xxx)", "cec(sop)", "cec(xxx)" );

  for ( auto const& benchmark : iscas_benchmarks())// ~(experiments::hyp) ) )//iscas_benchmarks( ) )//experiments::adder
  {
    fmt::print( "[i] processing {}\n", benchmark );
    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }
    depth_view daig{aig};

    lut_map_params ps_xxx;
    lut_map_stats st_xxx;
    ps_xxx.cut_enumeration_ps.cut_size = 7u;
    auto aig_xxx = xxx_balancing( aig, ps_xxx, &st_xxx );
    depth_view daig_xxx{aig_xxx};
    printf("xxx: ");
    abc_map(aig_xxx);
    const auto cec_xxx = abc_cec( daig_xxx, benchmark );

    lut_map_params ps_sop;
    lut_map_stats st_sop;
    ps_sop.cut_enumeration_ps.cut_size = 7u;
    auto aig_sop = sop_balancing( aig, ps_sop, &st_sop );
    depth_view daig_sop{aig_sop};
    printf("sop: ");
    abc_map(aig_sop);
    const auto cec_sop = abc_cec( daig_sop, benchmark );
    printf("           ======================           \n\n");

    exp( benchmark,
         aig_xxx.num_gates(), aig_sop.num_gates(),
         daig_xxx.depth(), daig_sop.depth(),
         cec_xxx, cec_sop
          );
  }

  exp.save();
  exp.table();

  return 0;
}

template<class Ntk>
void abc_map( Ntk const& ntk )
{
  using namespace mockturtle;

  aig_network res;
  write_blif( ntk, "/tmp/pre.blif" );

  std::string command = "abc -q \"read_library mcnc.genlib; r /tmp/pre.blif; st; dch; map -p; print_stats -p; print_stats -p;\"";
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
  printf("%s\n", result.c_str());
}
