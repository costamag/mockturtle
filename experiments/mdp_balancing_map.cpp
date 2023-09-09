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
 * included in all copies or substantial mdpions of the Software.
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
#include <mockturtle/algorithms/balancing/mdp_balancing.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <string>

#include <experiments.hpp>

using namespace mockturtle;

template<class Ntk> void abc_map( Ntk const& );
template<class Ntk> aig_network abc_sopbalancing( Ntk );

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, double, bool> exp( "mdp_balancing_map", "benchmark", "size", "depth", "size 4", "depth 4", "RT 4", "cec 4" );

  mdp_rebalancing<aig_network> mdp_balancing;

  for ( auto const& benchmark : epfl_benchmarks(~experiments::hyp))//epfl_benchmarks(~(experiments::hyp | experiments::div | experiments::adder | experiments::bar | experiments::log2 | experiments::max | experiments::multiplier | experiments::sin | experiments::sqrt | experiments::div | experiments::square ) ))//
  {
    fmt::print( "[i] processing {}\n", benchmark );
    aig_network xaig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xaig ) ) != lorina::return_code::success )
    {
      continue;
    }
    depth_view dxaig{xaig};

    balancing_params ps;
    balancing_stats st4;

    ps.progress = true;
    ps.only_on_critical_path = true;
    ps.cut_enumeration_ps.cut_size = 4u;

    auto xaig4 = abc_sopbalancing(xaig);
    depth_view dxaig4{xaig4};

    uint32_t depthNew = dxaig4.depth();
    uint32_t depthOld = dxaig.depth();
    
    printf("d(XAIG)=%d s(XAIG)=%d\n", depthNew, dxaig4.num_gates());

    uint32_t sizeBest{0xFFFFFFFF};
    depth_view<aig_network> dxaigBest {dxaig4};
    
    while( depthNew < depthOld )
    {

      xaig4 = abc_sopbalancing(xaig4);
      depth_view dxaigTmp0{xaig4}; 
      if( dxaigTmp0.depth() == depthNew )
      {
        xaig4 = balancing( xaig4, { mdp_balancing }, ps, &st4 );
      }
      depth_view dxaigTmp{xaig4}; 
      dxaig4 = dxaigTmp;
      depthOld = depthNew;
      depthNew = dxaigTmp.depth();
      printf("d(XAIG)=%d s(XAIG)=%d\n", depthNew, dxaig4.num_gates());

      if( depthNew < depthOld || ((depthNew < depthOld )&&(dxaig4.num_gates() < sizeBest) ))
      {
        sizeBest = dxaig4.num_gates();
        dxaigBest = dxaig4;
      }
    }
    abc_map( dxaigBest );

    ps.cut_enumeration_ps.cut_size = 4u;

    const auto cec4 = abc_cec( dxaigBest, benchmark );

    exp( benchmark,
         xaig.num_gates(), dxaig.depth(),
         dxaigBest.num_gates(), dxaigBest.depth(),
         to_seconds( st4.time_total ), cec4
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

  std::string command = "abc -q \"read_library mcnc.genlib; r /tmp/pre.blif; st; dch; map -p; print_stats;\"";
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


template<class Ntk>
aig_network abc_sopbalancing( Ntk ntk )
{
  aig_network res;
  write_blif( ntk, "/tmp/pre.blif" );

  std::string command = "abc -q \"r /tmp/pre.blif; if -g -K 6 -C 8; write_aiger /tmp/pre.aig\"";

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

  std::string string_path = ( "/tmp/pre.aig" );
  if( lorina::read_aiger( string_path, aiger_reader( res ) ) != lorina::return_code::success )
    std::cerr << "read_aiger failed" << std::endl;

  return res;
}