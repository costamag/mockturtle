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

#include <iostream>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/genlib.hpp>
#include <mockturtle/algorithms/emap.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/algorithms/cleanup.hpp>

#include <experiments.hpp>

using namespace experiments;
using namespace mockturtle;

std::pair<double, double> abc_map( aig_network const& aig, std::string const& library )
{
  write_aiger( aig, "/tmp/tmp.aig" );
  std::string command = fmt::format( "abc -q \"read /tmp/tmp.aig; read {}; &get; &nf -R 100; &put; print_stats;\"", library );

  std::array<char, 128> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  if ( !pipe )
  {
    throw std::runtime_error( "ABC: popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }

  /* parse the result */
  double area = -1, delay = -1;

  std::size_t pos = result.find( "area" );

  if ( pos != std::string::npos )
  {
    pos = result.find( "=", pos + 1 );
    std::string area_res = result.substr( pos + 1, result.find( " ", pos + 1 ) - pos - 1 );
    lorina::detail::trim( area_res );

    area = std::stod( area_res );

    pos = result.find( "=", pos + 1 );
    std::string delay_res = result.substr( pos + 1, result.find( " ", pos + 1 ) - pos - 1 );

    delay = std::stod( delay_res );
  }
  else
  {
    std::cout << "[e] failed to read the result\n";
  }

  return std::make_pair( area, delay );
}

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

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  fmt::print( "[i] processing technology library\n" );


  /* library to map to technology */
  std::vector<gate> gates;
  std::ifstream in( cell_libraries_path( "sky130" ) );

  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  tech_library_params tps;
  tech_library<5, classification_type::np_configurations> tech_lib( gates, tps );

  for ( auto const& benchmark : epfl_benchmarks( ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    std::vector<uint32_t> aig_size;
    std::vector<double> map_size;
    std::vector<uint32_t> aig_depth;
    std::vector<double> map_delay;
    std::vector<uint32_t> vheuristic;
    uint32_t heuristic=0;

    uint32_t num_old = aig.num_gates()+1;
    while( num_old > aig.num_gates() )
    {
        num_old = aig.num_gates();

        aig_network res2;

        aig = abc_opto( aig, benchmark, "rw" );
        aig = cleanup_dangling( aig );
        if( aig.num_gates() != num_old ){ vheuristic.push_back(0); }

        if( aig.num_gates() == num_old )
        {
            aig = abc_opto( aig, benchmark, "rs" );
            aig = cleanup_dangling( aig );
            if( aig.num_gates() != num_old ){ vheuristic.push_back(1); }

        }
        if( aig.num_gates() == num_old )
        {
            aig = abc_opto( aig, benchmark, "rf" );
            aig = cleanup_dangling( aig );
            if( aig.num_gates() != num_old ){ vheuristic.push_back(2); }

        }
        if( aig.num_gates() == num_old )
        {
            aig = abc_opto( aig, benchmark, "resyn2rs" );
            aig = cleanup_dangling( aig );
            if( aig.num_gates() != num_old ){ vheuristic.push_back(3); }

        }
        if( aig.num_gates() == num_old )
        {
            aig = abc_opto( aig, benchmark, "compress2rs" );
            aig = cleanup_dangling( aig );
            if( aig.num_gates() != num_old ){ vheuristic.push_back(4); }
        }

        if( aig.num_gates() == num_old )
        {
            vheuristic.push_back(5);
        }

        auto abc_res = abc_map( aig, cell_libraries_path( "sky130" ) );

        printf( "%d -> %f %f", aig.num_gates(), abc_res.first, abc_res.second );
        std::cout << std::endl;
        depth_view<aig_network> daig {aig};

        aig_size.push_back( aig.num_gates() );
        aig_depth.push_back( daig.depth() );
        map_size.push_back( abc_res.first );
        map_delay.push_back( abc_res.second );
    }

    write_aiger( aig, benchmark+"_optmap.aig" );

    printf("abc_aaig=np.array([");
    for( int i{0}; i< aig_size.size()-1; ++i )
    {
        if( i == aig_size.size()-2 )
            printf("%d])\n", aig_size[i]);
        else
            printf("%d, ", aig_size[i]);
    }

    printf("abc_amap=np.array([");
    for( int i{0}; i< map_size.size()-1; ++i )
    {
        if( i == map_size.size()-2 )
            printf("%f])\n", map_size[i]);
        else
            printf("%f, ", map_size[i]);
    }

    printf("abc_color=np.array([");
    for( int i{0}; i< vheuristic.size()-1; ++i )
    {
        if( i == vheuristic.size()-2 )
            printf("%d])\n", vheuristic[i]);
        else
            printf("%d, ", vheuristic[i]);
    }

    printf("d(aig)=[");
    for( auto x : aig_depth )
    {
        printf("%d, ", x);
    }
    printf("]\n");

    printf("d(map)=[");
    for( auto x : map_delay )
    {
        printf("%f, ", x);
    }
    printf("]\n");

    std::cout << std::endl;

    const auto cec = benchmark == "hyp" ? true : abc_cec( aig, benchmark );
    assert( cec && "[e] not equivalent" );

  }

  return 0;
}