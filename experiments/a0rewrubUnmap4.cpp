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
#include <mockturtle/algorithms/emap2.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/scg.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/boptimizer2.hpp>
#include <ctime>
#include <experiments.hpp>

using namespace experiments;
using namespace mockturtle;
using namespace std::chrono;

std::pair<double, double> abc_map( aig_network const& aig, std::string const& library, uint32_t nUnmap )
{
  write_aiger( aig, "/tmp/tmp.aig" );//&get; map; &put

  std::string unmap_string = "";
  for( int i{0}; i<nUnmap; ++i )
  {
    unmap_string += "unmap; mfs2 -a; strash; dfraig; resyn2rs; dfraig; compress2rs; dch; map -a; ";
  }

  std::string command = fmt::format( "abc -q \"read /tmp/tmp.aig; read {}; dch; map -a;" + unmap_string + "  print_stats;\"", library );

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

aig_network abc_if( aig_network const& ntk, std::string str_code, uint32_t K=4u )
{
  write_aiger( ntk, "/tmp/"+str_code+".aig" );
  std::string command = "abc -q \"r /tmp/"+str_code+".aig; dch -f; if -g; strash; dfraig; write_aiger /tmp/" + str_code + ".aig\"";

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

  aig_network res;

  std::string string_path = ("/tmp/"+str_code+".aig");
  if( lorina::read_aiger( string_path, aiger_reader( res ) ) != lorina::return_code::success )
    std::cerr << "read_blif failed" << std::endl;

  return res;
}

template<class Ntk>
Ntk abc_opto( Ntk const& ntk, std::string str_code, std::string abc_script = "resyn2rs" )
{
  write_aiger( ntk, "/tmp/"+str_code+".aig" );
  std::string command = "abc -q \"r /tmp/"+str_code+".aig; dfraig; " + abc_script + "; write_aiger /tmp/" + str_code + ".aig\"";

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

template<class Ntk>
Ntk abc_strash( Ntk const& ntk, std::string str_code )
{
  write_aiger( ntk, "/tmp/"+str_code+".aig" );
  std::string command = "abc -q \"r /tmp/"+str_code+".aig; mfs2 -a; mfs2 -a; mfs2 -a; strash; write_aiger /tmp/" + str_code + ".aig\"";

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

template<class Ntk>
Ntk abc_strash1( Ntk const& ntk, std::string str_code )
{
  write_aiger( ntk, "/tmp/"+str_code+".aig" );
  std::string command = "abc -q \"r /tmp/"+str_code+".aig; strash; write_aiger /tmp/" + str_code + ".aig\"";

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

  experiment<std::string, double, double, double, double, bool, bool> exp( "SCOPTA", "benchmark", "a(mfs)", "a(pmo)", "d(mfs)", "d(pmo)", "eq(mfs)", "eq(pmo)");

  fmt::print( "[i] processing technology library\n" );


  /* library to map to technology */
  std::vector<gate> gates;
  std::ifstream in( cell_libraries_path( "asap7" ) );//sky130

  pLibrary_t database( "asap7" );

  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  tech_library_params tps;
  tech_library<5, classification_type::np_configurations> tech_lib( gates, tps );


  double N{1};
  double rarea1{0};
  double rarea2{0};
  double rdept1{0};
  double rdept2{0};
  for ( auto const& benchmark : all_benchmarks( epfl & ~(experiments::ctrl | experiments::router) ) )
  {

    fmt::print( "[i] processing {}\n", benchmark );

    bool start=true;
    bool close=false;

    aig_network aig;
    if ( lorina::read_aiger( experiments::benchmark_path(benchmark), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    if( aig.num_gates()>300000 || benchmark == "hyp"  ) continue;

    aig = abc_opto( aig, benchmark, "resyn2rs" );
    aig = cleanup_dangling( aig );
    aig = abc_opto( aig, benchmark, "compress2rs" );
    aig = cleanup_dangling( aig );
    aig = abc_opto( aig, benchmark, "resyn2rs" );
    aig = cleanup_dangling( aig );
    aig = abc_opto( aig, benchmark, "compress2rs" );
    aig = cleanup_dangling( aig );
    aig = abc_opto( aig, benchmark, "resyn2rs" );
    aig = cleanup_dangling( aig );
    aig = abc_opto( aig, benchmark, "compress2rs" );
    aig = cleanup_dangling( aig );

    const auto cec = benchmark == "hyp" ? true : abc_cec( aig, benchmark );
    if(!cec)
      printf("ERROR\n");
    assert( cec && "[e] not equivalent" );

    uint32_t nUnmap{1};

    scopt::emap2_params ps;
    ps.required_time = std::numeric_limits<float>::max();
    ps.area_oriented_mapping = true;
    scopt::emap2_stats st;

    scopt::scg_network scg1 = emap2_klut( aig, tech_lib, ps, &st );
    scg1 = cleanup_scg( scg1 );

    scopt::scg_network scg2 = emap2_klut( aig, tech_lib, ps, &st );
    scg2 = cleanup_scg( scg2 );

    printf("A0[1]=%f\n", scg1.compute_area());
    aig_network dump1;
    for( int i{0}; i<nUnmap; ++i )
    {
      dump1 = scg1.unmap();
      dump1 = abc_strash( dump1, benchmark );
      dump1 = abc_opto( dump1, benchmark, "resyn2rs" );
      dump1 = abc_opto( dump1, benchmark, "compress2rs" );
      scg1 = emap2_klut( dump1, tech_lib, ps, &st );
      scg1 = cleanup_scg( scg1 );
    }

    rewrub_sc_params rps;
    rewrub_sc_stats rst_p1;
    printf("A0[2]=%f\n", scg2.compute_area());
    aig_network dump2;
    for( int i{0}; i<nUnmap; ++i )
    {
      dump2 = scg2.unmap();
      dump2 = abc_strash1( dump2, benchmark );
      dump2 = abc_opto( dump2, benchmark, "resyn2rs" );
      dump2 = abc_opto( dump2, benchmark, "compress2rs" );
      scg2 = emap2_klut( dump2, tech_lib, ps, &st );
      scg2 = cleanup_scg( scg2 );
      rewrub_sc( scg2, database, rps, &rst_p1 );
      rewrub_sc( scg2, database, rps, &rst_p1 );
      rewrub_sc( scg2, database, rps, &rst_p1 );
    }

    const auto cec1 = benchmark == "hyp" ? true : abc_cec( scg1, benchmark );
    if(!cec1)
      printf("ERROR 1\n");
    
    const auto cec2 = benchmark == "hyp" ? true : abc_cec( scg2, benchmark );
    if(!cec2)
      printf("ERROR 2\n");

    printf("a(mfs)=%6f\n", scg1.compute_area() );
    printf("a(pmo)=%6f\n", scg2.compute_area() );
    printf("d(mfs)=%6f\n", scg1.compute_worst_delay() );
    printf("d(pmo)=%6f\n", scg2.compute_worst_delay() );

    exp( benchmark, scg1.compute_area(), scg2.compute_area(), scg1.compute_worst_delay(), scg2.compute_worst_delay(), cec1, cec2 );

//    rewrub_sc_params rps;
//    rewrub_sc_stats rst_p1;
//
//    double aold1 = scg.compute_area();
//    double dold1 = scg.compute_worst_delay();
//
//
//  //  rps.max_trials = 10;
//  //  rps.max_pis = 8;
//  //  rps.max_divisors = 32;256u
//    double aold_N = scg.compute_area()+1;
//
//    double time_now=0;
//    std::clock_t begin = clock();
//
//    while( time_now<180 && aold_N > scg.compute_area() )
//    {
//      aold_N = scg.compute_area();
//      aold1 = scg.compute_area();
//      rewrub_sc( scg, database, rps, &rst_p1 );
//      
//      printf("1[a]%6f ", aold );
//      printf("-> %6f ", scg.compute_area() );
//      printf("1[d]%6f ", dold );
//      printf("-> %6f ", scg.compute_worst_delay() );
//
//      std::clock_t end_now = clock();
//      time_now = double(end_now - begin) / CLOCKS_PER_SEC;
//
//      std::cout << std::endl;
//    }
//    std::clock_t end1 = clock();
//
//    aig_network dump = scg.unmap();
//    aaig_old = dump.num_gates()+1;
//    aaig_new = dump.num_gates();
//    dump = abc_opto( dump, benchmark, "balance" );
//    //dump = cleanup_dangling( dump );
//    while( aaig_new < aaig_old )
//    {
//      dump = abc_opto( dump, benchmark, "resyn2rs" );
//      dump = cleanup_dangling( dump );
//      dump = abc_opto( dump, benchmark, "compress2rs" );
//      dump = cleanup_dangling( dump );
//////
//      aaig_old = aaig_new;
//      aaig_new = dump.num_gates();
//      printf("2 %d\n", aaig_new);
//////
//    }
//    auto const cec3 = benchmark == "hyp" ? true : abc_cec( dump, benchmark );
//    if( !cec3 )
//      printf("--------------------------------------------XXXX ERRORRRRRR\n");
//
//
//    scopt::scg_network scg2 = emap2_klut( dump, tech_lib, ps, &st );
//    scg2 = cleanup_scg( scg2 );
//
//    aold_N = scg2.compute_area()+1;
//
//    std::clock_t begin2 = clock();
//    while( time_now<180 && aold_N > scg2.compute_area() )
//    {
//      aold_N = scg2.compute_area();
//      rewrub_sc( scg2, database, rps, &rst_p1 );
//      
//      printf("2[a]%6f ", aold );
//      printf("-> %6f ", scg2.compute_area() );
//      printf("2[d]%6f ", dold );
//      printf("-> %6f ", scg2.compute_worst_delay() );
//
//      std::clock_t end_now = clock();
//      time_now = double(end_now - begin2) / CLOCKS_PER_SEC;
//
//      std::cout << std::endl;
//    }
//    std::clock_t end2 = clock();
//
//
//    printf("2)-> %6f ", scg2.compute_area() );
//    printf("{d}%6f ", dold );
//    printf("-> %6f ", scg2.compute_worst_delay() );
//
//    double time2 = double(end2 - begin) / CLOCKS_PER_SEC;
//    double time1 = double(end1 - begin) / CLOCKS_PER_SEC;
//
//    start=false;       
//
//    double const d_opt = scg.compute_worst_delay();
//


//    if(!cec2)
//      printf("ERROR2\n");
//
//
//    rarea1 = rarea1*(N-1)/N+(scg.compute_area()-aold)/(N*aold);
//    rdept1 = rdept1*(N-1)/N+(scg.compute_worst_delay()-dold)/(N*dold);
//
//    rarea2 = rarea2*(N-1)/N+(scg2.compute_area()-aold)/(N*aold);
//    rdept2 = rdept2*(N-1)/N+(scg2.compute_worst_delay()-dold)/(N*dold);
//
//    double dA1 = (scg.compute_area()-aold)/(aold);
//    double dA2 = (scg2.compute_area()-aold)/(aold);
//
//    double dD1 = (scg.compute_worst_delay()-dold)/(dold);
//    double dD2 = (scg2.compute_worst_delay()-dold)/(dold);
//
//    printf(" n1 =%f  aN =%f\n", dA1, dA2 );
//    printf("<a1>=%f <d1>=%f\n", rarea1, rdept1 );
//    printf("<aN>=%f <dN>=%f\n", rarea2, rdept2 );
//    std::cout << std::endl;
//    N+=1;
//
//
  }

  exp.save();
  exp.table();
  return 0;
}
