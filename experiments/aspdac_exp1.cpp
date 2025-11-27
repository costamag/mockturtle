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
#include <mockturtle/networks/mapped.hpp>
#include <mockturtle/utils/mapped_database.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/rewrub_mapped.hpp>
#include <ctime>
#include <experiments.hpp>

using namespace experiments;
using namespace mockturtle;
using namespace std::chrono;

template<class Ntk>
Ntk abc_strash( Ntk const& ntk, std::string str_code )
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

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, double,double,double,double,double,double, double,double,double,double,double,double, double,double,double,double,  bool,bool,bool,bool> exp( "DELAY_AWARE", "benchmark", "A(i)", "dA(i,1)", "dA(i,T)", "A(e)", "dA(e,1)", "dA(e,T)", "D(i)", "dD(i,1)", "dD(i,T)", "D(e)", "dD(e,1)", "dD(e,T)", "t(i,1)", "t(i,T)", "t(e,1)", "t(e,T)", "e(i,1)", "e(i,T)","e(e,1)","e(e,T)");

  fmt::print( "[i] processing technology library\n" );

  double RUNTIME_LIMIT=180;

  /* library to map to technology */
  std::vector<gate> gates;
  std::ifstream in( cell_libraries_path( "asap7" ) );

  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  tech_library_params tps;
  tech_library<9, classification_type::np_configurations> tech_lib( gates, tps );

  mapped_database_params ps;
  ps.num_vars = 4;
  ps.mps.area_oriented_mapping = true;
  mapped_database database( "../experiments/databases/mapped/asap7_area.db", tech_lib, ps );

  double dA_i_1{0};
  double dA_i_T{0};
  double dA_e_1{0};
  double dA_e_T{0};

  double dD_i_1{0};
  double dD_i_T{0};
  double dD_e_1{0};
  double dD_e_T{0};

  double DT_i_1{0};
  double DT_i_T{0};
  double DT_e_1{0};
  double DT_e_T{0};


  double raA_i_1{0};
  double raA_i_T{0};
  double raA_e_1{0};
  double raA_e_T{0};
  double raD_i_1{0};
  double raD_i_T{0};
  double raD_e_1{0};
  double raD_e_T{0};
  double raT_i_1{0};
  double raT_i_T{0};
  double raT_e_1{0};
  double raT_e_T{0};
  double nIteras{1};

  double area__Map_i{0};
  double delay_Map_i{0};
  double area__Map_e{0};
  double delay_Map_e{0};

  emap_params mps;
  mps.area_oriented_mapping = true;
  emap_stats st;

  rewrub_mapped_params rps;

  rps.delay_awareness = false;
  //rps.preserve_depth = true;
  rps.max_divisors = 20;
  rps.try_struct = true;
  rps.try_window = true;
  rps.try_simula = true;

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig_i;
    if ( lorina::read_aiger( experiments::benchmark_path(benchmark), aiger_reader( aig_i ) ) != lorina::return_code::success )
    {
      continue;
    }
    /* read the initial AIG size */
    if( aig_i.num_gates() > 200000 || benchmark == "hyp" )
      continue;

    binding_view<klut_network> luts = emap_klut( aig_i, tech_lib, mps, &st );

    printf("read as network type\n" );
    mapped_network scg_i( gates, luts );
    printf("mapped %f\n", scg_i.compute_area() );
    /* read the initial MAP size and delay*/
    area__Map_i = scg_i.compute_area();
    delay_Map_i = scg_i.compute_worst_delay();
    rps.required_time = delay_Map_i;

    /* optimize once */
    std::clock_t start_i_1 = clock();
    rewrub_mapped( scg_i, database, rps );
    std::clock_t time_i_1 = clock();
    DT_i_1 = double(time_i_1 - start_i_1) / CLOCKS_PER_SEC;
    dA_i_1 = (scg_i.compute_area()-area__Map_i)/(area__Map_i)*100;
    dD_i_1 = (scg_i.compute_worst_delay()-delay_Map_i)/(delay_Map_i)*100;

    printf("A[i, 1]= %6.2f -> dA[i,1] = %6.2f D[i, 1] = %6.2f -> dD[i,1] = %6.2f \n", scg_i.compute_area(), dA_i_1, scg_i.compute_worst_delay(), dD_i_1 );

    const auto cec_i_1 = benchmark == "hyp" ? true : abc_cec( scg_i, benchmark );
    if(!cec_i_1)
    {
      printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  ERROR i 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"); 
    }

    /* iteratively optimize */
    printf("PHASE [i,T]\n");
    std::clock_t start_i_T = clock();
    DT_i_T=DT_i_1;
    double aMap_i_T = scg_i.compute_area()+1;
    while( DT_i_T<RUNTIME_LIMIT && aMap_i_T > scg_i.compute_area() )
    {
      aMap_i_T = scg_i.compute_area();
      rewrub_mapped( scg_i, database, rps );
      
      printf("a[i, t]= %6.2f\n", scg_i.compute_area() );

      std::clock_t time_i_t = clock();
      DT_i_T = DT_i_1 + double(time_i_t - start_i_T) / CLOCKS_PER_SEC;
    }
    dA_i_T = (scg_i.compute_area()-area__Map_i)/(area__Map_i)*100;
    dD_i_T = (scg_i.compute_worst_delay()-delay_Map_i)/(delay_Map_i)*100;
    printf("A[i, T]= %6.2f -> dA[i,T] = %6.2f D[i, T] = %6.2f -> dD[i,T] = %6.2f \n", scg_i.compute_area(), dA_i_T, scg_i.compute_worst_delay(), dD_i_T );

    const auto cec_i_T = benchmark == "hyp" ? true : abc_cec( scg_i, benchmark );
    if(!cec_i_T)
    {
      printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  ERROR i T !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"); 
    }

    /* iteratively optimize the AIG */
    aig_network aig_e;
    if ( lorina::read_aiger( experiments::benchmark_path(benchmark), aiger_reader( aig_e ) ) != lorina::return_code::success )
    {
      continue;
    }
    /* read the initial AIG size */
    auto area_Aig = aig_e.num_gates() + 1;
    while( area_Aig > aig_e.num_gates() )
    {
      area_Aig = aig_e.num_gates();

      aig_e = abc_opto( aig_e, benchmark, "resyn2rs" );
      aig_e = cleanup_dangling( aig_e );
      printf("resyn2rs    -> %d\n", aig_e.num_gates() );
      aig_e = abc_opto( aig_e, benchmark, "compress2rs" );
      aig_e = cleanup_dangling( aig_e );
      printf("compress2rs -> %d\n", aig_e.num_gates() );
    }
    binding_view<klut_network> luts2 = emap_klut( aig_e, tech_lib, mps, &st );

    mapped_network scg_e( gates, luts2 );
    scg_e = cleanup_mapped( scg_e );
    /* read the initial MAP size and delay*/
    area__Map_e = scg_e.compute_area();
    delay_Map_e = scg_e.compute_worst_delay();
    rps.required_time = delay_Map_e;

    /* optimize once */
    std::clock_t start_e_1 = clock();
    rewrub_mapped( scg_e, database, rps );
    std::clock_t time_e_1 = clock();
    DT_e_1 = double(time_e_1 - start_e_1) / CLOCKS_PER_SEC;
    dA_e_1 = (scg_e.compute_area()-area__Map_e)/(area__Map_e)*100;
    dD_e_1 = (scg_e.compute_worst_delay()-delay_Map_e)/(delay_Map_e)*100;

    printf("A[e, 1]= %6.2f -> dA[e,1] = %6.2f D[e, 1] = %6.2f -> dD[e,1] = %6.2f \n", scg_e.compute_area(), dA_e_1, scg_e.compute_worst_delay(), dD_e_1 );

    const auto cec_e_1 = benchmark == "hyp" ? true : abc_cec( scg_e, benchmark );
    if(!cec_e_1)
    {
      printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  ERROR e 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"); 
    }

    /* iteratively optimize */
    printf("PHASE [e,T]\n");
    std::clock_t start_e_T = clock();

    double aMap_e_T = scg_e.compute_area()+1;

    DT_e_T = DT_e_1;
    while( DT_e_T<RUNTIME_LIMIT && aMap_e_T > scg_e.compute_area() )
    {
      aMap_e_T = scg_e.compute_area();
      rewrub_mapped( scg_e, database, rps );
      printf("a[e, t]= %6.2f\n", scg_e.compute_area() );
      std::clock_t time_e_t = clock();
      DT_e_T = DT_e_1 + double(time_e_t - start_e_T) / CLOCKS_PER_SEC;
    }

    dA_e_T = (scg_e.compute_area()-area__Map_e)/(area__Map_e)*100;
    dD_e_T = (scg_e.compute_worst_delay()-delay_Map_e)/(delay_Map_e)*100;
    printf("A[e, T]= %6.2f -> dA[e,T] = %6.2f D[e, T] = %6.2f -> dD[e,T] = %6.2f \n", scg_e.compute_area(), dA_e_T, scg_e.compute_worst_delay(), dD_e_T );

    const auto cec_e_T = benchmark == "hyp" ? true : abc_cec( scg_e, benchmark );
    if(!cec_e_T)
    {
      printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  ERROR e T !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"); 
    }


    printf("dA[i,1]=%6.2f\n", dA_i_1 );
    printf("dA[i,T]=%6.2f\n", dA_i_T );
    printf("dA[e,1]=%6.2f\n", dA_e_1 );
    printf("dA[e,T]=%6.2f\n\n", dA_e_T );

    printf("dD[i,1]=%6.2f\n", dD_i_1 );
    printf("dD[i,T]=%6.2f\n", dD_i_T );
    printf("dD[e,1]=%6.2f\n", dD_e_1 );
    printf("dD[e,T]=%6.2f\n\n", dD_e_T );

    printf("DT[i,1]=%6.2f\n", DT_i_1 );
    printf("DT[i,T]=%6.2f\n", DT_i_T );
    printf("DT[e,1]=%6.2f\n", DT_e_1 );
    printf("DT[e,T]=%6.2f\n\n", DT_e_T );

    if( cec_i_1 && cec_i_T && cec_e_1 && cec_e_T )
    {
      raA_i_1 = raA_i_1*(nIteras-1)/(nIteras) + dA_i_1/(nIteras);
      raA_i_T = raA_i_T*(nIteras-1)/(nIteras) + dA_i_T/(nIteras);
      raA_e_1 = raA_e_1*(nIteras-1)/(nIteras) + dA_e_1/(nIteras);
      raA_e_T = raA_e_T*(nIteras-1)/(nIteras) + dA_e_T/(nIteras);
      raD_i_1 = raD_i_1*(nIteras-1)/(nIteras) + dD_i_1/(nIteras);
      raD_i_T = raD_i_T*(nIteras-1)/(nIteras) + dD_i_T/(nIteras);
      raD_e_1 = raD_e_1*(nIteras-1)/(nIteras) + dD_e_1/(nIteras);
      raD_e_T = raD_e_T*(nIteras-1)/(nIteras) + dD_e_T/(nIteras);
      raT_i_1 = raT_i_1*(nIteras-1)/(nIteras) + DT_i_1/(nIteras);
      raT_i_T = raT_i_T*(nIteras-1)/(nIteras) + DT_i_T/(nIteras);
      raT_e_1 = raT_e_1*(nIteras-1)/(nIteras) + DT_e_1/(nIteras);
      raT_e_T = raT_e_T*(nIteras-1)/(nIteras) + DT_e_T/(nIteras);
      nIteras+=1;
    }

    printf( "<A_i_1> = %6.2f\n", raA_i_1);
    printf( "<A_i_T> = %6.2f\n", raA_i_T);
    printf( "<A_e_1> = %6.2f\n", raA_e_1);
    printf( "<A_e_T> = %6.2f\n", raA_e_T);
    printf( "<D_i_1> = %6.2f\n", raD_i_1);
    printf( "<D_i_T> = %6.2f\n", raD_i_T);
    printf( "<D_e_1> = %6.2f\n", raD_e_1);
    printf( "<D_e_T> = %6.2f\n", raD_e_T);
    printf( "<T_i_1> = %6.2f\n", raT_i_1);
    printf( "<T_i_T> = %6.2f\n", raT_i_T);
    printf( "<T_e_1> = %6.2f\n", raT_e_1);
    printf( "<T_e_T> = %6.2f\n", raT_e_T);

    exp( benchmark, area__Map_i, dA_i_1, dA_i_T, area__Map_e, dA_e_1, dA_e_T, delay_Map_i, dD_i_1, dD_i_T, delay_Map_e, dD_e_1, dD_e_T, DT_i_1, DT_i_T, DT_e_1, DT_e_T, cec_i_1, cec_i_T, cec_e_1, cec_e_T );

  }

  exp.save();
  exp.table();

  return 0;

}
