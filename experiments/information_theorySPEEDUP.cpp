#include <iostream>
//#include <catch.hpp>
//#include <execution>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <bitset>
#include <cmath>
#include <random>
#include <boost/dynamic_bitset.hpp>

#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/aig.hpp>

#include <mockturtle/networks/cover.hpp>
#include <mockturtle/algorithms/cover_to_graph.hpp>

#include <mockturtle/networks/pla2.hpp>// plaT for bottom, plaT0 for only greedy
//#include <mockturtle/networks/plaT0.hpp>
#include <lorina/blif.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/views/names_view.hpp>
#include <typeinfo>
#include <lorina/lorina.hpp>

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <fstream>
#include <string>
#include <omp.h>
#include <unistd.h>


//  #include <linux/perf_event.h>    /* Definition of PERF_* constants */
//  #include <linux/hw_breakpoint.h> /* Definition of HW_* constants */
//  #include <sys/syscall.h>         /* Definition of SYS_* constants */

//  int syscall(SYS_perf_event_open, struct perf_event_attr *attr,
//              pid_t pid, int cpu, int group_fd, unsigned long flags);

using namespace mockturtle;
using namespace kitty;

struct splitted_line{
  std::string first;
  std::string second;
};

splitted_line split_string_by_space( std::string line )
{
  splitted_line v_line;
  std::string delimiter = " ";
  size_t pos = 0; 
  std::string token;
  while ((pos = line.find(delimiter)) != std::string::npos) {
          token = line.substr(0, pos);
          v_line.first =  token ;
          line.erase(0, pos + delimiter.length());
          v_line.second = line;
        }
  return v_line;
}

struct XYdataset{
  std::vector<boost::dynamic_bitset<>> X;
  std::vector<boost::dynamic_bitset<>> Y;
  uint64_t nin;
  uint64_t nout;
  uint64_t ndata;
  
};

XYdataset dataset_loader( std::string file_name )
{
  using dyn_bitset = boost::dynamic_bitset<>;
  //using dbs_storage = std::vector<dyn_bitset>;
  XYdataset DS;

  std::string line;
  std::ifstream myfile ( file_name );
  if (myfile.is_open())
  {
    uint32_t r = 0;
    while ( getline (myfile,line) )
    {
      auto v_line = split_string_by_space( line );

      if ( line[0] == '.' )
      {
        if( v_line.first == ".i" )
        {
          DS.nin = std::stoi(v_line.second);
          //std::cout << "nin " << DS.nin << std::endl;
        }
        else if( v_line.first == ".o" )
        {
          DS.nout = std::stoi(v_line.second);
          //std::cout << "nout " << DS.nout << std::endl;

        }
        if( v_line.first == ".p" )
        {
          DS.ndata = std::stoi(v_line.second);
          //std::cout << "ndata " << DS.ndata << std::endl;

          boost::dynamic_bitset<>  empty_bitset( DS.ndata, 0u );
          for( uint32_t i{0u}; i < DS.nin; ++i )
            DS.X.push_back( empty_bitset );
          
          for( uint32_t i{0u}; i < 1; ++i )
            DS.Y.push_back( empty_bitset );
          
        }
      }
      else
      {
        dyn_bitset xline( v_line.first );
        dyn_bitset yline(v_line.second);

        for( uint32_t i{0u}; i < DS.nin; ++i )
            DS.X[i][r] = xline[i];

        for( uint32_t i{0u}; i < 1; ++i )
            DS.Y[i][r] = yline[i];
        r++;
      }
    }
    myfile.close();
  }
  else std::cout << "Unable to open file";
  return DS;
}



int main()
{
std::vector<uint32_t> BKS = { 74,43,21,30,20,40,50,69, 73, 75, 76, 77, 78,79
                              };//40:0
//perf_event_open();

std::cout << "NUM THREADS = " << omp_get_max_threads() << std::endl;
omp_set_num_threads( 8 );
//int num_threads = omp_get_num_threads();
#pragma omp parallel for 
for (uint32_t i = 0 ; i<100; i++) {
  uint32_t bsk = i;
  std::string str_code;
  if( bsk < 10 )
    str_code = "0"+std::to_string(bsk);
  else
    str_code = std::to_string(bsk);

  std::string path_TO_file = "/home/acostama/PhD/ADAPTIVE/IDSD/" + str_code + ".txt";
  std::string path_train = "/home/acostama/PhD/mockturtle/benchmarks/iwls2020-lsml-contest/benchmarks/train/train_txt/ex"+str_code+".train.txt";
  std::string path_test = "/home/acostama/PhD/mockturtle/benchmarks/iwls2020-lsml-contest/benchmarks/test/test_txt/ex"+str_code+".test.txt";
  std::string path_valid = "/home/acostama/PhD/mockturtle/benchmarks/iwls2020-lsml-contest/benchmarks/validation/validation_txt/ex"+str_code+".valid.txt";
    
  //using dyn_bitset = boost::dynamic_bitset<>;
  //using dbs_storage = std::vector<dyn_bitset>;

  std::cout << "condescending" << std::endl;
  auto DCl = dataset_loader( path_train );
  auto DCt = dataset_loader( path_test );
  auto DCv = dataset_loader( path_valid );
  pla2_network cpla( DCl.X, DCl.Y );
  cpla.add_output_file( path_TO_file, str_code );
  cpla.add_output_file( path_TO_file, str_code );

  std::cout << "informed" << std::endl;
  

  auto Dl = dataset_loader( path_train );
  auto Dt = dataset_loader( path_test );
  auto Dv = dataset_loader( path_valid );

  /*for( size_t i {0}; i<Dl.X.size();++i )
  {
    for( size_t j {0}; j<(int)(1*(float)Dv.X[0].size());++j )
      Dl.X[i].push_back(Dv.X[i][j]);
  }
  for( size_t j {0}; j<(int)(1*(float)Dv.X[0].size());++j )
    Dl.Y[0].push_back(Dv.Y[0][j]);*/
  
  pla2_network ipla( Dl.X, Dl.Y );
  bool top_decompose = true;
  bool bottom_decompose=true;
  bool dontknows = true;
  bool informed=true; 
  ipla.set_preferences( top_decompose , bottom_decompose,
                          dontknows, informed );
  ipla.add_output_file( path_TO_file, str_code );
  ipla.ME( Dl.X, Dl.Y, Dt.X, Dt.Y, Dv.X, Dv.Y );

}
  
  boost::dynamic_bitset<> i1(3);
  i1[0] = 0;
  i1[1] = 0;
  i1[2] = 1;
  boost::dynamic_bitset<> i2(3);
  i2[0] = 0;
  i2[1] = 1;
  i2[2] = 0;
  boost::dynamic_bitset<> i3(3);
  i3[0] = 0;
  i3[1] = 1;
  i3[2] = 1;
  pla2_network pla( {i1, i2}, {i3} );

  std::cout << pla.MI({i1},{i3})<<std::endl;
  std::cout << pla.MI({i2},{i3})<<std::endl;
  std::cout << pla.MI({i1,i2},{i3})<<std::endl;
  return 0;
}

  /*pla.print_pla();
  auto Vpr = pla.Pr( {Dl.X[0],Dl.X[1]});
  pla.print_Pr(Vpr);
  std::cout << pla.MI({Dl.X[0]},{Dl.Y[0]}) << std::endl;
  std::cout << pla.MI({Dl.X[1]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[2]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[3]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[4]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[0], Dl.X[1]},{Dl.Y[0]}) << std::endl;
  std::cout << pla.MI({Dl.X[0], Dl.X[2]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[0], Dl.X[3]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[0], Dl.X[4]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[1], Dl.X[2]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[1], Dl.X[3]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[1], Dl.X[4]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[2], Dl.X[3]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[2], Dl.X[4]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[3], Dl.X[4]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[0], Dl.X[1], Dl.X[2]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[0], Dl.X[1], Dl.X[3]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[0], Dl.X[1], Dl.X[4]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[0], Dl.X[2], Dl.X[3]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[0], Dl.X[2], Dl.X[4]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[0], Dl.X[3], Dl.X[4]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[1], Dl.X[2], Dl.X[3]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[1], Dl.X[2], Dl.X[4]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[2], Dl.X[3], Dl.X[4]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[0], Dl.X[1], Dl.X[2], Dl.X[3]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[0], Dl.X[1], Dl.X[2], Dl.X[4]},{Dl.Y[0]})<< std::endl;
  std::cout << pla.MI({Dl.X[0], Dl.X[1], Dl.X[2], Dl.X[3], Dl.X[4]},{Dl.Y[0]})<< std::endl;


  std::cout << "a" << std::endl;
  std::cout << "F41"<< std::endl;
  auto pair = pla.compute_cofactor( Dl.X, Dl.Y, 4, 1 );
  std::cout << "b"<< std::endl;
  pla.print_pla( std::make_pair(Dl.X, Dl.Y) );
  std::cout << "c"<< std::endl;
  pla.print_pla( pair );
  std::cout << "F40"<< std::endl;
  pair = pla.compute_cofactor( Dl.X, Dl.Y, 4, 0 );
  pla.print_pla( pair );

  std::cout << "d"<< std::endl;
  std::cout << "create function" << std::endl;
  std::cout << "TT= " << pla.create_function( Dl.X, Dl.Y, {0,1} ) << std::endl;
  pla.print_pla( std::make_pair(Dl.X, Dl.Y) );
  std::cout << "num nodes " << pla._num_nodes << std::endl;*/

  /*
  std::cout << "3,4,4,4: " <<cpla.Pk_f( 3,4,4,4 ) << std::endl;
  auto R = cpla.M1M2k( 4,4,4 );
  auto K = cpla.num_intersections( 4,4,4 );
  std::cout << "M1=" <<R.first<< " std=" <<R.second<< " num int="<< K << std::endl;
  std::cout << "4,4,4,4: " <<cpla.Pk_f( 4,4,4,4 ) << std::endl;
  R = cpla.M1M2k( 4,4,4 );
  K=cpla.num_intersections( 4,4,4 );
  std::cout << "M1=" <<R.first<< " std=" <<R.second<< " num int="<< K << std::endl;
  std::cout << "5,4,4,4: " <<cpla.Pk_f( 5,4,4,4 ) << std::endl;
  R = cpla.M1M2k( 4,4,4 );
  K=cpla.num_intersections( 4,4,4 );
  std::cout << "M1=" <<R.first<< " std=" <<R.second<< " num int="<< K << std::endl;
  std::cout << "2,0,4,4: " <<cpla.Pk_f( 2,0,4,4 ) << std::endl;
  R = cpla.M1M2k( 0,4,4 );
  K=cpla.num_intersections( 0,4,4 );
  std::cout << "M1=" <<R.first<< " std=" <<R.second<< " num int="<< K << std::endl;
  std::cout << "0,3,3,3: " <<cpla.Pk_f( 0,3,3,3 ) << std::endl;
  R = cpla.M1M2k( 3,3,3 );
  K=cpla.num_intersections( 3,3,3 );
  std::cout << "M1=" <<R.first<< " std=" <<R.second<< " num int="<< K << std::endl;
  std::cout << "3,5,8,10: " <<cpla.Pk_f( 3,5,8,10 ) << std::endl;
  R = cpla.M1M2k( 5,8,10 );
  K=cpla.num_intersections( 5,8,10 );
  std::cout << "M1=" <<R.first<< " std=" <<R.second<< " num int="<< K << std::endl;
  std::cout << "3,5,8,16: " <<cpla.Pk_f( 3,5,8,16 ) << std::endl;
  R = cpla.M1M2k( 5,8,16 );
  K=cpla.num_intersections( 5,8,16 );
  std::cout << "M1=" <<R.first<< " std=" <<R.second<< " num int="<< K << std::endl;
  std::cout << "0,5,8,16: " <<cpla.Pk_f( 0,5,8,16 ) << std::endl;
  std::cout << "3,5,8,20: " <<cpla.Pk_f( 3,5,8,20 ) << std::endl;
  R = cpla.M1M2k( 5,8,20 );
  K=cpla.num_intersections( 5,8,20 );
  std::cout << "M1=" <<R.first<< " std=" <<R.second<< " num int="<< K << std::endl;
  std::cout << "3,5,8,32: " <<cpla.Pk_f( 3,5,8,32 ) << std::endl;
  R = cpla.M1M2k( 5,8,32 );
  K=cpla.num_intersections( 5,8,32 );
  std::cout << "M1=" <<R.first<< " std=" <<R.second<< " num int="<< K << std::endl;
  std::cout << "3,5,8,35: " <<cpla.Pk_f( 3,5,8,35 ) << std::endl;
  R = cpla.M1M2k( 5,8,35 );
  K=cpla.num_intersections( 5,8,35 );
  std::cout << "M1=" <<R.first<< " std=" <<R.second<< " num int="<< K << std::endl;
  std::cout << "3,5,8,40: " <<cpla.Pk_f( 3,5,8,40 ) << std::endl;
  R = cpla.M1M2k( 5,8,40 );
  K=cpla.num_intersections( 5,8,40 );
  std::cout << "M1=" <<R.first<< " std=" <<R.second<< " num int="<< K << std::endl;
  std::cout << "3,5,8,100: " <<cpla.Pk_f( 3,5,8,100 ) << std::endl;
  R = cpla.M1M2k( 5,8,100 );
  K=cpla.num_intersections( 5,8,100 );
  std::cout << "M1=" <<R.first<< " std=" <<R.second<< " num int="<< K << std::endl;
  //cpla.cdsd( DCl.X, DCl.Y, DCt.X, DCt.Y, DCv.X, DCv.Y );
*/