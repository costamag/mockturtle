#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <bitset>
#include <cmath>
#include <random>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/aig.hpp>

#include <mockturtle/networks/cover.hpp>
#include <mockturtle/algorithms/cover_to_graph.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>

#include <mockturtle/algorithms/mi_decomposition_ptt.hpp>

#include <lorina/blif.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/views/names_view.hpp>
#include <typeinfo>
#include <lorina/lorina.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <fstream>
#include <string>
#include <omp.h>
#include <unistd.h>

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

XYdataset dataset_loader( std::string file_name )
{
  using dyn_bitset = kitty::partial_truth_table;
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
        }
        else if( v_line.first == ".o" )
        {
          DS.nout = std::stoi(v_line.second);
        }
        if( v_line.first == ".p" )
        {
          DS.ndata = std::stoi(v_line.second);
          kitty::partial_truth_table empty_bitset( DS.ndata );
          for( uint32_t i{0u}; i < DS.nin; ++i )
            DS.X.push_back( empty_bitset );
                    
        }
      }
      else
      {
        dyn_bitset xline( DS.nin );
        kitty::create_from_binary_string( xline, v_line.first );
        dyn_bitset yline( 1 );
        kitty::create_from_binary_string( yline, v_line.second );

        for( uint32_t i{0u}; i < DS.nin; ++i )
          kitty::get_bit(xline,i) == 1 ? kitty::set_bit(DS.X[i],r) : kitty::clear_bit(DS.X[i],r) ;

        
        DS.Y.add_bit( kitty::get_bit( yline,0 ) );
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
std::vector<uint32_t> BKS = { 94,95,96,97,98,99 };//40:0

std::cout << "*** simulations : iwls2020 ***" << std::endl;

std::cout << "NUM THREADS = " << omp_get_max_threads() << std::endl;
omp_set_num_threads( 1 );
//int num_threads = omp_get_num_threads();
#pragma omp parallel for 
for (uint32_t i = 0 ; i<100; i++) {
  uint32_t bsk = i;
  std::string str_code;
  if( bsk < 10 )
    str_code = "0"+std::to_string(bsk);
  else
    str_code = std::to_string(bsk);

  std::string path_TO_file = "/home/acostama/projects/EPFL/mockturtle/simulations/iwls20/" + str_code + ".txt";
  std::string path_train = "/home/acostama/projects/EPFL/mockturtle/benchmarks/iwls2020-lsml-contest/benchmarks/train/train_txt/ex"+str_code+".train.txt";
  std::string path_test = "/home/acostama/projects/EPFL/mockturtle/benchmarks/iwls2020-lsml-contest/benchmarks/test/test_txt/ex"+str_code+".test.txt";
  std::string path_valid = "/home/acostama/projects/EPFL/mockturtle/benchmarks/iwls2020-lsml-contest/benchmarks/validation/validation_txt/ex"+str_code+".valid.txt";
      

  auto Dl = dataset_loader( path_train );
  auto Dt = dataset_loader( path_test );
  auto Dv = dataset_loader( path_valid );

  mi_decomposition_params ps;
  ps.max_sup = 4;
  ps.is_informed = true;
  ps.try_top_decomposition = true;
  ps.try_creation = false;
  ps.try_xor_decomposition = true;
  ps.use_cumsum = false;
  ps.try_bottom_decomposition = true;
  ps.is_bottom_exact = true;

  ps.dontcares = true;

  klut_network klut;
  auto res = mi_decomposition_iwls20( Dl, klut, ps );

  auto aig = convert_klut_to_graph<aig_network>( klut );
  depth_view_params psD;
  psD.count_complements = false;
  depth_view depth_aig{aig, {}, psD};
  std::cout << ".bk " << bsk <<std::endl;
  std::cout << ".la " << compute_accuracy( Dl.X, Dl.Y, depth_aig )<<std::endl;
  std::cout << ".ta " << compute_accuracy( Dt.X, Dt.Y, depth_aig )<<std::endl;
  std::cout << ".va " << compute_accuracy( Dv.X, Dv.Y, depth_aig )<<std::endl;
  std::cout << ".ng " << depth_aig.num_gates() << std::endl;
  std::cout << ".sz " << depth_aig.size() << std::endl;
  std::cout << ".dt " << depth_aig.depth() << std::endl;
  std::cout << ".1t " << res._cnt._or << std::endl;
  std::cout << ".0t " << res._cnt._le << std::endl;
  std::cout << ".1c " << res._cnt._lt << std::endl;
  std::cout << ".0c " << res._cnt._and << std::endl;
  std::cout << ".ch " << res._cnt._ctj << std::endl;
  std::cout << ".bd " << res._cnt._btm << std::endl;
  std::cout << std::endl;




  std::ofstream myfile;
  myfile.open ( ("/home/acostama/projects/EPFL/mockturtle/simulations/iwls20/creation/"+str_code+".txt") ); 
  myfile << ".bk " << bsk <<std::endl;
  myfile << ".la " << compute_accuracy( Dl.X, Dl.Y, depth_aig )<<std::endl;
  myfile << ".ta " << compute_accuracy( Dt.X, Dt.Y, depth_aig )<<std::endl;
  myfile << ".va " << compute_accuracy( Dv.X, Dv.Y, depth_aig )<<std::endl;
  myfile << ".ng " << depth_aig.num_gates() << std::endl;
  myfile << ".sz " << depth_aig.size() << std::endl;
  myfile << ".dt " << depth_aig.depth() << std::endl;
  myfile << ".1t " << res._cnt._or << std::endl;
  myfile << ".0t " << res._cnt._le << std::endl;
  myfile << ".1c " << res._cnt._lt << std::endl;
  myfile << ".0c " << res._cnt._and << std::endl;
  myfile << ".ch " << res._cnt._ctj << std::endl;
  myfile << ".bd " << res._cnt._btm << std::endl;
  myfile.close();
  std::cout << std::endl;
}
  return 0;
}
