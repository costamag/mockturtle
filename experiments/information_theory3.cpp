#include <iostream>
//#include <catch.hpp>

#include <sstream>
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

#include <mockturtle/networks/plaT.hpp>
#include <lorina/blif.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/views/names_view.hpp>
#include <typeinfo>
#include <lorina/lorina.hpp>

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <fstream>
#include <string>

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
  uint32_t nin;
  uint32_t nout;
  uint32_t ndata;
  
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
        }
      }
      else
      {
        dyn_bitset xtrain(v_line.first);
        xtrain.push_back(0);
        DS.X.push_back(xtrain);        
        dyn_bitset ytrain(v_line.second);
        DS.Y.push_back(ytrain);

      }
      
    }
    myfile.close();
  }
  else std::cout << "Unable to open file";

  return DS;

}

double computeAcc(std::vector<boost::dynamic_bitset<>> inputs,std::vector<boost::dynamic_bitset<>> outputs, aig_network aig )
{
  //std::cout << aig.size() << std::endl;
  //std::cout << aig.num_gates() << std::endl;
  double acc = 0;
  //std::cout << "SZ=" << inputs[0].size() << std::endl;

  for (uint32_t d{0u}; d < inputs.size(); ++d )
  {
    std::vector<bool> inpt_v = {};
    //
    //for( int k = 0; k< (inputs[0].size()-1);++k )
    for( int k = (inputs[0].size()-2); k>= 0;--k )
    {
      //std::cout << "i:" << inputs[d][k] << "-";
      inpt_v.push_back( inputs[d][k] == 1 );
      //std::cout << inpt_v[k] << " ";
    }
    /*std::cout << "\nps:";
    for( uint32_t k{0u}; k<inpt_v.size();++k )
    {
      std::cout << inpt_v[k];
    }
    std::cout << std::endl;

    std::cout << "pg:" << inputs[d] << " ";*/
    bool sim_res = simulate<bool>( aig, default_simulator<bool>( inpt_v ) )[0];
    //std::cout << sim_res << ":" << outputs[d][0] << " ";
    //std::cout << std::endl;

    double deltaA = ( (sim_res == outputs[d][0]) ? (double)1.0/outputs.size() : 0.0 );
    acc += deltaA;
  }
  return 100*acc;
}


int main()
{  
  for( uint32_t it = 0; it < 10; ++it )
  {
  std::string str_code;
  if( it < 10 )
    str_code = "0"+std::to_string(it);
  else
    str_code = std::to_string(it);

  std::string path_train = "/home/acostama/PhD/mockturtle/benchmarks/iwls2020-lsml-contest/benchmarks/train/train_txt/ex"+str_code+".train.txt";
  std::string path_test = "/home/acostama/PhD/mockturtle/benchmarks/iwls2020-lsml-contest/benchmarks/test/test_txt/ex"+str_code+".test.txt";
  std::string path_valid = "/home/acostama/PhD/mockturtle/benchmarks/iwls2020-lsml-contest/benchmarks/validation/validation_txt/ex"+str_code+".valid.txt";

    
  using dyn_bitset = boost::dynamic_bitset<>;
  using dbs_storage = std::vector<dyn_bitset>;
  auto train_ds = dataset_loader( path_train );
  //std::cout << "nin = " << train_ds.nin << std::endl;
  //std::cout << "nout = " << train_ds.nout << std::endl;
  //std::cout << "ndata = " << train_ds.ndata << std::endl;

  auto test_ds = dataset_loader( path_test );
  auto valid_ds = dataset_loader( path_valid );

  //std::cout << "nin = " << test_ds.nin << std::endl;
  //std::cout << "nout = " << test_ds.nout << std::endl;
  //std::cout << "ndata = " << test_ds.ndata << std::endl;
  std::cout << "* * * * * * * * * * * * * * * *  " << std::endl;
  std::cout << "              " << str_code << "            " << std::endl;
  std::cout << "* * * * * * * * * * * * * * * *  " << std::endl;
  
  dyn_bitset N1 ( 5, 1u );
  dyn_bitset N2 ( 5, 3u );
  dyn_bitset N3 ( 5, 8u );


  std::cout << std::endl;
  std::cout << "INFORMED SHANNON + DSD " << std::endl; 
  plaT_network pla3( train_ds.X, train_ds.Y, 2, 4, 2 );
  //pla3.preprocess_muesli(0.3);

  //std::cout << "HD: " << pla3.HammingDistance(N1,N1)<< pla3.HammingDistance(N1,N2) << " " << pla3.HammingDistance(N1,N3) << std::endl ;
  pla3.it_dsd_shannon_decomposition(false, 0);
  std::cout << "test accuracy: " << pla3.compute_accuracy( test_ds.X, test_ds.Y ) << "%" << std::endl;
  std::cout << "valid accuracy: " << pla3.compute_accuracy( valid_ds.X, valid_ds.Y ) << "%" << std::endl;

  std::cout << std::endl;
  /*std::cout << "NOT INFORMED SHANNON " << std::endl; 
  plaT_network pla1( train_ds.X, train_ds.Y, 2, 4, 2 );
  pla1.it_shannon_decomposition(true, 0);
  std::cout << "test accuracy: " << pla1.compute_accuracy( test_ds.X, test_ds.Y ) << "%" << std::endl;
  std::cout << "valid accuracy: " << pla1.compute_accuracy( valid_ds.X, valid_ds.Y ) << "%" << std::endl;
  std::cout << std::endl;
  std::cout << "INFORMED SHANNON " << std::endl; 
  plaT_network pla2( train_ds.X, train_ds.Y, 2, 4, 2 );
  pla2.it_shannon_decomposition(false, 0);
  std::cout << "test accuracy: " << pla2.compute_accuracy( test_ds.X, test_ds.Y ) << "%" << std::endl;
  std::cout << "valid accuracy: " << pla2.compute_accuracy( valid_ds.X, valid_ds.Y ) << "%" << std::endl;
  
  std::cout << std::endl;*/

  }

  return 0;
}