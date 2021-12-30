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
#include <mockturtle/networks/pla.hpp>

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


int main()
{
  std::string path_train = "/home/acostama/PhD/mockturtle/benchmarks/iwls2020-lsml-contest/benchmarks/train/train_txt/ex00.train.txt";
  std::string path_test = "/home/acostama/PhD/mockturtle/benchmarks/iwls2020-lsml-contest/benchmarks/test/test_txt/ex00.test.txt";
  
  //std::ifstream myfile ( path_train + "/ex00.train.txt" );
  
  using dyn_bitset = boost::dynamic_bitset<>;
  using dbs_storage = std::vector<dyn_bitset>;
  auto train_ds = dataset_loader( path_train );
  std::cout << "nin = " << train_ds.nin << std::endl;
  std::cout << "nout = " << train_ds.nout << std::endl;
  std::cout << "ndata = " << train_ds.ndata << std::endl;

  auto test_ds = dataset_loader( path_test );
  std::cout << "nin = " << test_ds.nin << std::endl;
  std::cout << "nout = " << test_ds.nout << std::endl;
  std::cout << "ndata = " << test_ds.ndata << std::endl;

  pla_network pla1( train_ds.X, train_ds.Y, 3 );

  for ( uint32_t k{0u}; k<train_ds.nin; k++ )
  {
    std::cout << k << "[" << pla1.MI({k},{0}) << "] " << std::endl;
  }
  //pla_network pla1( train_ds.X, train_ds.Y, 3 );
  //pla1.muesli(2);


  std::vector<boost::dynamic_bitset<>> input_nodes;
  for ( uint32_t i {0u}; i < 16; ++i )
  {
    boost::dynamic_bitset<> idata( 5, i );
    input_nodes.push_back( idata );
  }

  std::vector<boost::dynamic_bitset<>> output_nodes;
  std::vector Voutput_nodes = {0,0,0,1,0,0,0,1,0,0,0,1,1,1,1,1}; /* ab + cde */



  for( uint32_t k {0u}; k < Voutput_nodes.size(); ++k )
  {
    boost::dynamic_bitset<> odata( 1, Voutput_nodes.at(k) );
    output_nodes.push_back( odata );
  }

  pla_network pla( input_nodes, output_nodes, 5 );
  pla.print_pla();
  pla.muesli(2);
  pla.print_pla();


  return 0;
}