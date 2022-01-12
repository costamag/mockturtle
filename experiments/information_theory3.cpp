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
  std::string str_code = "80";
  std::string path_train = "/home/acostama/PhD/mockturtle/benchmarks/iwls2020-lsml-contest/benchmarks/train/train_txt/ex"+str_code+".train.txt";
  std::string path_test = "/home/acostama/PhD/mockturtle/benchmarks/iwls2020-lsml-contest/benchmarks/test/test_txt/ex"+str_code+".test.txt";
  
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

  plaT_network pla1( train_ds.X, train_ds.Y, 4, 4 );

  for ( uint32_t k{0u}; k<train_ds.nin; k++ )
  {
    std::cout << k << "[" << pla1.MI({k},{0}) << "] " << std::endl;
  }

  /* Espresso STILL COVER TO ADJUST
  std::string path_train_blif = "../LUTs/ex"+str_code+"train.blif";//"../benchmarks/iwls2020-lsml-contest/benchmarks/train/train_txt/ex"+str_code+"train.blif";
  names_view<cover_network> cover_ntk;
  
  if ( lorina::read_blif( path_train_blif, blif_reader( cover_ntk ) ) != lorina::return_code::success )
  {
    std::cout << "read <testcase>.blif failed!\n";
    return -1;
  }

  aig_network aig_cv;
  convert_cover_to_graph( aig_cv, cover_ntk );
  double train_acc = computeAcc( train_ds.X, train_ds.Y, aig_cv );
  double test_acc = computeAcc( test_ds.X, test_ds.Y, aig_cv );
  std::cout << "TEST acc = " <<  test_acc << "%" << std::endl;
  std::cout << "TRAIN acc = " <<  train_acc << "%" << std::endl;
*/
/*#####################################################################
TEST WIT 3AND
std::vector<boost::dynamic_bitset<>> input_nodes;
  for ( uint32_t i {0u}; i < 8; ++i )
  {
    boost::dynamic_bitset<> idata( 3, i );
    input_nodes.push_back( idata );
    input_nodes[i].push_back(0);
  }


  std::vector<boost::dynamic_bitset<>> output_nodes;
  std::vector Voutput_nodes = {0,0,0,0,0,0,0,1};


  for( uint32_t k {0u}; k < Voutput_nodes.size(); ++k )
  {
    boost::dynamic_bitset<> odata( 1, Voutput_nodes.at(k) );
    output_nodes.push_back( odata );
  }
  str_code = "AA";
    std::string path_train_blif = "../LUTs/ex"+str_code+"train.blif";//"../benchmarks/iwls2020-lsml-contest/benchmarks/train/train_txt/ex"+str_code+"train.blif";
  cover_network cover_ntk;
  
  if ( lorina::read_blif( path_train_blif, blif_reader( cover_ntk ) ) != lorina::return_code::success )
  {
    std::cout << "read <testcase>.blif failed!\n";
    return -1;
  }

  aig_network aig_cv;

  //convert_cover_to_graph( aig_cv, cover_ntk );
  aig_cv = convert_cover_to_graph<aig_network>( cover_ntk );
  std::cout << aig_cv.size() << std::endl;
  std::cout << aig_cv.num_gates() << std::endl;

  double train_acc = computeAcc( input_nodes, output_nodes, aig_cv );
  //double test_acc = computeAcc( test_ds.X, test_ds.Y, aig_cv );
  //std::cout << "TEST acc = " <<  test_acc << "%" << std::endl;
  std::cout << "TRAIN acc = " <<  train_acc << "%" << std::endl;

//########################################################################*/
  std::cout << "not informed shannon " << std::endl; 
  plaT_network pla3( train_ds.X, train_ds.Y, 4, 4 );
  pla3.it_shannon_decomposition(true, 0);
  std::cout << "\n test accuracy: " << pla3.compute_accuracy( test_ds.X, test_ds.Y ) << "%" << std::endl;

  std::cout << "informed shannon " << std::endl; 
  pla1.it_shannon_decomposition(false, 0);
  std::cout << "\n test accuracy: " << pla1.compute_accuracy( test_ds.X, test_ds.Y ) << "%" << std::endl;


  std::cout << "informed shannon + dsd " << std::endl; 
  plaT_network pla2( train_ds.X, train_ds.Y, 4, 4 );
  pla2.it_dsd_shannon_decomposition(false, 0);
  std::cout << "\n test accuracy: " << pla2.compute_accuracy( test_ds.X, test_ds.Y ) << "%" << std::endl;



  /* informed SHANNON */

  /* not informed SHANNON TODO*/

  /* MUESLI */
  //pla1.preprocess_muesli(0.1);
  //pla1.muesli();

  /* MUESLI MODIFIED TODO*/

  



  //plaT_network pla_sh( inodes, onodes, 5, 3 );

/*
  std::vector<boost::dynamic_bitset<>> input_nodes;
  for ( uint32_t i {0u}; i < 16; ++i )
  {
    boost::dynamic_bitset<> idata( 5, i );
    input_nodes.push_back( idata );
  }

  std::vector<boost::dynamic_bitset<>> output_nodes;
  std::vector Voutput_nodes = {0,0,0,1,0,0,0,1,0,0,0,1,1,1,1,1}; /* ab + cde 


  for( uint32_t k {0u}; k < Voutput_nodes.size(); ++k )
  {
    boost::dynamic_bitset<> odata( 1, Voutput_nodes.at(k) );
    output_nodes.push_back( odata );
  }

  plaT_network pla( input_nodes, output_nodes, 5 );
  pla.print_pla();
  pla.muesli(2);
  pla.print_pla();*/


  return 0;
}