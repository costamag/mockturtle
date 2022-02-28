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

#include <mockturtle/networks/plaT.hpp>// plaT for bottom, plaT0 for only greedy
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
  return 100;
}


int main()
{ 

std::vector<uint32_t> BKSrd = { 
                              //10, 11, 15, 16, 
                              //30, 31, 32, 33, 34, 
                              //40, 41, 42, 43, 
                              //44, 45,
                              50, 51, 52, 53, 54, 55, 56, 57, 58,
                              60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
                              70, 71, 72, 73, 74, 75, 76, 77, 78, 79
                              };

std::vector<uint32_t> XOR = { 20, 21, 30, 40, 42, 43};

std::vector<uint32_t> BKS = { //0, 1, 2, 3, 5, 7, 
                              //9,
                              //10,
                              //4.  6.  7.  8.  9. 12. 13. 14. 17. 18. 19. 26. 27. 28. 29. 35. 36. 37.
                              //38. 39. 46. 47. 48. 49. 59. 80. 81. 82. 83. 84. 85. 86. 87. 88. 89. 90.
                              //91. 92. 93. 94. 95. 96. 97. 98. 99 
                              26, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89
                              //4,  6,  7,  8,  9, 12, 13, 14, 17, 18, 19, 26, 27, 28, 29, 35, 36, 37,
                              //38, 39, 46, 47, 48, 49, 59, , 90,
                              //91, 92, 93, 94, 95, 96, 97, 98, 99,
                              //21, 22, 23, 24,
                              //20,
                              //30, 31, 32, 33, 
                              //34, 
                              // 40, 41, 43, 44, 45,
                              //42,
                              //50, 51, 52, 56, 57, 58,
                              //53, 54, 55,  
                              // 60, 61, 62, 63, 64, 67, 68, 69,
                              // 70, 71, 72,73, 74, 75, 76, 
                              //77, 78, 79,                               
                              };
/*
for( uint32_t it = 74; it < 75; ++it )
{
  std::string str_code;
  if( it < 10 )
    str_code = "0"+std::to_string(it);
  else
    str_code = std::to_string(it);
*/

BKSrd = BKS;
std::cout << "NUM THREADS = " << omp_get_max_threads() << std::endl;
//omp_set_num_threads( 6 );
//for( uint32_t it = 0; it < BKSrd.size(); ++it )
//{
//int num_threads = omp_get_num_threads();
//#pragma omp parallel for 
for (uint32_t i = 0 ; i<1; i++) {
  // do something with i
  auto bsk =  95;//BKS[i];
        //do stuff with item
  bool is_dec_naive = false;
  bool try_bottom = false;
  bool is_bottom_greedy = false;
  bool only_shannon = true;
  bool try_top_xor = true;
  bool is_bottom_conservative = true;

  /*std::string str_code;
  if( BKSrd[it] < 10 )
    str_code = "0"+std::to_string(BKSrd[it]);
  else
    str_code = std::to_string(BKSrd[it]);*/

std::string str_code;
  if( bsk < 10 )
    str_code = "0"+std::to_string(bsk);
  else
    str_code = std::to_string(bsk);


  std::string path_TO_file = "/home/acostama/PhD/E3/" + str_code + ".txt";
  
  

  std::string path_train = "/home/acostama/PhD/mockturtle/benchmarks/iwls2020-lsml-contest/benchmarks/train/train_txt/ex"+str_code+".train.txt";
  std::string path_test = "/home/acostama/PhD/mockturtle/benchmarks/iwls2020-lsml-contest/benchmarks/test/test_txt/ex"+str_code+".test.txt";
  std::string path_valid = "/home/acostama/PhD/mockturtle/benchmarks/iwls2020-lsml-contest/benchmarks/validation/validation_txt/ex"+str_code+".valid.txt";
    
  //using dyn_bitset = boost::dynamic_bitset<>;
  //using dbs_storage = std::vector<dyn_bitset>;
  auto train_ds = dataset_loader( path_train );
  auto test_ds = dataset_loader( path_test );
  auto valid_ds = dataset_loader( path_valid );


  plaT_network pla3( train_ds.X, train_ds.Y, 2, 4, 2 ); // plaT_network for bottom, plaT0_network for only greedy
  pla3.add_test_set( test_ds.X, test_ds.Y );
  pla3.add_valid_set( valid_ds.X, valid_ds.Y );
  pla3.add_output_file( path_TO_file, str_code );
  uint64_t delta_supp = 3;
  pla3.it_dsd_shannon_decomposition(is_dec_naive, 0, try_bottom, is_bottom_greedy, only_shannon, try_top_xor, is_bottom_conservative, delta_supp );
  //uint32_t max_num_nodes = train_ds.X[0].size()*1.5;
  //pla3.BUDMA( 0.99, max_num_nodes );

  }


  return 0;
}