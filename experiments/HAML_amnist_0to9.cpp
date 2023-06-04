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
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/io/write_aiger.hpp>

#include <mockturtle/networks/cover.hpp>
#include "../include/mockturtle/algorithms/cover_to_graph.hpp"
#include "../include/mockturtle/algorithms/klut_to_graph.hpp"

//#include <mockturtle/algorithms/lfe/mi_decomposition.hpp>
//#include <mockturtle/algorithms/lfe/mi_decomposition.hpp>
#include "../include/mockturtle/algorithms/lfe/hyperdimensional_computing/model.hpp"
#include "../include/mockturtle/algorithms/lfe/hyperdimensional_computing/methods/accuracy_recovery.hpp"
#include "../include/mockturtle/algorithms/lfe/hyperdimensional_computing/methods/generators.hpp"
#include "../include/mockturtle/algorithms/lfe/hyperdimensional_computing/methods/selectors.hpp"
#include "../include/mockturtle/algorithms/lfe/hyperdimensional_computing/methods/selgenerators.hpp"
#include "../include/mockturtle/algorithms/lfe/projectors_in_hd.hpp"


#include <lorina/blif.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/views/names_view.hpp>
#include <typeinfo>
#include <lorina/lorina.hpp>
#include <lorina/truth.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/partial_truth_table.hpp>
#include <fstream>
#include <string>
//#include <omp.h>
#include <unistd.h>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <thread>
#include <mutex>
#include <algorithm>
#include <set>
#include "experiments.hpp"
#include <mockturtle/utils/stopwatch.hpp>
#include <random>
#include <mockturtle/algorithms/sim_resub.hpp>
#include <mockturtle/algorithms/simulation.hpp>

#include <mockturtle/algorithms/resubstitution.hpp>
//#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>

#include <iostream>

using namespace mockturtle;
using namespace hdc;
using namespace kitty;

struct XYdataset{
  std::vector<kitty::partial_truth_table> X;
  std::vector<kitty::partial_truth_table> Y;
  kitty::partial_truth_table M;
  uint64_t nin;
  uint64_t nout;
  uint64_t ndata;
};


template<typename Ntk>
double compute_accuracy( std::vector<kitty::partial_truth_table> const& X, std::vector<kitty::partial_truth_table> const& Y, kitty::partial_truth_table const& M, Ntk& ntk )
{
  std::cout << "acc" << std::endl;

  double acc = 0;
  partial_simulator sim( X );
  unordered_node_map<kitty::partial_truth_table, Ntk> node_to_value( ntk );
  simulate_nodes( ntk, node_to_value, sim );

  std::vector<kitty::partial_truth_table> V;
  std::cout << "#POs=" << ((ntk._storage->outputs).size()) << std::endl;
  for( auto i=0; i < (ntk._storage->outputs).size(); ++i ) 
  {
    V.push_back(node_to_value[ntk._storage->outputs[i]]);
    if( ntk.is_complemented( ntk._storage->outputs[i]) ) 
      V[i]=~V[i];
  }
  
  kitty::partial_truth_table diff = Y[0] | ~Y[0];
  diff = diff & (~(M&(V[0]^Y[0])));
  for( auto i=1; i < (ntk._storage->outputs).size(); ++i )
    diff = diff & (~(V[i]^Y[i]));

  /*for( auto i=0; i < (ntk._storage->outputs).size(); ++i )
  {
    kitty::print_binary( V[i] );
    std::cout << std::endl;
    kitty::print_binary( Y[i] );
    std::cout << std::endl;
  }*/

  acc = double(kitty::count_ones(diff))/diff.num_bits();
  std::cout << "end acc" << std::endl;

  return acc;
}

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

XYdataset dataset_loader( std::string file_name, uint32_t ndata )
{
  std::cout << "load" << std::endl;
  using dyn_bitset = kitty::partial_truth_table;
  XYdataset DS;
  uint32_t cnt_data = 0;
  std::string line;
  std::ifstream myfile ( file_name );
  if (myfile.is_open())
  {
    uint32_t r = 0;
    while ( getline (myfile,line) && ( cnt_data < ndata ) )
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
          DS.ndata = ndata;//std::stoi(v_line.second);
          kitty::partial_truth_table empty_bitset( DS.ndata );
          for( uint32_t i{0u}; i < DS.nin; ++i )
            DS.X.push_back( empty_bitset );
            
          for( uint32_t i{0u}; i < DS.nout; ++i )
            DS.Y.push_back( empty_bitset );
        }
      }
      else
      {
        cnt_data++;
        dyn_bitset xline( DS.nin );
        kitty::create_from_binary_string( xline, v_line.first );

        dyn_bitset yline( DS.nout+1 );
        kitty::create_from_binary_string( yline, v_line.second );

        for( uint32_t i{0u}; i < DS.nin; ++i )
          kitty::get_bit(xline,i) == 1 ? kitty::set_bit(DS.X[i],r) : kitty::clear_bit(DS.X[i],r) ;

        for( uint32_t i{0u}; i < DS.nout; ++i )
          kitty::get_bit(yline,i+1) == 1 ? kitty::set_bit(DS.Y[i],r) : kitty::clear_bit(DS.Y[i],r) ;

        DS.M.add_bit( kitty::get_bit( yline, 0 ) );
        r++;
      }
    }
    myfile.close();
  }
  else std::cout << "Unable to open file";
  
  DS.ndata = cnt_data;
  std::cout << "end load" << std::endl;
  return DS;
}

std::string DEC_ALGO{"DK_XTSD"};
using experiment_t = experiments::experiment<std::string, uint32_t, uint32_t, float, float, float, float>;
experiment_t exp_res( "/iwls2020/INTEGRATION/EX5/"+DEC_ALGO, "benchmark", "#gates", "depth", "train", "test", "valid", "runtime" );


#pragma region mutex
std::atomic<uint32_t> exp_id{0};
std::mutex exp_mutex;
#pragma endregion


#pragma region iwls2020 parameters
struct iwls2020_parameters
{
  public:
  uint32_t nImpurity{0};
  std::string dec_algo;
};
#pragma endregion iwls2020 parameters



#pragma region abc post processing
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

template<class Ntk>
xag_network abc_preprocess( Ntk const& ntk, std::string str_code, std::string abc_script = "share; resyn2rs" )
{
  xag_network res;
  write_blif( ntk, "/tmp/pre" + str_code + ".blif" );

  std::string command = "abc -q \"r /tmp/pre" + str_code + ".blif; " + abc_script + "; write_aiger /tmp/pre" + str_code + ".aig\"";

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

  std::string string_path = ( "/tmp/pre" + str_code + ".aig" );
  if( lorina::read_aiger( string_path, aiger_reader( res ) ) != lorina::return_code::success )
    std::cerr << "read_aiger failed" << std::endl;

  return res;
}

template<class Ntk>
void iterative_abc_opto( Ntk & ntk, std::string str_code, std::string abc_script = "resyn2rs"  )
{
  depth_view_params psD;
  psD.count_complements = true;
  depth_view depth_ntk_old{ntk, {}, psD};
  auto new_depth = depth_ntk_old.depth();
  auto old_depth = new_depth;
  auto new_num_gates = ntk.num_gates();
  auto old_num_gates = new_num_gates;
  bool is_first = true;
  while( is_first || (new_num_gates < old_num_gates) || (new_depth < old_depth) )
  {
    is_first = false;
    old_depth = new_depth;
    old_num_gates = new_num_gates;

    ntk = abc_opto( ntk, str_code, abc_script );
    ntk = cleanup_dangling( ntk );
    
    new_num_gates = ntk.num_gates();
    depth_view depth_ntk{ntk, {}, psD};
    new_depth = depth_ntk.depth();
  }
}
#pragma endregion abc post processing

#pragma region synthesis by high dimensional projection
template<class Ntk>
Ntk flow_hdp( std::vector<kitty::partial_truth_table>& X, std::vector<kitty::partial_truth_table>& Y, int const& topology = 1, project_params project_ps = {} )
{
  auto klut = project_in_hd( X, Y, topology, project_ps );
  Ntk ntk = convert_klut_to_graph<Ntk>( klut );
  ntk = cleanup_dangling( ntk );

  return ntk;
}
#pragma endregion region synthesis by high dimensional projection

void thread_run( iwls2020_parameters const& iwls2020_ps)
{
  std::string path_train = "../experiments/iwls2020/benchmarks/momnist/mnist60k_10_conv.txt"; //train2vDC.txt";
  std::string path_test = "../experiments/iwls2020/benchmarks/momnist/mnist10k_10_conv.txt"; //test2vDC.txt";
  //std::string path_valid = "../experiments/iwls2020/benchmarks/mnista/mnist_valid.txt";
  std::string output_path = "../experiments/iwls2020/results/MNIST/";

   std::cout << 101 << std::endl;
    auto Dl = dataset_loader( path_train, 10000 );
   std::cout << Dl.nin << " " << Dl.nout << " " << Dl.ndata << std::endl;
  
    auto Dt = dataset_loader( path_test, 10000 );
    //auto Dv = dataset_loader( path_valid );
    bool postprocess{false};

    xag_network xag;
    
    auto start = std::chrono::high_resolution_clock::now();
  
    project_params project_ps;
    project_ps.nImpurity = iwls2020_ps.nImpurity;
    xag = flow_hdp<xag_network>( Dl.X, Dl.Y, 3, project_ps );
    
    auto stop = std::chrono::high_resolution_clock::now();

    auto time_dec = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);

    depth_view d{xag};
    float la = 100*compute_accuracy( Dl.X, Dl.Y, Dl.M, d );
    float ta = 100*compute_accuracy( Dt.X, Dt.Y, Dt.M, d );
    //float va = 100*compute_accuracy( Dv.X, Dv.Y, d );
  //   fmt::print( "[i] obtained new result on {}: \n.g {}\n.d {} \n.l {} \n.t {} \n.v {}\n.c {}\n", 
  //              benchmark, xag.num_gates(), d.depth(), la, ta, va, to_seconds(res._cnt.time_dec) );
    fmt::print( "[i] obtained new result on mnist: \n.g {}\n.d {} \n.l {} \n.t {} \n.c {}", 
             xag.num_gates(), d.depth(), la, ta,to_seconds(time_dec)  );

    write_blif( xag, output_path +iwls2020_ps.dec_algo+ "mnist.blif" );

    std::ofstream myfile;
    myfile.open ( (output_path +"BLIFmnist"+iwls2020_ps.dec_algo+".txt") ); 
    myfile << ".l " << la <<std::endl;
    myfile << ".t " << ta <<std::endl;
    myfile << ".g " << xag.num_gates() << std::endl;
    myfile << ".d " << d.depth() << std::endl;
    myfile << ".c " << to_seconds(time_dec) << std::endl;

    myfile.close();
    std::cout << std::endl;

  }

int main( int argc, char* argv[] )
{
  using namespace experiments;

  iwls2020_parameters iwls2020_ps;
    if( argc == 1 )
    iwls2020_ps.nImpurity = 0;
  else if( argc == 2 )
    iwls2020_ps.nImpurity = atoi(argv[1]);

  iwls2020_ps.dec_algo = DEC_ALGO;

  std::vector<std::thread> threads;


  /* generate threads */
  for ( auto i = 0u; i < 1; ++i )
  {
    threads.emplace_back( thread_run, iwls2020_ps );
  }

  /* wait threads */
  for ( auto i = 0u; i < 1; ++i )
  {
    threads[i].join();
  }

  //exp_res.table( );
  //  exp_res.update( "best" );
  exp_res.save();
  exp_res.table();

  return 0;
}