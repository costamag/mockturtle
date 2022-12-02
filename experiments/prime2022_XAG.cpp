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

std::vector<uint32_t> IDS = { 90, 91,92,93,94,95,96,97,98,99,8,28,14,48 };

struct XYdataset{
  std::vector<kitty::partial_truth_table> X;
  kitty::partial_truth_table Y;
  uint64_t nin;
  uint64_t nout;
  uint64_t ndata;
  uint32_t conflicts_count{0};
};


template<typename Ntk>
bool simulate_input( kitty::partial_truth_table const& input_pattern, Ntk& ntk )
{
  std::vector<bool> inpt_v;

  for( uint64_t k{0u}; k<input_pattern.num_bits();++k )
  {
    inpt_v.push_back( ( ( kitty::get_bit( input_pattern, k ) == 1 ) ? true : false ) );
  }

  return simulate<bool>( ntk, default_simulator<bool>( inpt_v ) )[0];
}

template<typename Ntk>
double compute_accuracy( std::vector<kitty::partial_truth_table> const& X, kitty::partial_truth_table const& Y, Ntk& ntk )
{
  double acc = 0;
  partial_simulator sim( X );
  unordered_node_map<kitty::partial_truth_table, Ntk> node_to_value( ntk );
  simulate_nodes( ntk, node_to_value, sim );

  kitty::partial_truth_table V = node_to_value[ntk._storage->outputs[0]];
  if( ntk.is_complemented( ntk._storage->outputs[0]) ) 
    V=~V;
  acc = double(kitty::count_ones(~(V^Y)))/Y.num_bits();
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

XYdataset dataset_loader( std::string file_name )
{
  std::set<std::string> onset;
  std::set<std::string> offset;

  using dyn_bitset = kitty::partial_truth_table;
  XYdataset DS;
  DS.conflicts_count = 0;


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
        
        if( v_line.second == "0" )
        {
          if( onset.find(v_line.first) != onset.end() )
            DS.conflicts_count++;
          offset.insert( v_line.first );
        }
        else if( v_line.second == "1" )
        {
          if( offset.find(v_line.first) != offset.end() )
            DS.conflicts_count++;
          onset.insert( v_line.first );
        }
        else
          std::cerr << "[e] wrong label" << std::endl;

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

std::string DEC_ALGO{"sdec"};
using experiment_t = experiments::experiment<std::string, uint32_t, uint32_t, float, float, float, float>;
experiment_t exp_res( "/iwls2020/"+DEC_ALGO, "benchmark", "#gates", "depth", "train", "test", "valid", "runtime" );


#pragma region mutex
std::atomic<uint32_t> exp_id{0};
std::mutex exp_mutex;
#pragma endregion


#pragma region iwls2020 parameters
struct iwls2020_parameters
{
  std::string dec_algo;
  double frac_valid{0};
};
#pragma endregion iwls2020 parameters


#pragma region synthesis by high dimensional projection
template<class Ntk>
Ntk flow_hdp( std::vector<kitty::partial_truth_table>& X, std::vector<kitty::partial_truth_table>& Y, int const& topology = 1 )
{
  auto klut = project_in_hd( X, Y, topology );
  Ntk ntk = convert_klut_to_graph<Ntk>( klut );
  ntk = cleanup_dangling( ntk );

  return ntk;
}
#pragma endregion region synthesis by high dimensional projection




void thread_run( iwls2020_parameters const& iwls2020_ps, std::string const& run_only_one )
{
  std::string train_path = "../experiments/iwls2020/benchmarks/train/";
  std::string test_path = "../experiments/iwls2020/benchmarks/test/";
  std::string valid_path = "../experiments/iwls2020/benchmarks/validation/";
  std::string output_path = "../experiments/iwls2020/results/"+iwls2020_ps.dec_algo+"/";

  uint32_t id = exp_id++;


  while ( id < 100 )
  {
    //if(std::count(IDS.begin(), IDS.end(), id))
    //{
    /* read benchmark */
    std::string benchmark = fmt::format( "ex{:02}", id );
    if ( run_only_one != "" && benchmark != run_only_one )
    {
      id = exp_id++;
      continue;
    }
    std::cout << "[i] processing " << benchmark << "\n";
    
    std::string path_train = train_path+benchmark+".train.txt";
    std::string path_test = test_path+benchmark+".test.txt";
    std::string path_valid = valid_path+benchmark+".valid.txt";
      

    auto Dl = dataset_loader( path_train );
    auto Dt = dataset_loader( path_test );
    auto Dv = dataset_loader( path_valid );

  std::vector<kitty::partial_truth_table> X;
  kitty::partial_truth_table Y;

  if( iwls2020_ps.frac_valid != 0 )
  {
    for( uint32_t i = 0; i < Dl.X.size(); ++i )
    {
      for( uint32_t j = 0; j < uint32_t(iwls2020_ps.frac_valid*Dv.X[i].num_bits()); ++j )
      {
        Dl.X[i].add_bit( kitty::get_bit(Dv.X[i],j) );
      }
    }
    for( uint32_t j = 0; j < uint32_t(iwls2020_ps.frac_valid*Dv.Y.num_bits()); ++j )
    {
      Dl.Y.add_bit( kitty::get_bit(Dv.Y,j) );
    }
  }
    bool postprocess{false};

    xag_network xag;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    if( iwls2020_ps.dec_algo == "sdec" )
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 0 );
    }
    else if( iwls2020_ps.dec_algo == "isdec" )
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 1 );
    }
    else if( iwls2020_ps.dec_algo == "itsdec" )
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 2 );
    }
    else if( iwls2020_ps.dec_algo == "ixtsdec" )
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 3 );
    }
    else if( iwls2020_ps.dec_algo == "dcsdec" )
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 4 );
    }
    else if( iwls2020_ps.dec_algo == "dcxsdec" )
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 5 );
    }
    else if( iwls2020_ps.dec_algo == "muesli" )
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 6 );
    }
    else if( iwls2020_ps.dec_algo == "armuesli" )
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 7 );
    }
    else if( iwls2020_ps.dec_algo == "argmuesli" )
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 8 );
    }
    else if( iwls2020_ps.dec_algo == "fgen1024x1" )
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 9 );
    }
    else if( iwls2020_ps.dec_algo == "ifgen1024x10" ) 
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 11 );
    }
    else if( iwls2020_ps.dec_algo == "ifgen1024x10_S" )
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 12 );
    }
    else if( iwls2020_ps.dec_algo == "majgen1024x10" )
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 13 );
    }
    else if( iwls2020_ps.dec_algo == "forestS" )
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 14 );
    }
    else if( iwls2020_ps.dec_algo == "forestmuesli" )
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 15 );
    }
    else if( iwls2020_ps.dec_algo == "forestmuesli5" )
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 16 );
    }
    else if( iwls2020_ps.dec_algo == "ifgenS2048x1" ) 
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 18 );
    }
    else if( iwls2020_ps.dec_algo == "ifgenS4096x1" ) 
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 19 );
    }
    else if( iwls2020_ps.dec_algo == "ifgenS1024x2" ) 
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 20 );
    }
    else if( iwls2020_ps.dec_algo == "ifgenS1024x4" ) 
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 21 );
    }
    else if( iwls2020_ps.dec_algo == "idsdS" )
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 22 );
    }
    else if( iwls2020_ps.dec_algo == "vhds" )
    {
      auto Y = std::vector{Dl.Y};
      xag = flow_hdp<xag_network>( Dl.X, Y, 600 );
    }
    else
    {
      fmt::print( "[w] method named {} is not defined\n", iwls2020_ps.dec_algo );
    }
    
    auto stop = std::chrono::high_resolution_clock::now();

    auto time_dec = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);

    depth_view d{xag};
    float la = 100*compute_accuracy( Dl.X, Dl.Y, d );
    float ta = 100*compute_accuracy( Dt.X, Dt.Y, d );
    float va = 100*compute_accuracy( Dv.X, Dv.Y, d );
  //   fmt::print( "[i] obtained new result on {}: \n.g {}\n.d {} \n.l {} \n.t {} \n.v {}\n.c {}\n", 
  //              benchmark, xag.num_gates(), d.depth(), la, ta, va, to_seconds(res._cnt.time_dec) );
    fmt::print( "[i] obtained new result on {}: \n.g {}\n.d {} \n.l {} \n.w {} \n.t {} \n.v {}\n.c {}", 
            benchmark, xag.num_gates(), d.depth(), la, Dl.conflicts_count, ta, va,to_seconds(time_dec)  );

    exp_mutex.lock();
    exp_res( benchmark, xag.num_gates(), d.depth(), la, ta, va, to_seconds(time_dec) );
    exp_mutex.unlock();
    write_blif( xag, output_path +"BLIF/"+ benchmark + ".blif" );

    std::ofstream myfile;
    myfile.open ( (output_path +"RES/"+ benchmark + ".txt") ); 
    myfile << ".b " << fmt::format( "{:02}", id ) <<std::endl;
    myfile << ".l " << la <<std::endl;
    myfile << ".t " << ta <<std::endl;
    myfile << ".v " << va <<std::endl;
    myfile << ".g " << xag.num_gates() << std::endl;
    myfile << ".d " << d.depth() << std::endl;
    myfile << ".c " << to_seconds(time_dec) << std::endl;

    myfile.close();
    std::cout << std::endl;

    id = exp_id++;
  //}
  //else  
  //  id = exp_id++;
  }
}

int main( int argc, char* argv[] )
{
  using namespace experiments;

  iwls2020_parameters iwls2020_ps;
  iwls2020_ps.dec_algo = DEC_ALGO;
  iwls2020_ps.frac_valid = 0;

  std::string run_only_one = "";

  if ( argc == 2 )
    run_only_one = std::string( argv[1] );

  const auto processor_count = run_only_one != "" ? 1 : std::thread::hardware_concurrency();

  /* starting benchmark id */
  exp_id.store( 0 );

  std::vector<std::thread> threads;

  /* generate threads */
  fmt::print( "[i] Running on {} threads\n", processor_count );
  for ( auto i = 0u; i < processor_count; ++i )
  {
    threads.emplace_back( thread_run, iwls2020_ps, run_only_one );
  }

  /* wait threads */
  for ( auto i = 0u; i < processor_count; ++i )
  {
    threads[i].join();
  }

  //exp_res.table( );
  //  exp_res.update( "best" );
  exp_res.save();
  exp_res.table();

  return 0;
}