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
#include <omp.h>
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
  double delta_acc;
  for( uint64_t k {0u}; k < X[0].num_bits(); ++k )
  {
    kitty::partial_truth_table ipattern;
    for ( uint64_t j {0u}; j < X.size(); ++j )
      ipattern.add_bit( kitty::get_bit(X[j],k) );
            
    delta_acc = ( ( simulate_input( ipattern, ntk ) == kitty::get_bit(Y,k) ) ? (double)1.0/X[0].num_bits() : 0.0 );
    acc += delta_acc;

  }
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

std::string DEC_ALGO{"ifgen1024x1"};
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
aig_network abc_preprocess( Ntk const& ntk, std::string str_code, std::string abc_script = "share; resyn2rs" )
{
  aig_network res;
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
  uint64_t nin;
  uint64_t nout;
  uint64_t ndata;
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

    aig_network aig;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    if( iwls2020_ps.dec_algo == "sdec" )
    {
      auto Y = std::vector{Dl.Y};
      aig = flow_hdp<aig_network>( Dl.X, Y, 0 );
    }
    else if( iwls2020_ps.dec_algo == "isdec" )
    {
      auto Y = std::vector{Dl.Y};
      aig = flow_hdp<aig_network>( Dl.X, Y, 1 );
    }
    else if( iwls2020_ps.dec_algo == "itsdec" )
    {
      auto Y = std::vector{Dl.Y};
      aig = flow_hdp<aig_network>( Dl.X, Y, 2 );
    }
    else if( iwls2020_ps.dec_algo == "ixtsdec" )
    {
      auto Y = std::vector{Dl.Y};
      aig = flow_hdp<aig_network>( Dl.X, Y, 3 );
    }
    else if( iwls2020_ps.dec_algo == "dcsdec" )
    {
      auto Y = std::vector{Dl.Y};
      aig = flow_hdp<aig_network>( Dl.X, Y, 4 );
    }
    else if( iwls2020_ps.dec_algo == "dcxsdec" )
    {
      auto Y = std::vector{Dl.Y};
      aig = flow_hdp<aig_network>( Dl.X, Y, 5 );
    }
    else if( iwls2020_ps.dec_algo == "muesli" )
    {
      auto Y = std::vector{Dl.Y};
      aig = flow_hdp<aig_network>( Dl.X, Y, 6 );
    }
    else if( iwls2020_ps.dec_algo == "armuesli" )
    {
      auto Y = std::vector{Dl.Y};
      aig = flow_hdp<aig_network>( Dl.X, Y, 7 );
    }
    else if( iwls2020_ps.dec_algo == "argmuesli" )
    {
      auto Y = std::vector{Dl.Y};
      aig = flow_hdp<aig_network>( Dl.X, Y, 8 );
    }
    else if( iwls2020_ps.dec_algo == "fgen1024x1" )
    {
      auto Y = std::vector{Dl.Y};
      aig = flow_hdp<aig_network>( Dl.X, Y, 9 );
    }
    else if( iwls2020_ps.dec_algo == "ifgen1024x1" )
    {
      auto Y = std::vector{Dl.Y};
      aig = flow_hdp<aig_network>( Dl.X, Y, 10 );
    }
    else
    {
      fmt::print( "[w] method named {} is not defined\n", iwls2020_ps.dec_algo );
    }
    
    auto stop = std::chrono::high_resolution_clock::now();

    auto time_dec = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);

    if( postprocess )
      iterative_abc_opto( aig, benchmark, "resyn2rs" );

    depth_view d{aig};
    float la = 100*compute_accuracy( Dl.X, Dl.Y, d );
    float ta = 100*compute_accuracy( Dt.X, Dt.Y, d );
    float va = 100*compute_accuracy( Dv.X, Dv.Y, d );
  //   fmt::print( "[i] obtained new result on {}: \n.g {}\n.d {} \n.l {} \n.t {} \n.v {}\n.c {}\n", 
  //              benchmark, aig.num_gates(), d.depth(), la, ta, va, to_seconds(res._cnt.time_dec) );
    fmt::print( "[i] obtained new result on {}: \n.g {}\n.d {} \n.l {} \n.w {} \n.t {} \n.v {}\n.c {}", 
            benchmark, aig.num_gates(), d.depth(), la, Dl.conflicts_count, ta, va,to_seconds(time_dec)  );

    exp_mutex.lock();
    exp_res( benchmark, aig.num_gates(), d.depth(), la, ta, va, to_seconds(time_dec) );
    exp_mutex.unlock();
    write_aiger( aig, output_path +"AIG/"+ benchmark + ".aig" );

    std::ofstream myfile;
    myfile.open ( (output_path +"RES/"+ benchmark + ".txt") ); 
    myfile << ".b " << fmt::format( "{:02}", id ) <<std::endl;
    myfile << ".l " << la <<std::endl;
    myfile << ".t " << ta <<std::endl;
    myfile << ".v " << va <<std::endl;
    myfile << ".g " << aig.num_gates() << std::endl;
    myfile << ".d " << d.depth() << std::endl;
    myfile << ".c " << to_seconds(time_dec) << std::endl;

    myfile.close();
    std::cout << std::endl;

    id = exp_id++;
  }
}

int main( int argc, char* argv[] )
{
  using namespace experiments;

  iwls2020_parameters iwls2020_ps;
  iwls2020_ps.dec_algo = DEC_ALGO;
  iwls2020_ps.frac_valid = 1;

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