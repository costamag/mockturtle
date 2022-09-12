#include "experiments.hpp"
#include <mockturtle/traits.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/resyn_engines/xag_resyn.hpp>
//#include <mockturtle/algorithms/experimental/contest.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/detail/mffc_utils.hpp>
#include <mockturtle/algorithms/lfe/muesli.hpp>
#include <mockturtle/algorithms/lfe/sim_muesli.hpp>
//#include <mockturtle/algorithms/lfe/sim_decomposition.hpp>
//#include <mockturtle/algorithms/lfe/sim_decomposition_dp.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <kitty/print.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/partial_truth_table.hpp>
#include <lorina/aiger.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <set>

#include <sys/types.h>
#include <sys/stat.h>

using namespace mockturtle;
//using namespace mockturtle::experimental;
using namespace experiments;

int binomialCoefficients(int n, int k) {
   if (k == 0 || k == n)
   return 1;
   return binomialCoefficients(n - 1, k - 1) + binomialCoefficients(n - 1, k);
}

std::vector<int> combination;


std::vector<std::vector<int>> combinations; 

void go(std::vector<uint32_t> Bits, int offset, int k) {
  if (k == 0) {
    combinations.push_back( combination );
    //pretty_print(combination);
    return;
  }
  for (int i = offset; i <= Bits.size() - k; ++i) {
    combination.push_back(Bits[i]);
    go(Bits, i+1, k-1);
    combination.pop_back();
  }
}

double erase_and_print( std::vector<kitty::partial_truth_table> pats, kitty::partial_truth_table Y, uint32_t num_erase, std::string path )
{
  std::vector<sim_pattern<klut_network>> Xsp;
  for( auto tt : pats )
  {
    sim_pattern<klut_network> xsp( tt );
    Xsp.push_back(xsp);
  } 
  std::vector<uint32_t> reduced_support = {0,1,2};


  std::vector<uint32_t> erasable_bits;
  for (int i = 0; i < Y.num_bits(); ++i) { erasable_bits.push_back(i); }
  combinations = {};

  go( erasable_bits, 0, num_erase );
  
  uint32_t count = 0;
  uint32_t CONTA=0;
  std::string new_path = path+"/"+std::to_string(num_erase);
  
  int status1 = mkdir( (new_path).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH ); 

  int Bin = binomialCoefficients( 16 , num_erase );
  int delta = std::ceil( (double)Bin/20 );

  std::cout << Bin << " => delta=" << delta << std::endl;


  for( auto x : combinations )
  { 

    if( CONTA % delta != 0 )
    {
      CONTA++;
      continue;
    }
    std::cout << CONTA << " ";

    kitty::partial_truth_table on_f( Y.num_bits() ), amask1, amask0;
    on_f = Y;
    //kitty::print_binary( on_f ); std::cout << std::endl;
    amask1 = pats[3];
    //kitty::print_binary( amask1 ); std::cout << std::endl;
    amask0 = ~pats[3];
    //kitty::print_binary( amask0 ); std::cout << std::endl;

    for( uint32_t k = 0; k < x.size(); ++k )
    {
      kitty::clear_bit( amask0, x[k] );
      kitty::clear_bit( amask1, x[k] );
      kitty::clear_bit( on_f, x[k] );
    }

    std::reverse(x.begin(), x.end());

      kitty::partial_truth_table ytt = Y;
      std::vector<kitty::partial_truth_table> Xtt = pats;
      for( uint32_t j = 0; j < pats.size(); ++j )
      {
        for( uint32_t k = 0; k < x.size(); ++k )
          Xtt[j].erase_bit_shift( x[k] );
      }
  
      std::vector<kitty::partial_truth_table*> X;
      X.push_back( &Xtt[0] );
  

      for( uint32_t k = 0; k < x.size(); ++k )
        ytt.erase_bit_shift( x[k] );
      
      std::ofstream myfile;
      myfile.open (new_path+"/ex"+std::to_string(CONTA++)+".pla");
      myfile << ".i 4\n";
      myfile << ".o 1\n";
      std::string numpats_string = ".p " + std::to_string(ytt.num_bits()) + "\n";
      myfile << ".type fr\n";
      myfile << numpats_string;

      std::string row;

      for( uint32_t r=0; r < ytt.num_bits(); ++r )
      {
        row="";
        for( uint32_t c=0; c < Xtt.size(); ++c )
          row=std::to_string( kitty::get_bit( Xtt[c], r ) )+row;
        row = row + " " + std::to_string( kitty::get_bit( ytt, r ) ) + "\n";
      myfile << row;
      }
      myfile << ".e";
      myfile.close();
     
    }
  std::cout << std::endl;
  return std::round((double)count/(double)combinations.size() * 100.0) / 100.0;
}


template<class TT, class RewritingFn>
void print_pla( TT const& target, RewritingFn const& rewriting_fn_, std::string const& path  )
{
  using signal = typename aig_network::signal;
  std::optional<uint32_t> res;

  // first create a network 
  uint32_t num_pis = target.num_vars();
  aig_network ntk;
  std::vector<signal> pis( num_pis );
  std::generate( std::begin(pis), std::end(pis), [&](){ return ntk.create_pi(); } );
  
  signal osig;
  // run npn resynthesis
  auto const on_signal = [&]( auto const& s ) {
    uint32_t _num_nodes = mockturtle::detail::recursive_ref<aig_network>( ntk, ntk.get_node( s ) );
    mockturtle::detail::recursive_deref<aig_network>( ntk, ntk.get_node( s ) );
    if ( !res || *res > _num_nodes ) res = _num_nodes;
    osig = s;
    return true;
  };
  rewriting_fn_( ntk, target, std::begin(pis), std::end(pis), on_signal );

  std::vector<kitty::partial_truth_table> pats;
  kitty::partial_truth_table sim_pat( std::pow(2,pis.size()) );
  for( uint32_t i = 0; i < pis.size(); ++i )
  {
    kitty::create_nth_var( sim_pat, i );
    pats.push_back( sim_pat );
  }

  partial_simulator sim( pats );

  unordered_node_map<kitty::partial_truth_table, aig_network> node_to_value( ntk );
  simulate_nodes( ntk, node_to_value, sim );

  double fraction_right;
  std::vector<double> fractions;

  for( uint32_t i = 0; i < pats[0].num_bits(); ++i )
  {
    kitty::partial_truth_table Y = node_to_value[osig];
    
    erase_and_print( pats, Y, i, path );
  }
}


template<uint32_t num_vars>
void test_n_var_function( )
{
  using TT = kitty::dynamic_truth_table;
  using Ntk = aig_network;
  
  // prepare the database for lookup
  xag_npn_resynthesis<Ntk> resyn;
  xag_npn_resynthesis<Ntk, Ntk, xag_npn_db_kind::aig_complete > resyn_complete;
  

  std::vector<std::vector<double>> npn_fractions;
  TT target(num_vars);
  std::set<TT> reprs;
  do
  {
    const auto repr = kitty::exact_npn_canonization(target);

    if( reprs.find( std::get<0>(repr) ) != reprs.end() ) 
    {
      kitty::next_inplace( target );
      continue;
    }
    else
    {
      reprs.insert( std::get<0>(repr) );
      std::string TTref = kitty::to_hex( std::get<0>(repr) );
      //const char * direct = TTref;
      int status = mkdir( ("PLAS/"+TTref).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH ); 
      print_pla<TT>( target, resyn_complete, "PLAS/"+TTref );

      kitty::next_inplace( target );
    }
  } while ( !kitty::is_const0( target ) );
}

int main()
{
  int status = mkdir( "PLAS/", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH ); 
  
  test_n_var_function<4>();

  return 0;
}