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
#include <lorina/aiger.hpp>
#include <iostream>
#include <string>
#include <set>

using namespace mockturtle;
//using namespace mockturtle::experimental;
using namespace experiments;

template<class Ntk, class TT, class RewritingFn>
std::optional<uint32_t> database_lookup( TT const& target, RewritingFn const& rewriting_fn_, std::string const& func )
{
  using signal = typename Ntk::signal;
  std::optional<uint32_t> res;

  // first create a network 
  uint32_t num_pis = target.num_vars();
  Ntk ntk;
  std::vector<signal> pis( num_pis );
  std::generate( std::begin(pis), std::end(pis), [&](){ return ntk.create_pi(); } );
  
  signal osig;
  // run npn resynthesis
  auto const on_signal = [&]( auto const& s ) {
    uint32_t _num_nodes = mockturtle::detail::recursive_ref<Ntk>( ntk, ntk.get_node( s ) );
    mockturtle::detail::recursive_deref<Ntk>( ntk, ntk.get_node( s ) );
    if ( !res || *res > _num_nodes ) res = _num_nodes;
    osig = s;
    return true;
  };
  rewriting_fn_( ntk, target, std::begin(pis), std::end(pis), on_signal );
  
  //std::cout << ntk.num_gates() << " -> ";
  signal x = ntk.create_pi();
  signal f0;
  if( func == "xor" )
    f0 = ntk.create_xor( x, osig );
  else if( func == "and" )
    f0 = ntk.create_and( x, osig );
  else if( func == "or" )
    f0 = ntk.create_or( x, osig );
  else if( func == "lt" )
    f0 = ntk.create_lt( x, osig );
  else if( func == "le" )
    f0 = ntk.create_le( x, osig );
  else
    std::cerr << "error" << std::endl;
  ntk.create_po( f0 );
  
  //std::cout << ntk.num_gates() << std::endl;


  std::vector<kitty::partial_truth_table> pats;
  kitty::partial_truth_table sim_pat( std::pow(2,pis.size()+1) );
  for( uint32_t i = 0; i < pis.size()+1; ++i )
  {
    kitty::create_nth_var( sim_pat, i );
    pats.push_back( sim_pat );
  }

  partial_simulator sim( pats );

  unordered_node_map<kitty::partial_truth_table, Ntk> node_to_value( ntk );
  simulate_nodes( ntk, node_to_value, sim );

/*
  kitty::print_binary(node_to_value[f0]);
  std::cout << std::endl; 
*/
  std::vector<kitty::partial_truth_table*> X;
  X.push_back( &pats[0] );
  //std::cout << "inputs:" << std::endl;
  for( uint32_t i = 0; i < pats.size(); ++i )
  {
    X[0] = &pats[i];
    std::cout << kitty::mutual_information( X, &node_to_value[f0] ) << " ";

    /*std::cout << i << " I(Y" << i << ",F)=" << kitty::mutual_information( X, &node_to_value[f0] );
    std::cout << "  H(Y" << i << ")= " << kitty::entropy( X );    
    X[0] = &node_to_value[f0];
    std::cout << "  H(X^G)= " << kitty::entropy( X );
    X.push_back( &pats[i] );
    std::cout << "  H(Y" << i << ",X^G)= " << kitty::entropy( X );
    X.erase( X.begin() );
    X[0] = &node_to_value[osig];
    std::cout << "  H(G)= " << kitty::entropy( X ) << std::endl;*/

  }
  std::cout << std::endl;
  
  return res;
}


template<uint32_t num_vars>
void test_n_var_function( std::string func, bool only_npn = false )
{
  using TT = kitty::dynamic_truth_table;
  using Ntk = xag_network;
  
  // prepare the database for lookup
  xag_npn_resynthesis<Ntk> resyn;
  xag_npn_resynthesis<Ntk, Ntk, xag_npn_db_kind::aig_complete > resyn_complete;


  TT target(num_vars);
  std::set<TT> reprs;
  do
  {
    const auto repr = kitty::exact_npn_canonization(target);

    if( only_npn && ( reprs.find( std::get<0>(repr) ) != reprs.end() ) )
    {
      kitty::next_inplace( target );
      continue;
    }
    else
    {
      reprs.insert( std::get<0>(repr) );
      const auto database_res = database_lookup<Ntk, TT>( target, resyn_complete, func );
      uint32_t database_result = database_res? *database_res : std::numeric_limits<uint32_t>::max();

      
      //kitty::print_binary(target);
      //std::cout << " .r ";
      //kitty::print_hex(std::get<0>(repr) );
      
      
      //fmt::print( "\t: [db]:({}) ", 
      //            database_result);


      //std::cout << std::endl;


      //fmt::print( "{}\n", stdec_result );
      kitty::next_inplace( target );
    }
  } while ( !kitty::is_const0( target ) );
}

int main()
{

  // for ( auto const& benchmark : epfl_benchmarks() )
  // {
  //   // fmt::print( "[i] processing {}\n", benchmark );

  //   aig_network aig;
  //   auto const result = lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) );
  //   assert( result == lorina::return_code::success ); (void)result;
    
  //   uint32_t tot_mffc = 0;
  //   uint32_t big_mffc = 0;
  //   mockturtle::detail::initialize_values_with_fanout( aig );
  //   aig.foreach_node( [&]( auto n ){
  //     if ( mockturtle::detail::mffc_size( aig, n ) > 3 ) big_mffc++;
  //     tot_mffc++; 
  //   } );

  //   fmt::print( "{:.2f}\n", (double)big_mffc/tot_mffc );

  // }

  test_n_var_function<4>("xor", true);

  return 0;
}