#include "experiments.hpp"
#include <mockturtle/traits.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/resyn_engines/xag_resyn.hpp>
//#include <mockturtle/algorithms/experimental/contest.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/detail/mffc_utils.hpp>
#include <mockturtle/algorithms/lfe/sim_muesli.hpp>
#include <mockturtle/algorithms/lfe/hyperdimensional_computing/model.hpp>
//#include <mockturtle/algorithms/lfe/sim_decomposition.hpp>
//#include <mockturtle/algorithms/lfe/sim_decomposition_dp.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <lorina/aiger.hpp>
#include <kitty/statistics.hpp>
#include <iostream>
#include <string>
#include <set>

using namespace mockturtle;
using namespace hdc;
//using namespace mockturtle::experimental;
using namespace experiments;


template<class TT>
void create_bottomdec( TT const& oldtarget, TT const& target )
{
  std::vector<uint32_t> result;
  std::optional<uint32_t> res;
  // first create a network 
  klut_network ntk1;
  
  auto X0 = ntk1.create_pi();
  auto X1 = ntk1.create_pi();
  auto X2 = ntk1.create_pi();
  std::vector<klut_network::signal> X = { X0, X1, X2 };
  auto Xi = ntk1.create_pi();
  auto Xj = ntk1.create_pi();
  auto Xn = ntk1.create_and( Xi, Xj );
  auto f0 = ntk1.create_node( X, oldtarget );
  auto f1 = ntk1.create_node( X, target );
  auto fo = ntk1.create_ite( Xn, f1, f0 );
  ntk1.create_po( fo );

  std::vector<kitty::partial_truth_table> pats;
  kitty::partial_truth_table sim_pat( 32u );
  for( uint32_t i = 0; i < 5; ++i )
  {
    kitty::create_nth_var( sim_pat, i );
    pats.push_back( sim_pat );
  }

  partial_simulator sim( pats );

  unordered_node_map<kitty::partial_truth_table, klut_network> node_to_value( ntk1 );
  simulate_nodes( ntk1, node_to_value, sim );

  kitty::partial_truth_table Y = node_to_value[fo];

  std::vector<kitty::partial_truth_table*> Xptr = { &node_to_value[Xi] };
  std::cout << kitty::mutual_information( Xptr, &Y ) << " ";

  Xptr = { &node_to_value[Xj] };
  std::cout << kitty::mutual_information( Xptr, &Y ) << " - ";

  Xptr = { &node_to_value[Xi], &node_to_value[Xj] };
  std::cout << kitty::mutual_information( Xptr, &Y ) << " ";

  Xptr = { &node_to_value[Xn] };
  std::cout << kitty::mutual_information( Xptr, &Y ) << " ";

  Xptr = { &node_to_value[Xi], &node_to_value[Xn] };
  std::cout << kitty::mutual_information( Xptr, &Y ) << " ";

  Xptr = { &node_to_value[Xj], &node_to_value[Xn] };
  std::cout << kitty::mutual_information( Xptr, &Y ) << " ";

  Xptr = { &node_to_value[Xi], &node_to_value[Xj], &node_to_value[Xn] };
  std::cout << kitty::mutual_information( Xptr, &Y ) << std::endl;

}


void check_bottomdec()
{
  using TT = kitty::dynamic_truth_table;
  
  // prepare the database for lookup

  std::vector<std::vector<uint32_t>> num_gates;
  TT target(3);
  TT oldtarget(3);
  std::set<TT> reprs;
  do
  {
    const auto repr = kitty::exact_npn_canonization(target);
    if( ( reprs.find( std::get<0>(repr) ) != reprs.end() ) && ( reprs.find( std::get<0>(repr) ) != reprs.end() ) )
    {
      kitty::next_inplace( target );
      continue;
    }
    else
    {
      std::set<TT> reprsOLD;
    /////////////////////////
      do
      {
        const auto reprOLD = kitty::exact_npn_canonization(oldtarget);
        std::set<TT> reprsOLD;
    
      if( ( reprs.find( std::get<0>(reprOLD) ) != reprsOLD.end() ) && ( reprsOLD.find( std::get<0>(reprOLD) ) != reprsOLD.end() ) )
      {
        kitty::next_inplace( oldtarget );
        continue;
      }
      else
      {
        reprsOLD.insert( std::get<0>(reprOLD) );
        create_bottomdec<TT>( oldtarget, target );
        kitty::next_inplace( oldtarget );
      }
    } while ( !kitty::is_const0( target ) );
    /////////////////////////
      reprs.insert( std::get<0>(repr) );
      kitty::next_inplace( target );
    }

  } while ( !kitty::is_const0( target ) );
}

int main()
{
  check_bottomdec();
  return 0;
}