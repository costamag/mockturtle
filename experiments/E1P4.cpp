#include "experiments.hpp"
#include <mockturtle/traits.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/resyn_engines/xag_resyn.hpp>
//#include <mockturtle/algorithms/experimental/contest.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/detail/mffc_utils.hpp>
#include <mockturtle/algorithms/lfe/hyperdimensional_computing/model.hpp>
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
using namespace hdc;
//using namespace mockturtle::experimental;
using namespace experiments;


template<class Ntk, class TT, class RewritingFn>
std::vector<uint32_t> database_lookup( TT const& target, RewritingFn const& rewriting_fn_ )
{
  std::vector<uint32_t> result;
  std::optional<uint32_t> res;
  // first create a network 
  uint32_t num_pis = target.num_vars();
  Ntk ntk;
  std::vector<typename Ntk::signal> pis( num_pis );
  std::generate( std::begin(pis), std::end(pis), [&](){ return ntk.create_pi(); } );
  
  typename Ntk::signal osig;
  // run npn resynthesis
  auto const on_signal = [&]( auto const& s ) {
    uint32_t _num_nodes = mockturtle::detail::recursive_ref<Ntk>( ntk, ntk.get_node( s ) );
    mockturtle::detail::recursive_deref<Ntk>( ntk, ntk.get_node( s ) );
    if ( !res || *res > _num_nodes ) res = _num_nodes;
    osig = s;
    return true;
  };
  rewriting_fn_( ntk, target, std::begin(pis), std::end(pis), on_signal );
  
  ntk.create_po( osig );
  result.push_back(ntk.num_gates());
  
  //std::cout << ntk.num_gates() << std::endl;


  std::vector<kitty::partial_truth_table> pats;
  kitty::partial_truth_table sim_pat( std::pow(2,pis.size()) );
  for( uint32_t i = 0; i < pis.size(); ++i )
  {
    kitty::create_nth_var( sim_pat, i );
    pats.push_back( sim_pat );
  }

  partial_simulator sim( pats );

  unordered_node_map<kitty::partial_truth_table, Ntk> node_to_value( ntk );
  simulate_nodes( ntk, node_to_value, sim );

  std::vector<kitty::partial_truth_table> Y = {node_to_value[osig]};
  
  klut_network oklut;
  simulation_view oklut_sim{ oklut };

  model M( oklut_sim, pats, Y );
  std::vector<klut_network::signal> osignals;
  
  result.push_back(0);

  hdc::detail::selection_method selection_m = hdc::detail::selection_method::depth_selector;
  hdc::detail::selection_params selection_ps;
  selection_ps.max_new_supports = 6;
  selection_ps.max_selection_attempts = 30;
  selection_ps.support_size = 2;
  selection_ps.max_search_depth = std::numeric_limits<uint32_t>::max();

  std::vector<std::vector<klut_network::signal>> supports = select_variables( oklut_sim, selection_m, selection_ps );

  hdc::detail::creation_method creation_m = hdc::detail::creation_method::ifgenerator1;
  hdc::detail::creation_params creation_ps;
  creation_ps.max_nodes_total  = 6000;

  create_nodes( oklut_sim, supports, creation_m, creation_ps );
  hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::ixtsdec;
  hdc::detail::arecovery_params arecovery_ps;
  arecovery_ps.output=0;
  osignals.push_back( M.accuracy_recovery(arecovery_m, arecovery_ps) );

  M.ntk_.create_po( osignals[0] );

  result.push_back(M.ntk_.num_gates());
  M.print_summary();

  
  for( uint32_t i = 0; i < result.size(); ++i )
    std::cout << result[i] << " ";
  std::cout << std::endl; 
  return result;
}


template<uint32_t num_vars>
std::vector<std::vector<uint32_t>> synthesize_nf( bool only_npn = false )
{
  using TT = kitty::dynamic_truth_table;
  using Ntk = xag_network;
  
  // prepare the database for lookup
  xag_npn_resynthesis<Ntk> resyn;
  xag_npn_resynthesis<Ntk, Ntk, xag_npn_db_kind::xag_complete> resyn_complete;

  std::vector<std::vector<uint32_t>> num_gates;
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
      kitty::print_hex(target);
      std::cout << std::endl;
      reprs.insert( std::get<0>(repr) );
      num_gates.push_back(database_lookup<Ntk, TT>( target, resyn_complete ));

      kitty::next_inplace( target );
    }
  } while ( !kitty::is_const0( target ) );
  return num_gates;
}

int main()
{
  auto npn_fractions = synthesize_nf<3>(true);
  for( auto f : npn_fractions )
  {
    for( uint32_t i = 0; i < f.size(); ++i )
      std::cout << f[i] << " ";
    std::cout << std::endl;
  }

  return 0;
}