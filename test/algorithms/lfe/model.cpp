#include <catch.hpp>

#include <mockturtle/algorithms/lfe/simulation_view.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/algorithms/lfe/hyperdimensional_computing/model.hpp>
#include <mockturtle/algorithms/lfe/hyperdimensional_computing/methods/selgenerators.hpp>
#include <mockturtle/algorithms/lfe/hyperdimensional_computing/methods/selectors.hpp>
#include <mockturtle/algorithms/lfe/hyperdimensional_computing/methods/generators.hpp>
#include <mockturtle/algorithms/lfe/hyperdimensional_computing/methods/accuracy_recovery.hpp>
#include <kitty/partial_truth_table.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <kitty/constructors.hpp>

using namespace kitty;
using namespace mockturtle;
using namespace hdc;

TEST_CASE( "muesli f=ab+cde", "[model]" )
{
  klut_network oklut;
  simulation_view oklut_sim{ oklut };
  std::vector<partial_truth_table> ex;
  partial_truth_table tt(32u);
  for( uint32_t i = 0; i < 5; ++i )
  {
    create_nth_var(tt, i);
    ex.push_back(tt);
  }

  tt = (ex[0]&ex[1])|(ex[2]&(ex[3]&ex[4]));
  std::vector<partial_truth_table> targets = {tt};
  model M( oklut_sim, ex, targets );
  std::vector<signal<klut_network>> osignals;
  hdc::detail::selcreation_method selcreation_m = hdc::detail::selcreation_method::muesli;
  hdc::detail::selcreation_params selcreation_ps;
  selcreation_ps.re_initialize = false;
  selcreation_ps.verbose = true;
      
  selcreation_ps.output=0;
  M.add( selcreation_m, selcreation_ps );

  hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::none;
  hdc::detail::arecovery_params arecovery_ps;
  arecovery_ps.verbose = false;
  
  arecovery_ps.output=0;
  osignals.push_back( M.accuracy_recovery(arecovery_m, arecovery_ps) );

  M.ntk_.create_po( osignals[0] );

  klut_network fklut;
  fklut = M.ntk_;

  partial_simulator sim( ex );
  unordered_node_map<partial_truth_table, klut_network> node_to_value( fklut );
  simulate_nodes( fklut, node_to_value, sim );

  CHECK( node_to_value[osignals[0]]  == tt );
  M.print_summary();

}

TEST_CASE( "muesli f=ab+cd", "[model]" )
{
  klut_network oklut;
  simulation_view oklut_sim{ oklut };
  std::vector<partial_truth_table> ex;
  partial_truth_table tt(32u);
  for( uint32_t i = 0; i < 5; ++i )
  {
    create_nth_var(tt, i);
    ex.push_back(tt);
  }

  tt = (ex[0]&ex[1])|(ex[2]&ex[3]);
  std::vector<partial_truth_table> targets = {tt};
  model M( oklut_sim, ex, targets );
  std::vector<signal<klut_network>> osignals;
  hdc::detail::selcreation_method selcreation_m = hdc::detail::selcreation_method::muesli;
  hdc::detail::selcreation_params selcreation_ps;
  selcreation_ps.re_initialize = false;
  selcreation_ps.verbose = false;
      
  selcreation_ps.output=0;
  M.add( selcreation_m, selcreation_ps );

  hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::none;
  hdc::detail::arecovery_params arecovery_ps;
  arecovery_ps.verbose = false;
  
  arecovery_ps.output=0;
  osignals.push_back( M.accuracy_recovery(arecovery_m, arecovery_ps) );

  M.ntk_.create_po( osignals[0] );

  klut_network fklut;
  fklut = M.ntk_;

  partial_simulator sim( ex );
  unordered_node_map<partial_truth_table, klut_network> node_to_value( fklut );
  simulate_nodes( fklut, node_to_value, sim );

  CHECK( node_to_value[osignals[0]]  == tt );


}

TEST_CASE( "muesli f=(a^b)c", "[model]" )
{
  klut_network oklut;
  simulation_view oklut_sim{ oklut };
  std::vector<partial_truth_table> ex;
  partial_truth_table tt(8u);
  for( uint32_t i = 0; i < 3; ++i )
  {
    create_nth_var(tt, i);
    ex.push_back(tt);
  }

  tt = (ex[0]^ex[1])&ex[2];
  std::vector<partial_truth_table> targets = {tt};
  model M( oklut_sim, ex, targets );
  std::vector<signal<klut_network>> osignals;
  hdc::detail::selcreation_method selcreation_m = hdc::detail::selcreation_method::muesli;
  hdc::detail::selcreation_params selcreation_ps;
  selcreation_ps.re_initialize = false;
  selcreation_ps.verbose = true;
      
  selcreation_ps.output=0;
  M.add( selcreation_m, selcreation_ps );

  hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::none;
  hdc::detail::arecovery_params arecovery_ps;
  arecovery_ps.verbose = true;
  
  arecovery_ps.output=0;
  osignals.push_back( M.accuracy_recovery(arecovery_m, arecovery_ps) );

  M.ntk_.create_po( osignals[0] );

  klut_network fklut;
  fklut = M.ntk_;

  partial_simulator sim( ex );
  unordered_node_map<partial_truth_table, klut_network> node_to_value( fklut );
  simulate_nodes( fklut, node_to_value, sim );

  CHECK( node_to_value[osignals[0]]  == tt );

  std::cout << std::endl;
  for( uint32_t i = 0; i < M.ntk_.layer_to_signals.size(); ++i )
  {
    std::cout << "layer " << i << ": ";
    for( uint32_t j = 0; j < M.ntk_.layer_to_signals[i].size(); ++j )
      std::cout << M.ntk_.layer_to_signals[i][j] << " ";
    std::cout << std::endl;

  }
  M.print_summary();
}

TEST_CASE( "decomposition f=(a^b)c", "[model]" )
{
  klut_network oklut;
  simulation_view oklut_sim{ oklut };
  std::vector<partial_truth_table> ex;
  partial_truth_table tt(8u);
  for( uint32_t i = 0; i < 3; ++i )
  {
    create_nth_var(tt, i);
    ex.push_back(tt);
  }

  tt = (ex[0]^ex[1])&ex[2];
  std::vector<partial_truth_table> targets = {tt};
  model M( oklut_sim, ex, targets );
  std::vector<signal<klut_network>> osignals;

  hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::ixtsdecS;
  hdc::detail::arecovery_params arecovery_ps;
  arecovery_ps.verbose = true;
  
  arecovery_ps.output=0;
  osignals.push_back( M.accuracy_recovery(arecovery_m, arecovery_ps) );

  M.ntk_.create_po( osignals[0] );

  klut_network fklut;
  fklut = M.ntk_;

  partial_simulator sim( ex );
  unordered_node_map<partial_truth_table, klut_network> node_to_value( fklut );
  simulate_nodes( fklut, node_to_value, sim );

  CHECK( node_to_value[osignals[0]]  == tt );

  std::cout << std::endl;
  for( uint32_t i = 0; i < M.ntk_.layer_to_signals.size(); ++i )
  {
    std::cout << "layer " << i << ": ";
    for( uint32_t j = 0; j < M.ntk_.layer_to_signals[i].size(); ++j )
      std::cout << M.ntk_.layer_to_signals[i][j] << " ";
    std::cout << std::endl;

  }
  M.print_summary();
}


TEST_CASE( "selection", "[selector]" )
{
  klut_network oklut;
  simulation_view oklut_sim{ oklut };
  std::vector<partial_truth_table> ex;
  partial_truth_table tt(16u);
  for( uint32_t i = 0; i < 4; ++i )
  {
    create_nth_var(tt, i);
    ex.push_back(tt);
  }

  tt = (ex[0]^ex[1])|(ex[2]&ex[3]);
  std::vector<partial_truth_table> targets = {tt};
  model M( oklut_sim, ex, targets );
  std::vector<signal<klut_network>> osignals;

  hdc::detail::selection_method selection_m = hdc::detail::selection_method::depth_selector;
  hdc::detail::selection_params selection_ps;
  selection_ps.max_new_supports = 3;
  selection_ps.max_selection_attempts = 10;
  selection_ps.support_size = 2;
  selection_ps.max_search_depth = std::numeric_limits<uint32_t>::max();

  std::vector<std::vector<signal<klut_network>>> supports = select_variables( oklut_sim, selection_m, selection_ps );

  CHECK( supports.size() == 3 );

  selection_ps.max_new_supports = 10;
  selection_ps.max_selection_attempts = 50;
  supports = select_variables( oklut_sim, selection_m, selection_ps );

  CHECK( supports.size() == 6 );

  oklut_sim.create_and(supports[0][0],supports[0][1]);
  oklut_sim.create_and(supports[1][0],supports[1][1]);
  oklut_sim.create_and(supports[2][0],supports[2][1]);
  oklut_sim.create_and(supports[3][0],supports[3][1]);
  oklut_sim.create_and(supports[4][0],supports[4][1]);
  oklut_sim.create_and(supports[5][0],supports[5][1]);

  selection_ps.max_new_supports = 3;
  selection_ps.max_selection_attempts = 10;
  selection_ps.support_size = 2;
  selection_ps.max_search_depth = 1;
  
  supports = select_variables( oklut_sim, selection_m, selection_ps );
  bool is_correct_layer{true};
  for( uint32_t i = 0; i < supports.size(); ++i )
  {
    for( uint32_t j = 0; j < supports[i].size(); ++j )
      is_correct_layer &= M.ntk_.nodes_to_layer[M.ntk_.get_node(supports[i][j])] == ( M.ntk_.layer_to_signals.size() - 1 ) ;
  }
  CHECK( is_correct_layer == true );

  selection_ps.max_new_supports = 10;
  selection_ps.max_selection_attempts = 50;
  selection_ps.support_size = 2;
  selection_ps.max_search_depth = 2;
  
  supports = select_variables( oklut_sim, selection_m, selection_ps );
  is_correct_layer = true;
  for( uint32_t i = 0; i < supports.size(); ++i )
  {
    for( uint32_t j = 0; j < supports[i].size(); ++j )
      is_correct_layer &= M.ntk_.nodes_to_layer[M.ntk_.get_node(supports[i][j])] == ( M.ntk_.layer_to_signals.size() - 1 ) ;
  }
  CHECK( is_correct_layer == false );

  selection_m = hdc::detail::selection_method::layer_selector;
  selection_ps.max_new_supports = 4;
  selection_ps.max_selection_attempts = 30;
  selection_ps.support_size = 2;
  selection_ps.layer = 0;
    
  supports = select_variables( oklut_sim, selection_m, selection_ps );
  is_correct_layer = true;
  for( uint32_t i = 0; i < supports.size(); ++i )
  {
    for( uint32_t j = 0; j < supports[i].size(); ++j )
      is_correct_layer &= M.ntk_.nodes_to_layer[M.ntk_.get_node(supports[i][j])] == 0 ;
  }
  CHECK( is_correct_layer == true );
  CHECK( supports.size() == 4 );

  selection_ps.layer = 1;
    
  supports = select_variables( oklut_sim, selection_m, selection_ps );
  is_correct_layer = true;
  for( uint32_t i = 0; i < supports.size(); ++i )
  {
    for( uint32_t j = 0; j < supports[i].size(); ++j )
      is_correct_layer &= M.ntk_.nodes_to_layer[M.ntk_.get_node(supports[i][j])] == 1 ;
  }
  CHECK( is_correct_layer == true );
  CHECK( supports.size() == 4 );

}


TEST_CASE( "selection and creation", "[selector]" )
{
  klut_network oklut;
  simulation_view oklut_sim{ oklut };
  std::vector<partial_truth_table> ex;
  partial_truth_table tt(16u);
  for( uint32_t i = 0; i < 4; ++i )
  {
    create_nth_var(tt, i);
    ex.push_back(tt);
  }

  tt = (ex[0]^ex[1])|(ex[2]&ex[3]);
  std::vector<partial_truth_table> targets = {tt};
  model M( oklut_sim, ex, targets );
  std::vector<signal<klut_network>> osignals;

  hdc::detail::selection_method selection_m = hdc::detail::selection_method::depth_selector;
  hdc::detail::selection_params selection_ps;
  selection_ps.max_new_supports = 3;
  selection_ps.max_selection_attempts = 10;
  selection_ps.support_size = 2;
  selection_ps.max_search_depth = std::numeric_limits<uint32_t>::max();

  std::vector<std::vector<signal<klut_network>>> supports = select_variables( oklut_sim, selection_m, selection_ps );

  hdc::detail::creation_method creation_m = hdc::detail::creation_method::fgenerator1;
  hdc::detail::creation_params creation_ps;

  creation_ps.max_nodes_total  = 3;
  creation_ps.max_nodes_support = 1;

  create_nodes( oklut_sim, supports, creation_m, creation_ps );

  CHECK( M.ntk_.num_gates() == 3u );

}

TEST_CASE( "selection and creation more nodes", "[selector]" )
{
  klut_network oklut;
  simulation_view oklut_sim{ oklut };
  std::vector<partial_truth_table> ex;
  partial_truth_table tt(16u);
  for( uint32_t i = 0; i < 4; ++i )
  {
    create_nth_var(tt, i);
    ex.push_back(tt);
  }

  tt = (ex[0]^ex[1])|(ex[2]&ex[3]);
  std::vector<partial_truth_table> targets = {tt};
  model M( oklut_sim, ex, targets );
  std::vector<signal<klut_network>> osignals;

  hdc::detail::selection_method selection_m = hdc::detail::selection_method::depth_selector;
  hdc::detail::selection_params selection_ps;
  selection_ps.max_new_supports = 3;
  selection_ps.max_selection_attempts = 10;
  selection_ps.support_size = 2;
  selection_ps.max_search_depth = std::numeric_limits<uint32_t>::max();

  std::vector<std::vector<signal<klut_network>>> supports = select_variables( oklut_sim, selection_m, selection_ps );

  hdc::detail::creation_method creation_m = hdc::detail::creation_method::fgenerator1;
  hdc::detail::creation_params creation_ps;

  creation_ps.max_nodes_total  = 30;
  creation_ps.max_nodes_support = 20;

  create_nodes( oklut_sim, supports, creation_m, creation_ps );
  CHECK( M.ntk_.num_gates() == 23u );
  M.print_summary();
  hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::ixtsdec;
  hdc::detail::arecovery_params arecovery_ps;
  arecovery_ps.verbose = false;
  M.accuracy_recovery(arecovery_m, arecovery_ps);
  M.print_summary();
}


TEST_CASE( "selection and creation with functions sorting", "[selector]" )
{
  klut_network oklut;
  simulation_view oklut_sim{ oklut };
  std::vector<partial_truth_table> ex;
  partial_truth_table tt(16u);
  for( uint32_t i = 0; i < 4; ++i )
  {
    create_nth_var(tt, i);
    ex.push_back(tt);
  }

  tt = (ex[0]^ex[1])|(ex[2]&ex[3]);
  std::vector<partial_truth_table> targets = {tt};
  model M( oklut_sim, ex, targets );
  std::vector<signal<klut_network>> osignals;

  hdc::detail::selection_method selection_m = hdc::detail::selection_method::depth_selector;
  hdc::detail::selection_params selection_ps;
  selection_ps.max_new_supports = 4;
  selection_ps.max_selection_attempts = 20;
  selection_ps.support_size = 2;
  selection_ps.max_search_depth = std::numeric_limits<uint32_t>::max();

  std::vector<std::vector<signal<klut_network>>> supports = select_variables( oklut_sim, selection_m, selection_ps );

  hdc::detail::creation_method creation_m = hdc::detail::creation_method::ifgenerator1;
  hdc::detail::creation_params creation_ps;

  creation_ps.max_nodes_total  = 10;

  create_nodes( oklut_sim, supports, creation_m, creation_ps );
  CHECK( M.ntk_.num_gates() == 10u );
  M.print_summary();

}

