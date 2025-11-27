#include <catch.hpp>

#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>

#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/resynthesis/aig_syn.hpp>
#include <mockturtle/resynthesis/networks/bindings.hpp>
#include <mockturtle/resynthesis/networks/bound.hpp>
#include <mockturtle/resynthesis/list_simulator.hpp>
#include <mockturtle/resynthesis/power/node_simulator.hpp>
#include <mockturtle/resynthesis/power/analyzers/power_analyzer.hpp>

std::string const test_library = "GATE   inv1    1 O=!a;            PIN * INV 1 999 0.9 0.3 0.9 0.3\n"
                                 "GATE   inv2    2 O=!a;            PIN * INV 2 999 1.0 0.1 1.0 0.1\n"
                                 "GATE   nand2   2 O=!(a*b);        PIN * INV 1 999 1.0 0.2 1.0 0.2\n"
                                 "GATE   and2    3 O=a*b;           PIN * INV 1 999 1.7 0.2 1.7 0.2\n"
                                 "GATE   xor2    4 O=a^b;           PIN * UNKNOWN 2 999 1.9 0.5 1.9 0.5\n"
                                 "GATE   mig3    3 O=a*b+a*c+b*c;   PIN * INV 1 999 2.0 0.2 2.0 0.2\n"
                                 "GATE   xor3    5 O=a^b^c;         PIN * UNKNOWN 2 999 3.0 0.5 3.0 0.5\n"
                                 "GATE   buf     2 O=a;             PIN * NONINV 1 999 1.0 0.0 1.0 0.0\n"
                                 "GATE   zero    0 O=CONST0;\n"
                                 "GATE   one     0 O=CONST1;\n"
                                 "GATE   ha      5 C=a*b;           PIN * INV 1 999 1.7 0.4 1.7 0.4\n"
                                 "GATE   ha      5 S=!a*b+a*!b;     PIN * INV 1 999 2.1 0.4 2.1 0.4\n"
                                 "GATE   fa      6 C=a*b+a*c+b*c;   PIN * INV 1 999 2.1 0.4 2.1 0.4\n"
                                 "GATE   fa      6 S=a^b^c;         PIN * INV 1 999 3.0 0.4 3.0 0.4";

TEST_CASE( "Simulate nodes in mapped networks for power", "[node_simulator]" )
{
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bound_network bound( gates );
  /* network's definition */
  auto const x0 = bound.create_pi();
  auto const x1 = bound.create_pi();

  using signal = typename mockturtle::bound_network::signal;
  signal const f = bound.create_node( std::vector<signal>( { x0, x1 } ), "xor2" )[0];
  bound.create_po( f );

  uint32_t constexpr I = 2u;
  using Tt_t = kitty::static_truth_table<I>;
  std::array<Tt_t, 2> workload0;
  std::array<Tt_t, 2> workload1;
  kitty::create_from_binary_string( workload0[0], "1010" );
  kitty::create_from_binary_string( workload0[1], "0101" );
  kitty::create_from_binary_string( workload1[0], "1100" );
  kitty::create_from_binary_string( workload1[1], "0110" );

  uint32_t constexpr N = 16u;
  using Ntk = mockturtle::bound_network;
  using node_sim_t = mockturtle::node_simulator<Ntk, N, I>;

  node_sim_t sim0( bound.get_node( x0 ), workload0 , 2.5, 0, 5 );
  node_sim_t sim1( bound.get_node( x1 ), workload1 , 4, 1, 6 );
  sim0.print();
  sim1.print();

  double sensing = 0 + 0.5;
  double arrival = 6 + 1.9;
  node_sim_t simf( bound.get_node( f ), sensing, arrival );

  mockturtle::aig_list_simulator<Tt_t> list_simulator;
  std::vector<node_sim_t*> const sims_ptrs { &sim0, &sim1 };
  simf( bound, list_simulator, sims_ptrs );
  simf.print();

}

TEST_CASE( "Compute the power of a network", "[power_analyzer]" )
{
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bound_network bound( gates );
  /* network's definition */
  auto const x0 = bound.create_pi();
  auto const x1 = bound.create_pi();

  using signal = typename mockturtle::bound_network::signal;
  signal const f = bound.create_node( std::vector<signal>( { x0, x1 } ), "xor2" )[0];
  bound.create_po( f );

  uint32_t constexpr I = 2u;
  using Tt_t = kitty::static_truth_table<I>;
  std::array<Tt_t, 2> workload0;
  std::array<Tt_t, 2> workload1;
  kitty::create_from_binary_string( workload0[0], "1010" );
  kitty::create_from_binary_string( workload0[1], "0101" );
  kitty::create_from_binary_string( workload1[0], "1100" );
  kitty::create_from_binary_string( workload1[1], "0110" );

  uint32_t constexpr N = 16u;
  using Ntk = mockturtle::bound_network;
  using node_sim_t = mockturtle::node_simulator<Ntk, N, I>;

  node_sim_t sim0( bound.get_node( x0 ), workload0 , 2.5, 0, 5 );
  node_sim_t sim1( bound.get_node( x1 ), workload1 , 4, 1, 6 );
  sim0.print();
  sim1.print();

  double sensing = 0 + 0.5;
  double arrival = 6 + 1.9;
  node_sim_t simf( bound.get_node( f ), sensing, arrival );

  mockturtle::aig_list_simulator<Tt_t> list_simulator;
  std::vector<node_sim_t*> const sims_ptrs { &sim0, &sim1 };
  simf( bound, list_simulator, sims_ptrs );
  simf.print();
  std::cout << simf.get_dyn_switch() << std::endl;
  std::cout << simf.get_zero_delay() << std::endl;

  mockturtle::power_analyzer<Ntk, N, I> power( bound );
  std::vector<node_sim_t> pis_sims{ sim0, sim1 };
  power( pis_sims );
  std::cout << power.get_dyn_switch() << std::endl;
  std::cout << power.get_zero_delay() << std::endl;
}

TEST_CASE( "A glitchy network", "[power_analyzer]" )
{
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bound_network bound( gates );
  /* network's definition */
  auto const x0 = bound.create_pi();
  auto const x1 = bound.create_pi();
  auto const x2 = bound.create_pi();
  auto const x3 = bound.create_pi();

  using signal = typename mockturtle::bound_network::signal;
  signal const f0 = bound.create_node( std::vector<signal>( { x0, x1 } ), "xor2" )[0];
  signal const f1 = bound.create_node( std::vector<signal>( { f0, x2 } ), "xor2" )[0];
  signal const f2 = bound.create_node( std::vector<signal>( { f1, x3 } ), "xor2" )[0];
  bound.create_po( f2 );

  uint32_t constexpr I = 4u;
  using Tt_t = kitty::static_truth_table<I>;
  std::array<std::array<Tt_t, 2>, 4> workload;
  kitty::create_nth_var( workload[0][0], 0 );
  kitty::create_nth_var( workload[1][0], 1 );
  kitty::create_nth_var( workload[2][0], 2 );
  kitty::create_nth_var( workload[3][0], 3 );
  workload[0][1] = workload[1][0];
  workload[1][1] = workload[2][0];
  workload[2][1] = workload[3][0];
  workload[3][1] = workload[0][0];

  uint32_t constexpr N = 16u;
  using Ntk = mockturtle::bound_network;
  using node_sim_t = mockturtle::node_simulator<Ntk, N, I>;

  node_sim_t sim0( bound.get_node( x0 ), workload[0] , 0.0, -0.1, 0.01 );
  node_sim_t sim1( bound.get_node( x1 ), workload[1] , 0.0, -0.1, 0.01 );
  node_sim_t sim2( bound.get_node( x2 ), workload[2] , 0.0, -0.1, 0.01 );
  node_sim_t sim3( bound.get_node( x3 ), workload[3] , 0.0, -0.1, 0.01 );
  std::cout << "sim0.print()" << std::endl;
  sim0.print();
  std::cout << "sim1.print()" << std::endl;
  sim1.print();
  std::cout << "sim2.print()" << std::endl;
  sim2.print();
  std::cout << "sim3.print()" << std::endl;
  sim3.print();

  mockturtle::aig_list_simulator<Tt_t> list_simulator;
  std::vector<node_sim_t*> const sims_ptrs { &sim0, &sim1, &sim2, &sim3 };

  mockturtle::power_analyzer<Ntk, N, I> power( bound );
  std::vector<node_sim_t> pis_sims{ sim0, sim1, sim2, sim3 };
  power( pis_sims );
  power.print();
}

TEST_CASE( "Another glitchy network", "[power_analyzer]" )
{
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bound_network bound( gates );
  /* network's definition */
  auto const x0 = bound.create_pi();
  auto const x1 = bound.create_pi();
  auto const x2 = bound.create_pi();

  using signal = typename mockturtle::bound_network::signal;
  signal const f0 = bound.create_node( std::vector<signal>( { x0, x1 } ), "xor2" )[0];
  signal const f1 = bound.create_node( std::vector<signal>( { f0, x2 } ), "xor2" )[0];
  signal const f2 = bound.create_node( std::vector<signal>( { f1, f0 } ), "xor2" )[0];
  bound.create_po( f2 );

  uint32_t constexpr I = 4u;
  using Tt_t = kitty::static_truth_table<I>;
  std::array<std::array<Tt_t, 2>, 3> workload;
  kitty::create_nth_var( workload[0][0], 0 );
  kitty::create_nth_var( workload[1][0], 1 );
  kitty::create_nth_var( workload[2][0], 2 );
  workload[0][1] = workload[1][0];
  workload[1][1] = workload[2][0];
  workload[2][1] = workload[3][0];

  uint32_t constexpr N = 16u;
  using Ntk = mockturtle::bound_network;
  using node_sim_t = mockturtle::node_simulator<Ntk, N, I>;

  node_sim_t sim0( bound.get_node( x0 ), workload[0] , -0.05, -0.1, 0.0 );
  node_sim_t sim1( bound.get_node( x1 ), workload[1] , -0.05, -0.1, 0.0 );
  node_sim_t sim2( bound.get_node( x2 ), workload[2] , -0.05, -0.1, 0.0 );
  std::cout << "sim0.print()" << std::endl;
  sim0.print();
  std::cout << "sim1.print()" << std::endl;
  sim1.print();
  std::cout << "sim2.print()" << std::endl;
  sim2.print();

  mockturtle::aig_list_simulator<Tt_t> list_simulator;
  std::vector<node_sim_t*> const sims_ptrs { &sim0, &sim1, &sim2 };

  mockturtle::power_analyzer<Ntk, N, I> power( bound );
  std::vector<node_sim_t> pis_sims{ sim0, sim1, sim2 };
  power( pis_sims );
  power.print();
}