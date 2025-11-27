#include <catch.hpp>

#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>

#include <mockturtle/algorithms/emap.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/resynthesis/algorithms/resynthesize.hpp>
#include <mockturtle/resynthesis/power/node_simulator.hpp>
#include <mockturtle/resynthesis/database/homo_database.hpp>
#include <mockturtle/resynthesis/networks/bindings.hpp>
#include <mockturtle/resynthesis/networks/bound.hpp>
#include <mockturtle/utils/tech_library.hpp>


std::string const test_library = "GATE   inv1    1.0 O=!a ;         PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   and2    1.0 O=a*b;         PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   xor2    1.0 O=a^b;         PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   or2     1.0 O=a+b;         PIN * INV 1   999 1.0 0.0 1.0 0.0";

TEST_CASE( "Placement resynthesis - AND tree", "[placement_resynthesis]" )
{
  std::vector<mockturtle::gate> gates;
  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );
  using Ntk = mockturtle::bound_network;
  Ntk ntk( gates );
  mockturtle::bindings binds( gates );

  using signal = typename Ntk::signal;
  std::vector<signal> fs;
  /* Database construction from network */

  uint32_t constexpr num_vars = 4;
  uint32_t constexpr num_steps = 10;
  uint32_t constexpr cube_size = 4;
  uint32_t constexpr max_leaves = 10;
  uint32_t constexpr capacity = 10;
  using database_t = mockturtle::homo_database<num_vars, capacity>;
  database_t database( binds );
  
  Ntk db( gates );
  auto const db_a = db.create_pi();
  auto const db_na = db.get_node( db_a );
  auto const db_b = db.create_pi();
  auto const db_nb = db.get_node( db_b );
  auto const db_c = db.create_pi();
  auto const db_nc = db.get_node( db_c );
  auto const db_d = db.create_pi();
  auto const db_nd = db.get_node( db_d );

  auto const s1 = db.create_node( { db_a, db_b }, "and2" )[0];
  auto const s2 = db.create_node( { db_c, db_d }, "and2" )[0];
  auto const s3 = db.create_node( { s1, s2 }, "and2" )[0];
  std::array<typename Ntk::node, 4u> cut { db_na, db_nb, db_nc, db_nd };
  database.add_entry( db, cut, s3 );
  database.print();


  signal const a = ntk.create_pi();
  signal const b = ntk.create_pi();
  signal const c = ntk.create_pi();
  signal const d = ntk.create_pi();
  //signal const e = ntk.create_pi();
  //signal const f = ntk.create_pi();
  //signal const g = ntk.create_pi();
  //signal const h = ntk.create_pi();

  signal f1 = ntk.create_node( { a, b }, "and2" )[0];
  signal f2 = ntk.create_node( { f1, c }, "and2" )[0];
  signal f3 = ntk.create_node( { f2, d }, "and2" )[0];
  //signal f4 = ntk.create_node( { f3, e }, "and2" )[0];
  //signal f5 = ntk.create_node( { f4, f }, "and2" )[0];
  //signal f6 = ntk.create_node( { f5, g }, "and2" )[0];
  //signal f7 = ntk.create_node( { f6, h }, "and2" )[0];
  ntk.create_po( f3 );

  mockturtle::placement_params pps;
  pps.num_placers = 5;
  mockturtle::placement_analyzer<mockturtle::bound_network, 20, 20> analyzer( ntk, pps );
  analyzer.save_result( "/home/costamag/POST_MAPPING/mockturtle/build/graphs/test0_i0" );
  double cost_bef = analyzer.compute_cost();
  std::cout << "cost before " << cost_bef << std::endl;
  analyzer.fire();
  double cost_aft = analyzer.compute_cost();
  std::cout << "cost after  " << cost_aft << std::endl;
  analyzer.save_result( "/home/costamag/POST_MAPPING/mockturtle/build/graphs/test0_i1" );

  mockturtle::resynthesis_params rps;
  rps.preserve_depth = true;
  rps.dynamic_database = false;
  rps.try_struct = true;
  rps.try_window = false;
  rps.try_simula = false;
  rps.try_rewire = false;

  analyzer.print_state();
  mockturtle::place_resynthesize<Ntk, database_t, 4, 2, 4, 4, 20, 20>( ntk, database, rps );

  mockturtle::placement_analyzer<mockturtle::bound_network, 20, 20> analyzer2( ntk, pps );
  analyzer2.save_result( "/home/costamag/POST_MAPPING/mockturtle/build/graphs/test1_e0" );
  double cost_opt_0 = analyzer2.compute_cost();
  std::cout << "cost opt 0 " << cost_opt_0 << std::endl;
  analyzer2.fire();
  double cost_opt_1 = analyzer2.compute_cost();
  std::cout << "cost after  " << cost_opt_1 << std::endl;
  analyzer2.save_result( "/home/costamag/POST_MAPPING/mockturtle/build/graphs/test1_e1" );
}

TEST_CASE( "Placement resynthesis - large AND tree", "[placement_resynthesis]" )
{
  std::vector<mockturtle::gate> gates;
  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );
  using Ntk = mockturtle::bound_network;
  Ntk ntk( gates );
  mockturtle::bindings binds( gates );

  using signal = typename Ntk::signal;
  std::vector<signal> fs;
  /* Database construction from network */

  uint32_t constexpr num_vars = 4;
  uint32_t constexpr num_steps = 10;
  uint32_t constexpr cube_size = 4;
  uint32_t constexpr max_leaves = 10;
  uint32_t constexpr capacity = 10;
  using database_t = mockturtle::homo_database<num_vars, capacity>;
  database_t database( binds );
  
  Ntk db( gates );
  auto const db_a = db.create_pi();
  auto const db_na = db.get_node( db_a );
  auto const db_b = db.create_pi();
  auto const db_nb = db.get_node( db_b );
  auto const db_c = db.create_pi();
  auto const db_nc = db.get_node( db_c );
  auto const db_d = db.create_pi();
  auto const db_nd = db.get_node( db_d );

  auto const s1 = db.create_node( { db_a, db_b }, "and2" )[0];
  auto const s2 = db.create_node( { db_c, db_d }, "and2" )[0];
  auto const s3 = db.create_node( { s1, s2 }, "and2" )[0];
  std::array<typename Ntk::node, 4u> cut { db_na, db_nb, db_nc, db_nd };
  database.add_entry( db, cut, s3 );
  database.print();


  signal const a = ntk.create_pi();
  signal const b = ntk.create_pi();
  signal const c = ntk.create_pi();
  signal const d = ntk.create_pi();
  signal const e = ntk.create_pi();
  signal const f = ntk.create_pi();
  signal const g = ntk.create_pi();
  signal const h = ntk.create_pi();

  signal f1 = ntk.create_node( { a, b }, "and2" )[0];
  signal f2 = ntk.create_node( { f1, c }, "and2" )[0];
  signal f3 = ntk.create_node( { f2, d }, "and2" )[0];
  signal f4 = ntk.create_node( { f3, e }, "and2" )[0];
  signal f5 = ntk.create_node( { f4, f }, "and2" )[0];
  signal f6 = ntk.create_node( { f5, g }, "and2" )[0];
  signal f7 = ntk.create_node( { f6, h }, "and2" )[0];
  ntk.create_po( f7 );

  mockturtle::placement_params pps;
  pps.num_placers = 1;
  mockturtle::placement_analyzer<mockturtle::bound_network, 20, 20> analyzer( ntk, pps );
  analyzer.save_result( "/home/costamag/POST_MAPPING/mockturtle/build/graphs/test0_i0" );
  double cost_bef = analyzer.compute_cost();
  std::cout << "cost before " << cost_bef << std::endl;
  analyzer.fire();
  double cost_aft = analyzer.compute_cost();
  std::cout << "cost after  " << cost_aft << std::endl;
  analyzer.save_result( "/home/costamag/POST_MAPPING/mockturtle/build/graphs/test0_i1" );

  mockturtle::resynthesis_params rps;
  rps.preserve_depth = true;
  rps.dynamic_database = false;
  rps.try_struct = true;
  rps.try_window = false;
  rps.try_simula = false;
  rps.try_rewire = false;

  analyzer.print_state();
  mockturtle::place_resynthesize<Ntk, database_t, 4, 2, 4, 4, 20, 20>( ntk, database, rps );

  mockturtle::placement_analyzer<mockturtle::bound_network, 20, 20> analyzer2( ntk, pps );
  analyzer2.save_result( "/home/costamag/POST_MAPPING/mockturtle/build/graphs/test1_e0" );
  double cost_opt_0 = analyzer2.compute_cost();
  std::cout << "cost opt 0 " << cost_opt_0 << std::endl;
  analyzer2.fire();
  double cost_opt_1 = analyzer2.compute_cost();
  std::cout << "cost after  " << cost_opt_1 << std::endl;
  analyzer2.save_result( "/home/costamag/POST_MAPPING/mockturtle/build/graphs/test1_e1" );
}


#if 0
TEST_CASE( "Placement resynthesis - rewire EXOR", "[placement_resynthesis]" )
{
  std::vector<mockturtle::gate> gates;
  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );
  using Ntk = mockturtle::bound_network;
  Ntk ntk( gates );
  mockturtle::bindings binds( gates );

  using signal = typename Ntk::signal;
  std::vector<signal> fs;
  /* Database construction from network */

  uint32_t constexpr num_vars = 4;
  uint32_t constexpr num_steps = 10;
  uint32_t constexpr cube_size = 4;
  uint32_t constexpr max_leaves = 10;
  uint32_t constexpr capacity = 10;
  using database_t = mockturtle::homo_database<num_vars, capacity>;
  database_t database( binds );
  

  signal const a = ntk.create_pi();
  signal const b = ntk.create_pi();
  signal const c = ntk.create_pi();
  signal const d = ntk.create_pi();
  signal const e = ntk.create_pi();

  signal f1 = ntk.create_node( { a, b }, "xor2" )[0];
  signal f2 = ntk.create_node( { a }, "inv1" )[0];
  signal f3 = ntk.create_node( { b }, "inv1" )[0];
  signal f4 = ntk.create_node( { a, f3 }, "and2" )[0];
  signal f5 = ntk.create_node( { b, f2 }, "and2" )[0];
  signal f6 = ntk.create_node( { f4, f5 }, "or2" )[0];
  signal f7 = ntk.create_node( { f6, c }, "or2" )[0];
  signal f7n = ntk.create_node( { f7 }, "inv1" )[0];
  signal f8 = ntk.create_node( { c, d }, "and2" )[0];
  signal f9 = ntk.create_node( { c, e }, "and2" )[0];
  signal f10 = ntk.create_node( { f8, f9 }, "xor2" )[0];
  ntk.create_po( f1 );
  ntk.create_po( f7n );
  ntk.create_po( f10 );
  ntk.create_po( f10 );
  ntk.create_po( f10 );
  
  mockturtle::placement_params pps;
  pps.area_density = 0.1;
  pps.num_placers = 1;
  mockturtle::placement_analyzer<mockturtle::bound_network, 20, 20> analyzer( ntk, pps );
  analyzer.save_result( "/home/costamag/POST_MAPPING/mockturtle/build/graphs/test1_bef" );
  double cost_bef = analyzer.compute_cost();
  std::cout << "cost before " << cost_bef << std::endl;
  analyzer.fire();
  double cost_aft = analyzer.compute_cost();
  std::cout << "cost after  " << cost_aft << std::endl;
  analyzer.save_result( "/home/costamag/POST_MAPPING/mockturtle/build/graphs/test1_aft" );

  mockturtle::resynthesis_params rps;
  rps.preserve_depth = true;
  rps.dynamic_database = false;
  rps.try_struct = true;
  rps.try_window = true;
  rps.try_simula = true;
  rps.try_rewire = true;

  mockturtle::place_resynthesize<Ntk, database_t, 4, 2, 4, 4, 20, 20>( ntk, database, rps );

  mockturtle::placement_analyzer<mockturtle::bound_network, 20, 20> analyzer2( ntk, pps );
  analyzer2.save_result( "/home/costamag/POST_MAPPING/mockturtle/build/graphs/test1_opt_0" );
  double cost_opt_0 = analyzer2.compute_cost();
  std::cout << "cost opt 0 " << cost_opt_0 << std::endl;
  analyzer2.fire();
  double cost_opt_1 = analyzer2.compute_cost();
  std::cout << "cost after  " << cost_opt_1 << std::endl;
  analyzer2.save_result( "/home/costamag/POST_MAPPING/mockturtle/build/graphs/test1_opt_1" );  

}


TEST_CASE( "Power resynthesis inverter chain", "[power_resynthesis]" )
{
  std::vector<mockturtle::gate> gates;
  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );
  using Ntk = mockturtle::bound_network;
  Ntk ntk( gates );
  mockturtle::bindings binds( gates );

  using signal = typename Ntk::signal;
  std::vector<signal> fs;
  /* Database construction from network */
  signal const a = ntk.create_pi();

  uint32_t constexpr num_vars = 4;
  uint32_t constexpr num_steps = 10;
  uint32_t constexpr cube_size = 4;
  uint32_t constexpr max_leaves = 10;
  uint32_t constexpr capacity = 10;
  mockturtle::choice_order constexpr order = mockturtle::power_area_delay;
  using database_t = mockturtle::homo_database<num_vars, capacity, order>;
  database_t database( binds );
  
  mockturtle::mapped_list list( num_vars );
  list.create_po( num_vars/2 );
  database.add_entry( list );
  database.print();

  using node_simulator_t = mockturtle::node_simulator<Ntk, num_steps, cube_size>;

  /* Create the input simulations */
  std::vector<node_simulator_t> pis_sims;
  std::default_random_engine::result_type seed1 = 2;
  std::default_random_engine::result_type seed2 = 3;
  std::array<kitty::static_truth_table<cube_size>, 2> workload;
  int j = 1;
  ntk.foreach_pi([&]( auto n ){
    kitty::create_random( workload[0], seed1*j );
    kitty::create_random( workload[1], seed2*j );
    j += 2;
    node_simulator_t sim( n, workload, -0.1, 0 );
    pis_sims.push_back( sim );
  });

  for ( int i = 0; i < 4; ++i )
  {
    fs.clear();
    fs.push_back( ntk.create_node( { a }, "inv1" )[0] );
    fs.push_back( ntk.create_node( { fs[0] }, "inv1" )[0] );
    ntk.create_po( fs.back() );
    CHECK( ( ntk.fanout_size( ntk.get_node( fs.back() ) ) == 1 ) ); 
    CHECK( ( ntk.fanout_size( ntk.get_node( fs.back() ) ) == 1 ) ); 
  
  
    mockturtle::power_analyzer<Ntk, num_steps, cube_size> power_analyzer( ntk );
    power_analyzer( pis_sims );
  
    double const zero_delay0 = power_analyzer.get_zero_delay();
    double const dyn_switch0 = power_analyzer.get_dyn_switch();
    double const dyn_power0 = power_analyzer.get_dyn_power();
    double const glitching0 = dyn_switch0 - zero_delay0;
    mockturtle::resynthesis_params rps;
    rps.collect_mode = false;
    rps.try_rewire = i == 0;
    rps.try_struct = i == 1;
    rps.try_window = i == 2;
    rps.try_simula = i == 3;
    ntk.print();
    power_analyzer.print_power();
    rps.analyzer_ps.preserve_interval = true;
    rps.analyzer_ps.preserve_depth = false;
    mockturtle::glitch_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves>( ntk, database, rps );
    using node_simulator_t = mockturtle::node_simulator<Ntk, num_steps, cube_size>;
    ntk.print();
  
    power_analyzer( pis_sims );
    double const zero_delay1 = power_analyzer.get_zero_delay();
    double const dyn_switch1 = power_analyzer.get_dyn_switch();
    double const dyn_power1 = power_analyzer.get_dyn_power();
    double const glitching1 = dyn_switch1 - zero_delay1;
  
    power_analyzer.print_power();
    CHECK( zero_delay1 < zero_delay0 );
    CHECK( dyn_switch1 < dyn_switch0 );
    CHECK( dyn_power1 < dyn_power0 );
    CHECK( glitching1 == glitching0 );
  }

}

TEST_CASE( "Power resynthesis suboptimal EXOR 2", "[power_resynthesis]" )
{
  std::vector<mockturtle::gate> gates;
  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );
  using Ntk = mockturtle::bound_network;
  mockturtle::bindings binds( gates );
  
  using signal = typename Ntk::signal;
  std::vector<signal> fs;
  std::vector<signal> xs;
  /* Database construction from network */


  uint32_t constexpr num_vars = 4;
  uint32_t constexpr num_steps = 10;
  uint32_t constexpr cube_size = 4;
  uint32_t constexpr max_leaves = 10;
  uint32_t constexpr capacity = 10;
  mockturtle::choice_order constexpr order = mockturtle::power_area_delay;
  using database_t = mockturtle::homo_database<num_vars, capacity, order>;
  database_t database( binds );
  
  mockturtle::mapped_list list( num_vars );
  list.create_po( num_vars/2 );
  database.add_entry( list );
  database.print();

  using node_simulator_t = mockturtle::node_simulator<Ntk, num_steps, cube_size>;

  
  for ( int i = 1; i < 4; ++i )
  {
    fs.clear();
    xs.clear();
    Ntk ntk( gates );
    for ( int i = 0; i < 2; ++i )
    xs.push_back( std::move( ntk.create_pi() ) );
    /* Create the input simulations */
    std::vector<node_simulator_t> pis_sims;
    std::default_random_engine::result_type seed1 = 2;
    std::default_random_engine::result_type seed2 = 3;
    std::array<kitty::static_truth_table<cube_size>, 2> workload;
    int j = 1;
    ntk.foreach_pi([&]( auto n ){
      kitty::create_random( workload[0], seed1*j );
      kitty::create_random( workload[1], seed2*j );
      j += 2;
      node_simulator_t sim( n, workload, -0.1, 0 );
      pis_sims.push_back( sim );
    });
    fs.push_back( std::move( ntk.create_node( { xs[0] }, "inv1" )[0] ) );
    fs.push_back( std::move( ntk.create_node( { fs.back(), xs[1] }, "xor2" )[0] ) );
    fs.push_back( std::move( ntk.create_node( { fs.back() }, "inv1" )[0] ) );
    
    ntk.create_po( fs.back() );
    
    mockturtle::power_analyzer<Ntk, num_steps, cube_size> power_analyzer( ntk );
    power_analyzer( pis_sims );
    
    double const zero_delay0 = power_analyzer.get_zero_delay();
    double const dyn_switch0 = power_analyzer.get_dyn_switch();
    double const dyn_power0 = power_analyzer.get_dyn_power();
    double const glitching0 = dyn_switch0 - zero_delay0;
    mockturtle::resynthesis_params rps;
    rps.analyzer_ps.preserve_interval = true;
    rps.analyzer_ps.preserve_depth = false;
    rps.collect_mode = false;
    rps.try_rewire = false;
    rps.try_struct = i == 1;
    rps.try_window = i == 2;
    rps.try_simula = i == 3;
    ntk.print();
    power_analyzer.print_power();
    mockturtle::glitch_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves>( ntk, database, rps );
    using node_simulator_t = mockturtle::node_simulator<Ntk, num_steps, cube_size>;
    ntk.print();
  
    power_analyzer( pis_sims );
    double const zero_delay1 = power_analyzer.get_zero_delay();
    double const dyn_switch1 = power_analyzer.get_dyn_switch();
    double const dyn_power1 = power_analyzer.get_dyn_power();
    double const glitching1 = dyn_switch1 - zero_delay1;
  
    power_analyzer.print_power();
    CHECK( zero_delay1 < zero_delay0 );
    CHECK( dyn_switch1 < dyn_switch0 );
    CHECK( dyn_power1 < dyn_power0 );
    CHECK( glitching1 < glitching0 );
  }
}

TEST_CASE( "Power resynthesis EXOR2 to OR2", "[power_resynthesis]" )
{
  std::vector<mockturtle::gate> gates;
  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );
  using Ntk = mockturtle::bound_network;
  mockturtle::bindings binds( gates );
  
  using signal = typename Ntk::signal;
  std::vector<signal> fs;
  std::vector<signal> xs;
  /* Database construction from network */

  uint32_t constexpr num_vars = 4;
  uint32_t constexpr num_steps = 10;
  uint32_t constexpr cube_size = 4;
  uint32_t constexpr max_leaves = 8;
  uint32_t constexpr capacity = 10;
  mockturtle::choice_order constexpr order = mockturtle::power_area_delay;
  using database_t = mockturtle::homo_database<num_vars, capacity, order>;
  database_t database( binds );
  
  mockturtle::mapped_list list( num_vars );
  auto l4 = list.create_node( { 2, 3 }, 1 );
  auto l5 = list.create_node( { 2, 3 }, 1 );
  auto l6 = list.create_node( { l4, l5 }, 2 );
  list.create_po( l6 );
  database.add_entry( list );
  database.print();

  using node_simulator_t = mockturtle::node_simulator<Ntk, num_steps, cube_size>;

  
  for ( int i = 0; i < 4; ++i )
  {
    fs.clear();
    xs.clear();
    Ntk ntk( gates );
    for ( int i = 0; i < 2; ++i )
      xs.push_back( std::move( ntk.create_pi() ) );
    /* Create the input simulations */
    std::vector<node_simulator_t> pis_sims;
    std::default_random_engine::result_type seed1 = 2;
    std::default_random_engine::result_type seed2 = 3;
    std::array<kitty::static_truth_table<cube_size>, 2> workload;
    int j = 1;
    ntk.foreach_pi([&]( auto n ){
      kitty::create_random( workload[0], seed1*j );
      kitty::create_random( workload[1], seed2*j );
      j += 2;
      node_simulator_t sim( n, workload, -0.1, 0 );
      pis_sims.push_back( sim );
    });
    fs.push_back( std::move( ntk.create_node( { xs[0], xs[1] }, "and2" )[0] ) );
    fs.push_back( std::move( ntk.create_node( { xs[0], xs[1] }, "nor2" )[0] ) );
    fs.push_back( std::move( ntk.create_node( { xs[0] }, "inv1" )[0] ) );
    fs.push_back( std::move( ntk.create_node( { fs[2], xs[1] }, "xor2" )[0] ) );
    fs.push_back( std::move( ntk.create_node( { fs[0], fs[3] }, "or2" )[0] ) );
  
    ntk.create_po( fs[1] );
    ntk.create_po( fs[4] );
    
    mockturtle::power_analyzer<Ntk, num_steps, cube_size> power_analyzer( ntk );
    power_analyzer( pis_sims );
    
    double const zero_delay0 = power_analyzer.get_zero_delay();
    double const dyn_switch0 = power_analyzer.get_dyn_switch();
    double const dyn_power0 = power_analyzer.get_dyn_power();
    double const glitching0 = dyn_switch0 - zero_delay0;
    mockturtle::resynthesis_params rps;
    rps.analyzer_ps.preserve_interval = true;
    rps.analyzer_ps.preserve_depth = false;
    rps.collect_mode = false;
    rps.try_rewire = i == 0;
    rps.try_struct = i == 1;
    rps.try_window = i == 2;
    rps.try_simula = i == 3;
    ntk.print();
    power_analyzer.print_power();
    mockturtle::glitch_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves>( ntk, database, rps );
    using node_simulator_t = mockturtle::node_simulator<Ntk, num_steps, cube_size>;
    ntk.print();
  
    power_analyzer( pis_sims );
    double const zero_delay1 = power_analyzer.get_zero_delay();
    double const dyn_switch1 = power_analyzer.get_dyn_switch();
    double const dyn_power1 = power_analyzer.get_dyn_power();
    double const glitching1 = dyn_switch1 - zero_delay1;
  
    power_analyzer.print_power();
    CHECK( zero_delay1 < zero_delay0 );
    CHECK( dyn_switch1 < dyn_switch0 );
    CHECK( dyn_power1 < dyn_power0 );
    CHECK( glitching1 < glitching0 );
  }
}


TEST_CASE( "EXOR tree to EXOR tree", "[power_resynthesis]" )
{
  std::vector<mockturtle::gate> gates;
  std::istringstream in( asap7_subset );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );
  using Ntk = mockturtle::bound_network;
  mockturtle::bindings binds( gates );
  
  using signal = typename Ntk::signal;
  std::vector<signal> fs;
  std::vector<signal> xs;
  /* Database construction from network */

  uint32_t constexpr num_vars = 4;
  uint32_t constexpr num_steps = 20;
  uint32_t constexpr cube_size = 12;
  uint32_t constexpr max_leaves = 8;
  uint32_t constexpr capacity = 10;
  mockturtle::choice_order constexpr order = mockturtle::power_area_delay;
  using database_t = mockturtle::homo_database<num_vars, capacity, order>;
  database_t database( binds );
  
  mockturtle::mapped_list list( num_vars );
  auto l4 = list.create_node( { 2, 3 }, 1 );
  auto l5 = list.create_node( { 4, 5 }, 1 );
  auto l6 = list.create_node( { l4, l5 }, 2 );
  list.create_po( l6 );
  database.add_entry( list );
  database.print();

  using node_simulator_t = mockturtle::node_simulator<Ntk, num_steps, cube_size>;

  
  fs.clear();
  xs.clear();
  Ntk ntk( gates );
  for ( int i = 0; i < 4; ++i )
    xs.push_back( std::move( ntk.create_pi() ) );
  /* Create the input simulations */
  std::vector<node_simulator_t> pis_sims;
  std::default_random_engine::result_type seed = 1;
  std::array<kitty::static_truth_table<cube_size>, 2> workload;
  int j = 1;
  ntk.foreach_pi([&]( auto n ){
    kitty::create_random( workload[0], seed*j );
    kitty::create_random( workload[1], seed*( j + 1 ) );
    j += 2;
    node_simulator_t sim( n, workload, -1, 0 );
    pis_sims.push_back( sim );
  });
  fs.push_back( std::move( ntk.create_node( { xs[0], xs[1] }, "xor2" )[0] ) );
  fs.push_back( std::move( ntk.create_node( { xs[2], xs[3] }, "xor2" )[0] ) );
  fs.push_back( std::move( ntk.create_node( { fs[0], fs[1] }, "xor2" )[0] ) );

  ntk.create_po( fs[2] );
  
  mockturtle::power_analyzer<Ntk, num_steps, cube_size> power_analyzer( ntk );
  power_analyzer( pis_sims );
  
  double const zero_delay0 = power_analyzer.get_zero_delay();
  double const dyn_switch0 = power_analyzer.get_dyn_switch();
  double const dyn_power0 = power_analyzer.get_dyn_power();
  double const glitching0 = dyn_switch0 - zero_delay0;
  mockturtle::resynthesis_params rps;
  rps.collect_mode = false;
  rps.try_rewire = false;
  rps.try_struct = false;
  rps.try_window = true;
  rps.try_simula = false;
  rps.analyzer_ps.preserve_depth = false;
  rps.analyzer_ps.use_fanouts = true;
  ntk.print();
  power_analyzer.print_power();
  mockturtle::glitch_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves>( ntk, database, rps );
  using node_simulator_t = mockturtle::node_simulator<Ntk, num_steps, cube_size>;
  ntk.print();

  power_analyzer( pis_sims );
  double const zero_delay1 = power_analyzer.get_zero_delay();
  double const dyn_switch1 = power_analyzer.get_dyn_switch();
  double const dyn_power1 = power_analyzer.get_dyn_power();
  double const glitching1 = dyn_switch1 - zero_delay1;

  power_analyzer.print_power();
  CHECK( zero_delay1 < zero_delay0 );
  CHECK( dyn_switch1 < dyn_switch0 );
  CHECK( dyn_power1 < dyn_power0 );
  CHECK( glitching1 < glitching0 );
}
#endif