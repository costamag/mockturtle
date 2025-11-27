#include <catch.hpp>
#include <chrono>

#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>

#include <mockturtle/algorithms/emap.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/resynthesis/aig_syn.hpp>
#include <mockturtle/resynthesis/database/homo_database.hpp>
#include <mockturtle/resynthesis/database/mapped_list.hpp>
#include <mockturtle/resynthesis/list_simulator.hpp>
#include <mockturtle/resynthesis/networks/bindings.hpp>
#include <mockturtle/resynthesis/networks/bound.hpp>
#include <mockturtle/resynthesis/networks/mapped_to_bound.hpp>
#include <mockturtle/utils/tech_library.hpp>

std::string const test_library = "GATE   inv1    1 O=!a;            PIN * INV 1 999 0.9 0.3 0.9 0.3\n"
                                 "GATE   inv2    2 O=!a;            PIN * INV 2 999 1.0 0.1 1.0 0.1\n"
                                 "GATE   nand2   2 O=!(a*b);        PIN * INV 1 999 1.0 0.2 1.0 0.2\n"
                                 "GATE   and2    3 O=a*b;           PIN * INV 1 999 1.7 0.2 1.7 0.2\n"
                                 "GATE   or2     3 O=a+b;           PIN * INV 1 999 1.7 0.2 1.7 0.2\n"
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


std::string const mig_library = 
"GATE const0        0.00  z=CONST0;                           \n"                                  
"GATE const1        0.00  z=CONST1;                           \n"                                  
"GATE buf           0.50  Y=A;                                \n"                                  
"    PIN  A  UNKNOWN   1 999    1.00     0.00    1.00     0.00\n"
"GATE inv           0.50  Y=!A;                               \n"                                  
"    PIN  A  UNKNOWN   1 999    1.00     0.00    1.00     0.00\n"
"GATE maj           1.00  Y=(A * B) + (A * C) + (B * C);      \n"                                  
"    PIN  A  UNKNOWN   1 999    3.00     0.00    3.00     0.00\n"
"    PIN  B  UNKNOWN   1 999    1.00     0.00    1.00     0.00\n"
"    PIN  C  UNKNOWN   1 999    2.00     0.00    2.00     0.00";

std::string const sym_library = 
"GATE const0        0.00  z=CONST0;                           \n"                                  
"GATE const1        0.00  z=CONST1;                           \n"                                  
"GATE buf           0.50  Y=A;                                \n"                                  
"    PIN  A  UNKNOWN   1 999    1.00     0.00    1.00     0.00\n"
"GATE AO22x1_ASAP7_75t_R       0.13  Y=(A1 * A2) + (B1 * B2); \n"                                    
" PIN A1  UNKNOWN   1 999    22.74     0.00    22.74     0.00 \n"
" PIN A2  UNKNOWN   1 999    23.52     0.00    23.52     0.00 \n"
" PIN B1  UNKNOWN   1 999    22.75     0.00    22.75     0.00 \n"
" PIN B2  UNKNOWN   1 999    23.43     0.00    23.43     0.00";

TEST_CASE( "Database construction from MIG", "[bound_database]" )
{
  std::vector<mockturtle::gate> gates;

  std::istringstream in( mig_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bindings lib( gates );
  mockturtle::homo_database<6> db( lib );
  db.print();
}

TEST_CASE( "Database construction from SYM", "[bound_database]" )
{
  std::vector<mockturtle::gate> gates;

  std::istringstream in( sym_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bindings lib( gates );
  mockturtle::homo_database<6, 2> db( lib );
  db.print();

  using Ntk = mockturtle::bound_network;
  Ntk ntk( gates );

  /* network's definition */
  auto const a = ntk.create_pi();
  auto const na = ntk.get_node( a );
  auto const b = ntk.create_pi();
  auto const nb = ntk.get_node( b );
  auto const c = ntk.create_pi();
  auto const nc = ntk.get_node( c );
  auto const d = ntk.create_pi();
  auto const nd = ntk.get_node( d );
  auto const z = ntk.get_null();
  auto const f = ntk.create_node( { a, b, c, d }, "AO22x1_ASAP7_75t_R" )[0];
  std::array<typename Ntk::node, 6u> cut { nd, nc, nb, z, z, na };
  db.add_entry( ntk, cut, f );
  db.print();
}

#if 0
TEST_CASE( "Database construction", "[bound_database]" )
{
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bound_network bound( gates );
  /* network's definition */

  auto const bs = bound.get_bindings();
  mockturtle::homo_database<6> db( bs );
  using aig_network = mockturtle::aig_network;
  mockturtle::xag_npn_resynthesis<aig_network, aig_network, mockturtle::xag_npn_db_kind::aig_complete> resyn;
  mockturtle::exact_library_params eps;
  eps.np_classification = false;
  mockturtle::exact_library<aig_network> exact_lib( resyn, eps );

  using Tt4_t = kitty::static_truth_table<4>;
  using tt_hash4 = kitty::hash<Tt4_t>;
  std::unordered_set<Tt4_t, tt_hash4> classes;
  Tt4_t tt;
  do
  {
    const auto res = kitty::exact_p_canonization( tt );
    classes.insert( std::get<0>( res ) );

    auto perm = std::get<2>( res );

    kitty::print_binary( tt );
    std::cout << " ";

    kitty::print_binary( std::get<0>( res ) );
    std::cout << " ";
    for( auto p : perm )
    {
      std::cout << (uint32_t)p << " ";
    }
    std::cout << std::endl;

    kitty::next_inplace( tt );
  } while ( !kitty::is_const0( tt ) );

  using Tt6_t = kitty::static_truth_table<6>;
  Tt6_t tt6;
  mockturtle::klut_network klut;
  for ( auto const& entry : classes )
  {
    kitty::extend_to_inplace( tt6, entry );
    auto supp = kitty::min_base_inplace( tt6 );
  }
}

TEST_CASE( "Database construction for 3 inputs", "[bound_database]" )
{
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bound_network bound( gates );
  /* network's definition */

  auto const bs = bound.get_bindings();
  mockturtle::homo_database<3> db( bs );
  using aig_network = mockturtle::aig_network;
  mockturtle::xag_npn_resynthesis<aig_network, aig_network, mockturtle::xag_npn_db_kind::aig_complete> resyn;
  mockturtle::exact_library_params eps;
  eps.np_classification = false;
  mockturtle::exact_library<aig_network> exact_lib( resyn, eps );

  using Tt_t = kitty::dynamic_truth_table;
  using tt_hash = kitty::hash<Tt_t>;
  std::unordered_set<Tt_t, tt_hash> classes;
  Tt_t tt( 3u );
  do
  {
    const auto res = kitty::exact_p_canonization( tt );
    classes.insert( std::get<0>( res ) );
    auto perm = std::get<2>( res );
    kitty::print_binary( tt );
    std::cout << " ";

    kitty::print_binary( std::get<0>( res ) );
    std::cout << " ";
    for( auto p : perm )
    {
      std::cout << (uint32_t)p << " ";
    }
    std::cout << std::endl;

    kitty::next_inplace( tt );
  } while ( !kitty::is_const0( tt ) );

  using klut_network = mockturtle::klut_network;
  using aig_network = mockturtle::aig_network;
  klut_network klut;
  std::vector<typename mockturtle::klut_network::signal> pis;
  for ( int i = 0; i < 4; ++i )
  {
    pis.push_back( klut.create_pi() );
  }

  std::vector<Tt_t> tts;
  for ( auto const& entry : classes )
  {
    auto const f = klut.create_node( pis, entry );
    klut.create_po( f );
    tts.push_back( entry );
  }
  aig_network aig = mockturtle::convert_klut_to_graph<aig_network, klut_network>( klut );

  mockturtle::tech_library_params tps;
  tps.ignore_symmetries = false; // set to true to drastically speed-up mapping with minor delay increase
  tps.verbose = true;
  mockturtle::tech_library<9> tech_lib( gates, tps );
  mockturtle::emap_params mps;
  mps.matching_mode = mockturtle::emap_params::hybrid;
  mps.area_oriented_mapping = true;
  mps.map_multioutput = false;
  mps.relax_required = 0;
  mockturtle::emap_stats mst;
  mockturtle::cell_view<mockturtle::block_network> mapped = mockturtle::emap<9>( aig, tech_lib, mps, &mst );
  std::cout << mapped.compute_area() << " " << mapped.num_pos() << std::endl;

  int i = 0;
  using node_t = mockturtle::cell_view<mockturtle::block_network>::node;
  std::array<node_t, 3> vars;
  vars[0] = mapped.pi_at( 0 );
  vars[1] = mapped.pi_at( 1 );
  vars[2] = mapped.pi_at( 2 );
  mapped.foreach_po( [&]( auto fo ) {
    kitty::static_truth_table<3> func;
    kitty::create_from_binary_string( func, kitty::to_binary( tts[i++] ) );
    db.add_entry( mapped, vars, fo );
  } );
}

TEST_CASE( "Add one entry and extract it", "[bound_database]" )
{
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bound_network bound( gates );

  std::vector<kitty::static_truth_table<4u>> tts;
  for ( int i = 0; i < 4; ++i )
  {
    tts.emplace_back();
    kitty::create_nth_var( tts[i], i );
  }
 
  /* network's definition */
  auto const x0 = bound.create_pi();
  auto const x1 = bound.create_pi();
  auto const x2 = bound.create_pi();
  auto const n0 = bound.get_node( x0 );
  auto const n1 = bound.get_node( x1 );
  auto const n2 = bound.get_node( x2 );

  auto const func = tts[2] & ( tts[0] | tts[1] );

  using signal = typename mockturtle::bound_network::signal;
  signal const f0 = bound.create_node( std::vector<signal>( { x0, x1 } ), "or2" )[0];
  signal const f1 = bound.create_node( std::vector<signal>( { f0, x2 } ), "and2" )[0];
  bound.create_po( f1 );

  /* network's definition */
  auto const bs = bound.get_bindings();
  mockturtle::homo_database<4> db( bs );
  std::array<mockturtle::bound_network::node, 4> cut{ n0, n1, n2, n0 };
  db.add_entry( bound, cut, f1 );

  auto config = kitty::exact_p_canonization( func );
  auto perm = std::get<2>( config );
  auto res = db.get_entry( func );
  CHECK( res );
  auto perm2 = std::get<0>( *res );
  CHECK( perm == perm2 );
}

TEST_CASE( "Database construction for all 4 input functions", "[bound_database]" )
{
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bound_network bound( gates );
  /* network's definition */

  auto const bs = bound.get_bindings();
  mockturtle::homo_database<4> db( bs );
  using aig_network = mockturtle::aig_network;
  mockturtle::xag_npn_resynthesis<aig_network, aig_network, mockturtle::xag_npn_db_kind::aig_complete> resyn;
  mockturtle::exact_library_params eps;
  eps.np_classification = false;
  mockturtle::exact_library<aig_network> exact_lib( resyn, eps );

  using Tt_t = kitty::dynamic_truth_table;
  using tt_hash = kitty::hash<Tt_t>;
  std::unordered_set<Tt_t, tt_hash> classes;

  using klut_network = mockturtle::klut_network;
  using aig_network = mockturtle::aig_network;
  klut_network klut;
  std::vector<typename mockturtle::klut_network::signal> pis;
  for ( int i = 0; i < 4; ++i )
  {
    pis.push_back( klut.create_pi() );
  }

  std::vector<Tt_t> tts;
  kitty::dynamic_truth_table tt( 4u );
  do
  {
    kitty::next_inplace( tt );
    auto const f = klut.create_node( pis, tt );
    klut.create_po( f );
    tts.push_back( tt );
  } while ( !kitty::is_const0( tt ) );
  aig_network aig = mockturtle::convert_klut_to_graph<aig_network, klut_network>( klut );

  mockturtle::tech_library_params tps;
  tps.ignore_symmetries = false; // set to true to drastically speed-up mapping with minor delay increase
  tps.verbose = true;
  mockturtle::tech_library<9> tech_lib( gates, tps );
  mockturtle::emap_params mps;
  mps.matching_mode = mockturtle::emap_params::hybrid;
  mps.area_oriented_mapping = true;
  mps.map_multioutput = false;
  mps.relax_required = 0;
  mockturtle::emap_stats mst;
  mockturtle::cell_view<mockturtle::block_network> mapped = mockturtle::emap<9>( aig, tech_lib, mps, &mst );
  std::cout << mapped.compute_area() << " " << mapped.num_pos() << std::endl;

  int i = 0;
  using node_t = mockturtle::cell_view<mockturtle::block_network>::node;
  std::array<node_t, 4> vars;
  vars[0] = mapped.pi_at( 0 );
  vars[1] = mapped.pi_at( 1 );
  vars[2] = mapped.pi_at( 2 );
  vars[3] = mapped.pi_at( 3 );
  mapped.foreach_po( [&]( auto fo ) {
    kitty::static_truth_table<4> func;
    kitty::create_from_binary_string( func, kitty::to_binary( tts[i++] ) );
    db.add_entry( mapped, vars, fo );
  } );
}

TEST_CASE( "Database construction one network at a time", "[bound_database]" )
{
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bound_network bound( gates );
  /* network's definition */

  auto const bs = bound.get_bindings();
  mockturtle::homo_database<4> db( bs );
  
  using signal = mockturtle::bound_network::signal;
  using node = mockturtle::bound_network::node;
  std::array<signal, 4u> vars;
  vars[0] = bound.create_pi();
  vars[1] = bound.create_pi();
  vars[2] = bound.create_pi();
  vars[3] = bound.create_pi();

  auto f0 = bound.create_node( std::vector<signal>{ vars[3], vars[1] }, "and2" )[0];
  auto f1 = bound.create_node( std::vector<signal>{ vars[0], vars[2] }, "and2" )[0];
  auto f3 = bound.create_node( std::vector<signal>{ f0, f1 }, "xor2" )[0];
  bound.create_po( f3 );

  std::array<node, 4u> pis;
  pis[0] = bound.get_node( vars[0] );
  pis[1] = bound.get_node( vars[1] );
  pis[2] = bound.get_node( vars[2] );
  pis[3] = bound.get_node( vars[3] );
  db.add_entry( bound, pis, f3 );

}


TEST_CASE( "Database construction with T-canonization", "[bound_database]" )
{
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bound_network bound( gates );
  /* network's definition */

  auto const bs = bound.get_bindings();
  mockturtle::homo_database<6> db( bs );
  
  using signal = mockturtle::bound_network::signal;
  using node = mockturtle::bound_network::node;
  std::array<signal, 6u> vars;
  vars[0] = bound.create_pi();
  vars[1] = bound.create_pi();
  vars[2] = bound.create_pi();
  vars[3] = bound.create_pi();
  vars[4] = bound.create_pi();
  vars[5] = bound.create_pi();

  auto f0 = bound.create_node( std::vector<signal>{ vars[5], vars[4] }, "and2" )[0];
  auto f1 = bound.create_node( std::vector<signal>{ f0, vars[3] }, "and2" )[0];
  auto f2 = bound.create_node( std::vector<signal>{ f1, vars[2] }, "and2" )[0];
  auto f3 = bound.create_node( std::vector<signal>{ vars[1] }, "inv1" )[0];
  auto f4 = bound.create_node( std::vector<signal>{ vars[0], f3 }, "nand2" )[0];
  auto f5 = bound.create_node( std::vector<signal>{ f2, f4 }, "xor2" )[0];
  bound.create_po( f5 );

  std::array<node, 6u> pis;
  pis[0] = bound.get_node( vars[0] );
  pis[1] = bound.get_node( vars[1] );
  pis[2] = bound.get_node( vars[2] );
  pis[3] = bound.get_node( vars[3] );
  pis[4] = bound.get_node( vars[4] );
  pis[5] = bound.get_node( vars[5] );
  db.add_entry( bound, pis, f5 );

}
#endif