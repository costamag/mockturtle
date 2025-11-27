
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
#include <mockturtle/resynthesis/algorithms/resynthesize.hpp>
#include <mockturtle/resynthesis/analyzers/homo_xxa_analyzer.hpp>
#include <mockturtle/resynthesis/analyzers/homo_xpx_analyzer.hpp>
#include <mockturtle/resynthesis/analyzers/homo_pxx_analyzer.hpp>
#include <mockturtle/resynthesis/power/node_simulator.hpp>
#include <mockturtle/resynthesis/power/analyzers/power_analyzer.hpp>
//#include <mockturtle/resynthesis/power/optimization/power_rewrite.hpp>
#include <mockturtle/resynthesis/database/homo_database.hpp>
#include <mockturtle/resynthesis/database/mapped_list.hpp>
#include <mockturtle/resynthesis/list_simulator.hpp>
#include <mockturtle/resynthesis/networks/bindings.hpp>
#include <mockturtle/resynthesis/networks/bound.hpp>
#include <mockturtle/resynthesis/networks/mapped_to_bound.hpp>
#include <mockturtle/resynthesis/timing/arrival_analyzer.hpp>
#include <mockturtle/utils/tech_library.hpp>


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

int main()
{
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  using Ntk = mockturtle::bound_network;
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

 // mockturtle::arrival_analyzer<Ntk> arrival( bound );
 // arrival.compute_arrival_times();
 // double const delay_bef = arrival.worst_delay();
 // std::cout << "DELAY BEFORE " << delay_bef << std::endl;
  std::cout << bound.area() << std::endl;
  /* get the database */
  /* run resynthesis */
  /* build the analyzer */
  using database_t = mockturtle::homo_database<4u, 1u>;

  database_t db( bound.get_bindings() );
  //db.print();

  mockturtle::resynthesis_params rps;
  rps.dynamic_database = true;
  mockturtle::delay_resynthesize<Ntk, database_t, 4u>( bound, db, rps );

  //db.print();

  std::cout << bound.area() << std::endl;

 // double const delay_aft = arrival.worst_delay();

//  std::cout << "DELAY AFTER " << delay_aft << std::endl;

    return 0;
}