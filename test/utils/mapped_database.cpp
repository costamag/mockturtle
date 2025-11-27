#include <catch.hpp>

#include <cstdint>
#include <vector>

#include <lorina/genlib.hpp>
#include <lorina/super.hpp>
#include <mockturtle/algorithms/emap.hpp>
#include <mockturtle/generators/arithmetic.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/super_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/block.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mapped.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/utils/mapped_database.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/views/cell_view.hpp>
#include <mockturtle/views/dont_touch_view.hpp>

using namespace mockturtle;

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

TEST_CASE( "Mapped database read existing", "[mapped_database]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<3> lib( gates );


  mapped_network ntk( gates );
  std::vector<mapped_network::signal> pis( 4 );
  std::generate( std::begin( pis ), std::end( pis ),
                 [&]() { return ntk.create_pi(); } );
  /* constants */

  kitty::dynamic_truth_table tt( 2 );
  std::vector<kitty::dynamic_truth_table> tts;
  do
  {
    tts.push_back( tt );
    kitty::next_inplace( tt );
  } while ( !kitty::is_const0( tt ) );

  mapped_database_params ps;
  ps.num_vars = 2;
  mapped_database db( "../experiments/databases/mapped/test2.db", lib, ps );

  CHECK( tts.size() == 16 );
  
  std::vector<mapped_network::signal> inputs = { pis[2], pis[3] };

  /* 0000 */
  auto matches = db.get_matches( tts[0], inputs );
  CHECK( (*matches).first == std::vector<mapped_network::signal>{ pis[2], pis[3] } );

  /* 0001 */
  matches = db.get_matches( tts[1], inputs );
  CHECK( (*matches).first == std::vector<mapped_network::signal>{ pis[2], pis[3] } );

  /* 0010 */
  matches = db.get_matches( tts[2], inputs );
  CHECK( (*matches).first == std::vector<mapped_network::signal>{ pis[2], pis[3] } );

  /* 0011 */
  matches = db.get_matches( tts[3], inputs );
  CHECK( (*matches).first == std::vector<mapped_network::signal>{ pis[2], pis[3] } );

  /* 0100 */
  matches = db.get_matches( tts[4], inputs );
  CHECK( (*matches).first == std::vector<mapped_network::signal>{ pis[3], pis[2] } );

  /* 0101 */
  matches = db.get_matches( tts[5], inputs );
  CHECK( (*matches).first == std::vector<mapped_network::signal>{ pis[3], pis[2] } );

  /* 0110 */
  matches = db.get_matches( tts[6], inputs );
  CHECK( (*matches).first == std::vector<mapped_network::signal>{ pis[2], pis[3] } );

  /* 0111 */
  matches = db.get_matches( tts[7], inputs );
  CHECK( (*matches).first == std::vector<mapped_network::signal>{ pis[2], pis[3] } );

  /* 1000 */
  matches = db.get_matches( tts[8], inputs );
  CHECK( (*matches).first == std::vector<mapped_network::signal>{ pis[2], pis[3] } );

  /* 1001 */
  matches = db.get_matches( tts[9], inputs );
  CHECK( (*matches).first == std::vector<mapped_network::signal>{ pis[2], pis[3] } );

  /* 1010 */
  matches = db.get_matches( tts[10], inputs );
  CHECK( (*matches).first == std::vector<mapped_network::signal>{ pis[2], pis[3] } );

  /* 1011 */
  matches = db.get_matches( tts[11], inputs );
  CHECK( (*matches).first == std::vector<mapped_network::signal>{ pis[2], pis[3] } );

  /* 1100 */
  matches = db.get_matches( tts[12], inputs );
  CHECK( (*matches).first == std::vector<mapped_network::signal>{ pis[3], pis[2] } );

  /* 1101 */
  matches = db.get_matches( tts[13], inputs );
  CHECK( (*matches).first == std::vector<mapped_network::signal>{ pis[3], pis[2] } );

  /* 1110 */
  matches = db.get_matches( tts[14], inputs );
  CHECK( (*matches).first == std::vector<mapped_network::signal>{ pis[2], pis[3] } );

  /* 1111 */
  matches = db.get_matches( tts[15], inputs );
  CHECK( (*matches).first == std::vector<mapped_network::signal>{ pis[2], pis[3] } );

}



//1111 2 1 1 
//1110 2 1 1 2 0 1 3 0 2 4 5 2 6 
//1010 2 1 1 2 7 4 
//1001 2 1 1 2 0 2 4 3 4 5 
//1000 2 1 2 2 3 3 4 
//1011 2 1 1 2 0 2 4 3 2 5 
//0111 2 1 2 2 3 2 4 
//0110 2 1 2 2 3 4 4 
//0011 2 1 1 3 0 4 
//0010 2 1 1 3 0 2 4 2 3 5 
//0001 2 1 1 2 0 1 3 0 2 4 5 3 6 
//0000 2 1 0 