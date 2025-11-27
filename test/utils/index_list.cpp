#include <catch.hpp>

#include <kitty/static_truth_table.hpp>
#include <lorina/genlib.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/mapped.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/networks/muxig.hpp>
#include <mockturtle/utils/index_list.hpp>

using namespace mockturtle;

std::string const test_library = "GATE   inv1    1 O=!a;            PIN * INV 1 999 0.9 0.3 0.9 0.3\n"
                                 "GATE   inv2    2 O=!a;            PIN * INV 2 999 1.0 0.1 1.0 0.1\n"
                                 "GATE   nand2   2 O=!(a*b);        PIN * INV 1 999 1.0 0.2 1.0 0.2\n"
                                 "GATE   and2    3 O=a*b;           PIN * INV 1 999 1.7 0.2 1.7 0.2\n"
                                 "GATE   nor2    4 O=!(a+b);        PIN * INV 1 999 1.7 0.2 1.7 0.2\n"
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

TEST_CASE( "decode mig_index_list into mig_network", "[index_list]" )
{
  std::vector<uint32_t> const raw_list{ 4 | ( 1 << 8 ) | ( 2 << 16 ), 2, 4, 6, 10, 4, 8, 12 };
  mig_index_list mig_il{ raw_list };

  mig_network mig;
  decode( mig, mig_il );

  CHECK( mig.num_gates() == 2u );
  CHECK( mig.num_pis() == 4u );
  CHECK( mig.num_pos() == 1u );

  const auto tt = simulate<kitty::static_truth_table<4u>>( mig )[0];
  CHECK( tt._bits == 0xecc8 );
}

TEST_CASE( "encode mig_network into mig_index_list", "[index_list]" )
{
  mig_network mig;
  auto const a = mig.create_pi();
  auto const b = mig.create_pi();
  auto const c = mig.create_pi();
  auto const d = mig.create_pi();
  auto const t0 = mig.create_maj( a, b, c );
  auto const t1 = mig.create_maj( t0, b, d );
  mig.create_po( t1 );

  mig_index_list mig_il;
  encode( mig_il, mig );

  CHECK( mig_il.num_pis() == 4u );
  CHECK( mig_il.num_pos() == 1u );
  CHECK( mig_il.num_gates() == 2u );
  CHECK( mig_il.size() == 8u );
  CHECK( mig_il.raw() == std::vector<uint32_t>{ 4 | ( 1 << 8 ) | ( 2 << 16 ), 2, 4, 6, 4, 8, 10, 12 } );
  CHECK( to_index_list_string( mig_il ) == "{4 | 1 << 8 | 2 << 16, 2, 4, 6, 4, 8, 10, 12}" );
}

TEST_CASE( "decode xag_index_list into xag_network", "[index_list]" )
{
  std::vector<uint32_t> const raw_list{ 4 | ( 1 << 8 ) | ( 3 << 16 ), 2, 4, 6, 8, 12, 10, 14 };
  xag_index_list xag_il{ raw_list };

  xag_network xag;
  decode( xag, xag_il );

  CHECK( xag.num_gates() == 3u );
  CHECK( xag.num_pis() == 4u );
  CHECK( xag.num_pos() == 1u );

  const auto tt = simulate<kitty::static_truth_table<4u>>( xag )[0];
  CHECK( tt._bits == 0x7888 );
}

TEST_CASE( "encode xag_network into xag_index_list", "[index_list]" )
{
  xag_network xag;
  auto const a = xag.create_pi();
  auto const b = xag.create_pi();
  auto const c = xag.create_pi();
  auto const d = xag.create_pi();
  auto const t0 = xag.create_and( a, b );
  auto const t1 = xag.create_and( c, d );
  auto const t2 = xag.create_xor( t0, t1 );
  xag.create_po( t2 );

  xag_index_list xag_il;
  encode( xag_il, xag );

  CHECK( xag_il.num_pis() == 4u );
  CHECK( xag_il.num_pos() == 1u );
  CHECK( xag_il.num_gates() == 3u );
  CHECK( xag_il.size() == 8u );
  CHECK( xag_il.raw() == std::vector<uint32_t>{ 4 | ( 1 << 8 ) | ( 3 << 16 ), 2, 4, 6, 8, 12, 10, 14 } );
  CHECK( to_index_list_string( xag_il ) == "{4 | 1 << 8 | 3 << 16, 2, 4, 6, 8, 12, 10, 14}" );
}

TEST_CASE( "decode muxig_index_list into muxig_network", "[index_list]" )
{
  std::vector<uint32_t> const raw_list{ 4 | ( 1 << 8 ) | ( 2 << 16 ), 2, 4, 6, 10, 4, 8, 12 };
  muxig_index_list il{ raw_list };

  muxig_network ntk;
  decode( ntk, il );

  CHECK( ntk.num_gates() == 2u );
  CHECK( ntk.num_pis() == 4u );
  CHECK( ntk.num_pos() == 1u );

  const auto tt = simulate<kitty::static_truth_table<4u>>( ntk )[0];
  CHECK( tt._bits == 0xefc8 );
}

TEST_CASE( "decode mapped_index_list into mapped_network", "[index_list]" )
{
  std::vector<uint32_t> const raw_list{ 2, 1, 2, 2, 3, 2, 2, 2, 3, 4, 2, 4, 5, 5, 6 };
  mapped_index_list list{ raw_list };

  CHECK( raw_list == list.raw() );

  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );
  mapped_network ntk( gates );

  decode( ntk, list );

  CHECK( ntk.num_pis() == 2u );
  CHECK( ntk.num_pos() == 1u );
  CHECK( ntk.num_gates() == 3u );

  const auto tt = simulate<kitty::static_truth_table<2u>>( ntk )[0];
  CHECK( tt._bits == 0x6 );

}

TEST_CASE( "encode mapped_network into mapped_index_list", "[index_list]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );
  mapped_network ntk( gates );

  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const t0 = ntk.create_node( { a, b, c }, 6 );
  auto const t1 = ntk.create_node( { t0, b, d }, 6 );
  ntk.create_po( t1 );

  mapped_index_list list;
  encode( list, ntk );

  CHECK( list.num_pis() == 4u );
  CHECK( list.num_pos() == 1u );
  CHECK( list.size() == 13u );
  CHECK( list.raw() == std::vector<uint32_t>{ 4, 1, 3, 2, 3, 4, 6, 3, 6, 3, 5, 6, 7 } );
  CHECK( to_index_list_string( list ) == "{4, 1, 3, 2, 3, 4, 6, 3, 6, 3, 5, 6, 7}" );
}