#include <catch.hpp>

#include <sstream>
#include <string>

#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/networks/buffered.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/muxig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/networks/xmg.hpp>

#include <kitty/kitty.hpp>
#include <lorina/verilog.hpp>

using namespace mockturtle;

TEST_CASE( "read a VERILOG file into MIG network", "[verilog_reader]" )
{
  mig_network mig;

  std::string file{
      "module top( y1, y2, a, b, c ) ;\n"
      "  input a , b , c ;\n"
      "  output y1 , y2 ;\n"
      "  wire zero, g0, g1 , g2 , g3 , g4 ;\n"
      "  assign zero = 0 ;\n"
      "  assign g0 = a ;\n"
      "  assign g1 = ~c ;\n"
      "  assign g2 = g0 & g1 ;\n"
      "  assign g3 = a | g2 ;\n"
      "  assign g4 = ( ~a & b ) | ( ~a & c ) | ( b & c ) ;\n"
      "  assign g5 = g2 ^ g3 ^ g4;\n"
      "  assign g6 = ~( g4 & g5 );\n"
      "  assign y1 = g3 ;\n"
      "  assign y2 = g4 ;\n"
      "endmodule\n" };

  std::istringstream in( file );
  auto const result = lorina::read_verilog( in, verilog_reader( mig ) );

  /* structural checks */
  CHECK( result == lorina::return_code::success );
  CHECK( mig.size() == 11 );
  CHECK( mig.num_pis() == 3 );
  CHECK( mig.num_pos() == 2 );
  CHECK( mig.num_gates() == 7 );

  /* functional checks */
  default_simulator<kitty::dynamic_truth_table> sim( mig.num_pis() );
  const auto tts = simulate<kitty::dynamic_truth_table>( mig, sim );
  mig.foreach_po( [&]( auto const&, auto i ) {
    switch ( i )
    {
    case 0:
      CHECK( kitty::to_hex( tts[i] ) == "aa" );
      break;
    case 1:
      CHECK( kitty::to_hex( tts[i] ) == "d4" );
      break;
    }
  } );
}

TEST_CASE( "read a VERILOG file into XMG network", "[verilog_reader]" )
{
  xmg_network xmg;

  std::string file{
      "module top( y1, y2, a, b, c ) ;\n"
      "  input a , b , c ;\n"
      "  output y1 , y2 ;\n"
      "  wire zero, g0, g1 , g2 , g3 , g4 ;\n"
      "  assign zero = 0 ;\n"
      "  assign g0 = a ;\n"
      "  assign g1 = ~c ;\n"
      "  assign g2 = g0 & g1 ;\n"
      "  assign g3 = a | g2 ;\n"
      "  assign g4 = ( ~a & b ) | ( ~a & c ) | ( b & c ) ;\n"
      "  assign g5 = g2 ^ g3 ^ g4;\n"
      "  assign g6 = ~( g4 & g5 );\n"
      "  assign y1 = g3 ;\n"
      "  assign y2 = g4 ;\n"
      "endmodule\n" };

  std::istringstream in( file );
  auto const result = lorina::read_verilog( in, verilog_reader( xmg ) );

  /* structural checks */
  CHECK( result == lorina::return_code::success );
  CHECK( xmg.size() == 9 );
  CHECK( xmg.num_pis() == 3 );
  CHECK( xmg.num_pos() == 2 );
  CHECK( xmg.num_gates() == 5 );

  /* functional checks */
  default_simulator<kitty::dynamic_truth_table> sim( xmg.num_pis() );
  const auto tts = simulate<kitty::dynamic_truth_table>( xmg, sim );
  xmg.foreach_po( [&]( auto const&, auto i ) {
    switch ( i )
    {
    case 0:
      CHECK( kitty::to_hex( tts[i] ) == "aa" );
      break;
    case 1:
      CHECK( kitty::to_hex( tts[i] ) == "d4" );
      break;
    }
  } );
}

TEST_CASE( "read a VERILOG file into MuxIG network", "[verilog_reader]" )
{
  muxig_network ntk;

  std::string file{
      "module top( y1, a, b, c ) ;\n"
      "  input a , b , c ;\n"
      "  output y1 ;\n"
      "  wire zero, g1 , g2 , g3 , g4 ;\n"
      "  assign g1 = a & b ;\n"
      "  assign g2 = a | b ;\n"
      "  assign g3 = ~g2 ;\n"
      "  assign g4 = c ? g1 : g3 ;\n"
      "  assign y1 = g4 ;\n"
      "endmodule\n" };

  std::istringstream in( file );
  auto const result = lorina::read_verilog( in, verilog_reader( ntk ) );

  /* structural checks */
  CHECK( result == lorina::return_code::success );
  CHECK( ntk.size() == 7 );
  CHECK( ntk.num_pis() == 3 );
  CHECK( ntk.num_pos() == 1 );
  CHECK( ntk.num_gates() == 3 );

  /* functional checks */
  default_simulator<kitty::dynamic_truth_table> sim( ntk.num_pis() );
  const auto tts = simulate<kitty::dynamic_truth_table>( ntk, sim );
  CHECK( kitty::to_hex( tts[0] ) == "81" );
}

TEST_CASE( "read a VERILOG file with instances", "[verilog_reader]" )
{
  mig_network mig;

  std::string file{
      "module ripple_carry_adder( x1, x2, y );\n"
      "  input x1, x2;\n"
      "  output y;\n"
      "endmodule\n"
      "module top( a, b, c );\n"
      "  input [7:0] a, b ;\n"
      "  output [8:0] c;\n"
      "  ripple_carry_adder #(8) add1(.x1(a), .x2(b), .y(c));\n"
      "endmodule\n" };

  std::istringstream in( file );
  const auto result = lorina::read_verilog( in, verilog_reader( mig ) );
  mig = cleanup_dangling( mig );

  /* structural checks */
  CHECK( result == lorina::return_code::success );
  CHECK( mig.num_pis() == 16 );
  CHECK( mig.num_pos() == 9 );
  CHECK( mig.num_gates() == 32 );
}

TEST_CASE( "read a VERILOG file to create large Montgomery multiplier", "[verilog_reader]" )
{
  xag_network xag;

  std::string file{
      "module montgomery_multiplier( x1, x2, y );\n"
      "  input x1, x2;\n"
      "  output y;\n"
      "endmodule\n"
      "module top( a, b, c );\n"
      "  input [383:0] a, b;\n"
      "  output [383:0] c;\n"
      "  montgomery_multiplier #(384, 384'hfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000ffffffff, 384'h14000000140000000c00000002fffffffcfffffffafffffffbfffffffe00000000000000010000000100000001) mult(.x1(a), .x2(b), .y(c));\n"
      "endmodule\n" };

  std::istringstream in( file );
  verilog_reader reader( xag );
  const auto result = lorina::read_verilog( in, reader );
  xag = cleanup_dangling( xag );

  /* structural checks */
  CHECK( result == lorina::return_code::success );
  CHECK( xag.num_pis() == 768u );
  CHECK( xag.num_pos() == 384u );
  CHECK( xag.num_gates() == 909459u );

  /* name checks */
  CHECK( reader.name() == "top" );
  CHECK( reader.input_names() == std::vector<std::pair<std::string, uint32_t>>{ { { "a", 384 }, { "b", 384 } } } );
  CHECK( reader.output_names() == std::vector<std::pair<std::string, uint32_t>>{ { { "c", 384 } } } );
}

TEST_CASE( "read a VERILOG file with buffers", "[verilog_reader]" )
{
  buffered_mig_network mig;

  std::string file{
      "module buffer( i , o );\n"
      "  input i ;\n"
      "  output o ;\n"
      "endmodule\n"
      "module inverter( i , o );\n"
      "  input i ;\n"
      "  output o ;\n"
      "endmodule\n"
      "module top( x0 , x1 , y0 );\n"
      "  input x0 , x1 ;\n"
      "  output y0 ;\n"
      "  wire n3 , n4 , n5 , n6 ;\n"
      "  buffer  buf_n3( .i (x0), .o (n3) );\n"
      "  buffer  buf_n4( .i (n3), .o (n4) );\n"
      "  assign n5 = ~x1 & ~n4 ;\n"
      "  inverter  inv_n6( .i (n5), .o (n6) );\n"
      "  assign y0 = n6 ;\n"
      "endmodule\n" };

  std::istringstream in( file );
  const auto result = lorina::read_verilog( in, verilog_reader( mig ) );

  /* structural checks */
  CHECK( result == lorina::return_code::success );
  CHECK( mig.num_pis() == 2 );
  CHECK( mig.num_pos() == 1 );
  CHECK( mig.num_gates() == 1 );
  CHECK( mig.size() == 7 ); // 1 constant, 2 PIs, 1 gate, 3 buffers

  /* functional check */
  auto const po_values = simulate_buffered<2>( mig );
  CHECK( po_values[0]._bits == 0xe ); // or
}

TEST_CASE( "read VERILOG into buffered_crossed_klut", "[verilog_reader]" )
{
  buffered_crossed_klut_network ntk;

  std::string file{
      "module buffer( i , o );\n"
      "  input i ;\n"
      "  output o ;\n"
      "endmodule\n"
      "module inverter( i , o );\n"
      "  input i ;\n"
      "  output o ;\n"
      "endmodule\n"
      "module crossing( i1 , i2 , o1 , o2 );\n"
      "  input i1 , i2 ;\n"
      "  output o1 , o2 ;\n"
      "endmodule\n"
      "module top( x0 , x1 , y0 , y1 , y2 );\n"
      "  input x0 , x1 ;\n"
      "  output y0 , y1 , y2 ;\n"
      "  wire n4 , n5 , n6 , n7 , n8 , n9 , n10 , n11 , n12 , n13 ;\n"
      "  buffer buf_n4( .i (x0), .o (n4) );\n"
      "  crossing cross_n5( .i1 (x0), .i2 (x1), .o1 (n5_1), .o2 (n5_2) );\n"
      "  buffer buf_n6( .i (x1), .o (n6) );\n"
      "  buffer buf_n7( .i (n4), .o (n7) );\n"
      "  crossing cross_n8( .i1 (n4), .i2 (n5_2), .o1 (n8_1), .o2 (n8_2) );\n"
      "  crossing cross_n9( .i1 (n5_1), .i2 (n6), .o1 (n9_1), .o2 (n9_2) );\n"
      "  inverter inv_n10( .i (n6), .o (n10) );\n"
      "  assign n11 = ~n7 | ~n8_2 ;\n"
      "  assign n12 = n8_1 | n9_2 ;\n"
      "  assign n13 = n9_1 ^ n10 ;\n"
      "  assign y0 = n11 ;\n"
      "  assign y1 = n12 ;\n"
      "  assign y2 = n13 ;\n"
      "endmodule\n" };

  std::istringstream in( file );
  const auto result = lorina::read_verilog( in, verilog_reader( ntk ) );

  /* structural checks */
  CHECK( result == lorina::return_code::success );
  CHECK( ntk.num_pis() == 2 );
  CHECK( ntk.num_pos() == 3 );
  CHECK( ntk.size() == 14 ); // 2 constants, 2 PIs, 3 buffers, 1 inverter, 3 crossings, 3 gates

  /* functional check */
  auto const po_values = simulate_buffered<2>( ntk );
  CHECK( po_values[0]._bits == 0x7 ); // nand
  CHECK( po_values[1]._bits == 0xe ); // or
  CHECK( po_values[2]._bits == 0x9 ); // xnor
}

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

TEST_CASE( "Read structural verilog to mapped network", "[verilog_reader]" )
{

  using bound_network = mockturtle::bound_network<bound::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in_lib( test_library );
  auto result_lib = lorina::read_genlib( in_lib, genlib_reader( gates ) );
  CHECK( result_lib == lorina::return_code::success );

  std::string file{
      "module top( x0 , x1 , x2 , y0 , y1 , y2, y3 );\n"
      "  input x0 , x1, x2 ;\n"
      "  output y0 , y1 , y2, y3 ;\n"
      "  wire n4 , n5 , n6 ;\n"
      "  inv1 g0( .a (x0), .O (n4) );\n"
      "  fa   g1( .a (n4), .b (x1), .c (x2), .C (n5), .S (n6) );\n"
      "  inv1 g2( .a (n4), .O (y0) );\n"
      "  xor2 g3( .a (n6), .b (x2), .O (y1) );\n"
      "  buf g4( .a (n5), .O (y2) );\n"
      "  buf g5( .a (n6), .O (y3) );\n"
      "endmodule\n" };

  std::istringstream in_ntk( file );

  bound_network ntk( gates );
  const auto result_ntk = lorina::read_verilog( in_ntk, verilog_reader( ntk ) );

  /* structural checks */
  CHECK( result_ntk == lorina::return_code::success );
  CHECK( ntk.num_pis() == 3 );
  CHECK( ntk.num_pos() == 4 );
  CHECK( ntk.size() == 11 ); // 2 constants, 2 PIs, 3 buffers, 1 inverter, 3 crossings, 3 gates

  CHECK( ntk.num_gates() == 6 );

  std::ostringstream out;
  write_verilog( ntk, out );

  std::string expected = "module top( x0 , x1 , x2 , y0 , y1 , y2 , y3 );\n"
                         "  input x0 , x1 , x2 ;\n"
                         "  output y0 , y1 , y2 , y3 ;\n"
                         "  wire n5 , n6_0 , n6_1 ;\n"
                         "  inv1  g0( .a (x0), .O (n5) );\n"
                         "  inv1  g1( .a (n5), .O (y0) );\n"
                         "  fa    g2( .a (n5), .b (x1), .c (x2), .C (n6_0), .S (n6_1) );\n"
                         "  xor2  g3( .a (n6_1), .b (x2), .O (y1) );\n"
                         "  buf   g4( .a (n6_0), .O (y2) );\n"
                         "  buf   g5( .a (n6_1), .O (y3) );\n"
                         "endmodule\n";
  CHECK( out.str() == expected );
}

TEST_CASE( "Read structural verilog to mapped network  with the inputs permutated", "[verilog_reader]" )
{

  using bound_network = mockturtle::bound_network<bound::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in_lib( test_library );
  auto result_lib = lorina::read_genlib( in_lib, genlib_reader( gates ) );
  CHECK( result_lib == lorina::return_code::success );

  std::string file{
      "module top( x0 , x1 , x2 , y0 , y1 , y2, y3 );\n"
      "  input x0 , x1, x2 ;\n"
      "  output y0 , y1 , y2, y3 ;\n"
      "  wire n4 , n5 , n6 ;\n"
      "  inv1 g0( .a (x0), .O (n4) );\n"
      "  fa   g1( .b (x1), .C (n5), .a (n4), .c (x2), .S (n6) );\n"
      "  inv1 g2( .a (n4), .O (y0) );\n"
      "  xor2 g3( .a (n6), .O (y1), .b (x2) );\n"
      "  buf g4( .a (n5), .O (y2) );\n"
      "  buf g5( .O (y3), .a (n6) );\n"
      "endmodule\n" };

  std::istringstream in_ntk( file );

  bound_network ntk( gates );
  const auto result_ntk = lorina::read_verilog( in_ntk, verilog_reader( ntk ) );

  /* structural checks */
  CHECK( result_ntk == lorina::return_code::success );
  CHECK( ntk.num_pis() == 3 );
  CHECK( ntk.num_pos() == 4 );
  CHECK( ntk.size() == 11 ); // 2 constants, 2 PIs, 3 buffers, 1 inverter, 3 crossings, 3 gates

  CHECK( ntk.num_gates() == 6 );

  std::ostringstream out;
  write_verilog( ntk, out );

  std::string expected = "module top( x0 , x1 , x2 , y0 , y1 , y2 , y3 );\n"
                         "  input x0 , x1 , x2 ;\n"
                         "  output y0 , y1 , y2 , y3 ;\n"
                         "  wire n5 , n6_0 , n6_1 ;\n"
                         "  inv1  g0( .a (x0), .O (n5) );\n"
                         "  inv1  g1( .a (n5), .O (y0) );\n"
                         "  fa    g2( .a (n5), .b (x1), .c (x2), .C (n6_0), .S (n6_1) );\n"
                         "  xor2  g3( .a (n6_1), .b (x2), .O (y1) );\n"
                         "  buf   g4( .a (n6_0), .O (y2) );\n"
                         "  buf   g5( .a (n6_1), .O (y3) );\n"
                         "endmodule\n";
  CHECK( out.str() == expected );
}