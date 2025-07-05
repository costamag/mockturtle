#include <catch.hpp>

#include <cstdint>
#include <sstream>
#include <vector>

#include <lorina/genlib.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/super_reader.hpp>
#include <mockturtle/networks/mapped/bound_network.hpp>
#include <mockturtle/utils/databases/mapped/mapped_database.hpp>
#include <mockturtle/utils/tech_library.hpp>

using namespace mockturtle;

std::string const test_library = "GATE   zero    0 O=CONST0;\n"                                                     // 0
                                 "GATE   one     0 O=CONST1;\n"                                                     // 1
                                 "GATE   inv1    1 O=!a;                      PIN * INV 1 999 0.9 0.3 0.9 0.3\n"    // 2
                                 "GATE   inv2    2 O=!a;                      PIN * INV 2 999 1.0 0.1 1.0 0.1\n"    // 3
                                 "GATE   buf     2 O=a;                       PIN * NONINV 1 999 1.0 0.0 1.0 0.0\n" // 4
                                 "GATE   nand    2 O=!(a*b);                  PIN * INV 1 999 1.0 0.2 1.0 0.2\n"    // 5
                                 "GATE   maj3    8 O=(a*b)+(a*c)+(b*c);       PIN * INV 1 999 3.0 0.4 3.0 0.4\n";   // 6

TEST_CASE( "Adding lists implementing projection to the db", "[mapped_database]" )
{
  using bound_network = mockturtle::bound_network<bound::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  bound::augmented_library<bound::design_type_t::CELL_BASED> lib( gates );
  static constexpr uint32_t MaxNumVars = 4u;
  mapped_database<bound_network, MaxNumVars> db( lib );

  bound_list<bound::design_type_t::CELL_BASED> list0, list1, list2, list3;
  list0.add_inputs( MaxNumVars );
  auto const a = list0.pi_at( 0 );
  auto const b = list0.pi_at( 1 );
  auto const c = list0.pi_at( 2 );
  auto const d = list0.pi_at( 3 );
  list0.add_output( c );

  list1.add_inputs( MaxNumVars );
  list1.add_output( a );
  list2.add_inputs( MaxNumVars );
  list2.add_output( b );
  list3.add_inputs( MaxNumVars );
  list3.add_output( d );

  CHECK( db.add( list0 ) );
  CHECK( !db.add( list1 ) );
  CHECK( !db.add( list2 ) );
  CHECK( !db.add( list3 ) );
}

std::string const symmetric_library =
    "GATE INV                        1.00  Y=!A;                         \n"
    "    PIN  A  UNKNOWN   1 999    15.00     0.00    15.00     0.00     \n"
    "GATE AND2                       2.00  Y=(A * B);                    \n"
    "    PIN  A  UNKNOWN   1 999    25.00     0.00    25.00     0.00     \n"
    "    PIN  B  UNKNOWN   1 999    20.00     0.00    20.00     0.00     \n"
    "GATE MAJ3                       3.00  Y=(A * B) + (A * C) + (B * C);\n"
    "    PIN  A  UNKNOWN   1 999    35.00     0.00    35.00     0.00     \n"
    "    PIN  B  UNKNOWN   1 999    30.00     0.00    30.00     0.00     \n"
    "    PIN  C  UNKNOWN   1 999    25.00     0.00    25.00     0.00     \n"
    "GATE ASYM                       3.00  Y=((!A * B) + C);             \n"
    "    PIN  A  UNKNOWN   1 999    35.00     0.00    35.00     0.00     \n"
    "    PIN  B  UNKNOWN   1 999    30.00     0.00    30.00     0.00     \n"
    "    PIN  C  UNKNOWN   1 999    25.00     0.00    25.00     0.00     \n"
    "GATE AND4                       3.00  Y=((A * B) * (C * D));\n"
    "    PIN  A  UNKNOWN   1 999    35.00     0.00    35.00     0.00     \n"
    "    PIN  B  UNKNOWN   1 999    30.00     0.00    30.00     0.00     \n"
    "    PIN  C  UNKNOWN   1 999    25.00     0.00    25.00     0.00     \n"
    "    PIN  D  UNKNOWN   1 999    45.00     0.00    25.00     0.00     \n"
    "GATE RND4                       3.00  Y=(((!A * B) + C)^D);         \n"
    "    PIN  A  UNKNOWN   1 999    35.00     0.00    35.00     0.00     \n"
    "    PIN  B  UNKNOWN   1 999    30.00     0.00    30.00     0.00     \n"
    "    PIN  C  UNKNOWN   1 999    25.00     0.00    25.00     0.00     \n"
    "    PIN  D  UNKNOWN   1 999    65.00     0.00    25.00     0.00     \n"
    "GATE XOR2                       2.00  Y=(A ^ B);                    \n"
    "    PIN  A  UNKNOWN   1 999    25.00     0.00    25.00     0.00     \n"
    "    PIN  B  UNKNOWN   1 999    20.00     0.00    20.00     0.00     \n"
    "GATE FA                       3.00  C=(A * B) + (A * C) + (B * C);  \n"
    "    PIN  A  UNKNOWN   1 999    35.00     0.00    35.00     0.00     \n"
    "    PIN  B  UNKNOWN   1 999    30.00     0.00    30.00     0.00     \n"
    "    PIN  C  UNKNOWN   1 999    25.00     0.00    25.00     0.00     \n"
    "GATE FA                       3.00  S=( (A ^ B) ^ C );              \n"
    "    PIN  A  UNKNOWN   1 999    35.00     0.00    35.00     0.00     \n"
    "    PIN  B  UNKNOWN   1 999    30.00     0.00    30.00     0.00     \n"
    "    PIN  C  UNKNOWN   1 999    25.00     0.00    25.00     0.00     \n";

TEST_CASE( "Inserting lists with one-input node in mapped databases", "[mapped_database]" )
{
  using bound_network = mockturtle::bound_network<bound::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in( symmetric_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  bound::augmented_library<bound::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  mapped_database<bound_network, MaxNumVars> db( lib );

  bool first = true;
  for ( auto i = 0u; i < MaxNumVars; ++i )
  {
    bound_list<bound::design_type_t::CELL_BASED> list;
    list.add_inputs( MaxNumVars );
    list.add_output( list.add_gate( { i }, 0 ) );
    CHECK( !( first ^ db.add( list ) ) );
    first = false;
  }
}

TEST_CASE( "Inserting lists with two-input node in mapped databases", "[mapped_database]" )
{
  using bound_network = mockturtle::bound_network<bound::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in( symmetric_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  bound::augmented_library<bound::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  mapped_database<bound_network, MaxNumVars> db( lib );

  bool first = true;
  for ( auto i = 0u; i < MaxNumVars; ++i )
  {
    for ( auto j = 0u; j < MaxNumVars; ++j )
    {
      if ( i == j )
        continue;

      bound_list<bound::design_type_t::CELL_BASED> list;
      list.add_inputs( MaxNumVars );
      list.add_output( list.add_gate( { i, j }, 1 ) );
      CHECK( !( first ^ db.add( list ) ) );
      first = false;
    }
  }
}

TEST_CASE( "Inserting symmetric single-node lists with three inputs in mapped databases", "[mapped_database]" )
{
  using bound_network = mockturtle::bound_network<bound::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in( symmetric_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  bound::augmented_library<bound::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  mapped_database<bound_network, MaxNumVars> db( lib );

  bool first = true;
  for ( auto i = 0u; i < MaxNumVars; ++i )
  {
    for ( auto j = 0u; j < MaxNumVars; ++j )
    {
      if ( i == j )
        continue;

      for ( auto k = 0u; k < MaxNumVars; ++k )
      {
        if ( ( i == k ) || ( j == k ) )
          continue;

        bound_list<bound::design_type_t::CELL_BASED> list;
        list.add_inputs( MaxNumVars );
        list.add_output( list.add_gate( { i, j, k }, 2 ) );
        if ( first )
          CHECK( db.add( list ) );
        else
          CHECK( !db.add( list ) );
        first = false;
      }
    }
  }
}

TEST_CASE( "Inserting asymmetric single-node lists with three inputs in mapped databases", "[mapped_database]" )
{
  using bound_network = mockturtle::bound_network<bound::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in( symmetric_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  bound::augmented_library<bound::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  mapped_database<bound_network, MaxNumVars> db( lib );

  bool first = true;
  for ( auto i = 0u; i < MaxNumVars; ++i )
  {
    for ( auto j = 0u; j < MaxNumVars; ++j )
    {
      if ( i == j )
        continue;

      for ( auto k = 0u; k < MaxNumVars; ++k )
      {
        if ( ( i == k ) || ( j == k ) )
          continue;

        bound_list<bound::design_type_t::CELL_BASED> list;
        list.add_inputs( MaxNumVars );
        list.add_output( list.add_gate( { i, j, k }, 3 ) );
        CHECK( !( first ^ db.add( list ) ) );
        first = false;
      }
    }
  }
}

TEST_CASE( "Inserting symmetric single-node lists with 4 inputs in mapped databases", "[mapped_database]" )
{
  using bound_network = mockturtle::bound_network<bound::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in( symmetric_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  bound::augmented_library<bound::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  mapped_database<bound_network, MaxNumVars> db( lib );

  bool first = true;
  for ( auto i = 0u; i < MaxNumVars; ++i )
  {
    for ( auto j = 0u; j < MaxNumVars; ++j )
    {
      if ( i == j )
        continue;

      for ( auto k = 0u; k < MaxNumVars; ++k )
      {
        if ( ( i == k ) || ( j == k ) )
          continue;

        for ( auto l = 0u; l < MaxNumVars; ++l )
        {
          if ( ( i == l ) || ( j == l ) || ( k == l ) )
            continue;

          bound_list<bound::design_type_t::CELL_BASED> list;
          list.add_inputs( MaxNumVars );
          list.add_output( list.add_gate( { i, j, k, l }, 4 ) );
          CHECK( !( first ^ db.add( list ) ) );
          first = false;
        }
      }
    }
  }
}

TEST_CASE( "Inserting asymmetric single-node lists with 4 inputs in mapped databases", "[mapped_database]" )
{
  using bound_network = mockturtle::bound_network<bound::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in( symmetric_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  bound::augmented_library<bound::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  mapped_database<bound_network, MaxNumVars> db( lib );

  bool first = true;
  for ( auto i = 0u; i < MaxNumVars; ++i )
  {
    for ( auto j = 0u; j < MaxNumVars; ++j )
    {
      if ( i == j )
        continue;

      for ( auto k = 0u; k < MaxNumVars; ++k )
      {
        if ( ( i == k ) || ( j == k ) )
          continue;

        for ( auto l = 0u; l < MaxNumVars; ++l )
        {
          if ( ( i == l ) || ( j == l ) || ( k == l ) )
            continue;

          bound_list<bound::design_type_t::CELL_BASED> list;
          list.add_inputs( MaxNumVars );
          list.add_output( list.add_gate( { i, j, k, l }, 5 ) );
          CHECK( !( first ^ db.add( list ) ) );
          first = false;
        }
      }
    }
  }
}

TEST_CASE( "Inserting two nodes list in database", "[mapped_database]" )
{
  using bound_network = mockturtle::bound_network<bound::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in( symmetric_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  bound::augmented_library<bound::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  mapped_database<bound_network, MaxNumVars> db( lib );

  bool first = true;
  for ( auto i = 0u; i < MaxNumVars; ++i )
  {
    for ( auto j = 0u; j < MaxNumVars; ++j )
    {
      if ( i == j )
        continue;

      for ( auto k = 0u; k < MaxNumVars; ++k )
      {
        if ( ( i == k ) || ( j == k ) )
          continue;

        for ( auto l = 0u; l < MaxNumVars; ++l )
        {
          if ( ( i == l ) || ( j == l ) || ( k == l ) )
            continue;

          bound_list<bound::design_type_t::CELL_BASED> list;
          list.add_inputs( MaxNumVars );
          auto l0 = list.add_gate( { i, j }, 1 );
          auto l1 = list.add_gate( { l0, k, l }, 3 );
          list.add_output( l1 );
          CHECK( !( first ^ db.add( list ) ) );
          first = false;
        }
      }
    }
  }
}

TEST_CASE( "Dominant and dominated lists in mapped database", "[mapped_database]" )
{
  using bound_network = mockturtle::bound_network<bound::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in( symmetric_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  bound::augmented_library<bound::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  mapped_database<bound_network, MaxNumVars> db( lib );
  bound_list<bound::design_type_t::CELL_BASED> list1, list2, list3;
  list1.add_inputs( MaxNumVars );
  list2.add_inputs( MaxNumVars );
  list3.add_inputs( MaxNumVars );
  auto const l1_1 = list1.add_gate( { 1 }, 0 );
  auto const l1_2 = list1.add_gate( { 5 }, 0 );
  auto const l1_3 = list1.add_gate( { l1_1, 5 }, 1 );
  auto const l1_4 = list1.add_gate( { l1_2, 1 }, 1 );
  auto const l1_5 = list1.add_gate( { l1_3, l1_4 }, 6 );
  list1.add_output( l1_5 );

  auto const l2_1 = list2.add_gate( { 4, 0 }, 6 );
  list2.add_output( l2_1 );

  list3 = list1;

  CHECK( db.size() == 0 );
  CHECK( db.num_rows() == 0 );
  CHECK( db.add( list1 ) );
  CHECK( db.num_rows() == 1 );
  CHECK( db.size() == 1 );
  CHECK( db.add( list2 ) );
  CHECK( db.size() == 1 );
  CHECK( db.num_rows() == 1 );
  CHECK( !db.add( list3 ) );
  CHECK( db.num_rows() == 1 );
  CHECK( db.size() == 1 );
}

TEST_CASE( "Saving a mapped database", "[mapped_database]" )
{
  using bound_network = mockturtle::bound_network<bound::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in( symmetric_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  bound::augmented_library<bound::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  mapped_database<bound_network, MaxNumVars> db( lib );
  bound_list<bound::design_type_t::CELL_BASED> list1, list2, list3, list4;
  list1.add_inputs( MaxNumVars );
  list2.add_inputs( MaxNumVars );
  list3.add_inputs( MaxNumVars );
  list4.add_inputs( MaxNumVars );
  auto const l1_1 = list1.add_gate( { 1 }, 0 );
  auto const l1_2 = list1.add_gate( { 5 }, 0 );
  auto const l1_3 = list1.add_gate( { l1_1, 5 }, 1 );
  auto const l1_4 = list1.add_gate( { l1_2, 1 }, 1 );
  auto const l1_5 = list1.add_gate( { l1_3, l1_4 }, 6 );
  list1.add_output( l1_5 );

  auto const l2_1 = list2.add_gate( { 4, 0 }, 6 );
  list2.add_output( l2_1 );

  auto const l3_1 = list3.add_gate( { 1, 5, 2, 0 }, 4 );
  auto const l3_2 = list3.add_gate( { l3_1 }, 0 );
  auto const l3_3 = list3.add_gate( { 3, l3_2 }, 1 );
  list3.add_output( l3_3 );

  auto const l4_1 = list4.add_gate( { 2, 0, 3, 1 }, 4 );
  auto const l4_2 = list4.add_gate( { l4_1 }, 0 );
  auto const l4_3 = list4.add_gate( { 4, l4_2 }, 1 );
  list4.add_output( l4_3 );

  CHECK( db.add( list1 ) );
  CHECK( db.num_rows() == 1 );
  CHECK( db.size() == 1 );
  CHECK( db.add( list2 ) );
  CHECK( db.num_rows() == 1 );
  CHECK( db.size() == 1 );
  CHECK( db.add( list3 ) );
  CHECK( db.num_rows() == 2 );
  CHECK( !db.add( list4 ) );
  CHECK( db.num_rows() == 2 );
  CHECK( db.size() == 2 );
  std::stringstream out;
  db.commit( out );
  std::string const expected =
      "module top( x0 , x1 , x2 , x3 , x4 , x5 , y0 , y1 );\n"
      "  input x0 , x1 , x2 , x3 , x4 , x5 ;\n"
      "  output y0 , y1 ;\n"
      "  wire n10 , n12 ;\n"
      "  XOR2 g0( .A (x4), .B (x5), .Y (y0) );\n"
      "  AND4 g1( .A (x3), .B (x4), .C (x5), .D (x2), .Y (n12) );\n"
      "  INV  g2( .A (n12), .Y (n10) );\n"
      "  AND2 g3( .A (x0), .B (n10), .Y (y1) );\n"
      "endmodule\n";
  CHECK( out.str() == expected );
}
