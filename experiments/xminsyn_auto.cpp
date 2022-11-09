#include <kitty/constructors.hpp>
#include <mockturtle/algorithms/sfps/bottomup/xminsyn_auto.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/io/write_dot.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <iostream>
#include <fstream>

using namespace mockturtle;

int main()
{
  std::cout << " [1] for example, [2] for arbitrary truth table ";
  uint32_t inp;
  std::cin >> inp;
  if( inp == 1u )
  {
    kitty::dynamic_truth_table table( 3u );

    xag_network xag;
    const auto s1 = xag.create_pi();
    const auto s2 = xag.create_pi();
    const auto s3 = xag.create_pi();

    kitty::dynamic_truth_table x1(3u);
    kitty::dynamic_truth_table x2(3u);
    kitty::dynamic_truth_table x3(3u);
    kitty::create_nth_var( x1, 0 ); 
    kitty::create_nth_var( x2, 1 ); 
    kitty::create_nth_var( x3, 2 ); 

    table = (~x2 & ~x1 ) | ( x2 & x1 ) | ( x2 & ~x3 ) | ( x1 & x3 );

    auto f0 = xminsyn_auto( xag, table, { s1, s2, s3 } );
    xag.create_po(f0);

    xag = cleanup_dangling( xag );
    write_dot( xag, "tmp.dot" );
    write_verilog( xag, "tmp.v" );
    
  default_simulator<kitty::dynamic_truth_table> sim( 3 );
  const auto tt = simulate<kitty::dynamic_truth_table>( xag, sim )[0];
  kitty::print_binary( tt );
  std::cout << std::endl;
  kitty::print_binary( table );
  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << ( equal( tt, table ) ? " equal " : " different " ) << std::endl; 
  }
  else
  {
    std::string ttstr;
    std::cout << "enter your truth table" << std::endl;
    std::cin >> ttstr;
    uint32_t nin = uint32_t( log2(ttstr.length()) );
    std::cout << "num vars " << nin << std::endl;

    kitty::dynamic_truth_table table( nin );
    kitty::create_from_binary_string( table, ttstr );
    kitty::print_binary( table );
    std::cout << std::endl;

    std::vector<kitty::dynamic_truth_table> xs { nin, kitty::dynamic_truth_table( nin ) };
    xag_network xag;
    std::vector<xag_network::signal> pis;
    for( uint32_t i {0u}; i < nin; ++i )
    {
      pis.push_back( xag.create_pi() );
      kitty::create_nth_var( xs[i], i );
    }

    auto f0 = xminsyn_auto( xag, table, pis );
    xag.create_po(f0);

    xag = cleanup_dangling( xag );
    write_dot( xag, "tmp.dot" );
    write_verilog( xag, "tmp.v" );
    
    default_simulator<kitty::dynamic_truth_table> sim( nin );
    const auto tt = simulate<kitty::dynamic_truth_table>( xag, sim )[0];
    kitty::print_binary( tt );
    std::cout << std::endl;
    kitty::print_binary( table );
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << ( equal( tt, table ) ? " equal " : " different " ) << std::endl;

  }

  std::cout << std::endl;
  std::ifstream f("tmp.v");

  if (f.is_open())
      std::cout << f.rdbuf();

  /* structural checks */
/*
  kitty::dynamic_truth_table table1( 4u );

  xag_network xag1;
  const auto g1 = xag1.create_pi();
  const auto g2 = xag1.create_pi();
  const auto g3 = xag1.create_pi();
  const auto g4 = xag1.create_pi();

  kitty::dynamic_truth_table x4(4u);

  kitty::create_nth_var( x4, 3 ); 

  kitty::create_majority(table1);

  auto f1 = xminsyn_auto( xag1, table1, { g1, g2, g3, g4 } );
  xag1.create_po(f1);

  xag1 = cleanup_dangling( xag1 );
  write_dot( xag1, "test1.dot" );
    
  default_simulator<kitty::dynamic_truth_table> sim1( 4 );
  const auto tt1 = simulate<kitty::dynamic_truth_table>( xag1, sim1 )[0];
  kitty::print_binary( tt1 );
  std::cout << std::endl;
  kitty::print_binary( table1 );
  std::cout << std::endl;

  std::cout << ( equal( tt1, table1 ) ? " equal " : " different " ) << std::endl; */
}

