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

  kitty::dynamic_truth_table target(4u);
  xminsyn_auto_params ps;
  ps.verbose = false;
  ps.top2_decompose = true;
  ps.top_decompose = true;

  uint32_t ntot = 0;
  do
  {
    uint32_t nin = 4u;

    std::string name = kitty::to_hex( target );
    std::string dotname = "dot/" + name + ".dot";
    std::string verilogname = "verilog/" + name + ".v";

    std::vector<kitty::dynamic_truth_table> xs { nin, kitty::dynamic_truth_table( nin ) };
    xag_network xag;
    std::vector<xag_network::signal> pis;
    for( uint32_t i {0u}; i < nin; ++i )
    {
      pis.push_back( xag.create_pi() );
      kitty::create_nth_var( xs[i], i );
    }

    auto f0 = xminsyn_auto( xag, target, pis, ps );
    xag.create_po(f0);

    xag = cleanup_dangling( xag );
    
    default_simulator<kitty::dynamic_truth_table> sim( nin );
    const auto tt = simulate<kitty::dynamic_truth_table>( xag, sim )[0];

    if( equal( tt, target ) )
    {
      write_dot( xag, dotname );
      write_verilog( xag, verilogname );
      ntot += xag.num_gates();
      std::cout << name << " " << xag.num_gates() << std::endl;
    }
    else
    {
      //std::cout << "error" << std::endl;
      std::cout << "x " << name << " ";
      kitty::print_binary( target );
      std::cout << std::endl;
      //kitty::print_binary( tt );
      //std::cout << std::endl;
      //kitty::print_binary( target );
      //std::cout << std::endl;
      //std::cout << std::endl;
    }    
    kitty::next_inplace( target );

  }while ( !kitty::is_const0( target ) );

  std::cout << ntot << std::endl;
  //std::cout << std::endl;
  //std::ifstream f("tmp.v");

  //if (f.is_open())
  //    std::cout << f.rdbuf();


  /* structural checks */
/*
  kitty::dynamic_truth_table target1( 4u );

  xag_network xag1;
  const auto g1 = xag1.create_pi();
  const auto g2 = xag1.create_pi();
  const auto g3 = xag1.create_pi();
  const auto g4 = xag1.create_pi();

  kitty::dynamic_truth_table x4(4u);

  kitty::create_nth_var( x4, 3 ); 

  kitty::create_majority(target1);

  auto f1 = xminsyn_auto( xag1, target1, { g1, g2, g3, g4 } );
  xag1.create_po(f1);

  xag1 = cleanup_dangling( xag1 );
  write_dot( xag1, "test1.dot" );
    
  default_simulator<kitty::dynamic_truth_table> sim1( 4 );
  const auto tt1 = simulate<kitty::dynamic_truth_table>( xag1, sim1 )[0];
  kitty::print_binary( tt1 );
  std::cout << std::endl;
  kitty::print_binary( target1 );
  std::cout << std::endl;

  std::cout << ( equal( tt1, target1 ) ? " equal " : " different " ) << std::endl; */
}

