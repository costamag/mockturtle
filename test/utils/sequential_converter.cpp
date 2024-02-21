#include <catch.hpp>

#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/aqfp.hpp>
#include <mockturtle/networks/cover.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/sequential.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/networks/xmg.hpp>
#include <mockturtle/utils/sequential_converter.hpp>

using namespace mockturtle;

TEST_CASE( "create and use register in an AIG and convert to combinatorial", "[aig]" )
{
  sequential<aig_network> saig;

  CHECK( has_foreach_po_v<sequential<aig_network>> );
  CHECK( has_create_po_v<sequential<aig_network>> );
  CHECK( has_create_pi_v<sequential<aig_network>> );
  CHECK( has_create_ro_v<sequential<aig_network>> );
  CHECK( has_create_ri_v<sequential<aig_network>> );
  CHECK( has_create_and_v<sequential<aig_network>> );

  const auto x1 = saig.create_pi();
  const auto x2 = saig.create_pi();
  const auto x3 = saig.create_pi();

  CHECK( saig.size() == 4 );
  CHECK( saig.num_registers() == 0 );
  CHECK( saig.num_pis() == 3 );
  CHECK( saig.num_pos() == 0 );

  const auto f1 = saig.create_and( x1, x2 );
  saig.create_po( f1 );
  saig.create_po( !f1 );

  const auto f2 = saig.create_and( f1, x3 );
  saig.create_ri( f2 );

  const auto ro = saig.create_ro();
  saig.create_po( ro );

  CHECK( saig.num_pos() == 3 );
  CHECK( saig.num_registers() == 1 );

  saig.foreach_po( [&]( auto s, auto i ) {
    switch ( i )
    {
    case 0:
      CHECK( s == f1 );
      break;
    case 1:
      CHECK( s == !f1 );
      break;
    case 2:
      // Check if the output (connected to the register) data is the same as the node data being registered.
      CHECK( f2.data == saig.po_at( i ).data );
      break;
    default:
      CHECK( false );
      break;
    }
  } );

  network_converters_stats st;
  aig_network caig = sequential_to_combinatorial( saig, st );

  CHECK( st.num_pis == 3u );
  CHECK( st.num_pos == 3u );
  CHECK( caig.num_pis() == st.num_pis + saig.num_registers() );
  CHECK( caig.num_gates() == saig.num_gates() );
  CHECK( caig.num_pos() == st.num_pos + saig.num_registers() );

  sequential<aig_network> saig2 = combinatorial_to_sequential( caig, st );
  CHECK( saig2.num_pis() == 3u );
  CHECK( saig2.num_pos() == 3u );
  CHECK( saig2.num_registers() == 1u );


}
