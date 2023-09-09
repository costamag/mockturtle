#include <kitty/constructors.hpp>
#include <kitty/print.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/resyn_engines/xag_resyn.hpp>
#include <mockturtle/algorithms/node_resynthesis/exact.hpp>
#include <mockturtle/algorithms/resyn_engines/xag_resyn.hpp>
//#include "../../algorithms/ algorithms/techaware/sym_synthesis.hpp"
#include <mockturtle/algorithms/techaware/sym_synthesis.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/write_dot.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/algorithms/sim_resub.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/detail/switching_activity.hpp>
#include <string>
#include <bit>
#include <bitset>
#include <cstdint>

#include <iostream>
#include <fstream>
#include <string>
#include <set>

#include <sys/types.h>
#include <sys/stat.h>

using namespace mockturtle;
using namespace techaware;

template<class Ntk> void abc_sopmap( Ntk const& );
template<class Ntk> void abc_map( Ntk const& );

int main()
{
  using TT = kitty::dynamic_truth_table;
  using Ntk = aig_network;
  
  // prepare the database for lookup
  xag_npn_resynthesis<Ntk> resyn;
  xag_npn_resynthesis<Ntk, Ntk, xag_npn_db_kind::aig_complete> resyn_complete;
  
  std::string info;
  std::vector<std::vector<double>> npn_fractions;
  TT target(4u);
  std::set<TT> reprs;
  int id=0;

  std::vector<int> delta_covering;
  std::vector<int> delta_sym;
  int count{0};
  do
  {
    const auto repr = kitty::exact_npn_canonization(target);
   if(kitty::is_const0(target)) 
      reprs.insert( std::get<0>(repr) );
    if( kitty::is_const0(target) || reprs.find( std::get<0>(repr) ) != reprs.end() ) 
    {
      kitty::next_inplace( target );
      continue;
    }
    else
    {
      if (kitty::is_const0(target) || kitty::is_const0(~target) )
      {
        continue;
      }
      
      reprs.insert( std::get<0>(repr) );
      printf("[%2d]:", count);
      kitty::print_binary( target );
      //printf("\n");
      //aig_network aig_sat;
//
      //const auto a = aig_sat.create_pi();
      //const auto b = aig_sat.create_pi();
      //const auto c = aig_sat.create_pi();
      //const auto d = aig_sat.create_pi();
      //std::vector<aig_network::signal> pis = { a, b, c, d };
      //exact_aig_resynthesis<aig_network> resyn( false );
      //resyn( aig_sat, target, pis.begin(), pis.end(), [&]( auto const& f ) {
      //  aig_sat.create_po( f );
      //} );
      //printf("========================== EXACT SYNTHESIS ==========================\n");
      //printf("%u\n", aig_sat.num_gates() );
      //abc_sopmap(aig_sat);
      //std::vector<float> saSat = detail::switching_activity( aig_sat, 2048 );
      //printf("%d %d\n", saSat.size(), aig_sat.num_gates() );
      //float sumSat {0};
      //for( auto p : saSat )
      //  sumSat+=p;
//
      //printf("========================== SYMM SYNTHESIS ==========================\n");
//
      aig_network aig_sym;
      std::vector<uint32_t> T {0,0,0,0};
      sym_synthesis<aig_network> synt( target, T );
      std::vector<aig_network::signal> S;
      for( int i{0}; i<4; ++i )
        S.push_back(aig_sym.create_pi());
      if( !synt.net.error )
      {
        aig_sym.create_po(synt.rewrite( &aig_sym, S ));
        uint32_t out_level = synt.get_output_level();
        printf(" %u\n", aig_sym.num_gates() );
        
        default_simulator<kitty::dynamic_truth_table> sim( 4u );
        const auto tt = simulate<kitty::dynamic_truth_table>( aig_sym, sim )[0];

        kitty::print_binary(tt); printf("\n");        
        kitty::print_binary(target); printf("\n");
        assert(kitty::equal(tt, target));
//        abc_map(aig_sym);
//        std::vector<float> saSym = detail::switching_activity( aig_sym, 10000 );
//
//        printf("%d %d\n", saSym.size(), aig_sym.num_gates() );
//        float sumSym {0};
//        for( auto p : saSym )
//          sumSym+=p;
//        printf("ACTIVITY(SYM):%f\n", sumSym);

      }
      count++;
      //printf("====================================================================\n");

      kitty::next_inplace( target );

    }
  } while ( !kitty::is_const0( target ) );


  return 0;
}

template<class Ntk>
void abc_sopmap( Ntk const& ntk )
{
  using namespace mockturtle;

  xag_network res;
  write_blif( ntk, "/tmp/pre.blif" );

  std::string command = "abc -q \"read_library mcnc.genlib; r /tmp/pre.blif; if -g; st; dch; map ; print_stats -p; print_stats -p;\"";
  std::array<char, 128> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }
  printf("%s\n", result.c_str());
}

template<class Ntk>
void abc_map( Ntk const& ntk )
{
  using namespace mockturtle;

  xag_network res;
  write_blif( ntk, "/tmp/pre.blif" );

  std::string command = "abc -q \"read_library mcnc.genlib; r /tmp/pre.blif; st; dch; map; print_stats -p; print_stats -p;\"";
  std::array<char, 128> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }
  printf("%s\n", result.c_str());
}
