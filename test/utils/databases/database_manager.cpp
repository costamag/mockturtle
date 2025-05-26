#include <catch.hpp>

#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/utils/databases/database_manager.hpp>
#include <mockturtle/utils/index_list/index_list.hpp>
#include <mockturtle/utils/index_list/list_simulator.hpp>

using namespace mockturtle;

template<typename Ntk, typename List>
void test_npn_lookup()
{
  static constexpr uint32_t num_vars = 4u;
  using TT = kitty::static_truth_table<num_vars>;
  TT onset;
  database_manager<Ntk> mng;
  std::array<TT, num_vars> xs;
  std::vector<TT const*> xs_ptrs;
  for ( auto i = 0u; i < num_vars; ++i )
  {
    kitty::create_nth_var( xs[i], i );
    xs_ptrs.push_back( &xs[i] );
  }

  list_simulator<List, TT> sim_list;
  do
  {
    /* define the functionality */
    kitty::next_inplace( onset );
    kitty::ternary_truth_table<TT> tt( onset );

    /* boolean matching */
    auto info = mng.lookup_npn( tt );
    /* verify that at least a match is found */
    CHECK( info );
    Ntk ntk;
    std::vector<signal<Ntk>> pis;
    for ( auto i = 0u; i < num_vars; ++i )
    {
      pis.push_back( ntk.create_pi() );
    }
    /* consider all the sub-networks matching the functionality */
    info->foreach_entry( [&]( auto f ) {
      /* insert the sub-network into a network and simulate it */
      auto o = mng.insert( *info, ntk, f, pis.begin(), pis.end() );
      auto no = ntk.get_node( o );

      default_simulator<TT> sim;
      auto tts = simulate_nodes<TT, Ntk>( ntk, sim );
      auto res = ntk.is_complemented( o ) ? ~tts[no] : tts[no];
      CHECK( kitty::equal( res, onset ) );

      /* insert the sub-network into a list and simulate it */
      {
        List list( num_vars );
        std::vector<uint32_t> pis_list{ 2, 4, 6, 8 };
        auto lit_out = mng.insert( *info, list, f, pis_list.begin(), pis_list.end() );
        sim_list( list, xs_ptrs );
        auto [p_res_list, is_compl] = sim_list.get_simulation( list, xs_ptrs, lit_out );
        TT res_list = *p_res_list;
        if ( is_compl )
        {
          res_list = ~res_list;
        }
        CHECK( kitty::equal( res_list, onset ) );
      }
    } );
  } while ( !kitty::is_const0( onset ) );
}

TEST_CASE( "database for aig_network", "[database_manager]" )
{
  test_npn_lookup<aig_network, xag_index_list<true>>();
}

TEST_CASE( "database for xag_network", "[database_manager]" )
{
  test_npn_lookup<xag_network, xag_index_list<true>>();
}

TEST_CASE( "database for mig_network", "[database_manager]" )
{
  test_npn_lookup<mig_network, mig_index_list>();
}