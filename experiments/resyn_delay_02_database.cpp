/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2023  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>

#include <mockturtle/algorithms/aig_balancing.hpp>
#include <mockturtle/algorithms/emap.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/rewrite.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/resynthesis/aig_syn.hpp>
#include <mockturtle/resynthesis/algorithms/resynthesize.hpp>
#include <mockturtle/resynthesis/database/easy_database.hpp>
#include <mockturtle/resynthesis/networks/bindings.hpp>
#include <mockturtle/resynthesis/networks/bound.hpp>
#include <mockturtle/resynthesis/database/mapped_list.hpp>
#include <mockturtle/resynthesis/networks/mapped_to_bound.hpp>
#include <mockturtle/resynthesis/list_simulator.hpp>
#include <mockturtle/utils/tech_library.hpp>

#include <experiments.hpp>

using namespace mockturtle;

template<uint32_t num_vars, class Database_t, class Gates_t>
void fill_with_enum( Database_t&, Gates_t const& gates, std::vector<double> const& times );

template<uint32_t num_vars, class Database_t, class Gates_t>
void fill_with_map( Database_t&, Gates_t const& gates, std::vector<double> const& times );

int PERMS [24][4] = {
                      { 0, 1, 2, 3 },
                      { 3, 0, 1, 2 },
                      { 2, 3, 0, 1 },
                      { 1, 2, 3, 0 },
                      { 1, 0, 2, 3 },
                      { 3, 1, 0, 2 },
                      { 2, 3, 1, 0 },
                      { 0, 2, 3, 1 },
                      { 2, 1, 0, 3 },
                      { 3, 2, 1, 0 },
                      { 0, 3, 2, 1 },
                      { 1, 0, 3, 2 },
                      { 3, 1, 2, 0 },
                      { 0, 3, 1, 2 },
                      { 2, 0, 3, 1 },
                      { 1, 2, 0, 3 },
                      { 0, 2, 1, 3 },
                      { 3, 0, 2, 1 },
                      { 1, 3, 0, 2 },
                      { 2, 1, 3, 0 },
                      { 0, 1, 3, 2 },
                      { 2, 0, 1, 3 },
                      { 3, 2, 0, 1 },
                      { 1, 3, 2, 0 }
                    };

bool VERBOSE = false;
int main()
{
  using namespace experiments;
  using namespace mockturtle;

  /* setup the parameters */
  std::string library = "asap7";
  uint32_t constexpr num_vars = 4;
  uint32_t constexpr capacity = 11;
  double d_MUX = 37;
  std::vector<std::vector<double>> patterns { 
                                              //{ 0, 0, 3, 3 },
                                              //{ 0, 0, 0, 3 },
                                              //{ 0, 0, 0, 0 },
                                              //{ 0, 1, 2, 3 },
                                              //{ 0, 3, 3, 3 },
                                              //{ 0, 0, 0, 2 },
                                              //{ 0, 0, 2, 3 },
                                              //{ 0, 0, 0, 1 },
                                              //{ 0, 0, 1, 3 },
                                              //{ 0, 1, 3, 3 },
                                              //{ 0, 0, 1, 2 }
                                              { 0, 0, 0, 0 },
                                              { 111, 111, 0, 0 },
                                              { 111, 0, 111, 0 },
                                              { 111, 0, 0, 111 },
                                              { 0, 111, 111, 0 },
                                              { 0, 111, 0, 111 },
                                              { 0, 0, 111, 111 },
                                              { 0, 0, 0, 111 },
                                              { 0, 0, 111, 0 },
                                              { 0, 111, 0, 0 },
                                              { 111, 0, 0, 0 }
                                            };
//  for ( auto & pattern : patterns )
//  {
//    for( int i = 0; i < num_vars; ++i )
//      pattern[i] *= d_MUX;
//  }

  std::string file = library + "_enum_" + std::to_string( num_vars ) + ".db";

  /* parse the gates */
  std::vector<mockturtle::gate> gates;
  std::ifstream in( "/home/costamag/POST_MAPPING/mockturtle/experiments/cell_libraries/asap7.genlib" );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  if( result != lorina::return_code::success )
  {
    std::cerr << "[e] error in the library parsing" << std::endl;
    return 1;
  }

  /* build the basic database */
  bindings binds( gates );
  mockturtle::easy_database<num_vars, capacity> database( binds );
  for( auto pattern : patterns )
    fill_with_map<num_vars>( database, gates, pattern );
  
  database.dump( file );
  return 0;

}

typename aig_network::signal recursive_synt( aig_network & aig, int * perm, kitty::dynamic_truth_table func, int index )
{
  if( kitty::is_const0( func ) )
    return aig.get_constant( false );
  else if( kitty::is_const0( ~func ) )
    return aig.get_constant( true );
  
  auto x = aig.make_signal( aig.pi_at( perm[index] ) );
  auto tt0 = kitty::cofactor0( func, perm[index] );
  auto tt1 = kitty::cofactor1( func, perm[index] );
  if( kitty::is_const0( tt0 ) )
  {
    auto f1 = recursive_synt( aig, perm, tt1, index+1 );
    return aig.create_and( x, f1 );
  }
  else if( kitty::is_const0( tt1 ) )
  {
    auto f0 = recursive_synt( aig, perm, tt0, index+1 );
    return aig.create_and( !x, f0 );
  }
  else if( kitty::is_const0( ~tt0 ) )
  {
    auto f1 = recursive_synt( aig, perm, tt1, index+1 );
    return aig.create_or( !x, f1 );
  }
  else if( kitty::is_const0( ~tt1 ) )
  {
    auto f0 = recursive_synt( aig, perm, tt0, index+1 );
    return aig.create_or( x, f0 );
  }
  else
  {
    auto f0 = recursive_synt( aig, perm, tt0, index+1 );
    auto f1 = recursive_synt( aig, perm, tt1, index+1 );
    return aig.create_ite( x, f1, f0 );
  }
}

template<uint32_t num_vars, class Database_t, class Gates_t>
void fill_with_map( Database_t& database, Gates_t const& gates, std::vector<double> const& pattern )
{
  using namespace mockturtle;
  using Tt_t = kitty::dynamic_truth_table;
  using tt_hash = kitty::hash<Tt_t>;
  std::unordered_set<Tt_t, tt_hash> classes;
  Tt_t tt(4);
  do {
    const auto res = kitty::exact_p_canonization( tt );
    classes.insert( std::get<0>( res ) );
    kitty::next_inplace( tt );
  } while( !kitty::is_const0( tt ) );

int cnt = 0;
  tech_library_params tps;
  tps.ignore_symmetries = false;
  tps.verbose = true;
  tech_library<9> tech_lib( gates, tps );
  emap_params mps;
  mps.matching_mode = emap_params::hybrid;
  mps.area_oriented_mapping = false;
  mps.map_multioutput = false;
  aig_balancing_params bps;
  bps.minimize_levels = false;
  bps.fast_mode = true;


  for( auto tt : classes ) {
    if( cnt > 30 )
      break;
    std::cout << cnt++ << " / " << classes.size() << std::endl;

    if( !kitty::is_const0( tt ) && !kitty::is_const0( ~tt ) )
    {
if( VERBOSE )
{
  kitty::print_binary( tt ); std::cout << std::endl;
}
      klut_network klut;
      std::vector<typename klut_network::signal> pis;
      for( int i = 0; i < 4; ++i )
      {
        pis.push_back( klut.create_pi() );
      }

      auto const f = klut.create_node( pis, tt );
      klut.create_po( f );

      emap_stats mst;

      aig_network aig;

      for( auto p : PERMS )
      {
        auto fo = recursive_synt( aig, p, tt, 0 );
        aig.create_po( fo );
      }

      cell_view<block_network> mapped = emap<9>( aig, tech_lib, mps, &mst );
      auto bound = convert_mapped_to_bound( mapped, gates );

      arrival_analyzer<bound_network> arrival( bound );

      typename bound_network::signal fbest;
      //for( auto pattern : patterns )
      {
        double best_delay = std::numeric_limits<double>::max();
        arrival.compute_arrival_times( pattern );
        bound.foreach_po( [&]( auto fo ){
          auto no = bound.get_node( fo );
          if( arrival.get_time( no ) < best_delay )
          {
            best_delay = arrival.get_time( no );
            fbest = fo;
          }
        } );

        bound.foreach_pi( [&]( auto f ){
          std::cout << arrival.get_time( bound.get_node( f ) ) << " ";
        } );
        std::cout << best_delay << std::endl;
        
        std::array<typename bound_network::node, num_vars> xs;
        for( int i = 0; i < num_vars; ++i )
          xs[i] = bound.get_node( bound.pi_at( i ) );
        database.add_entry( bound, xs, fbest );
      }
    }

  //  kitty::next_inplace( tt );
  } //while( !kitty::is_const0( tt ) && cnt < 100 );

//  ntk.foreach_po( [&]( auto fo ){
//    database.add_entry( ntk, pis, fo );
//  } );
if( VERBOSE )
  database.print();
}

template<uint32_t num_vars, class Database_t, class Gates_t>
void fill_with_enum( Database_t& database, Gates_t const& gates, std::vector<double> const& pattern )
{
  using namespace mockturtle;
  bound_network ntk( gates );
  std::array<typename bound_network::node, 65536> funcs;
  for( int i = 0; i < 65536; ++i )
    funcs[i] = ntk.get_null();

  std::array<typename bound_network::node, num_vars> pis;
  for( int i = 0; i < num_vars; ++i )
  {
    pis[i] = ntk.get_node( ntk.create_pi() );
  }

  for( int i = 0; i < num_vars; ++i )
  {
    auto inv = ntk.create_node( { ntk.make_signal(ntk.pi_at( i )) }, "INVx1_ASAP7_75t_R")[0];
    ntk.create_po( ntk.make_signal(ntk.pi_at( i )) );
    ntk.create_po( inv );
  }

  std::vector<std::vector<uint64_t>> func_per_level;

  incomplete_node_map<kitty::static_truth_table<num_vars>, bound_network> node_to_sim( ntk );
  std::vector<uint64_t> new_funcs;
  for( int i = 0; i < num_vars; ++i )
  {

    auto n_pos = ntk.get_node( ntk.po_at( 2*i ) );
    auto n_neg = ntk.get_node( ntk.po_at( 2*i + 1 ) );

    kitty::create_nth_var( node_to_sim[n_pos], i );
    node_to_sim[ n_neg ] = ~node_to_sim[ n_pos ];

    uint64_t index_pos = node_to_sim[ n_pos ]._bits;
    funcs[index_pos] = n_pos;

    uint64_t index_neg = node_to_sim[ n_neg ]._bits;
    funcs[index_neg] = n_neg;
    new_funcs.push_back( index_pos );
    new_funcs.push_back( index_neg );
  }
  func_per_level.push_back( new_funcs );

  arrival_analyzer<bound_network> arrival( ntk );
  arrival.compute_arrival_times( pattern );
if( VERBOSE )
{
  for( int i = 0; i < num_vars; ++i )
    std::cout << arrival.get_time( pis[i] ) << " ";
  std::cout << std::endl;
}

  while( ntk.num_pos() < 65536 )
  {
    int level = func_per_level.size();
    new_funcs.clear();
    std::vector<uint32_t> gates_in_order { 4, 6, 7, 5, 9, 8 };
    std::vector<bool> fanin_swaps { false, true, true, true, false, true };
    for( int k=0; k < gates_in_order.size(); ++k )
    {
      if( ntk.num_pos() >= 65536 )
        break;

      auto i = gates_in_order[k];
      if( gates[i].num_vars == 2 )
      {
        std::cout << gates[i].name << std::endl;
        for( int i1 = 0; i1 < func_per_level[level - 1].size(); ++i1 )
        {
          if( ntk.num_pos() >= 65536 )
              break;
          std::cout << i1 << "/" << func_per_level[level - 1].size()-1 << std::endl; 
          for( int level2 = 0; level2 < level; ++level2 )
          {
            if( ntk.num_pos() >= 65536 )
                break;
            for( int i2 = 0; i2 < func_per_level[level2].size(); ++i2 )
            {
              if( ntk.num_pos() >= 65536 )
                break;
              auto index1 = func_per_level[level - 1][i1];
              auto index2 = func_per_level[level2][i2];
              if( index1 != index2 )
              {
                auto n1 = funcs[index1];
                auto n2 = funcs[index2];
                auto f1 = ntk.make_signal( n1 );
                auto f2 = ntk.make_signal( n2 );
                std::vector<typename bound_network::signal> P{ f1, f2 };
                if( fanin_swaps[k] )
                  P = { f2, f1 };
                auto f = ntk.create_node( P, gates[i].name )[0];
                std::vector<kitty::static_truth_table<num_vars>> xs { node_to_sim[n1], node_to_sim[n2] };
                if( fanin_swaps[k] )
                  xs = { node_to_sim[n2], node_to_sim[n1] };
                auto tt = ntk.compute( ntk.get_node( f ), xs.begin(), xs.end() );
//kitty::print_binary( tt );
//std::cout << std::endl;
                uint64_t index = tt._bits;
                if( ntk.is_null( funcs[index] ) )
                {
                  funcs[index] = ntk.get_node( f );
                  new_funcs.push_back( index );
                  node_to_sim.resize();
                  node_to_sim[ntk.get_node( f )] = tt;
                  //auto inv = ntk.create_node( { f }, "INVx1_ASAP7_75t_R")[0];
                  //auto tt_neg = ~tt;
                  //uint64_t index_neg = tt_neg._bits;
                  //funcs[index_neg] = ntk.get_node( inv );
                  //new_funcs.push_back( index_neg );
                  //node_to_sim.resize();
                  //node_to_sim[ntk.get_node( inv )] = tt_neg; 

                  std::cout << "add " << ntk.num_pos() << std::endl; 
                  ntk.create_po( f );

                  //std::cout << "add " << ntk.num_pos() << std::endl; 
                  //ntk.create_po( inv );

                }
                else
                {
                  if( arrival.get_time( ntk.get_node( f ) ) < arrival.get_time( funcs[index] ) )
                  {
                    std::cout << "sub " << arrival.get_time( ntk.get_node( f ) ) << " < " << arrival.get_time( funcs[index] ) << std::endl; 

                    ntk.substitute_node( ntk.get_node( funcs[index] ), f );
                    funcs[index] = ntk.get_node( f );
                    node_to_sim.resize();
                    node_to_sim[ntk.get_node( f )] = tt;

                    //auto inv = ntk.create_node( { f }, "INVx1_ASAP7_75t_R")[0];
                    //auto tt_neg = ~tt;
                    //uint64_t index_neg = tt_neg._bits;
                    //auto n_inv = ntk.get_node( inv );

                    //if( arrival.get_time( n_inv ) < arrival.get_time( funcs[index_neg] ) )
                    //{
                    //  funcs[index_neg] = ntk.get_node( inv );
                    //  node_to_sim.resize();
                    //  node_to_sim[ntk.get_node( inv )] = tt_neg; 
                    //  ntk.substitute_node( ntk.get_node( funcs[index_neg] ), inv );
                    //}


                  }
                }
              }
            }
          }
        }
      }
    }

    for( int i=0; i < gates.size(); ++i )
    {
                  if( ntk.num_pos() >= 65536 )
                    break;
      if( gates[i].num_vars == 3 )
      {
        std::cout << gates[i].name << std::endl;
        for( int i1 = 0; i1 < func_per_level[level - 1].size(); ++i1 )
        {
                  if( ntk.num_pos() >= 65536 )
                    break;
          std::cout << i1 << "/" << func_per_level[level - 1].size()-1 << std::endl; 
          for( int level2 = 0; level2 < level; ++level2 )
          {
                  if( ntk.num_pos() >= 65536 )
                    break;
            for( int i2 = 0; i2 < func_per_level[level2].size(); ++i2 )
            {
                  if( ntk.num_pos() >= 65536 )
                    break;
              for( int level3 = 0; level3 < level; ++level3 )
              {
                  if( ntk.num_pos() >= 65536 )
                    break;
                for( int i3 = 0; i3 < func_per_level[level3].size(); ++i3 )
                {
                  if( ntk.num_pos() >= 65536 )
                    break;

              auto index1 = func_per_level[level - 1][i1];
              auto index2 = func_per_level[level2][i2];
              auto index3 = func_per_level[level3][i3];
              {
                auto n1 = funcs[index1];
                auto n2 = funcs[index2];
                auto n3 = funcs[index3];
                auto f1 = ntk.make_signal( n1 );
                auto f2 = ntk.make_signal( n2 );
                auto f3 = ntk.make_signal( n3 );
                std::vector<typename bound_network::signal> P{ f1, f2, f3 };
                auto f = ntk.create_node( P, gates[i].name )[0];
                std::vector<kitty::static_truth_table<num_vars>> xs { node_to_sim[n1], node_to_sim[n2], node_to_sim[n3] };
                auto tt = ntk.compute( ntk.get_node( f ), xs.begin(), xs.end() );
//kitty::print_binary( tt );
//std::cout << std::endl;
                uint64_t index = tt._bits;
                if( ntk.is_null( funcs[index] ) )
                {
                  funcs[index] = ntk.get_node( f );
                  new_funcs.push_back( index );
                  node_to_sim.resize();
                  node_to_sim[ntk.get_node( f )] = tt;
                  //auto inv = ntk.create_node( { f }, "INVx1_ASAP7_75t_R")[0];
                  //auto tt_neg = ~tt;
                  //uint64_t index_neg = tt_neg._bits;
                  //funcs[index_neg] = ntk.get_node( inv );
                  //new_funcs.push_back( index_neg );
                  //node_to_sim.resize();
                  //node_to_sim[ntk.get_node( inv )] = tt_neg; 

                  std::cout << "add " << ntk.num_pos() << std::endl; 
                  ntk.create_po( f );
                  //std::cout << "add " << ntk.num_pos() << std::endl; 
                  //ntk.create_po( inv );

                }
                else
                {
                  if( arrival.get_time( ntk.get_node( f ) ) < arrival.get_time( funcs[index] ) )
                  {
                    std::cout << "sub " << arrival.get_time( ntk.get_node( f ) ) << " < " << arrival.get_time( funcs[index] ) << std::endl; 

                    ntk.substitute_node( ntk.get_node( funcs[index] ), f );
                    funcs[index] = ntk.get_node( f );
                    node_to_sim.resize();
                    node_to_sim[ntk.get_node( f )] = tt;

                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }


    func_per_level.push_back( new_funcs );

  }

  ntk.foreach_po( [&]( auto fo ){
    database.add_entry( ntk, pis, fo );
  } );

  database.print();
}
