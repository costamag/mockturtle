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
#include <mockturtle/resynthesis/networks/cleanup_bound.hpp>
#include <mockturtle/resynthesis/database/mapped_list.hpp>
#include <mockturtle/resynthesis/networks/mapped_to_bound.hpp>
#include <mockturtle/resynthesis/list_simulator.hpp>
#include <mockturtle/utils/tech_library.hpp>


#include <experiments.hpp>

using namespace mockturtle;

template<uint32_t num_vars, class Database_t, class Gates_t>
void fill_with_enum( Database_t&, Gates_t const& gates, std::vector<double> const& times );

bool VERBOSE = false;
int main()
{
  using namespace experiments;
  using namespace mockturtle;

  /* setup the parameters */
  std::string library = "asap7";
  uint32_t constexpr num_vars = 4;
  uint32_t constexpr capacity = 11;
  std::vector<std::vector<double>> patterns { 
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

  std::string file = library + "_enum3_" + std::to_string( num_vars ) + ".db";

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
    fill_with_enum<num_vars>( database, gates, pattern );
  
  database.dump( file );
  return 0;

}

template<class Tt_t, class index_list_t>
void simulate_list( Tt_t& result, index_list_t const& list, std::vector<Tt_t> & tts )
{
  list.foreach_gate( [&]( uint32_t lit0, uint32_t lit1 ){
    uint32_t const index0 = ( lit0 >> 1u ) - 1;
    uint32_t const index1 = ( lit1 >> 1u ) - 1;
    bool const is_compl0 = ( lit0 & 0x1 ) == 0x1;
    bool const is_compl1 = ( lit1 & 0x1 ) == 0x1;
    Tt_t const sim0 = ( lit0 & 0x1 ) == 0x1 ? ~tts[index0] : tts[index0];
    Tt_t const sim1 = ( lit1 & 0x1 ) == 0x1 ? ~tts[index1] : tts[index1];
    tts.emplace_back( sim0 & sim1 );
  } );
  /* the last entry of the list is the output */
  uint32_t const lit = list.po_at(0);
  uint32_t const index = ( lit >> 1u ) - 1;
  bool const is_compl = ( lit & 0x1 ) == 0x1;
  result = is_compl ? ~tts[index] : tts[index];
}

template<class Gate_t>
std::vector<Gate_t> sort_gates( std::vector<Gate_t> const& gates, uint32_t max_inputs, bool use_inv = false )
{
  std::vector<std::pair<Gate_t, double>> tgates;
  for ( int i = 0; i < gates.size(); ++i )
  {
    if ( gates[i].num_vars <= max_inputs )
    {
      if( !use_inv && gates[i].num_vars <= 1 )
        continue;

      double max_delay = std::numeric_limits<double>::max();
      for ( int j = 0; j < gates[i].num_vars; ++j )
      {
        double pin_delay = std::max( gates[i].pins[j].rise_block_delay, gates[i].pins[j].fall_block_delay );
        max_delay = std::min( pin_delay, max_delay );
      }
      tgates.push_back( std::make_pair( gates[i], max_delay ) );
    }
  }

  std::sort( tgates.begin(), tgates.end(), [&]( auto a, auto b ){
    return a.second < b.second;
  } );
  for( auto g : tgates )
  {
    std::cout << g.first.name << " " << g.second << std::endl;
  }
  std::vector<Gate_t> res;
  for ( auto g : tgates )
    res.push_back( g.first );
  return res;
}


template<uint32_t num_vars, class Database_t, class Gates_t>
void fill_with_enum( Database_t& database, Gates_t const& gates, std::vector<double> const& pattern )
{
  double d_MUX = 37;
  double Tmax = 500;
  double Tmin = 10;
  uint32_t steps = std::ceil( Tmax/Tmin );

  using namespace mockturtle;
  using signal = typename bound_network::signal;
  using node = typename bound_network::node;

  /* step 1 : find the P representatives */
  using Tt_t = kitty::static_truth_table<num_vars>;
  using tt_hash = kitty::hash<Tt_t>;
  std::unordered_map<Tt_t, uint32_t, tt_hash> tt_to_po;
  std::unordered_map<uint32_t, Tt_t> po_to_tt;

  /* store circuits for constant (0, 1) and projection functions xi */
  bound_network ntk( gates );

  std::vector<std::vector<uint32_t>> levels( steps );

  Tt_t const0;
  ntk.create_po( ntk.get_constant( false ) );
  ntk.create_po( ntk.get_constant( true ) );
  tt_to_po[const0] = 0;
  po_to_tt[0] = const0;
  tt_to_po[~const0] = 1;
  po_to_tt[1] = ~const0;
  std::array<node, num_vars> pis;
  std::vector<uint32_t> new_pos;
  for ( int i = 0; i < num_vars; ++i )
  {
    auto f = ntk.create_pi();
    pis[i] = ntk.get_node( f );
    Tt_t tt;
    kitty::create_nth_var( tt, i );
    auto index = ntk.create_po( f );
    new_pos.push_back( index );
    tt_to_po[tt] = index;
    po_to_tt[index] = tt;
    kitty::print_binary( tt ); std::cout << std::endl;
  }

  /* generate the arrival time analyzer */
  arrival_analyzer<bound_network> arrival( ntk );
  arrival.compute_arrival_times( pattern );

  for ( int i = 0; i < num_vars; ++i )
  {
    auto n = ntk.pi_at( i );
    auto s = ntk.make_signal( n );
    Tt_t tt;
    kitty::create_nth_var( tt, i );
    tt = ~tt;

    auto s_inv = ntk.create_node( { s }, "INVx1_ASAP7_75t_R")[0];

    auto n_inv = ntk.get_node( s_inv );
    auto index = ntk.create_po( s_inv );
    new_pos.push_back( index );
    tt_to_po[tt] = index;
    po_to_tt[index] = tt;
    kitty::print_binary( tt ); std::cout << std::endl;
  }

  for( int index : new_pos )
  {
    signal s = ntk.po_at( index );
    node n = ntk.get_node( s );
    double time = arrival.get_time( n );
    uint32_t lev = std::ceil( time / Tmin );
    levels[lev].push_back( index );
  }

//  std::sort( old_sigs.begin(), old_sigs.end(), [&]( auto fa, auto fb ){
//    node na = ntk.get_node( fa );
//    node nb = ntk.get_node( fb );
//    return arrival.get_time( na ) < arrival.get_time( nb );
//  } );

  /* Sort L by increasing delay */
  Gates_t ordered_gates = sort_gates( gates, 2 );

  bindings binds( ordered_gates );


  int num_stored = 2*num_vars + 2;
  int num_funcs = 65535;//1u << ( 1u << num_vars );
  double t = 0;
  // level, index
  //while (  )
  //{
    while( num_stored < num_funcs )
    {
      std::cout << num_stored << " / " << num_funcs << std::endl; 
      for ( int ig = 0; ig < ordered_gates.size(); ++ig )//ordered
      {
        auto g = ordered_gates[ig];//ordered
        uint32_t m = g.num_vars;
        if ( m == 2 )
        {
          double del_a = std::max( g.pins[0].rise_block_delay, g.pins[0].fall_block_delay );
          double del_b = std::max( g.pins[1].rise_block_delay, g.pins[1].fall_block_delay );

        for ( int lev1 = 0; lev1 < levels.size(); ++lev1 )
        {
          for ( int i1 = 0; i1 < levels[lev1].size(); ++i1 )
          {
            for ( int lev = 0; lev <= lev1; ++lev )
            {
              int max_i2 = lev == lev1 ? i1 : levels[lev].size();
              for( int i2 = 0; i2 < max_i2; ++i2 )
              {
                auto s1 = ntk.po_at( levels[lev1][i1] );
                auto s2 = ntk.po_at( levels[lev][i2] );
                auto n1 = ntk.get_node( s1 );
                auto n2 = ntk.get_node( s2 );
                node na = n1;
                node nb = n2;

                if ( del_a > del_b )
                {
                  if ( arrival.get_time( n1 ) > arrival.get_time( n2 ) )
                  {
                    na = n2;
                    nb = n1;
                  }
                  else
                  {
                    na = n1;
                    nb = n2;
                  }
                }
                else
                {
                  if ( arrival.get_time( n1 ) > arrival.get_time( n2 ) )
                  {
                    na = n1;
                    nb = n2;
                  }
                  else
                  {
                    na = n2;
                    nb = n1;
                  }
                }

                double ta = arrival.get_time( na );
                double tb = arrival.get_time( nb );
                double tA = std::max( ta + del_a, tb + del_b );
                uint32_t indexA = std::ceil( tA / Tmin );
                uint32_t indexT = std::ceil( t / Tmin );
                if ( indexT != indexA )
                  continue;

                signal sa = ntk.make_signal( na );
                signal sb = ntk.make_signal( nb );

                uint32_t ia = ntk.po_index( sa );
                uint32_t ib = ntk.po_index( sb );

                std::vector<kitty::static_truth_table<num_vars>> xs { po_to_tt[ia], po_to_tt[ib] };
                auto const& list = binds[ig].list;
                kitty::static_truth_table<num_vars> tt;
                simulate_list( tt, list, xs );

                if ( tt_to_po.find( tt ) == tt_to_po.end() )
                {
                  signal f = ntk.create_node( { sa, sb }, g.name )[0];
                  node n = ntk.get_node( f );
                  uint32_t index = ntk.create_po( f );
                  tt_to_po[tt] = index;
                  new_pos.push_back( index );
                  uint32_t lev = std::ceil( tA / Tmin );
                  levels[lev].push_back( index );
                  po_to_tt[index] = tt;
                  num_stored++;
                  kitty::print_binary( tt ); std::cout << " created # " << index << " node " << n << " " << arrival.get_time( n ) << std::endl;
                  database.add_entry( ntk, pis, f );
                  if ( ntk.num_pos() >= num_funcs )
                    return;
                }
              }
            }
          }
        }
        }
        else if( m == 3 )
      {
          double del_a = std::max( g.pins[0].rise_block_delay, g.pins[0].fall_block_delay );
          double del_b = std::max( g.pins[1].rise_block_delay, g.pins[1].fall_block_delay );
          double del_c = std::max( g.pins[2].rise_block_delay, g.pins[2].fall_block_delay );

        for ( int lev1 = 0; lev1 < levels.size(); ++lev1 )
        {
          for ( int i1 = 0; i1 < levels[lev1].size(); ++i1 )
          {
            for ( int lev2 = 0; lev2 <= lev1; ++lev2 )
            {
              int max_i2 = lev2 == lev1 ? i1 : levels[lev2].size();
              for( int i2 = 0; i2 < max_i2; ++i2 )
              {
                for ( int lev3 = 0; lev3 <= lev2; ++lev3 )
                {
                  int max_i3 = lev3 == lev2 ? i2 : levels[lev3].size();
                  for( int i3 = 0; i3 < max_i3; ++i3 )
                  {
                    auto s1 = ntk.po_at( levels[lev1][i1] );
                    auto s2 = ntk.po_at( levels[lev2][i2] );
                    auto s3 = ntk.po_at( levels[lev3][i3] );
                    auto n1 = ntk.get_node( s1 );
                    auto n2 = ntk.get_node( s2 );
                    auto n3 = ntk.get_node( s3 );
                    node na = n1;
                    node nb = n2;
                    node nc = n3;

                    std::vector<node> nds { n1, n2, n3 };
                    std::sort( nds.begin(), nds.end(), [&]( auto a, auto b ){
                      return arrival.get_time( a ) < arrival.get_time( b );
                    } );

                    if ( del_a > del_b && del_b > del_c )
                    {
                      na = nds[0];
                      nb = nds[1];
                      nc = nds[2];
                    }
                    else if ( del_a > del_c && del_c > del_b )
                    {
                      na = nds[0];
                      nb = nds[2];
                      nc = nds[1];
                    }
                    else if ( del_b > del_a && del_a > del_c )
                    {
                      nb = nds[0];
                      na = nds[1];
                      nc = nds[2];
                    }
                    else if ( del_b > del_c && del_c > del_a )
                    {
                      nb = nds[0];
                      nc = nds[1];
                      na = nds[2];
                    }
                    else if ( del_c > del_a && del_a > del_b )
                    {
                      nc = nds[0];
                      na = nds[1];
                      nb = nds[2];
                    }
                    else if ( del_c > del_b && del_b > del_a )
                    {
                      nc = nds[0];
                      nb = nds[1];
                      na = nds[2];
                    }

                    double ta = arrival.get_time( na );
                    double tb = arrival.get_time( nb );
                    double tc = arrival.get_time( nc );
                    double tA = std::max( std::max( ta + del_a, tb + del_b ), tc + del_c );
                    uint32_t indexA = std::ceil( tA / Tmin );
                    uint32_t indexT = std::ceil( t / Tmin );
                    if ( indexT != indexA )
                      continue;

                    signal sa = ntk.make_signal( na );
                    signal sb = ntk.make_signal( nb );
                    signal sc = ntk.make_signal( nc );

                    uint32_t ia = ntk.po_index( sa );
                    uint32_t ib = ntk.po_index( sb );
                    uint32_t ic = ntk.po_index( sc );

                    std::vector<kitty::static_truth_table<num_vars>> xs { po_to_tt[ia], po_to_tt[ib], po_to_tt[ic] };
                    auto const& list = binds[ig].list;
                    kitty::static_truth_table<num_vars> tt;
                    simulate_list( tt, list, xs );

                    if ( tt_to_po.find( tt ) == tt_to_po.end() )
                    {
                      signal f = ntk.create_node( { sa, sb, sc }, g.name )[0];
                      node n = ntk.get_node( f );
                      uint32_t index = ntk.create_po( f );
                      tt_to_po[tt] = index;
                      new_pos.push_back( index );
                      uint32_t lev = std::ceil( tA / Tmin );
                      levels[lev].push_back( index );
                      po_to_tt[index] = tt;
                      num_stored++;
                      kitty::print_binary( tt ); std::cout << " created # " << index << " node " << n << " " << arrival.get_time( n ) << std::endl;
                      database.add_entry( ntk, pis, f );
                      if ( ntk.num_pos() >= num_funcs )
                        return;
                    }
                  }
                }
              }
            }
          }
        }
        }
      }
      t += Tmin;
    }
  //}


//  std::cout << "#pos " << ntk.num_pos() << std::endl;
//  ntk.foreach_po( [&]( auto fo ){
//    std::cout << fo.index << std::endl;
//  } );

}
