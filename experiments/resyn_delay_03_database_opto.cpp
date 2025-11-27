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
void fill_with_opto( Database_t&, Gates_t const& gates, kitty::dynamic_truth_table tt, std::vector<double> const& times );

template<class Tt_t>
void simulate_tfi( mockturtle::bound_network & ntk, mockturtle::bound_network::signal f, std::unordered_map<mockturtle::bound_network::node, Tt_t> & map )
{
  auto n = ntk.get_node( f );
  if ( map.find( n ) != map.end() )
    return;
  
  std::vector<Tt_t> tts;
  ntk.foreach_fanin( n, [&]( auto fi, auto ii ){
    simulate_tfi( ntk, fi, map );
    auto ni = ntk.get_node( fi );
    tts.push_back( map[ni] );
  } );
  map[n] = ntk.compute( n, tts.begin(), tts.end() );
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

bool VERBOSE = false;
int main()
{
  using namespace experiments;
  using namespace mockturtle;

  /* setup the parameters */
  std::string library = "asap7";
  uint32_t constexpr num_vars = 4;
  uint32_t constexpr capacity = 12;
  double d_MUX = 37;
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
//  for ( auto & pattern : patterns )
//  {
//    for( int i = 0; i < num_vars; ++i )
//      pattern[i] *= d_MUX;
//  }

  std::string file = library + "_enum5_" + std::to_string( num_vars ) + ".db";
  std::string file_out = library + "_enum6_" + std::to_string( num_vars ) + ".db";

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
  mockturtle::easy_database<num_vars, capacity> database( binds, file );
  mockturtle::easy_database<num_vars, capacity> database_out( binds );

  auto ordered_gates = sort_gates( gates, 6 );
  bindings binds_o( ordered_gates );

  /* Enumerate all the P-classes */
  using Tt_t = kitty::static_truth_table<num_vars>;
  using tt_hash = kitty::hash<Tt_t>;

  for ( int p = 0; p < patterns.size(); ++p )
  {
  std::unordered_map<typename bound_network::node, Tt_t> node_to_tt;
  std::unordered_map<Tt_t, typename bound_network::node, tt_hash> tt_to_node;
  Tt_t tt;
  int cnt=0;
  using node = typename bound_network::node;
  bound_network ntk( gates );
  std::array<typename bound_network::node, num_vars> pis;
  for ( int i = 0; i < num_vars; ++i )
  {
    pis[i] = ntk.get_node( ntk.create_pi() );
    kitty::create_nth_var( node_to_tt[pis[i]], i );
  }

  uint32_t line  = 1;
  do {
    auto res = database.get_entry( tt );
    if( !res )
    {
      kitty::print_binary( tt );
      std::cout << std::endl;
    }
    else
    {
      uint32_t index = std::get<1>( *res );
      int pp = std::min( p, static_cast<int>(database[index].size() - 1) );
      auto list = database[index][pp];

      auto f = list.insert( ntk, pis );
      auto n = ntk.get_node( f );
      ntk.create_po( f );
      simulate_tfi( ntk, f, node_to_tt );
      tt_to_node[node_to_tt[n]] = n;
      //list.print();
      //kitty::print_binary( node_to_tt[ntk.get_node( f )] ); std::cout << std::endl;
    }

    //if( cnt++ > 20 )
    //  break;

      kitty::next_inplace( tt );
    } while( !kitty::is_const0( tt ) );

    arrival_analyzer<bound_network> arrival( ntk, patterns[p] );
    arrival.compute_arrival_times( patterns[p] );

    std::vector<typename bound_network::signal> ordered_sigs;
    using node = typename bound_network::node;
    using signal = typename bound_network::signal;
    std::unordered_map<node, uint32_t> node_to_index;
    ntk.foreach_node( [&]( auto ni, auto ii ){
      if ( !ntk.is_constant( ni ) )
      {
        node_to_index[ni] = ordered_sigs.size();
        ordered_sigs.push_back( ntk.make_signal( ni ) );
      }
    } ); 

//    std::stable_sort( ordered_sigs.begin(), ordered_sigs.end(), [&]( auto fa, auto fb ){
//      auto na = ntk.get_node( fa );
//      auto nb = ntk.get_node( fb );
//      return arrival.get_time( na ) < arrival.get_time( nb );
//    } );
//    
    int delta = 30;
    int imin = 0;
    //for ( int imin = 0; imin < 60; imin+=delta )
    {
    for ( int ig = 0; ig < ordered_gates.size(); ++ig )
    {
      auto g = ordered_gates[ig];
      uint32_t m = g.num_vars;
      if ( m == 3 )
      {
        auto const& list = binds_o[ig].list;

        int max_i = imin+delta;
        double del_a = std::max( g.pins[0].rise_block_delay, g.pins[0].fall_block_delay );
        double del_b = std::max( g.pins[1].rise_block_delay, g.pins[1].fall_block_delay );
        double del_c = std::max( g.pins[2].rise_block_delay, g.pins[2].fall_block_delay );

        for ( int i1 = imin; i1 < max_i; ++i1 )
        {
          auto f1 = ordered_sigs[i1];
          auto n1 = ntk.get_node( f1 );
          for ( int i2 = imin; i2 < i1; ++i2 )
          {
            if ( i1 == i2 )
              continue;
            auto f2 = ordered_sigs[i2]; 
            auto n2 = ntk.get_node( f2 );
            
            for ( int i3 = imin; i3 < i2; ++i3 )
            {
              if ( (i1 == i3) || (i2 == i3) )
                continue;
              if( (i1 + i2 + i3)%100 == 0  )
                std::cout << i1 << " " << i2 << " " << i3 << std::endl;
              auto f3 = ordered_sigs[i3]; 
              auto n3 = ntk.get_node( f3 );
              
              auto na = n1;
              auto nb = n2;
              auto nc = n3;

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
              auto fa = ntk.make_signal( na );
              auto fb = ntk.make_signal( nb );
              auto fc = ntk.make_signal( nc );

              double ta = arrival.get_time( na );
              double tb = arrival.get_time( nb );
              double tc = arrival.get_time( nc );
              double t_cand = std::max( std::max( ta + del_a, tb + del_b ), tc + del_c );

              std::vector<kitty::static_truth_table<num_vars>> xs { node_to_tt[na], node_to_tt[nb], node_to_tt[nc] };
              kitty::static_truth_table<num_vars> tt;
              simulate_list( tt, list, xs );
              auto n_curr = tt_to_node[ tt ];
              double t_curr = arrival.get_time( n_curr );

              if ( t_cand < (t_curr-1) )
              {
                std::cout << i1 << " " << i2 << " " << i3 << std::endl;
                std::cout << t_curr << "->" << t_cand << std::endl;
                auto f = ntk.create_node( { fa, fb, fc }, g.name )[0];
                auto index = node_to_index[n_curr];
                ordered_sigs[index] = f; 
                ntk.substitute_node( n_curr, f );
                auto n = ntk.get_node( f );
                node_to_index[n] = index;
                node_to_tt[n] = tt;
                tt_to_node[tt] = n;
              }
            }  
          } 
        }
      }//end m == 3
      else if ( m == 4 )
      {
        auto const& list = binds_o[ig].list;

        int max_i = imin+delta;
        std::vector<std::pair<int, double>> pin_del;
        for ( int i = 0; i < 4; ++i )
          pin_del.push_back( std::make_pair( i, std::max( g.pins[i].rise_block_delay, g.pins[i].fall_block_delay ) ) );

        for ( int i1 = imin; i1 < max_i; ++i1 )
        {
          auto f1 = ordered_sigs[i1];
          auto n1 = ntk.get_node( f1 );
          for ( int i2 = imin; i2 < i1; ++i2 )
          {
            if ( i1 == i2 )
              continue;
            auto f2 = ordered_sigs[i2]; 
            auto n2 = ntk.get_node( f2 );
            
            for ( int i3 = imin; i3 < i2; ++i3 )
            {
              if ( (i1 == i3) || (i2 == i3) )
                continue;
              if( (i1 + i2 + i3)%100 == 0  )
                std::cout << i1 << " " << i2 << " " << i3 << std::endl;
              auto f3 = ordered_sigs[i3]; 
              auto n3 = ntk.get_node( f3 );
              
            for ( int i4 = imin; i4 < i3; ++i4 )
            {
              if ( (i1 == i4) || (i2 == i4) || (i3 == i4) )
                continue;
              if( (i1 + i2 + i3 + i4)%100 == 0  )
                std::cout << i1 << " " << i2 << " " << i3  << " " << i4 << std::endl;
                auto f4 = ordered_sigs[i4]; 
                auto n4 = ntk.get_node( f4 );
                auto na = n1;
                auto nb = n2;
                auto nc = n3;
                auto nd = n4;

                std::vector<node> nds { n1, n2, n3, n4 };
                std::sort( nds.begin(), nds.end(), [&]( auto a, auto b ){
                  return arrival.get_time( a ) < arrival.get_time( b );
                } );

                std::sort( pin_del.begin(), pin_del.end(), [&]( auto a, auto b ){
                  return a.second > b.second;
                } );

                na = nds[pin_del[0].first];
                nb = nds[pin_del[1].first];
                nc = nds[pin_del[2].first];
                nd = nds[pin_del[3].first];

                auto fa = ntk.make_signal( na );
                auto fb = ntk.make_signal( nb );
                auto fc = ntk.make_signal( nc );
                auto fd = ntk.make_signal( nd );

                double ta = arrival.get_time( na );
                double tb = arrival.get_time( nb );
                double tc = arrival.get_time( nc );
                double td = arrival.get_time( nd );
                double t_cand = std::max( ta + pin_del[0].second, tb + pin_del[1].second );
                t_cand = std::max( t_cand, tc + pin_del[2].second );
                t_cand = std::max( t_cand, td + pin_del[3].second );

                std::vector<kitty::static_truth_table<num_vars>> xs { node_to_tt[na], node_to_tt[nb], node_to_tt[nc], node_to_tt[nd] };
                kitty::static_truth_table<num_vars> tt;
                simulate_list( tt, list, xs );
                auto n_curr = tt_to_node[ tt ];
                double t_curr = arrival.get_time( n_curr );

                if ( t_cand < (t_curr-1) )
                {
                  std::cout << i1 << " " << i2 << " " << i3 << " " << i4 << std::endl;
                  std::cout << t_curr << "->" << t_cand << std::endl;
                  auto f = ntk.create_node( { fa, fb, fc, fd }, g.name )[0];
                  auto index = node_to_index[n_curr];
                  ordered_sigs[index] = f; 
                  ntk.substitute_node( n_curr, f );
                  auto n = ntk.get_node( f );
                  node_to_index[n] = index;
                  node_to_tt[n] = tt;
                  tt_to_node[tt] = n;
                }
              }
            }  
          } 
        }
      }
      else if ( m == 5 )
      {
        auto const& list = binds_o[ig].list;

        int max_i = imin+delta;
        std::vector<std::pair<int, double>> pin_del;
        for ( int i = 0; i < 5; ++i )
          pin_del.push_back( std::make_pair( i, std::max( g.pins[i].rise_block_delay, g.pins[i].fall_block_delay ) ) );

        for ( int i1 = imin; i1 < max_i; ++i1 )
        {
          auto f1 = ordered_sigs[i1];
          auto n1 = ntk.get_node( f1 );
          for ( int i2 = imin; i2 < i1; ++i2 )
          {
            if ( i1 == i2 )
              continue;
            auto f2 = ordered_sigs[i2]; 
            auto n2 = ntk.get_node( f2 );
            
            for ( int i3 = imin; i3 < i2; ++i3 )
            {
              if ( (i1 == i3) || (i2 == i3) )
                continue;
              if( (i1 + i2 + i3)%100 == 0  )
                std::cout << i1 << " " << i2 << " " << i3 << std::endl;
              auto f3 = ordered_sigs[i3]; 
              auto n3 = ntk.get_node( f3 );
              
            for ( int i4 = imin; i4 < i3; ++i4 )
            {
              if ( (i1 == i4) || (i2 == i4) || (i3 == i4) )
                continue;
              if( (i1 + i2 + i3 + i4)%100 == 0  )
                std::cout << i1 << " " << i2 << " " << i3  << " " << i4 << std::endl;
                auto f4 = ordered_sigs[i4]; 
                auto n4 = ntk.get_node( f4 );

                for ( int i5 = imin; i5 < i4; ++i5 )
                {
                  if ( (i1 == i5) || (i2 == i5) || (i3 == i5) || (i4 == i5) )
                    continue;
                  if( (i1 + i2 + i3 + i4 + i5)%100 == 0  )
                    std::cout << i1 << " " << i2 << " " << i3  << " " << i4 << " " << i5 << std::endl;
                    auto f5 = ordered_sigs[i5]; 
                    auto n5 = ntk.get_node( f5 );
                    auto na = n1;
                    auto nb = n2;
                    auto nc = n3;
                    auto nd = n4;
                    auto ne = n5;

                    std::vector<node> nds { n1, n2, n3, n4, n5 };
                    std::sort( nds.begin(), nds.end(), [&]( auto a, auto b ){
                      return arrival.get_time( a ) < arrival.get_time( b );
                    } );

                    std::sort( pin_del.begin(), pin_del.end(), [&]( auto a, auto b ){
                      return a.second > b.second;
                    } );

                    na = nds[pin_del[0].first];
                    nb = nds[pin_del[1].first];
                    nc = nds[pin_del[2].first];
                    nd = nds[pin_del[3].first];
                    ne = nds[pin_del[4].first];

                    auto fa = ntk.make_signal( na );
                    auto fb = ntk.make_signal( nb );
                    auto fc = ntk.make_signal( nc );
                    auto fd = ntk.make_signal( nd );
                    auto fe = ntk.make_signal( ne );

                    double ta = arrival.get_time( na );
                    double tb = arrival.get_time( nb );
                    double tc = arrival.get_time( nc );
                    double td = arrival.get_time( nd );
                    double te = arrival.get_time( ne );

                    double t_cand = std::max( ta + pin_del[0].second, tb + pin_del[1].second );
                    t_cand = std::max( t_cand, tc + pin_del[2].second );
                    t_cand = std::max( t_cand, td + pin_del[3].second );
                    t_cand = std::max( t_cand, te + pin_del[4].second );

                    std::vector<kitty::static_truth_table<num_vars>> xs { node_to_tt[na], node_to_tt[nb], node_to_tt[nc], node_to_tt[nd], node_to_tt[ne] };
                    kitty::static_truth_table<num_vars> tt;
                    simulate_list( tt, list, xs );
                    auto n_curr = tt_to_node[ tt ];
                    double t_curr = arrival.get_time( n_curr );

                    if ( t_cand < (t_curr-1) )
                    {
                      std::cout << i1 << " " << i2 << " " << i3 << " " << i4 << " " << i5 << std::endl;
                      std::cout << t_curr << "->" << t_cand << std::endl;
                      auto f = ntk.create_node( { fa, fb, fc, fd, fe }, g.name )[0];
                      auto index = node_to_index[n_curr];
                      ordered_sigs[index] = f; 
                      ntk.substitute_node( n_curr, f );
                      auto n = ntk.get_node( f );
                      node_to_index[n] = index;
                      node_to_tt[n] = tt;
                      tt_to_node[tt] = n;
                    }
                  }

              }
            }  
          } 
        }
      }
      else if ( m == 6 )
      {
        auto const& list = binds_o[ig].list;

        int max_i = imin+delta;
        std::vector<std::pair<int, double>> pin_del;
        for ( int i = 0; i < 6; ++i )
          pin_del.push_back( std::make_pair( i, std::max( g.pins[i].rise_block_delay, g.pins[i].fall_block_delay ) ) );

        for ( int i1 = imin; i1 < max_i; ++i1 )
        {
          auto f1 = ordered_sigs[i1];
          auto n1 = ntk.get_node( f1 );
          for ( int i2 = imin; i2 < i1; ++i2 )
          {
            if ( i1 == i2 )
              continue;
            auto f2 = ordered_sigs[i2]; 
            auto n2 = ntk.get_node( f2 );
            
            for ( int i3 = imin; i3 < i2; ++i3 )
            {
              if ( (i1 == i3) || (i2 == i3) )
                continue;
              if( (i1 + i2 + i3)%100 == 0  )
                std::cout << i1 << " " << i2 << " " << i3 << std::endl;
              auto f3 = ordered_sigs[i3]; 
              auto n3 = ntk.get_node( f3 );
              
            for ( int i4 = imin; i4 < i3; ++i4 )
            {
              if ( (i1 == i4) || (i2 == i4) || (i3 == i4) )
                continue;
              if( (i1 + i2 + i3 + i4)%100 == 0  )
                std::cout << i1 << " " << i2 << " " << i3  << " " << i4 << std::endl;
                auto f4 = ordered_sigs[i4]; 
                auto n4 = ntk.get_node( f4 );

                for ( int i5 = imin; i5 < i4; ++i5 )
                {
                  if ( (i1 == i5) || (i2 == i5) || (i3 == i5) || (i4 == i5) )
                    continue;
                  if( (i1 + i2 + i3 + i4 + i5)%100 == 0  )
                    std::cout << i1 << " " << i2 << " " << i3  << " " << i4 << " " << i5 << std::endl;
                    auto f5 = ordered_sigs[i5]; 
                    auto n5 = ntk.get_node( f5 );
                    
                for ( int i6 = imin; i6 < i5; ++i6 )
                {
                  if ( (i1 == i6) || (i2 == i6) || (i3 == i6) || (i4 == i6)  || (i5 == i6) )
                    continue;
                  if( (i1 + i2 + i3 + i4 + i5 + i6)%100 == 0  )
                    std::cout << i1 << " " << i2 << " " << i3  << " " << i4 << " " << i5 << " " << i6 << std::endl;
                    auto f6 = ordered_sigs[i6]; 
                    auto n6 = ntk.get_node( f6 );
                    auto na = n1;
                    auto nb = n2;
                    auto nc = n3;
                    auto nd = n4;
                    auto ne = n5;
                    auto nf = n6;

                    std::vector<node> nds { n1, n2, n3, n4, n5, n6 };
                    std::sort( nds.begin(), nds.end(), [&]( auto a, auto b ){
                      return arrival.get_time( a ) < arrival.get_time( b );
                    } );

                    std::sort( pin_del.begin(), pin_del.end(), [&]( auto a, auto b ){
                      return a.second > b.second;
                    } );

                    na = nds[pin_del[0].first];
                    nb = nds[pin_del[1].first];
                    nc = nds[pin_del[2].first];
                    nd = nds[pin_del[3].first];
                    ne = nds[pin_del[4].first];
                    nf = nds[pin_del[5].first];

                    auto fa = ntk.make_signal( na );
                    auto fb = ntk.make_signal( nb );
                    auto fc = ntk.make_signal( nc );
                    auto fd = ntk.make_signal( nd );
                    auto fe = ntk.make_signal( ne );
                    auto ff = ntk.make_signal( nf );

                    double ta = arrival.get_time( na );
                    double tb = arrival.get_time( nb );
                    double tc = arrival.get_time( nc );
                    double td = arrival.get_time( nd );
                    double te = arrival.get_time( ne );
                    double tf = arrival.get_time( nf );

                    double t_cand = std::max( ta + pin_del[0].second, tb + pin_del[1].second );
                    t_cand = std::max( t_cand, tc + pin_del[2].second );
                    t_cand = std::max( t_cand, td + pin_del[3].second );
                    t_cand = std::max( t_cand, te + pin_del[4].second );
                    t_cand = std::max( t_cand, tf + pin_del[5].second );

                    std::vector<kitty::static_truth_table<num_vars>> xs { node_to_tt[na], node_to_tt[nb], node_to_tt[nc], node_to_tt[nd], node_to_tt[ne], node_to_tt[nf] };
                    kitty::static_truth_table<num_vars> tt;
                    simulate_list( tt, list, xs );
                    auto n_curr = tt_to_node[ tt ];
                    double t_curr = arrival.get_time( n_curr );

                    if ( t_cand < (t_curr-1) )
                    {
                      std::cout << i1 << " " << i2 << " " << i3 << " " << i4 << " " << i5 << " " << i6 << std::endl;
                      std::cout << t_curr << "->" << t_cand << std::endl;
                      auto f = ntk.create_node( { fa, fb, fc, fd, fe, ff }, g.name )[0];
                      auto index = node_to_index[n_curr];
                      ordered_sigs[index] = f; 
                      ntk.substitute_node( n_curr, f );
                      auto n = ntk.get_node( f );
                      node_to_index[n] = index;
                      node_to_tt[n] = tt;
                      tt_to_node[tt] = n;
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

    std::vector<double> tA_bef;
    ntk.foreach_po( [&]( auto fo ){
      auto no = ntk.get_node( fo );
      tA_bef.push_back( arrival.get_time( no ) );
    } );

    std::vector<double> tA_aft;
    ntk.foreach_po( [&]( auto fo ){
      auto no = ntk.get_node( fo );
      tA_aft.push_back( arrival.get_time( no ) );
    } );

    for ( int i = 0; i < tA_aft.size(); ++i )
    {
      if ( tA_aft[i] < tA_bef[i] )
      {
        auto f = ntk.po_at( i );
        simulate_tfi( ntk, f, node_to_tt );
        auto n = ntk.get_node( f );
        kitty::print_binary( node_to_tt[n] );
        std::cout << " " << tA_bef[i] << "->" << tA_aft[i] << std::endl;
      }
    }

    std::cout << " adding entry" << std::endl;
    ntk.foreach_po( [&]( auto fo ){
      database_out.add_entry( ntk, pis, fo );
    } );
  }
  database_out.dump( file_out );
  return 0;

}

template<uint32_t num_vars, class Database_t, class Gates_t>
void fill_with_opto( Database_t& database, Gates_t const& gates, kitty::dynamic_truth_table tt, std::vector<double> const& pattern )
{

}