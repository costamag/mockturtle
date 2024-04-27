/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
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

#include <iostream>
#include <string>
#include <fstream>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/genlib.hpp>
#include <mockturtle/algorithms/mapper2.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/scg.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/utils/node_map.hpp>


#include <experiments.hpp>

using namespace experiments;
using namespace mockturtle;

template<class TT>
TT compute( std::vector<TT> const& sims, TT const& function )
{
  uint32_t num_vars = sims.size();
  TT tmp(4u);
  TT sim = tmp.construct();
  for( uint32_t m{0}; m<(1u<<num_vars); ++m )
  {
    if( kitty::get_bit( function, m ) == 1 )
    {
      tmp = tmp | ~tmp;
      for( int i{0}; i<num_vars; ++i )
      {
        if( ( m >> i ) & 0x1 == 0x1 )
        {
          tmp &= sims[i];
        }
        else
        {
          tmp &= ~sims[i];
        }
      }
      sim |= tmp;
    }
  }
  return sim;
}

    struct comb_t
    {
      comb_t( int N, int K ): N(N), K(K), comb(K)
      {
        for( int j{0}; j<K; ++j )
        {
          comb[j]=j;
        }
        act=K-1;
      }

      std::optional<std::vector<uint32_t>> get()
      {
        bool finished=true;
        for( int i{0}; i<K; ++i )
        {
          if( comb[i]+K <= N+i )
          {
            finished=false;
            break;
          }
        }

        if( finished )
        {
          return std::nullopt;
        }
        auto res = comb;

        for( int j{act}; j>=0; --j )
        {
          comb[j]=comb[j]+1;
          if( comb[j]+K-1 < N+j )
          {
            for( int k{j+1}; k<=act; ++k )
            {
              comb[k]=comb[k-1]+1;
            }
            break;
          }
        }

        return res;
      }

      int N;
      int K;
      std::vector<uint32_t> comb;
      uint32_t act;
    };

template<class TT, class TTHASH>
void print_status( scg_network const& scg, std::vector<std::vector<scg_network::signal>> const& sigs_x_count, std::unordered_map<TT, scg_network::signal, TTHASH> const& existing, std::unordered_map<uint64_t, uint32_t> const& db_PClassMap, std::vector<double> const&  db_areas, std::unordered_map<uint64_t, uint32_t> const& db2_PClassMap, std::vector<double> const&  db2_areas )
{
  printf("#PIS=%3ld #NDS=%3ld\n", scg.num_pis(), scg.num_gates()); 
  for( uint32_t cnt{0}; cnt<sigs_x_count.size(); ++cnt )
  {
    printf("#sigs_x_count[%d]=%3ld\n", cnt, sigs_x_count[cnt].size() ); 
  }
  printf("|TT HASH TABLE|=%3ld\n", existing.size()); 
  for (const auto& [key, value] : existing )
  {
    kitty::print_binary(key);
    const auto res = kitty::exact_p_canonization( key );
    TT func_p = std::get<0>( res );
    uint64_t key64 = (func_p._bits)[0];
    if( db_PClassMap.find( key64 ) != db_PClassMap.end() )
      printf(" a(db):%f a(db2)%f\n", db_areas[db_PClassMap.at(key64)], db2_areas[db2_PClassMap.at(key64)]);
    
    //uint64_t pclass = db_PClassMap[ func_p._bits[0] ];
   // printf("->%d : db=%f\n", value.index, db_areas[db_PClassMap[pclass]] );
  }
}


int main()
{
  using namespace experiments;
  using namespace mockturtle;

  /* library to map to technology */
  std::vector<gate> gates;
  std::ifstream in( cell_libraries_path( "sky130" ) );

  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  tech_library_params tps;
  tech_library<5, classification_type::np_configurations> tech_lib( gates, tps );

  xag_npn_resynthesis<aig_network, aig_network, xag_npn_db_kind::aig_complete> resyn;
  exact_library_params eps;
  eps.np_classification = false;
  exact_library<aig_network> exact_lib( resyn, eps );

  /* read the database information */
  std::unordered_map<uint64_t, uint32_t> db_PClassMap;
  std::vector<double> db_areas;
  kitty::static_truth_table<4u> ttdb;
  int i{0};
  std::string line;
  std::ifstream fTts ("sky130.tts");
  if (fTts.is_open())
  {
    while ( std::getline (fTts,line) )
    {
      kitty::create_from_binary_string( ttdb, line );
      db_PClassMap[ttdb._bits]=i++;
    }
    fTts.close();
  }
  else
  {
    printf("not found\n");
  }

  std::ifstream fAreas ("sky130.area");
  if (fAreas.is_open())
  {
    while ( std::getline (fAreas,line) )
    {
      db_areas.push_back( std::stof( line ) );
    }
    fAreas.close();
  }
  else
  {
    printf("not found\n");
  }

  /* build the P-classes */
  using TT = kitty::dynamic_truth_table;
  using tt_hash = kitty::hash<TT>;

  std::unordered_set<TT, tt_hash> classes;
  TT tt(4u);
  int iter=0;
  do
  {
    const auto res = kitty::exact_p_canonization( tt );
    classes.insert( std::get<0>( res ) );
    kitty::next_inplace( tt );
    iter++;
  } while ( iter < 2 && !kitty::is_const0( tt ) );

  for( auto g : gates )
  {
    printf( "%d ", g.id );
    std::cout << g.name << " " << g.expression ;
    printf( " num_vars=%d ", g.num_vars );
    kitty::print_binary( g.function );
    printf(" area=%f \n", g.area );
  }

  std::sort( gates.begin(), gates.end(), [&]( auto g1, auto g2 ){ return g1.area < g2.area; } );

  printf("sorted\n");

  for( auto g : gates )
  {
    printf( "%d ", g.id );
    std::cout << g.name << " " << g.expression ;
    printf( " num_vars=%d ", g.num_vars );
    kitty::print_binary( g.function );
    printf(" area=%f \n", g.area );
  }

  /* the gates are sorted */
  std::unordered_map<TT, scg_network::signal, tt_hash> existing;

  scg_network scg;
  std::vector<std::vector<scg_network::signal>> sigs_x_count;
  std::unordered_map<uint64_t, uint32_t> db2_PClassMap;
  std::unordered_map<uint64_t, uint32_t> db2_node_to_idx;

  std::vector<double> db2_areas;  
  std::vector<TT> db2_tts;  
  uint64_t key64;

  /* fill in the 0 gates */
  sigs_x_count.emplace_back();
  TT tt0(4u);  
  
  std::vector<scg_network::signal> pis;
  TT tti(4u);
  for( int i{0}; i<4; ++i )
  {
    pis.push_back( scg.create_pi() );
    kitty::create_nth_var( tti, i );
    existing[tti]=pis[i];
    sigs_x_count[0].push_back( pis[i] );
    key64=(tti)._bits[0];
    db2_PClassMap[key64] = db2_areas.size();
    db2_node_to_idx[pis[i]] = db2_areas.size();
    db2_areas.push_back(0);
    db2_tts.push_back(tti);
  }

  sigs_x_count[0].push_back( scg.get_constant(false) );
  key64 = ((~tt0)._bits)[0];
  db2_PClassMap[key64] = db2_areas.size();
  db2_node_to_idx[scg.get_constant(false)] = db2_areas.size();

  db2_areas.push_back(0);
  db2_tts.push_back(~tt0);
  existing[~tt0]=(scg.get_constant(false));

  sigs_x_count[0].push_back( scg.get_constant(true) );

  key64 = ((tt0)._bits)[0];
  db2_PClassMap[key64] = db2_areas.size();
  db2_node_to_idx[scg.get_constant(true)] = db2_areas.size();
  db2_areas.push_back(0);
  db2_tts.push_back(tt0);
  existing[tt0]=(scg.get_constant(true));

  print_status( scg, sigs_x_count, existing, db_PClassMap, db_areas, db2_PClassMap, db2_areas );

  for( auto g : gates )
  {
    comb_t combs( 4, g.num_vars );
    
    int it{0};
    while( true )
    {
      auto comb = combs.get();
      if( !comb )
      {
        break;
      }
      std::vector<TT> sims;
      for( int i{0}; i<(*comb).size(); ++i )
      {
        printf("%d ", (*comb)[i] );
        sims.push_back( db2_tts[db2_node_to_idx[pis[(*comb)[i]]]] );
      }
      kitty::print_binary( g.function );
      printf("\n");
      auto sim = compute( sims, g.function );
      kitty::print_binary(sim); printf("<newsim\n");
      printf("\n");

      const auto res = kitty::exact_p_canonization( sim );
      auto ttcand = std::get<0>( res );
      key64 = ttcand._bits;
      if( db2_PClassMap.find( key64 )  ==  db2_PClassMap.end() )
      {
        db2_PClassMap[key64] = db2_areas.size();
        db2_node_to_idx[scg.get_constant(true)] = db2_areas.size();
        db2_areas.push_back( g.area );
        db2_tts.push_back( ttcand );
        
        existing[ttcand]=(scg.get_constant(true));
      }

    }
  }
  /* fill in the 1 gates */


  return 0;
}


