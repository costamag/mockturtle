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

/*!
  \file mapped_database.hpp
  \brief Implements utilities to store mapped sub-networks for resynthesis

  \author Andrea Costamagna
*/

#pragma once

#include <array>
#include <cassert>
#include <unordered_map>
#include <vector>

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/hash.hpp>
#include <kitty/npn.hpp>
#include <kitty/operators.hpp>
#include <kitty/print.hpp>
#include <kitty/static_truth_table.hpp>

#include <parallel_hashmap/phmap.h>

#include "../algorithms/emap.hpp"
#include "../algorithms/node_resynthesis/xag_npn.hpp"
#include "../io/genlib_reader.hpp"
#include "../io/super_reader.hpp"
#include "../networks/aig.hpp"
#include "../networks/block.hpp"
#include "../networks/klut.hpp"
#include "../utils/tech_library.hpp"
#include "../views/binding_view.hpp"
#include "../algorithms/cut_enumeration.hpp"
#include "../algorithms/cut_enumeration/rewrite_cut.hpp"
#include "include/supergate.hpp"
#include "index_list.hpp"
#include "standard_cell.hpp"
#include "struct_library.hpp"
#include "super_utils.hpp"

namespace mockturtle
{

struct mapped_database_params
{
  mapped_database_params()
  {
    /* 0 < Cut limit < 16 */
    cut_enumeration_ps.cut_limit = 8;
    cut_enumeration_ps.minimize_truth_table = true;
  }

  /*! \brief Load large gates with more than 6 inputs */
  emap_params mps;

  uint32_t num_vars{4};

  /*! \brief Cut enumeration parameters. */
  cut_enumeration_params cut_enumeration_ps{};

};

struct mapped_database_entry_t
{
  mapped_database_entry_t( mapped_index_list const& list, std::vector<gate> const& gates ):
  list( list )
  {
    area = list.compute_area( gates );
  }

  mapped_index_list list;
  double area{0.0};
};

struct mapped_database_row_t
{
  mapped_database_row_t()
  {
  }

  mapped_database_row_t( mapped_database_entry_t const& entry )
  {
    entries = {entry};
  }

  void insert( mapped_database_entry_t const& entry )
  {
    entries.push_back( entry );
  }

  std::vector<mapped_database_entry_t> entries;
};

/*! \brief Library of gates for Boolean matching
 *
 */
template<class TechLib>
class mapped_database
{
private:
  using TT = kitty::dynamic_truth_table;
  using tt_hash = kitty::hash<TT>;

public:

  mapped_database( std::string const& path, TechLib & lib, mapped_database_params const ps = {} ) :
  _tech_lib( lib ),
  _ps( ps )
  {
    std::ifstream ptr_file( path ); 
  
    /* Check if the file is present */
    if ( !ptr_file.is_open() ) 
    { 
      generate_database( path );
      return;
    }
    
    printf( "[w] database already exists, specifications ignored\n" );

    read_database( ptr_file );
    ptr_file.close(); 
  }

  mapped_database( std::string const& path, TechLib & lib, std::vector<TT> const& funcs, mapped_database_params const ps = {} ) :
  _tech_lib( lib ),
  _ps( ps )
  {
    std::ifstream ptr_file( path ); 
  
    /* Check if the file is present */
    if ( !ptr_file.is_open() ) 
    { 
      generate_database( path, funcs );
      return;
    }
    
    printf( "[w] database already exists, specifications ignored\n" );

    read_database( ptr_file );
    ptr_file.close(); 
  }

void generate_database( std::string const& path )
{
  auto resyn = xag_npn_resynthesis<aig_network, aig_network, xag_npn_db_kind::aig_complete>{};

  std::unordered_set<TT, tt_hash> classes;
  TT tt( _ps.num_vars );
  do
  {
    const auto res = kitty::exact_p_canonization( tt );
    classes.insert( std::get<0>( res ) );
    kitty::next_inplace( tt );
  } while ( !kitty::is_const0( tt ) );

  std::ofstream pfile;
  pfile.open ( path );

  int i{0};
  for ( auto const& tt : classes )
  {
    aig_network aig;
    std::vector<aig_network::signal> pis( _ps.num_vars );

    std::generate( std::begin( pis ), std::end( pis ),
      [&]() { return aig.create_pi(); } );

    std::vector<uint32_t> vars( _ps.num_vars );
    std::iota( vars.begin(), vars.end(), 0u );

  //  aig.create_po( shannon_decomposition( aig, tt, vars, pis, resyn ) );
    aig_network::signal out_sig;
    resyn( aig, tt, pis.begin(), pis.end(), [&]( auto const& f ) {
            out_sig = f;
            return false;
          } ); 
    aig.create_po( out_sig );

    binding_view<klut_network> luts = emap_klut( aig, _tech_lib, _ps.mps );

    mapped_index_list list;
    encode( list, luts );

    mapped_database_entry_t entry( list, _tech_lib.get_gates() );
    mapped_database_row_t data( entry );
    insert( tt, data );

    pfile << kitty::to_binary( tt ) << " ";
    pfile << "{ " << to_raw_index_list_string( list ) << "}\n";
  }
  pfile.close();
}

void generate_database( std::string const& path, std::vector<TT> const& funcs )
{
  auto resyn = xag_npn_resynthesis<aig_network, aig_network, xag_npn_db_kind::aig_complete>{};

  std::ofstream pfile;
  pfile.open ( path );

  int i{0};
  for ( auto const& tt : funcs )
  {
    aig_network aig;
    std::vector<aig_network::signal> pis( tt.num_vars() );

    std::generate( std::begin( pis ), std::end( pis ),
      [&]() { return aig.create_pi(); } );

    std::vector<uint32_t> vars( tt.num_vars() );
    std::iota( vars.begin(), vars.end(), 0u );

    //shannon_decomposition( aig, tt, vars, pis, resyn )
    aig_network::signal out_sig;
    resyn( aig, tt, pis.begin(), pis.end(), [&]( auto const& f ) {
            out_sig = f;
            return false;
          } ); 
    aig.create_po( out_sig );

    binding_view<klut_network> luts = emap_klut( aig, _tech_lib, _ps.mps );

    mapped_index_list list;
    encode( list, luts );
    mapped_database_entry_t entry( list, _tech_lib.get_gates() );
    mapped_database_row_t data( entry );
    insert( tt, data );

    pfile << kitty::to_binary( tt ) << " ";
    pfile << "{ " << to_raw_index_list_string( list ) << "}\n";
  }
  pfile.close();
}

void read_database( std::ifstream & pfile )
{
  std::vector<gate> const& gates = _tech_lib.get_gates();
  std::string line; 
  while ( std::getline( pfile, line ) ) 
  { 
    std::istringstream stream( line ); 
    std::vector<std::string> tokens; 
    std::string token; 
    std::vector<unsigned int> raw;
    bool first{ true };
    kitty::dynamic_truth_table tt;
    mapped_database_row_t data;
    while ( std::getline( stream, token, ' ' )) 
    { 
      if( first )
      {
        kitty::dynamic_truth_table func( log2( token.length() ) );
        kitty::create_from_binary_string( func, token );
        tt = func;
        first = false;
      }
      else
      {
        if( token == "{" )
        {
          raw.clear();
        }
        else if( token == "}" )
        {
          mapped_index_list list( raw );
          mapped_database_entry_t entry( list, _tech_lib.get_gates() );
          data.insert( entry );
        }
        else
        {
          raw.push_back( std::stoul( token ) );
        }
      }
    } 

    insert( tt, data );
  } 
}
 
template<class Ntk>
void refine_database( Ntk & ntk )
{
  static constexpr uint32_t num_vars = 4u;
  using network_cuts_t = dynamic_network_cuts<Ntk, num_vars, true, cut_enumeration_rewrite_cut>;
  using cut_manager_t = detail::dynamic_cut_enumeration_impl<Ntk, num_vars, true, cut_enumeration_rewrite_cut>;
  using cut_t = typename network_cuts_t::cut_t;
  using node_data = typename Ntk::storage::element_type::node_type;

  cut_enumeration_stats cst;
  network_cuts_t cuts( ntk.size() + ( ntk.size() >> 1 ) );
  cut_manager_t cut_manager( ntk, _ps.cut_enumeration_ps, cst, cuts );
  cut_manager.init_cuts();

  ntk.foreach_gate( [&]( auto const& n, auto i ) {
    cut_manager.clear_cuts( n );
    cut_manager.compute_cuts( n );

    for ( auto& cut : cuts.cuts( ntk.node_to_index( n ) ) )
    {
      /* skip trivial cut */
      if ( ( cut->size() == 1 && *cut->begin() == ntk.node_to_index( n ) ) )
      {
        continue;
      }

      std::vector<typename Ntk::node> leaves( *cut->size() );
      leaves.insert( leaves.begin(), *cut->begin(), *cut->end() );

      mapped_index_list list( ntk, n, leaves );
      double area_new = list.compute_area( _tech_lib.get_gates() );

      auto config = kitty::exact_p_canonization( cuts.truth_table( *cut ) );
      auto repr = std::get<0>( config );
      auto perm = std::get<2>( config );

      kitty::dynamic_truth_table tt( leaves.size() );
      kitty::create_from_binary_string( tt, kitty::to_binary( config ) );

      std::vector<uint32_t> perm2 ( leaves.size() );
      std::iota( perm2.begin(), perm2.end(), 0u );

      auto [perm3, index] = get_matches( tt, perm2 );
      printf( "%f %f\n", area_new, _data[index].area );
      
    }
  } );
}

/*! \brief Get matches.
 *
 * Given a truth table and the input data it 
 * - Extracts the representative
 * - Put in order the inputs to match the representative
 * - return a vector of matches
 * Future implementations will consider multiple orders and don't cares
 */
template<typename Input>
std::optional<std::pair<std::vector<Input>, uint32_t>> get_matches( TT const& tt, std::vector<Input> const& inputs )
{
  assert( tt.num_vars() == inputs.size() );

  const auto canon = kitty::exact_p_canonization( tt );
  auto func = std::get<0>( canon );
  auto perm = std::get<2>( canon );

  /* first match, no permutatons considered */
  std::vector<Input> ordered_inputs( inputs.size() );
  std::transform( perm.begin(), perm.end(), ordered_inputs.begin(), [&]( auto p ) { return inputs[p]; } );

  return std::make_pair( ordered_inputs, _tt_to_index[func] );
}

void insert( TT const& tt, mapped_database_row_t const& data )
{
  _tt_to_index[tt] = _data.size();
  _data.push_back( data );
}

template<typename Fn>
void foreach_entry( uint32_t index, Fn&& fn ) const
{
  for ( int i{0}; i<_data[index].entries.size(); ++i )
  {
    fn( _data[index].entries[i] );
  }
}

mapped_index_list get_list( uint32_t row, uint32_t column )
{
  return _data[row].entries[column].list;
}

private:
  phmap::flat_hash_map<TT, uint32_t, kitty::hash<TT>> _tt_to_index;
  std::vector<mapped_database_row_t> _data;
  mapped_database_params _ps;
  TechLib _tech_lib;

}; /* class exact_library */

} // namespace mockturtle
