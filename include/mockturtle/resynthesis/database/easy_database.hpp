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
  \file mapped_list.hpp
  \brief Index list for mapped networks

  \author Andrea Costamagna
*/

#pragma once

#include <fstream>
#include <vector>

#include <kitty/hash.hpp>
#include <kitty/constructors.hpp>

#include <parallel_hashmap/phmap.h>

#include "easy_storage.hpp"
#include "mapped_list.hpp"
#include "../list_simulator.hpp"
#include "../mapped_list_simulator.hpp"
#include "../networks/bindings.hpp"

namespace mockturtle
{
  /*! \brief Homogeneous database 
   * Homogeneous database, all the functionalities are represented
   * using static truth tables with num_vars inputs
   */
  template<uint32_t num_vars, uint32_t capacity = 1>
  class easy_database
  {
    /* each database entry  */
    using TT = kitty::static_truth_table<num_vars>;
    using tt_hash = kitty::hash<TT>;
    using simulator_t = mapped_list_simulator<bindings, TT>;
    using storage_t = easy_storage<num_vars, capacity>;
    using entry_t = typename storage_t::entry_t;
  
  public:

    easy_database( bindings const& gates, bool init = false )
      :  storage( std::make_shared<storage_t>() ),
        proj_fns(),
        simulator( gates ),
        lib( gates )
    {
      /* projection functions for functional simulation */
      proj_fns.reserve( num_vars );
      for( int i{0}; i < num_vars; ++i )
      {
        proj_fns.emplace_back();
        kitty::create_nth_var( proj_fns[i], i );
      }
      /* Initialize the database with the gates from the library */
      if ( init )
        _init();
    }

    std::vector<std::string> split( std::string str, char* delimiter )
    {
      std::vector<std::string> v;
      char *token = strtok( const_cast<char*>(str.c_str()), delimiter );
      while ( token != nullptr )
      {
        v.push_back( std::string( token ) );
        token = strtok( nullptr, delimiter );
      }
      
      return v;
    }

    easy_database( bindings const& gates, std::string const& path )
      :  lib( gates ), simulator( gates ), storage(std::make_shared<storage_t>())
    {
      /* projection functions for functional simulation */
      for( int i{0}; i < num_vars; ++i )
      {
        proj_fns.emplace_back();
        kitty::create_nth_var( proj_fns[i], i );
      }

      std::fstream file;
      file.open( path );
 
      mapped_list::element_type v;
      if ( file.is_open() ) {
        std::string line;
        while ( getline( file, line ) ) {
//std::cout << line << std::endl;
          std::vector<std::string> vstr = split( line, " " );
        
          mapped_list list( num_vars );
          uint32_t index = 0;
          uint32_t const num_pis = std::stoi( vstr[index++] );
          uint32_t const num_gates = std::stoi( vstr[index++] );
          uint32_t const num_pos = std::stoi( vstr[index++] );
          assert( num_pis == num_vars );
          double area = 0;
          for( int i = 0; i < num_gates; ++i )
          {
            uint32_t const fanin_size = std::stoi( vstr[index++] );
            std::vector<uint32_t> children( fanin_size );
            for( int j = 0; j < fanin_size; ++j )
            {
              children[j] = std::stoi( vstr[index++] );
            }
            uint32_t const id = std::stoi( vstr[index++] );
            v = list.create_node( children, id );
            area += lib[id].area;
          }
          
          list.create_po( std::stoi( vstr[index++] ) );
          list.compute_pin_delays( lib );
          list.set_area( area );
          add_entry( list );
        }

      }

      file.close();

    }

    /*! \brief Initialize the database using the gates from the library */
    void _init()
    {
      auto gates = lib.get_gates();
      for( auto const& g : gates )
      {
        if( g.num_vars > num_vars )
          continue;
        mapped_list list( num_vars );
        std::vector<mapped_list::element_type> children( g.num_vars );
        for( auto i = 0; i < g.num_vars; ++i )
        {
          children[i] = list.pi_at( i );
          std::cout << g.pins[i].input_load << std::endl;
        }
        if( children.size() > 0 )
        {
          mapped_list::element_type v = list.create_node( children, g.id );
          list.create_po( v );
          list.compute_pin_delays( lib );
          list.set_area( g.area );
        }
        else if( kitty::is_const0( g.function ) )
        {
          mapped_list::element_type v = list.get_constant( false );
          list.create_po( v );
          list.compute_pin_delays( lib );
          list.set_area( g.area );
        }
        else if( kitty::is_const0( ~g.function ) )
        {
          mapped_list::element_type v = list.get_constant( true );
          list.create_po( v );
          list.compute_pin_delays( lib );
          list.set_area( g.area );
        }
        else
        {
          std::cerr << "Unknown gate " << g.name << " no children and not constant" << std::endl;
        }
        add_entry( list );
      }
    }

    easy_database<num_vars> & get_database()
    {
      return (*this);
    }

    entry_t const& operator[]( uint32_t index ) const
    {
      return (*storage)[index];
    }

    bindings get_bindings() const
    {
      return lib;
    }

    

    /*! \brief add a database entry after p-canonization */
    template<class Ntk>
    void add_entry( Ntk & ntk, std::array<typename Ntk::node, num_vars> const& cut, typename Ntk::signal f )
    {
      /* construct the mapped list */
      mapped_list list( ntk, lib, cut, f );
      add_entry( list );
    }

    void add_entry( mapped_list & list )
    {
      TT res;
      /* get its functionaliy */
      /* get the new functionality */
      simulator( res, list, proj_fns );
      
//kitty::print_binary( res ); std::cout << std::endl;

      /* update the database */
      storage->add_entry( res, list );
    }

    /*! \brief extract the database entry matching the function */
    std::optional<std::tuple<std::vector<uint8_t>, uint32_t>> get_entry( TT const& tt )
    {
      auto match = storage->find( tt );
      if ( match == storage->end() )
      {
        return std::nullopt;
      }
      else // TODO: add the reward-based check
      {
        uint32_t index = match->second;
        std::vector<uint8_t> perm;
        for ( uint8_t i = 0; i < num_vars; ++i )
          perm.push_back( i );
        return std::tuple( perm, index );
      }
    }

    /*! \brief extract the database entry matching the function */
    std::optional<std::tuple<std::vector<uint8_t>, uint32_t>> get_entry( TT const& tt, TT const& care )
    {
      (void)care;

      auto match = storage->find( tt );
      if ( match == storage->end() )
      {
        return std::nullopt;
      }
      else // TODO: add the reward-based check
      {
        uint32_t index = match->second;
        std::vector<uint8_t> perm;
        for ( uint8_t i = 0; i < num_vars; ++i )
          perm.push_back( i );
        return std::tuple( perm, index );
      }
    }

    /*! \brief iterator to apply a lambda to all the nodes */
    template<typename Fn>
    void foreach_list( uint32_t index, Fn&& fn ) const
    {
      storage->foreach_list( index, fn );
    }

    /*! \brief iterator to apply a lambda to all the nodes */
    template<typename Fn>
    void foreach_entry( Fn&& fn ) const
    {
      storage->foreach_entry( fn );
    }

    void print()
    {
      foreach_entry( [&]( auto entry, auto ientry ){
        std::cout << "Entry " << ientry << std::endl;
        entry.foreach_list( [&]( auto list, auto ilist ){
          list.print();
        } );
      } );
    }

    void dump( std::string path )
    {
      std::ofstream file;
      file.open( path );
      foreach_entry( [&]( auto entry, auto ientry ){
        entry.foreach_list( [&]( auto list, auto ilist ){
          file << list.num_pis() << " ";
          file << list.num_gates() << " ";
          file << list.num_pos() << " ";
          list.foreach_gate( [&]( auto children, auto id ){
            file << children.size() << " ";
            for( auto child : children )
            {
              file << child << " ";
            }
            file << id << " ";
          } );
          list.foreach_po( [&]( auto lit ){
            file << lit << " ";
          } );
          file << "\n";
        } );
      } );
      file.close();
    }

    auto get_symm( uint32_t entry_index )
    {
      auto symm = storage->get_symm( entry_index );
      return symm;
    }


    /* store the nodes  */
    std::shared_ptr<storage_t> storage;
    std::vector<TT> proj_fns;
    simulator_t simulator;
    bindings lib;
  };

} // namespace mockturtle