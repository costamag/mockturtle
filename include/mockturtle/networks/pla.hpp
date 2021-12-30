/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2021  EPFL
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
  \file cover.hpp
  \brief single output cover logic network implementation
  \author Andrea Costamagna
*/

#pragma once

#include <iostream>
//#include <catch.hpp>

#include <sstream>
#include <string>
#include <vector>
#include <bitset>
#include <cmath>
#include <random>
#include <boost/dynamic_bitset.hpp>

#include <mockturtle/networks/klut.hpp>
#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>

namespace mockturtle
{

struct index_to_signal
{
  index_to_signal()
  {
    storage.reserve( 10000u );
  }
  void insert( uint64_t pla_index, uint64_t klut_signal )
  {
    storage[pla_index] = klut_signal;
  }

  std::unordered_map<uint64_t, uint64_t> storage;
};

  class pla_network
  {
    #pragma region Types and constructors
    using dyn_bitset = boost::dynamic_bitset<>;
    using dbs_storage = std::vector<dyn_bitset>;


    public:
      pla_network( dbs_storage input_nodes, dbs_storage output_nodes, uint32_t max_act )
      : _nodes( input_nodes ),
        _outputs( output_nodes ),
        _num_nodes( input_nodes.at(0).size() - 1 ),
        _num_outputs( output_nodes.at(0).size() ),
        _num_data( input_nodes.size() ),
        _max_act( max_act )
        {
          _init();
        }

    protected:
      inline void _init()
      {

        for ( uint32_t i {1u}; i < _num_nodes ; ++i )
        {
          auto pi = klut.create_pi();
          _itos.insert( i, pi );
        }
        _act = 0;

      }
    #pragma endregion
    
    #pragma region visual
    public:
      void print_pla()
      {
        for ( uint32_t i {0u}; i < _num_data; ++i )
            std::cout << _outputs.at(i) << ":" << _nodes.at(i)  << std::endl ;
      }

      void print_probabilities( std::vector<double> probabilities )
      {
        uint32_t num_vars = probabilities.size();
        std::cout << std::endl;
        for ( uint32_t mask {0u}; mask < num_vars; ++mask )
        {
          dyn_bitset BS ( (uint32_t)log2(num_vars), mask );
          std::cout << "|P(" << BS << ") = " << probabilities.at(mask) << std::endl ;
        }
        std::cout << std::endl;

      }

      void print_active_list()
      {
        std::cout << "\nactive list:";
        for(uint32_t k{0u}; k<_active_list.size(); ++k)
        {
          std::cout << _active_list[k] << " ";
        }
        std::cout << "\n";
      }
    #pragma endregion
      
    #pragma region basic_function
    std::vector<double> Pr( std::vector<uint64_t> indeces_nodes, std::vector<uint64_t> indeces_outputs )
    {
      uint32_t size_P_space = std::pow( 2, indeces_nodes.size() + indeces_outputs.size() );
      std::vector<double> probabilities;
      double_t proba;
      double eq_flag_nodes, eq_flag_outputs;
      dyn_bitset b1_nodes ( _num_nodes+1, 1u );
      dyn_bitset b1_outputs ( _num_outputs, 1u );


      for ( uint64_t x_u64_t {0u}; x_u64_t < size_P_space; ++x_u64_t )
      {
        dyn_bitset xin ( _num_nodes+1, x_u64_t );
        dyn_bitset mask_nodes ( _num_nodes+1, 0u );
        dyn_bitset X_nodes ( _num_nodes+1, 0u );
        uint32_t jeff;
        for ( uint32_t j {0u}; j < indeces_nodes.size(); ++j )
        {
          jeff = indeces_outputs.size() + j;
          mask_nodes |= ( b1_nodes << indeces_nodes.at(j) );
          X_nodes |= ( ( ( ( b1_nodes << jeff ) & xin ) >> jeff ) << indeces_nodes.at(j) );
        }

        dyn_bitset mask_outputs ( _num_outputs, 0u );
        uint64_t u64_t_X_outputs {0u};

        for ( uint32_t j {0u}; j < indeces_outputs.size(); ++j )
        {
          mask_outputs |= ( b1_outputs << indeces_outputs.at(j) );
          u64_t_X_outputs |= ( ( ( ( 1u << j ) & x_u64_t ) >> j ) << indeces_outputs.at(j) );
        }
        dyn_bitset X_outputs( _num_outputs, u64_t_X_outputs );

        proba = 0;

        for ( uint32_t i {0u}; i < _num_data; ++i )
        {

          if ( ( indeces_nodes.size() != 0 ) && ( indeces_outputs.size() != 0 ) )
          {
            eq_flag_nodes = ( X_nodes == ( mask_nodes & _nodes.at(i) ) ) ? 1 : 0;
            eq_flag_outputs = ( X_outputs == ( mask_outputs & _outputs.at(i) ) ) ? 1 : 0;
          }
          else if ( indeces_nodes.size() == 0 )
          {
            eq_flag_nodes = 1;
            eq_flag_outputs = ( X_outputs == ( mask_outputs & _outputs.at(i) ) ) ? 1 : 0;

          }
          else if ( indeces_outputs.size() == 0 )
          {
            eq_flag_nodes = ( X_nodes == ( mask_nodes & _nodes.at(i) ) ) ? 1 : 0;
            eq_flag_outputs = 1;
          }
          proba += eq_flag_outputs*eq_flag_nodes/_num_data;

        }
        probabilities.push_back( proba );
      }

      return probabilities;

    }

    double H( std::vector<uint64_t> indeces_nodes, std::vector<uint64_t> indeces_outputs )
    {
      auto proba = Pr( indeces_nodes, indeces_outputs );
      uint32_t size_P_space = proba.size();
      double entropy { 0 };
      double deltaH { 0 };

      for ( uint64_t xin {0u}; xin < size_P_space; ++xin )
      {
        deltaH = ( proba[xin] == 0 ) ? 0 : -1*proba[xin]*log2( proba[xin] );
        entropy += deltaH;
      } 
      return entropy;
    }

    double MI ( std::vector<uint64_t> Xindeces, std::vector<uint64_t> Yindeces )
    {
      auto Hx = H( Xindeces, {} );
      auto Hy = H( {}, Yindeces );
      auto Hxy = H( Xindeces, Yindeces ); 

      return ( Hx + Hy - Hxy );
    }

    #pragma endregion

    #pragma region new_node
    void fill_active_list( )
    {
      double mi_loc;
      double mi_max = 0;
      uint32_t idx;
      /* first active variable */
      for( uint32_t i {0u}; i < _num_nodes; ++i )
      {
        mi_loc = MI( {i}, {0} );
        uint32_t idx; 

        if ( mi_loc >= mi_max )
        {
          mi_max = mi_loc;
          idx = i;
        }
        _active_list = {idx};
      }

      std::vector<uint64_t> inv_indeces;

      for ( uint32_t i {1u}; i < _max_act; ++i )
      {
        mi_max = 0;
        inv_indeces = _active_list;
        inv_indeces.emplace_back(0);
        for ( uint32_t j {0u}; j < _num_nodes; ++j )
        { 

          if ( std::find( _active_list.begin(), _active_list.end(), j ) != _active_list.end() )
            continue;
          
          inv_indeces.at(i) = j;

          mi_loc = MI( inv_indeces, {0} );
          if ( mi_loc >= mi_max )
          {
            mi_max = mi_loc;
            idx = j;
          }
        }
      _active_list.push_back(idx);

      }
    }

    std::string create_fn( std::vector<uint32_t> support )
    {

      uint32_t nin_node = support.size();
      std::cout << "supp size = " << nin_node << std::endl;
      uint32_t domain_size = pow( 2, nin_node );
      uint32_t Ci0, Ci1;
      dyn_bitset mask ( _num_nodes + 1, 0u );
      dyn_bitset X ( _num_nodes + 1, 0u );
      dyn_bitset Bit1 ( _num_nodes + 1, 1u );
      dyn_bitset Bit0 ( _num_nodes + 1, 0u );

      dyn_bitset Bit1_outputs ( _num_outputs, 1u );

      std::default_random_engine generator;
      std::bernoulli_distribution distribution(0.5);
      std::string tt_str;
      
      auto mask0 = ~( Bit1 << _num_nodes );

      for ( uint32_t j {0u}; j < _num_data; ++j )
      {
        _nodes.at(j) &= mask0; 
      }

      for ( uint32_t x_u64_t {0u}; x_u64_t < domain_size; ++x_u64_t )
      {
        dyn_bitset xin ( _num_nodes+1, x_u64_t );
        Ci0 = 0;
        Ci1 = 0;
        mask = Bit0;
        X = Bit0;

        for ( uint32_t j {0u}; j < support.size(); ++j )
        {
          mask |= ( Bit1 << support.at(j) );
          X |= ( ( ( ( Bit1 << j ) & xin ) >> j ) << support.at(j) );
        }

        for ( uint32_t j {0u}; j < _num_data; ++j )
        {
          if ( X == ( mask & _nodes.at(j) ) )
            ( ( _outputs.at(j) & Bit1_outputs ) == Bit1_outputs ) ? Ci1++ : Ci0++;
        }

        auto new_val = Bit0;
        if( Ci1 > Ci0 )
        {
          new_val = Bit1 << ( _num_nodes );
          tt_str = "1" + tt_str;
        }
        else if( Ci1 == Ci0 )
        {
          if (distribution(generator))
          {
            new_val = Bit1 << ( _num_nodes );
            tt_str = "1" + tt_str;
          }
          else
          {
            tt_str = "0" + tt_str;
          }
        }
        else
        {
          tt_str = "0" + tt_str;
        }

        for ( uint32_t j {0u}; j < _num_data; ++j )
        {
          if ( X == ( mask & _nodes.at(j) ) )
            _nodes.at(j) |= new_val;
        }
      }
      
      return tt_str;

    }

    void create_klut_node( std::vector<uint32_t> support, std::string tt_str )
    {
      auto supp_size = support.size();
      kitty::dynamic_truth_table tt( supp_size );
      create_from_binary_string( tt, tt_str );
      std::vector<uint64_t> klut_signals;
      for ( uint32_t i {0u}; i < supp_size; ++i )
        {
          klut_signals.push_back( _itos.storage[support[i]] );
        }
      auto f0 = klut.create_node( klut_signals, tt );
      _itos.insert( _num_nodes ,f0 );
      _num_nodes++;

      dbs_storage dbs_nodes;
      dbs_storage dbs_outputs;
      
      print_pla();

      for ( uint32_t k {0u}; k < _num_data; ++k )
      {
        _nodes.at(k).push_back(0);
      }


    }

    bool improve_fn( uint32_t nact )
    {
      std::vector<uint32_t> support;
      std::string tt_str;
      

      fill_active_list( );
      print_active_list( );
      support = {};
      std::cout << "support: ";
      for ( uint32_t k{0u}; k < nact; ++k )
      {
        std::cout << _active_list.at( _act + k ) << " ";
        support.push_back( _active_list.at( _act + k ));
      }
      std::cout << std::endl;
      //
      std::vector<uint64_t> first_act;
      std::cout << "a ";
      for (uint32_t k {0u}; k <= _act; ++k )
      {
        std::cout << k << std::endl;
        first_act.push_back(_active_list.at(k));
      }
      //
      //auto mi_old = MI( {support.at( 0 )}, {0} );
      auto mi_old = MI( first_act, {0} );
      
      tt_str = create_fn( support );
      std::cout << "truth table: " <<tt_str << std::endl;
      std::cout << _act << std::endl;
      first_act.at(_act) = _num_nodes;
     // auto mi_new = MI( {_num_nodes}, {0} );
     auto mi_new = MI( first_act , {0} );
      
      std::cout << "mi_new " << mi_new << std::endl;
      std::cout << "mi_old " << mi_old << std::endl;

      if ( mi_new > mi_old )
      {
        std::cout << "new node created. Stored at: " << _num_nodes << std::endl;
        create_klut_node( support, tt_str );
        return true;
      }
      std::cout << "Node is not kept. Remove new function " << std::endl;
      return false;
    }

    void muesli( uint32_t nact )
    {
      double eps_I_H = MI( {_num_nodes}, {0} )/H({},{0});
      double eps_th = 0.99;
      bool success;
      std::string tt_str;
      std::vector<uint32_t> support;


      while( ( eps_I_H < eps_th ) && ( _sup < _max_sup ) ) 
      {
        _act = 0;
        success = false;
        while( ( success == false ) && ( _act < _max_act ) )
        {
          success = improve_fn(2);
          if (success)
          {
            eps_I_H =  MI( {_num_nodes-1}, {0} )/H({},{0});
          }
          else
          {
            eps_I_H =  MI( {_num_nodes}, {0} )/H({},{0});
            std::cout << "true: mi(f^;f)/H(f)=" << eps_I_H << std::endl;
            _act++;
          }
        }
        while( success == true )
        {
          success = improve_fn(2);
          if (success)
            eps_I_H =  MI( {_num_nodes-1}, {0} )/H({},{0});
          else
            eps_I_H =  MI( {_num_nodes}, {0} )/H({},{0});

          std::cout << "still true: mi(f^;f)/H(f)=" << eps_I_H << std::endl;

        }
      }
      fill_active_list( );
      auto f0 = klut.create_po(_itos.storage[_active_list[0]]);
    }
    #pragma endregion


    public:
      dbs_storage _nodes;   /* storage element: value of the output at each example */
      dbs_storage _outputs; /* storage element: value of the output at each example */
      uint32_t _num_data;   /* number of examples */
      uint32_t _num_nodes;
      uint32_t _num_outputs;
      klut_network klut;
      std::vector<uint64_t> _active_list;
      index_to_signal _itos;
      uint32_t _act;
      uint32_t _sup = 0;
      uint32_t _max_act;
      uint32_t _max_sup = 2;


  };

} // namespace mockturtle