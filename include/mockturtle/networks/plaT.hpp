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
      pla_network( dbs_storage input_nodes, dbs_storage output_nodes, uint32_t max_act, uint32_t max_sup = 2, uint32_t init_sup = 2 )
      : _nodes( input_nodes ),
        _outputs( output_nodes ),
        _num_nodes( input_nodes.at(0).size() - 1 ),
        _num_outputs( output_nodes.at(0).size() ),
        _num_data( input_nodes.size() ),
        _max_act( max_act ),
        _max_sup( max_sup ),
        _init_sup( init_sup )
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

    #pragma region basic_function_given_data
    std::vector<double> Pr_gd( std::vector<uint64_t> indeces_nodes, std::vector<uint64_t> indeces_outputs, 
                            dbs_storage nodes, dbs_storage outputs, 
                            uint32_t num_nodes )
    {
      uint32_t size_P_space = std::pow( 2, indeces_nodes.size() + indeces_outputs.size() );
      std::vector<double> probabilities;
      double_t proba;
      double eq_flag_nodes, eq_flag_outputs;
      dyn_bitset b1_nodes ( num_nodes+1, 1u );
      dyn_bitset b1_outputs ( _num_outputs, 1u );
      auto num_data = nodes.size();


      for ( uint64_t x_u64_t {0u}; x_u64_t < size_P_space; ++x_u64_t )
      {
        dyn_bitset xin ( num_nodes+1, x_u64_t );
        dyn_bitset mask_nodes ( num_nodes+1, 0u );
        dyn_bitset X_nodes ( num_nodes+1, 0u );
        uint32_t jeff;
        //std::cout << "z1 "<<std::endl;
        for ( uint32_t j {0u}; j < indeces_nodes.size(); ++j )
        {
          jeff = indeces_outputs.size() + j;
          mask_nodes |= ( b1_nodes << indeces_nodes.at(j) );
          X_nodes |= ( ( ( ( b1_nodes << jeff ) & xin ) >> jeff ) << indeces_nodes.at(j) );
        }
        //std::cout << "z2 "<<std::endl;

        dyn_bitset mask_outputs ( _num_outputs, 0u );
        uint64_t u64_t_X_outputs {0u};

        //std::cout << "z3 "<<std::endl;
        for ( uint32_t j {0u}; j < indeces_outputs.size(); ++j )
        {
          mask_outputs |= ( b1_outputs << indeces_outputs.at(j) );
          u64_t_X_outputs |= ( ( ( ( 1u << j ) & x_u64_t ) >> j ) << indeces_outputs.at(j) );
        }
        dyn_bitset X_outputs( _num_outputs, u64_t_X_outputs );

        proba = 0;
        //std::cout << "z4 "<<std::endl;
        for ( uint32_t i {0u}; i < num_data; ++i )
        {

          if ( ( indeces_nodes.size() != 0 ) && ( indeces_outputs.size() != 0 ) )
          {
            eq_flag_nodes = ( X_nodes == ( mask_nodes & nodes.at(i) ) ) ? 1 : 0;
            eq_flag_outputs = ( X_outputs == ( mask_outputs & outputs.at(i) ) ) ? 1 : 0;
          }
          else if ( indeces_nodes.size() == 0 )
          {
            eq_flag_nodes = 1;
            eq_flag_outputs = ( X_outputs == ( mask_outputs & outputs.at(i) ) ) ? 1 : 0;

          }
          else if ( indeces_outputs.size() == 0 )
          {
            eq_flag_nodes = ( X_nodes == ( mask_nodes & nodes.at(i) ) ) ? 1 : 0;
            eq_flag_outputs = 1;
          }
          proba += eq_flag_outputs*eq_flag_nodes/num_data;

        }
        //std::cout << "z5 "<<std::endl;
        probabilities.push_back( proba );
      }

      return probabilities;

    }

    double H_gd( std::vector<uint64_t> indeces_nodes, std::vector<uint64_t> indeces_outputs, 
                            dbs_storage nodes, dbs_storage outputs, 
                            uint32_t num_nodes )
    {
      /*std::cout << "a1: " << num_nodes << std::endl;
      for( uint32_t k{0u}; k<indeces_nodes.size(); ++k )
        std::cout << indeces_nodes[k] << std::endl;
      for( uint32_t k{0u}; k<indeces_outputs.size(); ++k )
        std::cout << indeces_outputs[k] << std::endl;
      for( uint32_t k{0u}; k<nodes.size(); ++k )
        std::cout << nodes[k] << ": " << outputs[k] << std::endl;*/

      auto proba = Pr_gd( indeces_nodes, indeces_outputs, nodes, outputs, num_nodes );
      //std::cout << "a2" << std::endl;
      uint32_t size_P_space = proba.size();
      //std::cout << "a2" << std::endl;
      double entropy { 0 };
      double deltaH { 0 };

      for ( uint64_t xin {0u}; xin < size_P_space; ++xin )
      {
        deltaH = ( proba[xin] == 0 ) ? 0 : -1*proba[xin]*log2( proba[xin] );
        entropy += deltaH;
      } 
      return entropy;
    }

    double MI_gd ( std::vector<uint64_t> Xindeces, std::vector<uint64_t> Yindeces,
                   dbs_storage nodes, dbs_storage outputs, uint32_t num_nodes )
    {
      //std::cout << "a" << std::endl;
      auto Hx = H_gd( Xindeces, {}, nodes, outputs, num_nodes );
      //std::cout << "b" << std::endl;
      auto Hy = H_gd( {}, Yindeces, nodes, outputs, num_nodes );
      //std::cout << "c" << std::endl;
      auto Hxy = H_gd( Xindeces, Yindeces, nodes, outputs, num_nodes ); 

      return ( Hx + Hy - Hxy );
    }

    std::string create_fn_gd( std::vector<uint32_t> support, dbs_storage nodes, dbs_storage outputs )
    {
      uint32_t num_nodes = support.size();
      uint32_t num_data = nodes.size();
      uint32_t nin_node = support.size();
      std::cout << "supp size = " << nin_node << std::endl;
      uint32_t domain_size = pow( 2, nin_node );
      uint32_t Ci0, Ci1;
      dyn_bitset mask ( num_nodes + 1, 0u );
      dyn_bitset X ( num_nodes + 1, 0u );
      dyn_bitset Bit1 ( num_nodes + 1, 1u );
      dyn_bitset Bit0 ( num_nodes + 1, 0u );

      dyn_bitset Bit1_outputs ( _num_outputs, 1u );

      std::default_random_engine generator;
      std::bernoulli_distribution distribution(0.5);
      std::string tt_str;
      
      auto mask0 = ~( Bit1 << num_nodes );
      //std::cout << "mask " << mask0 << " domain size " << domain_size << std::endl;

      for ( uint32_t j {0u}; j < nodes.size(); ++j )
      {
        nodes.at(j) &= mask0; 
      }
      //std::cout << "nodes " << nodes[0] << std::endl;


      for ( uint32_t x_u64_t {0u}; x_u64_t < domain_size; ++x_u64_t )
      {
        dyn_bitset xin ( num_nodes+1, x_u64_t );
        Ci0 = 0;
        Ci1 = 0;
        mask = Bit0;
        X = Bit0;
        //std::cout << "pre " << std::endl;
        for ( uint32_t j {0u}; j < support.size(); ++j )
        {
          mask |= ( Bit1 << support.at(j) );
          X |= ( ( ( ( Bit1 << j ) & xin ) >> j ) << support.at(j) );
        }
        //std::cout << "mask " << mask << " X "<< X << std::endl;

        for ( uint32_t j {0u}; j < nodes.size(); ++j )
        {
          if ( X == ( mask & nodes.at(j) ) )
            ( ( outputs.at(j) & Bit1_outputs ) == Bit1_outputs ) ? Ci1++ : Ci0++;
        }
        //std::cout << "C1 " << Ci1<< " C0 " << Ci0 << std::endl; 

        auto new_val = Bit0;
        if( Ci1 > Ci0 )
        {
          new_val = Bit1 << ( num_nodes );
          tt_str = "1" + tt_str;
        }
        else if( Ci1 == Ci0 )
        {
          if (distribution(generator))
          {
            new_val = Bit1 << ( num_nodes );
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

        for ( uint32_t j {0u}; j < num_data; ++j )
        {
          if ( X == ( mask & nodes.at(j) ) )
            nodes.at(j) |= new_val;
        }
      }
      //std::cout << "str: " << tt_str << std::endl;
      
      return tt_str;

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
      //std::cout << "supp size = " << nin_node << std::endl;
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

    bool improve_fn( )
    {
      std::vector<uint32_t> support;
      std::string tt_str;
      

      fill_active_list( );
      print_active_list( );
      support = {};
      std::cout << "support: ";
      for ( uint32_t k{0u}; k < _sup; ++k )
      {
        std::cout << _active_list.at( _act + k ) << " ";
        support.push_back( _active_list.at( _act + k ));
      }
      std::cout << std::endl;
      //
      std::vector<uint64_t> first_act;
      //std::cout << "a ";
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

    bool not_done( uint32_t best_idx )
    {
      double eps_I_H = MI( {best_idx}, {0} )/H({},{0});
      if ( eps_I_H >= _eps_th )
        return false;
      else
        return true;
      
      if ( eps_I_H > _eps_best )
      {
        _eps_best = eps_I_H;
        _idx_fn = best_idx;
      }
    }

    void muesli( double eps_th = 0.99 )
    {

      bool success; /* true if found a function improving the mi */
      std::string tt_str; /* contains the tt of the new node */
      std::vector<uint32_t> support; /* contains the support of the new function */
      uint32_t best_idx; /* used to keep track of the new best approximation */

      /* nlist <- sort_nlist_by_I(nlist,1) */
      _eps_th = eps_th;
      double eps_best = 0;
      auto max_act_tmp = _max_act; /* store temporarily */
      _max_act = 1; /* set to one to speed up the search of the max mi */
      fill_active_list(); 
      _idx_fn = _active_list[0];
      best_idx = _idx_fn;
      _max_act = max_act_tmp;
      /* sup <- 2 */
      _sup = _init_sup;
      /* while ( not_done(nlist) AND sup < max_sup )  < turned into <= else max_sup 2 gives false*/
      while( not_done(best_idx) && ( _sup <= _max_sup ) ) 
      {
        _act = 0; /* select which node in the list is active */
        success = false;
        do {

          success = improve_fn();
          best_idx = _num_nodes;
          if (success)
            best_idx -= 1;
               
          std::cout << "best idx = " << best_idx << std::endl;
          std::cout << MI( {best_idx}, {0} )/H({},{0}) << std::endl;
          if ( not_done( best_idx ) == false )
            break;

          _act++;       

        } while( ( success == false ) && ( _act <= _max_act ) );
        if ( success == true )
        {
          if ( not_done( best_idx ) == false )
            break;
          _sup = _init_sup;
          while( success == true )
          {
            success = improve_fn();
            best_idx = _num_nodes;
            if (success)
              best_idx -= 1;
            std::cout << "best idx = " << best_idx << std::endl;
            std::cout << MI( {best_idx}, {0} )/H({},{0}) << std::endl;
          }
        }
        else
        {
          _sup++;
        }
      }
      
      max_act_tmp = _max_act; /* store temporarily */
      _max_act = 1; /* set to one to speed up the search of the max mi */
      fill_active_list(); 
      best_idx = _active_list[0];
      _max_act = max_act_tmp;
      
      auto f0 = klut.create_po(_itos.storage[_active_list[0]]);
    }

    void muesli_modified( double eps_th = 0.99 )
    {

      bool success; /* true if found a function improving the mi */
      std::string tt_str; /* contains the tt of the new node */
      std::vector<uint32_t> support; /* contains the support of the new function */
      uint32_t best_idx; /* used to keep track of the new best approximation */

      /* nlist <- sort_nlist_by_I(nlist,1) */
      _eps_th = eps_th;
      double eps_best = 0;
      auto max_act_tmp = _max_act; /* store temporarily */
      _max_act = 1; /* set to one to speed up the search of the max mi */
      fill_active_list(); 
      _idx_fn = _active_list[0];
      best_idx = _idx_fn;
      _max_act = max_act_tmp;
      /* sup <- 2 */
      _sup = _init_sup;
      /* while ( not_done(nlist) AND sup < max_sup )  < turned into <= else max_sup 2 gives false*/
      while( not_done(best_idx) && ( _sup <= _max_sup ) ) 
      {
        _act = 0; /* select which node in the list is active */
        success = false;
        do {
          success = improve_fn();
          best_idx = _num_nodes;
          if (success)
            best_idx -= 1;
               
          std::cout << "best idx = " << best_idx << std::endl;
          std::cout << MI( {best_idx}, {0} )/H({},{0}) << std::endl;
          if ( not_done( best_idx ) == false )
            break;

          _act++;       

        } while( ( success == false ) && ( _act <= _max_act ) );
        if ( success == true )
        {
          if ( not_done( best_idx ) == false )
            break;
          _sup = _init_sup;
        }
        else
        {
          _sup++;
        }
      }
      
      max_act_tmp = _max_act; /* store temporarily */
      _max_act = 1; /* set to one to speed up the search of the max mi */
      fill_active_list(); 
      best_idx = _active_list[0];
      _max_act = max_act_tmp;
      
      auto f0 = klut.create_po(_itos.storage[_active_list[0]]);
    }


    #pragma endregion

    #pragma region details_muesli_preprocessing
    template<typename T>
    void swap( T& a, T& b)
    {
      T t = a;
      a = b;
      b = t;
    }

    int partition ( std::vector<uint32_t>& support, std::vector<double>& attribute , uint32_t low, uint32_t high )
    {
    double pivot = attribute[high];    // pivot
    int i = (low-1);  // Index of smaller element
 
    for (int j = low; j < high; j++)
      {
        // If current element is smaller than or
        // equal to pivot
        if ( attribute[j] >= pivot)
        {
            i++;    // increment index of smaller element
            swap(attribute[i], attribute[j]);
            swap(support[i], support[j]);
        }
      }
      swap(attribute[i + 1], attribute[high]);
      swap(support[i + 1], support[high]);
      return (i + 1);
    }

    void quicksort_by_attribute( std::vector<uint32_t>& support, std::vector<double>& attribute,  int low, int high )
    {
      if (low < high)
      {
        /* pi is partitioning index, arr[p] is now
           at right place */
        auto pi = partition( support, attribute, low, high);
 
        // Separately sort elements before
        // partition and after partition
        quicksort_by_attribute( support, attribute, low, pi - 1);
        quicksort_by_attribute( support, attribute, pi + 1, high);
      }
    }

    std::vector<std::vector<uint32_t>> group_by_mi( std::vector<uint32_t> const& support, std::vector<double> const& mi_v, double dI = 0 )
    {
      std::vector<std::vector<uint32_t>> Pi;
      std::vector<double> miP;
      uint32_t idxPi = 0;
      Pi.push_back({support[0]});
      miP.push_back(mi_v[0]);
      for( uint32_t k{1u}; k<support.size(); ++k )
      {
        if( mi_v[k] >= ( miP[idxPi] - miP[idxPi]*dI) )
        {
          Pi[idxPi].push_back( support[k] );
        }
        else
        {
          idxPi++;
          Pi.push_back( {support[k]} );
          miP.push_back(mi_v[k]);
        }
      }
      return Pi;
    }

    uint32_t r_create_fn_from_support( std::vector<uint32_t> p, std::vector<uint32_t> given_klg = {}, uint64_t o_idx = 0 )
    {
      std::cout << "\n{ ";
        for( uint32_t j{0u}; j <p.size(); ++j )
        {
          std::cout << p[j] << " ";
        }
        std::cout << "}\n";
      if( given_klg.size() == 0 )
      {
        if ( p.size() == 1 )
        {
          return p[0];
        }
        else if( p.size() <= _max_sup )
        {
          auto tt_new = create_fn( p );
          create_klut_node( p, tt_new );
          return ( _num_nodes - 1 );
        }
        else // |pi| > _max_sup
        {
          auto x = p[0];
          std::vector<double> mi_v;
          std::vector<uint32_t> p1;

          for( uint32_t k{1u}; k < p.size(); ++k )
          {
            mi_v.emplace_back( MI( { p[k], x }, { o_idx }) );
            p1.emplace_back( p[k] );
          }
          auto P1 = group_by_mi( p1, mi_v, 0 );
          std::vector<uint32_t> Fns;
          std::vector<double> mi_Fns;

          Fns.push_back( r_create_fn_from_support( P1[0], {x}, 0 ) );
          mi_Fns.push_back( MI( {Fns[0]}, {o_idx} ) );
          for ( uint32_t k {1u}; k<P1.size(); ++k )
          {
            Fns.push_back( r_create_fn_from_support( P1[k], {}, 0 ) );
            mi_Fns.push_back( MI( {Fns[k]}, {o_idx} ) );
          }
          quicksort_by_attribute( Fns, mi_Fns, 0, (Fns.size()-1) );
          if( Fns.size() == 1 )
          {
            return Fns[0];
          }
          else
          {
            auto Fold = Fns[0];
            for( uint32_t j{1u}; j < Fns.size(); ++j )
            {
              auto tt_new = create_fn( {Fold, Fns[j]} );
              create_klut_node( {Fold, Fns[j]}, tt_new );
              Fold = ( _num_nodes - 1 );
            }
            return Fold;
          }
        }

      }
      else
      {
        if ( p.size()+given_klg.size() <= _max_sup )
        {
          for ( uint32_t k{0u}; k<given_klg.size(); ++k )
            p.push_back( given_klg[k] );

          return r_create_fn_from_support( p, {}, 0 );
        }
        else // |p|+|given| > kmax
        {
          auto y = p[0];
          p.erase(p.begin());
          auto f0 = r_create_fn_from_support( p, {y}, 0 );
          return r_create_fn_from_support( {f0}, given_klg, 0 );
        }
      }
    }

    void group_by_symmetry( std::vector<uint32_t>& support, uint32_t o_idx = 0 )
    {
      /* compute the MI of all the nodes */
      std::vector<double> mi_v;
      for( uint32_t k{0u}; k < support.size(); ++k )
        mi_v.emplace_back( MI( { support[k] }, { o_idx }) );

      for( uint32_t i{0u}; i<support.size(); ++i )
        std::cout << support[i] << " " << mi_v[i] << std::endl;

      quicksort_by_attribute( support, mi_v, 0, (support.size()-1) );

      for( uint32_t i{0u}; i<support.size(); ++i )
        std::cout << support[i] << " " << mi_v[i] << std::endl;

      auto Pi = group_by_mi( support, mi_v );
      for (uint32_t k{0u}; k<Pi.size(); ++k )
      {
        auto p = Pi[k];
        std::cout << "\n{ ";
        for( uint32_t j{0u}; j <p.size(); ++j )
        {
          std::cout << p[j] << " ";
        }
        std::cout << "}\n";
        if ( p.size() > 1 )
        {
          if( p.size() <= _max_sup )
          {
            auto tt_new = create_fn( p );
            create_klut_node( p, tt_new );
          }
          else
          {
            r_create_fn_from_support( p, {}, 0 );
          }
        }
      }
    }
    #pragma endregion

    #pragma region preprocess_muesli
    void preprocess_muesli( )
    {
      std::vector<uint32_t> support;
      for( auto k = 0u; k < _num_nodes; ++k )
        support.push_back(k);
      group_by_symmetry( support );


    }
    #pragma endregion

    /* return klut signal */
    #pragma region it_shannon_decomposition
    uint64_t it_shannon_decomposition_step( std::vector<uint32_t> support, dbs_storage nodes_remaining, dbs_storage outputs_remaining, uint32_t o_idx = 0 )
    {
      //std::cout << 1 << std::endl;
      uint32_t num_nodes = nodes_remaining[0].size() - 1;
      double mi_max = 0;
      double mi_new;
      uint32_t x_s;

      bool all_ones = true;
      bool all_zeros = true;
      //std::cout << 2 << std::endl;

      for( uint32_t k{0u}; k < outputs_remaining.size(); ++k )
      {
        if( outputs_remaining[k][o_idx] == 0 )
          all_ones = false;
        else if ( outputs_remaining[k][o_idx] == 1 )
          all_zeros = false;
        else
          std::cerr << "none valid " << std::endl;
      }
      //std::cout << 3 << std::endl;

      if( all_ones == true )
        return klut.get_constant( true );
      
      if( all_zeros == true )
        return klut.get_constant( false );

      //std::cout << 4 << std::endl;


      for( uint32_t k{0u}; k<nodes_remaining.size(); ++k )
        //std::cout << nodes_remaining[k] << ": " << outputs_remaining[k] << std::endl;

      if( support.size() < _max_sup )
      {
        auto tt_tmp = create_fn_gd( support, nodes_remaining, outputs_remaining );
        //std::cout << tt_tmp << std::endl;
        create_klut_node( support, tt_tmp );
        return _itos.storage[_num_nodes];
      }

      //std::cout << 5 << std::endl;

      for( uint32_t k{0u}; k < support.size(); ++k )
      {
        /* print */
        //std::cout << "supp el: " << support[k] << " " << k << "/" << support.size() << std::endl;
        //std::cout << nodes_remaining.size() << " " << outputs_remaining.size() << std::endl;
        //for( uint32_t jj {0u}; jj<nodes_remaining.size(); ++jj )
          //std::cout << nodes_remaining[jj] << ":" << outputs_remaining[jj]<< std::endl;
        /* print */
        mi_new = MI_gd( { support[k] },{o_idx}, nodes_remaining, outputs_remaining, support.size() );
        //std::cout << 5.1 << std::endl;

        if( mi_new >= mi_max )
        {
          mi_max = mi_new;
          x_s = support[k];
        }
      }
      dbs_storage nodes0, nodes1, outputs0, outputs1;   /* storage element: value of the output at each example */

      //std::cout << 6 << std::endl;

      std::vector<uint32_t> new_support;
      dyn_bitset mask ( num_nodes + 1, 1u );
      mask = mask << x_s;
      for (uint32_t k {0u}; k < nodes_remaining.size(); ++k )
      {
        //std::cout << 7 << std::endl;
        //std::cout << mask << std::endl;
        //std::cout << nodes_remaining[k] << std::endl;


        if ( ( mask & nodes_remaining[k] ) == mask ) /* f1 */
        {
          //std::cout << 8 << std::endl;
          boost::dynamic_bitset<> new_bs;
          for ( uint32_t j{0u}; j < support.size(); ++j )
          {
            //std::cout << 9 << std::endl;
            if( support[j] != support[x_s] )
            {
              new_bs.push_back( nodes_remaining[k][j] );
            } 
          }
          new_bs.push_back( 0 );
          nodes1.push_back( new_bs );
          outputs1.push_back( outputs_remaining[k] );
        }
        else /* f0 */
        {
          //std::cout << 10 << std::endl;
          boost::dynamic_bitset<> new_bs;
          for ( uint32_t j{0u}; j < support.size(); ++j )
          {
            //std::cout << 11 << std::endl;
            if( support[j] != support[x_s] )
            {
              new_bs.push_back( nodes_remaining[k][j] );
            } 
          }
          new_bs.push_back( 0 );

          nodes0.push_back( new_bs );
          outputs0.push_back( outputs_remaining[k] );
        }
      }
      for ( uint32_t j{0u}; j < support.size(); ++j )
      {
        //std::cout << 12 << std::endl;
        if( support[j] != support[x_s] )
          new_support.push_back( support[j] );
      }
      //std::cout << 13 << std::endl;
      auto f1 = klut.create_and( _itos.storage[support[x_s]], it_shannon_decomposition_step( new_support, nodes1, outputs1, 0) );
      auto f0 = klut.create_and( klut.create_not(_itos.storage[support[x_s]]), it_shannon_decomposition_step( new_support, nodes0, outputs0, 0) );

      return klut.create_or( f1, f0 );
      /* construct the substorage blocks */
    }

    void it_shannon_decomposition( uint32_t o_idx = 0 )
    {
      std::vector<uint32_t> initial_support;
      for( uint32_t k{0u}; k < _num_nodes; ++k )
        initial_support.push_back( k );
      
      auto f0 = it_shannon_decomposition_step( initial_support, _nodes, _outputs, 0 );
      klut.create_po( f0 );
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
      uint32_t _max_sup;
      uint32_t _init_sup;
      double _eps_th;
      double _eps_best;
      uint32_t _idx_fn;



  };

} // namespace mockturtle