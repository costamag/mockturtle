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
#include <algorithm>
#include <string>
#include <vector>
#include <bitset>
#include <cmath>
#include <random>
#include <boost/dynamic_bitset.hpp>
#include <unordered_map>

#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>


#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>

#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/algorithms/aig_algebraic_rewriting.hpp>


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


  class plaT0_network
  {
    #pragma region Types and constructors
    using dyn_bitset = boost::dynamic_bitset<>;
    using dbs_storage = std::vector<dyn_bitset>;


    public:
      plaT0_network( dbs_storage input_nodes, dbs_storage output_nodes, uint64_t max_act, uint64_t max_sup = 2, uint64_t init_sup = 2 )
      : _input_nodes( input_nodes ),
        _nodes( input_nodes ),
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
        std::vector<uint64_t> Xindeces;
        for ( uint64_t i {0u}; i < _num_nodes ; ++i )
        {
          auto pi = klut.create_pi();
          _itos.insert( i, pi );

          MI({i},{0});
        }
        _act = 0;

      }
    #pragma endregion
    
    #pragma region visual
    public:
      void print_pla()
      {
        for ( uint64_t i {0u}; i < _num_data; ++i )
            std::cout << _outputs.at(i) << ":" << _nodes.at(i)  << std::endl ;
      }
      void print_pla_gd( dbs_storage const& nodes, dbs_storage const& outputs )
      {
        for ( uint64_t i {0u}; i < nodes.size(); ++i )
            std::cout << outputs.at(i) << ":" << nodes.at(i)  << std::endl ;
      }

      void print_probabilities( std::vector<double> probabilities )
      {
        uint64_t num_vars = probabilities.size();
        std::cout << std::endl;
        for ( uint64_t mask {0u}; mask < num_vars; ++mask )
        {
          dyn_bitset BS ( (uint64_t)log2(num_vars), mask );
          std::cout << "|P(" << BS << ") = " << probabilities.at(mask) << std::endl ;
        }
        std::cout << std::endl;

      }

      void print_active_list()
      {
        std::cout << "\nactive list:";
        for(uint64_t k{0u}; k<_active_list.size(); ++k)
        {
          std::cout << "k:";
          std::cout << _active_list[k] << "] ";
        }
        std::cout << "\n";
      }
    #pragma endregion
      
    #pragma region basic_function
    std::vector<double> Pr( std::vector<uint64_t> indeces_nodes, std::vector<uint64_t> indeces_outputs )
    {
      uint64_t size_P_space = std::pow( 2, indeces_nodes.size() + indeces_outputs.size() );
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
        uint64_t jeff;
        for ( uint64_t j {0u}; j < indeces_nodes.size(); ++j )
        {
          jeff = indeces_outputs.size() + j;
          mask_nodes |= ( b1_nodes << indeces_nodes.at(j) );
          X_nodes |= ( ( ( ( b1_nodes << jeff ) & xin ) >> jeff ) << indeces_nodes.at(j) );
        }

        dyn_bitset mask_outputs ( _num_outputs, 0u );
        uint64_t u64_t_X_outputs {0u};

        for ( uint64_t j {0u}; j < indeces_outputs.size(); ++j )
        {
          mask_outputs |= ( b1_outputs << indeces_outputs.at(j) );
          u64_t_X_outputs |= ( ( ( ( 1u << j ) & x_u64_t ) >> j ) << indeces_outputs.at(j) );
        }
        dyn_bitset X_outputs( _num_outputs, u64_t_X_outputs );

        proba = 0;

        for ( uint64_t i {0u}; i < _num_data; ++i )
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
      uint64_t size_P_space = proba.size();
      double entropy { 0 };
      double deltaH { 0 };

      for ( uint64_t xin {0u}; xin < size_P_space; ++xin )
      {
        deltaH = ( proba[xin] == 0 ) ? 0 : -1*proba[xin]*log2( proba[xin] );
        entropy += deltaH;
      } 
      return entropy;
    }

    double MI ( std::vector<uint64_t> Xindeces, std::vector<uint64_t> Yindeces, bool overwrite = false )
    {
      std::stringstream ss;
      std::copy( Xindeces.begin(), Xindeces.end(), std::ostream_iterator<int>(ss, " "));
      std::string s = ss.str();
      s = s.substr(0, s.length()-1);
      bool not_present = (_mi_storage.find(s) == _mi_storage.end());
      if ( not_present || overwrite )
      {
        auto Hx = H( Xindeces, {} );
        auto Hy = H( {}, Yindeces );
        auto Hxy = H( Xindeces, Yindeces ); 

        if( not_present )
          _mi_storage.insert(std::make_pair(s,(Hx + Hy - Hxy)));
        else
          _mi_storage.at(s)=Hx + Hy - Hxy;
          
        //std::cout << "I(" << s << ";f)=" << _mi_storage.at(s) << std::endl;
        return _mi_storage.at(s);
      }
      else
      {
        //std::cout << "I(" << s << ";f)=" << _mi_storage.at(s) << std::endl;
        return _mi_storage.at(s);
      }
    }

    #pragma endregion

    #pragma region basic_function_given_data
    std::vector<double> Pr_gd( std::vector<uint64_t> const& indeces_nodes, std::vector<uint64_t> const& indeces_outputs, 
                            dbs_storage const& nodes, dbs_storage const& outputs, 
                            uint64_t const& num_nodes )
    {
      uint64_t size_P_space = std::pow( 2, indeces_nodes.size() + indeces_outputs.size() );
      std::vector<double> probabilities;
      double_t proba;
      double eq_flag_nodes, eq_flag_outputs;
      dyn_bitset b1_nodes ( num_nodes, 1u );
      dyn_bitset b1_outputs ( _num_outputs, 1u );
      auto num_data = nodes.size();


      for ( uint64_t x_u64_t {0u}; x_u64_t < size_P_space; ++x_u64_t )
      {
        dyn_bitset xin ( num_nodes, x_u64_t );
        dyn_bitset mask_nodes ( num_nodes, 0u );
        dyn_bitset X_nodes ( num_nodes, 0u );
        uint64_t jeff;
        for ( uint64_t j {0u}; j < indeces_nodes.size(); ++j )
        {
          jeff = indeces_outputs.size() + j;
          mask_nodes |= ( b1_nodes << indeces_nodes.at(j) );
          X_nodes |= ( ( ( ( b1_nodes << jeff ) & xin ) >> jeff ) << indeces_nodes.at(j) );
        }

        dyn_bitset mask_outputs ( _num_outputs, 0u );
        uint64_t u64_t_X_outputs {0u};

        for ( uint64_t j {0u}; j < indeces_outputs.size(); ++j )
        {
          mask_outputs |= ( b1_outputs << indeces_outputs.at(j) );
          u64_t_X_outputs |= ( ( ( ( 1u << j ) & x_u64_t ) >> j ) << indeces_outputs.at(j) );
        }
        dyn_bitset X_outputs( _num_outputs, u64_t_X_outputs );

        proba = 0;
        for ( uint64_t i {0u}; i < num_data; ++i )
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
        probabilities.push_back( proba );
      }

      return probabilities;

    }

    double H_gd( std::vector<uint64_t>  const& indeces_nodes, std::vector<uint64_t> const& indeces_outputs, 
                            dbs_storage const& nodes, dbs_storage const& outputs, 
                            uint64_t const& num_nodes )
    {

      auto proba = Pr_gd( indeces_nodes, indeces_outputs, nodes, outputs, num_nodes );
      uint64_t size_P_space = proba.size();
      double entropy { 0 };
      double deltaH { 0 };

      for ( uint64_t xin {0u}; xin < size_P_space; ++xin )
      {
        deltaH = ( proba[xin] == 0 ) ? 0 : -1*proba[xin]*log2( proba[xin] );
        entropy += deltaH;
      } 
      return entropy;
    }

    double MI_gd ( std::vector<uint64_t> const& Xindeces, std::vector<uint64_t> const& Yindeces,
                   dbs_storage const& nodes, dbs_storage const& outputs, uint64_t num_nodes )
    {
      auto Hx = H_gd( Xindeces, {}, nodes, outputs, num_nodes );
      auto Hy = H_gd( {}, Yindeces, nodes, outputs, num_nodes );
      auto Hxy = H_gd( Xindeces, Yindeces, nodes, outputs, num_nodes ); 

      return ( Hx + Hy - Hxy );
    }

    std::string create_fn_gd( std::vector<uint64_t> support, dbs_storage& nodes, dbs_storage outputs )
    {
      uint64_t num_nodes = nodes[0].size();
      uint64_t num_data = nodes.size();
      uint64_t nin_node = support.size();
      //std::cout << "supp size = " << nin_node << std::endl;
      uint64_t domain_size = pow( 2, nin_node );
      uint64_t Ci0, Ci1;
      for( uint64_t k{0u}; k < nodes.size(); ++k )
      {
        nodes[k].push_back(0);
      }

      dyn_bitset mask ( num_nodes+1, 0u );
      dyn_bitset X ( num_nodes+1, 0u );
      dyn_bitset Bit1 ( num_nodes+1, 1u );
      dyn_bitset Bit0 ( num_nodes+1, 0u );

      dyn_bitset Bit1_outputs ( _num_outputs, 1u );

      std::default_random_engine generator;
      std::bernoulli_distribution distribution(0.5);
      std::string tt_str;
      
      auto mask0 = ~( Bit1 << num_nodes );

      for ( uint64_t j {0u}; j < nodes.size(); ++j )
      {
        nodes.at(j) &= mask0; 
      }


      for ( uint64_t x_u64_t {0u}; x_u64_t < domain_size; ++x_u64_t )
      {
        dyn_bitset xin ( num_nodes+1, x_u64_t );
        Ci0 = 0;
        Ci1 = 0;
        mask = Bit0;
        X = Bit0;
        for ( uint64_t j {0u}; j < support.size(); ++j )
        {
          mask |= ( Bit1 << support.at(j) );
          X |= ( ( ( ( Bit1 << j ) & xin ) >> j ) << support.at(j) );
        }

        for ( uint64_t j {0u}; j < nodes.size(); ++j )
        {
          if ( X == ( mask & nodes.at(j) ) )
            ( ( outputs.at(j) & Bit1_outputs ) == Bit1_outputs ) ? Ci1++ : Ci0++;
        }

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

        for ( uint64_t j {0u}; j < num_data; ++j )
        {
          if ( X == ( mask & nodes.at(j) ) )
            nodes[j] |= new_val ;

        }
      }
      //std::cout << "str: " << tt_str << std::endl;
      
      return tt_str;

    }

    #pragma endregion
//########################################################################################################
    #pragma region new_node_given_data
    std::vector<uint64_t> active_list_gd( dbs_storage const& nodes_remaining, dbs_storage const& outputs_remaining )
    {
      double mi_loc;
      double mi_max = 0;
      uint64_t idx;
      std::vector<uint64_t> active_list;
      /* first active variable */
      for( uint64_t i {0u}; i < nodes_remaining[0].size(); ++i )
      {
        mi_loc = MI_gd( {i}, {0}, nodes_remaining, outputs_remaining, nodes_remaining[0].size() );
        uint64_t idx; 

        if ( mi_loc >= mi_max )
        {
          mi_max = mi_loc;
          idx = i;
          active_list = {idx};
        }
      }

      std::vector<uint64_t> inv_indeces;

      for ( uint64_t i {1u}; i < _max_act; ++i )
      {
        mi_max = 0;
        inv_indeces = active_list;
        inv_indeces.emplace_back(0);
        for ( uint64_t j {0u}; j < nodes_remaining[0].size(); ++j )
        { 
          if ( std::find( active_list.begin(), active_list.end(), j ) != active_list.end() )
            continue;
          
          inv_indeces.at(i) = j;

          mi_loc = MI_gd( inv_indeces, {0}, nodes_remaining, outputs_remaining, nodes_remaining[0].size() );
          if ( mi_loc >= mi_max )
          {
            mi_max = mi_loc;
            idx = j;
          }
        }
      active_list.push_back(idx);
      }
      return active_list;
    }
    #pragma endregion
//########################################################################################################


    #pragma region new_node
    void fill_active_list( )
    {
      double mi_loc;
      double mi_max = 0;
      uint64_t idx;
      /* first active variable */
      for( uint64_t i {0u}; i < _num_nodes; ++i )
      {
        mi_loc = MI( {i}, {0} );
        uint64_t idx; 

        if ( mi_loc >= mi_max )
        {
          mi_max = mi_loc;
          idx = i;
          _active_list = {idx};
        }
      }

      std::vector<uint64_t> inv_indeces;

      for ( uint64_t i {1u}; i < _max_act; ++i )
      {
        mi_max = 0;
        inv_indeces = _active_list;
        inv_indeces.emplace_back(0);
        for ( uint64_t j {0u}; j < _num_nodes; ++j )
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

    std::string create_fn( std::vector<uint64_t> support )
    {

      uint64_t nin_node = support.size();
      uint64_t domain_size = pow( 2, nin_node );
      uint64_t Ci0, Ci1;
      dyn_bitset mask ( _num_nodes + 1, 0u );
      dyn_bitset X ( _num_nodes + 1, 0u );
      dyn_bitset Bit1 ( _num_nodes + 1, 1u );
      dyn_bitset Bit0 ( _num_nodes + 1, 0u );

      dyn_bitset Bit1_outputs ( _num_outputs, 1u );

      std::default_random_engine generator;
      std::bernoulli_distribution distribution(0.5);
      std::string tt_str;
      
      auto mask0 = ~( Bit1 << _num_nodes );

      for ( uint64_t j {0u}; j < _num_data; ++j )
      {
        _nodes.at(j) &= mask0; 
      }

      for ( uint64_t x_u64_t {0u}; x_u64_t < domain_size; ++x_u64_t )
      {
        dyn_bitset xin ( _num_nodes+1, x_u64_t );
        Ci0 = 0;
        Ci1 = 0;
        mask = Bit0;
        X = Bit0;

        for ( uint64_t j {0u}; j < support.size(); ++j )
        {
          mask |= ( Bit1 << support.at(j) );
          X |= ( ( ( ( Bit1 << j ) & xin ) >> j ) << support.at(j) );
        }

        for ( uint64_t j {0u}; j < _num_data; ++j )
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

        for ( uint64_t j {0u}; j < _num_data; ++j )
        {
          if ( X == ( mask & _nodes.at(j) ) )
            _nodes.at(j) |= new_val;
        }
      }
      
      return tt_str;

    }

    void create_klut_node( std::vector<uint64_t> support, std::string tt_str )
    {
      auto supp_size = support.size();
      kitty::dynamic_truth_table tt( supp_size );
      create_from_binary_string( tt, tt_str );
      //std::cout << "TT: " ;
      //kitty::print_binary(tt);
      std::vector<uint64_t> klut_signals;
      for ( uint64_t i {0u}; i < supp_size; ++i )
      {
        klut_signals.push_back( _itos.storage[support[i]] );
      }
      auto f0 = klut.create_node( klut_signals, tt );
      //std::cout <<"new signal:" << f0 << std::endl;
      _itos.insert( _num_nodes ,f0 );
      _num_nodes++;

      dbs_storage dbs_nodes;
      dbs_storage dbs_outputs;
      
      //print_pla();

      //for ( uint64_t k {0u}; k < _num_data; ++k )
      //{
        //_nodes.at(k).push_back(0);
      //}
    }


    bool improve_fn( )
    {
      std::vector<uint64_t> support;
      std::string tt_str;
      

      fill_active_list( );
      print_active_list( );
      support = {};
      if( _act+_sup > _active_list.size() )
      {
        return false;
      }
      //std::cout << "support: ";
      for ( uint64_t k{0u}; k < _sup; ++k )
      {
        //std::cout << _active_list.at( _act + k ) << " ";
        support.push_back( _active_list.at( _act + k ));
      }
      //std::cout << std::endl;
      
      std::vector<uint64_t> first_act;
      for (uint64_t k {0u}; k <= _act; ++k )
      {
        first_act.push_back(_active_list.at(k));
      }
      auto mi_old = MI( first_act, {0} ); //#########################################################################
      
      tt_str = create_fn( support );
      //std::cout << "truth table: " <<tt_str << std::endl;
      first_act.at(_act) = _num_nodes;
      //print_pla();
      auto mi_new = MI( first_act , {0}, true );
      
      //std::cout << "mi_new " << mi_new << std::endl;
      //std::cout << "mi_old " << mi_old << std::endl;

      if ( mi_new > mi_old )
      {
        //std::cout << "new node created. Stored at: " << _num_nodes << std::endl;
        create_klut_node( support, tt_str );
        return true;
      }
      else
      {
      std::stringstream ss;
      std::copy( first_act.begin(), first_act.end(), std::ostream_iterator<int>(ss, " "));
      std::string s = ss.str();
      s = s.substr(0, s.length()-1);
        _mi_storage.at(s) = mi_old;
      }
      //std::cout << "Node is not kept. Remove new function " << std::endl;
      return false;
    }

    bool not_done( uint64_t best_idx )
    {
      double eps_I_H = MI( {best_idx}, {0} )/H({},{0});
      //std::cout << "I(n*;f)/H(f)= " << eps_I_H << std::endl;
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
      std::vector<uint64_t> support; /* contains the support of the new function */
      uint64_t best_idx; /* used to keep track of the new best approximation */

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
          {
            best_idx -= 1;
          }
               
          //std::cout << "best idx = " << best_idx << std::endl;
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
            //std::cout << "best idx = " << best_idx << std::endl;
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
      std::cout << "node with maximum mutual information is n*=" << _active_list[0] << std::endl;
      std::cout << "maximum mutual information is I(n*;f)=" << MI({best_idx},{0}) << std::endl;
      auto f0 = klut.create_po(_itos.storage[_active_list[0]]);
      _training_accuracy = compute_accuracy( _input_nodes, _outputs );
      std::cout << "training accuracy: " << _training_accuracy << "%" << std::endl;
    }

    void muesli_modified( double eps_th = 0.99 )
    {

      bool success; /* true if found a function improving the mi */
      std::string tt_str; /* contains the tt of the new node */
      std::vector<uint64_t> support; /* contains the support of the new function */
      uint64_t best_idx; /* used to keep track of the new best approximation */

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

    int partition ( std::vector<uint64_t>& support, std::vector<double>& attribute , uint64_t low, uint64_t high )
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

    void quicksort_by_attribute( std::vector<uint64_t>& support, std::vector<double>& attribute,  int low, int high )
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

    std::vector<std::vector<uint64_t>> group_by_mi( std::vector<uint64_t> const& support, std::vector<double> const& mi_v )
    {
      std::vector<std::vector<uint64_t>> Pi;
      std::vector<double> miP;
      uint64_t idxPi = 0;
      Pi.push_back({support[0]});
      miP.push_back(mi_v[0]);
      for( uint64_t k{1u}; k<support.size(); ++k )
      {
        if( mi_v[k] >= ( miP[idxPi] - miP[idxPi]*_dI) )
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

    uint64_t r_create_fn_from_support( std::vector<uint64_t> p, std::vector<uint64_t> given_klg = {}, uint64_t o_idx = 0 )
    {
      std::cout << "\n{ ";
        for( uint64_t j{0u}; j <p.size(); ++j )
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
          std::vector<uint64_t> p1;

          for( uint64_t k{1u}; k < p.size(); ++k )
          {
            mi_v.emplace_back( MI( { p[k], x }, { o_idx }) );
            p1.emplace_back( p[k] );
          }
          quicksort_by_attribute( p1, mi_v, 0, (p1.size()-1) );
          auto P1 = group_by_mi( p1, mi_v );
          std::vector<uint64_t> Fns;
          std::vector<double> mi_Fns;

          Fns.push_back( r_create_fn_from_support( P1[0], {x}, 0 ) );
          mi_Fns.push_back( MI( {Fns[0]}, {o_idx} ) );
          for ( uint64_t k {1u}; k<P1.size(); ++k )
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
            for( uint64_t j{1u}; j < Fns.size(); ++j )
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
          for ( uint64_t k{0u}; k<given_klg.size(); ++k )
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

    void group_by_symmetry( std::vector<uint64_t>& support, uint64_t o_idx = 0 )
    {
      /* compute the MI of all the nodes */
      std::vector<double> mi_v;
      for( uint64_t k{0u}; k < support.size(); ++k )
        mi_v.emplace_back( MI( { support[k] }, { o_idx }) );

      for( uint64_t i{0u}; i<support.size(); ++i )
        std::cout << support[i] << " " << mi_v[i] << std::endl;

      quicksort_by_attribute( support, mi_v, 0, (support.size()-1) );

      for( uint64_t i{0u}; i<support.size(); ++i )
        std::cout << support[i] << " " << mi_v[i] << std::endl;

      auto Pi = group_by_mi( support, mi_v );
      for (uint64_t k{0u}; k<Pi.size(); ++k )
      {
        auto p = Pi[k];
        std::cout << "\n{ ";
        for( uint64_t j{0u}; j <p.size(); ++j )
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
    void preprocess_muesli( double dI = 0 )
    {
      _dI = dI;
      std::vector<uint64_t> support;
      for( auto k = 0u; k < _num_nodes; ++k )
        support.push_back(k);
      group_by_symmetry( support );


    }
    #pragma endregion

    /* return klut signal */
    #pragma region it_shannon_decomposition
    uint64_t it_shannon_decomposition_step( std::vector<uint64_t> support, dbs_storage nodes_remaining, dbs_storage outputs_remaining, bool is_dec_naive = false, uint64_t o_idx = 0 )
    {

      if( nodes_remaining.size() == 0 )
        return klut.get_constant( false );
      
      if( nodes_remaining[0].size() == 0 )
        return klut.get_constant( false );

      if( nodes_remaining[0].size() != support.size() )
        std::cerr << "not same size" << std::endl;

      uint64_t num_nodes = nodes_remaining[0].size();
      double mi_max = 0;
      double mi_new;
      uint64_t x_s;

      bool all_ones = true;
      bool all_zeros = true;

      for( uint64_t k{0u}; k < outputs_remaining.size(); ++k )
      {
        if( outputs_remaining[k][o_idx] == 0 )
          all_ones = false;
        else if ( outputs_remaining[k][o_idx] == 1 )
          all_zeros = false;
        else
          std::cerr << "none valid " << std::endl;
      }

      if( all_ones == true )
        return klut.get_constant( true );
      
      if( all_zeros == true )
        return klut.get_constant( false );

      if( support.size() <= _max_sup )
      {
        auto tt_tmp = create_fn_gd( support, nodes_remaining, outputs_remaining );
        //std::cout << tt_tmp << std::endl;
        create_klut_node( support, tt_tmp );
        return _itos.storage[_num_nodes-1];
      }

      if ( is_dec_naive )
      {
        x_s = 0;//support[0];
      }
      else
      {
        for( uint64_t k{0u}; k < support.size(); ++k )
        {
          mi_new = MI_gd( { k },{o_idx}, nodes_remaining, outputs_remaining, support.size() ); // chanched support[k] -> k
          if( mi_new >= mi_max )
          {
            mi_max = mi_new;
            x_s = k;//support[k];
          }
        }
      }

      dbs_storage nodes0, nodes1, outputs0, outputs1;   /* storage element: value of the output at each example */

      std::vector<uint64_t> new_support;
      dyn_bitset mask ( num_nodes, 1u );
      mask = mask << x_s;
      for (uint64_t k {0u}; k < nodes_remaining.size(); ++k )
      {

        if ( ( mask & nodes_remaining[k] ) == mask ) /* f1 */
        {
          boost::dynamic_bitset<> new_bs;
          for ( uint64_t j{0u}; j < support.size(); ++j )
          {
            if( support[j] != support[x_s] )
            {
              new_bs.push_back( nodes_remaining[k][j] );
            } 
          }
          nodes1.push_back( new_bs );
          outputs1.push_back( outputs_remaining[k] );
        }
        else /* f0 */
        {
          boost::dynamic_bitset<> new_bs;
          for ( uint64_t j{0u}; j < support.size(); ++j )
          {
            if( support[j] != support[x_s] )
            {
              new_bs.push_back( nodes_remaining[k][j] );
            } 
          }

          nodes0.push_back( new_bs );
          outputs0.push_back( outputs_remaining[k] );
        }
      }
      for ( uint64_t j{0u}; j < support.size(); ++j )
      {
        if( support[j] != support[x_s] )
          new_support.push_back( support[j] );
      }
      auto fa1 = it_shannon_decomposition_step( new_support, nodes1, outputs1,is_dec_naive, 0);
      auto f1 = klut.create_and( _itos.storage[support[x_s]], fa1 );
      auto fa0 = it_shannon_decomposition_step( new_support, nodes0, outputs0, is_dec_naive, 0);
      auto f0 = klut.create_and( klut.create_not(_itos.storage[support[x_s]]), fa0 );

      auto fn = klut.create_or( f1, f0 );
      
      return fn;
      /* construct the substorage blocks */
    }

    void it_shannon_decomposition( bool is_dec_naive = false, uint64_t o_idx = 0 )
    {

      std::vector<uint64_t> initial_support;
      dbs_storage nodes;
      for( uint64_t k{0u}; k < _num_nodes; ++k )
        initial_support.push_back( k );
      
      for( uint64_t d{0u}; d < _nodes.size(); ++d )
      {
        boost::dynamic_bitset<> dtmp;
        for( uint64_t k{0u}; k < _num_nodes; ++k )
        {
          dtmp.push_back(_nodes[d][k]);
        } 
        nodes.push_back(dtmp); 
      }

      auto f0 = it_shannon_decomposition_step( initial_support, nodes, _outputs, is_dec_naive, 0 );
      klut.create_po( f0 );
      _training_accuracy = compute_accuracy( _input_nodes, _outputs );
      std::cout << "training accuracy: " << _training_accuracy << "%" << std::endl;
      std::cout << "number of gates: " << aig.num_gates()<< std::endl;
      std::cout << "size: " << aig.size()<< std::endl;
      depth_view_params ps;
      ps.count_complements = true;
      depth_view depth_aig{aig, {}, ps};
      std::cout << "num levels: " << depth_aig.depth()<< std::endl;
    }
    #pragma endregion

    // ##############################################################################################

#pragma region dsd_shannon


struct res_BD_type
{
  bool is_created;
  uint64_t signal;
  std::vector<uint64_t> Supp;
  std::vector<uint64_t> A;
  uint64_t idx_node;
  uint64_t idx_newFn;
  double mi;
  bool r_del;
  bool c_del;
  bool rc_del;
  std::string tt;

};


res_BD_type try_bottom_decomposition_EXP( std::vector<uint64_t>& support, dbs_storage& nodes_remaining, dbs_storage& outputs_remaining, double& MImax )
{
  // is_created = true if a new node is worth being added due to bottom decomposition
  // signal = signal of the klut node created 
  //std::cout << "in try bottom Imax =" << MImax << std::endl;
  res_BD_type res_BD;
  res_BD.is_created = false;

  dbs_storage new_nodes;
  dbs_storage nodes_tmp;
  std::vector<uint64_t> Apart;
  std::vector<uint64_t> Spart;  
  
// START ################################
  // consider all pairs and store if mi_supp == mi_Fnew AND mi_supp > mi_max
  nodes_tmp = nodes_remaining;
  for( uint64_t r = 0; r < nodes_remaining[0].size() ; ++r )
  {
    for( uint64_t c = r+1; c < nodes_remaining[0].size() ; ++c )
    {
      Apart = {r,c};
      Spart = {support[r], support[c]};
      //std::cout << "r:" << r << " c:" << c << std::endl; //X1
      //print_pla_gd( nodes_tmp, outputs_remaining );
      auto tt_tmp = create_fn_gd( Apart, nodes_tmp, outputs_remaining );

      //std::cout << tt_tmp << std::endl;
      //print_pla_gd( nodes_tmp, outputs_remaining );

      auto mi_supp = MI_gd( Apart, { 0 }, nodes_tmp, outputs_remaining, nodes_tmp[0].size() ); // support[k] -> k
      auto mi_Fnew = MI_gd( { nodes_tmp[0].size() - 1 }, { 0 }, nodes_tmp, outputs_remaining, nodes_tmp[0].size() ); // support[k] -> k
      auto mi_Fr = MI_gd( { (nodes_tmp[0].size() - 1), r }, { 0 }, nodes_tmp, outputs_remaining, nodes_tmp[0].size() );
      auto mi_Fc = MI_gd( { (nodes_tmp[0].size() - 1), c }, { 0 }, nodes_tmp, outputs_remaining, nodes_tmp[0].size() );
      auto mi_Frc = MI_gd( { (nodes_tmp[0].size() - 1), r, c }, { 0 }, nodes_tmp, outputs_remaining, nodes_tmp[0].size() );

      if ( mi_Fnew > MImax )//( ( mi_supp == mi_Fnew ) && ( mi_supp > MImax ) )
      {
        MImax = mi_Fnew;
        res_BD.is_created = true;
        res_BD.tt = tt_tmp;
        res_BD.Supp = Spart;
        res_BD.A = Apart;
        res_BD.mi = mi_Fnew;
        new_nodes = nodes_tmp;

        if( mi_Frc == mi_Fnew )
        {
          res_BD.rc_del = true; 
          res_BD.r_del = false; 
          res_BD.c_del = false; 
        }
        else if( mi_Fr == mi_Fnew )
        {
          res_BD.rc_del = false; 
          res_BD.r_del = true; 
          res_BD.c_del = false; 
        }
        else if( mi_Fc == mi_Fnew )
        {
          res_BD.rc_del = false; 
          res_BD.r_del = false; 
          res_BD.c_del = true; 
        }
        else
        {
          res_BD.rc_del = false; 
          res_BD.r_del = false; 
          res_BD.c_del = false; 
        }
      }
      std::vector<uint64_t> fls_support = support;
      remove_column( fls_support, nodes_tmp, nodes_remaining[0].size() );
    }
  }
  // modify
  if( res_BD.is_created )
  {
    std::cout << "created f(A[" << res_BD.A[0] << "],A[" << res_BD.A[1] << "])=f(" << 
                  res_BD.Supp[0] << "," << res_BD.Supp[1] << ")=" << res_BD.tt;

    //std::cout << "At the moment the number of nodes is "<< _num_nodes << " will bw index of new fn" << std::endl;
    //std::cout << "support has size " << support.size() << std::endl;
    res_BD.idx_node = _num_nodes;
    support.push_back(_num_nodes);
    //std::cout << "after inserting" << _num_nodes << " in the back, support has size " << support.size() << std::endl;
    create_klut_node( res_BD.Supp, res_BD.tt );
    //std::cout << "now. created node. _num_nodes is now " << _num_nodes << std::endl;
    //std::cout << "Indeed, signal stored at _num_nodes-1 in _itos is " << _itos.storage[_num_nodes-1] << std::endl;
    res_BD.signal = _itos.storage[_num_nodes-1];
    nodes_remaining = new_nodes;
    
    if ( res_BD.rc_del )
    {
      std::cout << " --> del {r,c}" << std::endl;
      remove_column( support, nodes_remaining, std::max(res_BD.A[0], res_BD.A[1]) );
      remove_column( support, nodes_remaining, std::min(res_BD.A[0], res_BD.A[1]) );
    }
    else if ( res_BD.r_del )
    {
      std::cout << " --> del {r}" << std::endl;
      remove_column( support, nodes_remaining, res_BD.A[0] );
    }
    else if ( res_BD.c_del )
    {
      std::cout << " --> del {c}" << std::endl;
      remove_column( support, nodes_remaining, res_BD.A[1] );
    }
    else
    {
      std::cout << " --> del {}" << std::endl;
    }
    res_BD.idx_newFn = nodes_remaining[0].size()-1;      

  }
// END ################################
  return res_BD;
}


void prepare_cofactor( dbs_storage const& nodes_remaining, dbs_storage const& outputs_remaining, 
                        uint64_t x_idx, bool Id, 
                        dbs_storage& nodesId, dbs_storage& outputsId )
{
  //std::cout << "Id=" << Id << std::endl;
  for( uint64_t dt {0u}; dt < nodes_remaining.size(); ++dt )
  {
    if( nodes_remaining[dt][x_idx] == Id )
    {
      outputsId.push_back(outputs_remaining[dt]);
      dyn_bitset dbset;
      for( uint64_t k{0u}; k < nodes_remaining[0].size(); ++k )
      {
        if( k != x_idx )
          dbset.push_back( nodes_remaining[dt][k] );
      }
      nodesId.push_back(dbset);
    }
  } // filled the cofactor

}

uint64_t HammingDistance( dyn_bitset& a, dyn_bitset& b  )
{
  auto a_xor_b = a^b;
  return a_xor_b.count();
}

bool is_f1_eqto_not_f0( dbs_storage const& nodes_remaining, dbs_storage const& outputs_remaining, uint64_t& count_max, uint64_t x_idx )
{
  uint64_t hd_max = 1;
  double Rt = 1.0;
  uint64_t count_x = 0;
  uint64_t count_neg = 0;

  dbs_storage nodes0, nodes1, outputs0, outputs1;
  prepare_cofactor( nodes_remaining, outputs_remaining, x_idx, 0, nodes0, outputs0 );
  prepare_cofactor( nodes_remaining, outputs_remaining, x_idx, 1, nodes1, outputs1 );
  
  for ( uint64_t n {0u}; n < nodes0.size(); ++n )
  {
    for ( uint64_t m {0u}; m < nodes1.size(); ++m )
    {
      
      if( HammingDistance( nodes0[n], nodes1[m] ) <= hd_max )
      {
        std::cout << nodes0[n] << std::endl;
        std::cout << nodes1[m] << std::endl;
        count_x++;
        count_neg += ( outputs0[n] == ~outputs1[m] );
      }
    } // explored N1
  }// explored N0

  if( ( count_neg >= count_max ) && ( count_neg >= Rt*count_x ) )
  {
    count_max = count_neg;
    return true;
  }
  return false;
}

bool is_f1_eqto_not_f0_hash_gd( dbs_storage const& nodes0, dbs_storage const& nodes1, 
                                dbs_storage const& outputs0, dbs_storage const& outputs1 )
{
  uint64_t count_neg = 0;
  std::unordered_map<std::string, double> str_nodes0;

  /* fill hash table */
  for( uint64_t k {0u}; k < nodes0.size(); ++k )
  {
    std::string s;
    to_string( nodes0[k], s );
    str_nodes0.insert(std::make_pair(s,outputs0[k][0]));
  }


  for ( uint64_t m {0u}; m < nodes1.size(); ++m )
  {
    std::string s ;
    to_string( nodes1[m], s );

    if( str_nodes0.find(s) != str_nodes0.end() )
    {
      //std::cout << s << ":" << str_nodes0.at(s) << std::endl;
      //std::cout << nodes1[m] << ":" << outputs1[m] << std::endl;
      if( str_nodes0.at(s) == outputs1[m][0] )
      {
        //std::cout << "F" << std::endl;
        return false;
      }
      else
      {
        //std::cout << "T" << std::endl;
        count_neg++;
      }
    }
  } 

  if( count_neg > 0 ) // CONSIDER CHANGING TO >= 0 if motivated
  {
    //std::cout << "x";
    return true;
  }
  return false;
}
          
void remove_column( std::vector<uint64_t>& support, dbs_storage& nodes_remaining, uint64_t x_s )
{
  dbs_storage new_nodes_remaining;
  std::vector<uint64_t> new_support;
  for( uint64_t dt {0u}; dt < nodes_remaining.size(); ++dt )
  {
    dyn_bitset dbset;
    for( uint64_t k {0u}; k < nodes_remaining[0].size(); ++k )
    {
      if( k != x_s )
        dbset.push_back( nodes_remaining[dt][k] );
    } 
    new_nodes_remaining.push_back( dbset );
  }

  for( uint64_t k {0u}; k < nodes_remaining[0].size(); ++k )
  {
    if( k != x_s )
      new_support.push_back( support[k] );
  }

  support = new_support;
  nodes_remaining = new_nodes_remaining;
}

void remove_column_and_invert( std::vector<uint64_t>& support, dbs_storage& nodes_remaining, dbs_storage& outputs_remaining, uint64_t x_s )
{
  dbs_storage new_nodes_remaining;
  std::vector<uint64_t> new_support;
  for( uint64_t dt {0u}; dt < nodes_remaining.size(); ++dt )
  {
    dyn_bitset dbset;
    for( uint64_t k {0u}; k < nodes_remaining[0].size(); ++k )
    {
      if( k != x_s )
        dbset.push_back( nodes_remaining[dt][k] );
    } 
    if ( nodes_remaining[dt][x_s] == 1 )
      outputs_remaining[dt][0] = ~outputs_remaining[dt][0];
    new_nodes_remaining.push_back( dbset );
  }

  for( uint64_t k {0u}; k < nodes_remaining[0].size(); ++k )
  {
    if( k != x_s )
      new_support.push_back( support[k] );
  }

  support = new_support;
  nodes_remaining = new_nodes_remaining;
}

bool check_if_all( dbs_storage const& outputs, bool const& val )
{

}
bool cec_all_val( dbs_storage const& outputs_remaining, bool val )
{
  bool ans = true;

  for( uint64_t k{0u}; k < outputs_remaining.size(); ++k )
  {
    if( outputs_remaining[k][0] != val )
      ans = false;
  }
  return ans;
}

uint64_t it_dsd_shannon_decomposition_step( std::vector<uint64_t> support, dbs_storage nodes_remaining, dbs_storage outputs_remaining, bool is_dec_naive = false, uint64_t o_idx = 0 )
    {
      std::cout << "\nN---------------------------------N" << std::endl;
      std::cout << "|s|[" << support.size() <<"]|s|.";
      std::cout << "|N|[" << nodes_remaining.size() <<"]|N|.";
      std::cout << "|No|[" << nodes_remaining[0].size() <<"]|oN|.";
      std::cout << "|O|[" << outputs_remaining.size() <<"]|O|."<< std::endl;

      assert( ( nodes_remaining.size() == outputs_remaining.size() ) ); // check nodes and outputs have the same length

      if( nodes_remaining.size() == 0 )
        return klut.get_constant( false );
      if( nodes_remaining[0].size() == 0 )
        return klut.get_constant( false );

      uint64_t num_nodes = nodes_remaining[0].size();
      assert( ( num_nodes == support.size() ) ); // check support length

      double mi_max = 0;
      double mi_new;
      uint64_t x_s=0;

      bool all_ones = cec_all_val(outputs_remaining, 1 );
      bool all_zeros = cec_all_val(outputs_remaining, 0 );

      if( all_ones == true )
      {
        //std::cout << "return true " << std::endl; //X 
        return klut.get_constant( true );
      }
      if( all_zeros == true )
      {
        //std::cout << "return false " << std::endl; //X 
        return klut.get_constant( false );
      }
  
      /* If support sufficiently small create function */
      if( support.size() <= _max_sup )
      {
        std::vector<uint64_t> supp_alt;
        for( uint64_t el = 0; el < support.size(); ++el )
          supp_alt.push_back(el);
        
        //print_pla_gd( nodes_remaining, outputs_remaining );
        auto tt_tmp = create_fn_gd( supp_alt, nodes_remaining, outputs_remaining );
        //std::cout << "TT=" << tt_tmp << std::endl;
        //std::cout << "_num_nodes=" << _num_nodes << std::endl;
        //std::cout << "storage size=" << _itos.storage.size() << std::endl;
        create_klut_node( support, tt_tmp );
        //std::cout << "_num_nodes=" << _num_nodes << std::endl;
        //std::cout << "storage size=" <<_itos.storage.size() << std::endl;
        return _itos.storage[_num_nodes-1];
      }

      for( uint64_t k{0u}; k < support.size(); ++k )
      {
        mi_new = MI_gd( { k }, { o_idx }, nodes_remaining, outputs_remaining, support.size() ); // support[k] -> k
        //std::cout << "new MI : " << mi_new << std::endl;
        //std::cout << "max MI : " << mi_max << std::endl;

        if( mi_new > mi_max )
        {
          //std::cout << "updated " << std::endl;
          //std::cout <<  std::endl;
          mi_max = mi_new;
          x_s = k;//support[k];
        }
      }
      Hnew = H_gd( { x_s }, { o_idx }, nodes_remaining, outputs_remaining, support.size() );
      Inew = MI_gd( { x_s }, { o_idx }, nodes_remaining, outputs_remaining, support.size() );
      std::cout << "eps[" << support(x_s) << "]=" << Inew/Hnew << std::endl; // support[k] -> k
      //std::cout << "xs= " << x_s << std::endl;
      //std::cout << "supp[xs]= " << support[x_s] << std::endl;

      dbs_storage nodes0, nodes1, outputs0, outputs1;   /* storage element: value of the output at each example */

      std::vector<uint64_t> new_support;
      std::vector<uint64_t> support_UP;

      /* fill cofactors */
      //print_pla_gd(nodes_remaining, outputs_remaining);

      prepare_cofactor( nodes_remaining, outputs_remaining, x_s, 0,  nodes0, outputs0 );
      prepare_cofactor( nodes_remaining, outputs_remaining, x_s, 1,  nodes1, outputs1 );

      //print_pla_gd(nodes0, outputs0);

      for ( uint64_t j{0u}; j < support.size(); ++j )
      {
        if( j != x_s )
          new_support.push_back( support[j] );
      }

      /* START checK
      std::cout << "supp size = " << support.size() << std::endl;
      std::cout << "nodes size = " << nodes_remaining.size() << std::endl;
      std::cout << "nodes[0] size = " << nodes_remaining[0].size() << std::endl;
      // END checK*/

      /* TOP DECOMPOSITION */
      
      bool is_F0_taut = cec_all_val(outputs0, 1 );
      bool is_F1_taut = cec_all_val(outputs1, 1 );
      bool is_F0_cont = cec_all_val(outputs0, 0 );
      bool is_F1_cont = cec_all_val(outputs1, 0 );

      if ( is_F1_taut )
      {
        std::cout << "F1=1 " << std::endl; //X
        auto F0 = it_dsd_shannon_decomposition_step( new_support, nodes0, outputs0, is_dec_naive, 0 );
        return klut.create_or( _itos.storage[support[x_s]], F0 );
      }
      else if ( is_F0_taut )
      {
        auto F1 = it_dsd_shannon_decomposition_step( new_support, nodes1, outputs1, is_dec_naive, 0 );
        std::cout << "F0=1 " << std::endl; //X
        return klut.create_le( _itos.storage[support[x_s]], F1 );
      }
      else if ( is_F1_cont )
      {
        auto F0 = it_dsd_shannon_decomposition_step( new_support, nodes0, outputs0, is_dec_naive, 0 );
        std::cout << "F1=0 " << std::endl; //X
        return klut.create_lt( _itos.storage[support[x_s]], F0 );
      }
      else if ( is_F0_cont )
      {
        auto F1 = it_dsd_shannon_decomposition_step( new_support, nodes1, outputs1, is_dec_naive, 0 );
        std::cout << "F0=0 " << std::endl; //X
        return klut.create_and( _itos.storage[support[x_s]], F1 );
      }
      else if( is_f1_eqto_not_f0_hash_gd( nodes0, nodes1, outputs0, outputs1 ) )
      {
          std::cout << "xXORf0' " << std::endl; //X
        auto pi_sig = _itos.storage[support[x_s]];
          //print_pla_gd(nodes_remaining, outputs_remaining); //X
          //std::cout << "RM+INV " << x_s << std::endl;//X
        remove_column_and_invert( support, nodes_remaining, outputs_remaining, x_s ); // checked correct
          //print_pla_gd(nodes_remaining, outputs_remaining); //X

        auto f0bar = it_dsd_shannon_decomposition_step( support, nodes_remaining, outputs_remaining, is_dec_naive, 0 );
          
        return klut.create_xor( pi_sig , f0bar );
      }
      else // xor
      {
        //print_pla_gd(nodes_remaining, outputs_remaining);
        res_BD_type res_BD = try_bottom_decomposition_EXP( support, nodes_remaining, outputs_remaining, mi_max ); 
        //print_pla_gd(nodes_remaining, outputs_remaining);


        if ( res_BD.is_created )
        {
          /* START checK
          std::cout << "BTM" << std::endl;
          std::cout << "supp size = " << support.size() << std::endl;
          std::cout << "nodes size = " << nodes_remaining.size() << std::endl;
          std::cout << "nodes[0] size = " << nodes_remaining[0].size() << std::endl;
          // END checK */
          //std::cout << "x xor f0'" << std::endl;
          std::cout << "BTM " << std::endl; //X
          return it_dsd_shannon_decomposition_step( support, nodes_remaining, outputs_remaining, is_dec_naive, 0 );
        }
        
      }

      // NO TOP DECOMPOSTION - BOTTOM DECOMPOSITION
      //for now only shannon
      std::cout << "SH1 " << std::endl; //X
      auto F1 = it_dsd_shannon_decomposition_step( new_support, nodes1, outputs1, is_dec_naive, 0 );
      std::cout << "SH0 " << std::endl; //X

      auto F0 = it_dsd_shannon_decomposition_step( new_support, nodes0, outputs0, is_dec_naive, 0 );
      auto f0 = klut.create_and( klut.create_not(_itos.storage[support[x_s]]), F0 );
      auto f1 = klut.create_and(_itos.storage[support[x_s]], F1 );

      return klut.create_or( f1, f0 );
      /* construct the substorage blocks */
    }

    void it_dsd_shannon_decomposition( bool is_dec_naive = false, uint64_t o_idx = 0 )
    {

      std::vector<uint64_t> initial_support;
      dbs_storage nodes;
      for( uint64_t k{0u}; k < _num_nodes; ++k )
        initial_support.push_back( k );
      
      for( uint64_t d{0u}; d < _nodes.size(); ++d )
      {
        boost::dynamic_bitset<> dtmp;
        for( uint64_t k{0u}; k < _num_nodes; ++k )
        {
          dtmp.push_back(_nodes[d][k]);
        } 
        nodes.push_back(dtmp); 
      }

      auto f0 = it_dsd_shannon_decomposition_step( initial_support, nodes, _outputs, is_dec_naive, 0 );
      std::cout << std::endl;
      klut.create_po( f0 );
      _training_accuracy = compute_accuracy( _input_nodes, _outputs );
      std::cout << "training accuracy: " << _training_accuracy << "%" << std::endl;
      std::cout << "number of gates: " << aig.num_gates()<< std::endl;
      std::cout << "size: " << aig.size()<< std::endl;
      depth_view_params ps;
      ps.count_complements = true;
      depth_view depth_aig{aig, {}, ps};
      std::cout << "num levels: " << depth_aig.depth()<< std::endl;
    }
    #pragma endregion // top dsd + shannon

    // ##############################################################################################

    #pragma region simulate
    bool simulate_input( dyn_bitset const& input_pattern, bool convertToAig = true )
    {
      if( convertToAig )
        aig = convert_klut_to_graph<aig_network>( klut );

      std::vector<bool> inpt_v;
      for( uint64_t k{0u}; k<input_pattern.size();++k )
      {
        inpt_v.push_back( ( ( input_pattern[k] == 1 ) ? true : false ) );
      }

      return simulate<bool>( aig, default_simulator<bool>( inpt_v ) )[0];
    }

    double compute_accuracy( dbs_storage const& nodes, dbs_storage const& outputs )
    {
      aig = convert_klut_to_graph<aig_network>( klut );
      //aig = cleanup_dangling( aig1 );
      //std::cout << "pre number of gates: " << aig.num_gates()<< std::endl;
      //std::cout << "pre size: " << aig.size()<< std::endl;
      //aig_algebraic_rewriting( aig );

      double acc = 0;
      double delta_acc;
      for( uint64_t k {0u}; k < nodes.size(); ++k )
      {
        dyn_bitset ipattern;
        for ( uint64_t j {0u}; j < (nodes[0].size()-1); ++j )
          ipattern.push_back( nodes[k][j] );
        
        delta_acc = ( ( simulate_input( ipattern, false ) == outputs[k][0] ) ? (double)1.0/nodes.size() : 0.0 );
        acc += delta_acc;
      }
      return 100*acc;
    }

    bool simulate_at_node( dyn_bitset const& input_pattern, klut_network const& klut_node )
    {
      aig_network aig_node = convert_klut_to_graph<aig_network>( klut_node );

      std::vector<bool> inpt_v;
      for( uint64_t k{0u}; k<input_pattern.size();++k )
      {
        inpt_v.push_back( ( ( input_pattern[k] == 1 ) ? true : false ) );
      }

      return simulate<bool>( aig_node, default_simulator<bool>( inpt_v ) )[0];
    }

    bool compare_nodes( dbs_storage const& nodes, klut_network const& klut_cec, bool is_same_sign )
    {
      bool equal = true;
      std::cout << "same? " << is_same_sign << std::endl;
      
      for( uint64_t k {0u}; k < nodes.size(); ++k )
      {
        dyn_bitset ipattern;
        for ( uint64_t j {0u}; j < nodes[0].size(); ++j )
          ipattern.push_back( nodes[k][j] );
        
        std::cout << simulate_at_node( ipattern, klut_cec ) << " ";

        if( is_same_sign )
        {
          if ( simulate_at_node( ipattern, klut_cec ) == 1 )
            return false;
        }
        else
        {
          
          if ( simulate_at_node( ipattern, klut_cec ) == 0 )
            return false;
        }
      }
      return true;
    }
    #pragma endregion

    public:
      dbs_storage _input_nodes;
      dbs_storage _nodes;   /* storage element: value of the output at each example */
      dbs_storage _outputs; /* storage element: value of the output at each example */
      uint64_t _num_data;   /* number of examples */
      uint64_t _num_nodes;
      uint64_t _num_outputs;
      klut_network klut;
      aig_network aig;
      std::vector<uint64_t> _active_list;
      index_to_signal _itos;
      uint64_t _act;
      uint64_t _sup = 0;
      uint64_t _max_act;
      uint64_t _max_sup;
      uint64_t _init_sup;
      double _eps_th;
      double _eps_best;
      uint64_t _idx_fn;
      double _training_accuracy;
      double _dI = 0;
      std::unordered_map<std::string, double> _mi_storage;

  };

} // namespace mockturtle