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
  \file simulation_view.hpp
  \brief View that attaches simulation patterns and allows for their simulation
  \author Andrea Costamagna
*/

#pragma once

#include "../../networks/klut.hpp"
#include "../../utils/node_map.hpp"
#include "sim_patterns.hpp"
#include <kitty/partial_truth_table.hpp>
#include <kitty/print.hpp>
#include <utility>
#include <map>
#include <set>


namespace mockturtle
{

template<typename Ntk>
class simulation_view : public Ntk
{
public:
  using storage = typename Ntk::storage;
  using node    = typename Ntk::node;
  using signal  = typename Ntk::signal;

public:
  simulation_view( Ntk & ntk )
    : Ntk(ntk),
    nodes_to_patterns( *this ),
    nodes_to_size_fanin( *this ),
    nodes_to_layer( *this )
  {
    static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
    static_assert( has_compute_v<Ntk, kitty::dynamic_truth_table>, "Ntk does not implement the compute function" );
    static_assert( has_create_pi_v<Ntk>, "Ntk does not implement the create_pi function" );
    static_assert( has_create_po_v<Ntk>, "Ntk does not implement the create_po function" );
    static_assert( has_create_node_v<Ntk>, "Ntk does not implement the create_node function" );
    static_assert( has_create_and_v<Ntk>, "Ntk does not implement the create_and function" );
    static_assert( has_create_or_v<Ntk>, "Ntk does not implement the create_or function" );
    static_assert( has_create_lt_v<Ntk>, "Ntk does not implement the create_lt function" );
    static_assert( has_create_le_v<Ntk>, "Ntk does not implement the create_le function" );
    static_assert( has_create_xor_v<Ntk>, "Ntk does not implement the create_xor function" );
    static_assert( has_create_maj_v<Ntk>, "Ntk does not implement the create_maj function" );
    static_assert( has_create_ite_v<Ntk>, "Ntk does not implement the create_ite function" );
    static_assert( has_get_constant_v<Ntk>, "Ntk does not implement the get_constant function" );
    static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node function" );
    static_assert( has_foreach_pi_v<Ntk>, "Ntk does not implement the foreach_pi function" );
    static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po function" );
    static_assert( has_foreach_gate_v<Ntk>, "Ntk does not implement the foreach_gate function" );
    static_assert( has_make_signal_v<Ntk>, "Ntk does not implement the make_signal function" );
    static_assert( has_num_pis_v<Ntk>, "Ntk does not implement the num_pis function" );
    
    signal s0 = Ntk::get_constant( false );
    signal s1 = Ntk::get_constant( true );
    kitty::partial_truth_table tt0;
    sim_pattern<Ntk> pat0 = sim_pattern<Ntk>( tt0, s0 );
    sim_pattern<Ntk> pat1 = sim_pattern<Ntk>( tt0, s1 );
    nodes_to_patterns[ Ntk::get_node(s0) ] = sim_patterns.size();
    sim_patterns.push_back( pat0 );
    nodes_to_patterns[ Ntk::get_node(s1) ] = sim_patterns.size();
    sim_patterns.push_back( pat1 );

  }

  signal create_pi( kitty::partial_truth_table pat, std::string const& name = std::string() )
  {
    signal s = Ntk::create_pi( name );
    sim_pattern<Ntk> input_pat = sim_pattern<Ntk>( pat, s, true );
    nodes_to_patterns[Ntk::get_node(s)] = sim_patterns.size();
    sim_patterns.push_back( input_pat );
    input_patterns.push_back( input_pat );
    
    if( layer_to_signals.size() == 0 )
    {
      layer_to_signals.push_back(std::vector{s});
      summary.push_back(std::to_string( s ));
    }
    else
    {
      layer_to_signals[0].push_back(s);
      summary[0] = summary[0] + " " + std::to_string( s ); 
    }
    
    nodes_to_size_fanin[Ntk::get_node(s)] = 0;
    sim_patterns[get_node_pattern( s )].flag_sized = true;

    nodes_to_layer[Ntk::get_node(s)] = 0;
    
    return s;
  }

  std::vector<sim_pattern<Ntk>> get_input_patterns()
  {
    return input_patterns;
  }

#pragma region Create unary functions
  signal create_not( signal const& a )
  {
    return _create_node( {a}, 3 );
  }
  #pragma endregion

#pragma region Create binary functions
  signal create_and( signal a, signal b )
  {
    return _create_node( {a, b}, 4 );
  }

  signal create_nand( signal a, signal b )
  {
    return _create_node( { a, b }, 5 );
  }

  signal create_or( signal a, signal b )
  {
    return _create_node( {a, b}, 6 );
  }

  signal create_lt( signal a, signal b )
  {
    return _create_node( {a, b}, 8 );
  }

  signal create_le( signal a, signal b )
  {
    return _create_node( {a, b}, 11 );
  }

  signal create_xor( signal a, signal b )
  {
    return _create_node( {a, b}, 12 );
  }
#pragma endregion

#pragma region Create ternary functions
signal create_maj( signal a, signal b, signal c )
  {
    return _create_node( {a, b, c}, 14 );
  }

  signal create_ite( signal a, signal b, signal c )
  {
    return _create_node( {a, b, c}, 16 );
  }

  signal create_xor3( signal a, signal b, signal c )
  {
    return _create_node( {a, b, c}, 18 );
  }
#pragma endregion

#pragma region Create arbitrary functions
  signal create_node( std::vector<signal> const& children, kitty::dynamic_truth_table const& function )
  {
    if ( children.size() == 0u )
    {
      assert( function.num_vars() == 0u );
      return Ntk::get_constant( !kitty::is_const0( function ) );
    }

    return _create_node( children, Ntk::_storage->data.cache.insert( function ) );
  }

  signal _create_node( std::vector<signal> const& children, uint32_t literal )
  {
    signal f0 = Ntk::_create_node( children, literal );

    std::vector<kitty::partial_truth_table> fanins_sim;
    uint count_fanin_gates = 0;
    std::string summary_entry = "{ " + std::to_string( f0 ) + ": ";
    for( size_t i = 0; i < children.size(); ++i )
    {
      if( Ntk::is_pi(children[i]) )
        fanins_sim.push_back( input_patterns[ get_input_pattern(children[i]) ].pat );
      else
        fanins_sim.push_back( sim_patterns[ get_node_pattern(children[i]) ].pat );
      
      count_fanin_gates += nodes_to_size_fanin[get_node(children[i])];
      summary_entry = summary_entry + std::to_string( children[i] ) + " ";
    }
    summary_entry = summary_entry + kitty::to_binary( Ntk::_storage->data.cache[literal] ) + " } ";
        
    nodes_to_size_fanin[get_node(f0)] = count_fanin_gates+1;

    kitty::partial_truth_table new_pat = Ntk::compute( Ntk::get_node(f0), fanins_sim.begin(), fanins_sim.end() );
    nodes_to_patterns[f0] = sim_patterns.size();

    sim_pattern<Ntk> new_sim_pat = sim_pattern<Ntk>( new_pat, f0 );
    
    sim_patterns.push_back( new_sim_pat );
    
    uint32_t new_layer = 0;
    for( uint32_t i = 0; i < children.size(); ++i )
    {
      if( nodes_to_layer[children[i]] > new_layer )
        new_layer =  nodes_to_layer[children[i]];
    }
    new_layer++;
    
    if( layer_to_signals.size() == new_layer )
    {
      layer_to_signals.push_back(std::vector{f0});
      summary.push_back( summary_entry );
    }
    else
    {
      layer_to_signals[new_layer].push_back(f0);
      summary[new_layer]+= summary_entry;
    }

    nodes_to_layer[f0]=new_layer;
    return f0;
  }

  signal clone_node( klut_network const& other, node const& source, std::vector<signal> const& children )
  {
    assert( !children.empty() );
    const auto tt = other._storage->data.cache[other._storage->nodes[source].data[1].h1];
    return create_node( children, tt );
  }
#pragma endregion

#pragma region initialization
void initialize_network( std::vector<kitty::partial_truth_table> const& examples )
{
  clear_simulated();
  std::vector<sim_pattern<Ntk>> void_patterns;

  input_patterns = void_patterns;
  sim_patterns = void_patterns;

  signal s0 = Ntk::get_constant( false );
  signal s1 = Ntk::get_constant( true );
  kitty::partial_truth_table tt0(examples[0].num_bits());
  sim_pattern<Ntk> pat0 = sim_pattern<Ntk>( tt0, s0, true );
  sim_pattern<Ntk> pat1 = sim_pattern<Ntk>( ~tt0, s1, true );
  nodes_to_patterns[ Ntk::get_node(s0) ] = sim_patterns.size();
  sim_patterns.push_back( pat0 );
  nodes_to_patterns[ Ntk::get_node(s1) ] = sim_patterns.size();
  sim_patterns.push_back( pat1 );

  if(examples.size() == Ntk::num_pis())
  {
    Ntk::foreach_pi( [&]( auto const& n, auto i ) {
      signal s = Ntk::make_signal(n);
      sim_pattern<Ntk> input_pat = sim_pattern<Ntk>( examples[i], s, true );
      nodes_to_patterns[n] = sim_patterns.size();
      sim_patterns.push_back( input_pat );
      input_patterns.push_back( input_pat );
    } );
  }
  else
  {
    for( size_t i = 0; i < examples.size(); ++i )
      create_pi(examples[i]);
  }

  Ntk::foreach_gate( [&]( auto const& n ) {
    if( !Ntk::is_pi(n) )
      simulate_fanin_cone( n );
} );
}
#pragma endregion

#pragma region simulation
void simulate_fanin_cone( node const& n )
{
  uint32_t child_layer;
  uint32_t max_child_layer = 0;

  if( ! Ntk::is_pi( n ) )
  {
  if( !( sim_patterns[nodes_to_patterns[Ntk::get_node(Ntk::get_constant(false))]].simulated ) )
  {
    uint32_t num_bits = input_patterns[0].pat.num_bits();
    kitty::partial_truth_table tt0( num_bits );
    sim_patterns[nodes_to_patterns[Ntk::get_node(Ntk::get_constant(false))]].pat = tt0;
    sim_patterns[nodes_to_patterns[Ntk::get_node(Ntk::get_constant(true))]].pat = ~tt0;

    sim_patterns[nodes_to_patterns[Ntk::get_node(Ntk::get_constant(false))]].simulated = true;
    sim_patterns[nodes_to_patterns[Ntk::get_node(Ntk::get_constant(true))]].simulated = true;
  }

  if( !nodes_to_patterns.has( n ) || ( nodes_to_patterns.has( n ) && !sim_patterns[nodes_to_patterns[n]].simulated ) )
  {

    std::vector<kitty::partial_truth_table> fanin_values( Ntk::fanin_size( n ) );
    Ntk::foreach_fanin( n, [&]( auto const& f, auto i ) {
      node nchild = Ntk::get_node( f );

      if ( nodes_to_patterns.has( nchild ) && sim_patterns[get_node_pattern(nchild)].simulated )
      {
        uint32_t npat_idx = get_node_pattern(f);
        fanin_values[i] = sim_patterns[ npat_idx ].pat;
      }
      else if ( !nodes_to_patterns.has( nchild ) || (nodes_to_patterns.has( nchild ) && !sim_patterns[get_node_pattern(nchild)].simulated ) )
      {
        simulate_fanin_cone( nchild );
        uint32_t npat_idx = get_node_pattern(f);
        fanin_values[i] = sim_patterns[npat_idx].pat;
      }
      else
      {
        std::cerr << "not valid " << std::endl; 
      }

      child_layer = sim_patterns[get_node_pattern(nchild)].layer;
      if( child_layer > max_child_layer )
        max_child_layer = child_layer;
    } );


    kitty::partial_truth_table pat = Ntk::compute( n, fanin_values.begin(), fanin_values.end() );

    if( !nodes_to_patterns.has( n ) || (nodes_to_patterns.has( n ) && !sim_patterns[get_node_pattern(n)].simulated) )
    {
      nodes_to_patterns[n] = sim_patterns.size();
      sim_pattern<Ntk> spat = sim_pattern<Ntk>( pat, Ntk::make_signal( n ) ); 
      sim_patterns.push_back( spat );
    }
    else
    {
      sim_pattern<Ntk> spat = sim_pattern<Ntk>( pat, Ntk::make_signal( n ) ); 
      sim_patterns[nodes_to_patterns[n]] = spat;
    }
  }
  }
  sim_patterns[nodes_to_patterns[n]].simulated = true;
  sim_patterns[nodes_to_patterns[n]].layer = max_child_layer+1;
  
  if( layer_to_signals.size() == max_child_layer+1 )
    layer_to_signals.push_back( std::vector{ Ntk::make_signal(n) } );
  else if( layer_to_signals.size() > max_child_layer+1 )
    layer_to_signals[max_child_layer+1].push_back(Ntk::make_signal(n));
  else
    std::cerr << "[e] problems in layers filling" << std::endl;

}

void simulate_network()
{
  Ntk::foreach_po( [&]( auto f ) {
    node n = Ntk::get_node(f);
    simulate_fanin_cone( n );
  } );
}

void simulate_fanin_cone_explicit(  std::vector<sim_pattern<Ntk>> & dest_patterns, std::vector<sim_pattern<Ntk>> & examples, 
                                    unordered_node_map<uint32_t, Ntk> & nodes_to_tmp_patterns ,node const& n )
{

  if( !Ntk::is_pi( n ) )
  {
    std::vector<kitty::partial_truth_table> fanin_values( Ntk::fanin_size( n ) );

    Ntk::foreach_fanin( n, [&]( auto const& f, auto i ) {
      node nchild = Ntk::get_node( f );
      
      if ( Ntk::is_pi( nchild ) || ( f == Ntk::get_constant(false)) || ( f == Ntk::get_constant(true)) )
        fanin_values[i] = examples[ nodes_to_tmp_patterns[f] ].pat;
      else
      {
        simulate_fanin_cone_explicit( dest_patterns, examples, nodes_to_tmp_patterns , nchild );
        fanin_values[i] = examples[ nodes_to_tmp_patterns[f] ].pat;
      }
    } );

    kitty::partial_truth_table pat = Ntk::compute( n, fanin_values.begin(), fanin_values.end() );

    if ( nodes_to_tmp_patterns.has( n ) )
    {
      examples[nodes_to_tmp_patterns[n]].pat = pat;
      examples[nodes_to_tmp_patterns[n]].sig = Ntk::make_signal( n );
      if( (examples[nodes_to_tmp_patterns[n]].pat != examples[nodes_to_tmp_patterns[0]].pat)&&(examples[nodes_to_tmp_patterns[n]].pat != examples[nodes_to_tmp_patterns[1]].pat) )
        dest_patterns.push_back(examples[nodes_to_tmp_patterns[n]]);
    }
    else
    {
      nodes_to_tmp_patterns[n] = examples.size();
      sim_pattern<Ntk> spat = sim_pattern<Ntk>( pat, Ntk::make_signal( n ) ); 
      examples.push_back( spat );
      if( (examples[nodes_to_tmp_patterns[n]].pat != examples[nodes_to_tmp_patterns[0]].pat)&&(examples[nodes_to_tmp_patterns[n]].pat != examples[nodes_to_tmp_patterns[1]].pat) )
        dest_patterns.push_back(examples[nodes_to_tmp_patterns[n]]);
    }
  }
  examples[nodes_to_tmp_patterns[n]].simulated = true;
}
#pragma endregion simulation

#pragma region clear flags
void clear_flag()
{
  for( size_t i = 0; i < sim_patterns.size(); ++i )
    sim_patterns[i].flag = false;
}
void clear_simulated()
{
  for( size_t i = 0; i < sim_patterns.size(); ++i )
    sim_patterns[i].simulated = false;
}
void clear_weight()
{
  for( size_t i = 0; i < sim_patterns.size(); ++i )
    sim_patterns[i].weight = -1;
}
#pragma endregion flags

#pragma region count fanin nodes
void clear_network_fanin_size_from_node( node n )
{
  if( sim_patterns[ n ].flag_sized == false )
  {
    Ntk::foreach_fanin( n, [&]( auto const& f ) {
      node nchild = get_node( f );
      if( nodes_to_size_fanin[nchild] != 0 )
      {
        //nodes_to_size_fanin[nchild] = 0;
        clear_network_fanin_size_from_node( nchild );
      }
    } );
    nodes_to_size_fanin[n] = 0;
    sim_patterns[ n ].flag_sized = true;
  }
}

void update_network_fanin_size()
{
  Ntk::foreach_gate( [&]( auto const& n) {
    if( !sim_patterns[n].flag_sized )
    {
      uint32_t counter = 1;
      Ntk::foreach_fanin( n, [&]( auto const& f) {
        node nchild = get_node( f );
        counter+=nodes_to_size_fanin[nchild]; 
      } );
      nodes_to_size_fanin[n] = counter;
    }
  } );
}

#pragma endregion count fanin nodes


uint32_t get_node_pattern( signal s )
{ 
  node n = Ntk::get_node( s );
  return nodes_to_patterns[n];
}

uint32_t get_node( signal s )
{ 
  return Ntk::get_node( s );
}

uint32_t get_input_pattern( signal s )
{ 
  node n = Ntk::get_node( s );
  return nodes_to_patterns[n]-2;
}

signal get_constant( bool bool_val = false )
{ 
  return Ntk::get_constant( bool_val );
}

uint32_t is_pi( signal s )
{ 
  return Ntk::is_pi( s );
}

public:
  std::vector<sim_pattern<Ntk>> input_patterns;
  unordered_node_map<uint32_t, Ntk> nodes_to_patterns;
  unordered_node_map<uint32_t, Ntk> nodes_to_size_fanin;
  std::vector<sim_pattern<Ntk>> sim_patterns;

  std::vector<std::vector<signal>> layer_to_signals;
  unordered_node_map<uint32_t, Ntk> nodes_to_layer;

public:
  uint32_t layer_pointer = 0;
  unsigned seed = 336;
  std::vector<std::string> summary;
  std::vector<kitty::partial_truth_table> targets;
  std::set<std::pair<std::vector<signal>,std::string>> available_nodes;

};

} /* namespace mockturtle */
