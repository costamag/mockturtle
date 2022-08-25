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
  \file graph_to_lfe.hpp
  \brief Convert a network into a vector of input-output simulation patterns.
  \author Andrea Costamagna
*/

#pragma once

#include "../../traits.hpp"

#include <vector>
#include <utility>      // std::pair, std::make_pair
#include <string>       // std::string
#include <iostream>     // std::cout

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>
#include <kitty/bit_operations.hpp>
#include <kitty/operators.hpp>
#include <kitty/print.hpp>

#include "mockturtle/algorithms/simulation.hpp"


namespace mockturtle
{
using ClfeNtk = std::pair<std::vector<kitty::dynamic_truth_table>, std::vector<kitty::dynamic_truth_table>>;
using PlfeNtk = std::pair<std::vector<kitty::partial_truth_table>, kitty::partial_truth_table>;
using PlfeNtk_MO = std::pair<std::vector<kitty::partial_truth_table>, std::vector<kitty::partial_truth_table>>;


template<class Ntk>
struct lfeNtk
{
  ClfeNtk complete;
  PlfeNtk partial;
  PlfeNtk_MO partial_MO;
  std::vector<signal<Ntk>> signals;
  bool use_MO{false};
  uint64_t oidx{0};
};
namespace detail
{

template<class Ntk>
class graph_to_lfe_impl
{
public:
  graph_to_lfe_impl( Ntk& ntk, int oidx )
      : _ntk( ntk ),
      _oidx( oidx )
  {
  }

  lfeNtk<Ntk> run( )
  {
    
    default_simulator<kitty::dynamic_truth_table> sim( _ntk.num_pis() );
    unordered_node_map<kitty::dynamic_truth_table, Ntk> node_to_value( _ntk );
    simulate_nodes<kitty::dynamic_truth_table>( _ntk, node_to_value, sim );
    
    lfeNtk<Ntk> result;
    
    #pragma region complete
    std::vector<kitty::dynamic_truth_table> xs;
    size_t i = 0;
    _ntk.foreach_pi( [&]( auto const& node, auto index  ) {
      xs.push_back( node_to_value[ node ] );
      result.signals.push_back( _ntk._storage->inputs[i++] );
    } );

    std::vector<kitty::dynamic_truth_table> gs;
    std::vector<kitty::dynamic_truth_table> fs;

    _ntk.foreach_po( [&]( auto const& node, auto index ) {
      _ntk._storage->nodes[node].data[0].h2 = 1;
    } );

    _ntk.foreach_gate( [&]( auto const& node, auto index  ) {
      if( _ntk._storage->nodes[node].data[0].h2 == 0 )
      {
        gs.push_back( node_to_value[ node ] );
        result.signals.push_back( _ntk._storage->hash[_ntk._storage->nodes[node]] );
      }
    } );

    _ntk.foreach_po( [&]( auto const& node, auto index ) {
      fs.push_back( node_to_value[ node ] );
    } );

    for( const auto& x : gs )
      xs.push_back(x);

    result.complete = std::make_pair(xs,fs);
    #pragma endregion complete

    #pragma region partial_MO
    std::vector<kitty::partial_truth_table> xs_MO;
    std::vector<kitty::partial_truth_table> fs_MO;
    for( auto x : xs )
    {
      kitty::partial_truth_table xs_tmp( xs[0].num_bits() );
      kitty::create_from_binary_string( xs_tmp, kitty::to_binary( x ) );
      xs_MO.push_back( xs_tmp );
    }

    for( auto f : fs )
    {
      kitty::partial_truth_table fs_tmp( fs[0].num_bits() );
      kitty::create_from_binary_string( fs_tmp, kitty::to_binary( f ) );
      fs_MO.push_back( fs_tmp );
    }

    result.partial_MO = std::make_pair(xs_MO,fs_MO);

    #pragma endregion partial_MO

    #pragma region partial
    std::vector<kitty::partial_truth_table> xsP;
    if( _oidx < 0 )
    {
      for( auto i = 0; i < xs.size(); ++i )
      {
        xsP.emplace_back();
        for( auto j = 0; j < fs.size(); ++j )
        {
          for( auto k = 0; k < xs[i].num_bits(); ++k )
          {
            xsP[i].add_bit( kitty::get_bit( xs[i], k) ); // xsP[size0 + i]
          }
        }
      }
    }
    else
    {
      for( auto i = 0; i < xs.size(); ++i )
      {
        xsP.emplace_back();
        for( auto k = 0; k < xs[i].num_bits(); ++k )
        {
          xsP[i].add_bit( kitty::get_bit( xs[i], k) ); // xsP[size0 + i]
        }
      }      
    }

    kitty::partial_truth_table fsP;
    if( _oidx < 0 )
    {
      for( auto i = 0; i < fs.size(); ++i )
      {
        for( auto j = 0; j < fs[0].num_bits(); ++j )
        {
          fsP.add_bit( kitty::get_bit( fs[i], j) );
        }
      }
    }
    else
    {
      for( auto j = 0; j < fs[0].num_bits(); ++j )
      {
        fsP.add_bit( kitty::get_bit( fs[_oidx], j) );
      }
    }
    result.partial = std::make_pair(xsP,fsP);

    #pragma endregion partial

    return result;
  }

private:
  Ntk& _ntk;
  int _oidx;

};
} // namespace detail

/*! \brief Convert a network into LFE, i.e., a set of input-output simulation patterns
 *
 * comments
 *
 * \tparam Ntk Type of the origin network..
 * \return An equivalent LFE simulation patterns
 */
template<class Ntk>
lfeNtk<Ntk> graph_to_lfe( Ntk& ntk, int oidx = -1 )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_get_constant_v<Ntk>, "Ntk does not implement the get_constant method" );
  static_assert( has_constant_value_v<Ntk>, "Ntk does not implement the constant_value method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_pi_v<Ntk>, "Ntk does not implement the foreach_pi method" );
  static_assert( has_foreach_gate_v<Ntk>, "Ntk does not implement the foreach_gate method" );
  static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin method" );
  static_assert( has_fanin_size_v<Ntk>, "Ntk does not implement the fanin_size method" );
  static_assert( has_num_pos_v<Ntk>, "Ntk does not implement the num_pos method" );
  static_assert( has_compute_v<Ntk, kitty::dynamic_truth_table>, "Ntk does not implement the compute method for SimulationType" );

  detail::graph_to_lfe_impl<Ntk> impl( ntk, oidx );
  return impl.run();
}

} // namespace mockturtle