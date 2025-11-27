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
  \file struct_optimizer.hpp
  \brief Inplace rewrite

  \author Andrea Costamagna
*/

#pragma once

#include <fmt/format.h>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>
#include "../../algorithms/cnf.hpp"

namespace mockturtle
{

template<uint32_t num_vars, uint32_t num_vars_l>
class truth_table_list_t
{
public:
  using lit_t = uint32_t;
  uint32_t null_lit = std::numeric_limits<uint32_t>::max();
  static constexpr uint32_t offset = 2;

  struct node_t
  {
    node_t( lit_t lit ):
      lit( lit )
    {}

    node_t( std::array<lit_t, num_vars> children, kitty::static_truth_table<num_vars> tt, lit_t lit ):
      children( children ), func( tt ), care( tt | ~tt ), lit( lit )
    {}

    node_t( std::array<lit_t, num_vars> children, kitty::static_truth_table<num_vars> tt, kitty::static_truth_table<num_vars> mk, lit_t lit ):
      children( children ), func( tt ), care( mk ), lit( lit )
    {}

    lit_t lit;
    std::array<lit_t, num_vars> children; 
    kitty::static_truth_table<num_vars> func;
    kitty::static_truth_table<num_vars> care;
  };

  truth_table_list_t()
    : gen( 5 )
  {
    for ( uint32_t i = 0; i < num_vars_l; ++i )
    {
      nodes.emplace_back( i + offset );
    }
    for ( int i = 0; i < num_vars_l; ++i )
      kitty::create_nth_var( _xl[i], i );
    for ( int i = 0; i < num_vars; ++i )
      kitty::create_nth_var( _xs[i], i );
  }

  /*! \brief iterator to apply a lambda to all the nodes */
  template<typename Fn>
  void foreach_gate( Fn&& fn ) const
  {
    for ( uint64_t i = num_vars_l; i < nodes.size(); i++ )
    {
      fn( nodes[i] );
    }
  }

  bool is_null( lit_t lit )
  {
    return lit == null_lit;
  }


  lit_t delay_oriented_shannon( kitty::static_truth_table<num_vars_l> tt, std::vector<double> times, std::vector<uint32_t> indices )
  {
    if ( kitty::is_const0( tt ) )
      return 0;
    if ( kitty::is_const0( ~tt ) )
      return 1;

    auto supp = kitty::min_base_inplace( tt );
    if ( (supp.size() > indices.size()) || ( supp.back() >= indices.size() ) )
      return null_lit;

    int j = supp.size()-1;
    for ( int i = indices.size() - 1; i >= 0; --i )
    {
      if ( (j >= 0) && (i == static_cast<int>( supp[j])) )
      {
        j -= 1;
      }
      else
      {
        indices.erase( indices.begin() + i );
        times.erase( times.begin() + i );
      }
    }
//std::cout << indices.size() << std::endl;
    /* here the truth table is reduced and the support minimized */
    if ( indices.size() <= num_vars )
    {
      kitty::static_truth_table<num_vars> func;
      kitty::shrink_to_inplace( func, tt );
      std::array<lit_t, num_vars> children;
      for ( int i = 0; i < num_vars; ++i )
      {
        if ( i < indices.size() )
        {
          children[i] = indices[i] + offset;
        }
        else
        {
          children[i] = null_lit;
        }
      }
      lit_t lit = nodes.size() + offset;
      nodes.emplace_back( children, func, lit );
      return lit;
    }
    /* find the index with the worst delay */
    std::vector<uint32_t> sorted;
    for ( int i = 0; i < indices.size(); ++i )
      sorted.push_back( i );

    std::sort( sorted.begin(), sorted.end(), [&]( auto a, auto b ){
      return times[a] > times[b];
    } );
    std::array<uint32_t, num_vars> children;
    std::array<kitty::static_truth_table<num_vars_l>, num_vars> tts;
    auto tt0 = kitty::cofactor0( tt, sorted[0] );
    auto tt1 = kitty::cofactor1( tt, sorted[0] );
    if ( kitty::is_const0( tt0 ) || kitty::is_const0( ~tt0 ) )
    {
      auto lit_new = delay_oriented_shannon( tt1, times, indices );
      if ( lit_new == null_lit )
        return null_lit;
      children = { indices[sorted[0]] + offset, lit_new, null_lit, null_lit };
      tts[0] = _xl[sorted[0]];
      tts[1] = tt1;
      auto [fun, sim, ex] = chatterjee_synthesis( tts, 2, tt );
      auto lit = create_node( children, fun );
      return lit;
    }
    if ( kitty::is_const0( tt1 ) || kitty::is_const0( ~tt1 ) )
    {
      auto lit_new = delay_oriented_shannon( tt0, times, indices );
      if ( lit_new == null_lit )
        return null_lit;
      children = { indices[sorted[0]] + offset, lit_new, null_lit, null_lit };
      tts[0] = _xl[sorted[0]];
      tts[1] = tt0;
      auto [fun, sim, ex] = chatterjee_synthesis( tts, 2, tt );
      auto lit = create_node( children, fun );
      return lit;
    }
    if ( kitty::equal( tt0, ~tt1 ) )
    {
      auto lit_new = delay_oriented_shannon( tt0, times, indices );
      if ( lit_new == null_lit )
        return null_lit;
      children = { indices[sorted[0]] + offset, lit_new, null_lit, null_lit };
      tts[0] = _xl[sorted[0]];
      tts[1] = tt0;
      auto [fun, sim, ex] = chatterjee_synthesis( tts, 2, tt );
      auto lit = create_node( children, fun );
      return lit;
    }

    tts[0] = _xl[sorted[0]];
    tts[1] = tt0;
    tts[2] = tt1;
    children = { indices[sorted[0]] + offset, null_lit, null_lit, null_lit };
    children[1] = delay_oriented_shannon( tt0, times, indices );
    children[2] = delay_oriented_shannon( tt1, times, indices );
    auto [fun, sim, ex] = chatterjee_synthesis( tts, 3, tt );
    auto lit = create_node( children, fun );
    return lit; 
}

  bool minimize_support( kitty::static_truth_table<num_vars_l>& tt, std::vector<double>& times, std::vector<uint32_t>& indices )
  {
    auto supp = kitty::min_base_inplace( tt );
    if ( (supp.size() > indices.size()) || ( supp.back() >= indices.size() ) )
      return false;

    int j = supp.size()-1;
    for ( int i = indices.size() - 1; i >= 0; --i )
    {
      if ( (j >= 0) && (i == static_cast<int>( supp[j])) )
      {
        j -= 1;
      }
      else
      {
        indices.erase( indices.begin() + i );
        times.erase( times.begin() + i );
      }
    }
    return true;
  }

  lit_t create_node( std::array<lit_t, num_vars> const& children, kitty::static_truth_table<num_vars> const& fun )
  {
    if ( kitty::is_const0( fun ) )
      return 0;
    if ( kitty::is_const0( ~fun ) )
      return 1;

    lit_t lit = nodes.size() + offset;
    nodes.emplace_back( children, fun, lit );
    return lit;
  }

  lit_t create_node( std::array<lit_t, num_vars> const& children, kitty::static_truth_table<num_vars> const& func, kitty::static_truth_table<num_vars> const& care )
  {
    if ( kitty::is_const0( care & func ) )
      return 0;
    if ( kitty::equal( care & func, care ) )
      return 1;

    lit_t lit = nodes.size() + offset;
    nodes.emplace_back( children, func, care, lit );
    return lit;
  }


  std::tuple<kitty::static_truth_table<num_vars>, kitty::static_truth_table<num_vars_l>, bool> chatterjee_synthesis( std::array<uint32_t, num_vars> const& vars, kitty::static_truth_table<num_vars_l> const& tt, kitty::static_truth_table<num_vars_l> const& mk )
  {
    kitty::static_truth_table<num_vars_l> sim;
    kitty::static_truth_table<num_vars> fun;

    uint32_t nvars = 0;
    for ( auto v : vars )
    {
      nvars += ( v != null_lit ) ? 1 : 0;
    }

    uint32_t num_minterms = 1u << nvars;
    kitty::static_truth_table<num_vars_l> tmp_l;
    kitty::static_truth_table<num_vars_l> flip_l0;
    kitty::static_truth_table<num_vars_l> flip_l1;
    kitty::static_truth_table<num_vars> tmp_s;
    kitty::static_truth_table<num_vars> flip_s0;
    kitty::static_truth_table<num_vars> flip_s1;
    uint32_t num0_best = 0;
    uint32_t num1_best = 0;
    std::bernoulli_distribution d( 0.5 );
    for ( int m = 0; m < num_minterms; ++m )
    {
      tmp_l |= ~tmp_l;
      tmp_s |= ~tmp_s;
      int k = 0;
      for ( int v = 0; v < num_vars; ++v )
      {
        if ( vars[v] != null_lit )
        {
          tmp_l &= ( ( ( m >> k ) & 0x1 ) == 0x1 ) ? _xl[vars[v]] : ~_xl[vars[v]];
          tmp_s &= ( ( ( m >> k ) & 0x1 ) == 0x1 ) ? _xs[k] : ~_xs[k];
          k++;
        }
      }

      uint32_t num_1 = kitty::count_ones( mk & tmp_l & tt );
      uint32_t num_0 = kitty::count_ones( mk & tmp_l & ~tt );

      if ( num_1 > num1_best )
      {
        flip_l1 = tmp_l;
        flip_s1 = tmp_s;
        num1_best = num_1;
      }
      if ( num_0 > num1_best )
      {
        flip_l0 = tmp_l;
        flip_s0 = tmp_s;
        num1_best = num_0;
      }

      bool set_to_one = ( num_1 > num_0 ) | ( ( num_1 == num_0 ) && d( gen ) );

      if ( set_to_one ) 
      {
        sim |= tmp_l;
        fun |= tmp_s;
      }
      else
      {
        sim &= ~tmp_l;
        fun &= ~tmp_s;
      }
    }

    if ( kitty::is_const0( fun ) )
    {
      sim |= flip_l1;
      fun |= flip_s1;
    }
    else if ( kitty::is_const0( ~fun ) )
    {
      sim &= ~flip_l0;
      fun &= ~flip_s0;
    }

    bool is_same = kitty::equal( mk&sim, mk&tt );

    return std::make_tuple( fun, sim, is_same );
  }

  std::tuple<kitty::static_truth_table<num_vars>, kitty::static_truth_table<num_vars_l>, bool> close_synthesis( std::array<kitty::static_truth_table<num_vars_l>, num_vars> const& itts, uint32_t nvars, spfd_manager_t<num_vars_l, num_vars_l> const& spfd )
  {
    int j = 0;
  //  for ( auto t : itts )
  //  {
  //    kitty::print_binary( t ); std::cout << " " << j++ << std::endl;
  //  }
    auto const& func = spfd.func[1];
    auto const& cares = spfd.masks;

    kitty::static_truth_table<num_vars> fun;
    kitty::static_truth_table<num_vars_l> sim;
    bool is_same = true;

    uint32_t constexpr num_minterms = static_cast<uint32_t>( 1u << num_vars );
    uint32_t constexpr num_patterns = static_cast<uint32_t>( 1u << num_vars_l );

    std::array<uint32_t, num_patterns> simbit_to_min;
    std::array<kitty::static_truth_table<num_vars_l>, num_minterms> min_to_sim;
    kitty::static_truth_table<num_vars_l> tmp_l;
    kitty::static_truth_table<num_vars> tmp_s;

    for ( int m {0}; m < ( 1u << nvars ); ++m )
    {
      tmp_l |= ~tmp_l;
      tmp_s |= ~tmp_s;
      uint32_t mint = m;
      for ( uint32_t v = 0; v < nvars; ++v )
      {
        tmp_l &= ( ( ( mint >> v ) & 0x1 ) == 0x1 ) ? itts[v] : ~itts[v];
        tmp_s &= ( ( ( mint >> v ) & 0x1 ) == 0x1 ) ? _xs[v] : ~_xs[v];
      }
      min_to_sim[m] |= tmp_l;
      for ( int b = 0; b < num_patterns; ++b )
      {
        if ( kitty::get_bit( tmp_l, b ) == 1 )
          simbit_to_min[b] = mint;
      }
    }

    kitty::static_truth_table<num_vars> fun_s;
    kitty::static_truth_table<num_vars> done;
    
    for ( uint32_t m = 0; m < (1u << nvars); ++m )
    {
      tmp_l |= ~tmp_l;
      tmp_s |= ~tmp_s;
      for ( int v = 0; v < nvars; ++v )
      {
        tmp_l &= ( ( ( m >> v ) & 0x1 ) == 0x1 ) ? itts[v] : ~itts[v];
        tmp_s &= ( ( ( m >> v ) & 0x1 ) == 0x1 ) ? _xs[v] : ~_xs[v];
      }

      for ( int k = 0; k < spfd.num_masks; ++k )
      {
        auto const& mk = spfd.masks[k];

        if ( kitty::count_ones( mk & min_to_sim[m] ) > 1 )
        {
          for ( int b = 0; b < num_patterns; ++b )
          {
            if ( ( kitty::get_bit( mk, b ) == 1 ) && ( kitty::get_bit( done, simbit_to_min[b] ) == 1 ) && ( simbit_to_min[b] != m ) )
            {
              int v1 = kitty::get_bit( mk & func, b );
              int v2 = kitty::count_ones( mk & func & tmp_l ) > 0 ? 1 : 0;
              if ( v1 != v2 )
              {
                if ( kitty::get_bit( fun_s, simbit_to_min[b] ) == 0 )
                {
                  kitty::set_bit( fun_s, m );
                  fun |= tmp_s;
                  sim |= tmp_l;
                }
                else
                {
                  kitty::clear_bit( fun_s, m );
                  fun &= ~tmp_s;
                  sim &= ~tmp_l;
                }
              }
              else
              {
                if ( kitty::get_bit( fun, simbit_to_min[b] ) == 1 )
                {
                  kitty::set_bit( fun_s, m );
                  fun |= tmp_s;
                  sim |= tmp_l;
                }
                else
                {
                  kitty::clear_bit( fun_s, m );
                  fun &= ~tmp_s;
                  sim &= ~tmp_l;
                }
              }
            }
          }
        }
      }
      kitty::set_bit( done, m );
    }

    auto res = std::make_tuple( fun, sim, is_same );
    return res;
  }

  std::tuple<kitty::static_truth_table<num_vars>, kitty::static_truth_table<num_vars>, bool> chatterjee_synthesis( std::array<kitty::static_truth_table<num_vars_l>, num_vars> const& tts, uint32_t nvars, kitty::static_truth_table<num_vars_l> const& tt, kitty::static_truth_table<num_vars_l> const& mk )
  {
    kitty::static_truth_table<num_vars> mas;
    kitty::static_truth_table<num_vars> fun;

    uint32_t num_minterms = 1u << nvars;
    kitty::static_truth_table<num_vars_l> tmp_l;
    kitty::static_truth_table<num_vars_l> flip_l0;
    kitty::static_truth_table<num_vars_l> flip_l1;
    uint32_t num0_best = 0;
    uint32_t num1_best = 0;
    
    kitty::static_truth_table<num_vars> tmp_s;
    kitty::static_truth_table<num_vars> flip_s0;
    kitty::static_truth_table<num_vars> flip_s1;
    std::bernoulli_distribution d( 0.5 );
    for ( int m = 0; m < num_minterms; ++m )
    {
      tmp_l |= ~tmp_l;
      tmp_s |= ~tmp_s;
      for ( int v = 0; v < nvars; ++v )
      {
        tmp_l &= ( ( ( m >> v ) & 0x1 ) == 0x1 ) ? tts[v] : ~tts[v];
        tmp_s &= ( ( ( m >> v ) & 0x1 ) == 0x1 ) ? _xs[v] : ~_xs[v];
      }

      uint32_t num_1 = kitty::count_ones( tmp_l & tt & mk );
      uint32_t num_0 = kitty::count_ones( tmp_l & ~tt & mk );

      bool set_to_one = ( num_1 > num_0 ) || ( ( num_1 == num_0 ) && d( gen ) );
//std::cout << "M=" << m << " " << num_0 << " " << num_1 << std::endl;
      if ( num_0 > num0_best )
      {
        flip_l0 = tmp_l;
        flip_s0 = tmp_s;
        num0_best = num_0;
      }
      if ( num_1 > num1_best )
      {
        flip_l1 = tmp_l;
        flip_s1 = tmp_s;
        num1_best = num_1;
      }

      if ( kitty::count_ones( tmp_l ) > 0 )
        mas |= tmp_s;

      if ( set_to_one ) 
      {
        fun |= tmp_s;
      }
      else
      {
        fun &= ~tmp_s;
      }
    }

    if ( kitty::is_const0( fun ) )
    {
      fun |= flip_s1;
    }
    else if ( kitty::is_const0( ~fun ) )
    {
      fun &= ~flip_s0;
    }


    return std::make_tuple( fun, mas, true );
  }


  std::tuple<kitty::static_truth_table<num_vars>, kitty::static_truth_table<num_vars_l>> chatterjee_synthesis_with_sim( std::array<kitty::static_truth_table<num_vars_l>, num_vars> const& tts, uint32_t nvars, kitty::static_truth_table<num_vars_l> const& tt, kitty::static_truth_table<num_vars_l> const& mk )
  {
    kitty::static_truth_table<num_vars_l> sim;
    kitty::static_truth_table<num_vars> fun;

    uint32_t num_minterms = 1u << nvars;
    kitty::static_truth_table<num_vars_l> tmp_l;
    kitty::static_truth_table<num_vars_l> flip_l0;
    kitty::static_truth_table<num_vars_l> flip_l1;
    uint32_t num0_best = 0;
    uint32_t num1_best = 0;
    
    kitty::static_truth_table<num_vars> tmp_s;
    kitty::static_truth_table<num_vars> flip_s0;
    kitty::static_truth_table<num_vars> flip_s1;
    std::bernoulli_distribution d( 0.5 );
    for ( int m = 0; m < num_minterms; ++m )
    {
      tmp_l |= ~tmp_l;
      tmp_s |= ~tmp_s;
      for ( int v = 0; v < nvars; ++v )
      {
        tmp_l &= ( ( ( m >> v ) & 0x1 ) == 0x1 ) ? tts[v] : ~tts[v];
        tmp_s &= ( ( ( m >> v ) & 0x1 ) == 0x1 ) ? _xs[v] : ~_xs[v];
      }

      uint32_t num_1 = kitty::count_ones( tmp_l & tt & mk );
      uint32_t num_0 = kitty::count_ones( tmp_l & ~tt & mk );
      if ( num_1*num_0 > 0 )
        std::cerr << "[e] Illigal in CH" << std::endl;
      bool set_to_one = ( num_1 > num_0 ) || ( ( num_1 == num_0 ) && d( gen ) );
//std::cout << "M=" << m << " " << num_0 << " " << num_1 << std::endl;
      if ( num_0 > num0_best )
      {
        flip_l0 = tmp_l;
        flip_s0 = tmp_s;
        num0_best = num_0;
      }
      if ( num_1 > num1_best )
      {
        flip_l1 = tmp_l;
        flip_s1 = tmp_s;
        num1_best = num_1;
      }


      if ( set_to_one ) 
      {
        fun |= tmp_s;
        sim |= tmp_l;
      }
      else
      {
        fun &= ~tmp_s;
        sim &= ~tmp_l;
      }
    }

    if ( kitty::is_const0( fun ) )
    {
      fun |= flip_s1;
      sim |= flip_l1;
    }
    else if ( kitty::is_const0( ~fun ) )
    {
      fun &= ~flip_s0;
      sim &= ~flip_l0;
    }


    return std::make_tuple( fun, sim );
  }



  std::tuple<kitty::static_truth_table<num_vars>, kitty::static_truth_table<num_vars_l>, bool> chatterjee_synthesis( std::array<kitty::static_truth_table<num_vars_l>, num_vars> const& tts, uint32_t nvars, kitty::static_truth_table<num_vars_l> const& tt )
  {
    kitty::static_truth_table<num_vars_l> sim;
    kitty::static_truth_table<num_vars> fun;

    uint32_t num_minterms = 1u << nvars;
    kitty::static_truth_table<num_vars_l> tmp_l;
    kitty::static_truth_table<num_vars_l> flip_l0;
    kitty::static_truth_table<num_vars_l> flip_l1;
    uint32_t num0_best = 0;
    uint32_t num1_best = 0;
    
    kitty::static_truth_table<num_vars> tmp_s;
    kitty::static_truth_table<num_vars> flip_s0;
    kitty::static_truth_table<num_vars> flip_s1;
    std::bernoulli_distribution d( 0.5 );
    for ( int m = 0; m < num_minterms; ++m )
    {
      tmp_l |= ~tmp_l;
      tmp_s |= ~tmp_s;
      for ( int v = 0; v < nvars; ++v )
      {
        tmp_l &= ( ( ( m >> v ) & 0x1 ) == 0x1 ) ? tts[v] : ~tts[v];
        tmp_s &= ( ( ( m >> v ) & 0x1 ) == 0x1 ) ? _xs[v] : ~_xs[v];
      }

      uint32_t num_1 = kitty::count_ones( tmp_l & tt );
      uint32_t num_0 = kitty::count_ones( tmp_l & ~tt );

      bool set_to_one = ( num_1 > num_0 ) || ( ( num_1 == num_0 ) && d( gen ) );
//std::cout << "M=" << m << " " << num_0 << " " << num_1 << std::endl;
      if ( num_0 > num0_best )
      {
        flip_l0 = tmp_l;
        flip_s0 = tmp_s;
        num0_best = num_0;
      }
      if ( num_1 > num1_best )
      {
        flip_l1 = tmp_l;
        flip_s1 = tmp_s;
        num1_best = num_1;
      }

      if ( set_to_one ) 
      {
        sim |= tmp_l;
        fun |= tmp_s;
      }
      else
      {
        sim &= ~tmp_l;
        fun &= ~tmp_s;
      }
    }

    if ( kitty::is_const0( fun ) )
    {
      sim |= flip_l1;
      fun |= flip_s1;
    }
    else if ( kitty::is_const0( ~fun ) )
    {
      sim &= ~flip_l0;
      fun &= ~flip_s0;
    }

    bool is_same = kitty::equal( sim, tt );

    return std::make_tuple( fun, sim, is_same );
  }

  std::array<uint32_t, num_vars> find_best_support( std::vector<int> const& indices, kitty::static_truth_table<num_vars_l> const& tt, kitty::static_truth_table<num_vars_l> const& mk )
  {
    std::array<uint32_t, num_vars> supp;
    spfd_manager.init( tt, mk );
    double best_cost = std::numeric_limits<double>::max();
    for ( int i = 0; i < num_vars; ++i )
    {
      supp[i] = null_lit;
      if ( spfd_manager.is_covered() )
        continue;
      for ( int idx : indices )
      {
        double cost = spfd_manager.evaluate( _xl[idx] );
        if ( cost < best_cost )
        {
          best_cost = cost;
          supp[i] = idx;
        }
      }
      if ( supp[i] != null_lit )
        spfd_manager.update( _xl[supp[i]] );
    }
    if ( spfd_manager.is_covered() )
      return supp;
  }

  lit_t bottom_up_synthesis( kitty::static_truth_table<num_vars_l> tt, std::vector<double> times, std::vector<uint32_t> indices )
  {
//std::cout << "A" << std::endl;
    if ( kitty::is_const0( tt ) )
      return 0;
    if ( kitty::is_const0( ~tt ) )
      return 1;

    if( !minimize_support( tt, times, indices ) )
      return null_lit;
//std::cout << "B" << std::endl;

    std::vector<int> sorted_vars;
    for ( int i = 0; i < indices.size(); ++i )
      sorted_vars.push_back( i );
    std::sort( sorted_vars.begin(), sorted_vars.end(), [&]( auto ia, auto ib ){
      return times[ia] < times[ib];
    } );

    /* phase 1 : build the bottom layer */
    /* build the first LUT with the earliest arriving variables */
    std::array<uint32_t, num_vars> children;
    std::array<uint32_t, num_vars> children_indices;
    std::array<kitty::static_truth_table<num_vars_l>, num_vars> xs;
    for ( int i = 0; i < num_vars; ++i )
    {
      children[i] = i < indices.size() ? indices[sorted_vars[i]] + offset : null_lit;
      if ( i < indices.size() )
        xs[i] = _xl[sorted_vars[i]];
    }
    auto [fun1, sim1, ex1] = chatterjee_synthesis( xs, std::min( num_vars, static_cast<uint32_t>(indices.size()) ), tt );
    auto lit1 = create_node( children, fun1 );
    if ( ex1 )
      return lit1;
//std::cout << "C" << std::endl;

    /* cover the information graph with the information */
    spfd_manager.init( tt );
    spfd_manager.update( sim1 );

    kitty::static_truth_table<num_vars_l> target_2;
    target_2 = ( spfd_manager.masks[0] & tt ) | ( spfd_manager.masks[1] & ~tt );

    std::array<std::array<uint32_t, num_vars>, 3u> supp2;
    std::array<kitty::static_truth_table<num_vars_l>, 3u> funcs2 { target_2, tt, tt };
    std::array<kitty::static_truth_table<num_vars_l>, 3u> masks2 { spfd_manager.masks[0]|spfd_manager.masks[1],
                                                                   spfd_manager.masks[0],
                                                                   spfd_manager.masks[1] };

    supp2[0] = find_best_support( sorted_vars, funcs2[0], masks2[0] );
    supp2[1] = find_best_support( sorted_vars, funcs2[1], masks2[1] );
    supp2[2] = find_best_support( sorted_vars, funcs2[2], masks2[2] );
    std::array<kitty::static_truth_table<num_vars>, 3u> funs2;
    std::array<kitty::static_truth_table<num_vars_l>, 3u> sims2;
    std::array<bool, 3u> exs2;
    std::array<std::tuple<kitty::static_truth_table<num_vars>, kitty::static_truth_table<num_vars_l>, bool>, 3> cands2;
    cands2[0] = chatterjee_synthesis( supp2[0], funcs2[0], masks2[0] );
    cands2[1] = chatterjee_synthesis( supp2[1], funcs2[1], masks2[1] );
    cands2[2] = chatterjee_synthesis( supp2[2], funcs2[2], masks2[2] );
    for ( int i = 0; i < 3; ++i )
    {
      funs2[i] = std::get<0>( cands2[i] );
      sims2[i] = std::get<1>( cands2[i] );
      exs2[i] = std::get<2>( cands2[i] );
    }
//std::cout << "D" << std::endl;

    spfd_manager.init( tt );
    spfd_manager.update( sim1 );
    std::array<double, 3u> cost_2;
    for ( int i = 0; i < 3; ++i )
      cost_2[i] = spfd_manager.evaluate( sims2[i] );
// std::cout << cost_2[0] << " " << cost_2[1] << " " << cost_2[2] << std::endl;
    int choice_2 = 0;
//  kitty::print_binary( funs2[0] ); std::cout << std::endl;

    for ( int i = 1; i < 3; ++i )
    {
//    kitty::print_binary( funs2[i] ); std::cout << std::endl;
      if ( !kitty::is_const0(~funs2[i]) && !kitty::is_const0(funs2[i]) && (cost_2[i] < cost_2[choice_2]) )
        choice_2 = i;
    }
//std::cout << "chosen " << choice_2 << std::endl;
    for ( int i = 0; i < num_vars; ++i )
    {
      auto idx = supp2[choice_2][i];
      children[i] = idx != null_lit ? indices[idx] + offset : null_lit;
    }
    auto lit2 = create_node( children, funs2[choice_2] );
    auto sim2 = sims2[choice_2];

    spfd_manager.init( tt );
    spfd_manager.update( sim1 );
    spfd_manager.update( sim2 );

//  std::cout << spfd_manager.num_edges << std::endl;

    if ( spfd_manager.is_covered() )
    {
      std::array<uint32_t, num_vars> children_top { lit1, lit2, null_lit, null_lit };
      std::array<kitty::static_truth_table<num_vars_l>, num_vars> tts;
      tts[0] = sim1;
      tts[1] = sim2;
      auto [funx, simx, exx] = chatterjee_synthesis( tts, 2, tt );
      if ( !exx )
        std::cerr << "[e] simulation not matching" << std::endl;
      auto litx = create_node( children_top, funx );
      return litx;
    }

    /* CONSTRUCTING 3rd TT */
    spfd_manager.init( tt );
    spfd_manager.update( sim1 );
    spfd_manager.update( sim2 );

    std::array<std::array<uint32_t, num_vars>, 22u> supp3;
    auto const& masks_3 = spfd_manager.masks;
    std::array<kitty::static_truth_table<num_vars_l>, 22u> funcs3 { tt,
                                                                   tt, 
                                                                   tt,
                                                                   tt,
                                                                   ( masks_3[0] & ~tt ) | ( masks_3[1] & tt ),
                                                                   ( masks_3[0] & ~tt ) | ( masks_3[2] & tt ),
                                                                   ( masks_3[0] & ~tt ) | ( masks_3[3] & tt ),
                                                                   ( masks_3[1] & ~tt ) | ( masks_3[2] & tt ),
                                                                   ( masks_3[1] & ~tt ) | ( masks_3[3] & tt ),
                                                                   ( masks_3[2] & ~tt ) | ( masks_3[3] & tt ),
                                                                   ( masks_3[0] & ~tt ) | ( masks_3[1] & tt ) | ( masks_3[2] & tt ),
                                                                   ( masks_3[0] & tt ) | ( masks_3[1] & ~tt ) | ( masks_3[2] & tt ),
                                                                   ( masks_3[0] & tt ) | ( masks_3[1] & tt ) | ( masks_3[2] & ~tt ),
                                                                   ( masks_3[0] & ~tt ) | ( masks_3[1] & tt )  | ( masks_3[3] & tt ),
                                                                   ( masks_3[0] & tt )  | ( masks_3[1] & ~tt ) | ( masks_3[3] & tt ),
                                                                   ( masks_3[0] & tt )  | ( masks_3[1] & tt )  | ( masks_3[3] & ~tt ),
                                                                   ( masks_3[0] & ~tt ) | ( masks_3[2] & tt )  | ( masks_3[3] & tt ),
                                                                   ( masks_3[0] & tt )  | ( masks_3[2] & ~tt ) | ( masks_3[3] & tt ),
                                                                   ( masks_3[0] & tt )  | ( masks_3[2] & tt )  | ( masks_3[3] & ~tt ),
                                                                   ( masks_3[1] & ~tt ) | ( masks_3[2] & tt )  | ( masks_3[3] & tt ),
                                                                   ( masks_3[1] & tt )  | ( masks_3[2] & ~tt ) | ( masks_3[3] & tt ),
                                                                   ( masks_3[1] & tt )  | ( masks_3[2] & tt )  | ( masks_3[3] & ~tt )
                                                                    };
    std::array<kitty::static_truth_table<num_vars_l>, 22u> masks3 { masks_3[0],
                                                                   masks_3[1],
                                                                   masks_3[2],
                                                                   masks_3[3],
                                                                   masks_3[0] | masks_3[1],
                                                                   masks_3[0] | masks_3[2],
                                                                   masks_3[0] | masks_3[3],
                                                                   masks_3[1] | masks_3[2],
                                                                   masks_3[1] | masks_3[3],
                                                                   masks_3[2] | masks_3[3],
                                                                   masks_3[0] | masks_3[1] | masks_3[2],
                                                                   masks_3[0] | masks_3[1] | masks_3[2],
                                                                   masks_3[0] | masks_3[1] | masks_3[2],
                                                                   masks_3[0] | masks_3[1] | masks_3[3],
                                                                   masks_3[0] | masks_3[1] | masks_3[3],
                                                                   masks_3[0] | masks_3[1] | masks_3[3],
                                                                   masks_3[0] | masks_3[2] | masks_3[3],
                                                                   masks_3[0] | masks_3[2] | masks_3[3],
                                                                   masks_3[0] | masks_3[2] | masks_3[3],
                                                                   masks_3[1] | masks_3[2] | masks_3[3],
                                                                   masks_3[1] | masks_3[2] | masks_3[3],
                                                                   masks_3[1] | masks_3[2] | masks_3[3]
                                                                  };
    for ( int i = 0; i < 22; ++i )
    {
      supp3[i] = find_best_support( sorted_vars, funcs3[i], masks3[i] );
    }
    std::array<kitty::static_truth_table<num_vars>, 22u> funs3;
    std::array<kitty::static_truth_table<num_vars_l>, 22u> sims3;
    std::array<bool, 22u> exs3;
    std::array<std::tuple<kitty::static_truth_table<num_vars>, kitty::static_truth_table<num_vars_l>, bool>, 22u> cands3;
    for ( int i = 0; i < 22; ++i )
    {
      cands3[i] = chatterjee_synthesis( supp3[i], funcs3[i], masks3[i] );
    }
    for ( int i = 0; i < 22u; ++i )
    {
      funs3[i] = std::get<0>( cands3[i] );
      sims3[i] = std::get<1>( cands3[i] );
      exs3[i] = std::get<2>( cands3[i] );
    }
//std::cout << "D" << std::endl;

    spfd_manager.init( tt );
    spfd_manager.update( sim1 );
    spfd_manager.update( sim2 );
    std::array<double, 22u> cost_3;
    for ( int i = 0; i < 22; ++i )
    {
      cost_3[i] = spfd_manager.evaluate( sims3[i] );
//  std::cout << cost_3[i] << " ";
    }
//  std::cout << std::endl;
    int choice_3 = 0;
    for ( int i = 1; i < 22; ++i )
    {
      if ( cost_3[i] < cost_3[choice_3] )
        choice_3 = i;
    }
//  std::cout << "chosen " << choice_2 << std::endl;
    for ( int i = 0; i < num_vars; ++i )
    {
      auto idx = supp3[choice_3][i];
      children[i] = idx != null_lit ? indices[idx] + offset : null_lit;
    }
    auto lit3 = create_node( children, funs3[choice_3] );
    auto sim3 = sims3[choice_3];

    spfd_manager.init( tt );
    spfd_manager.update( sim1 );
    spfd_manager.update( sim2 );
    spfd_manager.update( sim3 );

    if ( spfd_manager.is_covered() )
    {
      std::array<uint32_t, num_vars> children_top { lit1, lit2, lit3, null_lit };
      std::array<kitty::static_truth_table<num_vars_l>, num_vars> tts;
      tts[0] = sim1;
      tts[1] = sim2;
      tts[2] = sim3;
      auto [funx, simx, exx] = chatterjee_synthesis( tts, 3, tt );
      if ( !exx )
        std::cerr << "[e] simulation not matching" << std::endl;
      auto litx = create_node( children_top, funx );
      return litx;
    }

    spfd_manager.init( tt );
    spfd_manager.update( sim1 );
    spfd_manager.update( sim2 );
    spfd_manager.update( sim3 );
    int nvars = 0;
    double best_cost = spfd_manager.num_edges;
    std::array<kitty::static_truth_table<num_vars_l>, num_vars> tts_end;
    uint32_t best_index = 0;
    while ( ( nvars < num_vars ) && !spfd_manager.is_covered() )
    {
      children[nvars] = null_lit;
      for ( int i = 0; i < sorted_vars.size(); ++i )
      {
        auto const& cost = spfd_manager.evaluate( _xl[sorted_vars[i]] );
        if ( cost < best_cost )
        {
          best_cost = cost;
          children[nvars] = indices[sorted_vars[i]]+offset;
          best_index = sorted_vars[i];
        }
      }
      if ( children[nvars] != null_lit )
      {
        spfd_manager.update( _xl[best_index] );
        tts_end[nvars] = _xl[best_index];
//kitty::print_binary( tts_end[nvars] ); std::cout << " " << nvars << std::endl;
        nvars++;
      }
      else
      {
        break;
      }
    }

    if ( spfd_manager.is_covered() )
    {
      spfd_manager.init( tt );
      spfd_manager.update( sim1 );
      spfd_manager.update( sim2 );
      spfd_manager.update( sim3 );
      
      kitty::static_truth_table<num_vars_l> sim4 = tts_end[0];
      uint32_t lit4 = children[0];


      if ( nvars > 1 )
      {
        return null_lit;
        auto [fun4, sim4_l, ex4] = close_synthesis( tts_end, nvars, spfd_manager );
        sim4 = sim4_l;
        lit4 = create_node( children, fun4 );
      }
      children = { lit1, lit2, lit3, lit4 };
      std::array<kitty::static_truth_table<num_vars_l>, num_vars> tts;
      tts[0] = sim1;
      tts[1] = sim2;
      tts[2] = sim3;
      tts[3] = sim4;
      auto [funx, simx, exx] = chatterjee_synthesis( tts, 4, tt );
//kitty::print_binary( simx ); std::cout << std::endl;
//kitty::print_binary( tt ); std::cout << std::endl;

      if ( !exx )
        std::cerr << "[e] simulation not matching" << std::endl;
      auto litx = create_node( children, funx );
      return litx;
    }
    //std::cout << spfd_manager.num_edges << std::endl;

    return null_lit;
  }

  lit_t delay_oriented_shannon3( kitty::static_truth_table<num_vars_l> tt, std::vector<double> times, std::vector<uint32_t> indices )
  {
    if ( kitty::is_const0( tt ) )
      return 0;
    if ( kitty::is_const0( ~tt ) )
      return 1;

    auto supp = kitty::min_base_inplace( tt );
    if ( (supp.size() > indices.size()) || ( supp.back() >= indices.size() ) )
      return null_lit;

    int j = supp.size()-1;
    for ( int i = indices.size() - 1; i >= 0; --i )
    {
      if ( (j >= 0) && (i == static_cast<int>( supp[j])) )
      {
        j -= 1;
      }
      else
      {
        indices.erase( indices.begin() + i );
        times.erase( times.begin() + i );
      }
    }
//std::cout << indices.size() << std::endl;
    /* here the truth table is reduced and the support minimized */
    if ( indices.size() <= num_vars )
    {
      kitty::static_truth_table<num_vars> func;
      kitty::shrink_to_inplace( func, tt );
      std::array<lit_t, num_vars> children;
      for ( int i = 0; i < num_vars; ++i )
      {
        if ( i < indices.size() )
        {
          children[i] = indices[i] + offset;
        }
        else
        {
          children[i] = null_lit;
        }
      }
      lit_t lit = nodes.size() + offset;
      nodes.emplace_back( children, func, lit );
      return lit;
    }
    /* find the index with the worst delay */
    std::vector<uint32_t> sorted;
    for ( int i = 0; i < indices.size(); ++i )
      sorted.push_back( i );

    std::sort( sorted.begin(), sorted.end(), [&]( auto a, auto b ){
      return times[a] > times[b];
    } );
    std::array<uint32_t, num_vars> children;
    std::array<kitty::static_truth_table<num_vars_l>, num_vars> tts;
    auto tt0 = kitty::cofactor0( tt, sorted[0] );
    auto tt1 = kitty::cofactor1( tt, sorted[0] );
    if ( kitty::is_const0( tt0 ) || kitty::is_const0( ~tt0 ) )
    {
      auto lit_new = delay_oriented_shannon3( tt1, times, indices );
      if ( lit_new == null_lit )
        return null_lit;
      children = { indices[sorted[0]] + offset, lit_new, null_lit, null_lit };
      tts[0] = _xl[sorted[0]];
      tts[1] = tt1;
      auto [fun, sim, ex] = chatterjee_synthesis( tts, 2, tt );
      auto lit = create_node( children, fun );
      return lit;
    }
    if ( kitty::is_const0( tt1 ) || kitty::is_const0( ~tt1 ) )
    {
      auto lit_new = delay_oriented_shannon3( tt0, times, indices );
      if ( lit_new == null_lit )
        return null_lit;
      children = { indices[sorted[0]] + offset, lit_new, null_lit, null_lit };
      tts[0] = _xl[sorted[0]];
      tts[1] = tt0;
      auto [fun, sim, ex] = chatterjee_synthesis( tts, 2, tt );
      auto lit = create_node( children, fun );
      return lit;
    }
    if ( kitty::equal( tt0, ~tt1 ) )
    {
      auto lit_new = delay_oriented_shannon3( tt0, times, indices );
      if ( lit_new == null_lit )
        return null_lit;
      children = { indices[sorted[0]] + offset, lit_new, null_lit, null_lit };
      tts[0] = _xl[sorted[0]];
      tts[1] = tt0;
      auto [fun, sim, ex] = chatterjee_synthesis( tts, 2, tt );
      auto lit = create_node( children, fun );
      return lit;
    }

    auto tt00 = kitty::cofactor0( tt0, sorted[1] );
    auto tt01 = kitty::cofactor1( tt0, sorted[1] );
    auto tt10 = kitty::cofactor0( tt1, sorted[1] );
    auto tt11 = kitty::cofactor1( tt1, sorted[1] );
    auto tt000 = kitty::cofactor0( tt00, sorted[2] );
    auto tt001 = kitty::cofactor0( tt01, sorted[2] );
    auto tt010 = kitty::cofactor0( tt10, sorted[2] );
    auto tt011 = kitty::cofactor0( tt11, sorted[2] );
    auto tt100 = kitty::cofactor1( tt00, sorted[2] );
    auto tt101 = kitty::cofactor1( tt01, sorted[2] );
    auto tt110 = kitty::cofactor1( tt10, sorted[2] );
    auto tt111 = kitty::cofactor1( tt11, sorted[2] );


    std::vector<kitty::static_truth_table<num_vars_l>> reminders3 { tt000, tt001, tt010, tt011, tt100, tt101, tt110, tt111 };
    for ( int i = reminders3.size()-1; i >= 0; --i )
    {
      if ( kitty::is_const0( reminders3[i] ) || kitty::is_const0( ~reminders3[i] ) )
        reminders3.erase( reminders3.begin() + i );
    }
    int i = 0;
    while ( i < reminders3.size() )
    {
      for ( int j = reminders3.size()-1; j > i ; --j )
      {
        if ( kitty::equal( reminders3[i], reminders3[j] ) || kitty::equal( reminders3[i], ~reminders3[j] ) )
          reminders3.erase( reminders3.begin() + j );
      }
      i++;
    }
    std::vector<kitty::static_truth_table<num_vars_l>> reminders2 { tt00, tt01, tt10, tt11 };
    for ( int i = reminders2.size()-1; i >= 0; --i )
    {
      if ( kitty::is_const0( reminders2[i] ) || kitty::is_const0( ~reminders2[i] ) )
        reminders2.erase( reminders2.begin() + i );
    }
    i = 0;
    while ( i < reminders2.size() )
    {
      for ( int j = reminders2.size()-1; j > i ; --j )
      {
        if ( kitty::equal( reminders2[i], reminders2[j] ) || kitty::equal( reminders2[i], ~reminders2[j] ) )
          reminders2.erase( reminders2.begin() + j );
      }
      i++;
    }
#if 1
    if ( reminders2.size() == 1 )
    {
      tts[0] = _xl[sorted[0]];
      tts[1] = _xl[sorted[1]];
      children = { indices[sorted[0]] + offset, indices[sorted[1]] + offset, null_lit, null_lit };
      children[2] = delay_oriented_shannon3( reminders2[0], times, indices );
      if ( children[2] == null_lit )
        return null_lit;
      tts[2] = reminders2[0];
      auto [fun, sim, ex] = chatterjee_synthesis( tts, 3, tt );
      if ( !ex )
      {
        kitty::print_binary( sim ); std::cout << std::endl;
        kitty::print_binary( tt ); std::cout << std::endl;
      }
      auto lit = create_node( children, fun );
      return lit; 
    }
#endif
#if 1
    if ( reminders3.size() == 1 )
    {
      tts[0] = _xl[sorted[0]];
      tts[1] = _xl[sorted[1]];
      tts[2] = _xl[sorted[2]];
      children = { indices[sorted[0]] + offset, indices[sorted[1]] + offset, indices[sorted[2]] + offset, null_lit };
      children[3] = delay_oriented_shannon3( reminders3[0], times, indices );
      if ( children[3] == null_lit )
        return null_lit;
      tts[3] = reminders3[0];
      auto [fun, sim, ex] = chatterjee_synthesis( tts, 4, tt );
      if ( !ex )
      {
        kitty::print_binary( sim ); std::cout << std::endl;
        kitty::print_binary( tt ); std::cout << std::endl;
      }
      auto lit = create_node( children, fun );
      return lit; 
    }
#endif
#if 0
    {
      spfd_manager.init( tt );
      spfd_manager.update( _xl[sorted[0]] );
      spfd_manager.update( _xl[sorted[1]] );
      spfd_manager.update( _xl[sorted[2]] );
      uint32_t trans = std::numeric_limits<uint32_t>::max();
      double best_time = std::numeric_limits<double>::max();
      kitty::static_truth_table<num_vars_l> ttnew;
      for ( int m = 1; m < ( 1u << spfd_manager.num_masks )/2; ++m  )
      {
        ttnew = tt;
        for ( int k = 0; k < spfd_manager.num_masks; ++k )
        {
          if ( ( ( m >> k ) & 0x1 ) == 0x1 )
            ttnew ^= spfd_manager.masks[k];
        }
        auto supp = kitty::min_base_inplace( ttnew );
        if ( supp.size() < indices.size() )
        {
          double worst_time = 0;
          for ( auto s : supp )
          {
            worst_time = std::max( worst_time, times[s] );
          }
          if ( worst_time < best_time )
          {
            best_time = worst_time;
            trans = m;
          }
        }
      }
      if ( ( trans != std::numeric_limits<uint32_t>::max() ) && ( best_time < times[sorted[0]] ) )
      {
        ttnew = tt;
        for ( int k = 0; k < spfd_manager.num_masks; ++k )
        {
          if ( ( ( trans >> k ) & 0x1 ) == 0x1 )
            ttnew ^= spfd_manager.masks[k];
        }
        children[0] = indices[sorted[0]] + offset;
        children[1] = indices[sorted[1]] + offset;
        children[2] = indices[sorted[2]] + offset;
        children[3] = delay_oriented_shannon3( ttnew, times, indices );
        
        tts[0] = _xl[sorted[0]];
        tts[1] = _xl[sorted[1]];
        tts[2] = _xl[sorted[2]];
        tts[3] = ttnew;

        auto [fun, sim, ex] = chatterjee_synthesis( tts, 4, tt );
        auto lit = create_node( children, fun );
        return lit;
      }
    }
#endif
#if 1
    if ( reminders2.size() == 2 )
    {
      tts[0] = _xl[sorted[0]];
      tts[1] = _xl[sorted[1]];
      children = { indices[sorted[0]] + offset, indices[sorted[1]] + offset, null_lit, null_lit };
      for ( int i = 0; i < reminders2.size(); ++i )
      {
        children[2+i] = delay_oriented_shannon3( reminders2[i], times, indices );
        if ( children[2+i] == null_lit )
          return null_lit;
        tts[2+i] = reminders2[i];
      }
      auto [fun, sim, ex] = chatterjee_synthesis( tts, 4, tt );
      if ( !ex )
      {
        kitty::print_binary( sim ); std::cout << std::endl;
        kitty::print_binary( tt ); std::cout << std::endl;
      }
      auto lit = create_node( children, fun );
      return lit; 
    }
#endif

    tts[0] = _xl[sorted[0]];
    tts[1] = tt0;
    tts[2] = tt1;
    children = { indices[sorted[0]] + offset, null_lit, null_lit, null_lit };
    children[1] = delay_oriented_shannon3( tt0, times, indices );
    children[2] = delay_oriented_shannon3( tt1, times, indices );
    auto [fun, sim, ex] = chatterjee_synthesis( tts, 3, tt );
    auto lit = create_node( children, fun );
    return lit; 


#if 0
    lit_t lit_idx = indices[worst_index] + offset;
    auto tt0 = kitty::cofactor0( tt, worst_index );
    auto tt1 = kitty::cofactor1( tt, worst_index );
    if ( kitty::equal( tt0, ~tt1 ) )
    {
      lit_t l0 = delay_oriented_shannon( tt0, times, indices );

      kitty::static_truth_table<num_vars> x0;
      kitty::static_truth_table<num_vars> x1;
      kitty::static_truth_table<num_vars> fn;
      kitty::create_nth_var( x0, 0 );
      kitty::create_nth_var( x1, 1 );
      fn = x0 ^ x1;
      lit_t lit = nodes.size() + offset;
      std::array<lit_t, num_vars> children{ lit_idx, l0, null_lit, null_lit };
      nodes.emplace_back( children, fn, lit );
      return lit;
    }
    lit_t l0 = delay_oriented_shannon( tt0, times, indices );
    lit_t l1 = delay_oriented_shannon( tt1, times, indices );

    if ( l0 == 0 ) // F = x F1
    {
      kitty::static_truth_table<num_vars> x0;
      kitty::static_truth_table<num_vars> x1;
      kitty::static_truth_table<num_vars> fn;
      kitty::create_nth_var( x0, 0 );
      kitty::create_nth_var( x1, 1 );
      fn = x0 & x1;
      lit_t lit = nodes.size() + offset;
      std::array<lit_t, num_vars> children{ lit_idx, l1, null_lit, null_lit };
      nodes.emplace_back( children, fn, lit );
      return lit;
    }
    else if ( l0 == 1 ) // F = x' + F1
    {
      kitty::static_truth_table<num_vars> x0;
      kitty::static_truth_table<num_vars> x1;
      kitty::static_truth_table<num_vars> fn;
      kitty::create_nth_var( x0, 0 );
      kitty::create_nth_var( x1, 1 );
      fn = ( ~x0 ) | x1;
      lit_t lit = nodes.size() + offset;
      std::array<lit_t, num_vars> children{ lit_idx, l1, null_lit, null_lit };
      nodes.emplace_back( children, fn, lit );
      return lit;
    }
    else if ( l1 == 0 ) // F = x'F0
    {
      kitty::static_truth_table<num_vars> x0;
      kitty::static_truth_table<num_vars> x1;
      kitty::static_truth_table<num_vars> fn;
      kitty::create_nth_var( x0, 0 );
      kitty::create_nth_var( x1, 1 );
      fn = ( ~x0 ) & x1;
      lit_t lit = nodes.size() + offset;
      std::array<lit_t, num_vars> children{ lit_idx, l0, null_lit, null_lit };
      nodes.emplace_back( children, fn, lit );
      return lit;
    }
    else if ( l1 == 1 ) // F = x + F0
    {
      kitty::static_truth_table<num_vars> x0;
      kitty::static_truth_table<num_vars> x1;
      kitty::static_truth_table<num_vars> fn;
      kitty::create_nth_var( x0, 0 );
      kitty::create_nth_var( x1, 1 );
      fn = x0 | x1;
      lit_t lit = nodes.size() + offset;
      std::array<lit_t, num_vars> children{ lit_idx, l0, null_lit, null_lit };
      nodes.emplace_back( children, fn, lit );
      return lit;
    }
    else
    {
      kitty::static_truth_table<num_vars> x0;
      kitty::static_truth_table<num_vars> x1;
      kitty::static_truth_table<num_vars> x2;
      kitty::static_truth_table<num_vars> fn;
      kitty::create_nth_var( x0, 0 );
      kitty::create_nth_var( x1, 1 );
      kitty::create_nth_var( x2, 2 );
      fn = ( x0 & x1 ) | ( ( ~x0 ) & x2 );
      lit_t lit = nodes.size() + offset;
      std::array<lit_t, num_vars> children{ lit_idx, l1, l0, null_lit };
      nodes.emplace_back( children, fn, lit );
      return lit;
    }
#endif
  }

  std::vector<uint32_t> reduce_support( std::vector<uint32_t> & vars, kitty::static_truth_table<num_vars_l> func, kitty::static_truth_table<num_vars_l> care )
  {
    std::vector<uint32_t> supp;
    for ( auto v : vars )
    {
      auto func0 = kitty::cofactor0( func, v );
      auto func1 = kitty::cofactor1( func, v );
      auto care0 = kitty::cofactor0( care, v );
      auto care1 = kitty::cofactor1( care, v );
      if ( kitty::equal( func0&care0&care1, func1&care0&care1 ) )
      {
        func = (func0&care0) | (func1&care1);
        care = care0 | care1;
      }
      else
        supp.push_back( v );
    }

    return supp;
  }

  std::tuple<kitty::static_truth_table<num_vars>, kitty::static_truth_table<num_vars>> get_function( std::vector<uint32_t> const& vars, kitty::static_truth_table<num_vars_l> const& func, kitty::static_truth_table<num_vars_l> const& care )
  {
    kitty::static_truth_table<num_vars> tt;
    kitty::static_truth_table<num_vars> mk;

    uint32_t num_minterms = 1u << vars.size();

    kitty::static_truth_table<num_vars_l> tmp_l;
    kitty::static_truth_table<num_vars> tmp_s;
    for ( uint32_t m = 0; m < num_minterms; ++m )
    {
      tmp_l |= ~tmp_l;
      tmp_s |= ~tmp_s;
      for ( int v = 0; v < vars.size(); ++v )
      {
        tmp_l &= ( ( ( m >> v ) & 0x1 ) == 0x1 ) ? _xl[vars[v]] : ~_xl[vars[v]];
        tmp_s &= ( ( ( m >> v ) & 0x1 ) == 0x1 ) ? _xs[v] : ~_xs[v];
      }
      if ( kitty::count_ones( func & tmp_l ) > 0 )
        tt |= tmp_s;
      else
        tt &= ~tmp_s;
      if ( kitty::count_ones( care & tmp_l ) > 0 )
        mk |= tmp_s;
      else
        mk &= ~tmp_s;
    }
    return std::make_tuple( tt, mk );
  }

  std::tuple<kitty::static_truth_table<num_vars>, kitty::static_truth_table<num_vars>, kitty::static_truth_table<num_vars_l>> get_function_with_sim( std::vector<uint32_t> const& vars, kitty::static_truth_table<num_vars_l> const& func, kitty::static_truth_table<num_vars_l> const& care )
  {
    //for ( auto v : vars )
    //  std::cout << v << " ";
    //std::cout << " <= vars" << std::endl;
    kitty::static_truth_table<num_vars> tt;
    kitty::static_truth_table<num_vars> mk;
    kitty::static_truth_table<num_vars_l> sim;

    uint32_t num_minterms = 1u << vars.size();

    kitty::static_truth_table<num_vars_l> tmp_l;
    kitty::static_truth_table<num_vars> tmp_s;
    for ( uint32_t m = 0; m < num_minterms; ++m )
    {
      //std::cout << "M=" << m << std::endl;
      tmp_l |= ~tmp_l;
      tmp_s |= ~tmp_s;
      for ( uint32_t v = 0; v < vars.size(); ++v )
      {
        tmp_l &= ( ( ( m >> v ) & 0x1 ) == 0x1 ) ? _xl[vars[v]] : ~_xl[vars[v]];
        //if ( kitty::is_const0( _xs[v] )||kitty::is_const0( ~_xs[v] ) )
        {
          kitty::create_nth_var( _xs[v], v );
          //kitty::print_binary( _xs[v] ); std::cout << " xs[" << v << "]" << std::endl;
        }
     //   kitty::print_binary( _xs[v] ); std::cout << " xs[" << v << "]" << std::endl;
        tmp_s &= ( ( ( m >> v ) & 0x1 ) == 0x1 ) ? _xs[v] : ~_xs[v];
      }
      //kitty::print_binary( tmp_s ); std::cout << " tm" <<std::endl;

      uint32_t num_ones = kitty::count_ones( care & func & tmp_l );
      uint32_t num_zeros = kitty::count_ones( care & ~func & tmp_l );
      if ( (num_ones > 0)&&(num_zeros > 0) )
      std::cerr << "[e] illegal termination condition" << std::endl;
      if ( num_ones > 0 )
      {
        tt |= tmp_s;
        sim |= tmp_l;
      }
      else
      {
        tt &= ~tmp_s;
        sim &= ~tmp_l;
      }
      //std::cout << "#1=" << num_ones << " #0=" << num_zeros << std::endl;
      //kitty::print_binary( tmp_s ); std::cout << " tm" <<std::endl;
      //kitty::print_binary( tt ); std::cout << " tt" <<std::endl;
      if ( kitty::count_ones( care & tmp_l ) > 0 )
        mk |= tmp_s;
      else
        mk &= ~tmp_s;
    }
    return std::make_tuple( tt, mk, sim );
  }

  std::array<uint32_t, num_vars> get_leaves( std::vector<uint32_t> const& vars, std::vector<uint32_t> const& indices )
  {
    std::array<uint32_t, num_vars> leaves;
    for ( int i = 0; i < num_vars; ++i )
      leaves[i] = i < vars.size() ? indices[vars[i]] + offset : null_lit;

    return leaves;
  }

  lit_t delay_oriented_shannon_dc( kitty::static_truth_table<num_vars_l> func, kitty::static_truth_table<num_vars_l> care, std::vector<double> times, std::vector<uint32_t> indices )
  {
bool VERBOSE=false;
if ( VERBOSE )
{
std::cout << "here 0" << std::endl;
for ( int i = 0; i < func.num_bits(); ++i )
{
  int c = kitty::get_bit( care, i );
  int b = kitty::get_bit( func, i );
  if ( c == 0 )
    std::cout << "*";
  else
    std::cout << b;
}
std::cout << std::endl;
}
    if ( kitty::is_const0( func & care ) )
      return 0;
    if ( kitty::equal( func & care, care ) )
      return 1;
    std::vector<uint32_t> vars;
    for ( uint32_t i = 0; i < indices.size(); ++i )
      vars.push_back( i );
    
    auto supp = reduce_support( vars, func, care );
if ( VERBOSE )
{
for ( auto s : supp )
  std::cout << s << " ";
std::cout << std::endl;
}
    if ( supp.size() <= num_vars )
    {
if ( VERBOSE )
{
std::cout << "done" << std::endl;
}
      auto [ tt, mk, sim ] = get_function_with_sim( supp, func, care );

if (VERBOSE)
{
for ( int i = 0; i < mk.num_bits(); ++i )
{
  int c = kitty::get_bit( mk, i );
  int b = kitty::get_bit( tt, i );
  if ( c == 0 )
    std::cout << "*";
  else
    std::cout << b;
}
std::cout << std::endl;
}
      auto children = get_leaves( supp, indices );
if ( VERBOSE )
{
for ( auto c : children )
  std::cout << c << " ";
std::cout << std::endl;
}
      auto lit = create_node( children, tt );
      if( lit <= 1 )
        return lit;
      _xl[lit-offset] = sim;
if ( VERBOSE )
std::cout << lit << std::endl;
      return lit;
    }

    /* try top decomposition */
    std::sort( supp.begin(), supp.end(), [&]( auto va, auto vb ){
      return times[va] > times[vb];
    } );

    auto const& func0 = kitty::cofactor0( func, supp[0] );
    auto const& func1 = kitty::cofactor1( func, supp[0] );
    auto care0 = kitty::cofactor0( care, supp[0] );
    auto care1 = kitty::cofactor1( care, supp[0] );
    std::array<uint32_t, num_vars> children;
    std::array<kitty::static_truth_table<num_vars_l>, num_vars> tts;
    if ( ( kitty::count_ones( func0 & care ) == 0 ) || ( kitty::count_ones( (~func0)&care ) == 0 ) )
    {
if ( VERBOSE )
  std::cout << "topp 1" << std::endl;
      auto new_care = ( kitty::count_ones( func0 & care ) == 0 ) ? care & _xl[supp[0]] : care & ~_xl[supp[0]];
      auto lit_new = delay_oriented_shannon_dc( func1, new_care, times, indices );
      if ( lit_new == null_lit )
        return null_lit;
      if( lit_new <= 1 )
        return lit_new;
      children = { indices[supp[0]] + offset, lit_new, null_lit, null_lit };
      tts[0] = _xl[supp[0]];
      tts[1] = _xl[lit_new-offset];

      auto [fun, sim] = chatterjee_synthesis_with_sim( tts, 2, func, care );
      auto lit = create_node( children, fun );
      if( lit <= 1 )
        return lit;
      _xl[lit-offset] = sim;

      return lit;
    }
    if ( ( kitty::count_ones( func1 & care ) == 0 ) || ( kitty::count_ones( (~func1)&care ) == 0 ) )
    {
if ( VERBOSE )
  std::cout << "top 0" << std::endl;
      auto new_care = ( kitty::count_ones( func1 & care ) == 0 ) ? care & ~_xl[supp[0]] : care & _xl[supp[0]];

      auto lit_new = delay_oriented_shannon_dc( func0, new_care, times, indices );
      if ( lit_new == null_lit )
        return null_lit;
      if( lit_new <= 1 )
        return lit_new;
      children = { indices[supp[0]] + offset, lit_new, null_lit, null_lit };
      tts[0] = _xl[supp[0]];
      tts[1] = _xl[lit_new-offset];
      auto [fun, sim] = chatterjee_synthesis_with_sim( tts, 2, func, care );
      auto lit = create_node( children, fun );
      if( lit_new <= 1 )
        return lit;
      _xl[lit-offset] = sim;

      return lit;
    }
    if ( ( kitty::equal( func0 & care, (~func1) & care ) ) )
    {
if ( VERBOSE )
  std::cout << "top ^" << std::endl;

      auto new_func = func0;
      auto lit_new = delay_oriented_shannon_dc( new_func, care, times, indices );
      if ( lit_new == null_lit )
        return null_lit;
      if( lit_new <= 1 )
        return lit_new;
      children = { indices[supp[0]] + offset, lit_new, null_lit, null_lit };
      tts[0] = _xl[supp[0]];
      tts[1] = new_func;
      auto [fun, sim] = chatterjee_synthesis_with_sim( tts, 2, func, care );
      auto lit = create_node( children, fun );
      if( lit <= 1 )
        return lit;
      _xl[lit-offset] = sim;

      return lit;
    }
    int szA, szB;

    {
      spfd_manager.init( func, care );
      spfd_manager.update( _xl[supp[0]] );
      spfd_manager.update( _xl[supp[1]] );
if ( VERBOSE )
  std::cout << "#masks " << spfd_manager.num_masks << std::endl;
      children = { indices[supp[0]] + offset, indices[supp[1]] + offset, null_lit, null_lit };
      auto const& mks = spfd_manager.masks;
      std::array<std::array<kitty::static_truth_table<num_vars_l>, 4u>, 28u> cands;
      cands[0] =  { mks[0], mks[1]|mks[2]|mks[3], func & mks[0], func & ( mks[3]|mks[1]|mks[2]) };
      cands[1] =  { mks[0], mks[1]|mks[2]|mks[3], func & mks[0], ( mks[1] & ( func ^ mks[1] ) ) | ( func & ( mks[2] | mks[3] ) ) };
      cands[2] =  { mks[0], mks[1]|mks[2]|mks[3], func & mks[0], ( mks[2] & ( func ^ mks[2] ) ) | ( func & ( mks[1] | mks[3] ) ) };
      cands[3] =  { mks[0], mks[1]|mks[2]|mks[3], func & mks[0], ( mks[3] & ( func ^ mks[3] ) ) | ( func & ( mks[1] | mks[2] ) ) };
      cands[4] =  { mks[1], mks[0]|mks[2]|mks[3], func & mks[1], func & ( mks[0]|mks[3]|mks[2]) };
      cands[5] =  { mks[1], mks[0]|mks[2]|mks[3], func & mks[1], ( mks[0] & ( func ^ mks[0] ) ) | ( func & ( mks[2] | mks[3] ) ) };
      cands[6] =  { mks[1], mks[0]|mks[2]|mks[3], func & mks[1], ( mks[2] & ( func ^ mks[2] ) ) | ( func & ( mks[0] | mks[3] ) ) };
      cands[7] =  { mks[1], mks[0]|mks[2]|mks[3], func & mks[1], ( mks[3] & ( func ^ mks[3] ) ) | ( func & ( mks[0] | mks[2] ) ) };
      cands[8] =  { mks[2], mks[0]|mks[1]|mks[3], func & mks[2], func & ( mks[0]|mks[1]|mks[3]) };
      cands[9] =  { mks[2], mks[0]|mks[1]|mks[3], func & mks[2], ( mks[0] & ( func ^ mks[0] ) ) | ( func & ( mks[1] | mks[3] ) ) };
      cands[10] = { mks[2], mks[0]|mks[1]|mks[3], func & mks[2], ( mks[1] & ( func ^ mks[1] ) ) | ( func & ( mks[0] | mks[3] ) ) };
      cands[11] = { mks[2], mks[0]|mks[1]|mks[3], func & mks[2], ( mks[3] & ( func ^ mks[3] ) ) | ( func & ( mks[0] | mks[1] ) ) };
      cands[12] = { mks[3], mks[1]|mks[2]|mks[0], func & mks[3], func & ( mks[0]|mks[1]|mks[2]) };
      cands[13] = { mks[3], mks[1]|mks[2]|mks[0], func & mks[3], ( mks[0] & ( func ^ mks[0] ) ) | ( func & ( mks[1] | mks[2] ) ) };
      cands[14] = { mks[3], mks[1]|mks[2]|mks[0], func & mks[3], ( mks[1] & ( func ^ mks[1] ) ) | ( func & ( mks[0] | mks[2] ) ) };
      cands[15] = { mks[3], mks[1]|mks[2]|mks[0], func & mks[3], ( mks[2] & ( func ^ mks[2] ) ) | ( func & ( mks[0] | mks[1] ) ) };
      cands[16] = { mks[0] | mks[1], mks[2] | mks[3], func & ( mks[0] | mks[1] ), func & ( mks[2] | mks[3] ) };
      cands[17] = { mks[0] | mks[1], mks[2] | mks[3], ( func & mks[1] ) | ( ( func ^ mks[0] ) & mks[0] ) , func & ( mks[2] | mks[3] ) };
      cands[18] = { mks[0] | mks[1], mks[2] | mks[3], func & ( mks[0] | mks[1] ), ( func & mks[3] ) | ( ( func ^ mks[2] ) & mks[2] ) };
      cands[19] = { mks[0] | mks[1], mks[2] | mks[3], ( func & mks[1] ) | ( ( func ^ mks[0] ) & mks[0] ), ( func & mks[3] ) | ( ( func ^ mks[2] ) & mks[2] ) };
      cands[20] = { mks[0] | mks[3], mks[2] | mks[1], func & ( mks[0] | mks[3] ), func & ( mks[1] | mks[2] ) };
      cands[21] = { mks[0] | mks[3], mks[2] | mks[1], ( func & mks[3] ) | ( ( func ^ mks[0] ) & mks[0] ) , func & ( mks[1] | mks[2] ) };
      cands[22] = { mks[0] | mks[3], mks[2] | mks[1], func & ( mks[0] | mks[3] ), ( func & mks[2] ) | ( ( func ^ mks[1] ) & mks[1] ) };
      cands[23] = { mks[0] | mks[3], mks[2] | mks[1], ( func & mks[3] ) | ( ( func ^ mks[0] ) & mks[0] ), ( func & mks[2] ) | ( ( func ^ mks[1] ) & mks[1] ) };
      cands[24] = { mks[0] | mks[2], mks[3] | mks[1], func & ( mks[0] | mks[2] ), func & ( mks[1] | mks[3] ) };
      cands[25] = { mks[0] | mks[2], mks[3] | mks[1], ( func & mks[2] ) | ( ( func ^ mks[0] ) & mks[0] ) , func & ( mks[1] | mks[3] ) };
      cands[26] = { mks[0] | mks[2], mks[3] | mks[1], func & ( mks[0] | mks[2] ), ( func & mks[3] ) | ( ( func ^ mks[1] ) & mks[1] ) };
      cands[27] = { mks[0] | mks[2], mks[3] | mks[1], ( func & mks[2] ) | ( ( func ^ mks[0] ) & mks[0] ), ( func & mks[3] ) | ( ( func ^ mks[1] ) & mks[1] ) };
      int best_i = -1;
      double best_t = times[supp[1]];
      double best_p = ( ( vars.size()-1) * ( vars.size()-1) );
      for ( int i = 0; i < 28; ++i )
      {
        auto const ttA = cands[i][2];
        auto const ttB = cands[i][3];
        auto const mkA = cands[i][0];
        auto const mkB = cands[i][1];
      //  if ( kitty::is_const0( ttA & mkA ) || kitty::is_const0( ttB & mkB )|| kitty::is_const0( ~ttA & mkA ) || kitty::is_const0( ~ttB & mkB ) )
      //    continue;
if( false && VERBOSE )
{
 kitty::print_binary( ttA ); std::cout << std::endl;
 kitty::print_binary( ttB ); std::cout << std::endl;
 kitty::print_binary( mkA ); std::cout << std::endl;
 kitty::print_binary( mkB ); std::cout << std::endl;
}
        auto suppA = reduce_support( vars, ttA, mkA );
        auto suppB = reduce_support( vars, ttB, mkB );
        if ( (suppA.size() >= supp.size())||(suppB.size() >= supp.size()) )
          continue;
if( VERBOSE )
{
std::cout << "I" << i << std::endl;
 for ( auto s : suppA )
   std::cout << s << "." << times[s] << " ";
 std::cout << std::endl;
 for ( auto s : suppB )
   std::cout << s << "." << times[s] << " ";
 std::cout << std::endl;
}
        std::sort( suppA.begin(), suppA.end(), [&]( auto va, auto vb ){
          return times[va] > times[vb];
        } );
        std::sort( suppB.begin(), suppB.end(), [&]( auto va, auto vb ){
          return times[va] > times[vb];
        } );

        if ( (suppA.size() > 0) && ( times[suppA[0]] > times[supp[0]] ))
          continue;
        if ( (suppB.size() > 0) && (times[suppB[0]] > times[supp[0]]) )
          continue;
        double tA = suppA.size() > 0 ? times[suppA[0]] : times[supp[1]];
        double tB = suppB.size() > 0 ? times[suppB[0]] : times[supp[1]];

        double worst = std::max( tA, tB );
        if ( worst < best_t )
        {
          best_t = worst;
          best_i = i;
          best_p = std::max( 1, (int)suppA.size())*std::max( 1, (int)suppB.size());
          szA = suppA.size();
          szB = suppB.size();
        }
        else if ( ( worst == best_t ) && 
                  ( ( suppA.size() == 0 ) || ( suppA[0] != supp[0] ) ) && 
                  ( ( suppB.size() == 0 ) || ( suppB[0] != supp[0] ) ) &&
                  ( std::max( 1, (int)suppA.size())*std::max( 1, (int)suppB.size()) < best_p ) )
        {
          best_t = worst;
          best_i = i;
          best_p = std::max( 1, (int)suppA.size())*std::max( 1, (int)suppB.size());
          szA = suppA.size();
          szB = suppB.size();
        }
      }
      if ( best_i >= 0 )
      {
if ( VERBOSE )
{
 std::cout << "*" << best_i << " ";
 std::cout << best_t << " ";  
 std::cout << best_p << std::endl;
}
        auto const ttA = cands[best_i][2];
        auto const ttB = cands[best_i][3];
        auto const mkA = cands[best_i][0];
        auto const mkB = cands[best_i][1];
if( VERBOSE )
{
for ( int b = 0; b < ttA.num_bits(); ++b )
{
  int bit = kitty::get_bit( ttA, b );
  if ( kitty::get_bit( mkA, b ) == 1 )
    std::cout << bit;
  else
    std::cout << "*";
}
std::cout << std::endl;
for ( int b = 0; b < ttB.num_bits(); ++b )
{
  int bit = kitty::get_bit( ttB, b );
  if ( kitty::get_bit( mkB, b ) == 1 )
    std::cout << bit;
  else
    std::cout << "*";
}
std::cout << std::endl;
}
        children = { indices[supp[0]] + offset, indices[supp[1]] + offset, null_lit, null_lit };

        tts[0] = _xl[supp[0]];
        tts[1] = _xl[supp[1]];
        int num_v = 3;
        if ( szA == 0 )
        {
          children[2] = delay_oriented_shannon_dc( ttB, mkB, times, indices );
          if( children[2] > 1 )
            tts[2] = _xl[children[2]-offset];
          else
            tts[2] ^= children[2] == 0 ? tts[2] : ~tts[2];

        }
        else if ( szB == 0 )
        {
          children[2] = delay_oriented_shannon_dc( ttA, mkA, times, indices );
          if( children[2] > 1 )
            tts[2] = _xl[children[2]-offset];
          else
            tts[2] ^= children[2] == 0 ? tts[2] : ~tts[2];
        }
        else
        {
          children[2] = delay_oriented_shannon_dc( ttA, mkA, times, indices );
          children[3] = delay_oriented_shannon_dc( ttB, mkB, times, indices );
          if ( children[3] <= 1 )
          {
            num_v = 3;
            if( children[2] > 1 )
              tts[2] = _xl[children[2]-offset];
            else
              tts[2] ^= children[2] == 0 ? tts[2] : ~tts[2];
          }
          else if ( children[2] <= 1 )
          {
            num_v = 3;
            children[2]=children[3];
            if( children[2] > 1 )
              tts[2] = _xl[children[2]-offset];
            else
              tts[2] ^= children[2] == 0 ? tts[2] : ~tts[2];
          }
          else
          {
            if( children[2] > 1 )
              tts[2] = _xl[children[2]-offset];
            else
              tts[2] ^= children[2] == 0 ? tts[2] : ~tts[2];
            if( children[3] > 1 )
              tts[3] = _xl[children[3]-offset];
            else
              tts[3] ^= children[3] == 0 ? tts[3] : ~tts[3];
            num_v = 4;
          }
        }
if( VERBOSE )
 std::cout << "ch" << std::endl;
        auto [fun, sim] = chatterjee_synthesis_with_sim( tts, num_v, func, care );
if( VERBOSE )
  std::cout << "cr" << std::endl;

        auto lit = create_node( children, fun );
        if ( lit <= 1 )
          return lit;
if ( VERBOSE )
  std::cout << "as" << lit << " " << lit-offset << std::endl;
        _xl[lit-offset] = sim;
if ( VERBOSE )
  std::cout << "re" << std::endl;

        return lit; 
      }

    }
// std::cout << "shannon" << std::endl;

    tts[0] = _xl[supp[0]];
    children = { indices[supp[0]] + offset, null_lit, null_lit, null_lit };
    children[1] = delay_oriented_shannon_dc( func0, care & ~_xl[supp[0]], times, indices );
    children[2] = delay_oriented_shannon_dc( func1, care & _xl[supp[0]], times, indices );
    if( children[1] > 1 )
      tts[1] = _xl[children[1]-offset];
    else
      tts[1] ^= children[1] == 0 ? tts[1] : ~tts[1];
    if( children[2] > 1 )
      tts[2] = _xl[children[2]-offset];
    else
      tts[2] ^= children[2] == 0 ? tts[2] : ~tts[2];
    auto [fun, sim] = chatterjee_synthesis_with_sim( tts, 3, func, care );
    auto lit = create_node( children, fun );
    _xl[lit-offset] = sim;

    return lit; 
  }



  lit_t dc_decompose( kitty::static_truth_table<num_vars_l> func, std::vector<double> times, std::vector<uint32_t> indices )
  {
    /* check for constant */
    if ( kitty::is_const0( func ) )
      return 0;
    if ( kitty::is_const0( ~func ) )
      return 1;
    /* step 1 : reduce the support */
    auto supp = kitty::min_base_inplace( func );
    if ( ( supp.size() > indices.size() ) || ( supp.back() >= indices.size() ) )
      return null_lit;

    int j = supp.size()-1;
    for ( int i = indices.size() - 1; i >= 0; --i )
    {
      if ( (j >= 0) && (i == static_cast<int>( supp[j])) )
      {
        j -= 1;
      }
      else
      {
        indices.erase( indices.begin() + i );
        times.erase( times.begin() + i );
      }
    }
//std::cout << "1" << std::endl;
    /* here the truth table is reduced and the support minimized */
    if ( indices.size() <= num_vars )
    {
      kitty::static_truth_table<num_vars> tt;
      kitty::shrink_to_inplace( tt, func );
      std::array<lit_t, num_vars> children;
      for ( int i = 0; i < num_vars; ++i )
      {
        if ( i < indices.size() )
        {
          children[i] = indices[i] + offset;
        }
        else
        {
          children[i] = null_lit;
        }
      }
      lit_t lit = nodes.size() + offset;
      nodes.emplace_back( children, tt, lit );
      return lit;
    }
    std::vector<std::pair<int, double>> var_time;
    for ( int i = 0; i < indices.size(); ++i )
      var_time.push_back( std::make_pair( i, times[i] ) );
    std::sort( var_time.begin(), var_time.end(), [&]( auto a, auto b ){
      return a.second > b.second;
    } );
//std::cout << "2" << std::endl;

    if ( indices.size() <= ( 2*num_vars - 1 ) )
    {
      spfd_manager.init( func );
      std::array<kitty::static_truth_table<num_vars_l>, num_vars> tts;
      std::array<uint32_t, num_vars> children_free;
//std::cout << "3" << std::endl;
      for ( int i = 0; i < num_vars - 1; ++i )
      {
        spfd_manager.update( _xl[var_time[i].first] );
        tts[i] = _xl[var_time[i].first];
        children_free[i] = indices[var_time[i].first] + offset;
      }
      std::vector<uint32_t> tmp_indices;
      std::vector<double> tmp_times;
//std::cout << "4" << std::endl;

      for ( int m = 0; m < ( 1u << spfd_manager.num_masks )/2; ++m )
      {
        tmp_indices = indices;
        tmp_times = times;
        kitty::static_truth_table<num_vars_l> func_bound;
//std::cout << "5" << std::endl;
        for ( int r = 0; r < spfd_manager.num_masks; ++r )
        {
          if ( ( ( m >> r ) & 0x1 ) == 0x1 )
            func_bound |= spfd_manager.masks[r] & func;
          else
            func_bound |= spfd_manager.masks[r] & ~func;
        }

        auto supp = kitty::min_base_inplace( func_bound );

        int j = supp.size()-1;
        bool not_good = false;
//std::cout << "6" << std::endl;

        for ( int i = tmp_indices.size() - 1; i >= 0; --i )
        {
          //std::cout << supp.size() << "=supsize | " << j << std::endl;
          //std::cout << tmp_times.size() << "=tmptimes | " << i << std::endl;
          if ( (j >= 0) && (i == static_cast<int>( supp[j])) )
          {
            j -= 1;
            if ( supp[j] == var_time[0].first )
              not_good = true;
          }
          else
          {
            tmp_indices.erase( tmp_indices.begin() + i );
            tmp_times.erase( tmp_times.begin() + i );
          }
        }
        if ( not_good )
        {
          //std::cout << "not good" << std::endl;
        }
        else
        {

//        for ( int v = tmp_indices.size() - 1; v >= 0; --v )
//        {
//          auto func0 = kitty::cofactor0( func_bound, v );
//          auto func1 = kitty::cofactor1( func_bound, v );
//          if ( kitty::equal( func0, func1 ) )
//          {
//            tmp_indices.erase( tmp_indices.begin() + v );
//            tmp_times.erase( tmp_times.begin() + v );
//          }
//        }
          if ( tmp_indices.size() <= num_vars )
          {
            auto supp = kitty::min_base_inplace( func_bound );
            kitty::static_truth_table<num_vars> tt_bound;
            kitty::shrink_to_inplace( tt_bound, func_bound );
            std::array<lit_t, num_vars> children_bound;
//std::cout << "7" << std::endl;

            for ( int i = 0; i < num_vars; ++i )
            {
              if ( i < tmp_indices.size() )
              {
                children_bound[i] = tmp_indices[i] + offset;
              }
              else
              {
                children_bound[i] = null_lit;
              }
            }
            lit_t lit_bound = nodes.size() + offset;
            nodes.emplace_back( children_bound, tt_bound, lit_bound );
            tts[num_vars-1] = func_bound;
            children_free[num_vars-1] = lit_bound;
            kitty::static_truth_table<num_vars> tt_free;
            kitty::static_truth_table<num_vars_l> tmp;
//std::cout << "8" << std::endl;

            for ( int min = 0; min < ( 1u << num_vars ); ++min )
            {
              tmp |= ~tmp;
//std::cout << "9" << std::endl;

              for ( int v = 0; v < num_vars; ++v )
              {
                if ( ( ( min >> v ) & 0x1 ) == 0x1 )
                  tmp &= tts[v];
                else
                  tmp &= ~tts[v];
              }
              if ( kitty::count_ones( tmp & func ) > 0 )
                kitty::set_bit( tt_free, min );
              else
                kitty::clear_bit( tt_free, min );
            }

            lit_t lit_free = nodes.size() + offset;
            nodes.emplace_back( children_free, tt_free, lit_free );

            return lit_free;
          }
        }
      }

    }


    if ( indices.size() <= ( 3*num_vars - 2 ) )
    {
      spfd_manager.init( func );
      std::array<kitty::static_truth_table<num_vars_l>, num_vars> tts;
      std::array<uint32_t, num_vars> children_free;
      for ( int i = 0; i < num_vars - 2; ++i )
      {
        spfd_manager.update( _xl[var_time[i].first] );
        tts[i] = _xl[var_time[i].first];
        children_free[i] = indices[var_time[i].first] + offset;
      }
      std::vector<uint32_t> tmp_indicesA;
      std::vector<double> tmp_timesA;
      std::vector<uint32_t> tmp_indicesB;
      std::vector<double> tmp_timesB;
      for ( int t = 1; t < 8; ++t )
      {
        std::vector<uint32_t> A;
        std::vector<uint32_t> B;
        for ( int k = 0; k < 4; ++k )
        {
          if ( ( ( t >> k ) & 0x1 ) == 0x1 )
            A.push_back( k );
          else
            B.push_back( k );
        }
        kitty::static_truth_table<num_vars_l> funcA;
        kitty::static_truth_table<num_vars_l> careA;
        kitty::static_truth_table<num_vars_l> funcB;
        kitty::static_truth_table<num_vars_l> careB;
        auto const& masks = spfd_manager.masks;
        for ( auto m : A )
          careA |= masks[m];
        for ( auto m : B )
          careB |= masks[m];

        for ( int mA = 0; mA < ( 1u << A.size() ); ++mA )
        {
          for ( int mB = 0; mB < ( 1u << B.size() ); ++mB )
          {
            funcA ^= funcA;
            funcB ^= funcB;
            for ( int iA = 0; iA < A.size(); ++iA )
            {
              if ( ((mA >> iA)&0x1)==0x1 )
                funcA |= func & masks[A[iA]];
              else
                funcA |= ~func & masks[A[iA]];
            }
            for ( int iB = 0; iB < B.size(); ++iB )
            {
              if ( ((mB >> iB)&0x1)==0x1 )
                funcB |= func & masks[B[iB]];
              else
                funcB |= ~func & masks[B[iB]];
            }
            /* check for support reduction */

            auto indicesA = indices;
            auto indicesB = indices;
            std::vector<int> suppA;
            for ( int i = 0; i < indicesA.size(); ++i )
            {
              auto tt0 = kitty::cofactor0( funcA, i );
              auto tt1 = kitty::cofactor1( funcA, i );
              if ( !kitty::equal( tt0 & careA, tt1 & careA ) )
                suppA.push_back( i );
            }
            std::vector<int> suppB;
            for ( int i = 0; i < indicesB.size(); ++i )
            {
              auto tt0 = kitty::cofactor0( funcB, i );
              auto tt1 = kitty::cofactor1( funcB, i );
              if ( !kitty::equal( tt0 & careB, tt1 & careB ) )
                suppB.push_back( i );
            }
            if ( (suppA.size() <= num_vars)&&(suppB.size() <= num_vars) )
            {
              std::cout << "3DEC" << std::endl;
              /* extract the functionality and the simulation*/
              kitty::static_truth_table<num_vars_l> simA;
              kitty::static_truth_table<num_vars_l> simB;
              kitty::static_truth_table<num_vars> ttA;
              kitty::static_truth_table<num_vars> ttB;
              kitty::static_truth_table<num_vars_l> tmp_l;
              kitty::static_truth_table<num_vars> tmp;
              for ( int m = 0; m < ( 1u << num_vars ); ++m )
              {
                tmp_l ^= tmp_l;
                tmp |= ~tmp;
                for ( int v = 0; v < suppA.size(); ++v )
                {
                  if ( ((m >> v)&0x1)==0x1 )
                    tmp_l &= _xl[suppA[v]];
                  else
                    tmp_l &= ~_xl[suppA[v]];
                }
                if ( kitty::count_ones( tmp_l & careA ) > 0 )
                {
                  if ( kitty::count_ones( tmp_l & funcA ) > 0 )
                  {
                    simA |= tmp_l;
                    kitty::set_bit( ttA, m );
                  }
                  else
                    kitty::clear_bit( ttA, m );
                }
              }
std::cout << "q0" << std::endl;
              for ( int m = 0; m < ( 1u << num_vars ); ++m )
              {
                tmp_l ^= tmp_l;
                tmp |= ~tmp;
                for ( int v = 0; v < suppB.size(); ++v )
                {
                  if ( ((m >> v)&0x1)==0x1 )
                    tmp_l &= _xl[suppB[v]];
                  else
                    tmp_l &= ~_xl[suppB[v]];
                }
                if ( kitty::count_ones( tmp_l & careB ) > 0 )
                {
                  if ( kitty::count_ones( tmp_l & funcB ) > 0 )
                  {
                    simB |= tmp_l;
                    kitty::set_bit( ttB, m );
                  }
                  else
                    kitty::clear_bit( ttB, m );
                }
              }
std::cout << "q1" << std::endl;

              /* define the support literals */
              std::array<uint32_t, num_vars> childrenA;
              std::array<uint32_t, num_vars> childrenB;
              for ( int v = 0; v < num_vars; ++v )
              {
                childrenA[v] = ( v < suppA.size() ) ? indices[suppA[v]] + offset : null_lit; 
                childrenB[v] = ( v < suppB.size() ) ? indices[suppB[v]] + offset : null_lit; 
              }
              lit_t litA = nodes.size() + offset;
              nodes.emplace_back( childrenA, ttA, litA );
              lit_t litB = nodes.size() + offset;
              nodes.emplace_back( childrenB, ttB, litB );
              children_free[num_vars - 2] = litA;
              children_free[num_vars - 1] = litB;
              tts[num_vars-2] = simA;
              tts[num_vars-1] = simB;
std::cout << "q2" << std::endl;

              kitty::static_truth_table<num_vars> ttF;
              for ( int m = 0; m < ( 1u << num_vars ); ++m )
              {
                tmp_l ^= tmp_l;
                tmp |= ~tmp;
                for ( int v = 0; v < num_vars; ++v )
                {
                  if ( ((m >> v)&0x1)==0x1 )
                    tmp_l &= tts[v];
                  else
                    tmp_l &= ~tts[v];
                }
                if ( kitty::count_ones( tmp_l & func ) > 0 )
                {
                  kitty::set_bit( ttF, m );
                }
                else
                  kitty::clear_bit( ttF, m );
              }
std::cout << "q3" << std::endl;
              lit_t litF = nodes.size() + offset;
              nodes.emplace_back( children_free, ttF, litF );
              return litF;
          }
        }
      }
    }
    }

//std::cout << "get out " << std::endl;
    return null_lit;
  }

  lit_t delay_oriented_satdec( kitty::static_truth_table<num_vars_l> tt, std::vector<double> times, std::vector<uint32_t> indices )
  {
    spfd_manager.reset();
    spfd_manager.init( tt );

    if ( kitty::is_const0( tt ) )
      return 0;
    if ( kitty::is_const0( ~tt ) )
      return 1;

    auto supp = kitty::min_base_inplace( tt );
    if ( (supp.size() > indices.size()) || ( supp.back() >= indices.size() ) )
      return null_lit;
    int j = supp.size()-1;
    for ( int i = indices.size() - 1; i >= 0; --i )
    {
      if ( (j >= 0) && (i == static_cast<int>( supp[j])) )
      {
        j -= 1;
      }
      else
      {
        indices.erase( indices.begin() + i );
        times.erase( times.begin() + i );
      }
    }

    for ( int i = indices.size() - 1; i >= 0; --i )
    {
      auto tt0 = kitty::cofactor0( tt, i );
      auto tt1 = kitty::cofactor1( tt, i );
      if ( kitty::equal( tt0, tt1 ) )
      {
        indices.erase( indices.begin() + i );
        times.erase( times.begin() + i );
      }
    }

    /* here the truth table is reduced and the support minimized */
    if ( indices.size() <= num_vars )
    {
      kitty::static_truth_table<num_vars> func;
      kitty::shrink_to_inplace( func, tt );
      std::array<lit_t, num_vars> children;
      for ( int i = 0; i < num_vars; ++i )
      {
        if ( i < indices.size() )
        {
          children[i] = indices[i] + offset;
        }
        else
        {
          children[i] = null_lit;
        }
      }
      lit_t lit = nodes.size() + offset;
      nodes.emplace_back( children, func, lit );
      if ( lit > ( num_vars_l + 2 + 5 ) )
        return null_lit;
      return lit;
    }
    else if ( indices.size() == ( num_vars + 1 ) )
    {
      /* find the index with the worst delay */
      spfd_manager.reset();
      spfd_manager.init( tt );
      std::vector<std::pair<uint32_t, double>> idx_time;
      std::array<uint32_t, num_vars> bound_indices;
      std::array<uint32_t, num_vars> free_indices;
      for ( uint32_t i = 0; i < indices.size(); ++i )
      {
        idx_time.push_back( std::make_pair( i, times[i] ) );
      }
      std::sort( idx_time.begin(), idx_time.end(), [&]( auto a, auto b ){
        return a.second > b.second;
      } );
      std::array<kitty::static_truth_table<num_vars_l>, num_vars_l> xs;
      for ( int i = 0; i < num_vars_l; ++i )
      {
        kitty::create_nth_var( xs[i], i );
      }
      for ( int i = 0; i < num_vars - 1; ++i )
      {
        free_indices[i] = indices[idx_time[i].first] + offset;
        spfd_manager.update( xs[idx_time[i].first] );
      }
      for ( int a = 0; a < ( 1u << spfd_manager.num_masks )/2; ++a )
      {
        kitty::static_truth_table<num_vars_l> tt_new;
        std::cout << "a " << a << std::endl;
        for ( int m = 0; m < spfd_manager.num_masks; ++m )
        {
          auto const& mask = spfd_manager.masks[m];
          if ( ( ( a >> m ) & 0x1 ) == 0x1 )
            tt_new |= mask & tt;
          else
            tt_new |= mask & ~tt;
        }
        std::cout << "b " << a << std::endl;

        auto supp_b = kitty::min_base_inplace( tt_new );
        for ( int k = supp_b.size()-1; k >= 0; --k )
        {
          auto ttn0 = kitty::cofactor0( tt_new, k );
          auto ttn1 = kitty::cofactor1( tt_new, k );
          if( kitty::equal( ttn0, ttn1 ) )
          {
            supp_b.erase( supp_b.begin() + k );
          }
        }
        std::cout << "c " << a << std::endl;

        if ( supp_b.size() <= num_vars )
        {
        std::cout << "d " << a << std::endl;

          for ( int i = 0; i < supp_b.size(); ++i )
            bound_indices[i] = indices[supp_b[i]] + offset;
          for ( int i = supp_b.size(); i < num_vars; ++i )
            bound_indices[i] = null_lit;

          lit_t lit = nodes.size() + offset;
          kitty::static_truth_table<num_vars> tt_bound;
          kitty::shrink_to_inplace( tt_bound, tt_new );
          nodes.emplace_back( bound_indices, tt_bound, lit );
          free_indices[num_vars-1]=lit;
          std::array<kitty::static_truth_table<num_vars_l>, num_vars> tts;
          for( int i = 0; i < num_vars-1; ++i )
            tts[i] = xs[idx_time[i].first];
        std::cout << "e " << a << std::endl;

          tts[num_vars-1] = tt_new;
          kitty::static_truth_table<num_vars> res;
          kitty::static_truth_table<num_vars_l> tmp;
          kitty::static_truth_table<num_vars> tt_free;
        std::cout << "f " << a << std::endl;

          for ( int min = 0; min < ( 1u << num_vars ); ++min )
          {
            tmp |= ~tmp;
            for ( int i = 0; i < num_vars; ++i )
            {
              if ( ( ( min >> i ) & 0x1 ) == 0x1 )
                tmp &= tts[i];
              else
                tmp &= ~tts[i];
            }
            if ( kitty::count_ones( tmp & tt ) > 0 )
              kitty::set_bit( tt_free, min );
            else
              kitty::clear_bit( tt_free, min );
            lit = nodes.size() + offset;
            nodes.emplace_back( free_indices, tt_free, lit );
          }
        std::cout << "g " << lit << std::endl;
          if ( lit > ( num_vars_l + 2 + 5 ) )
            return null_lit;
          return lit;
        }
      }
        std::cout << "h " << std::endl;

      return null_lit;
    }
    else if ( ( indices.size() > ( num_vars + 1 ) ) && ( indices.size() <= ( 2*num_vars + 2 ) ) )
    {
      /* find the index with the worst delay */
      spfd_manager.reset();
      spfd_manager.init( tt );
      std::vector<std::pair<uint32_t, double>> idx_time;
      std::array<uint32_t, num_vars> bound_indices;
      std::array<uint32_t, num_vars> free_indices;
      for ( uint32_t i = 0; i < indices.size(); ++i )
      {
        idx_time.push_back( std::make_pair( i, times[i] ) );
      }
      std::sort( idx_time.begin(), idx_time.end(), [&]( auto a, auto b ){
        return a.second > b.second;
      } );
      std::array<kitty::static_truth_table<num_vars_l>, num_vars_l> xs;
      for ( int i = 0; i < num_vars_l; ++i )
      {
        kitty::create_nth_var( xs[i], i );
      }
      free_indices[0] = indices[idx_time[0].first] + offset;
      spfd_manager.update( xs[idx_time[0].first] );
      free_indices[1] = indices[idx_time[1].first] + offset;
      spfd_manager.update( xs[idx_time[1].first] );

      auto tt00 = kitty::cofactor0( kitty::cofactor0( tt, idx_time[0].first ), idx_time[1].first );
      auto tt01 = kitty::cofactor0( kitty::cofactor1( tt, idx_time[0].first ), idx_time[1].first );
      auto tt10 = kitty::cofactor1( kitty::cofactor0( tt, idx_time[0].first ), idx_time[1].first );
      auto tt11 = kitty::cofactor1( kitty::cofactor1( tt, idx_time[0].first ), idx_time[1].first );
      kitty::static_truth_table<num_vars_l> ttA;
      kitty::static_truth_table<num_vars_l> ttB;
      std::vector<uint32_t> indicesA;
      std::vector<uint32_t> indicesB;
      std::vector<uint8_t> suppA;
      std::vector<uint8_t> suppB;
      std::vector<double> timesA;
      std::vector<double> timesB;
      auto const& masks = spfd_manager.masks;
      for ( int a = 0; a < 1; a++ )
      {
        if ( a == 0 )
        {
          ttA = ( tt00 & ( masks[0] | masks[1] ) ) | ( tt01 & ( masks[2] | masks[3] ) );
          ttB = ( tt10 & ( masks[0] | masks[1] ) ) | ( tt11 & ( masks[2] | masks[3] ) );
          suppA = kitty::min_base_inplace( ttA );
          suppB = kitty::min_base_inplace( ttB );
        }

        for ( int i = 0; i < suppA.size(); ++i )
        {
          indicesA.push_back( indices[suppA[i]] );
          timesA.push_back( times[suppA[i]] );
        }
        for ( int i = 0; i < suppB.size(); ++i )
        {
          indicesB.push_back( indices[suppB[i]] );
          timesB.push_back( times[suppB[i]] );
        }

        auto litA = delay_oriented_satdec( ttA, timesA, indicesA );
        auto litB = delay_oriented_satdec( ttB, timesB, indicesB );
        if ( (litA == null_lit)||(litB == null_lit) )
          return null_lit;
        free_indices[2] = litA;
        free_indices[3] = litB;
        std::array<kitty::static_truth_table<num_vars_l>, num_vars> xs;
        kitty::create_nth_var( xs[0], idx_time[0].first );
        kitty::create_nth_var( xs[1], idx_time[1].first );
        xs[2] = ttA;
        xs[3] = ttB;
        kitty::static_truth_table<num_vars_l> tmp;
        kitty::static_truth_table<num_vars> func;
        for ( int m = 0; m < ( 1u << num_vars ); ++m )
        {
          tmp |= ~tmp;
          for ( int v = 0; v < num_vars; ++v )
          {
            if ( ( ( m >> v ) & 0x1 ) == 0x1 )
              tmp &= xs[v];
            else
              tmp &= ~xs[v];
          }
          if ( kitty::count_ones( tmp ) > 0 )
            kitty::set_bit( func, m );
          else
            kitty::clear_bit( func, m );
        }
        uint32_t lit = nodes.size() + offset;
        nodes.emplace_back( free_indices, func, lit );
        if ( lit > ( num_vars_l + 2 + 5 ) )
            return null_lit;
        return lit;
      }
      return null_lit;
    }
    /* find the index with the worst delay */
    double worst_time = std::numeric_limits<double>::min();
    uint32_t worst_index = 0;
    for ( int i = 0; i < indices.size(); ++i )
    {
      if ( times[i] > worst_time )
      {
        worst_time = times[i];
        worst_index = i;
      }
    }

    lit_t lit_idx = indices[worst_index] + offset;
    auto tt0 = kitty::cofactor0( tt, worst_index );
    auto tt1 = kitty::cofactor1( tt, worst_index );
    if ( kitty::equal( tt0, ~tt1 ) )
    {
      lit_t l0 = delay_oriented_satdec( tt0, times, indices );

      kitty::static_truth_table<num_vars> x0;
      kitty::static_truth_table<num_vars> x1;
      kitty::static_truth_table<num_vars> fn;
      kitty::create_nth_var( x0, 0 );
      kitty::create_nth_var( x1, 1 );
      fn = x0 ^ x1;
      lit_t lit = nodes.size() + offset;
      std::array<lit_t, num_vars> children{ lit_idx, l0, null_lit, null_lit };
      nodes.emplace_back( children, fn, lit );
        if ( lit > ( num_vars_l + 2 + 5 ) )
            return null_lit;
      return lit;
    }
    lit_t l0 = delay_oriented_satdec( tt0, times, indices );
    lit_t l1 = delay_oriented_satdec( tt1, times, indices );
    if ( ( l0 == null_lit ) || ( l1 == null_lit ) )
      return null_lit;

    if ( l0 == 0 ) // F = x F1
    {
      kitty::static_truth_table<num_vars> x0;
      kitty::static_truth_table<num_vars> x1;
      kitty::static_truth_table<num_vars> fn;
      kitty::create_nth_var( x0, 0 );
      kitty::create_nth_var( x1, 1 );
      fn = x0 & x1;
      lit_t lit = nodes.size() + offset;
      std::array<lit_t, num_vars> children{ lit_idx, l1, null_lit, null_lit };
      nodes.emplace_back( children, fn, lit );
        if ( lit > ( num_vars_l + 2 + 5 ) )
            return null_lit;
      return lit;
    }
    else if ( l0 == 1 ) // F = x' + F1
    {
      kitty::static_truth_table<num_vars> x0;
      kitty::static_truth_table<num_vars> x1;
      kitty::static_truth_table<num_vars> fn;
      kitty::create_nth_var( x0, 0 );
      kitty::create_nth_var( x1, 1 );
      fn = ( ~x0 ) | x1;
      lit_t lit = nodes.size() + offset;
      std::array<lit_t, num_vars> children{ lit_idx, l1, null_lit, null_lit };
      nodes.emplace_back( children, fn, lit );
        if ( lit > ( num_vars_l + 2 + 5 ) )
            return null_lit;
      return lit;
    }
    else if ( l1 == 0 ) // F = x'F0
    {
      kitty::static_truth_table<num_vars> x0;
      kitty::static_truth_table<num_vars> x1;
      kitty::static_truth_table<num_vars> fn;
      kitty::create_nth_var( x0, 0 );
      kitty::create_nth_var( x1, 1 );
      fn = ( ~x0 ) & x1;
      lit_t lit = nodes.size() + offset;
      std::array<lit_t, num_vars> children{ lit_idx, l0, null_lit, null_lit };
      nodes.emplace_back( children, fn, lit );
        if ( lit > ( num_vars_l + 2 + 5 ) )
            return null_lit;
      return lit;
    }
    else if ( l1 == 1 ) // F = x + F0
    {
      kitty::static_truth_table<num_vars> x0;
      kitty::static_truth_table<num_vars> x1;
      kitty::static_truth_table<num_vars> fn;
      kitty::create_nth_var( x0, 0 );
      kitty::create_nth_var( x1, 1 );
      fn = x0 | x1;
      lit_t lit = nodes.size() + offset;
      std::array<lit_t, num_vars> children{ lit_idx, l0, null_lit, null_lit };
      nodes.emplace_back( children, fn, lit );
        if ( lit > ( num_vars_l + 2 + 5 ) )
            return null_lit;
      return lit;
    }
    else
    {
      kitty::static_truth_table<num_vars> x0;
      kitty::static_truth_table<num_vars> x1;
      kitty::static_truth_table<num_vars> x2;
      kitty::static_truth_table<num_vars> fn;
      kitty::create_nth_var( x0, 0 );
      kitty::create_nth_var( x1, 1 );
      kitty::create_nth_var( x2, 2 );
      fn = ( x0 & x1 ) | ( ( ~x0 ) & x2 );
      lit_t lit = nodes.size() + offset;
      std::array<lit_t, num_vars> children{ lit_idx, l1, l0, null_lit };
      nodes.emplace_back( children, fn, lit );
      if ( lit > ( num_vars_l + 2 + 5 ) )
          return null_lit;
      return lit;
    }  
  }

  void print()
  {
    for ( auto nd : nodes )
    {
      for ( int i = 0; i < num_vars; ++i )
        std::cout << nd.children[i] << " ";
      kitty::print_binary( nd.func );
      std::cout << " : " << nd.lit << std::endl;
    }
    std::cout << std::endl;
  }

private:
  std::vector<node_t> nodes;
  std::set<uint32_t> done;
  spfd_manager_t<num_vars_l, num_vars_l> spfd_manager;
  std::array<kitty::static_truth_table<num_vars_l>, num_vars_l+100> _xl;
  std::array<kitty::static_truth_table<num_vars>, num_vars> _xs;

  std::random_device rd;
  std::mt19937 gen;

};

} /* namespace mockturtle */