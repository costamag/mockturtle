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

/*!
  \file esop_balancing.hpp
  \brief ESOP-based balancing engine for `balancing` algorithm

  \author Alessandro Tempia Calvino
  \author Heinz Riener
  \author Mathias Soeken
*/

#pragma once

#include <algorithm>
#include <cstdint>
#include <queue>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <kitty/cube.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/esop.hpp>
#include <kitty/hash.hpp>
#include <kitty/operations.hpp>
#include <kitty/spp.hpp>

#include "../../traits.hpp"
#include "../../utils/stopwatch.hpp"
#include "../balancing.hpp"
#include "../exorcism.hpp"
#include "../../algorithms/techaware/sym_synthesis.hpp"
#include "utils.hpp"

namespace mockturtle
{

using namespace techaware;

template<class Ntk>
struct sym_rebalancing
{
  void operator()( Ntk& dest, kitty::dynamic_truth_table const& function, std::vector<arrival_time_pair<Ntk>> const& inputs, uint32_t best_level, uint32_t best_cost, rebalancing_function_callback_t<Ntk> const& callback ) const
  {
    bool inverted = false;
    auto [and_terms, max_level, num_and_gates] = create_esop_function( dest, function, inputs, inverted );
    arrival_time_pair<Ntk> cand = balanced_xor_tree( dest, and_terms );

    auto [csym, num_sym_gates, symError] = create_symm_function( dest, function, inputs );

    //printf("d(SYM)=%d a(SYM)=%d, d(ESOP)=%d a(ESOP)=%d\n", csym.level, num_sym_gates, cand.level, num_and_gates );

    //if ( symError && ( cand.level < best_level || ( cand.level == best_level && num_and_gates < best_cost ) ) )
    //{
    //  cand.f = cand.f ^ inverted;
    //  callback( cand, num_and_gates );
    //}
    if ( !symError && (csym.level < best_level || ( csym.level == best_level && num_sym_gates < best_cost ) ))
    {
      callback( csym, num_sym_gates );
    }
  }

private:
  std::tuple<arrival_time_queue<Ntk>, uint32_t, uint32_t> create_esop_function( Ntk& dest, kitty::dynamic_truth_table const& func, std::vector<arrival_time_pair<Ntk>> const& arrival_times, bool& inverted ) const
  {
    if ( spp_optimization )
    {
      return create_function_from_spp( dest, func, arrival_times, inverted );
    }
    else
    {
      return create_function_from_esop( dest, func, arrival_times, inverted );
    }
  }

  std::tuple<arrival_time_pair<Ntk>, uint32_t, bool> create_symm_function( Ntk& dest, kitty::dynamic_truth_table const& func, std::vector<arrival_time_pair<Ntk>> const& arrival_times ) const
  {

    std::vector<uint32_t> T;
    std::vector<signal<Ntk>> S;
    for( auto time : arrival_times )
    {
      T.push_back( time.level );
      S.push_back( time.f );
    }

    sym_synthesis<Ntk> synt( func, T );
    
    signal<Ntk> sigOut = synt.rewrite( &dest, S );
    uint32_t out_level = synt.get_output_level();
    arrival_time_pair<Ntk> csym {sigOut, out_level};

    return { csym, synt.get_num_nodes(), synt.net.error };
  }

  std::tuple<arrival_time_queue<Ntk>, uint32_t, uint32_t> create_function_from_esop( Ntk& dest, kitty::dynamic_truth_table const& func, std::vector<arrival_time_pair<Ntk>> const& arrival_times, bool& inverted ) const
  {
    const auto esop = create_sop_form( func, inverted );

    stopwatch<> t_tree( time_tree_balancing );
    arrival_time_queue<Ntk> and_terms;
    uint32_t max_level{};
    uint32_t num_and_gates{};
    for ( auto const& cube : esop )
    {
      arrival_time_queue<Ntk> product_queue;
      for ( auto i = 0u; i < func.num_vars(); ++i )
      {
        if ( cube.get_mask( i ) )
        {
          const auto [f, l] = arrival_times[i];
          product_queue.push( { cube.get_bit( i ) ? f : dest.create_not( f ), l } );
        }
      }
      if ( product_queue.size() )
      {
        num_and_gates += static_cast<uint32_t>( product_queue.size() ) - 1u;
      }
      arrival_time_pair<Ntk> and_res = balanced_and_tree( dest, product_queue );
      and_terms.push( and_res );
      max_level = std::max( max_level, and_res.level );
    }
    return { and_terms, max_level, num_and_gates };
  }

  std::tuple<arrival_time_queue<Ntk>, uint32_t, uint32_t> create_function_from_spp( Ntk& dest, kitty::dynamic_truth_table const& func, std::vector<arrival_time_pair<Ntk>> const& arrival_times, bool& inverted ) const
  {
    const auto esop = create_sop_form( func, inverted );
    const auto [spp, sums] = kitty::simple_spp( esop, func.num_vars() );

    stopwatch<> t_tree( time_tree_balancing );
    arrival_time_queue<Ntk> and_terms;
    uint32_t max_level{};
    uint32_t num_and_gates{};
    for ( auto const& cube : spp )
    {
      arrival_time_queue<Ntk> product_queue;
      for ( auto i = 0u; i < func.num_vars(); ++i )
      {
        if ( cube.get_mask( i ) )
        {
          const auto [f, l] = arrival_times[i];
          product_queue.push( { cube.get_bit( i ) ? f : dest.create_not( f ), l } );
        }
      }
      for ( auto i = 0u; i < sums.size(); ++i )
      {
        if ( cube.get_mask( func.num_vars() + i ) )
        {
          std::vector<signal<Ntk>> xor_terms;
          uint32_t xor_level{};
          for ( auto j = 0u; j < func.num_vars(); ++j )
          {
            if ( ( sums[i] >> j ) & 1 )
            {
              const auto [f, l] = arrival_times[j];
              xor_terms.push_back( f );
              xor_level = std::max( xor_level, l );
            }
          }
          const auto f = dest.create_nary_xor( xor_terms );
          product_queue.push( { cube.get_bit( func.num_vars() + i ) ? f : dest.create_not( f ), xor_level } );
        }
      }
      if ( product_queue.size() )
      {
        num_and_gates += static_cast<uint32_t>( product_queue.size() ) - 1u;
      }
      arrival_time_pair<Ntk> and_res = balanced_and_tree( dest, product_queue );
      and_terms.push( and_res );
      max_level = std::max( max_level, and_res.level );
    }
    return { and_terms, max_level, num_and_gates };
  }

  arrival_time_pair<Ntk> balanced_and_tree( Ntk& dest, arrival_time_queue<Ntk>& queue ) const
  {
    if ( queue.empty() )
    {
      return { dest.get_constant( true ), 0u };
    }

    while ( queue.size() > 1u )
    {
      auto [s1, l1] = queue.top();
      queue.pop();
      auto [s2, l2] = queue.top();
      queue.pop();
      const auto s = dest.create_and( s1, s2 );
      const auto l = std::max( l1, l2 ) + 1;
      queue.push( { s, l } );
    }
    return queue.top();
  }

  arrival_time_pair<Ntk> balanced_xor_tree( Ntk& dest, arrival_time_queue<Ntk>& queue ) const
  {
    if ( queue.empty() )
    {
      return { dest.get_constant( true ), 0u };
    }

    while ( queue.size() > 1u )
    {
      auto [s1, l1] = queue.top();
      queue.pop();
      auto [s2, l2] = queue.top();
      queue.pop();
      const auto s = dest.create_xor( s1, s2 );
      const auto l = std::max( l1, l2 ) + 1;
      queue.push( { s, l } );
    }
    return queue.top();
  }

  std::vector<kitty::cube> create_sop_form( kitty::dynamic_truth_table const& func, bool& inverted ) const
  {
    stopwatch<> t( time_sop );
    inverted = false;

    if ( auto it = sop_hash_.find( func ); it != sop_hash_.end() )
    {
      sop_cache_hits++;
      return it->second;
    }

    if ( both_phases )
    {
      if ( auto it = sop_hash_.find( ~func ); it != sop_hash_.end() )
      {
        inverted = true;
        sop_cache_hits++;
        return it->second;
      }
    }

    sop_cache_misses++;
    std::vector<kitty::cube> sop = mockturtle::exorcism( func );

    if ( both_phases )
    {
      std::vector<kitty::cube> n_sop = mockturtle::exorcism( ~func );

      if ( n_sop.size() < sop.size() )
      {
        inverted = true;
        return sop_hash_[~func] = n_sop;
      }
      else if ( n_sop.size() == sop.size() )
      {
        /* compute literal cost */
        uint32_t lit = 0, n_lit = 0;
        for ( auto const& c : sop )
        {
          lit += c.num_literals();
        }
        for ( auto const& c : n_sop )
        {
          n_lit += c.num_literals();
        }

        if ( n_lit < lit )
        {
          inverted = true;
          return sop_hash_[~func] = n_sop;
        }
      }
    }

    return sop_hash_[func] = sop;
  }

private:
  mutable std::unordered_map<kitty::dynamic_truth_table, std::vector<kitty::cube>, kitty::hash<kitty::dynamic_truth_table>> sop_hash_;

public:
  bool both_phases{ false };
  bool spp_optimization{ false };
  bool mux_optimization{ false };

public:
  mutable uint32_t sop_cache_hits{};
  mutable uint32_t sop_cache_misses{};

  mutable stopwatch<>::duration time_sop{};
  mutable stopwatch<>::duration time_tree_balancing{};
};

} // namespace mockturtle