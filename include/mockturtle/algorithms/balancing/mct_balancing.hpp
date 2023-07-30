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
  \file mct_balancing.hpp
  \brief SOP-based balancing engine for `balancing` algorithm

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
#include <kitty/hash.hpp>
#include <kitty/isop.hpp>
#include <kitty/operations.hpp>

#include <mockturtle/algorithms/decompose/DecSolver.hpp>
#include <mockturtle/algorithms/ccgame/solvers/cusco.hpp>

#include "../../traits.hpp"
#include "../../utils/stopwatch.hpp"
#include "../balancing.hpp"
#include "utils.hpp"

namespace mockturtle
{

/*! \brief SOP rebalancing function
 *
 * This class can be used together with the generic `balancing` function.  It
 * converts each cut function into an SOP and then performs weight-oriented
 * tree balancing on the AND terms and the outer OR function.
 */
template<class Ntk>
struct mct_rebalancing
{
  void operator()( Ntk& dest, kitty::dynamic_truth_table const& function, std::vector<arrival_time_pair<Ntk>> const& inputs, uint32_t best_level, uint32_t best_cost, rebalancing_function_callback_t<Ntk> const& callback ) const
  {
    auto [and_terms, num_and_gates] = create_function( dest, function, inputs );
    const auto num_gates = num_and_gates + ( and_terms.empty() ? 0u : static_cast<uint32_t>( and_terms.size() ) - 1u );
    const auto cand = balanced_tree( dest, and_terms, false );

  //  ADDON

  ccgame::report_t<Ntk> rep;
  typedef DecSolver<kitty::dynamic_truth_table, Ntk> solver_t;

  kitty::dynamic_truth_table mask = function.construct();
  mask = mask | ~mask;
  solver_t solver( {function}, {mask} );

    std::vector<kitty::dynamic_truth_table> xs;
    for( int i{0}; i<function.num_vars(); ++i )
    {
      xs.emplace_back( function.num_vars() );
      kitty::create_nth_var( xs[i], i );
    }
    /* define the parameters */
    int nIters = 10;
    ccgame::cusco_ps cusco_ps( ccgame::solver_t::_SYM_RDE, nIters );
    std::vector<signal<Ntk>> S;
    for( int i{0}; i < inputs.size(); ++i )
    {
      cusco_ps.T.push_back( inputs[i].level );
      S.push_back( inputs[i].f );
    }
    /* solve */
    ccgame::cusco<Ntk> csolver( xs, {function} );
//printf("w0\n");
    rep = csolver.solve( cusco_ps, S, &dest );
//printf("w1\n");
    arrival_time_pair<Ntk> csym{ rep.osig, rep.levels };
//printf("w2\n");

    //printf("%d SOP = %d SYM = %d\n", rep.Esl, cand.level, csym.level );
    if( !rep.Esl || cand.level < csym.level )
    {
      if ( cand.level < best_level || ( cand.level == best_level && num_gates < best_cost ) )
      {
        callback( cand, num_gates );
      }
    }
    else if( cand.level > csym.level )
    {
      if ( csym.level < best_level || ( csym.level == best_level && rep.ntk.num_gates() < best_cost ) )
      {
        callback( csym, rep.ntk.num_gates() );
      }
    }
    else if( cand.level == csym.level && ( num_gates < rep.ntk.num_gates() ) )
    {
      if ( cand.level < best_level || ( cand.level == best_level && num_gates < best_cost ) )
      {
        callback( cand, num_gates );
      }
    }
    else if( cand.level == csym.level && ( num_gates == rep.ntk.num_gates() ) )
    {
      if ( csym.level < best_level || ( csym.level == best_level && rep.ntk.num_gates() < best_cost ) )
      {
        callback( csym, rep.ntk.num_gates() );
      }
    }
  //  END ADDON

  /*  if( !rep.Esl )//|| cand.level < csym.level )
    {
      if ( cand.level < best_level || ( cand.level == best_level && num_gates < best_cost ) )
      {
        callback( cand, num_gates );
      }
    }
    //else if( cand.level > csym.level )
    else
    {
      if ( csym.level < best_level || ( csym.level == best_level && rep.ntk.num_gates() < best_cost ) )
      {
        callback( csym, rep.ntk.num_gates() );
      }
    }
//    else if( cand.level == csym.level && ( num_gates < rep.ntk.num_gates() ) )
//    {
//      if ( cand.level < best_level || ( cand.level == best_level && num_gates < best_cost ) )
//      {
//        callback( cand, num_gates );
//      }
//    }
//    else if( cand.level == csym.level && ( num_gates == rep.ntk.num_gates() ) )
//    {
//      if ( csym.level < best_level || ( csym.level == best_level && rep.ntk.num_gates() < best_cost ) )
//      {
//        callback( csym, rep.ntk.num_gates() );
//      }
//    }
*/


  }

private:
  std::pair<arrival_time_queue<Ntk>, uint32_t> create_function( Ntk& dest, kitty::dynamic_truth_table const& func, std::vector<arrival_time_pair<Ntk>> const& arrival_times ) const
  {
    const auto sop = create_sop_form( func );

    stopwatch<> t_tree( time_tree_balancing );
    arrival_time_queue<Ntk> and_terms;
    uint32_t num_and_gates{};
    for ( auto const& cube : sop )
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
      if ( !product_queue.empty() )
      {
        num_and_gates += static_cast<uint32_t>( product_queue.size() ) - 1u;
      }
      and_terms.push( balanced_tree( dest, product_queue ) );
    }
    return { and_terms, num_and_gates };
  }

  arrival_time_pair<Ntk> balanced_tree( Ntk& dest, arrival_time_queue<Ntk>& queue, bool _and = true ) const
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
      const auto s = _and ? dest.create_and( s1, s2 ) : dest.create_or( s1, s2 );
      const auto l = std::max( l1, l2 ) + 1;
      queue.push( { s, l } );
    }
    return queue.top();
  }

  std::vector<kitty::cube> create_sop_form( kitty::dynamic_truth_table const& func ) const
  {
    stopwatch<> t( time_sop );
    if ( auto it = sop_hash_.find( func ); it != sop_hash_.end() )
    {
      sop_cache_hits++;
      return it->second;
    }
    else
    {
      sop_cache_misses++;
      return sop_hash_[func] = kitty::isop( func ); 
    }
  }

private:
  mutable std::unordered_map<kitty::dynamic_truth_table, std::vector<kitty::cube>, kitty::hash<kitty::dynamic_truth_table>> sop_hash_;

public:
  mutable uint32_t sop_cache_hits{};
  mutable uint32_t sop_cache_misses{};

  mutable stopwatch<>::duration time_sop{};
  mutable stopwatch<>::duration time_tree_balancing{};
};

} // namespace mockturtle