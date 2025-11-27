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
  \file timing.hpp
  \brief Compute timing information of a network.

  This engine can be used for power analysis of mapped network.
  For each node, the following information is stored:
  - The sensing time : the first time at which a transition can happen
  - The sensing time : the first time at which the output is stable
  - A vector of simulation patterns identifying quantized timesteps in this interval
  The class has the following template parameters:
  \param N : number of timesteps used for simulating each signal
  \param I : number of inputs to be used for simulating the static truth table ( 2^I input pairs )

  \author Andrea Costamagna
*/

#pragma once

#include "../../utils/node_map.hpp"
#include "../../algorithms/detail/resub_utils.hpp"
#include "../../algorithms/reconv_cut.hpp"
#include "../../algorithms/cut_enumeration.hpp"
#include "../../algorithms/cut_enumeration/rewrite_cut.hpp"

namespace mockturtle
{

struct window_manager_params
{
  window_manager_params()
  {
    /* 0 < Cut limit < 16 */
    cut_enumeration_ps.cut_limit = 8;//12
    /* otherwise creates problem in the presence of reconvergence (A^B)^(B^C) becomes {A,C} */
    cut_enumeration_ps.minimize_truth_table = true;
  }

  uint32_t skip_fanout_limit_for_divisors{ 100 };
  uint32_t max_num_divisors{ 128 };
  bool preserve_depth{ false };
  cut_enumeration_params cut_enumeration_ps{};
};

template<class Ntk, uint32_t max_pis = 12>
class window_manager
{
public:
  using node = typename Ntk::node;
  using cut_comp = detail::reconvergence_driven_cut_impl<Ntk>;
  using cut_comp_parameters_type = typename cut_comp::parameters_type;
  using cut_comp_statistics_type = typename cut_comp::statistics_type;
  using network_cuts_t = dynamic_network_cuts<Ntk, max_pis, true, cut_enumeration_rewrite_cut>;
  using cut_manager_t = detail::dynamic_cut_enumeration_impl<Ntk, max_pis, true, cut_enumeration_rewrite_cut>;
  using cut_t = typename network_cuts_t::cut_t;

public:
  window_manager( Ntk & ntk, window_manager_params ps = {} )
    : ntk( ntk ),
      ps( ps ),
      cuts_st(),
      cuts( ntk, cut_comp_parameters_type{ max_pis }, cuts_st ),
      cuts2( ntk.size() + ( ntk.size() >> 1 ) ),
      cut_manager( ntk, ps.cut_enumeration_ps, cst2, cuts2 ),
      valid( false )
  {
    cut_manager.init_cuts2();
  }

  void reset()
  {
    root = ntk.get_null();
    divs.clear();
    mffc.clear();
    leaves.clear();
    valid = false;
  }

  bool run( node const& n )
  {
    root = n;
    divs.clear();
    mffc.clear();
    leaves.clear();
    valid = true;

    /* set the flags for storage: visited should be trav_id(),
     * value should be
     * leaves : 0 | ( trav_id() << 2u )
     * mffc : 1 | ( trav_id() << 2u )
     * divs : 2 | ( trav_id() << 2u )
     */
    ntk.incr_trav_id();
    /* find a reconvergence-driven cut */
    if ( ntk.fanin_size( n ) == 1 )
    {
      auto const& children = ntk.get_children( n );
      auto ni = ntk.get_node( children[0] );
      if ( ntk.is_pi( ni ) )
      {
        leaves = { ni };
        divs = { ni };
        mffc = { n };
        return true;
      }
      else
      {
        leaves = cuts.run( { ni } ).first;
      }
    }
    else
    {
      leaves = cuts.run( { n } ).first;
      if ( leaves.size() == 1 )
        leaves = ntk.get_nodes_fanin( n );
    }
    for( auto const& l : leaves )
    {
      make_leaf( l );
      // std::cout << l << "_l ";
      divs.emplace_back( l );
    }
    // std::cout << std::endl;
    
    
    /* collect the maximum fanout free cone */
    if( !collect_mffc( n ) )
    {
      valid = false;
      std::cout << "failed collecting mffc" << std::endl;
      return false;
    }
    
    /* collect the divisors in the tfi */
    if( !collect_divs_tfi( n ) )
    {
      valid = false;
      std::cout << "failed collecting divs tfi" << std::endl;
      
      return false;
    }
    /* TODO: expand divisors */
    if( !collect_divs_sides( n ) )
    {
      valid = false;
      std::cout << "failed collecting divs sides" << std::endl;

      return false;
    }
    return true;
  }

  std::vector<node> const & get_divisors() const
  {
    std::vector<node> const & res = divs;
    return res;
  }

  bool is_leaf( node const& n ) const
  {
    return ntk.value( n ) == ( 3u | ( ntk.trav_id() << 2u ) );
  }

  void make_leaf( node const& n ) const
  {
    ntk.set_value( n, ( 3u | ( ntk.trav_id() << 2u ) ));
  }

  bool is_mffc( node const& n ) const
  {
    return ntk.value( n ) == ( 1u | ( ntk.trav_id() << 2u ) );
  }

  void make_mffc( node const& n ) const
  {
    ntk.set_value( n, ( 1u | ( ntk.trav_id() << 2u ) ));
  }

  bool is_divisor( node const& n ) const
  {
    return ntk.value( n ) == ( 2u | ( ntk.trav_id() << 2u ) );
  }

  void make_divisor( node const& n ) const
  {
    ntk.set_value( n, ( 2u | ( ntk.trav_id() << 2u ) ));
  }

  std::vector<node> const& get_leaves() const
  {
    return leaves;
  }

  std::vector<node> const& get_divs() const
  {
    return divs;
  }

  std::vector<node> const& get_mffc() const
  {
    return mffc;
  }

  node const& get_root() const
  {
    return root;
  }

#pragma region Mffc
private:

  bool collect_mffc( node const& n )
  {
    for ( const auto& l : leaves )
      ntk.incr_fanout_size( l );
    /* dereference the node */
    mffc.emplace_back( n );
    make_mffc( n );
    node_deref_rec( n );
    std::reverse( mffc.begin(), mffc.end() );
    for( auto m : mffc )
    {
      if( ntk.is_pi( m ) )
        return false;
    }
    /* reference it back */
    node_ref_rec( n );
    for ( const auto& l : leaves )
      ntk.decr_fanout_size( l );

    return true;
  }

  /* ! \brief Dereference the node's MFFC */
  void node_deref_rec( node const& n )
  {
    if ( ntk.is_pi( n ) )
      return;

//std::cout << n << "_m ";
//std::cout << std::endl;

    ntk.foreach_fanin( n, [&]( const auto& f ) {
      auto const& p = ntk.get_node( f );

      ntk.decr_fanout_size( p );
      if ( ntk.fanout_size( p ) == 0 )
      {
        mffc.emplace_back( p );
        make_mffc( p );
        node_deref_rec( p );
      }
    } );
//for( auto m : mffc )
//  std::cout << m << ".. ";
//std::cout << std::endl;
  }

  /* ! \brief Reference the node's MFFC */
  void node_ref_rec( node const& n )
  {
    if ( ntk.is_pi( n ) )
      return;

    ntk.foreach_fanin( n, [&]( const auto& f ) {
      auto const& p = ntk.get_node( f );

      auto v = ntk.fanout_size( p );
      ntk.incr_fanout_size( p );
      if ( v == 0 )
      {
        node_ref_rec( p );
      }
    } );
  }
#pragma endregion

#pragma region Divisors
  bool collect_divs_tfi( node const& n )
  {    
    if( is_leaf( n ) )
    {
      return true;
    }
    if( is_divisor( n ) )
    {
      return true;
    }
    /* PI and not leaf means error */
    if( ntk.is_pi( n ) )
    {
      std::cout << "pi not leaf" << std::endl;
      return false;
    }
    bool is_good = true;
    ntk.foreach_fanin( n, [&]( auto fi ){
      node const ni = ntk.get_node( fi );
      is_good &= collect_divs_tfi( ni );
    } );
    if( !is_mffc( n ) && !is_leaf( n ) )
    {
      divs.emplace_back( n );
      make_divisor( n );
    }
    return is_good;
  }

  bool collect_divs_sides( node const& n )
  {
    bool quit = false;
    bool is_good = true;
    for ( auto i = 0u; i < divs.size(); ++i )
    {
      auto const d = divs.at( i );

      if ( ntk.fanout_size( d ) > ps.skip_fanout_limit_for_divisors )
      {
        continue;
      }
      if ( divs.size() >= ps.max_num_divisors )
      {
        break;
      }

      /* if the fanout has all fanins in the set, add it */
      bool add = true;
      ntk.foreach_fanout( d, [&]( node const& p ) {
        if ( is_leaf( p ) || is_divisor( p ) || is_mffc( p ) )
        {
          return true; /* next fanout */
        }
        if ( ntk.fanout_size( p ) <= 0 )
        {
          return true; /* next fanout */
        }

        bool all_fanins_visited = true;
        ntk.foreach_fanin( p, [&]( const auto& g ) {
          if ( !is_leaf( ntk.get_node(g) ) && !is_divisor( ntk.get_node(g) ) )
          {
            all_fanins_visited = false;
            return false; /* terminate fanin-loop */
          }
          return true; /* next fanin */
        } );

        if ( !all_fanins_visited )
        {
          add = false;
          return true; /* next fanout */
        }

        bool has_root_as_child = false;
        ntk.foreach_fanin( p, [&]( const auto& g ) {
          if ( ntk.get_node( g ) == root )
          {
            has_root_as_child = true;
            return false; /* terminate fanin-loop */
          }
          return true; /* next fanin */
        } );

        if ( has_root_as_child )
        {
          std::cout << "has root as child " << root << std::endl;
          return true; /* next fanout */
        }
        is_good &= all_fanins_visited;
        divs.emplace_back( p );
        make_divisor( p );

        /* quit computing divisors if there are too many of them */
        if ( divs.size() >= ps.max_num_divisors )
        {
          quit = true;
          return false; /* terminate fanout-loop */
        }

        return true; /* next fanout */
      } );

      if ( quit )
      {
        break;
      }
    }
    return is_good;
  }
#pragma endregion

public:
  /*! \brief iterator to apply a lambda to all the nodes */
  template<typename Fn>
  void foreach_leaf( Fn&& fn ) const
  {
    for ( uint32_t i = 0; i < leaves.size(); i ++ )
    {
      fn( leaves[i], i );
    }
  }

  /*! \brief iterator to apply a lambda to all the nodes */
  template<typename Fn>
  void foreach_divisor( Fn&& fn ) const
  {
    for ( uint32_t i = 0; i < divs.size(); i ++ )
    {
      fn( divs[i], i );
    }
  }

  /*! \brief iterator to apply a lambda to all the nodes */
  template<typename Fn>
  void foreach_mffc( Fn&& fn ) const
  {
    for ( uint32_t i = 0; i < mffc.size(); i++ )
    {
      fn( mffc[i], i );
    }
  }

  void print() const
  {
    std::cout << "leaves" << std::endl;
    foreach_leaf( [&]( auto l, auto i ){
      std::cout << l << " ";
    } );
    std::cout << std::endl;
    std::cout << "divs" << std::endl;
    foreach_divisor( [&]( auto d, auto i ){
      std::cout << d << " ";
    } );
    std::cout << std::endl;
    std::cout << "mffc" << std::endl;
    foreach_mffc( [&]( auto m, auto i ){
      std::cout << m << " ";
    } );
    std::cout << std::endl;
  }

  bool is_valid() const
  {
    return valid;
  }

private:
  Ntk & ntk;
  network_cuts_t cuts2;
  window_manager_params ps;
  cut_comp_statistics_type cuts_st;
  cut_comp cuts;
  cut_enumeration_stats cst2;
  cut_manager_t cut_manager;



public:
  std::vector<node> leaves;
  std::vector<node> divs;
  std::vector<node> mffc;
  node root;
  bool valid; 

};

} // namespace mockturtle