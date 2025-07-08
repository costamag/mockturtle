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
  \file window_manager.hpp
  \brief Manager of windows for peephole optimization.

  This engine can be used to manage
  \param N : number of timesteps used for simulating each signal
  \param I : number of inputs to be used for simulating the static truth table ( 2^I input pairs )

  \author Andrea Costamagna
*/

#pragma once

#include "../../../utils/node_map.hpp"
#include "../../../views/depth_view.hpp"
#include "../../cut_enumeration.hpp"
#include "../../cut_enumeration/rewrite_cut.hpp"
#include "../../reconv_cut.hpp"

namespace mockturtle
{

struct window_manager_params
{
  window_manager_params()
  {
    /* 0 < Cut limit < 16 */
    cut_enumeration_ps.cut_limit = 8; // 12
    /* otherwise creates problem in the presence of reconvergence (A^B)^(B^C) becomes {A,C} */
    cut_enumeration_ps.minimize_truth_table = true;
  }

  uint32_t skip_fanout_limit_for_divisors{ 100 };
  uint32_t max_num_divisors{ 128 };
  bool preserve_depth{ false };
  cut_enumeration_params cut_enumeration_ps{};
};

template<typename Ntk>
struct window_manager_stats
{
  /* stats of the cut computation engine */
  typename detail::reconvergence_driven_cut_impl<Ntk>::statistics_type cuts_st;
};

template<class Ntk, uint32_t MaxNumLeaves = 12>
class window_manager
{
public:
  using node_index_t = typename Ntk::node;
  using cut_comp = detail::reconvergence_driven_cut_impl<Ntk>;
  using cut_comp_parameters_type = typename cut_comp::parameters_type;

  struct window_t
  {
    using node_index_t = typename Ntk::node;
    static constexpr auto null_node = std::numeric_limits<node_index_t>::max();

    /* leaves of the window */
    std::vector<node_index_t> leaves;
    /* divisors of the window */
    std::vector<node_index_t> divs;
    /* portion of the network removed when the root is removed */
    std::vector<node_index_t> mffc;
    /* root node for the window's construction */
    node_index_t root;
    /* flag set to true when a window is available */
    bool valid = false;

    void reset( node_index_t n = null_node )
    {
      root = n;
      divs.clear();
      mffc.clear();
      leaves.clear();
      valid = n != null_node;
    }

    void init( node_index_t const& n )
    {
      reset( n );
    }
  };

public:
  window_manager( Ntk& ntk, window_manager_params const& ps, window_manager_stats<Ntk>& st )
      : ntk_( ntk ),
        ps_( ps ),
        st_( st ),
        cuts_( ntk, cut_comp_parameters_type{ MaxNumLeaves }, st.cuts_st )
  {
  }

  void reset()
  {
    window_.reset();
  }

  bool run( node_index_t const& n )
  {
    window_.init( n );

    /* set the flags for storage: visited should be trav_id(),
     * value should be
     * leaves : 0 | ( trav_id() << 2u )
     * mffc   : 1 | ( trav_id() << 2u )
     * divs   : 2 | ( trav_id() << 2u )
     */
    ntk_.incr_trav_id();
    // find a reconvergence-driven cut
    if ( ntk_.fanin_size( n ) == 1 )
    {
      auto const& children = ntk_.get_children( n );
      auto ni = ntk_.get_node( children[0] );
      if ( ntk_.is_pi( ni ) )
      {
        window_.leaves = { ni };
        window_.divs = { ni };
        window_.mffc = { n };
        return true;
      }
      else
      {
        window_.leaves = cuts_.run( { ni } ).first;
      }
    }
    else
    {
      window_.leaves = cuts_.run( { n } ).first;
      while ( ( window_.leaves.size() == 1 ) && ( !ntk_.is_pi( n ) ) )
        window_.leaves = ntk_.get_fanins( n );
    }

    for ( auto const& l : window_.leaves )
    {
      make_leaf( l );
      window_.divs.emplace_back( l );
    }

    /* collect the maximum fanout free cone */
    if ( !collect_mffc( n ) )
    {
      window_.valid = false;
      std::cout << "failed collecting mffc" << std::endl;
      return false;
    }

    /* collect the divisors in the tfi */
    if ( !collect_divs_tfi( n ) )
    {
      window_.valid = false;
      std::cout << "failed collecting divs tfi" << std::endl;

      return false;
    }
    /* TODO: expand divisors */
    if ( !collect_divs_sides( n ) )
    {
      window_.valid = false;
      std::cout << "failed collecting divs sides" << std::endl;

      return false;
    }
    return true;
  }

  std::vector<node_index_t> const& get_divisors() const
  {
    std::vector<node_index_t> const& res = window_.divs;
    return res;
  }

  bool is_leaf( node_index_t const& n ) const
  {
    return ntk_.value( n ) == ( 3u | ( ntk_.trav_id() << 2u ) );
  }

  void make_leaf( node_index_t const& n ) const
  {
    ntk_.set_value( n, ( 3u | ( ntk_.trav_id() << 2u ) ) );
  }

  bool is_mffc( node_index_t const& n ) const
  {
    return ntk_.value( n ) == ( 1u | ( ntk_.trav_id() << 2u ) );
  }

  void make_mffc( node_index_t const& n ) const
  {
    ntk_.set_value( n, ( 1u | ( ntk_.trav_id() << 2u ) ) );
  }

  bool is_divisor( node_index_t const& n ) const
  {
    return ntk_.value( n ) == ( 2u | ( ntk_.trav_id() << 2u ) );
  }

  void make_divisor( node_index_t const& n ) const
  {
    ntk_.set_value( n, ( 2u | ( ntk_.trav_id() << 2u ) ) );
  }

  std::vector<node_index_t> const& get_leaves() const
  {
    return window_.leaves;
  }

  std::vector<node_index_t> const& get_divs() const
  {
    return window_.divs;
  }

  std::vector<node_index_t> const& get_mffc() const
  {
    return window_.mffc;
  }

  node_index_t const& get_root() const
  {
    return window_.root;
  }

#pragma region Mffc
private:
  bool collect_mffc( node_index_t const& n )
  {
    for ( const auto& l : window_.leaves )
      ntk_.incr_fanout_size( l );
    /* dereference the node */
    window_.mffc.emplace_back( n );
    make_mffc( n );
    node_deref_rec( n );
    std::reverse( window_.mffc.begin(), window_.mffc.end() );
    for ( auto m : window_.mffc )
    {
      if ( ntk_.is_pi( m ) )
        return false;
    }
    /* reference it back */
    node_ref_rec( n );
    for ( const auto& l : window_.leaves )
      ntk_.decr_fanout_size( l );

    return true;
  }

  /* ! \brief Dereference the node's MFFC */
  void node_deref_rec( node_index_t const& n )
  {
    if ( ntk_.is_pi( n ) )
      return;

    ntk_.foreach_fanin( n, [&]( const auto& f ) {
      auto const& p = ntk_.get_node( f );

      ntk_.decr_fanout_size( p );
      if ( ntk_.fanout_size( p ) == 0 )
      {
        window_.mffc.emplace_back( p );
        make_mffc( p );
        node_deref_rec( p );
      }
    } );
  }

  /* ! \brief Reference the node's MFFC */
  void node_ref_rec( node_index_t const& n )
  {
    if ( ntk_.is_pi( n ) )
      return;

    ntk_.foreach_fanin( n, [&]( const auto& f ) {
      auto const& p = ntk_.get_node( f );

      auto v = ntk_.fanout_size( p );
      ntk_.incr_fanout_size( p );
      if ( v == 0 )
      {
        node_ref_rec( p );
      }
    } );
  }
#pragma endregion

#pragma region Divisors
  bool collect_divs_tfi( node_index_t const& n )
  {
    if ( is_leaf( n ) )
    {
      return true;
    }
    if ( is_divisor( n ) )
    {
      return true;
    }
    /* PI and not leaf means error */
    if ( ntk_.is_pi( n ) )
    {
      std::cout << "pi not leaf" << std::endl;
      return false;
    }
    bool is_good = true;
    ntk_.foreach_fanin( n, [&]( auto fi ) {
      node_index_t const ni = ntk_.get_node( fi );
      is_good &= collect_divs_tfi( ni );
    } );
    if ( !is_mffc( n ) && !is_leaf( n ) )
    {
      window_.divs.emplace_back( n );
      make_divisor( n );
    }
    return is_good;
  }

  bool collect_divs_sides( node_index_t const& n )
  {
    bool quit = false;
    bool is_good = true;
    for ( auto i = 0u; i < window_.divs.size(); ++i )
    {
      auto const d = window_.divs.at( i );

      if ( ntk_.fanout_size( d ) > ps_.skip_fanout_limit_for_divisors )
      {
        continue;
      }
      if ( window_.divs.size() >= ps_.max_num_divisors )
      {
        break;
      }

      /* if the fanout has all fanins in the set, add it */
      bool add = true;
      ntk_.foreach_fanout( d, [&]( node_index_t const& p ) {
        if ( is_leaf( p ) || is_divisor( p ) || is_mffc( p ) )
        {
          return true; /* next fanout */
        }
        if ( ntk_.fanout_size( p ) <= 0 )
        {
          return true; /* next fanout */
        }

        bool all_fanins_visited = true;
        ntk_.foreach_fanin( p, [&]( const auto& g ) {
          if ( !is_leaf( ntk_.get_node( g ) ) && !is_divisor( ntk_.get_node( g ) ) )
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
        ntk_.foreach_fanin( p, [&]( const auto& g ) {
          if ( ntk_.get_node( g ) == window_.root )
          {
            has_root_as_child = true;
            return false; /* terminate fanin-loop */
          }
          return true; /* next fanin */
        } );

        if ( has_root_as_child )
        {
          std::cout << "has root as child " << window_.root << std::endl;
          return true; /* next fanout */
        }
        is_good &= all_fanins_visited;
        window_.divs.emplace_back( p );
        make_divisor( p );

        /* quit computing divisors if there are too many of them */
        if ( window_.divs.size() >= ps_.max_num_divisors )
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
    for ( uint32_t i = 0; i < window_.leaves.size(); i++ )
    {
      fn( window_.leaves[i], i );
    }
  }

  /*! \brief iterator to apply a lambda to all the nodes */
  template<typename Fn>
  void foreach_divisor( Fn&& fn ) const
  {
    for ( uint32_t i = 0; i < window_.divs.size(); i++ )
    {
      fn( window_.divs[i], i );
    }
  }

  /*! \brief iterator to apply a lambda to all the nodes */
  template<typename Fn>
  void foreach_mffc( Fn&& fn ) const
  {
    for ( uint32_t i = 0; i < window_.mffc.size(); i++ )
    {
      fn( window_.mffc[i], i );
    }
  }

  void print() const
  {
    std::cout << "leaves" << std::endl;
    foreach_leaf( [&]( auto l, auto i ) {
      std::cout << l << " ";
    } );
    std::cout << std::endl;
    std::cout << "divs" << std::endl;
    foreach_divisor( [&]( auto d, auto i ) {
      std::cout << d << " ";
    } );
    std::cout << std::endl;
    std::cout << "mffc" << std::endl;
    foreach_mffc( [&]( auto m, auto i ) {
      std::cout << m << " ";
    } );
    std::cout << std::endl;
  }

  bool is_valid() const
  {
    return window_.valid;
  }

private:
  Ntk& ntk_;
  window_manager_params const& ps_;
  window_manager_stats<Ntk>& st_;
  cut_comp cuts_;
  window_t window_;
};

} // namespace mockturtle