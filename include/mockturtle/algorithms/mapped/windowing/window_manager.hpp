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
  \brief Construction of windows for mapped networks.

  \author Andrea Costamagna
*/

#pragma once

namespace mockturtle
{

struct window_manager_params
{
  bool preserve_depth = true;
  uint32_t odc_levels = 0;
  uint32_t cut_limit = 8;
  uint32_t skip_fanout_limit_for_divisors{ 100 };
  uint32_t max_num_divisors{ 128 };
};

struct window_manager_stats
{
};

template<class Ntk>
class window_manager
{
public:
  using node_index_t = typename Ntk::node;
  using signal_t = typename Ntk::signal;

  struct window_t
  {
    node_index_t pivot;
    std::vector<node_index_t> tfos;
    std::vector<node_index_t> mffc;
    std::vector<node_index_t> divs;
    std::vector<node_index_t> outputs;
    std::vector<node_index_t> inputs;
  };

public:
  window_manager( Ntk& ntk, window_manager_params const& ps, window_manager_stats& st )
      : ntk_( ntk ),
        ps_( ps ),
        st_( st )
  {
  }

  [[nodiscard]] bool run( node_index_t const& n )
  {
    init( n );

    // label the MFFC nodes
    collect_mffc_nodes();
    for ( auto const& m : window_.mffc )
    {
      ntk_.set_visited( m, ntk_.trav_id() );
      make_alien( m );
    }
    window_.mffc.clear();

    window_.inputs = { window_.pivot };
    // expand the leaves until reaching the boundary of the MFFC
    expand_leaves( [&]( node_index_t const& v ) { return ( ntk_.visited( v ) == ntk_.trav_id() ) && !ntk_.is_pi( v ); },
                   [&]( node_index_t const& v ) { make_mffc( v ); } );

    collect_mffc_nodes();

    /* expand toward the tfo */
    if ( ps_.odc_levels > 0 )
    {
      collect_nodes_tfo();
      // expand the divisors set to try removing upper-leaves with reconvergence
      collect_side_divisors();
    }

    // expand the leaves to find reconvergences
    expand_leaves( [&]( node_index_t const& n ) { return !ntk_.is_pi( n ); },
                   [&]( node_index_t const& v ) { make_divisor( v ); } );

    // expand the divisors set to try removing upper-leaves with reconvergence
    collect_side_divisors();

    if ( ps_.odc_levels > 0 )
    {
      topological_sort( window_.tfos );
      topological_sort( window_.outputs );
    }
    topological_sort( window_.mffc );
    topological_sort( window_.divs );
    topological_sort( window_.inputs );

    return true;
  }

  std::vector<node_index_t> const& get_divisors() const
  {
    std::vector<node_index_t> const& res = window_.divs;
    return res;
  }

  std::vector<node_index_t> const& get_tfos() const
  {
    return window_.tfos;
  }

  std::vector<node_index_t> const& get_outputs() const
  {
    return window_.outputs;
  }

  std::vector<node_index_t> const& get_leaves() const
  {
    return window_.inputs;
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
    return window_.pivot;
  }

#pragma region Miscellanea
private:
  void init( node_index_t const& n )
  {
    ntk_.incr_trav_id();
    window_.pivot = n;
    make_mffc( n );
    window_.outputs.clear();
    window_.tfos.clear();
    window_.mffc.clear();
    window_.divs.clear();
    window_.inputs.clear();
  }

  void topological_sort( std::vector<node_index_t>& nodes )
  {
    std::stable_sort( nodes.begin(), nodes.end(), [&]( auto const& a, auto const& b ) {
      return ntk_.level( a ) < ntk_.level( b );
    } );
  }

  bool is_tfo( node_index_t const& n ) const
  {
    return ntk_.value( n ) == ( 4u | ( ntk_.trav_id() << 3u ) );
  }

  void make_tfo( node_index_t const& n ) const
  {
    ntk_.set_value( n, ( 4u | ( ntk_.trav_id() << 3u ) ) );
  }

  bool is_contained( node_index_t const& n ) const
  {
    return ( ( ntk_.value( n ) >> 3u ) == ntk_.trav_id() );
  }

  bool is_output( node_index_t const& n ) const
  {
    return ntk_.value( n ) == ( 5u | ( ntk_.trav_id() << 3u ) );
  }

  void make_output( node_index_t const& n ) const
  {
    ntk_.set_value( n, ( 5u | ( ntk_.trav_id() << 3u ) ) );
  }

  void make_alien( node_index_t const& n ) const
  {
    ntk_.set_value( n, 0u );
  }

  bool is_leaf( node_index_t const& n ) const
  {
    return ntk_.value( n ) == ( 3u | ( ntk_.trav_id() << 3u ) );
  }

  void make_leaf( node_index_t const& n ) const
  {
    ntk_.set_value( n, ( 3u | ( ntk_.trav_id() << 3u ) ) );
  }

  bool is_mffc( node_index_t const& n ) const
  {
    return ntk_.value( n ) == ( 1u | ( ntk_.trav_id() << 3u ) );
  }

  void make_mffc( node_index_t const& n ) const
  {
    ntk_.set_value( n, ( 1u | ( ntk_.trav_id() << 3u ) ) );
  }

  bool is_divisor( node_index_t const& n ) const
  {
    return ntk_.value( n ) == ( 2u | ( ntk_.trav_id() << 3u ) );
  }

  void make_divisor( node_index_t const& n ) const
  {
    ntk_.set_value( n, ( 2u | ( ntk_.trav_id() << 3u ) ) );
  }
#pragma endregion

#pragma region TFO
  void collect_nodes_tfo()
  {
    std::vector<node_index_t> outputs;
    std::vector<node_index_t> inputs;
    std::vector<node_index_t> tfos;
    window_.outputs = { window_.pivot };
    window_.tfos = {};

    uint32_t num_leaves = window_.inputs.size();

    for ( auto lev = 0; lev < ps_.odc_levels; ++lev )
    {
      outputs.clear();
      inputs.clear();
      tfos.clear();
      for ( auto const& n : window_.outputs )
      {
        int cnt = 0;
        ntk_.foreach_fanout( n, [&]( auto no ) {
          if ( !is_tfo( no ) && !is_output( no ) )
          {
            outputs.push_back( no );
            make_output( no );
          }
          cnt++;
        } );
        if ( cnt == 0 )
        {
          make_output( n );
          outputs.push_back( n );
        }
        else if ( n != window_.pivot )
        {
          make_tfo( n );
          tfos.push_back( n );
        }
      }
      for ( node_index_t const& n : outputs )
      {
        ntk_.foreach_fanin( n, [&]( auto const& fi, auto ii ) {
          auto const ni = ntk_.get_node( fi );
          if ( !is_contained( ni ) )
          {
            inputs.push_back( ni );
            make_leaf( ni );
            num_leaves += 1;
          }
        } );
      }
      if ( num_leaves <= ps_.cut_limit )
      {
        window_.outputs = outputs;
        for ( auto const& n : inputs )
          window_.inputs.push_back( n );
        for ( auto const& n : tfos )
          window_.tfos.push_back( n );
      }
      else
      {
        for ( auto const& n : inputs )
          make_alien( n );
        for ( auto const& n : window_.outputs )
          make_output( n );
      }
    }
  }
#pragma endregion

#pragma region Leaves
  template<typename DoExpandFn, typename Apply>
  void expand_leaves( DoExpandFn&& do_expand, Apply&& apply )
  {
    while ( true )
    {
      int best_cost = ps_.cut_limit - window_.inputs.size() + 1;
      std::optional<node_index_t> best_leaf;
      for ( auto const& l : window_.inputs )
      {
        if ( do_expand( l ) )
        {
          int cost = compute_leaf_cost( l );
          if ( cost < best_cost )
          {
            best_cost = cost;
            best_leaf = std::make_optional( l );
          }
        }
      }
      if ( best_leaf )
      {
        ntk_.foreach_fanin( *best_leaf, [&]( auto const& fi, auto ii ) {
          auto const ni = ntk_.get_node( fi );
          if ( !is_contained( ni ) )
          {
            window_.inputs.push_back( ni );
            make_leaf( ni );
          }
        } );
        auto& leaves = window_.inputs;
        apply( *best_leaf );
        leaves.erase( std::remove( leaves.begin(), leaves.end(), *best_leaf ), leaves.end() );
      }
      else
      {
        return;
      }
    }
  }

  int compute_leaf_cost( node_index_t const& n )
  {
    if ( ntk_.is_pi( n ) )
      return ps_.cut_limit;

    int cost = -1;
    ntk_.foreach_fanin( n, [&]( auto const& fi, auto ii ) {
      auto const ni = ntk_.get_node( fi );
      if ( !is_contained( ni ) )
        cost += 1;
    } );
    return cost;
  }
#pragma endregion Leaves

#pragma region Mffc
  void collect_mffc_nodes()
  {
    /* dereference the node_index_t */
    window_.mffc.clear();
    make_mffc( window_.pivot );
    window_.mffc = { window_.pivot };

    for ( auto const& l : window_.inputs )
      ntk_.incr_fanout_size( l );

    node_deref_rec( window_.pivot );
    node_ref_rec( window_.pivot );

    for ( auto const& l : window_.inputs )
      ntk_.decr_fanout_size( l );
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
        make_mffc( p );
        window_.mffc.push_back( p );
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
  void collect_side_divisors()
  {
    uint32_t max_level = std::numeric_limits<uint32_t>::min();
    for ( auto const& n : window_.outputs )
    {
      max_level = std::max( max_level, ntk_.level( n ) );
    }

    window_.divs = window_.inputs;

    bool done = false;
    while ( !done )
    {
      done = true;
      uint32_t num_divs = window_.divs.size();
      for ( auto i = 0u; i < num_divs; ++i )
      {
        auto const d = window_.divs[i];
        if ( ntk_.fanout_size( d ) > ps_.skip_fanout_limit_for_divisors )
          continue;
        if ( window_.divs.size() >= ps_.max_num_divisors )
        {
          done = true;
          break;
        }
        /* if the fanout has all fanins in the set, add it */
        bool add = true;
        ntk_.foreach_fanout( d, [&]( node_index_t const& n ) {
          if ( window_.divs.size() >= ps_.max_num_divisors )
          {
            done = true;
            return;
          }

          if ( ps_.preserve_depth && ( ntk_.level( n ) >= max_level ) )
            return;

          ntk_.foreach_fanin( n, [&]( const auto& f ) {
            if ( !is_contained( ntk_.get_node( f ) ) )
            {
              add = false;
              return;
            }
          } );

          if ( add && !is_tfo( n ) && !is_output( n ) )
          {
            if ( is_leaf( n ) )
            {
              auto& leaves = window_.inputs;
              leaves.erase( std::remove( leaves.begin(), leaves.end(), n ), leaves.end() );
            }
            window_.divs.push_back( n );
            make_divisor( n );
            done = false;
          }
        } );
      }
    }
  }

private:
  Ntk& ntk_;
  window_t window_;
  window_manager_params const& ps_;
  window_manager_stats& st_;
};

} // namespace mockturtle