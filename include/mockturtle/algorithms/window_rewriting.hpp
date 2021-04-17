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
  \file window_rewriting.hpp
  \brief Window rewriting

  \author Heinz Riener
*/

#include "../utils/index_list.hpp"
#include "../utils/stopwatch.hpp"
#include "../utils/window_utils.hpp"
#include "../views/topo_view.hpp"
#include "../views/window_view.hpp"
#include "../utils/debugging_utils.hpp"

#include <abcresub/abcresub2.hpp>
#include <fmt/format.h>
#include <stack>

#pragma once

namespace mockturtle
{

struct window_rewriting_params
{
  uint64_t cut_size{6};
  uint64_t num_levels{5};

  bool filter_cyclic_substitutions{false};
}; /* window_rewriting_params */

struct window_rewriting_stats
{
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{0};

  /*! \brief Time for constructing windows. */
  stopwatch<>::duration time_window{0};

  /*! \brief Time for optimizing windows. */
  stopwatch<>::duration time_optimize{0};

  /*! \brief Time for substituting. */
  stopwatch<>::duration time_substitute{0};

  /*! \brief Time for updating level information. */
  stopwatch<>::duration time_levels{0};

  /*! \brief Time for updating window outputs. */
  stopwatch<>::duration time_update_vector{0};

  /*! \brief Time for topological sorting. */
  stopwatch<>::duration time_topo_sort{0};

  /*! \brief Time for encoding index_list. */
  stopwatch<>::duration time_encode{0};

  /*! \brief Total number of calls to the resub. engine. */
  uint64_t num_substitutions{0};
  uint64_t num_restrashes{0};
  uint64_t num_windows{0};
  uint64_t gain{0};

  window_rewriting_stats operator+=( window_rewriting_stats const& other )
  {
    time_total += other.time_total;
    time_window += other.time_window;
    time_optimize += other.time_optimize;
    time_substitute += other.time_substitute;
    time_levels += other.time_levels;
    time_update_vector += other.time_update_vector;
    time_topo_sort += other.time_topo_sort;
    time_encode += other.time_encode;
    num_substitutions += other.num_substitutions;
    num_windows += other.num_windows;
    gain += other.gain;
    return *this;
  }

  void report() const
  {
    stopwatch<>::duration time_other =
      time_total - time_window - time_topo_sort - time_optimize - time_substitute - time_levels - time_update_vector;

    fmt::print( "===========================================================================\n" );
    fmt::print( "[i] Windowing =  {:7.2f} ({:5.2f}%) (#win = {})\n",
                to_seconds( time_window ), to_seconds( time_window ) / to_seconds( time_total ) * 100, num_windows );
    fmt::print( "[i] Top.sort =   {:7.2f} ({:5.2f}%)\n", to_seconds( time_topo_sort ), to_seconds( time_topo_sort ) / to_seconds( time_total ) * 100 );
    fmt::print( "[i] Enc.list =   {:7.2f} ({:5.2f}%)\n", to_seconds( time_encode ), to_seconds( time_encode ) / to_seconds( time_total ) * 100 );
    fmt::print( "[i] Optimize =   {:7.2f} ({:5.2f}%) (#resubs = {}, est. gain = {})\n",
                to_seconds( time_optimize ), to_seconds( time_optimize ) / to_seconds( time_total ) * 100, num_substitutions, gain );
    fmt::print( "[i] Substitute = {:7.2f} ({:5.2f}%) (#hash upd. = {})\n",
                to_seconds( time_substitute ),
                to_seconds( time_substitute ) / to_seconds( time_total ) * 100,
                num_restrashes );
    fmt::print( "[i] Upd.levels = {:7.2f} ({:5.2f}%)\n", to_seconds( time_levels ), to_seconds( time_levels ) / to_seconds( time_total ) * 100 );
    fmt::print( "[i] Upd.win =    {:7.2f} ({:5.2f}%)\n", to_seconds( time_update_vector ), to_seconds( time_update_vector ) / to_seconds( time_total ) * 100 );
    fmt::print( "[i] Other =      {:7.2f} ({:5.2f}%)\n", to_seconds( time_other ), to_seconds( time_other ) / to_seconds( time_total ) * 100 );
    fmt::print( "---------------------------------------------------------------------------\n" );
    fmt::print( "[i] TOTAL =      {:7.2f}\n", to_seconds( time_total ) );
    fmt::print( "===========================================================================\n" );
  }
}; /* window_rewriting_stats */

namespace detail
{

template<typename Ntk>
bool is_contained_in_tfi_recursive( Ntk const& ntk, typename Ntk::node const& node, typename Ntk::node const& n )
{
  if ( ntk.color( node ) == ntk.current_color() )
  {
    return false;
  }
  ntk.paint( node );

  if ( n == node )
  {
    return true;
  }

  bool found = false;
  ntk.foreach_fanin( node, [&]( typename Ntk::signal const& fi ){
    if ( is_contained_in_tfi_recursive( ntk, ntk.get_node( fi ), n ) )
    {
      found = true;
      return false;
    }
    return true;
  });

  return found;
}

} /* namespace detail */

template<typename Ntk>
bool is_contained_in_tfi( Ntk const& ntk, typename Ntk::node const& node, typename Ntk::node const& n )
{
  /* do not even build the TFI, but just search for the node */
  ntk.new_color();
  return detail::is_contained_in_tfi_recursive( ntk, node, n );
}

namespace detail
{

template<class Ntk>
class window_rewriting_impl
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  explicit window_rewriting_impl( Ntk& ntk, window_rewriting_params const& ps, window_rewriting_stats& st )
    : ntk( ntk )
    , ps( ps )
    , st( st )
  {
    auto const update_level_of_new_node = [&]( const auto& n ) {
      ntk.resize_levels();
      update_node_level( n );
    };

    auto const update_level_of_existing_node = [&]( node const& n, const auto& old_children ) {
      (void)old_children;
      ntk.resize_levels();
      update_node_level( n );
    };

    auto const update_level_of_deleted_node = [&]( node const& n ) {
      assert( ntk.fanout_size( n ) == 0u );
      ntk.set_level( n, -1 );
    };

    ntk._events->on_add.emplace_back( update_level_of_new_node );
    ntk._events->on_modified.emplace_back( update_level_of_existing_node );
    ntk._events->on_delete.emplace_back( update_level_of_deleted_node );
  }

  void run()
  {
    stopwatch t( st.time_total );

    create_window_impl windowing( ntk );
    uint32_t const size = ntk.size();
    for ( uint32_t n = 0u; n < std::min( size, ntk.size() ); ++n )
    {
      if ( ntk.is_constant( n ) || ntk.is_ci( n ) || ntk.is_dead( n ) )
      {
        continue;
      }

      if ( const auto w = call_with_stopwatch( st.time_window, [&]() { return windowing.run( n, ps.cut_size, ps.num_levels ); } ) )
      {
        ++st.num_windows;

        auto topo_win = call_with_stopwatch( st.time_topo_sort, ( [&](){
          window_view win( ntk, w->inputs, w->outputs, w->nodes );
          topo_view topo_win{win};
          return topo_win;
        }) );

        abc_index_list il;
        call_with_stopwatch( st.time_encode, [&]() {
          encode( il, topo_win );
        } );

        auto il_opt = optimize( il );
        if ( !il_opt )
        {
          continue;
        }

        std::vector<signal> signals;
        for ( auto const& i : w->inputs )
        {
          signals.push_back( ntk.make_signal( i ) );
        }

        std::vector<signal> outputs;
        topo_win.foreach_co( [&]( signal const& o ){
          outputs.push_back( o );
        });

        std::vector<signal> new_outputs;
        uint32_t counter{0};

        ++st.num_substitutions;
        insert( ntk, std::begin( signals ), std::end( signals ), *il_opt,
                [&]( signal const& _new )
                {
                  auto const _old = outputs.at( counter++ );
                  if ( _old == _new )
                  {
                    return true;
                  }

                  /* ensure that _old is not in the TFI of _new */
                  // assert( !is_contained_in_tfi( ntk, ntk.get_node( _new ), ntk.get_node( _old ) ) );
                  if ( ps.filter_cyclic_substitutions && is_contained_in_tfi( ntk, ntk.get_node( _new ), ntk.get_node( _old ) ) )
                  {
                    std::cout << "undo resubstitution " << ntk.get_node( _old ) << std::endl;
                    if ( ntk.fanout_size( ntk.get_node( _new ) ) == 0u )
                    {
                      ntk.take_out_node( ntk.get_node( _new ) );
                    }
                    return false;
                  }

                  auto const updates = substitute_node( ntk.get_node( _old ), topo_win.is_complemented( _old ) ? !_new : _new );
                  update_vector( outputs, updates );
                  return true;
                });

        /* recompute levels and depth */
        // call_with_stopwatch( st.time_levels, [&]() { ntk.update_levels(); } );

        /* ensure that no dead nodes are reachable */
        assert( count_reachable_dead_nodes( ntk ) == 0u );

        /* ensure that the network structure is still acyclic */
        assert( network_is_acylic( ntk ) );

        /* ensure that the levels and depth is correct */
        assert( check_network_levels( ntk ) );

        /* update internal data structures in windowing */
        windowing.resize( ntk.size() );
      }
    }

    /* ensure that no dead nodes are reachable */
    assert( count_reachable_dead_nodes( ntk ) == 0u );
  }

private:
  /* optimize an index_list and return the new list */
  std::optional<abc_index_list> optimize( abc_index_list const& il, bool verbose = false )
  {
    stopwatch t( st.time_optimize );

    int *raw = ABC_CALLOC( int, il.size() + 1u );
    uint64_t i = 0;
    for ( auto const& v : il.raw() )
    {
      raw[i++] = v;
    }
    raw[1] = 0; /* fix encoding */

    abcresub::Abc_ResubPrepareManager( 1 );
    int *new_raw = nullptr;
    int num_resubs = 0;
    uint64_t new_entries = abcresub::Abc_ResubComputeWindow( raw, ( il.size() / 2u ), 1000, -1, 0, 0, 0, 0, &new_raw, &num_resubs );
    abcresub::Abc_ResubPrepareManager( 0 );

    if ( verbose )
    {
      fmt::print( "Performed resub {} times.  Reduced {} nodes.\n",
                  num_resubs, new_entries > 0 ? ( ( il.size() / 2u ) - new_entries ) : 0 );
    }
    st.gain += new_entries > 0 ? ( ( il.size() / 2u ) - new_entries ) : 0;

    if ( raw )
    {
      ABC_FREE( raw );
    }

    if ( new_entries > 0 )
    {
      std::vector<uint32_t> values;
      for ( uint32_t i = 0; i < 2*new_entries; ++i )
      {
        values.push_back( new_raw[i] );
      }
      values[1u] = 1; /* fix encoding */
      if ( new_raw )
      {
        ABC_FREE( new_raw );
      }
      return abc_index_list( values, il.num_pis() );
    }
    else
    {
      assert( new_raw == nullptr );
      return std::nullopt;
    }
  }

  /* substitute the node with a signal and return all strashing updates */
  std::vector<std::pair<node, signal>> substitute_node( node const& old_node, signal const& new_signal )
  {
    stopwatch t( st.time_substitute );

    std::vector<std::pair<node, signal>> updates;
    std::stack<std::pair<node, signal>> to_substitute;
    to_substitute.push( {old_node, new_signal} );
    while ( !to_substitute.empty() )
    {
      const auto [_old, _new] = to_substitute.top();
      to_substitute.pop();

      auto const p = std::make_pair( _old, _new );
      if ( std::find( std::begin( updates ), std::end( updates ), p ) == std::end( updates ) )
      {
        updates.push_back( p );
      }

      const auto parents = ntk.fanout( _old );
      for ( auto n : parents )
      {
        if ( const auto repl = ntk.replace_in_node( n, _old, _new ); repl )
        {
          to_substitute.push( *repl );
          ++st.num_restrashes;
        }
      }

      /* check outputs */
      ntk.replace_in_outputs( _old, _new );

      /* reset fan-in of old node */
      ntk.take_out_node( _old );
    }

    return updates;
  }

  void update_vector( std::vector<signal>& vs, std::vector<std::pair<node, signal>> const& updates )
  {
    stopwatch t( st.time_update_vector );

    for ( auto it = std::begin( vs ); it != std::end( vs ); ++it )
    {
      for ( auto it2 = std::begin( updates ); it2 != std::end( updates ); ++it2 )
      {
        if ( ntk.get_node( *it ) == it2->first )
        {
          *it = ntk.is_complemented( *it ) ? !it2->second : it2->second;
        }
      }
    }
  }

  /* recursively update the node levels and the depth of the critical path */
  void update_node_level( node const& n )
  {
    uint32_t const curr_level = ntk.level( n );

    uint32_t max_level = 0;
    ntk.foreach_fanin( n, [&]( const auto& f ) {
      auto const p = ntk.get_node( f );
      auto const fanin_level = ntk.level( p );
      if ( fanin_level > max_level )
      {
        max_level = fanin_level;
      }
    } );
    ++max_level;

    if ( ntk.depth() < max_level )
    {
      ntk.set_depth( max_level );
    }

    if ( curr_level != max_level )
    {
      ntk.set_level( n, max_level );

      ntk.foreach_fanout( n, [&]( const auto& p ) {
        update_node_level( p );
      } );
    }
  }

private:
  Ntk& ntk;
  window_rewriting_params ps;
  window_rewriting_stats& st;
}; /* window_rewriting_impl */

} /* detail */

template<class Ntk>
void window_rewriting( Ntk& ntk, window_rewriting_params const& ps = {}, window_rewriting_stats* pst = nullptr )
{
  window_rewriting_stats st;
  detail::window_rewriting_impl<Ntk>( ntk, ps, st ).run();
  if ( pst )
  {
    *pst = st;
  }
}

} /* namespace mockturtle */