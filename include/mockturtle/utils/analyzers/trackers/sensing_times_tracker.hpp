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
 * \file sensing_times_tracker.hpp
 * \brief Compute sensing time information of a network and updates upon change
 *
 * \author Andrea Costamagna
 */

#pragma once

#include "../../network_exploration/tfo_manager.hpp"
#include "../../signal_map.hpp"
#include <limits>

namespace mockturtle
{

/*! \brief Engine to evaluate the sensing times of a network.
 *
 * This engine computes the sensing times of a network and keeps them up-to-date.
 * During construction it is possible to specify a vector of sensing times. If
 * not provided, the engine assumes 0 sensing time at all the PIs. At construction,
 * the sensing times are propagated in the network.
 *
 * The engine is equipped with a transitive fanout (TFO) manager to handle the
 * timing updates. Two events trigger the update:
 * - Node addition: The sensing time of the node is computed from the fanins
 * - Node modification: The sensing time of the TFO of the affected nodes is updated.
 *
 * \tparam Ntk the network type to be analyzed.
 *
 * \verbatim embed:rst

   Example

   .. code-block:: c++

      bound_network ntk( gates );
      auto const a = ntk.create_pi();
      auto const f1 = ntk.create_node( { a }, 0 );
      sensing_times_tracker sensing( ntk );
      auto const f2 = ntk.create_node( { f1 }, 0 );
   \endverbatim
 */
template<class Ntk>
class sensing_times_tracker
{
public:
  using node_index_t = typename Ntk::node;
  using signal_t = typename Ntk::signal;

public:
  sensing_times_tracker( Ntk& ntk )
      : ntk_( ntk ),
        times_( ntk ),
        tfo_( ntk )
  {
    init();
  }

  sensing_times_tracker( Ntk& ntk, std::vector<double> const& input_sensings )
      : ntk_( ntk ),
        times_( ntk ),
        input_( input_sensings ),
        tfo_( ntk )
  {
    init();
  }

  void init()
  {
    /* check if the network implements the needed functionalities */
    static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
    static_assert( has_incr_trav_id_v<Ntk>, "Ntk does not implement the incr_trav_id method" );
    static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin method" );
    static_assert( has_foreach_output_v<Ntk>, "Ntk does not implement the foreach_output method" );
    static_assert( has_foreach_fanout_v<Ntk>, "Ntk does not implement the foreach_fanout method" );
    static_assert( has_foreach_pi_v<Ntk>, "Ntk does not implement the foreach_pi method" );
    static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
    static_assert( has_set_value_v<Ntk>, "Ntk does not implement the set_visited method" );
    static_assert( has_value_v<Ntk>, "Ntk does not implement the set_visited method" );
    static_assert( has_set_visited_v<Ntk>, "Ntk does not implement the set_visited method" );
    static_assert( has_visited_v<Ntk>, "Ntk does not implement the set_visited method" );
    static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
    static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
    compute_sensing_times();

    add_event_ = ntk_.events().register_add_event( [&]( const node_index_t& n ) {
      times_.resize();
      tfo_.resize();
      compute_sensing_time( n );
    } );

    modified_event_ = ntk_.events().register_modified_event( [&]( const auto& n, auto old_children ) {
      /* The tfo of the new node is modified as it acquires new fanouts */
      update_sensing_times_tfo( n );

      /* The tfo of the old children is modified since they loose the old node */
      for ( auto const& f : old_children )
      {
        update_sensing_times_tfo( ntk_.get_node( f ) );
      }
    } );
  }

  ~sensing_times_tracker()
  {
    if ( add_event_ )
    {
      ntk_.events().release_add_event( add_event_ );
    }

    if ( modified_event_ )
    {
      ntk_.events().release_modified_event( modified_event_ );
    }
  }

#pragma region Interface methods
public:
  [[nodiscard]] double get_time( signal_t const f ) const
  {
    return times_[f];
  }
#pragma endregion

#pragma region Implementation details
private:
  bool is_marked_ready( node_index_t const& n )
  {
    return ntk_.value( n ) == ntk_.trav_id();
  };

  void make_ready( node_index_t const& n )
  {
    ntk_.set_value( n, ntk_.trav_id() );
  };

  void compute_sensing_times()
  {
    times_.reset();
    if ( input_.size() < ntk_.num_pis() )
    {
      input_ = std::vector<double>( ntk_.num_pis(), 0 );
    }

    ntk_.incr_trav_id();
    ntk_.foreach_pi( [&]( auto n, auto index ) {
      times_[ntk_.make_signal( n )] = input_[index];
      make_ready( n );
    } );

    ntk_.foreach_po( [&]( auto f ) {
      compute_sensing_times_tfi( f );
    } );
  }

  /*! \brief Compute the sensing times of the nodes in the TFI of a signal's node*/
  void compute_sensing_times_tfi( signal_t const& f )
  {
    node_index_t n = ntk_.get_node( f );
    if ( is_marked_ready( n ) || ntk_.is_pi( n ) )
    {
      return;
    }

    ntk_.foreach_fanin( n, [&]( auto fi, auto ii ) {
      compute_sensing_times_tfi( fi );
    } );
    compute_sensing_time( n );

    make_ready( n );
  }

  /*! \brief Efficient update of the sensing times in the TFO of a node.\
   */
  void update_sensing_times_tfo( node_index_t const& n )
  {
    /* mark all the nodes that must be considered */
    tfo_.init( n );

    std::vector<node_index_t> old_worklist{ n };
    std::vector<node_index_t> new_worklist;
    old_worklist.reserve( 100 );
    new_worklist.reserve( 100 );
    tfo_.mark_seen( n );

    bool progress{ true };
    /* iterate over all the nodes having some inputs ready */
    while ( progress && !old_worklist.empty() )
    {
      progress = false;
      for ( node_index_t const& u : old_worklist )
      {
        bool ready = true;
        ntk_.foreach_fanin( u, [&]( auto f, auto i ) {
          node_index_t ni = ntk_.get_node( f );
          ready &= ntk_.is_pi( ni ) || !tfo_.belongs_to_tfo( ni ) || tfo_.is_marked_ready( ni );
        } );
        tfo_.mark_seen( u );
        if ( ready )
        {
          progress = true;
          tfo_.mark_ready( u );
          ntk_.foreach_output( u, [&]( auto const& fu ) {
            double old_sensing = times_[fu];
            compute_sensing_time_at_pin( fu );
            if ( std::abs( times_[fu] - old_sensing ) > std::numeric_limits<double>::epsilon() )
            {
              ntk_.foreach_fanout( fu, [&]( node_index_t const& o ) {
                if ( !tfo_.is_marked_seen( o ) ) // only insert if not seen
                {
                  new_worklist.push_back( o );
                  tfo_.mark_seen( o );
                }
              } );
            }
          } );
        }
        else
        {
          new_worklist.push_back( u );
        }
      }
      old_worklist.swap( new_worklist );
      new_worklist.clear();
    }
    if ( !progress )
    {
      std::cerr << "[e] Infinite loop in sensing times analyzer" << std::endl;
    }
  }

  void compute_sensing_time_at_pin( signal_t const& f )
  {
    node_index_t n = ntk_.get_node( f );
    if ( ntk_.is_pi( n ) )
    {
      times_[f] = input_[ntk_.pi_index( n )];
    }
    else
    {
      auto time = std::numeric_limits<double>::max();
      ntk_.foreach_fanin( n, [&]( auto const& fi, auto const ii ) {
        time = std::min( time, times_[fi] + ntk_.get_min_pin_delay( f, ii ) );
      } );
      times_[f] = time;
    }
  }

  void compute_sensing_time( node_index_t const& n )
  {
    ntk_.foreach_output( n, [&]( auto const& f ) {
      compute_sensing_time_at_pin( f );
    } );
  }
#pragma endregion

private:
  Ntk& ntk_;
  incomplete_signal_map<double, Ntk> times_;
  tfo_manager<Ntk> tfo_;
  std::vector<double> input_;
  /* events */
  std::shared_ptr<typename network_events<Ntk>::add_event_type> add_event_;
  std::shared_ptr<typename network_events<Ntk>::modified_event_type> modified_event_;
};

} // namespace mockturtle