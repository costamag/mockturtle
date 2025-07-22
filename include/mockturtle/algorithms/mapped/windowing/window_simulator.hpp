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

#include "../../../utils/signal_map.hpp"

namespace mockturtle
{

template<class Ntk, uint32_t CubeSizeLeaves = 12>
class window_simulator
{
public:
  using signal_t = typename Ntk::signal;
  using node_index_t = typename Ntk::node;
  using signature_t = kitty::static_truth_table<CubeSizeLeaves>;
  static constexpr uint32_t num_bits = 1u << CubeSizeLeaves;

public:
  window_simulator( Ntk& ntk )
      : ntk_( ntk ),
        sig_to_sim_( ntk )
  {
    sims_.reserve( 1000u );
    init();
  }

  void run( window_manager<Ntk> const& window )
  {
    sims_.reserve( window.size() );
    sims_.resize( CubeSizeLeaves );
    sig_to_sim_.reset();

    assign_inputs( window );

    window.foreach_divisor( [&]( auto const& f, auto i ) {
      auto const n = ntk_.get_node( f );
      if ( f.output == 0 && !window.is_leaf( n ) )
      {
        compute( window, n );
      }
    } );

    window.foreach_mffc( [&]( auto const& n, auto i ) {
      compute( window, n );
    } );

    window.foreach_tfo( [&]( auto const& n, auto i ) {
      compute( window, n );
    } );

    window.foreach_output( [&]( auto const& f, auto i ) {
      if ( f.output == 0 )
        compute( window, ntk_.get_node( f ) );
    } );
  }

  signature_t const& get( signal_t const& f ) const
  {
    return sims_[sig_to_sim_[f]];
  }

  signature_t const compute_observability_careset( window_manager<Ntk> const& window )
  {
    signature_t care;
    auto n = window.get_pivot();
    auto const& outputs = window.get_outputs();
    if ( ntk_.num_outputs( n ) == outputs.size() )
    {
      bool no_odcs = true;
      for ( auto i = 0u; i < outputs.size(); ++i )
      {
        no_odcs &= ( ntk_.get_node( outputs[i] ) == window.get_pivot() );
      }
      if ( no_odcs )
      {
        return ~care;
      }
    }

    for ( uint32_t m = 1u; m < ( 1u << ntk_.num_outputs( n ) ); ++m )
    {
      int i = 0;
      ntk_.foreach_output( n, [&]( auto const& f ) {
        if ( ( m >> i ) & 0x1 > 0 )
          sims_[sig_to_sim_[f]] = ~sims_[sig_to_sim_[f]];
        i++;
      } );
      window.foreach_tfo( [&]( auto no, auto io ) {
        re_compute( window, no );
      } );

      window.foreach_output( [&]( auto f, auto io ) {
        if ( f.output == 0 )
        {
          auto const n = ntk_.get_node( f );
          std::vector<signature_t> old;
          ntk_.foreach_output( n, [&]( auto const& fo ) {
            old.push_back( sims_[sig_to_sim_[fo]] );
          } );

          re_compute( window, n );

          ntk_.foreach_output( n, [&]( auto const& fo ) {
            care |= ( old[io] ^ sims_[sig_to_sim_[fo]] );
          } );
        }
      } );

      ntk_.foreach_output( n, [&]( auto const& f ) {
        sims_[sig_to_sim_[f]] = ~sims_[sig_to_sim_[f]];
      } );
      window.foreach_tfo( [&]( auto no, auto io ) {
        re_compute( window, no );
      } );

      window.foreach_output( [&]( auto f, auto io ) {
        if ( f.output == 0 )
        {
          auto const n = ntk_.get_node( f );
          re_compute( window, n );
        }
      } );
    }

    return care;
  }

  void compute( window_manager<Ntk> const& window, node_index_t const& n )
  {
    if ( window.is_leaf( n ) )
      return;

    std::vector<signature_t const*> sim_ptrs;
    ntk_.foreach_fanin( n, [&]( auto const& fi, auto ii ) {
      sim_ptrs.push_back( &sims_[sig_to_sim_[fi]] );
    } );
    auto const tts = ntk_.compute( n, sim_ptrs );
    uint32_t io = 0;
    ntk_.foreach_output( n, [&]( auto const& fo ) {
      if ( !sig_to_sim_.has( fo ) )
      {
        sig_to_sim_[fo] = sims_.size();
        sims_.push_back( tts[io++] );
      }
    } );
    return;
  }

  void re_compute( window_manager<Ntk> const& window, node_index_t const& n )
  {
    if ( window.is_leaf( n ) )
      return;

    std::vector<signature_t const*> sim_ptrs;
    ntk_.foreach_fanin( n, [&]( auto const& fi, auto ii ) {
      sim_ptrs.push_back( &sims_[sig_to_sim_[fi]] );
    } );
    auto const tts = ntk_.compute( n, sim_ptrs );
    uint32_t io = 0;
    ntk_.foreach_output( n, [&]( auto const& fo ) {
      sims_[sig_to_sim_[fo]] = tts[io++];
    } );
    return;
  }

private:
  void init()
  {
    sims_.resize( CubeSizeLeaves );
    for ( auto i = 0u; i < CubeSizeLeaves; ++i )
      kitty::create_nth_var( sims_[i], i );
  }

  void assign_inputs( window_manager<Ntk> const& window )
  {
    sims_.resize( CubeSizeLeaves );
    window.foreach_input( [&]( auto const& f, auto i ) {
      sig_to_sim_[f] = i;
    } );
  }

private:
  Ntk& ntk_;
  std::vector<signature_t> sims_;
  incomplete_signal_map<uint32_t, Ntk> sig_to_sim_;
};

} // namespace mockturtle