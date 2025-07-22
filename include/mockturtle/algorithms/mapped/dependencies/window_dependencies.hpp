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
  \file window_dependencies.hpp
  \brief Compute dependencies allowing to window the pivot node.

  \author Andrea Costamagna
*/

#pragma once

#include "../../../utils/constexpr_functions.hpp"
#include "../boolean/spfd_manager.hpp"
#include "dependency_cut.hpp"

namespace mockturtle
{

template<class Ntk, uint32_t CubeSizeLeaves = 6u, uint32_t MaxNumVars = 6u, uint32_t MaxCubeSize = 12u>
class window_dependencies
{

public:
  using signal_t = typename Ntk::signal;
  using node_index_t = typename Ntk::node;
  static constexpr uint32_t NumBits = 1u << CubeSizeLeaves;
  static constexpr uint32_t ExaNumPairs = log2_ceil( NumBits * ( NumBits - 1 ) / 2 );
  static constexpr uint32_t NumPairs = std::min( MaxCubeSize, ExaNumPairs );
  using information_t = kitty::static_truth_table<NumPairs>;
  using signature_t = kitty::static_truth_table<CubeSizeLeaves>;

public:
  window_dependencies( Ntk& ntk )
      : ntk_( ntk )
  {
    if constexpr ( CubeSizeLeaves > 6u )
    {
      signature_t tmp;
      kitty::simd::test_avx2_advantage( tmp, CubeSizeLeaves );
    }
  }

  template<typename WinMng, typename WinSim>
  void run( WinMng const& window, WinSim& simulator )
  {
    assert( ( std::is_same<signature_t, typename WinSim::signature_t>::value && "signatures have different type" ) );
    cuts_.clear();

    // initialize the approximate information manager
    load_information( window, simulator );
    // identify candidates through branch and bound
    identify_candidates( window, simulator );
    // check if the candidates are valid cuts. If not, complete them.
    exactify_candidates( window, simulator );
  }

  template<typename Fn>
  void foreach_cut( Fn&& fn )
  {
    for ( auto i = 0u; i < cuts_.size(); ++i )
    {
      fn( cuts_[i], i );
    }
  }

#pragma region Information Loading
private:
  void load_information( information_t& info, signature_t const& sign, signature_t const& care )
  {
    uint32_t const num_bits_sign = sign.num_bits();
    uint32_t const num_bits_info = info.num_bits();

    information_t sign_info_pos;
    information_t sign_info_neg;
    auto const sign_pos = sign;
    auto const sign_neg = kitty::simd::unary_not( sign );

    for ( auto block = 0u; block < sign.num_blocks(); ++block )
    {
      *( sign_info_pos.begin() + block ) = *( sign_pos.begin() + block );
      *( sign_info_neg.begin() + block ) = *( sign_neg.begin() + block );
    }

    uint32_t bit_info = 0, bit_sign = 0;
    while ( ( bit_sign < num_bits_sign ) && ( ( bit_info + num_bits_sign - bit_sign ) < num_bits_info ) )
    {
      if ( kitty::get_bit( care, bit_sign ) )
      {
        auto const tmp = kitty::get_bit( sign, bit_sign ) ? kitty::shift_right( sign_info_neg, bit_sign + 1 ) : kitty::shift_right( sign_info_pos, bit_sign + 1 );
        info = kitty::simd::binary_or( info, kitty::shift_left( tmp, bit_info ) );
        bit_info += num_bits_sign - bit_sign - 1;
      }
      bit_sign += 1;
    }
  }

  template<typename WinMng, typename WinSim>
  void load_information( WinMng const& window, WinSim& simulator )
  {
    assert( ( std::is_same<signature_t, typename WinSim::signature_t>::value && "signatures have different type" ) );
    signature_t const& care = simulator.get_careset();

    auto const n = window.get_pivot();

    divs_info_.clear();
    divs_info_.reserve( window.num_divisors() );
    window.foreach_divisor( [&]( auto const& f, auto i ) {
      divs_info_.emplace_back();
      information_t& info = divs_info_.back();
      auto const& sign = simulator.get( f );

      load_information( info, sign, care );
    } );

    root_info_.clear();
    root_info_.reserve( ntk_.num_outputs( n ) );
    ntk_.foreach_output( n, [&]( auto const& f ) {
      root_info_.emplace_back();
      information_t& info = root_info_.back();
      auto const& sign = simulator.get( f );
      load_information( info, sign, care );
    } );

    info_from_.resize( divs_info_.size() );
    info_from_.back() = divs_info_.back();
    for ( int i = info_from_.size() - 2; i >= 0; --i )
    {
      info_from_[i] = info_from_[i + 1] | divs_info_[i];
    }

    certain_from_.resize( divs_info_.size() );
    int i = divs_info_.size() - 1;
    while ( i >= 0 )
    {
      if ( ( i < ( divs_info_.size() - 1 ) ) && certain_from_[i + 1] )
        certain_from_[i] = true;
      else
      {
        certain_from_[i] = true;
        for ( auto const& info : root_info_ )
        {
          if ( !kitty::equal( info_from_[i] & info, info ) )
            certain_from_[i] = false;
        }
      }
      i--;
    }
  }
#pragma endregion

#pragma region Candidates Identification

  void update_information( information_t const& remove, std::vector<information_t>& todos )
  {
    for ( auto& todo : todos )
    {
      todo = kitty::simd::binary_and( todo, kitty::simd::unary_not( remove ) );
    }
  }

  bool is_done( std::vector<information_t>& todos )
  {
    for ( auto const& todo : todos )
    {
      if ( !kitty::is_const0( todo ) )
      {
        return false;
      }
    }
    return true;
  }

  bool is_possible_from( uint32_t index, std::vector<information_t> const& todos )
  {
    if ( certain_from_[index] )
    {
      return true;
    }
    for ( auto const& todo : todos )
    {
      if ( !kitty::equal( kitty::simd::binary_and( info_from_[index], todo ), todo ) )
      {
        return false;
      }
    }
    return true;
  }

  bool contains( std::vector<uint32_t> const& cut1, std::vector<uint32_t> const& cut2 )
  {
    if ( cut2.size() > cut1.size() )
      return false;
    bool contains = true;
    auto i = 0;
    auto j = 0;
    while ( contains && ( i < cut1.size() ) && ( j < cut2.size() ) )
    {
      if ( cut1[i] < cut2[j] )
        i++;
      else if ( cut1[i] == cut2[j] )
      {
        i++;
        j++;
      }
      else if ( cut1[i] > cut2[j] )
        contains = false;
    }
    if ( contains && ( j == cut2.size() ) )
      return true;
    return false;
  }

  bool contains_previous( std::vector<uint32_t> const& cut, std::vector<std::vector<uint32_t>> const& cuts )
  {
    for ( auto const& other : cuts )
    {
      if ( contains( cut, other ) )
        return true;
    }
    return false;
  }

  void add_to_cuts( std::vector<uint32_t> const& cut, std::vector<std::vector<uint32_t>>& cuts )
  {
    std::vector<uint32_t> erase;
    for ( auto i = 0u; i < cuts.size(); ++i )
    {
      auto const other = cuts[i];
      if ( contains( other, cut ) )
        erase.push_back( i );
    }
    for ( auto it = erase.rbegin(); it < erase.rend(); ++it )
    {
      cuts.erase( cuts.begin() + *it );
    }
    cuts.push_back( cut );
  }

  void identify_candidates_recursive( std::vector<std::vector<uint32_t>>& cuts,
                                      std::vector<uint32_t> cut,
                                      uint32_t begin,
                                      std::vector<information_t> todos )
  {
    if ( contains_previous( cut, cuts ) )
      return;
    if ( is_done( todos ) )
    {
      add_to_cuts( cut, cuts );
      return;
    }
    if ( ( begin >= divs_info_.size() ) )
      return;
    if ( cut.size() >= MaxNumVars )
      return;

    if ( ( begin < ( divs_info_.size() - 1 ) ) && is_possible_from( begin + 1, todos ) )
      identify_candidates_recursive( cuts, cut, begin + 1, todos );

    if ( is_possible_from( begin, todos ) )
    {
      cut.push_back( begin );
      update_information( divs_info_[begin], todos );
      identify_candidates_recursive( cuts, cut, begin + 1, todos );
    }
  }

  template<typename WinMng, typename WinSim>
  void identify_candidates( WinMng const& window, WinSim& simulator )
  {
    auto todos = root_info_;
    uint32_t const begin = 0;
    std::vector<uint32_t> cut;
    std::vector<std::vector<uint32_t>> cuts;
    identify_candidates_recursive( cuts, cut, begin, todos );

    auto const n = window.get_pivot();
    cuts_.clear();
    for ( auto const& cut : cuts )
    {
      std::vector<signal_t> leaves;
      for ( auto const& index : cut )
      {
        leaves.push_back( window.get_divisor( index ) );
      }
      cuts_.emplace_back( dependency_t::WINDOW_DEP, n, leaves );
    }
  }

#pragma endregion

#pragma region Exactify Candidates
  template<typename WinMng, typename WinSim>
  void exactify_candidates( WinMng const& window, WinSim& simulator )
  {
    assert( ( std::is_same<signature_t, typename WinSim::signature_t>::value && "signatures have different type" ) );
    signature_t const& care = simulator.get_careset();

    static constexpr uint32_t capacity = ( 1u << MaxNumVars );

    auto const n = window.get_pivot();
    std::vector<signature_t const*> funcs_ptrs;
    ntk_.foreach_output( n, [&]( auto const& f ) {
      funcs_ptrs.push_back( &simulator.get( f ) );
    } );

    spfds_.init( funcs_ptrs, care );

    std::vector<uint32_t> erase;
    for ( auto i = 0u; i < cuts_.size(); ++i )
    {
      auto& cut = cuts_[i];
      bool const add = exactify_candidates_greedy( cut, window, simulator );
      if ( add )
      {
        std::vector<signature_t const*> in_ptrs;
        for ( auto const& l : cut )
          in_ptrs.push_back( &simulator.get( l ) );
        ntk_.foreach_output( n, [&]( auto const& f ) {
          auto const func = extract_function<signature_t, MaxNumVars>( in_ptrs, simulator.get( f ), care );
          cut.add_func( func );
        } );
      }
      else
      {
        erase.push_back( i );
      }
    }
  }

  template<typename WinMng, typename WinSim>
  bool exactify_candidates_greedy( dependency_cut_t<Ntk, MaxNumVars>& cut, WinMng const& window, WinSim& simulator )
  {
    spfds_.reset();
    for ( auto const& l : cut )
      spfds_.update( simulator.get( l ) );

    auto cnt = cut.size();
    uint32_t best_num_edges = spfds_.get_num_edges();
    while ( !spfds_.is_covered() && !spfds_.is_saturated() && ( cnt < MaxNumVars ) )
    {
      std::optional<uint32_t> best_div;
      window.foreach_divisor( [&]( auto const& d, auto i ) {
        auto const& sign = simulator.get( d );
        auto const num_edges = spfds_.evaluate( sign );
        if ( num_edges < best_num_edges )
        {
          best_div = std::make_optional( d );
          best_num_edges = num_edges;
        }
      } );
      if ( best_div )
        cut.add_leaf( *best_div );
      else
        return false;
      cnt++;
    }

    return spfds_.is_covered();
  }
#pragma endregion

private:
  Ntk& ntk_;
  std::vector<dependency_cut_t<Ntk, MaxNumVars>> cuts_;
  std::vector<information_t> divs_info_;
  std::vector<information_t> info_from_;
  std::vector<bool> certain_from_;
  std::vector<information_t> root_info_;
  spfd_manager<signature_t, ( 1u << MaxNumVars )> spfds_;
};

} // namespace mockturtle