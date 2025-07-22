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
  \file spfd_manager.hpp
  \brief Classes to perform functional analysis using SPFDs

  \author Andrea Costamagna
*/

#pragma once

#include "../../../utils/constexpr_functions.hpp"
#include <kitty/simd_operations.hpp>

namespace mockturtle
{
template<typename Tt_t, uint32_t MaxNumMasks>
class spfd
{
public:
  spfd( Tt_t const& func, Tt_t const& care )
  {
    init( func, care );
  }

  void init( Tt_t const& func, Tt_t const& care )
  {
    care_ = care;
    func_[0] = kitty::simd::binary_and( care, func );
    func_[0] = kitty::simd::binary_and( care, kitty::simd::unary_not( func ) );
    masks_[0] = care;
    num_masks_ = 1;
    num_edges_ = count_edges( 0 );

    bool const is_killed = num_edges_ <= 0;

    kills_[0] = is_killed;
    num_kills_ = static_cast<uint32_t>( is_killed );
  }

  void reset()
  {
    masks_[0] = care_;
    num_masks_ = 1;
    num_edges_ = count_edges( 0 );

    bool const is_killed = num_edges_ <= 0;

    kills_[0] = is_killed;
    num_kills_ = static_cast<uint32_t>( is_killed );
  }

  bool update( Tt_t const& tt )
  {
    num_edges_ = 0;
    for ( uint32_t iMask{ 0 }; iMask < num_masks_; ++iMask )
    {
      if ( kills_[iMask] )
      {
        kills_[num_masks_ + iMask] = true;
        num_kills_++;
      }
      else
      {
        masks_[num_masks_ + iMask] = kitty::simd::binary_and( masks_[iMask], tt );
        masks_[iMask] &= kitty::simd::unary_not( tt );

        if ( is_killed( num_masks_ + iMask ) )
        {
          kills_[num_masks_ + iMask] = true;
          num_kills_++;
        }
        else
        {
          kills_[num_masks_ + iMask] = false;
          num_edges_ += count_edges( num_masks_ + iMask );
        }

        if ( is_killed( iMask ) )
        {
          kills_[iMask] = true;
          num_kills_++;
        }
        else
        {
          kills_[iMask] = false;
          num_edges_ += count_edges( iMask );
        }
      }
    }
    num_masks_ = num_masks_ << 1u;
    return true;
  }

  uint32_t evaluate( Tt_t const& tt ) const
  {
    uint32_t res = 0;
    for ( auto iMask{ 0 }; iMask < num_masks_; ++iMask )
    {
      if ( !kills_[iMask] )
      {
        Tt_t const mask1 = kitty::simd::binary_and( masks_[iMask], tt );
        Tt_t const mask0 = kitty::simd::binary_and( masks_[iMask], kitty::simd::unary_not( tt ) );
        res += kitty::count_ones( kitty::simd::binary_and( func_[1], mask0 ) ) * kitty::count_ones( kitty::simd::binary_and( func_[0], mask0 ) );
        res += kitty::count_ones( kitty::simd::binary_and( func_[1], mask1 ) ) * kitty::count_ones( kitty::simd::binary_and( func_[0], mask1 ) );
      }
    }
    return res;
  }

  bool is_covered() const
  {
    return num_masks_ <= num_kills_;
  }

  bool is_saturated() const
  {
    return num_masks_ >= MaxNumMasks;
  }

  uint32_t get_num_edges() const
  {
    return num_edges_;
  }

private:
  bool is_killed( uint32_t i ) const
  {
    bool const only_onset = ( kitty::count_ones( kitty::simd::binary_and( masks_[i], func_[0] ) ) == 0 );
    bool const only_ofset = ( kitty::count_ones( kitty::simd::binary_and( masks_[i], func_[1] ) ) == 0 );
    return only_onset || only_ofset;
  }

  uint32_t count_edges( uint32_t i ) const
  {
    Tt_t const onset = kitty::simd::binary_and( masks_[i], func_[0] );
    Tt_t const ofset = kitty::simd::binary_and( masks_[i], func_[1] );
    return kitty::count_ones( onset ) * kitty::count_ones( ofset );
  }

  uint32_t count_edges() const
  {
    uint32_t res{ 0 };
    for ( auto i = 0u; i < num_masks_; ++i )
    {
      res += count_edges( i );
    }
    return res;
  }

private:
  std::array<Tt_t, MaxNumMasks> masks_;
  std::array<bool, MaxNumMasks> kills_;
  uint32_t num_masks_;
  uint32_t num_edges_;
  uint32_t num_kills_;
  Tt_t care_;
  std::array<Tt_t, 2u> func_;
};

template<typename Tt_t, uint32_t MaxNumMasks>
class spfd_manager
{
public:
  spfd_manager()
  {
    Tt_t tmp;
    if constexpr ( log2_ceil( MaxNumMasks ) > 6u )
    {
      kitty::simd::test_avx2_advantage( tmp, tmp.num_vars() );
    }
  }

  void reset()
  {
    for ( auto& spfd : spfds_ )
      spfd.reset();
  }

  void init( std::vector<Tt_t const*> const& targets, Tt_t const& care )
  {
    spfds_.reserve( targets.size() );
    for ( Tt_t const* ptt : targets )
      spfds_.emplace_back( *ptt, care );
  }

  void update( Tt_t const& tt )
  {
    for ( auto& spfd : spfds_ )
      spfd.update( tt );
  }

  bool is_covered()
  {
    for ( auto const& spfd : spfds_ )
      if ( !spfd.is_covered() )
        return false;
    return true;
  }

  bool is_saturated()
  {
    for ( auto const& spfd : spfds_ )
      if ( spfd.is_saturated() )
        return true;
    return false;
  }

  uint32_t get_num_edges() const
  {
    uint32_t res = 0;
    for ( auto const& spfd : spfds_ )
      res += spfd.get_num_edges();
    return res;
  }

  uint32_t evaluate( Tt_t const& tt ) const
  {
    uint32_t res = 0;
    for ( auto const& spfd : spfds_ )
      res += spfd.evaluate( tt );
    return res;
  }

private:
  std::vector<spfd<Tt_t, MaxNumMasks>> spfds_;
};

} // namespace mockturtle