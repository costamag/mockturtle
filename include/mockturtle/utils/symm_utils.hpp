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
  \file symm_utils.hpp
  \brief Implements utilities to handle functional symmetries

  \author Andrea Costamagna
*/

#pragma once

#include <kitty/kitty.hpp>

namespace mockturtle
{

struct symmetries_t
{
  uint64_t data = 0u;

  symmetries_t() = default;

  template<typename TT>
  symmetries_t( TT const& tt )
  {
    data = 0u;
    assert( tt.num_vars() <= 8 );

    for ( auto i = 0u; i < tt.num_vars(); ++i )
    {
      auto const tt1 = kitty::cofactor1( tt, i );
      auto const tt0 = kitty::cofactor0( tt, i );
      if ( kitty::equal( tt0, tt1 ) )
        continue;
      for ( auto j = i + 1; j < tt.num_vars(); ++j )
      {
        auto const tt1 = kitty::cofactor1( tt, j );
        auto const tt0 = kitty::cofactor0( tt, j );
        if ( kitty::equal( tt0, tt1 ) )
          continue;

        auto const tt01 = kitty::cofactor0( tt1, i );
        auto const tt10 = kitty::cofactor1( tt0, i );
        if ( kitty::equal( tt01, tt10 ) )
        {
          set( i, j );
        }
      }
    }
  }

  constexpr void set( uint8_t i, uint8_t j ) noexcept
  {
    uint64_t const mask = ( 1lu << j ) | ( 1lu << i );
    data |= ( mask << ( 8 * i ) );
    data |= ( mask << ( 8 * j ) );
  }

  constexpr bool symmetric( uint8_t i, uint8_t j ) const noexcept
  {
    return ( ( ( data >> ( 8 * i ) ) >> j ) & ( ( data >> ( 8 * j ) ) >> i ) & 0x1 ) > 0;
  }

  constexpr bool has_symmetries( uint8_t i ) const noexcept
  {
    return ( ( data >> ( 8 * i ) ) & 0xFF ) > 0;
  }
};

/*! \brief Permutes the input variables to have the ones closest to the output last
 *
 *  if fn = a < b and all variables are symmetric
 *    3
 *   2
 *  1
 * 0
 * */
template<typename E, typename V, typename CompFn>
void sort_symmetric( std::vector<E>& entries, std::vector<V>& values, symmetries_t const& symm, CompFn&& fn )
{
  assert( ( entries.size() == values.size() ) &&
          "[e] entries and values should have the same size" );
  std::vector<uint8_t> inputs( entries.size() );
  std::iota( inputs.begin(), inputs.end(), 0 );

  for ( uint8_t i = 0; i < entries.size(); ++i )
  {
    if ( symm.has_symmetries( i ) )
    {
      uint8_t k = i;
      int j = i - 1;
      V value = values[inputs[i]];
      bool swapped = true;
      while ( swapped && ( j >= 0 ) )
      {
        if ( symm.symmetric( inputs[k], inputs[j] ) )
        {
          if ( fn( value, values[j] ) )
          {
            std::swap( entries[k], entries[j] );
            std::swap( values[k], values[j] );
            std::swap( inputs[k], inputs[j] );
            k = j;
            swapped = true;
          }
          else
          {
            swapped = false;
          }
        }
        j--;
      }
    }
  }
}

/*! \brief Store permutation transformations for up to 16 inputs */
struct permutation_t
{
  uint64_t fmap = 0;
  uint64_t imap = 0;
  uint8_t num_vars = 0;

  permutation_t() = default;
  permutation_t( permutation_t&& ) = default;
  permutation_t( permutation_t const& ) = default;
  permutation_t( std::vector<uint8_t> const& perm )
  {
    num_vars = perm.size();
    for ( uint8_t i{ 0u }; i < perm.size(); ++i )
    {
      set( i, perm[i] );
    }
  }

  permutation_t& operator=( permutation_t const& ) = default;
  permutation_t& operator=( permutation_t&& ) = default;

  constexpr uint8_t forward( uint8_t i ) const noexcept
  {
    return ( fmap >> ( 4 * i ) ) & 0xF;
  }

  constexpr uint8_t inverse( uint8_t i ) const noexcept
  {
    return ( imap >> ( 4 * i ) ) & 0xF;
  }

  constexpr void set( uint8_t i, uint8_t v ) noexcept
  {
    fmap |= static_cast<uint64_t>( v ) << ( 4 * i );
    imap |= static_cast<uint64_t>( i ) << ( 4 * v );
  }

  constexpr bool operator==( const permutation_t& other ) const noexcept
  {
    return ( fmap == other.fmap ) && ( imap == other.imap );
  }
};

template<typename T>
std::vector<T> inverse_permute( permutation_t const& perm, std::vector<T> const& vec )
{
  assert( vec.size() < 16 );
  std::vector<T> res( vec.size() );
  for ( auto i = 0u; i < vec.size(); ++i )
  {
    res[i] = vec[perm.inverse( i )];
  }
  return res;
}

template<typename T>
std::vector<T> forward_permute( permutation_t const& perm, std::vector<T> const& vec )
{
  assert( vec.size() < 16 );
  std::vector<T> res( vec.size() );
  for ( auto i = 0u; i < vec.size(); ++i )
  {
    res[i] = vec[perm.forward( i )];
  }
  return res;
}

template<typename T>
void inverse_permute_inplace_single( permutation_t const& perm, std::vector<T>& vec )
{
  vec = inverse_permute( perm, vec );
}

template<typename... Vecs>
void inverse_permute_inplace( permutation_t const& perm, Vecs&... vecs )
{
  ( inverse_permute_inplace_single( perm, vecs ), ... );
}

template<typename T>
void forward_permute_inplace_single( permutation_t const& perm, std::vector<T>& vec )
{
  vec = forward_permute( perm, vec );
}

template<typename... Vecs>
void forward_permute_inplace( permutation_t const& perm, Vecs&... vecs )
{
  ( forward_permute_inplace_single( perm, vecs ), ... );
}

template<typename T, typename Fn>
void foreach_permutation( std::vector<T>& vec, Fn&& fn )
{
  auto const n = vec.size();
  auto perms = [&]( int i, auto&& perms_ref ) -> void {
    if ( i >= n )
    {
      fn( vec );
      return;
    }

    for ( int j = i; j < n; ++j )
    {
      std::swap( vec[i], vec[j] );
      perms_ref( i + 1, perms_ref );
      std::swap( vec[i], vec[j] ); // backtrack
    }
  };

  perms( 0, perms );
}

} // namespace mockturtle