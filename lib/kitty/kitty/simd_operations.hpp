/* kitty: C++ truth table library
 * Copyright (C) 2017-2022  EPFL
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
  \file simd_operations.hpp
  \brief Implements efficient alternatives of some common operations using SIMD

  \author Andrea Costamagna
*/

#pragma once

#include <cstdint>
#include <cassert>
#include <numeric>
#include <algorithm>

#include "static_truth_table.hpp"
#include "partial_truth_table.hpp"
#include "dynamic_truth_table.hpp"
#include "detail/mscfix.hpp"

#ifndef KITTY_HAS_AVX2
#if defined( __x86_64__ ) || defined( _M_X64 )
#define KITTY_HAS_AVX2 1
#else
#define KITTY_HAS_AVX2 0
#endif
#endif

#if KITTY_HAS_AVX2
#include <immintrin.h>
#include <cpuid.h>
#endif

namespace kitty
{

namespace simd
{
enum class BinaryOp
{
  AND,
  OR,
  XOR,
  LT // Less Than (ANDNOT)
};

inline bool has_avx2_cached()
{
#if KITTY_HAS_AVX2
  static const bool cached = []
  {
    unsigned int eax, ebx, ecx, edx;
    if ( !__get_cpuid_max( 0, nullptr ) || __get_cpuid_max( 0, nullptr ) < 7 )
      return false;
    __cpuid_count( 7, 0, eax, ebx, ecx, edx );
    return ( ebx & ( 1 << 5 ) ) != 0; // Bit 5 of EBX indicates AVX2
  }();
  return cached;
#else
  return false;
#endif
}

#if KITTY_HAS_AVX2
inline __m256i get_simd_op( BinaryOp op, __m256i vr, __m256i v2 )
{
  switch ( op )
  {
  case BinaryOp::AND:
    return _mm256_and_si256( vr, v2 );
  case BinaryOp::OR:
    return _mm256_or_si256( vr, v2 );
  case BinaryOp::XOR:
    return _mm256_xor_si256( vr, v2 );
  case BinaryOp::LT:
    return _mm256_andnot_si256( vr, v2 ); // Using ANDNOT for LT
  default:
    throw std::invalid_argument( "Unsupported operation" );
  }
}
#endif

template<typename TT>
inline TT bitwise_binop( const TT& tta, const TT& ttb, BinaryOp op, std::function<uint64_t( uint64_t, uint64_t )> scalar_op )
{
  TT result = tta;
  size_t size = tta.num_blocks();
  auto& datar = result._bits;
  const auto& data2 = ttb._bits;

  size_t i = 0;
#if KITTY_HAS_AVX2
  if ( has_avx2_cached() )
  {
    for ( ; i + 4 <= size; i += 4 )
    {
      __m256i vr = _mm256_loadu_si256( reinterpret_cast<const __m256i*>( &datar[i] ) );
      __m256i v2 = _mm256_loadu_si256( reinterpret_cast<const __m256i*>( &data2[i] ) );
      vr = get_simd_op( op, vr, v2 );
      _mm256_storeu_si256( reinterpret_cast<__m256i*>( &datar[i] ), vr );
    }
  }
#endif

  // Handle remaining elements with scalar operation
  for ( ; i < size; ++i )
  {
    datar[i] = scalar_op( datar[i], data2[i] );
  }

  result.mask_bits();

  return result;
}

template<typename TT>
inline TT bitwise_and( const TT& tta, const TT& ttb )
{
  return bitwise_binop(
      tta, ttb,
      BinaryOp::AND, // SIMD part for OR
      []( uint64_t a, uint64_t b )
      { return a & b; } // Scalar part for OR
  );
}

template<typename TT>
inline TT bitwise_or( const TT& tta, const TT& ttb )
{
  return bitwise_binop(
      tta, ttb,
      BinaryOp::OR, // SIMD part for OR
      []( uint64_t a, uint64_t b )
      { return a | b; } // Scalar part for OR
  );
}

template<typename TT>
inline TT bitwise_xor( const TT& tta, const TT& ttb )
{
  return bitwise_binop(
      tta, ttb,
      BinaryOp::XOR, // SIMD part for OR
      []( uint64_t a, uint64_t b )
      { return a ^ b; } // Scalar part for OR
  );
}

template<typename TT>
inline TT bitwise_lt( const TT& tta, const TT& ttb )
{
  return bitwise_binop(
      tta, ttb,
      BinaryOp::LT, // SIMD part for OR
      []( uint64_t a, uint64_t b )
      { return ~a & b; } // Scalar part for OR
  );
}

template<typename TT>
inline TT unary_not( const TT& tt )
{
  TT result = tt;
  size_t size = tt.num_blocks();
  auto& data = result._bits;

  size_t i = 0;
#if KITTY_HAS_AVX2
  if ( has_avx2_cached() )
  {
    const __m256i all_ones = _mm256_set1_epi64x( -1 );
    for ( ; i + 4 <= size; i += 4 )
    {
      __m256i v = _mm256_loadu_si256( reinterpret_cast<const __m256i*>( &data[i] ) );
      v = _mm256_xor_si256( v, all_ones ); // ~v == v ^ all_ones
      _mm256_storeu_si256( reinterpret_cast<__m256i*>( &data[i] ), v );
    }
  }
#endif

  for ( ; i < size; ++i )
  {
    data[i] = ~data[i];
  }
  return result;
}

template<typename TT>
inline void set_zero( TT& tt )
{
  size_t size = tt.num_blocks();
  auto& data = tt._bits;

  size_t i = 0;
#if KITTY_HAS_AVX2
  if ( has_avx2_cached() )
  {
    for ( ; i + 4 <= size; i += 4 )
    {
      __m256i v = _mm256_setzero_si256(); // Set the vector to zero
      _mm256_storeu_si256( reinterpret_cast<__m256i*>( &data[i] ), v );
    }
  }
#endif

  // Handle remaining elements
  for ( ; i < size; ++i )
  {
    data[i] = 0;
  }
}

template<typename TT>
inline void set_ones( TT& tt )
{
  size_t size = tt.num_blocks();
  auto& data = tt._bits;

  size_t i = 0;
#if KITTY_HAS_AVX2
  if ( has_avx2_cached() )
  {
    for ( ; i + 4 <= size; i += 4 )
    {
      __m256i v = _mm256_set1_epi64x( -1 );
      _mm256_storeu_si256( reinterpret_cast<__m256i*>( &data[i] ), v );
    }
  }
#endif

  // Handle remaining elements
  for ( ; i < size; ++i )
  {
    data[i] = (uint64_t)( -1 );
  }
}

template<uint32_t NumVars>
inline void set_ones( kitty::static_truth_table<NumVars, true> & tt )
{
  tt |= ~tt;
}

template<uint32_t NumVars>
inline void set_zero( kitty::static_truth_table<NumVars, true> & tt )
{
  tt ^= tt;
}

template<uint32_t NumVars>
inline kitty::static_truth_table<NumVars> unary_not( kitty::static_truth_table<NumVars, true> const& tt )
{
  return ~tt;
}

template<uint32_t NumVars>
inline kitty::static_truth_table<NumVars> bitwise_and( const kitty::static_truth_table<NumVars, true> & tta, const kitty::static_truth_table<NumVars, true> & ttb )
{
  return tta & ttb;
}

template<uint32_t NumVars>
inline kitty::static_truth_table<NumVars> bitwise_or( const kitty::static_truth_table<NumVars, true> & tta, const kitty::static_truth_table<NumVars, true> & ttb )
{
  return tta | ttb;
}

template<uint32_t NumVars>
inline kitty::static_truth_table<NumVars> bitwise_xor( const kitty::static_truth_table<NumVars, true> & tta, const kitty::static_truth_table<NumVars, true> & ttb )
{
  return tta ^ ttb;
}

template<uint32_t NumVars>
inline kitty::static_truth_table<NumVars> bitwise_lt( const kitty::static_truth_table<NumVars, true> & tta, const kitty::static_truth_table<NumVars, true> & ttb )
{
  return ~tta & ttb;
}

} // namespace simd

} // namespace kitty
