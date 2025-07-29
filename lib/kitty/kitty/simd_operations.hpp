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
  \brief Implements efficient alternatives of common truth-table operations.

  \author Andrea Costamagna
*/

#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <mutex>
#include <numeric>
#include <unordered_map>

#include "detail/mscfix.hpp"
#include "dynamic_truth_table.hpp"
#include "operators.hpp"
#include "partial_truth_table.hpp"
#include "static_truth_table.hpp"

/*! Check if the code is compiled on a platform that might support AVX2 */
#ifndef KITTY_HAS_AVX2
#if defined( __x86_64__ ) || defined( _M_X64 ) || defined( _MSC_VER )
#define KITTY_HAS_AVX2 1
#else
#define KITTY_HAS_AVX2 0
#endif
#endif

#if KITTY_HAS_AVX2
#if defined( _MSC_VER )
#include <intrin.h>
#else
#include <cpuid.h>
#include <immintrin.h>
#endif
#endif

namespace kitty
{

namespace simd
{

/*! \brief Enumeration for compile-time dispatch and minimize code redundancies. */
enum class Operation : uint32_t
{
  AND,
  OR,
  XOR,
  LT,
  NOT,
  CONST0,
  CONST1,
  SIZE
};

/*! Check if AVX2 is supported on this machine. */
inline bool has_avx2_cached()
{
#if KITTY_HAS_AVX2
#if defined( _MSC_VER )
  static const bool cached = [] {
    int cpuInfo[4];
    __cpuid( cpuInfo, 0 );
    if ( cpuInfo[0] < 7 )
      return false;
    __cpuidex( cpuInfo, 7, 0 );
    return ( cpuInfo[1] & ( 1 << 5 ) ) != 0;
  }();
  return cached;
#else
  static const bool cached = [] {
    unsigned int eax, ebx, ecx, edx;
    unsigned int max_leaf = __get_cpuid_max( 0, nullptr );
    if ( !max_leaf || max_leaf < 7 )
      return false;
    __cpuid_count( 7, 0, eax, ebx, ecx, edx );
    return ( ebx & ( 1 << 5 ) ) != 0;
  }();
  return cached;
#endif
#else
  return false;
#endif
}

/*! Check if AVX2 is supported on this machine. */
template<Operation Op, typename TT>
inline bool use_avx2_cached( TT const&, uint32_t );

/*! Compile-time dispatch of the vector operations. */
#if KITTY_HAS_AVX2
template<Operation Op>
constexpr auto vector_operation()
{
  if constexpr ( Op == Operation::AND )
    return []( __m256i a, __m256i b ) { return _mm256_and_si256( a, b ); };
  else if constexpr ( Op == Operation::OR )
    return []( __m256i a, __m256i b ) { return _mm256_or_si256( a, b ); };
  else if constexpr ( Op == Operation::XOR )
    return []( __m256i a, __m256i b ) { return _mm256_xor_si256( a, b ); };
  else if constexpr ( Op == Operation::LT )
    return []( __m256i a, __m256i b ) { return _mm256_andnot_si256( a, b ); };
}
#endif

/*! Compile-time dispatch of the scalar operations. */
template<Operation Op, typename T>
constexpr auto scalar_operation()
{
  if constexpr ( Op == Operation::AND )
    return std::bit_and<T>{};
  else if constexpr ( Op == Operation::OR )
    return std::bit_or<T>{};
  else if constexpr ( Op == Operation::XOR )
    return std::bit_xor<T>{};
  else if constexpr ( Op == Operation::LT )
    return []( T a, T b ) { return ~a & b; };
}

/*! \brief Universal function for vectorized operations between two truth tables.
 *
 * Computes the bitwise operation \f$tt_a Op tt_b\f$ using 256-bit AVX2 registers.
 * Each register processes four 64-bit words in parallel, enabling efficient
 * parallel computation across the truth tables. When the number of bits is not
 * a multiple of four, the function fallsback to the scalar version on the tailing
 * words.
 *
 * \tparam Op Binary operation.
 * \tparam TT Truth table type.
 * \param tta First truth table.
 * \param ttb Second truth table.
 * \return The truth table obtained by applying Op.
 */
template<Operation Op, bool UseCache, typename TT>
inline TT binary_operation( const TT& tta, const TT& ttb )
{
  assert( tta.num_blocks() == ttb.num_blocks() );

  TT result = tta;
  auto& datar = result._bits;
  const auto& data2 = ttb._bits;
  using T = typename std::decay_t<decltype( datar[0] )>;

  size_t i = 0;
#if KITTY_HAS_AVX2
  size_t size = tta.num_blocks();
  if ( has_avx2_cached() && size >= 4 && ( !UseCache || use_avx2_cached<Op, TT>( tta, tta.num_vars() ) ) )
  {
    auto binary_op = vector_operation<Op>();
    for ( ; i + 3 < size; i += 4 )
    {
      __m256i vr = _mm256_loadu_si256( reinterpret_cast<const __m256i*>( &datar[i] ) );
      __m256i v2 = _mm256_loadu_si256( reinterpret_cast<const __m256i*>( &data2[i] ) );
      vr = binary_op( vr, v2 );
      _mm256_storeu_si256( reinterpret_cast<__m256i*>( &datar[i] ), vr );
    }
  }
#endif
  /* Fallback to the scalar version for the remaining words. */
  auto scalar_op = scalar_operation<Op, T>();
  assert( datar.size() == data2.size() );
  assert( datar.size() >= i );
  std::transform( datar.cbegin() + i, datar.cend(), data2.cbegin() + i, datar.begin() + i, scalar_op );

  result.mask_bits();
  return result;
}

/*! \brief Perform a vectorized bitwise AND between two truth tables.
 *
 * Computes the bitwise AND \f$tt_a \land tt_b\f$ using 256-bit AVX2 registers.
 * Each register processes four 64-bit words in parallel, enabling efficient
 * parallel computation across the truth tables.
 *
 * \tparam TT Truth table type.
 * \param tta First truth table.
 * \param ttb Second truth table.
 * \return The truth table obtained by applying the binary AND.
 */
template<typename TT, bool UseCache = true>
inline TT binary_and( const TT& tta, const TT& ttb )
{
  return binary_operation<Operation::AND, UseCache>( tta, ttb );
}

/*! \brief Perform a vectorized bitwise OR between two truth tables.
 *
 * Computes the bitwise OR \f$tt_a \lor tt_b\f$ using 256-bit AVX2 registers.
 * Each register processes four 64-bit words in parallel, enabling efficient
 * parallel computation across the truth tables.
 *
 * \tparam TT Truth table type.
 * \param tta First truth table.
 * \param ttb Second truth table.
 * \return The truth table obtained by applying the binary OR.
 */
template<typename TT, bool UseCache = true>
inline TT binary_or( const TT& tta, const TT& ttb )
{
  return binary_operation<Operation::OR, UseCache>( tta, ttb );
}

/*! \brief Perform a vectorized bitwise XOR between two truth tables.
 *
 * Computes the bitwise XOR \f$tt_a \lxor tt_b\f$ using 256-bit AVX2 registers.
 * Each register processes four 64-bit words in parallel, enabling efficient
 * parallel computation across the truth tables.
 *
 * \tparam TT Truth table type.
 * \param tta First truth table.
 * \param ttb Second truth table.
 * \return The truth table obtained by applying the binary XOR.
 */
template<typename TT, bool UseCache = true>
inline TT binary_xor( const TT& tta, const TT& ttb )
{
  return binary_operation<Operation::XOR, UseCache>( tta, ttb );
}

/*! \brief Perform a vectorized bitwise LT ( Lower Than ) between two truth tables.
 *
 * Computes the bitwise LT \f$ ~tt_a \land tt_b\f$ using 256-bit AVX2 registers.
 * Each register processes four 64-bit words in parallel, enabling efficient
 * parallel computation across the truth tables.
 *
 * \tparam TT Truth table type.
 * \param tta First truth table.
 * \param ttb Second truth table.
 * \return The truth table obtained by applying the binary LT.
 */
template<typename TT, bool UseCache = true>
inline TT binary_lt( const TT& tta, const TT& ttb )
{
  return binary_operation<Operation::LT, UseCache>( tta, ttb );
}

/*! \brief Perform a vectorized inversion of a truth tables.
 *
 * Computes the inverse \f$ ~tt \f$ using 256-bit AVX2 registers.
 * Each register processes four 64-bit words in parallel, enabling efficient
 * parallel inversion across the truth tables.
 *
 * \tparam TT Truth table type.
 * \param tt Truth table.
 * \return The negated truth table.
 */
template<typename TT, bool UseCache = true>
inline TT unary_not( const TT& tt )
{
  TT result = tt;
  size_t size = tt.num_blocks();
  auto& data = result._bits;

  size_t i = 0;
#if KITTY_HAS_AVX2
  if ( has_avx2_cached() && ( !UseCache || use_avx2_cached<Operation::NOT, TT>( tt, tt.num_vars() ) ) )
  {
    const __m256i all_ones = _mm256_set1_epi64x( -1 );
    for ( ; i + 3 < size; i += 4 )
    {
      __m256i v = _mm256_loadu_si256( reinterpret_cast<const __m256i*>( &data[i] ) );
      v = _mm256_xor_si256( v, all_ones );
      _mm256_storeu_si256( reinterpret_cast<__m256i*>( &data[i] ), v );
    }
  }
#endif
  std::transform( data.cbegin() + i, data.cend(), result.begin() + i, std::bit_not<>() );

  result.mask_bits();
  return result;
}

/*! \brief Vectorized set of a truth-table to a constant value.
 *
 * Assign the bits of a truth table to a constant using 256-bit AVX2 registers.
 * Each register processes four 64-bit words in parallel, enabling efficient
 * parallel assignment across the truth tables.
 *
 * \tparam TT Truth table type.
 * \tparam Const Constant value to be assigned ( 0 for contradiction, anything for tautology ).
 * \param tt Truth table.
 */
template<typename TT, Operation Op, bool UseCache = true>
inline void set_const( TT& tt )
{
  size_t size = tt.num_blocks();
  auto& data = tt._bits;

  size_t i = 0;
#if KITTY_HAS_AVX2

  if ( has_avx2_cached() && ( !UseCache || use_avx2_cached<Op, TT>( tt, tt.num_vars() ) ) )
  {
    __m256i v;
    if constexpr ( Op == Operation::CONST0 )
    {
      v = _mm256_setzero_si256();
    }
    else
    {
      v = _mm256_set1_epi64x( -1 );
    }

    for ( ; i + 3 < size; i += 4 )
    {
      _mm256_storeu_si256( reinterpret_cast<__m256i*>( &data[i] ), v );
    }
  }
#endif
  uint64_t v;

  if constexpr ( Op == Operation::CONST0 )
  {
    v = 0;
  }
  else
  {
    v = (uint64_t)( -1 );
  }

  // Handle remaining elements
  if ( tt.num_blocks() > i )
    std::fill( data.begin() + i, data.end(), v );
}

/*! \brief Reset all the bits of a truth table to 0 through vectorization.
 *
 * Set all the bits of a truth table to 0 using 256-bit AVX2 registers.
 * Each register processes four 64-bit words in parallel, enabling efficiently
 * setting the truth table to the desird value.
 *
 * \tparam TT Truth table type.
 * \param tt Truth table.
 */
template<typename TT, bool UseCache = true>
inline void set_zero( TT& tt )
{
  set_const<TT, Operation::CONST0, UseCache>( tt );
}

/*! \brief Reset all the bits of a truth table to 1 through vectorization.
 *
 * Set all the bits of a truth table to 1 using 256-bit AVX2 registers.
 * Each register processes four 64-bit words in parallel, enabling efficiently
 * setting the truth table to the desird value.
 *
 * \tparam TT truth table type.
 * \param tt Truth table.
 */
template<typename TT, bool UseCache = true>
inline void set_ones( TT& tt )
{
  set_const<TT, Operation::CONST1, UseCache>( tt );
}

class benchmarking
{
  static constexpr auto num_cases = 100u;
  static constexpr double eps = 0.1;

public:
  template<typename FnSisd, typename FnSimd, typename TT>
  bool test_noreturn( FnSisd fn_sisd, FnSimd fn_simd, TT const& tt ) const
  {
    TT tt1 = tt.construct();
    TT tt2 = tt.construct();
    double time_diff = 0;
    double time_sisd = 0;
    double time_simd = 0;
    for ( auto i = 0u; i < num_cases; ++i )
    {
      create_random( tt1, i );
      create_random( tt2, i );

      run_noreturn_with_time<FnSisd, TT>( fn_sisd, tt1, time_sisd );
      run_noreturn_with_time<FnSimd, TT>( fn_simd, tt2, time_simd );

      time_diff += ( time_simd - time_sisd ) / time_sisd / static_cast<double>( num_cases );
    }
#if KITTY_HAS_AVX2
    bool const has_avx2 = has_avx2_cached();
    if ( has_avx2 )
    {
      return time_diff < -eps;
    }
#endif
    return false;
  }

  template<typename FnSisd, typename FnSimd, typename TT>
  bool test_unary( FnSisd fn_sisd, FnSimd fn_simd, TT tt ) const
  {
    double time_diff = 0;
    double time_sisd = 0;
    double time_simd = 0;
    TT tt1 = tt.construct();
    for ( auto i = 0u; i < num_cases; ++i )
    {
      create_random( tt1, i );

      run_with_time<FnSisd, TT>( fn_sisd, tt1, time_sisd );
      run_with_time<FnSimd, TT>( fn_simd, tt1, time_simd );

      time_diff += ( time_simd - time_sisd ) / time_sisd / static_cast<double>( num_cases );
    }
#if KITTY_HAS_AVX2
    bool const has_avx2 = has_avx2_cached();
    if ( has_avx2 )
    {
      return time_diff < -eps;
    }
#endif
    return false;
  }

  template<typename FnSisd, typename FnSimd, typename TT>
  bool test_binary( FnSisd fn_sisd, FnSimd fn_simd, TT tt ) const
  {
    double time_diff = 0;
    double time_sisd = 0;
    double time_simd = 0;
    TT tt1 = tt.construct();
    TT tt2 = tt.construct();
    for ( auto i = 0u; i < num_cases; ++i )
    {
      create_random( tt1, i );

      run_with_time<FnSisd, TT>( fn_sisd, tt1, tt2, time_sisd );
      run_with_time<FnSimd, TT>( fn_simd, tt1, tt2, time_simd );

      time_diff += ( time_simd - time_sisd ) / time_sisd / static_cast<double>( num_cases );
    }
#if KITTY_HAS_AVX2
    bool const has_avx2 = has_avx2_cached();
    if ( has_avx2 )
    {
      return time_diff < -eps;
    }
#endif
    return false;
  }

  template<typename F, typename TT>
  auto run_with_time( F func, TT& tt, double& t ) const
  {
    auto start = std::chrono::high_resolution_clock::now();
    auto const res = func( tt );
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    t = elapsed.count();
    return res;
  }

  template<typename F, typename TT>
  auto run_with_time( F func, TT& tt1, TT& tt2, double& t ) const
  {
    auto start = std::chrono::high_resolution_clock::now();
    auto const res = func( tt1, tt2 );
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    t = elapsed.count();
    return res;
  }

  template<typename F, typename TT>
  void run_noreturn_with_time( F func, TT& tt, double& t ) const
  {
    auto start = std::chrono::high_resolution_clock::now();
    func( tt );
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    t = elapsed.count();
  }
};

template<Operation Op, typename TT, typename FnSisd, typename FnSimd>
inline bool use_avx2_cached( TT const& tt, FnSisd fn_sisd, FnSimd fn_simd, uint32_t num_vars )
{
  static std::unordered_map<int, bool> cache;
  static std::mutex mutex;

  {
    std::lock_guard<std::mutex> lock( mutex );
    auto it = cache.find( num_vars );
    if ( it != cache.end() )
    {
      return it->second;
    }
  }

  // Not cached yet compute it
  benchmarking b;
  bool result = false;

  if constexpr ( Op == Operation::CONST0 || Op == Operation::CONST1 )
  {
    result = b.test_noreturn(
        fn_sisd,
        fn_simd,
        tt );
  }
  else if constexpr ( Op == Operation::NOT )
  {
    result = b.test_unary(
        fn_sisd,
        fn_simd,
        tt );
  }
  else
  {
    result = b.test_binary(
        fn_sisd,
        fn_simd,
        tt );
  }

  {
    std::lock_guard<std::mutex> lock( mutex );
    cache[num_vars] = result;
  }

  return result;
}

template<Operation Op, typename TT>
struct use_avx2_cached_impl
{
  bool eval( TT const& tt, uint32_t num_vars )
  {
    static_assert( Op != Op, "use_avx2_cached_impl not specialized for this Operation" );
    return false;
  }
};

template<typename TT>
struct use_avx2_cached_impl<Operation::AND, TT>
{
  bool eval( TT const& tt, uint32_t num_vars )
  {
    return use_avx2_cached<Operation::AND>( tt, []( const TT& t1, const TT& t2 ) { return kitty::binary_and( t1, t2 ); }, []( const TT& t1, const TT& t2 ) { return simd::binary_and<TT, false>( t1, t2 ); }, num_vars );
  }
};

template<typename TT>
struct use_avx2_cached_impl<Operation::OR, TT>
{
  bool eval( const TT& tt, uint32_t num_vars )
  {
    return use_avx2_cached<Operation::OR>( tt, []( const TT& t1, const TT& t2 ) { return kitty::binary_or( t1, t2 ); }, []( const TT& t1, const TT& t2 ) { return simd::binary_or<TT, false>( t1, t2 ); }, num_vars );
  }
};

template<typename TT>
struct use_avx2_cached_impl<Operation::XOR, TT>
{
  bool eval( const TT& tt, uint32_t num_vars )
  {
    return use_avx2_cached<Operation::XOR>( tt, []( const TT& t1, const TT& t2 ) { return kitty::binary_xor( t1, t2 ); }, []( const TT& t1, const TT& t2 ) { return simd::binary_xor<TT, false>( t1, t2 ); }, num_vars );
  }
};

template<typename TT>
struct use_avx2_cached_impl<Operation::LT, TT>
{
  bool eval( const TT& tt, uint32_t num_vars )
  {
    return use_avx2_cached<Operation::LT>( tt, []( const TT& t1, const TT& t2 ) { return kitty::binary_and( ~t1, t2 ); }, []( const TT& t1, const TT& t2 ) { return simd::binary_lt<TT, false>( t1, t2 ); }, num_vars );
  }
};

template<typename TT>
struct use_avx2_cached_impl<Operation::NOT, TT>
{
  bool eval( const TT& tt, uint32_t num_vars )
  {
    return use_avx2_cached<Operation::NOT>( tt, []( const TT& t ) { return kitty::unary_not( t ); }, []( const TT& t ) { return simd::unary_not<TT, false>( t ); }, num_vars );
  }
};

template<typename TT>
struct use_avx2_cached_impl<Operation::CONST0, TT>
{
  bool eval( const TT& tt, uint32_t num_vars )
  {
    return use_avx2_cached<Operation::CONST0>( tt, []( TT& t ) { t = t ^ t; }, []( TT& t ) { simd::set_zero<TT, false>( t ); }, num_vars );
  }
};

template<typename TT>
struct use_avx2_cached_impl<Operation::CONST1, TT>
{
  bool eval( const TT& tt, uint32_t num_vars )
  {
    return use_avx2_cached<Operation::CONST1>( tt, []( TT& t ) { t = t ^ ~t; }, []( TT& t ) { simd::set_ones<TT, false>( t ); }, num_vars );
  }
};

template<Operation Op, typename TT>
inline bool use_avx2_cached( TT const& tt, uint32_t num_vars )
{
  if ( tt.num_vars() <= 6u )
    return false;
  use_avx2_cached_impl<Op, TT> impl;
  return impl.eval( tt, num_vars );
}

/*! \brief Test and caches if it is better to use the scalar or the vector version.
 *
 * Each operation is tested once for the specified truth table type and for the
 * specified number of variables. The benchmarking determines if the vectorized ( AVX2 )
 * implementation should be preferred for this machin, truth table size, and truth table
 * type.
 *
 * \tparam TT Truth table type.
 * \param tt Reference truth table needed for construction purposes.
 * \param num_vars Number of variables in the truth table.
 */
template<typename TT>
void test_avx2_advantage( TT const& tt, uint32_t num_vars )
{
  use_avx2_cached<Operation::AND, TT>( tt, num_vars );
  use_avx2_cached<Operation::OR, TT>( tt, num_vars );
  use_avx2_cached<Operation::XOR, TT>( tt, num_vars );
  use_avx2_cached<Operation::LT, TT>( tt, num_vars );
  use_avx2_cached<Operation::NOT, TT>( tt, num_vars );
  use_avx2_cached<Operation::CONST0, TT>( tt, num_vars );
  use_avx2_cached<Operation::CONST1, TT>( tt, num_vars );
}

/*! Fallback to the default unary NOT for small static truth tables. */
template<uint32_t NumVars>
inline kitty::static_truth_table<NumVars> unary_not( kitty::static_truth_table<NumVars, true> const& tt )
{
  return ~tt;
}

/*! Fallback to the default bitwise AND for small static truth tables. */
template<uint32_t NumVars>
inline kitty::static_truth_table<NumVars> binary_and( const kitty::static_truth_table<NumVars, true>& tta, const kitty::static_truth_table<NumVars, true>& ttb )
{
  return tta & ttb;
}

/*! Fallback to the default bitwise OR for small static truth tables. */
template<uint32_t NumVars>
inline kitty::static_truth_table<NumVars> binary_or( const kitty::static_truth_table<NumVars, true>& tta, const kitty::static_truth_table<NumVars, true>& ttb )
{
  return tta | ttb;
}

/*! Fallback to the default bitwise XOR for small static truth tables. */
template<uint32_t NumVars>
inline kitty::static_truth_table<NumVars> binary_xor( const kitty::static_truth_table<NumVars, true>& tta, const kitty::static_truth_table<NumVars, true>& ttb )
{
  return tta ^ ttb;
}

/*! Fallback to the default bitwise LT for small static truth tables. */
template<uint32_t NumVars>
inline kitty::static_truth_table<NumVars> binary_lt( const kitty::static_truth_table<NumVars, true>& tta, const kitty::static_truth_table<NumVars, true>& ttb )
{
  return ~tta & ttb;
}

/*! Implementation set to constant 1 for small static truth tables. */
template<uint32_t NumVars>
inline void set_ones( kitty::static_truth_table<NumVars, true>& tt )
{
  tt |= ~tt;
}

/*! Implementation set to constant 0 for small static truth tables. */
template<uint32_t NumVars>
inline void set_zero( kitty::static_truth_table<NumVars, true>& tt )
{
  tt ^= tt;
}

} // namespace simd

} // namespace kitty