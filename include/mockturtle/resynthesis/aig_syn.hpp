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
  \file aig_syn.hpp
  \brief Synthesis by recursive decomposition for AIGs.

  \author Andrea Costamagna
*/

#pragma once

#include "../utils/index_list.hpp"
#include "../utils/node_map.hpp"
#include "../utils/stopwatch.hpp"

#include <fmt/format.h>
#include <kitty/kitty.hpp>

#include <algorithm>
#include <optional>
#include <type_traits>
#include <vector>

namespace mockturtle
{

template<class Tt_t>
class aig_syn_decompose
{
public:
  using index_list_t = large_xag_index_list;

  explicit aig_syn_decompose() noexcept
  {
    static_assert( std::is_same_v<typename index_list_t::element_type, uint32_t>, "Incompatible element type with the synthesis implementation" );
  }

  /*! \brief Perform AIG synthesis.
   * 
   * Synthesize an AIG index list using a decomposition approach.
   * Differently from the resynthesis engine, this engine guarantees the return
   * of an index list.
   * 
   * TODO: For now decomposition is limited to Shannon's decomposition.
   * Add unateness-guided heuristics and Boolean matching to reduce the list size.
   * 
   * \param target Truth table of the care set to be synthesized.
   */
  index_list_t operator()( Tt_t const& target )
  {
    care = target;
    return synthesize();
  }

private:
  /*! \brief Synthesize the function as an AIG list.
   * 
   * Since this class is intended to be called to synthesize more
   * than one function, the `index_list` argument must be cleared
   * from previous information. Next, decomposition-based synthesis
   * can be performed.
   * 
   */
  index_list_t synthesize()
  {
    index_list.clear();
    int num_vars = care.num_vars();
    index_list.add_inputs( num_vars );
    auto const lit = shannon_decompose( care, num_vars - 1 );
    index_list.add_output( lit );
    return index_list;
  }

  /*! \brief Perform the Shannon decomposition of a truth table.
   * 
   * 
   * 
   */
  uint32_t shannon_decompose( Tt_t const& func, int i )
  {
    if( kitty::is_const0( func ) )
    {
      return 0;
    }
    else if ( kitty::is_const0( ~func ) )
    {
      return 1;
    }
    const Tt_t tt0 = kitty::cofactor0( func, i );
    const Tt_t tt1 = kitty::cofactor1( func, i );
    
    if( kitty::equal( tt0, tt1 ) )
    {
      return shannon_decompose( tt1, i - 1 );
    }
    /* recursively synthesize the cofactors */
    const uint32_t lit0 = shannon_decompose( tt0, i - 1 );
    const uint32_t lit1 = shannon_decompose( tt1, i - 1 );
    const uint32_t lit = ( i + 1 ) << 1;
    /* case-by-case synthesis */
    if( ( lit0 == 0 ) && ( lit1 == 0 ) ) // F = x'0 +x0 = 0
    {
      return 0;
    }
    else if( ( lit0 == 1 ) && ( lit1 == 1 ) ) // F = x'1 +x1 = 1
    {
      return 1;
    }
    else if( ( lit0 == 0 ) && ( lit1 == 1 ) ) // F = x'0 +x1 = x
    {
      return lit;
    }
    else if( ( lit0 == 1 ) && ( lit1 == 0 ) ) // F = x'1 +x0 = x'
    {
      return lit ^ 0x1;
    }
    else if( ( lit0 > 1 ) && ( lit1 == 1 ) ) // F = x'F0 +x1 = F0 + x
    {
      const uint32_t inv = index_list.add_and( lit ^ 0x1, lit0 ^ 0x1 );
      return inv ^ 0x1;
    }
    else if( ( lit0 == 1 ) && ( lit1 > 1 ) ) // F = x'1 +xF1 = F1 + x'
    {
      const uint32_t inv = index_list.add_and( lit, lit1 ^ 0x1 );
      return inv ^ 0x1;
    }
    else if( ( lit0 > 1 ) && ( lit1 == 0 ) ) // F = x'F0 +x0 = x'F0
    {
      return index_list.add_and( lit ^ 0x1, lit0 );
    }
    else if( ( lit0 == 0 ) && ( lit1 > 1 ) ) // F = x'0 +xF1 = xF1
    {
      return index_list.add_and( lit, lit1 );
    }
    const uint32_t and0 = index_list.add_and( lit ^ 0x1, lit0 );
    const uint32_t and1 = index_list.add_and( lit, lit1 );
    const uint32_t inv = index_list.add_and( and1 ^ 0x1, and0 ^ 0x1 );
    const uint32_t res = inv ^ 0x1;
    return res;
  }

private:
  Tt_t care;
  index_list_t index_list;
}; /* aig_syn_decompose */

} /* namespace mockturtle */