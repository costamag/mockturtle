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
  \file homo_xpx_analyzer.hpp
  \brief Inplace rewrite

  \author Andrea Costamagna
*/

#pragma once

#include "../database/mapped_list.hpp"

#include "../windowing/window_manager.hpp"
#include "../../algorithms/simulation.hpp"
#include "../list_simulator.hpp"
#include "../timing/sensing_analyzer.hpp"
#include "../timing/arrival_analyzer.hpp"
#include "../power/node_simulator.hpp"
#include "../../algorithms/detail/resub_utils.hpp"
#include "../../algorithms/reconv_cut.hpp"

#include <fmt/format.h>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>

namespace mockturtle
{

  template<uint32_t cube_size, uint32_t num_vars>
  class spfd_manager_t
  {
public:
    using truth_table_t = kitty::static_truth_table<cube_size>;
    static constexpr uint32_t capacity = 1u << num_vars;

    spfd_manager_t()
      :  masks(),
         killed(),
         num_masks( 0 ),
         num_killed( 0 ),
         num_edges( 0 ),
         care(),
         safe_care(),
         func()
    {}
    
    void init( truth_table_t const& target )
    {
      func[1] =  target;
      func[0] = ~target;
      care = func[1]|~func[1];
      safe_care = care;
      reset();
    }

    void init( truth_table_t const& target, truth_table_t const& careset )
    {
      care = careset;
      safe_care = careset;
      func[1] =  target & careset;
      func[0] = ~target & careset;
      reset();
    }

    void reset()
    {
      masks[0] = safe_care;
      num_masks = 1;
      num_edges = kitty::count_ones( func[1] ) * kitty::count_ones( func[0] );  
      killed[0] = num_edges > 0 ? false : true;
      num_killed = num_edges > 0 ? 0 : 1;
    }

    void reset( truth_table_t const& modified_care, bool complement )
    {
      masks[0] = complement ? safe_care & ~modified_care : safe_care & modified_care;
      num_masks = 1;
      num_edges = kitty::count_ones( func[1] & masks[0] ) * kitty::count_ones( func[0] & masks[0] );  
      killed[0] = num_edges > 0 ? false : true;
      num_killed = num_edges > 0 ? 0 : 1;
    }

    bool update( truth_table_t const& tt )
    {
      num_edges = 0;
      for( uint32_t iMask{0}; iMask < num_masks; ++iMask )
      {
        if( killed[iMask] )
        {
          killed[num_masks+iMask] = true;
          num_killed++;
        }
        else
        {
          masks[num_masks+iMask] = masks[iMask] & tt;
          masks[iMask] &= ~tt;

          if( ( kitty::count_ones( masks[num_masks+iMask] & func[1] ) == 0 ) || ( kitty::count_ones( masks[num_masks+iMask] & func[0] ) == 0 ) )
          {
            killed[num_masks+iMask] = true;
            num_killed++;
          }
          else
          {
            killed[num_masks+iMask] = false;
            num_edges += kitty::count_ones( func[1] & masks[num_masks+iMask] ) * kitty::count_ones( func[0] & masks[num_masks+iMask] );  
          }

          if( kitty::count_ones( masks[iMask] & func[1] ) == 0 || kitty::count_ones( masks[iMask] & func[0] ) == 0 )
          {
            killed[iMask] = true;
            num_killed++;
          }
          else
          {
            killed[iMask] = false;
            num_edges += kitty::count_ones( func[1] & masks[iMask] ) * kitty::count_ones( func[0] & masks[iMask] );  
          }
        }
      }
      num_masks = num_masks << 1u;
      return true;
    }

    uint32_t evaluate( truth_table_t const& tt )
    {
      uint32_t res=0;
      for( auto iMask{0}; iMask<num_masks; ++iMask )
      {
        if( !killed[iMask] )
        {
          res += kitty::count_ones( func[1] & masks[iMask] &  tt ) * kitty::count_ones( func[0] & masks[iMask] & tt );  
          res += kitty::count_ones( func[1] & masks[iMask] & ~tt ) * kitty::count_ones( func[0] & masks[iMask] & ~tt );  
        }
      }
      return res;
    } 

    uint32_t evaluate( truth_table_t const& tt1, truth_table_t const& tt2 )
    {
      uint32_t res=0;
      for( auto iMask{0}; iMask<num_masks; ++iMask )
      {
        if( !killed[iMask] )
        {
          res+= kitty::count_ones( func[1] & masks[iMask] &  tt1 & tt2 ) * kitty::count_ones( func[0] & masks[iMask] & tt1 & tt2 );  
          res+= kitty::count_ones( func[1] & masks[iMask] & ~tt1 & tt2 ) * kitty::count_ones( func[0] & masks[iMask] & ~tt1 & tt2 );  
          res+= kitty::count_ones( func[1] & masks[iMask] & ~tt1 &~tt2 ) * kitty::count_ones( func[0] & masks[iMask] & ~tt1 & ~tt2 );  
          res+= kitty::count_ones( func[1] & masks[iMask] & tt1 & ~tt2 ) * kitty::count_ones( func[0] & masks[iMask] & tt1 &~tt2 );  
        }
      }
      return res;
    } 

    uint32_t evaluate( truth_table_t const& tt1, truth_table_t const& tt2, truth_table_t const& tt3 )
    {
      uint32_t res=0;
      for( auto iMask{0}; iMask<num_masks; ++iMask )
      {
        if( !killed[iMask] )
        {
          res+= kitty::count_ones( func[1] & masks[iMask] &  tt1 &  tt2 &  tt3 ) * kitty::count_ones( func[0] & masks[iMask] &  tt1 &  tt2 &  tt3 );  
          res+= kitty::count_ones( func[1] & masks[iMask] & ~tt1 &  tt2 &  tt3 ) * kitty::count_ones( func[0] & masks[iMask] & ~tt1 &  tt2 &  tt3 );  
          res+= kitty::count_ones( func[1] & masks[iMask] & ~tt1 & ~tt2 &  tt3 ) * kitty::count_ones( func[0] & masks[iMask] & ~tt1 & ~tt2 &  tt3 );  
          res+= kitty::count_ones( func[1] & masks[iMask] &  tt1 & ~tt2 &  tt3 ) * kitty::count_ones( func[0] & masks[iMask] &  tt1 & ~tt2 &  tt3 );  
          res+= kitty::count_ones( func[1] & masks[iMask] &  tt1 &  tt2 & ~tt3 ) * kitty::count_ones( func[0] & masks[iMask] &  tt1 &  tt2 & ~tt3 );  
          res+= kitty::count_ones( func[1] & masks[iMask] & ~tt1 &  tt2 & ~tt3 ) * kitty::count_ones( func[0] & masks[iMask] & ~tt1 &  tt2 & ~tt3 );  
          res+= kitty::count_ones( func[1] & masks[iMask] & ~tt1 & ~tt2 & ~tt3 ) * kitty::count_ones( func[0] & masks[iMask] & ~tt1 & ~tt2 & ~tt3 );  
          res+= kitty::count_ones( func[1] & masks[iMask] &  tt1 & ~tt2 & ~tt3 ) * kitty::count_ones( func[0] & masks[iMask] &  tt1 & ~tt2 & ~tt3 );  
        }
      }
      return res;
    } 

    bool is_covered()
    {
      return num_masks <= num_killed;
    }

    bool is_saturated()
    {
      return num_masks >= capacity;
    }

    uint32_t get_num_edges() const
    {
      return num_edges;
    }

    std::array<truth_table_t, capacity> masks;
    std::array<bool, capacity> killed;
    uint32_t num_masks;
    uint32_t num_killed;
    uint32_t num_edges;
    truth_table_t care;
    truth_table_t safe_care;
    std::array<truth_table_t, 2u> func;
  };

} // namespace mockturtle