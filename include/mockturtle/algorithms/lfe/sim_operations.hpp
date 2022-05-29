/* kitty: C++ truth table library
 * Copyright (C) 2017-2021  EPFL
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
  \file operations.hpp
  \brief Implements several operations on simulation patterns

  \author Andrea Costamagna
*/

#pragma once

#include <algorithm>
#include <cassert>

#include "sim_patterns.hpp"

namespace mockturtle
{

  #pragma region cofactors
  template<typename Ntk>
  std::pair< std::vector<sim_pattern<Ntk>>, std::vector<sim_pattern<Ntk>> > compute_cofactor( std::vector<sim_pattern<Ntk>> const& X, std::vector<sim_pattern<Ntk>> const& F, 
                                                                                      bool const& id, size_t idx )
  {
    if( X.size()==0 )
      return std::make_pair( X, F );
    assert( X[0].pat.num_bits() == F[0].pat.num_bits() );
    assert( idx < X.size() );

    kitty::partial_truth_table M;
    M = id ?  X[idx].pat : ~X[idx].pat;
    std::vector<sim_pattern<Ntk>> Xid;
    std::vector<sim_pattern<Ntk>> Fid;

    uint32_t ref_idx = M.num_bits();
    for( uint32_t i{0u}; i < M.num_bits(); ++i )
    {
      if( kitty::get_bit(M,i) == 1 )
      {
        ref_idx = i;
        break;
      }
    }

    if (ref_idx != M.num_bits() )
    {
      kitty::partial_truth_table rM ( kitty::count_ones(M) );

      for( size_t i{0}; i < X.size(); ++i )
      {
        Xid.push_back( rM );
        Xid[i].sig = X[i].sig;
      }
      for( size_t i{0}; i < F.size(); ++i )
        Fid.push_back( rM );

      uint64_t current_idx = ref_idx;
      size_t k = 0;
      bool new1_found = true;
      do
      {
        for( size_t i = 0; i < F.size(); ++i )
          kitty::get_bit( F[i].pat, current_idx ) == 1 ? kitty::set_bit( Fid[i].pat, k ) : kitty::clear_bit( Fid[i].pat, k );
        for( size_t i{0}; i < X.size(); ++i )
          kitty::get_bit(X[i].pat,current_idx) == 1 ? kitty::set_bit(Xid[i].pat,k) : kitty::clear_bit(Xid[i].pat,k);
        k++;  
        new1_found = false;

        for( auto i{current_idx+1}; i < M.num_bits(); ++i )
        {
          if( kitty::get_bit(M,i) == 1 )
          {
            current_idx = i;
            new1_found = true;
            break;
          }
        }
        if( current_idx >= M.num_bits() )
          new1_found = false;

      } while ( new1_found );
      Xid.erase( Xid.begin()+idx );
    }

    return std::make_pair( Xid, Fid );
  }

  template<typename Ntk>
  std::pair< std::vector<sim_pattern<Ntk>>, sim_pattern<Ntk> > compute_cofactor( std::vector<sim_pattern<Ntk>> const& X, sim_pattern<Ntk> const& F, 
                                                                                      bool const& id, size_t idx )
  {
    if( X.size()==0 )
      return std::make_pair( X, F );
    assert( X[0].pat.num_bits() == F.pat.num_bits() );
    assert( idx < X.size() );

    kitty::partial_truth_table M;
    M = id ?  X[idx].pat : ~X[idx].pat;
    std::vector<sim_pattern<Ntk>> Xid;

    uint32_t ref_idx = M.num_bits();
    for( uint32_t i{0u}; i < M.num_bits(); ++i )
    {
      if( kitty::get_bit(M,i) == 1 )
      {
        ref_idx = i;
        break;
      }
    }

    kitty::partial_truth_table rM ( kitty::count_ones(M) );
    sim_pattern<Ntk> Fid = sim_pattern<Ntk>( rM );

    if (ref_idx != M.num_bits() )
    {

      for( size_t i{0}; i < X.size(); ++i )
      {
        Xid.push_back( rM );
        Xid[i].sig = X[i].sig;
      }

      uint64_t current_idx = ref_idx;
      size_t k = 0;
      bool new1_found = true;
      do
      {
        kitty::get_bit( F.pat, current_idx ) == 1 ? kitty::set_bit( Fid.pat, k ) : kitty::clear_bit( Fid.pat, k );
        for( size_t i{0}; i < X.size(); ++i )
          kitty::get_bit(X[i].pat,current_idx) == 1 ? kitty::set_bit(Xid[i].pat,k) : kitty::clear_bit(Xid[i].pat,k);
        k++;  
        new1_found = false;

        for( auto i{current_idx+1}; i < M.num_bits(); ++i )
        {
          if( kitty::get_bit(M,i) == 1 )
          {
            current_idx = i;
            new1_found = true;
            break;
          }
        }
        if( current_idx >= M.num_bits() )
          new1_found = false;

      } while ( new1_found );
      Xid.erase( Xid.begin()+idx );
    }

    return std::make_pair( Xid, Fid );
  }

  template<typename Ntk>
  std::pair< std::vector<sim_pattern<Ntk>>, std::vector<sim_pattern<Ntk>> > compute_cofactor0( std::vector<sim_pattern<Ntk>> const& X, std::vector<sim_pattern<Ntk>> const& F, size_t idx = 0 )
  {
    return compute_cofactor( X, F, 0, idx );
  }

  template<typename Ntk>
  std::pair< std::vector<sim_pattern<Ntk>>, std::vector<sim_pattern<Ntk>> > compute_cofactor1( std::vector<sim_pattern<Ntk>> const& X, std::vector<sim_pattern<Ntk>> const& F, size_t idx = 0 )
  {
    return compute_cofactor( X, F, 1, idx );
  }

  template<typename Ntk>
  std::pair< std::vector<sim_pattern<Ntk>>, sim_pattern<Ntk> > compute_cofactor0( std::vector<sim_pattern<Ntk>> const& X, sim_pattern<Ntk> const& F, size_t idx = 0 )
  {
    return compute_cofactor( X, F, 0, idx );
  }

  template<typename Ntk>
  std::pair< std::vector<sim_pattern<Ntk>>, sim_pattern<Ntk> > compute_cofactor1( std::vector<sim_pattern<Ntk>> const& X, sim_pattern<Ntk> const& F, size_t idx = 0 )
  {
    return compute_cofactor( X, F, 1, idx );
  }
  #pragma end region cofactors
  
  template<typename Ntk>
  void remove_column_and_invert( std::vector<sim_pattern<Ntk>>& X, sim_pattern<Ntk>& Y, uint32_t idx )
  {
    Y.pat ^= X[idx].pat;
    X.erase( X.begin() + idx );
  }

} // namespace mockturtle