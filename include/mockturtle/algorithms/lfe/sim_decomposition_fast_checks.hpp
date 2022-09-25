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
  \file decomposition.hpp
  \brief Check decomposition properties (and perform decomposition) of a function

  \author Andrea Costamagna
*/

#pragma once

#include "sim_utils.hpp"
#include "sim_patterns.hpp"
#include "../../lib/kitty/kitty/statistics.hpp"
#include <set>

namespace mockturtle
{

enum class sim_top_decomposition_fast
{
  none,
  and_,
  or_,
  lt_,
  le_,
  xor_
};

#pragma region is_xor_decomposable_fast
template<typename Ntk, typename TT>
bool is_xor_decomposable_fast( std::vector<sim_pattern<Ntk>>& X, std::vector<uint32_t> & support, TT & onset, TT const& amask1, TT const& amask0 )
{
  uint32_t count_neg=0;
  
  std::set<std::string> onset_Fo; // contains the minterms F0(minterm) = 1
  std::set<std::string> offset_Fo; // contains the minterms F0(minterm) = 0
  std::set<std::string> already; // contains the minterms of F1 that have already been considered

  uint32_t N0 = 0;
  for( uint32_t i = 0; i < amask0.num_bits(); ++i )
  {
    if( kitty::get_bit( amask0, i ) == 1 )
    {
      std::string minterm = "";
      for( uint32_t j{0}; j<support.size(); ++j  )
        minterm += std::to_string( kitty::get_bit( X[support[j]].pat, i ) ); 

      if( kitty::get_bit( onset, i ) == 1 ) // F0(minterm) == 1
      {
        if( offset_Fo.find(minterm) != offset_Fo.end() )
          return false;
        else if( onset_Fo.find(minterm) == onset_Fo.end() )
        {
          N0++;
          onset_Fo.insert(minterm);
        }
      }
      else // F0(minterm)=0
      {
        if( onset_Fo.find(minterm) != onset_Fo.end() )
          return false;
        else if( offset_Fo.find(minterm) == offset_Fo.end() )
        {
          N0++;
          offset_Fo.insert(minterm);
        }
      }
    }
  }

  uint32_t N1 = 0;

  for( uint32_t i = 0; i < amask1.num_bits(); ++i )
  {
    //bool is_new = false;
    if( kitty::get_bit( amask1, i ) == 1 ) // F1(minterm)==1
    {
      std::string minterm = "";
      for( uint32_t j{0}; j<support.size(); ++j  )
        minterm += std::to_string( kitty::get_bit( X[support[j]].pat, i ) ); 

      if( already.find(minterm) == already.end() )
        N1++;

      if( kitty::get_bit( onset, i ) == 1 )
      {
        if( onset_Fo.find( minterm ) != onset_Fo.end() )
          return false;
        else if( (offset_Fo.find( minterm ) != offset_Fo.end())&& ( already.find(minterm) == already.end() )   )
          count_neg++;
      }
      else
      {
        if( offset_Fo.find( minterm ) != offset_Fo.end() )
          return false;
        else
        {
          if( (onset_Fo.find( minterm ) != onset_Fo.end())&& ( already.find(minterm) == already.end() ) ) 
            count_neg++;
        }
      }
      if( already.find(minterm) == already.end() )
        already.insert(minterm);
    }
  }
      
  uint32_t n = support.size()+1;
  std::pair<double,double> R = M1M2k(N0,N1,n);
  //int Min = std::max(1,(int)floor(R.first-R.second));
  //int Max = std::max(1,(int)ceil(R.first+R.second));

  bool is_satisfied = ( ( CumSum( count_neg+(int)ceil(R.second), N0, N1, n ) >= 1-0.001 ) && ( count_neg >1 ) ); // CHANGE !!!

  if( is_satisfied )
    return true;

  return false;
  }  
#pragma endregion is_xor_decomposable_fast

#pragma region is_xor_decomposable_fast
template<typename Ntk, typename TT>
bool is_dc_xor_decomposable_fast( std::vector<sim_pattern<Ntk>>& X, std::vector<uint32_t> & support, TT & onset, TT const& amask1, TT const& amask0 )
{
  uint32_t count_neg=0;
  
  std::set<std::string> onset_Fo; // contains the minterms F0(minterm) = 1
  std::set<std::string> offset_Fo; // contains the minterms F0(minterm) = 0
  std::set<std::string> already; // contains the minterms of F1 that have already been considered

  uint32_t N0 = 0;
  for( uint32_t i = 0; i < amask0.num_bits(); ++i )
  {
    if( kitty::get_bit( amask0, i ) == 1 )
    {
      std::string minterm = "";
      for( uint32_t j{0}; j<support.size(); ++j  )
        minterm += std::to_string( kitty::get_bit( X[support[j]].pat, i ) ); 

      if( kitty::get_bit( onset, i ) == 1 ) // F0(minterm) == 1
      {
        if( offset_Fo.find(minterm) != offset_Fo.end() )
          return false;
        else if( onset_Fo.find(minterm) == onset_Fo.end() )
        {
          N0++;
          onset_Fo.insert(minterm);
        }
      }
      else // F0(minterm)=0
      {
        if( onset_Fo.find(minterm) != onset_Fo.end() )
          return false;
        else if( offset_Fo.find(minterm) == offset_Fo.end() )
        {
          N0++;
          offset_Fo.insert(minterm);
        }
      }
    }
  }

  uint32_t N1 = 0;

  for( uint32_t i = 0; i < amask1.num_bits(); ++i )
  {
    //bool is_new = false;
    if( kitty::get_bit( amask1, i ) == 1 ) // F1(minterm)==1
    {
      std::string minterm = "";
      for( uint32_t j{0}; j<support.size(); ++j  )
        minterm += std::to_string( kitty::get_bit( X[support[j]].pat, i ) ); 

      if( already.find(minterm) == already.end() )
        N1++;

      if( kitty::get_bit( onset, i ) == 1 )
      {
        if( onset_Fo.find( minterm ) != onset_Fo.end() )
          return false;
        else if( (offset_Fo.find( minterm ) != offset_Fo.end())&& ( already.find(minterm) == already.end() )   )
          count_neg++;
      }
      else
      {
        if( offset_Fo.find( minterm ) != offset_Fo.end() )
          return false;
        else
        {
          if( (onset_Fo.find( minterm ) != onset_Fo.end())&& ( already.find(minterm) == already.end() ) ) 
            count_neg++;
        }
      }
      if( already.find(minterm) == already.end() )
        already.insert(minterm);
    }
  }
  
  return true;

  }  
#pragma endregion is_xor_decomposable_fast

#pragma region is_top_decomposable_fast

  template<typename Ntk, typename TT>
  sim_top_decomposition_fast is_top_decomposable_fast( std::vector<sim_pattern<Ntk>> X, std::vector<uint32_t> & support, TT & onset, TT & amask1, TT & amask0, bool try_xor = false )
  {

    if( (kitty::count_ones( onset & amask0 ) == 0) ) // F0 = 0
      return sim_top_decomposition_fast::and_ ;
    else if( (onset & amask1) == amask1 ) // F1 = 1
      return sim_top_decomposition_fast::or_ ;
    else if( kitty::count_ones( onset &amask1 ) == 0 ) // F1 = 0
      return sim_top_decomposition_fast::lt_ ;
    else if( ( onset &amask0 ) == amask0 ) // F0 = 1
      return sim_top_decomposition_fast::le_ ;
    else if( try_xor && is_xor_decomposable_fast( X, support, onset, amask1, amask0 ) )
      return sim_top_decomposition_fast::xor_ ;
    else
      return sim_top_decomposition_fast::none ;
  }
#pragma endregion is_top_decomposable_fast

#pragma region is_dc_top_decomposable_fast

  template<typename Ntk, typename TT>
  sim_top_decomposition_fast is_dc_top_decomposable_fast( std::vector<sim_pattern<Ntk>> X, std::vector<uint32_t> & support, TT & onset, TT & amask1, TT & amask0, bool try_xor = false )
  {
    if( (kitty::count_ones( onset & amask0 ) == 0) ) // F0 = 0
      return sim_top_decomposition_fast::and_ ;
    else if( (onset & amask1) == amask1 ) // F1 = 1
      return sim_top_decomposition_fast::or_ ;
    else if( kitty::count_ones( onset &amask1 ) == 0 ) // F1 = 0
      return sim_top_decomposition_fast::lt_ ;
    else if( ( onset &amask0 ) == amask0 ) // F0 = 1
      return sim_top_decomposition_fast::le_ ;
    else if( try_xor && is_dc_xor_decomposable_fast( X, support, onset, amask1, amask0 ) )
      return sim_top_decomposition_fast::xor_ ;
    else
      return sim_top_decomposition_fast::none ;
  }
#pragma endregion is_top_decomposable_fast

#pragma endregion is_bottom_decomposable

} // namespace mockturtle