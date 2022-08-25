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
#include "../../lib/kitty/kitty/statistics.hpp"

namespace mockturtle
{

enum class sim_top_decomposition
{
  none,
  and_,
  or_,
  lt_,
  le_,
  xor_
};

#pragma region is_xor_decomposable
template<typename Ntk>
bool is_xor_decomposable(  std::pair<std::vector<sim_pattern<Ntk>>, sim_pattern<Ntk>> const& XY0, std::pair<std::vector<sim_pattern<Ntk>>, sim_pattern<Ntk>> const& XY1 )
{
  size_t count_neg = 0;
  std::unordered_map<std::string, bool> str_nodes0;
  std::unordered_map<std::string, bool> already;

  size_t N0 = 0;
  for( size_t k {0u}; k < XY0.first[0].pat.num_bits(); ++k )
  {
    kitty::partial_truth_table pattern;
    for( size_t j{0}; j < XY0.first.size(); ++j )
      pattern.add_bit( kitty::get_bit(XY0.first[j].pat,k) );
    std::string s;
    s = kitty::to_binary( pattern );
    if( str_nodes0.find(s) == str_nodes0.end() )
    {
      N0++;
      str_nodes0.insert( std::make_pair(s, kitty::get_bit(XY0.second.pat,k)));
    }
    else if ( str_nodes0.at(s) != kitty::get_bit(XY0.second.pat, k) )
    {
      return false;
    }
  }

  uint32_t N1 = 0;
  for ( size_t k {0u}; k < XY1.first[0].pat.num_bits(); ++k )
  {
    kitty::partial_truth_table pattern;
    for( size_t j{0}; j < XY1.first.size(); ++j )
      pattern.add_bit( kitty::get_bit(XY1.first[j].pat, k) );
    std::string s;
    s = kitty::to_binary( pattern );
    if( already.find(s) == already.end() )
      N1++;

    if( str_nodes0.find(s) != str_nodes0.end() )
    {   
      if( str_nodes0.at(s) == kitty::get_bit(XY1.second.pat,k) )
        return false;
      else
      {
        if(already.find(s) == already.end())
          count_neg++;
      }
    }
    if(already.find(s)==already.end())
      already.insert(std::make_pair(s, kitty::get_bit( XY1.second.pat, k )));
  }
  
  size_t n = XY0.first.size()+1;
  std::pair<double,double> R = M1M2k(N0,N1,n);
  int Min = std::max(1,(int)floor(R.first-R.second));
  int Max = std::max(1,(int)ceil(R.first+R.second));

  bool is_satisfied = ( ( CumSum( count_neg+(int)ceil(R.second), N0, N1, n ) >= 1-0.001 ) && ( count_neg >=2 ) );

  if( is_satisfied )
    return true;

  return false;
  }

template<typename Ntk>
bool is_xor_decomposable(  std::pair<std::vector<sim_pattern<Ntk>>, std::vector<sim_pattern<Ntk>>> const& XY0, std::pair<std::vector<sim_pattern<Ntk>>, std::vector<sim_pattern<Ntk>>> const& XY1 )
{
  size_t count_neg = 0;
  std::unordered_map<std::string, bool> str_nodes0;
  std::unordered_map<std::string, bool> already;

  size_t N0 = 0;
  for( size_t k {0u}; k < XY0.first[0].pat.num_bits(); ++k )
  {
    kitty::partial_truth_table pattern;
    for( size_t j{0}; j < XY0.first.size(); ++j )
      pattern.add_bit( kitty::get_bit(XY0.first[j].pat,k) );
    std::string s;
    s = kitty::to_binary( pattern );
    if( str_nodes0.find(s) == str_nodes0.end() )
    {
      N0++;
      str_nodes0.insert( std::make_pair(s, kitty::get_bit(XY0.second[0].pat,k)));
    }
    else if ( str_nodes0.at(s) != kitty::get_bit(XY0.second[0].pat, k) )
    {
      return false;
    }
  }

  uint32_t N1 = 0;
  for ( size_t k {0u}; k < XY1.first[0].pat.num_bits(); ++k )
  {
    kitty::partial_truth_table pattern;
    for( size_t j{0}; j < XY1.first.size(); ++j )
      pattern.add_bit( kitty::get_bit(XY1.first[j].pat, k) );
    std::string s;
    s = kitty::to_binary( pattern );
    if( already.find(s) == already.end() )
      N1++;

    if( str_nodes0.find(s) != str_nodes0.end() )
    {   
      if( str_nodes0.at(s) == kitty::get_bit(XY1.second[0].pat,k) )
        return false;
      else
      {
        if(already.find(s) == already.end())
          count_neg++;
      }
    }
    if(already.find(s)==already.end())
      already.insert(std::make_pair(s, kitty::get_bit( XY1.second[0].pat, k )));
  }
  
  size_t n = XY0.first.size()+1;
  std::pair<double,double> R = M1M2k(N0,N1,n);
  int Min = std::max(1,(int)floor(R.first-R.second));
  int Max = std::max(1,(int)ceil(R.first+R.second));

  bool is_satisfied = ( ( CumSum( count_neg+(int)ceil(R.second), N0, N1, n ) >= 1-0.001 ) && ( count_neg >=2 ) );

  if( is_satisfied )
    return true;

  return false;
  }  
#pragma endregion is_xor_decomposable

#pragma region is_top_decomposable
template<typename Ntk>
  sim_top_decomposition is_top_decomposable( std::pair<std::vector<sim_pattern<Ntk>>, sim_pattern<Ntk>> const& XY0, std::pair<std::vector<sim_pattern<Ntk>>, sim_pattern<Ntk>> const& XY1 )
  {
    if( (XY0.second.pat.num_bits() == 0) || (kitty::count_ones(XY0.second.pat) == 0) ) // F0 = 0
      return sim_top_decomposition::and_ ;
    else if( (XY1.second.pat.num_bits() != 0) && (kitty::count_ones(XY1.second.pat) == XY1.second.pat.num_bits()) ) // F1 = 1
      return sim_top_decomposition::or_ ;
    else if( (XY1.second.pat.num_bits() == 0) || (kitty::count_ones(XY1.second.pat) == 0) ) // F1 = 0
      return sim_top_decomposition::lt_ ;
    else if( (XY0.second.pat.num_bits() != 0) && (kitty::count_ones(XY0.second.pat) == XY0.second.pat.num_bits()) ) // F0 = 1
      return sim_top_decomposition::le_ ;
    else if( is_xor_decomposable( XY0, XY1 ) )
      return sim_top_decomposition::xor_ ;
    else
      return sim_top_decomposition::none ;
  }

  template<typename Ntk>
  sim_top_decomposition is_top_decomposable( std::pair<std::vector<sim_pattern<Ntk>>, std::vector<sim_pattern<Ntk>>> const& XY0, std::pair<std::vector<sim_pattern<Ntk>>, std::vector<sim_pattern<Ntk>>> const& XY1 )
  {
    if( (XY0.second[0].pat.num_bits() == 0) || (kitty::count_ones(XY0.second[0].pat) == 0) ) // F0 = 0
      return sim_top_decomposition::and_ ;
    else if( (XY1.second[0].pat.num_bits() != 0) && (kitty::count_ones(XY1.second[0].pat) == XY1.second[0].pat.num_bits()) ) // F1 = 1
      return sim_top_decomposition::or_ ;
    else if( (XY1.second[0].pat.num_bits() == 0) || (kitty::count_ones(XY1.second[0].pat) == 0) ) // F1 = 0
      return sim_top_decomposition::lt_ ;
    else if( (XY0.second[0].pat.num_bits() != 0) && (kitty::count_ones(XY0.second[0].pat) == XY0.second[0].pat.num_bits()) ) // F0 = 1
      return sim_top_decomposition::le_ ;
    else if( is_xor_decomposable( XY0, XY1 ) )
      return sim_top_decomposition::xor_ ;
    else
      return sim_top_decomposition::none ;
  }
#pragma endregion is_top_decomposable

#pragma region is_bottom_decomposable

 bool filter4_chatterjee_result( std::string& tt )
  {
    assert( tt.length() == 4 );
    if( ( tt == "1111") || ( tt == "0000") || ( tt == "1100") || ( tt == "0011") || ( tt == "1010") || ( tt == "0101") )
      return false;
    else
      return true;
  }

template<typename Ntk>
struct bottom_res
{
  bool found{false};
  chj_result chj;
  std::vector<signal<Ntk>> children;
};

template<typename Ntk>
bottom_res<Ntk> is_bottom_decomposable(  std::vector<sim_pattern<Ntk>>& X, sim_pattern<Ntk> & Y, double Imax,
                                      std::vector<double>& mi_vect, std::vector<uint32_t>& idx_vect )
{
  quicksort_by_attribute( idx_vect, mi_vect );

  std::vector<size_t> indeces2;

  bottom_res<Ntk> res;

  double Isupp, Ifnew, Ifr, Ifc, Ifrc;  
  chj_result chj_res;

  for( size_t i=0;  i < idx_vect.size()-1 ; ++i )
  {
    size_t r = idx_vect[i];

    double Ir = mi_vect[i];

    for( size_t j=i+1;  j < idx_vect.size() ; ++j )
    {
      if( mi_vect[i] != mi_vect[j] )
        break;

      size_t c = idx_vect[j];
      indeces2 = { r, c };
      std::vector<kitty::partial_truth_table * > Xtmp = { &(X[r].pat), &(X[c].pat) };
      chj_res = chatterjee_method( Xtmp, &(Y.pat) );
      
      Isupp = kitty::mutual_information( Xtmp, &(Y.pat) ); 
      Ifnew = kitty::mutual_information( std::vector{&chj_res.pat}, &(Y.pat) );
      std::vector<kitty::partial_truth_table *> V2 = { &(chj_res.pat), &(X[r].pat) };
      Ifr = kitty::mutual_information( V2, &(Y.pat) );
      std::vector<kitty::partial_truth_table *> V3 = { &(chj_res.pat), &(X[c].pat) };
      Ifc = kitty::mutual_information( V3, &(Y.pat) );
      std::vector<kitty::partial_truth_table *> V4 = { &(chj_res.pat), &(X[r].pat), &(X[c].pat) };
      Ifrc = kitty::mutual_information( V4, &(Y.pat) );

      bool is_informative = filter4_chatterjee_result( chj_res.tt );
      bool is_better = ( Ifnew > mi_vect[i] ) && ( Ifnew > mi_vect[j] ) && ( Ifnew >= Imax );
      bool filter = is_informative && is_better;

      if( filter && ( Isupp == Ifnew ) && (Ifrc == Ifnew) && ( Ifr == Ifnew ) && ( Ifc == Ifnew) )
      {
        res.found = true;
        res.chj = chj_res;
        res.children = { X[r].sig, X[c].sig };
        X.erase(X.begin()+std::max( r, c ) );
        X.erase(X.begin()+std::min( r, c ) );
        return res;
      }
    }
  }
  return res;
}
#pragma endregion is_bottom_decomposable

} // namespace mockturtle