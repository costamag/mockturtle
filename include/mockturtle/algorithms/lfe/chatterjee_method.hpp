/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2021  EPFL
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
  \file chatterjee_method.hpp
  \brief Statistically optimal truth-table learning from examples

  \author Andrea Costamagna
*/
#pragma once


#include "../../traits.hpp"

#include <kitty/partial_truth_table.hpp>
#include <kitty/operators.hpp>
#include <kitty/constructors.hpp>
#include <vector>

namespace mockturtle
{

  struct chj_result
  {
    std::string tt;
    bool both_not0_and_eq {false};
    bool both_not0 {false};
  };

  namespace detail
  {
    using TT = kitty::partial_truth_table;

    std::string decToBinary( int n, size_t N )
    {
      // array to store binary number
      std::vector<int> binaryNum;
      std::string S = "";
 
      // counter for binary array
      int i = 0;
      while (n > 0) {
        // storing remainder in binary array
        binaryNum.push_back( n % 2 );
        S = (n % 2 == 1 ? "1"+S : "0"+S );
        n = n / 2;
        i++;
      }
    
      if( binaryNum.size() < N )
      { 
        for( auto j{0u}; j < N - binaryNum.size(); ++j )
        {
          S = "0"+S;
        }
      }
      return S;
    }

    chj_result chatterjee_method_impl( std::vector<TT>& X, TT& Y, std::vector<size_t> indeces = {} )
    {
      chj_result chj_res;

      bool both_not0_and_eq = true; // all input patterns => same output
      bool both_not0 = true; // all input patterns => same output
      size_t num_vars = indeces.size(); // number of variables
      size_t num_patterns = pow( 2,num_vars ); // number of patterns to consider
      TT signal0( X[0].num_bits() ); // 0 for all examples
      TT mask_examples;  // mask selecting the examples for which the pattern under analysis is present
      TT new_values; // new values at each example
      new_values = signal0; // initialize to 0
      std::string tt = "";
      std::default_random_engine generator;
      std::bernoulli_distribution distribution(0.5);
      uint64_t C0, C1;
      for( uint64_t k {0u}; k < num_patterns; ++k )
      {
        mask_examples = ~signal0;
        TT mask_pattern( num_vars );
        std::string Sk = decToBinary( k, num_vars );

        kitty::create_from_binary_string( mask_pattern, Sk );

        for( uint64_t j {0u}; j < num_vars; ++j )
        {
          if( kitty::get_bit( mask_pattern, j) == 1 )
            mask_examples &= X[indeces[j]];
          else if( kitty::get_bit(mask_pattern, j ) == 0 )
            mask_examples &= ~X[indeces[j]];
          else
            std::cerr << "invalid" << std::endl;
        }
        C1 = kitty::count_ones(( mask_examples & Y ));
        C0 = kitty::count_ones(( mask_examples & ~Y ));
        auto r = distribution(generator);

        if( ( C1 > C0 ) || ( ( C1 == C0 ) && ( r >= 0.5 ) ) )
        {
          new_values |= mask_examples;
          tt = "1"+tt;
        }
        else if( ( C1 < C0 ) || ( ( C1 == C0 ) && ( r < 0.5 ) ) )
        {
          tt = "0"+tt;
        }
        if( (C1 != 0 && C0 != 0) )
          both_not0 &= false;
        if( (C1 != 0 && C0 != 0) && ( C0 == C1 ) )
          both_not0_and_eq &= false;
      }
      X.push_back( new_values );
      
      chj_res.tt = tt;
      chj_res.both_not0_and_eq = both_not0_and_eq ;
      chj_res.both_not0 = both_not0;
      
      return chj_res;
    }

  template<typename Ntk>
  signal<Ntk> apply_chatterjee_impl( Ntk& ntk, std::vector<signal<Ntk>>& support, std::vector<TT>& X, TT & Y, std::vector<uint64_t> indeces )
  {
    std::string tt_str = chatterjee_method_impl( X, Y, indeces ).tt;
    kitty::dynamic_truth_table tt( support.size() );
    kitty::create_from_binary_string( tt, tt_str );
    return ntk.create_node( support, tt );
  }

  } /* namespace detail */

  chj_result chatterjee_method( std::vector<kitty::partial_truth_table>& X, kitty::partial_truth_table& Y, std::vector<uint64_t> indeces = {} )
  {
    if( indeces.size() == 0 )
    {
      for( size_t i{0}; i < X.size(); ++i )
        indeces.push_back(i);
    }
    return detail::chatterjee_method_impl( X, Y, indeces );
  }

  chj_result chatterjee_method( std::vector<kitty::partial_truth_table>& X, std::vector<kitty::partial_truth_table>& Y, size_t oidx, std::vector<uint64_t> indeces = {} )
  {
    return chatterjee_method( X, Y[oidx], indeces );
  }

  template<typename Ntk>
  signal<Ntk> apply_chatterjee( Ntk& ntk, std::vector<signal<Ntk>>& support, std::vector<kitty::partial_truth_table>& X, kitty::partial_truth_table & Y, std::vector<uint64_t> indeces = {} )
  {
    static_assert( is_network_type_v<Ntk>, "NtkDest is not a network type" );
    static_assert( has_create_node_v<Ntk>, "NtkDest does not implement the create_node method" );
    if( indeces.size() == 0 )
    {
      assert( support.size() == X.size() );
      for( size_t i{0}; i < X.size(); ++i )
        indeces.push_back(i);
    }
    else
    {
      assert( indeces.size() == support.size() );
    }
    return detail::apply_chatterjee_impl( ntk, support, X, Y, indeces );
  }
  
  template<typename Ntk>
  signal<Ntk> apply_chatterjee( Ntk& ntk, std::vector<signal<Ntk>>& support, std::vector<kitty::partial_truth_table>& X, std::vector<kitty::partial_truth_table>& Y, size_t oidx, std::vector<uint64_t> indeces = {} )
  {
    return apply_chatterjee( ntk, support, X, Y[oidx], indeces );
  }

  } /* namespace mockturtle */