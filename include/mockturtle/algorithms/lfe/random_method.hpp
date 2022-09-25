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
  \file random_method.hpp
  \brief Bayes-optimal truth-table learning from examples
  \author Andrea Costamagna
*/

#pragma once

#include "../../traits.hpp"
#include <kitty/partial_truth_table.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operators.hpp>
#include <kitty/constructors.hpp>
#include <kitty/print.hpp>
#include <vector>
#include <random>
#include <kitty/properties.hpp>


namespace mockturtle
{

  struct chj_result
  {
    std::string tt;
    kitty::partial_truth_table pat;
    bool both_not0_and_eq {false};
    bool both_not0 {false};
    kitty::dynamic_truth_table dtt;
  };

  namespace detail
  {

    template<typename TT>
    struct random_method_impl
    {
      public:
        random_method_impl( std::vector<TT * >& X, TT * Y, int& seed = 123 ) 
          : X(X), 
            Y(Y),
            seed(seed)
        {
        }

      std::string dec2bin( int n, size_t N )
      {
        std::vector<int> binaryNum;
        std::string S = "";

        int i = 0;
        while (n > 0) 
        {
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

      chj_result run()
      {
        chj_result chj_res;

        bool both_not0_and_eq = true; // all input patterns => same output
        bool both_not0 = true; // all input patterns => same output
        size_t num_vars = X.size(); // number of variables
        size_t num_patterns = pow( 2, num_vars ); // number of patterns to consider
        size_t nconstr = std::is_same<TT, kitty::dynamic_truth_table>::value ? (int)log2( (*X[0]).num_bits() ) : (*X[0]).num_bits();
        TT signal0( nconstr ); // 0 for all examples
        TT mask_examples( nconstr );  // mask selecting the examples for which the pattern under analysis is present
        TT new_values; // new values at each example
        new_values = signal0; // initialize to 0
        std::string tt = "";

        uint32_t C0, C1;
        std::vector<uint32_t> flippable_bits;
        for( uint32_t k {0u}; k < num_patterns; ++k )
        {
          mask_examples = ~signal0;
          kitty::partial_truth_table mask_pattern( num_vars );

          std::string Sk = dec2bin( k, num_vars );

          kitty::create_from_binary_string( mask_pattern, Sk );


          for( uint32_t j {0u}; j < num_vars; ++j )
          {
            if( kitty::get_bit( mask_pattern, j) == 1 )
              mask_examples &= *X[j];
            else if( kitty::get_bit(mask_pattern, j ) == 0 )
              mask_examples &= ~(*X[j]);
            else
              std::cerr << "invalid" << std::endl;
          }
        C1 = kitty::count_ones(( mask_examples & *Y ));
        C0 = kitty::count_ones(( mask_examples & ~(*Y) ));
        srand(seed++);
        double r = rand()%2;
        if( C1 == C0 )
          flippable_bits.push_back( k );

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

      kitty::dynamic_truth_table dtt_attempt(X.size());
      kitty::create_from_binary_string( dtt_attempt, tt );

      if( !kitty::is_trivial( dtt_attempt ) && !kitty::is_trivial( ~dtt_attempt ) )
      {
        chj_res.tt = tt;
        chj_res.pat = new_values;
        chj_res.both_not0_and_eq = both_not0_and_eq ;
        chj_res.both_not0 = both_not0;
        chj_res.dtt = dtt_attempt;
      }
      else
      {
        uint32_t num_flips = std::ceil( (double)flippable_bits.size()/2 );
        while( flippable_bits.size() > num_flips )
        {
          auto eidx = rand()%( flippable_bits.size() );
          flippable_bits.erase( flippable_bits.begin() + eidx );
        }
        for( uint32_t i{0}; i < flippable_bits.size(); ++i )
          tt[flippable_bits[i]] = ( tt[flippable_bits[i]]==char('0') ) ? char('1') : char('0');

        new_values = signal0; // initialize to 0
        for( uint32_t k {0u}; k < num_patterns; ++k )
        {
          mask_examples = ~signal0;
          kitty::partial_truth_table mask_pattern( num_vars );

          std::string Sk = dec2bin( k, num_vars );

          kitty::create_from_binary_string( mask_pattern, Sk );


          for( uint32_t j {0u}; j < num_vars; ++j )
          {
            if( kitty::get_bit( mask_pattern, j) == 1 )
              mask_examples &= *X[j];
            else if( kitty::get_bit(mask_pattern, j ) == 0 )
              mask_examples &= ~(*X[j]);
            else
              std::cerr << "invalid" << std::endl;
          }
          if(tt[k]==char('1'))
            new_values |= mask_examples;
        }

        chj_res.tt = tt;
        chj_res.pat = new_values;
        chj_res.both_not0_and_eq = false ;
        chj_res.both_not0 = false;
        kitty::dynamic_truth_table dtt(X.size());
        kitty::create_from_binary_string( dtt, tt );
        chj_res.dtt = dtt_attempt;
      }      
      return chj_res;
    }

    private:
      std::vector<TT * >& X;
      TT * Y;
      int & seed;
    };

  } /* namespace detail */

  template<typename TT>
  chj_result random_method( std::vector<TT *>& X, TT * Y, int& seed )
  {
    return detail::random_method_impl<TT>( X, Y, seed ).run();
  }

  template<typename TT>
  chj_result random_method( std::vector<TT * >& X, std::vector<TT * >& Y, size_t oidx, int& seed )
  {
    return random_method( X, Y[oidx], seed );
  }

    template<typename TT>
  chj_result random_method( std::vector<TT *>& X, TT * Y, int seed = 123 )
  {
    return detail::random_method_impl<TT>( X, Y, seed ).run();
  }

  template<typename TT>
  chj_result random_method( std::vector<TT * >& X, std::vector<TT * >& Y, size_t oidx, int seed = 123 )
  {
    return random_method( X, Y[oidx], seed );
  }

  } /* namespace mockturtle */