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
  \file sim_create_nodes.hpp
  \brief Statistically optimal truth-table learning from examples

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

namespace mockturtle
{
  template< typename TT >
  struct sim_create_nodes_result
  {
    std::vector<std::string> tt_v;
    std::vector<TT> pat_v;
    std::vector<kitty::dynamic_truth_table> dtt_v;
  };

  namespace detail
  {

    template<typename TT>
    struct sim_create_nodes_impl
    {
      public:
        sim_create_nodes_impl( std::vector<TT * >& X, TT * Y ) 
          : X(X), 
            Y(Y) 
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

      sim_create_nodes_result<TT> run()
      {
        sim_create_nodes_result<TT> chj_res;

        size_t num_vars = X.size(); // number of variables
        size_t num_patterns = pow( 2, num_vars ); // number of patterns to consider
        size_t nconstr = std::is_same<TT, kitty::dynamic_truth_table>::value ? (int)log2( (*X[0]).num_bits() ) : (*X[0]).num_bits();
        TT signal0( nconstr ); // 0 for all examples
        TT mask_examples;  // mask selecting the examples for which the pattern under analysis is present
        std::vector<TT> new_values_v = { signal0 }; // new values at each example

        std::vector<std::string> tt_v = {""};
        std::default_random_engine generator;
        std::bernoulli_distribution distribution(0.5);
        uint64_t C0, C1;
        for( uint64_t k {0u}; k < num_patterns; ++k )
        {
          mask_examples = ~signal0;
          kitty::partial_truth_table mask_pattern( num_vars );
          std::string Sk = dec2bin( k, num_vars );

          kitty::create_from_binary_string( mask_pattern, Sk );

          for( uint64_t j {0u}; j < num_vars; ++j )
          {
            if( kitty::get_bit( mask_pattern, j) == 1 )
              mask_examples &= *X[j];
            else if( kitty::get_bit(mask_pattern, j ) == 0 )
              mask_examples &= ~*X[j];
            else
              std::cerr << "invalid" << std::endl;
        }
        C1 = kitty::count_ones(( mask_examples & *Y ));
        C0 = kitty::count_ones(( mask_examples & ~*Y ));
        //auto r = distribution(generator);

        if( ( C0 == 0 ) && ( C1 != 0 ) )
        {
          for( size_t j = 0; j < new_values_v.size(); ++j )
          {
            new_values_v[j] |= mask_examples;
            tt_v[j] = "1"+tt_v[j];
          }
        }
        else if( ( C1 == 0 ) && ( C0 != 0 ) )
        {
          for( size_t j = 0; j < new_values_v.size(); ++j )
            tt_v[j] = "0"+tt_v[j];
        }
        else
        {
          size_t M = new_values_v.size();
          for( size_t j = 0; j < M; ++j )
          {
            new_values_v.push_back(new_values_v[j]);
            new_values_v[j] |= mask_examples;
            tt_v.push_back(tt_v[j]);
            tt_v[j] = "1"+tt_v[j];
            tt_v[ tt_v.size() - 1 ] = "0"+tt_v[ tt_v.size() - 1 ];
          }
        }
      }
      
      chj_res.tt_v = tt_v;
      chj_res.pat_v = new_values_v;

      kitty::dynamic_truth_table dtt(X.size());
      for( size_t i = 0; i < tt_v.size(); ++i )
      {
        kitty::create_from_binary_string( dtt, tt_v[i] );
        chj_res.dtt_v.push_back( dtt );
      }

      return chj_res;
    }

    private:
      std::vector<TT * >& X;
      TT * Y;
    };

  } /* namespace detail */

  template<typename TT>
  sim_create_nodes_result<TT> sim_create_nodes_method( std::vector<TT * >& X, TT * Y )
  {
    return detail::sim_create_nodes_impl( X, Y ).run();
  }

  template<typename TT>
  sim_create_nodes_result<TT> sim_create_nodes_method( std::vector<TT * >& X, std::vector<TT * >& Y, size_t oidx )
  {
    return sim_create_nodes_method( X, Y[oidx] );
  }

  } /* namespace mockturtle */