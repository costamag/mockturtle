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
  \file node_creation.hpp
  \brief Statistics-based truth-table learning from examples
  \author Andrea Costamagna
*/

#pragma once

#include "../../traits.hpp"
#include <kitty/properties.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/partial_truth_table.hpp>

#include <vector>
#include <random>


namespace mockturtle
{
  std::string decimal_to_binary( uint32_t dec_number, uint32_t num_bits )
  {
    std::vector<uint32_t> bin_number;
    std::string bin_string = "";

    uint32_t i = 0;
    while (dec_number > 0) 
    {
      bin_number.push_back( dec_number % 2 );
      bin_string = ( dec_number % 2 == 1 ? "1" + bin_string : "0" + bin_string );
      dec_number = dec_number / 2;
      i++;
    }

    if( bin_number.size() < num_bits )
    { 
      for( auto j{0u}; j < num_bits - bin_number.size(); ++j )
      {
        bin_string = "0"+bin_string;
      }
    }
    return bin_string;
  }

  template<typename TT>
  struct chatterjee_result
  {
    std::string tt;
    TT pat;
  };

  /*! \brief Parameters for Chatterjee's node function creation. */
  struct chatterjee_method_params
  {
   /*! \brief Make non-trivial if this value is false. */
    bool detrivialize{true};
   /*! \brief Seed for the Bernoulli sampling. */
    uint32_t seed{123};
  };

  template<typename TT>
  struct nodes_enumeration_result
  {
    std::vector<std::string> tt_v;
    std::vector<TT> pat_v;
    std::vector<uint32_t> degree_correlation_v;
  };

  /*! \brief Parameters for Chatterjee's node function creation. */
  struct nodes_enumeration_params
  {
   /*! \brief Make non-trivial if this value is false. */
    bool detrivialize{true};
   /*! \brief Seed for the Bernoulli sampling. */
    uint32_t seed{123};
  };

  namespace detail
  {

    template<typename TT>
    struct chatterjee_method_impl
    {
      public:
        chatterjee_method_impl( std::vector<TT*>& X, TT * Y, chatterjee_method_params & ps = {} ) 
          : X(X), Y(Y), ps(ps)
        {
        }

      public:
        chatterjee_result<TT> run()
        {
          chatterjee_result<TT> chj_res;

          uint32_t num_vars = X.size();
          uint32_t num_patterns = pow( 2, num_vars );
          
          TT const0 = (*X[0]).construct();
          TT mask_examples = (*X[0]).construct();
          TT new_values = (*X[0]).construct();
          std::string tt = "";

          uint32_t C0, C1;

          for( uint32_t k {0u}; k < num_patterns; ++k )
          {
            mask_examples = ~const0;
            kitty::partial_truth_table mask_pattern( num_vars );

            std::string binary_k = decimal_to_binary( k, num_vars );
            kitty::create_from_binary_string( mask_pattern, binary_k );

            for( uint32_t j {0u}; j < num_vars; ++j )
            {
              if( kitty::get_bit( mask_pattern, j) == 1 )
                mask_examples &= *X[j];
              else if( kitty::get_bit(mask_pattern, j ) == 0 )
                mask_examples &= ~(*X[j]);
            }

            C1 = kitty::count_ones(( mask_examples & *Y ));
            C0 = kitty::count_ones(( mask_examples & ~(*Y) ));

            srand(ps.seed++);
            double r = rand() % 2;

            if( ( C1 > C0 ) || ( ( C1 == C0 ) && ( r >= 0.5 ) ) )
            {
              new_values |= mask_examples;
              tt = "1"+tt;
            }
            else if( ( C1 < C0 ) || ( ( C1 == C0 ) && ( r < 0.5 ) ) )
            {
              tt = "0"+tt;
            }
          }

          kitty::dynamic_truth_table dtt_attempt(X.size());
          kitty::create_from_binary_string( dtt_attempt, tt );

          chj_res.tt = tt;
          chj_res.pat = new_values;  
        
          return chj_res;
        }

      private:
        std::vector<TT*> & X;
        TT * Y;
        chatterjee_method_params & ps;
    };

    template<typename TT>
    struct nodes_enumeration_impl
    {
      public:
        nodes_enumeration_impl( std::vector<TT*>& X, TT * Y, nodes_enumeration_params & ps = {} ) 
          : X(X), Y(Y), ps(ps)
        {
        }

      public:
        nodes_enumeration_result<TT> run()
        {
          nodes_enumeration_result<TT> nen_res;

          uint32_t num_vars = X.size();
          uint32_t num_patterns = pow( 2, num_vars );
          
          TT const0 = (*X[0]).construct();
          TT mask_examples = (*X[0]).construct();
          std::vector<TT> new_values_v = std::vector{(*X[0]).construct()};
          std::vector<std::string> tt_v = {""};

          uint32_t C0, C1;
          std::vector<uint32_t> degree_correlation_v = {0};

          for( uint64_t k {0u}; k < num_patterns; ++k )
          {
            mask_examples = ~const0;
            kitty::partial_truth_table mask_pattern( num_vars );

            std::string binary_k = decimal_to_binary( k, num_vars );
            kitty::create_from_binary_string( mask_pattern, binary_k );

            for( uint64_t j {0u}; j < num_vars; ++j )
            {
              if( kitty::get_bit( mask_pattern, j) == 1 )
                mask_examples &= *X[j];
              else if( kitty::get_bit(mask_pattern, j ) == 0 )
                mask_examples &= ~*X[j];
            }
            
            C1 = kitty::count_ones(( mask_examples & *Y ));
            C0 = kitty::count_ones(( mask_examples & ~*Y ));

            if( ( C0 == 0 ) && ( C1 != 0 ) )
            {
              for( size_t j = 0; j < new_values_v.size(); ++j )
              {
                new_values_v[j] |= mask_examples;
                tt_v[j] = "1"+tt_v[j];
                degree_correlation_v[j] += C1;
              }
            }
            else if( ( C1 == 0 ) && ( C0 != 0 ) )
            {
              for( size_t j = 0; j < new_values_v.size(); ++j )
              {
                tt_v[j] = "0"+tt_v[j];
                degree_correlation_v[j] += C0;
              }
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
                degree_correlation_v.push_back(degree_correlation_v[j]);
                degree_correlation_v[j] += C1;
                degree_correlation_v[degree_correlation_v.size() - 1 ] += C0;
              }
            }
          }
      
          kitty::dynamic_truth_table dtt(X.size());
          for( size_t i = 0; i < tt_v.size(); ++i )
          {
            kitty::create_from_binary_string( dtt, tt_v[i] );

            if( !kitty::is_trivial( dtt ) )
            {
              if( nen_res.degree_correlation_v.size()==0 || nen_res.degree_correlation_v[nen_res.degree_correlation_v.size()-1] >= degree_correlation_v[i] )
              {
                nen_res.tt_v.push_back( tt_v[i] );
                nen_res.pat_v.push_back( new_values_v[i] );
                nen_res.degree_correlation_v.push_back( degree_correlation_v[i] );
              }
              else
              {
                for( uint32_t j = 0; j < nen_res.degree_correlation_v.size(); ++j )
                {
                  if( degree_correlation_v[i] > nen_res.degree_correlation_v[j] )
                  {
                    nen_res.tt_v.insert( nen_res.tt_v.begin()+j, tt_v[i] );
                    nen_res.pat_v.insert( nen_res.pat_v.begin()+j, new_values_v[i] );
                    nen_res.degree_correlation_v.insert( nen_res.degree_correlation_v.begin()+j, degree_correlation_v[i] );
                    break;
                  }
                }
              }
            }
          } 
        
          return nen_res;
        }

      private:
        std::vector<TT*> & X;
        TT * Y;
        nodes_enumeration_params & ps;
    };

  } /* namespace detail */

  template<typename TT>
  chatterjee_result<TT> chatterjee_method( std::vector<TT*>& X, TT * Y, chatterjee_method_params ps = {} )
  {
    return detail::chatterjee_method_impl<TT>( X, Y, ps ).run();
  }

  template<typename TT>
  chatterjee_result<TT> chatterjee_method( std::vector<TT*>& X, std::vector<TT*>& Y, uint32_t oidx, chatterjee_method_params ps = {} )
  {
    return chatterjee_method( X, Y[oidx], ps );
  }

  template<typename TT>
  nodes_enumeration_result<TT> nodes_enumeration( std::vector<TT*>& X, TT * Y, nodes_enumeration_params ps = {} )
  {
    return detail::nodes_enumeration_impl<TT>( X, Y, ps ).run();
  }

  template<typename TT>
  nodes_enumeration_result<TT> nodes_enumeration( std::vector<TT*>& X, std::vector<TT*>& Y, uint32_t oidx, nodes_enumeration_params ps = {} )
  {
    return nodes_enumeration( X, Y[oidx], ps );
  }

  } /* namespace mockturtle */