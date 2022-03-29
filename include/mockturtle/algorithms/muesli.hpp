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
  \file dsd_decomposition.hpp
  \brief DSD decomposition

  \author Andrea Costamagna
*/

#pragma once

#include <cstdint>
#include <vector>

#include "../traits.hpp"

#include <fmt/format.h>
#include <kitty/constructors.hpp>
#include <kitty/decomposition.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/dbs_statistical_bit_operations.hpp>
#include <kitty/operators.hpp>
#include <kitty/print.hpp>
#include <bitset>


#include <mockturtle/algorithms/graph_to_lfe.hpp>


namespace mockturtle
{

/*! \brief Parameters for dsd_decomposition */
struct muesli_params
{
  /* input parameters */
  size_t max_act{10};
  size_t max_sup{3};
  double eps_th{1};
  size_t init_sup{2};
  bool is_po_only_if_exact{true};
  /* output parameters */
  bool is_exact_fn{false};
};

using dbitset = boost::dynamic_bitset<>;
using dbitset_vector = std::vector<dbitset>;

namespace detail
{

class muesli_impl
{

public:
  muesli_impl( klut_network& ntk, lfeNtk<klut_network>& examples, muesli_params& ps )
      : _klut( ntk ),
        _examples( examples ),
        _ps( ps )
  {
    _num_out = examples.complete.second.size();
  }

  void remove_column_and_invert( dbitset_vector& X, dbitset& Y, uint64_t idx )
  {
    Y ^= X[idx];
    X.erase( X.begin() + idx );
  }

  #pragma region storage units
    struct Istorage
  {
    std::unordered_map<std::string, double> Fnew;
    std::unordered_map<std::string, double> Fr;
    std::unordered_map<std::string, double> Fc;
    std::unordered_map<std::string, double> Frc;
    std::unordered_map<std::string, double> supp;

    void clear()
    {
      Fnew.clear();
      Fr.clear();
      Fc.clear();
      Frc.clear();
      supp.clear();
    }
  };
  struct muesli_vars
  {
    size_t act{0};
    size_t sup{2}; 
  };

  #pragma endregion sotorage units

  #pragma region quicksort
    template<typename T>
    void swap( T& a, T& b)
    {
      T t = a;
      a = b;
      b = t;
    }

    int partition ( std::vector<uint64_t>& support, std::vector<double>& attribute , uint64_t low, uint64_t high )
    {
    double pivot = attribute[high];    // pivot
    int i = (low-1);  // Index of smaller element
 
    for (int j = low; j < high; j++)
      {
        if ( attribute[j] >= pivot)
        {
            i++;    // increment index of smaller element
            swap(attribute[i], attribute[j]);
            swap(support[i], support[j]);
        }
      }
      swap(attribute[i + 1], attribute[high]);
      swap(support[i + 1], support[high]);
      return (i + 1);
    }

    void quicksort_by_attribute( std::vector<uint64_t>& support, std::vector<double>& attribute,  int low, int high )
    {
      if (low > high)
      {
        auto pi = partition( support, attribute, low, high);
        quicksort_by_attribute( support, attribute, low, pi - 1);
        quicksort_by_attribute( support, attribute, pi + 1, high);
      }
    }
    #pragma endregion

  #pragma region chatterjee
  struct chj_info
  {
    bool is_exact{false};
    std::string tt_str;
    dbitset x;
  };

  chj_info chatterjee_method( dbitset_vector const& X, dbitset const& Y )
  {
    chj_info ninfo;
    bool only_zeros = true;
    ninfo.is_exact = true;
    uint64_t N = X.size();
    uint64_t pow2N = pow(2,N);
    dbitset bit02N( X[0].size(), 0 );
    dbitset Kmask, new_values;
    new_values = bit02N;

    ninfo.tt_str = "";
    std::default_random_engine generator;
    std::bernoulli_distribution distribution(0.5);
    uint64_t C0, C1;
    std::vector<uint64_t> C1s;
    std::vector<uint64_t> C1s_r;
    std::vector<size_t> random_indeces;
    for( uint64_t k {0u}; k < pow2N; ++k )
    {
      Kmask = ~bit02N;
      dbitset maskN( N, k );
      for( uint64_t j {0u}; j < N; ++j )
      {

        if( maskN[j] == 1 )
          Kmask &= X[j];
        else if( maskN[j] == 0 )
          Kmask &= ~X[j];
        else
          std::cerr << "invalid" << std::endl;
      }

      C1 = ( Kmask & Y ).count();
      C0 = ( Kmask & ~Y ).count();
      auto r = distribution( generator );


      if( ( C1 > C0 ) || ( ( C1 == C0 ) && ( r >= 0.5 ) ) )
      {
        only_zeros = false;
        new_values |= Kmask;
        ninfo.tt_str = "1"+ninfo.tt_str;
      }
      else if( ( C1 < C0 ) || ( ( C1 == C0 ) && ( r < 0.5 ) ) )
      {
        ninfo.tt_str = "0"+ninfo.tt_str;
      }
      if( C1 != 0 && C0 != 0 )
        ninfo.is_exact &= false;
      if( C1 == C0 )
      {
        random_indeces.push_back( k );
        C1s.push_back(C1);
        C1s_r.push_back(C1);
      }
      else if( C1 > C0 )
      {
        C1s.push_back(C1);
        C1s_r.push_back(0);
      }
      else
      {
        C1s.push_back(0);
        C1s_r.push_back(0);
      }
    }
    ninfo.x = new_values;
    // don't want contraddictions to appear
    if( only_zeros && ( random_indeces.size() != 0 ) )
    {
      ninfo.tt_str = "";
      uint64_t maxC1_r = 0;
      uint64_t idx_max = 0;
      for( auto i = 0; i < C1s_r.size(); ++i )
        idx_max = ( C1s_r[i] > maxC1_r ) ? i : idx_max;

      for( auto i = 0; i < C1s_r.size(); ++i )
          ninfo.tt_str = ( i == idx_max ) ? "1"+ninfo.tt_str : "0"+ninfo.tt_str ;
      
      dbitset bit1( C1s.size(), 1u );
      std::cout << bit1 <<" " <<  ( bit1 << idx_max ) <<std::endl;
      
      Kmask = ~bit02N;
      dbitset maskN( N, idx_max );
      for( uint64_t j {0u}; j < N; ++j )
      {
        if( maskN[j] == 1 )
          Kmask &= X[j];
        else if( maskN[j] == 0 )
          Kmask &= ~X[j];
        else
          std::cerr << "invalid" << std::endl;
      }
      ninfo.x |= Kmask;
    }

    return ninfo;
  }
  #pragma endregion chatterjee

  #pragma region active list
  
  struct alist_info
  {
    dbitset_vector X;
    std::vector<size_t> indeces;
    std::vector<signal<klut_network>> support;
  };

  alist_info fill_active_list( std::vector<signal<klut_network>> const& support, dbitset_vector const& X, dbitset const& Y, int sizeA = -1 )
  {
    alist_info ninfo;
    std::vector<size_t> indeces_orig;
    for( size_t i = 0; i < X.size(); ++i )
      indeces_orig.push_back( i );

    if( sizeA < 0 )
      sizeA = ( sizeA < 1 ) ? _ps.max_act : sizeA;
    double Imax = 0;
    double Inew = 0;
    size_t active = 0;
    
    for( size_t i=0; i < X.size(); ++i )
    {
      Inew = kitty::mutual_information( X[i], Y );
      if( Inew > Imax )
      {
        active = i;
        Imax = Inew;
      }
    }
    ninfo.indeces.push_back( active );      
    ninfo.X.push_back( X[active] );
    ninfo.support.push_back( support[active] );

    if( sizeA > 1 )
    {
      for( size_t i=1; i<std::min( sizeA, (int)X.size() ); ++i )
      {
        active = 0;
        Imax = 0;
        Inew = 0;
        for( size_t j=0; j<X.size(); ++j )
          {
            if ( std::find(ninfo.indeces.begin(), ninfo.indeces.end(), j) == ninfo.indeces.end() )
            {
              ninfo.X.push_back(X[j]);
              Inew = kitty::mutual_information( ninfo.X, Y );
              if( Inew > Imax )
              {
                active = j;
                Imax = Inew;
              }
              ninfo.X.erase( ninfo.X.begin() + ninfo.X.size() - 1 );
            }
          }
          ninfo.X.push_back( X[active] );
          ninfo.indeces.push_back( active );
          ninfo.support.push_back( support[active] );

        }
      }
    return ninfo;  
  }
  #pragma endregion active list

  #pragma region not done
  struct nd_info
  {
    bool truth;
    signal<klut_network> sig;
  };

  nd_info not_done( std::vector<signal<klut_network>> const& support, dbitset_vector const& X, dbitset const& Y )
  {
    double Inew = 0;
    double Imax = 0;
    double Hy = kitty::entropy( Y );
    double best_eps = 0;
    size_t best_idx = 0;
    for( size_t i=0; i< X.size(); ++i )
    {
      Inew = kitty::mutual_information( X[i], Y );
      if( Inew > Imax )
      {
        best_idx = i;
        best_eps = Inew / Hy;
        Imax = Inew;
      }
    }
    nd_info ninfo;

    if( _vars.sup > _ps.max_sup )
    {
      ninfo.truth = false;
      return ninfo;
    }

    ninfo.sig = support[best_idx];
    
    ninfo.truth = best_eps < _ps.eps_th;
    if( ninfo.truth )
    {
      is_exact_functionality = true;
      _ps.is_exact_fn = true;
    }
    return ninfo;
  }
#pragma endregion not done

#pragma region improve function
bool is_not_trivial( std::string tt_str )
{
  bool ans = true;
  ans &= ( tt_str != "1111" ); 
  ans &= ( tt_str != "0000" ); 
  ans &= ( tt_str != "1100" ); 
  ans &= ( tt_str != "0011" ); 
  ans &= ( tt_str != "1010" ); 
  ans &= ( tt_str != "0101" );
  return ans; 
}

struct improve_info
{
  std::string tt_str;
  bool success;
  dbitset x;
  std::vector<signal<klut_network>> support;
};
  improve_info improve_fn( std::vector<signal<klut_network>> const& support, dbitset_vector const & X, dbitset const& Y )
    {
      improve_info ninfo;
      auto A_info = fill_active_list( support, X, Y );

      if( _vars.act+_vars.sup > A_info.support.size() )
      {
        ninfo.success = false;
        return ninfo;
      }
      
      std::vector<signal<klut_network>> new_support;
      dbitset_vector Xnew;

      for ( size_t k{0}; k < _vars.sup; ++k )
      {
        Xnew.push_back( A_info.X[ _vars.act + k ]);
        new_support.push_back( A_info.support[ _vars.act + k ] );
        //std::cout << A_info.X[ _vars.act + k ] << std::endl;

      }

      dbitset_vector first_act;
      for (uint64_t k {0u}; k <= _vars.act; ++k )
      {
        first_act.push_back(A_info.X[k]);
      }
      
      
      auto mi_old = kitty::mutual_information( first_act, Y ); //#########################################################################
      
      auto chj_info = chatterjee_method( Xnew, Y );

      //std::cout << "> " << chj_info.x << std::endl;
      std::cout << chj_info.tt_str << std::endl;

      ninfo.tt_str = chj_info.tt_str;
      ninfo.support = new_support;
      first_act[ first_act.size() - 1 ] = chj_info.x;

      auto mi_new = kitty::mutual_information( first_act , Y );
      std::cout << mi_new << " > "<< mi_old << " ? " << std::endl; 
      if (  mi_new > mi_old && is_not_trivial( ninfo.tt_str ) )
      {
        ninfo.x = chj_info.x;
        ninfo.success = true;
        return ninfo;

      }
      ninfo.success = false;
      return ninfo;
    }
  #pragma endregion improve function

#pragma region muesli
  signal<klut_network> muesli_step( std::vector<signal<klut_network>>& support, dbitset_vector & X, dbitset & Y )
  {
    bool success; /* true if found a function improving the mi */
    std::string tt_str; /* contains the tt of the new node */
    size_t best_idx; /* used to keep track of the new best approximation */
    double eps_best = 0;
    double Inew;
    double Imax = 0;
    auto support_active = support;
    auto X_active = X;
    auto support_copy = support;
    auto X_copy = X;

    _vars.sup = _ps.init_sup; // sup <- 2
    while( not_done( support, X, Y ).truth && ( _vars.sup <= _ps.max_sup ) ) 
    {

      success = false;
      dbitset_vector X_new;
      std::vector<size_t> support_new;
      _vars.act = 0;
      do {
          if( _vars.sup > _ps.max_sup )
            break;
          auto impr_info = improve_fn( support, X, Y ); // improve mi
          success = impr_info.success;
          
          if ( success )
          {
            X.push_back( impr_info.x );
            for( auto s : impr_info.support )
              std::cout << s << std::endl;
            //std::cout << impr_info.x << std::endl;

            kitty::dynamic_truth_table tt( impr_info.support.size() );
            kitty::create_from_binary_string( tt, impr_info.tt_str );
            support.push_back( _klut.create_node( impr_info.support, tt ) );
          }
          if ( not_done( support, X, Y ).truth == false )
            break;
          _vars.act++;
        } while( ( success == false ) && ( _vars.act <= _ps.max_act ) );
        if ( success == true )
        {
          if ( not_done( support, X, Y ).truth == false )
            break;
          
          _vars.sup = _ps.init_sup;

          while( success == true )
          {
            if( _vars.sup > _ps.max_sup )
              break;
            auto impr_info = improve_fn( support, X, Y );            
            success = impr_info.success;

            if ( success )
            {
              for( auto s : impr_info.support )
                std::cout << s << std::endl;
              //std::cout << impr_info.x << std::endl;
              
              X.push_back( impr_info.x );
              kitty::dynamic_truth_table tt( impr_info.support.size() );
              kitty::create_from_binary_string( tt, impr_info.tt_str );
              support.push_back( _klut.create_node( impr_info.support, tt ) );
            }
          }
        }
        else
        {
          _vars.sup++;
          if( _vars.sup > _ps.max_sup )
            break;
        }
      }

      Inew = 0;
      Imax = 0;
      best_idx = 0;
      for( size_t i=0; i< X.size(); ++i )
      {
        Inew = kitty::mutual_information( X[i], Y );
        if( Inew > Imax )
          best_idx = i;
      }
      return support[best_idx];

  }
  #pragma endregion muesli

  signal<klut_network> run()
  {
    return muesli_step( _examples.signals, _examples.partial.first, _examples.partial.second );
  }
private:
  klut_network& _klut;
  lfeNtk<klut_network> _examples;
  signal<klut_network> _support;
  uint64_t _num_out;
  Istorage _Icoll;
  muesli_params _ps;
  muesli_vars _vars;
public:
  bool is_exact_functionality = false;
};
} // namespace detail

#pragma region iwls2020

klut_network muesli( dbitset_vector& X, dbitset& Y, muesli_params& ps )
{
  lfeNtk<klut_network> examples;
  examples.partial.first = X;
  examples.partial.second = Y;
  klut_network klut;
  for( auto x : X )
    examples.signals.push_back( klut.create_pi() );
  detail::muesli_impl impl( klut, examples, ps );

  auto osignal = impl.run();
  if( (ps.is_po_only_if_exact && impl.is_exact_functionality) || !ps.is_po_only_if_exact )
  {
    klut.create_po( osignal );
  }

  return klut;
}
#pragma endregion iwls2020
} // namespace mockturtle