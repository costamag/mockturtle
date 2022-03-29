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
#include <mockturtle/algorithms/muesli.hpp>

#include <kitty/operators.hpp>
#include <kitty/print.hpp>
#include <bitset>


#include <mockturtle/algorithms/graph_to_lfe.hpp>


namespace mockturtle
{

/*! \brief Parameters for dsd_decomposition */
struct it_decomposition_params
{
  /*! \brief Apply XOR decomposition. */
  bool is_informed{true};
  size_t max_sup{3};
  bool try_top_decomposition{true};
  bool try_bottom_decomposition{true};
  bool try_creation{false};
  bool try_xor_decomposition{false};
  bool use_cumsum{false};
  bool is_bottom_exact{false};
  bool is_trivial{true};
};

struct detection_counter
{
  size_t _and;
  size_t _or;
  size_t _lt;
  size_t _le;
  size_t _xor;
  size_t _btm;
  size_t _ctj;
  size_t _cre;

  detection_counter():
    _and(0),
    _or(0),
    _lt(0),
    _le(0),
    _xor(0),
    _btm(0),
    _ctj(0),
    _cre(0)
  {}

};

struct XYdataset{
  std::vector<boost::dynamic_bitset<>> X;
  boost::dynamic_bitset<> Y;
  uint64_t nin;
  uint64_t nout;
  uint64_t ndata;
};

using dbitset = boost::dynamic_bitset<>;
using dbitset_vector = std::vector<dbitset>;

namespace detail
{

enum class it_top_decomposition
{
  none,
  and_,
  or_,
  lt_,
  le_,
  xor_
};

class it_decomposition_impl
{

public:
  it_decomposition_impl( klut_network& ntk, lfeNtk<klut_network>& examples, it_decomposition_params& ps )
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

  #pragma endregion sotorage units

  #pragma region cofactors
  std::pair<dbitset_vector, dbitset> compute_cofactor( dbitset_vector const& X, dbitset const& Y, 
                                                              uint64_t const& id, uint64_t idx = 0 )
  {
    if( X.size()==0 )
      return std::make_pair(X,Y);
    assert( X[0].size() == Y.size() );
    assert( idx < X.size() );

    dbitset M;
    M = (id == 1) ?  X[idx] : ~X[idx];
    dbitset_vector Xid;
    dbitset Yid;

    typedef boost::dynamic_bitset<>::size_type size_type;
    const size_type npos = boost::dynamic_bitset<>::npos;
    size_type first_idx = M.find_first();

    if (first_idx != npos )
    {
      dbitset rM ( M.count(), 0 );
      Yid = rM;

      for( size_t i{0}; i < X.size(); ++i )
        Xid.push_back( rM );

      size_type current_idx = first_idx;
      size_t k = 0;
      do
      {
        Yid[k] = Y[current_idx];
        for( size_t i{0}; i < X.size(); ++i )
          Xid[i][k] = X[i][current_idx]; 
        k++;  
        current_idx = M.find_next( current_idx );
      } while ( current_idx != boost::dynamic_bitset<>::npos );
      Xid.erase( Xid.begin()+idx );
    }

    return std::make_pair( Xid, Yid );
  }
  #pragma endregion cofactors

  #pragma region cum sum
    double Pk_f( uint64_t const& k, uint64_t const& N0, uint64_t const& N1, uint const& n )
    {
      
      uint64_t Nh = std::max(N0,N1);
      uint64_t Nl = std::min(N0,N1);
      uint64_t n_inf = 10;
      if( n>n_inf || Nl==0 || Nh == 0)
        return (k==0 ? (double)1:(double)0);
      double Pk = 1;
      if(k>Nl)
      {
        "k-Nl check";
        return 0;
      }
      if( pow(2,n-1)+k < Nh+Nl) 
      {
        "pow check";
        return 0;
      }
      if(Nh==pow(2,n-1))
      {
        if( k==Nl )
        {
          return 1;
        } 
      }

      for(auto j=0;j<Nl-k;++j)
      {
        Pk*=(1-(double)Nh/(pow(2,n-1)-j));
      }
      
      if(k!=0)
      {
        for( auto j = 0; j<k ; ++j )
        {
          double Ak = (double)(Nl-j)/(j+1);
          double Bk = (double)(Nh-j)/(pow(2,n-1)-Nl+j+1);
          Pk *= Ak*Bk;
        }
      }

      return Pk;
    }
    std::pair<double,double> M1M2k(uint64_t const& N0, uint64_t const& N1, uint64_t const& n )
    {
      uint64_t Nh = std::max(N0,N1);
      uint64_t Nl = std::min(N0,N1);
      uint64_t n_inf = 32;
      if( n>n_inf )
        return std::make_pair(0,0);
      uint64_t kmin = std::max(1, (int)((Nh+Nl)-pow(2,n-1)));
      double Pk = Pk_f( kmin, N0, N1, n );
      double M1 = kmin*Pk;
      double M2 = kmin*kmin*Pk;
      for( uint64_t k=kmin+1; k< Nl+1 ;++k)
      {
        
        //double Ak = (double)(Nh-k+1)/(pow(2,n-1)-(Nh+Nl)+k);
        //Ak = Pk*(Nl-k+1);
        double Ak = k*Pk_f( k, N0, N1, n );
        M1 += Ak;
        M2 += Ak*k;
      }
      return std::make_pair(M1,sqrt(M2-M1*M1));
    }
    uint64_t num_intersections(uint64_t const& N0, uint64_t const& N1, uint64_t const& n)
    {
      auto R = M1M2k(N0,N1,n);
      return std::max((int)floor(R.first-3*R.second),1);
    }

    double CumSum(uint64_t const& kmax, uint64_t const& N0, uint64_t const& N1, uint64_t const& n )
    {
      double CS =0;
      for(uint64_t k=0; k<kmax+1; ++k)
      {
        double dP=Pk_f( k, N0, N1, n );
        //std::cout << "Pk_f("<< k<<","<< N0<<"," <<N1 <<","<< n<<" )=" << dP<<std::endl;
        CS+=dP;
      }
      //std::cout << CS << std::endl;
      return CS;
    }
  #pragma endregion cum sum

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

  #pragma region is_xor_decomposable
  bool is_xor_decomposable(  std::pair<dbitset_vector, dbitset> const XY0, std::pair<dbitset_vector, dbitset> const XY1 )
  {
    size_t count_neg = 0;
    std::unordered_map<std::string, bool> str_nodes0;
    std::unordered_map<std::string, bool> already;

    size_t N0 = 0;
    /* fill hash table */
    if( XY0.first.size() == 0 )
      return false;
    for( size_t k {0u}; k < XY0.first[0].size(); ++k )
    {
      dbitset pattern;
      for( size_t j{0}; j < XY0.first.size(); ++j )
        pattern.push_back( XY0.first[j][k] );
      std::string s;
      to_string( pattern, s );
      if( str_nodes0.find(s) == str_nodes0.end() )
      {
        N0++;
        str_nodes0.insert( std::make_pair(s, XY0.second[k]));
      }
      else if ( str_nodes0.at(s) != XY0.second[k] )
      {
        return false;
      }
    }

    uint64_t N1 = 0;
    for ( size_t k {0u}; k < XY1.first[0].size(); ++k )
    {
      dbitset pattern;
      for( size_t j{0}; j < XY1.first.size(); ++j )
        pattern.push_back( XY1.first[j][k] );
      std::string s;
      to_string( pattern, s );
      if( already.find(s) == already.end() )
        N1++;

      if( str_nodes0.find(s) != str_nodes0.end() )
      {   
        if( str_nodes0.at(s) == XY1.second[k] )
          return false;
        else
        {
          if(already.find(s) == already.end())
            count_neg++;
        }
      }
      if(already.find(s)==already.end())
        already.insert(std::make_pair(s,XY1.second[k]));
    }
    
    auto n = XY0.first.size()+1;
    auto R = M1M2k(N0,N1,n);
    auto Min = std::max(1,(int)floor(R.first-R.second));
    auto Max = std::max(1,(int)ceil(R.first+R.second));
      
    bool condition1 = ( ( CumSum( count_neg+(int)ceil(R.second), N0, N1, n ) >= 1-0.001 ) && ( count_neg >=2 ) );
    bool condition2 = ( count_neg >= 1 ) ;
    bool is_satisfied = _ps.use_cumsum ? condition1 : condition2;
    /*std::cout << "c1 " << condition1 << std::endl; 
    std::cout << "c2 " << condition2 << std::endl; 
    std::cout << "cs " << is_satisfied << std::endl; */
    if( is_satisfied ) // modifications are possible
      return true;

    return false;
    }
  #pragma endregion is_xor_decomposable

  #pragma region chatterjee
  std::pair<std::string, bool> chatterjee_method( dbitset_vector& X, dbitset Y )
  {
      bool is_exact = true;
      uint64_t N = X.size();
      uint64_t pow2N = pow(2,N);
      dbitset bit02N( X[0].size(), 0 );
      dbitset Kmask, new_values;
      new_values = bit02N;
      std::string tt = "";
      std::default_random_engine generator;
      std::bernoulli_distribution distribution(0.5);
      uint64_t C0, C1;
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
        auto r = distribution(generator);

        if( ( C1 > C0 ) || ( ( C1 == C0 ) && ( r >= 0.5 ) ) )
        {
          new_values |= Kmask;
          tt = "1"+tt;
        }
        else if( ( C1 < C0 ) || ( ( C1 == C0 ) && ( r < 0.5 ) ) )
        {
          tt = "0"+tt;
        }
        if( C1 != 0 && C0 != 0 )
          is_exact &= false;
      }
      X.push_back( new_values );
      return std::make_pair( tt, is_exact );
  }
  signal<klut_network> apply_chatterjee( std::vector<signal<klut_network>>& support, dbitset_vector& X, dbitset const& Y)
  {
    auto tt_str = chatterjee_method( X, Y ).first;
    kitty::dynamic_truth_table tt( support.size() );
    kitty::create_from_binary_string( tt, tt_str );
    return _klut.create_node( support, tt );
  }
  #pragma endregion chatterjee

  #pragma region decomposition_checks
  it_top_decomposition is_top_decomposable( auto& XY0, auto& XY1 )
  {
    if( (XY0.second.size() == 0) || (XY0.second.count() == 0) ) // F0 = 0
    {
      _cnt._and++;
      return it_top_decomposition::and_ ;
    }
    else if( (XY1.second.size() != 0) && (XY1.second.count() == XY1.second.size()) ) // F1 = 1
    {
      _cnt._or++;
      return it_top_decomposition::or_ ;
    }
    else if( (XY1.second.size() == 0) || (XY1.second.count() == 0) ) // F1 = 0
    {
      _cnt._lt++;
      return it_top_decomposition::lt_ ;
    }
    else if( (XY0.second.size() != 0) && (XY0.second.count() == XY0.second.size()) ) // F0 = 1
    {
      _cnt._le++;
      return it_top_decomposition::le_ ;
    }
    else if( _ps.try_xor_decomposition && is_xor_decomposable( XY0, XY1 ) )
    {
      _cnt._xor++;
      return it_top_decomposition::xor_ ;
    }
    else
      return it_top_decomposition::none ;
  }

  bool is_bottom_decomposable(  std::vector<signal<klut_network>>& support, dbitset_vector& X, dbitset const& Y, double Imax,
                                      std::vector<double>& Ivect, std::vector<size_t>& IDXvect )
    {
      quicksort_by_attribute( IDXvect, Ivect, 0, Ivect.size()-1 );

      auto original_support = support;

      std::vector<size_t> indeces2;
      std::vector<signal<klut_network>> support2;
      bool flag = false;

      double Isupp, Ifnew, Ifr, Ifc, Ifrc;  
      dbitset_vector Xtmp;

      for( size_t i=0;  i < IDXvect.size()-1 ; ++i )
      {
          auto r = IDXvect[i];
          auto c = IDXvect[i+1];
          indeces2 = { r, c };
          support2 = { original_support[r], original_support[c] };

          std::string Sr, Sc;
          uint64_t Sr_64t = original_support[r];
          uint64_t Sc_64t = original_support[c];
          Sr = std::to_string( Sr_64t );
          Sc = std::to_string( Sc_64t );
          std::string support_key = Sr + " " + Sc;

          Xtmp = { X[r], X[c] };
          auto chj = chatterjee_method( Xtmp, Y );
          auto tt_str = chj.first;
          assert( Xtmp.size() == 3 );

          if( (_Icoll.Frc).find( support_key ) == (_Icoll.Frc).end() )
          {
            Isupp = kitty::mutual_information( std::vector{ X[r], X[c] }, Y ); 
            Ifnew = kitty::mutual_information( Xtmp[2], Y );
            Ifr = kitty::mutual_information( std::vector{ Xtmp[2], X[r] }, Y );
            Ifc = kitty::mutual_information( std::vector{ Xtmp[2], X[c] }, Y );
            Ifrc = kitty::mutual_information( std::vector{ Xtmp[2], X[r], X[c] }, Y );
            (_Icoll.Fnew).insert(std::make_pair( support_key, Ifnew ));
            (_Icoll.Frc).insert(std::make_pair( support_key, Ifrc ));
            (_Icoll.Fr).insert(std::make_pair( support_key, Ifr ));
            (_Icoll.Fc).insert(std::make_pair( support_key, Ifc ));
            (_Icoll.supp).insert(std::make_pair( support_key, Isupp ));
          } 
          else
          {
            Isupp = (_Icoll.supp).at( support_key );
            Ifnew = (_Icoll.Fnew).at( support_key );
            Ifr = (_Icoll.Fr).at( support_key );
            Ifc = (_Icoll.Fc).at( support_key );
            Ifrc = (_Icoll.Frc).at( support_key );
          }
 
          bool exact_flag = _ps.is_bottom_exact ? chj.second : true;
          /*std::cout << "is bottom exact " << _ps.is_bottom_exact << std::endl;
          std::cout << "chj.second " << chj.second << std::endl;
          std::cout << "exact flag " << exact_flag << std::endl;
*/
          if( ( Isupp == Ifnew ) && (Ifrc == Ifnew) && ( Ifr == Ifnew ) && ( Ifc == Ifnew) && exact_flag )
          {

            kitty::dynamic_truth_table tt(2u);
            create_from_binary_string( tt, tt_str );
            support.push_back( _klut.create_node( support2, tt ) );
            X.push_back( Xtmp[2] );
            X.erase(X.begin()+std::max( r, c ) );
            X.erase(X.begin()+std::min( r, c ) );
            support.erase(support.begin()+std::max( r, c ));
            support.erase(support.begin()+std::min( r, c ));
            return true;
          }
        //}
      }
      return false;
    }
  #pragma endregion decomposition_checks

  #pragma region creation procedures
  bool is_new_created(  std::vector<signal<klut_network>>& support, dbitset_vector& X, dbitset const& Y, double Imax,
                                      std::vector<double>& Ivect, std::vector<size_t>& IDXvect )
    {
      quicksort_by_attribute( IDXvect, Ivect, 0, Ivect.size()-1 );

      auto original_support = support;

      std::vector<size_t> indeces2;
      std::vector<signal<klut_network>> support2;
      bool flag = false;

      double Isupp, Ifnew, Ifr, Ifc, Ifrc;  
      dbitset_vector Xtmp;

      for( size_t i=0;  i < IDXvect.size()-1 ; ++i )
      {
          auto r = IDXvect[i];
          auto c = IDXvect[i+1];
          indeces2 = { r, c };
          support2 = { original_support[r], original_support[c] };

          std::string Sr, Sc;
          uint64_t Sr_64t = original_support[r];
          uint64_t Sc_64t = original_support[c];
          Sr = std::to_string( Sr_64t );
          Sc = std::to_string( Sc_64t );
          std::string support_key = Sr + " " + Sc;

          Xtmp = { X[r], X[c] };
          auto chj = chatterjee_method( Xtmp, Y );
          auto tt_str = chj.first;
          assert( Xtmp.size() == 3 );

          if( (_Icoll.Fnew).find( support_key ) == (_Icoll.Fnew).end() )
          { 
            Isupp = kitty::mutual_information( std::vector{ X[r], X[c] }, Y ); 
            Ifnew = kitty::mutual_information( Xtmp[2], Y );
            (_Icoll.Fnew).insert(std::make_pair( support_key, Ifnew ));
            if( Ifnew >= Imax )
            {
              kitty::dynamic_truth_table tt(2u);
              create_from_binary_string( tt, tt_str );
              support.push_back( _klut.create_node( support2, tt ) );
              X.push_back( Xtmp[2] );
                return true;
            }
          } 
          else
          {
            Ifnew = (_Icoll.Fnew).at( support_key );
          }
 
          //bool exact_flag = _ps.is_bottom_exact ? chj.second : true;
          /*std::cout << "is bottom exact " << _ps.is_bottom_exact << std::endl;
          std::cout << "chj.second " << chj.second << std::endl;
          std::cout << "exact flag " << exact_flag << std::endl;
*/
        //}
      }
      return false;
    }
  #pragma endregion creation procedures

  #pragma region idsd
  signal<klut_network> idsd_step( std::vector<signal<klut_network>> support, dbitset_vector& X, dbitset & Y )
  {
    if( X.size() == 0 )
      return _klut.get_constant( false );    

    if( X[0].size() == 0 )
      return _klut.get_constant( false );

    assert( ( support.size() == X.size() ) );
    assert( ( X[0].size() == Y.size() ) );

    if( Y.count() == 0 ) // contraddiction
      return _klut.get_constant( false );
    else if( Y.count() == Y.size() ) // tautology
      return _klut.get_constant( true );


    if( support.size() <= _ps.max_sup )
    {
      _cnt._ctj++;
      return apply_chatterjee( support, X, Y );
    }
    uint64_t idx = 0;
    double Inew = 0;
    double Imax = 0;
    std::vector<double> Ivect;
    std::vector<size_t> IDXvect;
    
    if( _ps.is_informed )
    {
      size_t i = 0;
      for( auto x : X )
      {
        Inew = kitty::mutual_information( x, Y );
        if( Inew >= Imax )
        {
          idx = i;
          Imax = Inew;
        }
        IDXvect.push_back( i++ );
        Ivect.push_back( Inew );
      }
    }
    auto XY0 = compute_cofactor( X, Y, 0, idx );
    auto XY1 = compute_cofactor( X, Y, 1, idx ); 

    auto reduced_support = support;
    reduced_support.erase( reduced_support.begin() + idx );

    if( _ps.try_top_decomposition )
    {
      auto res = is_top_decomposable( XY0, XY1 );
      if ( res != it_top_decomposition::none )
      {
        _Icoll.clear();
        switch ( res )
        {
        default:
          assert( false );
        case it_top_decomposition::and_:
        {
          auto F1 = idsd_step( reduced_support, XY1.first, XY1.second );
          return _klut.create_and( support[idx], F1 );
        }
        case it_top_decomposition::or_:
        {
          auto F0 = idsd_step( reduced_support, XY0.first, XY0.second );
          return _klut.create_or( support[idx], F0 );
        }
        case it_top_decomposition::lt_:
        {
          auto F0 = idsd_step( reduced_support, XY0.first, XY0.second );
          return _klut.create_lt( support[idx], F0 );
        }
        case it_top_decomposition::le_:
        {  
          auto F1 = idsd_step( reduced_support, XY1.first, XY1.second );
          return _klut.create_le( support[idx], F1 );
        }
        case it_top_decomposition::xor_:
        {
          remove_column_and_invert( X, Y, idx ); 
          return _klut.create_xor( support[idx] , idsd_step( reduced_support, X, Y ) );
        }
        }
      }
    }
    if( _ps.try_bottom_decomposition && is_bottom_decomposable( support, X, Y, Imax, Ivect, IDXvect ) )
    {
      _cnt._btm++;
      return idsd_step( support, X, Y );
    }
    if( _ps.try_creation && is_new_created( support, X, Y, Imax, Ivect, IDXvect ) )
    {
      _cnt._cre++;
      return idsd_step( support, X, Y );
    } 
    _Icoll.clear();
    auto F0 = idsd_step( reduced_support, XY0.first, XY0.second );
    _Icoll.clear();
    auto F1 = idsd_step( reduced_support, XY1.first, XY1.second );

    auto f0 = _klut.create_and( _klut.create_not( support[idx] ), F0 );
    auto f1 = _klut.create_and( support[idx], F1 );

    return _klut.create_or( f1, f0 );

  }
  #pragma endregion idsd

  signal<klut_network> run()
  {
    return idsd_step( _examples.signals, _examples.partial.first, _examples.partial.second );
  }
private:
  klut_network& _klut;
  lfeNtk<klut_network> _examples;
  signal<klut_network> _support;
  uint64_t _num_out;
  Istorage _Icoll;
  it_decomposition_params _ps;
public:
  detection_counter _cnt;
};

} // namespace detail

#pragma region iwls2022
void print_LFE( auto LFE, bool only_complete = false )  
{
  std::cout << "complete:" << std::endl;
  for( auto x : LFE.complete.first)
  {
    kitty::print_binary(x);std::cout<<std::endl;
  }
  auto n = LFE.complete.first[0].num_bits();
  for( auto i = 0u; i < n; ++i )
    std::cout << "-";
  std::cout << std::endl;
  for( auto x : LFE.complete.second)
  {
    kitty::print_binary(x);std::cout<<std::endl;
  } 
  if( !only_complete )
  {
    std::cout << "partial:" << std::endl;
    for( auto x : LFE.partial.first)
    {
      std::cout << x <<std::endl;
    }
    n = LFE.partial.first[0].size();
    for( auto i = 0u; i < n; ++i )
      std::cout << "-";
    std::cout << std::endl;

    std::cout << LFE.partial.second <<std::endl; 
  }
}

/* Algorithms:
 * it_decomposition      : multi output given an original klut network
 * trivial_decomposition : decompose all the outputs
*/
void it_decomposition( klut_network& klut, it_decomposition_params& ps )
{
  if( ps.is_trivial )
  {
    std::vector<uint64_t> output_nodes;
    klut.foreach_po( [&]( auto const& node, auto index ) {
      output_nodes.push_back( node );
    } );
    auto num_pos = klut.num_pos();
    auto examples = graph_to_lfe( klut );
    auto tt0 = examples.partial.second;
    std::vector<signal<klut_network>> outputs;
  
    for( auto i = 0; i < examples.complete.second.size(); ++i ) // you could work on the order
    {
      examples = graph_to_lfe( klut, i );
      detail::it_decomposition_impl impl( klut, examples, ps );
      klut.substitute_node( output_nodes[i], impl.run() );
    }
  }
  else
  {
    std::vector<uint64_t> output_nodes;
    klut.foreach_po( [&]( auto const& node, auto index ) {
      output_nodes.push_back( node );
    } );
    auto num_pos = klut.num_pos();
    auto examples = graph_to_lfe( klut );
    std::vector<signal<klut_network>> outputs;

    detail::it_decomposition_impl impl( klut, examples, ps );
    outputs = {impl.run()};
    if( num_pos > 1 )
    {
      outputs = {};
      auto examples = graph_to_lfe( klut );   
      for( auto i = 0; i < examples.complete.second.size(); ++i ) // you could work on the order
      {
        auto examplesN = graph_to_lfe( klut, i );
        detail::it_decomposition_impl impl( klut, examplesN, ps );
        klut.substitute_node( output_nodes[i], impl.run() );
      }
    }
    else
    {
      klut.substitute_node( output_nodes[0], outputs[0] );
    }
  }
}
#pragma endregion iwls2022

#pragma region iwls2020
template<typename Ntk>
bool simulate_input( mockturtle::dbitset const& input_pattern, Ntk& ntk )
{
  std::vector<bool> inpt_v;
  for( uint64_t k{0u}; k<input_pattern.size();++k )
  {
    inpt_v.push_back( ( ( input_pattern[k] == 1 ) ? true : false ) );
  }

  return simulate<bool>( ntk, default_simulator<bool>( inpt_v ) )[0];
}

template<typename Ntk>
double compute_accuracy( mockturtle::dbitset_vector const& X, mockturtle::dbitset const& Y, Ntk& ntk )
{
  double acc = 0;
  double delta_acc;
  for( uint64_t k {0u}; k < X[0].size(); ++k )
  {
    dbitset ipattern;
    for ( uint64_t j {0u}; j < X.size(); ++j )
      ipattern.push_back( X[j][k] );
        
    delta_acc = ( ( simulate_input( ipattern, ntk ) == Y[k] ) ? (double)1.0/X[0].size() : 0.0 );
    acc += delta_acc;
  }
  return acc;
}

detail::it_decomposition_impl it_decomposition_iwls20( XYdataset& dt, klut_network& klut, it_decomposition_params& ps )
{
  lfeNtk<klut_network> examples;
  examples.partial = std::make_pair( dt.X, dt.Y );

  muesli_params mps;
  ps.max_sup = 2;
  klut = muesli( examples.partial.first, examples.partial.second, mps );

  //for( size_t i = 0; i < dt.nin; ++i )
  //  examples.signals.push_back( klut.create_pi() );
  
  detail::it_decomposition_impl impl( klut, examples, ps );
  if( !mps.is_exact_fn )
    klut.create_po( impl.run() );
  return impl;
}
#pragma endregion iwls2020
} // namespace mockturtle