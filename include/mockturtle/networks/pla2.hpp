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
  \file cover.hpp
  \brief single output cover logic network implementation
  \author Andrea Costamagna
*/

#pragma once

#include <iostream>
//#include <catch.hpp>
#include <ctime>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include <bitset>
#include <cmath>
#include <random>
#include <boost/dynamic_bitset.hpp>
#include <unordered_map>

#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>

#include <fstream>
#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>

#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/algorithms/aig_algebraic_rewriting.hpp>


namespace mockturtle
{
  struct index_to_signal
  {
    index_to_signal()
    {
      storage.reserve( 10000u );
    }
    void insert( uint64_t pla_index, uint64_t klut_signal )
    {
      storage[pla_index] = klut_signal;
    }
    std::unordered_map<uint64_t, uint64_t> storage;
  };

  struct calls_counter
  {
    uint32_t F1T = 0;
    uint32_t F1C = 0;
    uint32_t F0T = 0;
    uint32_t F0C = 0;
    uint32_t Chj = 0;
    uint32_t XOR = 0;
    uint32_t Fo = 0;
    uint32_t Frc = 0;
    double train_acc = 0;
    double test_acc = 0;
    double valid_acc = 0;
  };

  class pla2_network
  {
    #pragma region Types and constructors
    using dbitset = boost::dynamic_bitset<>;
    using dbitset_vector = std::vector<dbitset>;
    
    public:
    pla2_network( dbitset_vector input_nodes, dbitset_vector output_nodes, uint64_t max_sup = 4 )
      : _nodes( input_nodes ),
        _outputs( output_nodes ),
        _num_nodes( input_nodes.size() ),
        _max_sup( max_sup )
        {
          _init();
        }

    protected:
      inline void _init()
      {
        for ( uint64_t i {0u}; i < _nodes.size() ; ++i )
        {
          auto pi = _klut.create_pi();
          _itos.insert( i, pi );
        }
      }
    #pragma endregion

    #pragma region visual
    public:
    void print_pla()
    {

      std::cout << "===";
      for ( uint64_t i {0u}; i < _nodes[0].size(); ++i )
        std::cout << "=";
      std::cout << std::endl;
      for ( uint64_t i {0u}; i < _nodes.size(); ++i )
        std::cout << "X" << i << ":" << _nodes.at(i)  << std::endl ;
      std::cout << "---";
      for ( uint64_t i {0u}; i < _nodes[0].size(); ++i )
        std::cout << "-";
      std::cout << std::endl; 
      for ( uint64_t i {0u}; i < _outputs.size(); ++i )
        std::cout << "Y" << i << ":" << _outputs.at(i)  << std::endl ;
      std::cout << "===";
      for ( uint64_t i {0u}; i < _nodes[0].size(); ++i )
        std::cout << "=";
      std::cout << std::endl;
    }

    void print_pla( std::pair<dbitset_vector,dbitset_vector> pair )
    {

      std::cout << "===";
      for ( size_t i {0u}; i < pair.first[0].size(); ++i )
        std::cout << "=";
      std::cout << std::endl;
      for ( uint64_t i {0u}; i < pair.first.size(); ++i )
        std::cout << "X" << i << ":" << pair.first.at(i)  << std::endl ;
      std::cout << "---";
      for ( uint64_t i {0u}; i < pair.first[0].size(); ++i )
        std::cout << "-";
      std::cout << std::endl; 
      for ( uint64_t i {0u}; i < pair.second.size(); ++i )
        std::cout << "Y" << i << ":" << pair.second.at(i)  << std::endl ;
      std::cout << "===";
      for ( uint64_t i {0u}; i < pair.first[0].size(); ++i )
        std::cout << "=";
      std::cout << std::endl;
    }

    void print_Pr( std::vector<double> Vpr )
    {
      for( uint32_t k {0u}; k < Vpr.size(); ++k )
      {
        boost::dynamic_bitset<> A (sqrt(Vpr.size()),k);
        std::cout << "|P(" << A << ")=" << Vpr[k] << std::endl; 
      }
    }

    void print_features()
    {
      depth_view_params ps;
      ps.count_complements = true;
      depth_view depth_aig{_aig, {}, ps};
      std::cout << ".F1T : " << _ID << _cnt.F1T << std::endl;
      std::cout << ".F1C : " << _ID<< _cnt.F1C << std::endl;
      std::cout << ".F0T : " << _ID<< _cnt.F1T << std::endl;
      std::cout << ".F0C : " << _ID<< _cnt.F1C << std::endl;
      std::cout << ".XOR : " << _ID<< _cnt.XOR << std::endl;
      std::cout << ".2-OR : " << _ID<< _cntOR << std::endl;
      std::cout << ".2-LT : " << _ID<< _cntLT << std::endl;
      std::cout << ".2-LE : " << _ID<< _cntLE << std::endl;
      std::cout << ".2-AND : " << _ID<< _cntAND << std::endl;
      std::cout << ".2-XOR : " << _ID<< _cntXOR << std::endl;

      std::cout << ".Fo : " << _ID<< _cnt.Fo << std::endl;
      std::cout << ".Frc : "<< _ID << _cnt.Frc << std::endl;
      std::cout << ".c   : " << _ID<< _cnt.Chj << std::endl;
      std::cout << ".g   : " << _ID<< depth_aig.num_gates() << std::endl;
      std::cout << ".s   : "<< _ID << depth_aig.size() << std::endl;
      std::cout << ".d   : " << _ID<< depth_aig.depth() << std::endl;
      std::cout << ".l   : " << _ID<< _cnt.train_acc << std::endl;
      std::cout << ".t   : " << _ID<< _cnt.test_acc << std::endl;
      std::cout << ".v   : "<< _ID << _cnt.valid_acc << std::endl;
      std::cout << ".a "<< _ID<< _duration <<'\n';


      if( _has_file )
      {
        std::ofstream _myfile;
        _myfile.open ( _pathTOfile );
        _myfile << ".b " << _IDs << std::endl; 
        _myfile << ".e muesli enhanced" << std::endl;
        _myfile << ".F1T " << _cnt.F1T << std::endl;
        _myfile << ".F1C " << _cnt.F1C << std::endl;
        _myfile << ".F0T " << _cnt.F1T << std::endl;
        _myfile << ".F0C " << _cnt.F1C << std::endl;
        _myfile << ".XOR " << _cnt.XOR << std::endl;
        _myfile << ".cntOR " << _cntOR << std::endl;
        _myfile << ".cntLT " << _cntLT << std::endl;
        _myfile << ".cntLE " << _cntLE << std::endl;
        _myfile << ".cntAND " << _cntAND << std::endl;
        _myfile << ".cntXOR " << _cntXOR << std::endl;
        _myfile << ".Fo " << _cnt.Fo << std::endl;
        _myfile << ".Frc " << _cnt.Frc << std::endl;
        _myfile << ".c " << _cnt.Chj << std::endl;
        _myfile << ".g " << depth_aig.num_gates() << std::endl;
        _myfile << ".s " << depth_aig.size() << std::endl;
        _myfile << ".d " << depth_aig.depth() << std::endl;
        _myfile << ".l " << _cnt.train_acc << std::endl;
        _myfile << ".t " << _cnt.test_acc << std::endl;
        _myfile << ".v " << _cnt.valid_acc << std::endl;
        _myfile <<".a "<< _duration <<'\n';
        _myfile.close();
      }
      
    }

    void add_output_file( std::string pathTOfile, std::string ID = "" )
    {
      _ID = "[" + ID + "] ";
      _IDs = ID;
      _pathTOfile = pathTOfile;
      _has_file = true;
    }
    void set_preferences( bool top_decompose, bool bottom_decompose,bool dontknows, bool informed )
    {
      _top_decompose = top_decompose;
      _bottom_decompose = bottom_decompose;
      _dontknows = dontknows;
      _informed = informed;
    }
    #pragma endregion

    #pragma region statistics
    std::vector<double> Pr( dbitset_vector const& nodes )
    {
      std::vector<double> Vpr = {};
      uint64_t N = nodes.size();
      uint64_t pow2N = pow(2,N);
      dbitset bit2N0( nodes[0].size(), 0 );
      for( uint64_t k {0u}; k < pow2N; ++k )
      {
        auto Vk = ~bit2N0;
        dbitset maskN( N, k );
        for( uint64_t j {0u}; j < N; ++j )
        {
          if( maskN[j] == 1 )
            Vk &= nodes[j];
          else if( maskN[j] == 0 )
            Vk &= ~nodes[j];
          else
            std::cerr << "invalid" << std::endl;
        }
        Vpr.push_back( (double)Vk.count()/nodes[0].size() );
      }
      return Vpr;
    }

    double H( dbitset_vector const& nodes )
    {
      auto P =  Pr( nodes );
      double h = 0;
      for( uint64_t k{0u}; k < P.size(); ++k )
      {
        if( P[k] != 0 )
          h += -P[k]*log2(P[k]);
      }
      return h;
    }

    double MI( dbitset_vector const& X, dbitset_vector const& Y, std::vector<uint64_t> support = {} )
    {
      auto N = X[0].size();
      assert( N == Y[0].size() );
      auto XY = X;
      double Hx, Hy, Hxy;

      for( uint32_t k {0u}; k<Y.size(); ++k )
        XY.push_back(Y[k]); 
      Hx = H(X);
      Hy = H(Y);
      Hxy = H(XY);
      return Hx+Hy-Hxy;
    }

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
    #pragma endregion

    #pragma region cofactors_manipulation
    std::pair<dbitset_vector, dbitset_vector> compute_cofactor( dbitset_vector const& X, dbitset_vector const& Y, 
                                                                 uint64_t const idx, uint64_t const& id )
    {
      
      //std::cout << "idx: " << idx << std::endl;
      if( X.size()==0 )
        return std::make_pair(X,Y);
      assert( X[0].size() == Y[0].size() );
      assert( idx < X.size() );

      dbitset M;
      M = (id == 1) ?  X[idx] : ~X[idx];
      dbitset_vector Xid, Yid;

      typedef boost::dynamic_bitset<>::size_type size_type;
      const size_type npos = boost::dynamic_bitset<>::npos;
      size_type first_idx = M.find_first();

      if (first_idx != npos )
      {
        dbitset rM ( M.count(), 0 );
        Yid = { rM };
        for( size_t i{0}; i < X.size(); ++i )
          Xid.push_back( { rM } );
  
        size_type current_idx = first_idx;
        size_t k = 0;
        do
        {
          Yid[0][k] = Y[0][current_idx];
          for( size_t i{0}; i < X.size(); ++i )
            Xid[i][k] = X[i][current_idx]; 
          
          current_idx = M.find_next(current_idx);
          k++;
        } while ( current_idx != boost::dynamic_bitset<>::npos);
        Xid.erase( Xid.begin()+idx );
      }

    return std::make_pair( Xid, Yid );
    }
    #pragma endregion

    #pragma region function and node creation
    std::string create_function( dbitset_vector& X, dbitset_vector const& Y, std::vector<uint64_t> indeces = {} )
    {
      if( indeces.size() == 0 )
      {
        for( size_t i{0}; i<X.size(); ++i )
          indeces.push_back( i );
      }
      uint64_t N = indeces.size();
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
        C1 = ( Kmask & Y[0] ).count();
        C0 = ( Kmask & ~Y[0] ).count();
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
      }
      X.push_back( new_values );
      return tt;
    }

    uint64_t create_klut_node( std::vector<uint64_t> support, std::string const& tt_str )
    {
      kitty::dynamic_truth_table tt( support.size() );
      create_from_binary_string( tt, tt_str );
      std::vector<uint64_t> klut_signals;
      for ( uint64_t i {0u}; i < support.size(); ++i )
        klut_signals.push_back( _itos.storage[support[i]] );

      auto f0 = _klut.create_node( klut_signals, tt );
      _itos.insert( _num_nodes ,f0 );
      _num_nodes++;
      return f0;
    }
    #pragma endregion

    uint64_t show_max(uint64_t const& N0, uint64_t const& N1, uint64_t const& n )
    {
      uint64_t kmax = 0;
      double Pmax = 0;
      double Pnew = 0;
      for(uint64_t k = 0; k<std::min(N0,N1);++k)
      {
        Pnew = Pk_f( k, N0, N1, n );
        if( Pnew > Pmax )
        {
          Pmax = Pnew;
          kmax = k;
        }
      }
      return kmax;
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
    #pragma region DSD
    bool is_F1_not_F0(  std::pair<dbitset_vector, dbitset_vector> const XY0, std::pair<dbitset_vector, dbitset_vector> const XY1, 
                        uint64_t min_intersection = 0 )
    {

      uint64_t count_neg = 0;
      std::unordered_map<std::string, double> str_nodes0;
      std::unordered_map<std::string, double> already;

      uint64_t N0 = 0;
      /* fill hash table */
      for( uint64_t k {0u}; k < XY0.first[0].size(); ++k )
      {
        dbitset pattern;
        for( size_t j{0}; j < XY0.first.size(); ++j )
          pattern.push_back( XY0.first[j][k] );
        std::string s;
        to_string( pattern, s );
        if( str_nodes0.find(s) == str_nodes0.end() )
          N0++;
        str_nodes0.insert(std::make_pair(s,XY0.second[0][k]));
      }

      uint64_t N1 = 0;

      for ( uint64_t k {0u}; k < XY1.first[0].size(); ++k )
      {
        dbitset pattern;
        for( size_t j{0}; j < XY1.first.size(); ++j )
          pattern.push_back( XY1.first[j][k] );
        std::string s;
        to_string( pattern, s );
        if(already.find(s)==already.end())
        {
          N1++;
        }

        if( str_nodes0.find(s) != str_nodes0.end() )
        {   

          if( str_nodes0.at(s) == XY1.second[0][k] )
          {
            return false;
          }
          else
          {
            if(already.find(s) == already.end())
            {
              count_neg++;
            }
          }
        }
        if(already.find(s)==already.end())
        {
          already.insert(std::make_pair(s,XY1.second[0][k]));
        }
      }
      auto n = XY0.first.size()+1;
      auto R = M1M2k(N0,N1,n);
      auto Min = std::max(1,(int)floor(R.first-R.second));
      auto Max = std::max(1,(int)ceil(R.first+R.second));
      
//
      if(  (CumSum( count_neg+(int)ceil(R.second), N0, N1, n ) >= 1-0.001 ) &&(count_neg >=2) )//CumSum( count_neg, N0, N1, n ) >= 1 ) // CONSIDER CHANGING TO >= 0 if motivated
      {
        std::cout << "cum sum = "<< CumSum( count_neg, N0, N1, n ) <<std::endl;
        return true;
      }
      return false;
    }

    bool is_F1_F0(  std::pair<dbitset_vector, dbitset_vector> const XY0, 
                    std::pair<dbitset_vector, dbitset_vector> const XY1, std::vector<uint64_t>& where1,
                    uint64_t min_intersection = 0 )
    {
      uint64_t count = 0;
      std::unordered_map<std::string, double> str_nodes0;
      if( ( XY0.first.size() == 0 ) || ( XY1.first.size() == 0 ) )
        return true;
      if( (XY0.first[0].size() == 0) || (XY1.first[0].size() == 0) )
        return true;

      /* fill hash table */
      for( uint64_t k {0u}; k < XY0.first[0].size(); ++k )
      {
        dbitset pattern;
        for( size_t j{0}; j < XY0.first.size(); ++j )
          pattern.push_back( XY0.first[j][k] );
        std::string s;
        to_string( pattern, s );
        str_nodes0.insert(std::make_pair(s,XY0.second[0][k]));
      }

      for ( uint64_t k {0u}; k < XY1.first[0].size(); ++k )
      {
        dbitset pattern;
        for( size_t j{0}; j < XY1.first.size(); ++j )
          pattern.push_back( XY1.first[j][k] );
        std::string s;
        to_string( pattern, s );

        if( str_nodes0.find(s) != str_nodes0.end() )
        {
          if( str_nodes0.at(s) == XY1.second[0][k] )
          {
            where1.push_back(k);
            count++;
          }
          else
          {
            return false;
          }
        }
      }

      if( count >= min_intersection ) // CONSIDER CHANGING TO >= 0 if motivated
      {
        return true;
      }
      return false;
    }

    void remove_column_and_invert( dbitset_vector& X, dbitset_vector& Y, uint64_t idx_max )
    {
      Y[0] ^= X[idx_max];
      X.erase( X.begin() + idx_max );
    }

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

    struct new_nodes_storage
    {
    bool is_created;
    std::vector<uint64_t> support;
    std::vector<uint64_t> indeces;
    double I;
    std::string tt;
    bool rc_del;
    };

void add_tt_to_hash( std::string const& tt_new )
{
  if( _tt_counter.find( tt_new ) == _tt_counter.end() )
  {
    _tt_counter.insert(std::make_pair( tt_new, 1 ));
  }
  else
  {
    _tt_counter[tt_new] += 1;
  }
}

    bool try_ME_step( std::vector<uint64_t>& support, dbitset_vector& X, dbitset_vector& Y, double& Imax )
    {
      new_nodes_storage nns;
      nns.is_created = false;

      dbitset new_node;
      dbitset_vector Xtmp;
      std::vector<uint64_t> original_support = support;

      std::vector<uint64_t> indeces2;
      std::vector<uint64_t> support2;
      double Isupp, Ifnew, Ifr, Ifc, Ifrc;  

      for( uint64_t r = 0; r < (X.size()-1) ; ++r )
      {
        for( uint64_t c = r+1; c < X.size() ; ++c )
        {
          indeces2 = {r,c};
          support2 = {original_support[r], original_support[c]};

          std::string Sr, Sc;
          uint64_t Sr_64t = original_support[r];
          uint64_t Sc_64t = original_support[c];
          Sr = std::to_string( Sr_64t );
          Sc = std::to_string( Sc_64t );
          std::string support_key = Sr + " " + Sc;
      
          Xtmp = { X[r], X[c] };
          auto tt_tmp = create_function( Xtmp, Y, indeces2 );

          if( (_Icoll.Frc).find( support_key ) == (_Icoll.Frc).end() )
          {
            Isupp = MI( {X[r], X[c]}, Y, support2 ); 
            Ifnew = MI( { Xtmp[Xtmp.size()-1] }, Y ); // support[k] -> k
            Ifr = MI( {  Xtmp[Xtmp.size()-1], X[r] }, Y );
            Ifc = MI( {  Xtmp[Xtmp.size()-1], X[c] }, Y );
            Ifrc = MI( {  Xtmp[Xtmp.size()-1], X[r], X[c] }, Y );
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

          if ( Ifnew > Imax )//( ( mi_supp == mi_Fnew ) && ( mi_supp > MImax ) )
          {
            Imax = Ifnew;
            nns.is_created = true;
            nns.tt = tt_tmp;
            nns.support = support2;
            nns.indeces = indeces2;
            nns.I = Ifnew;
            new_node = Xtmp[Xtmp.size()-1];

            if( ( Isupp == Ifnew ) && (Ifrc == Ifnew) && ( Ifr == Ifnew ) && ( Ifc == Ifnew) )
              nns.rc_del = true; 
            else
              nns.rc_del = false; 
          }
    }
  }
  // modify

  if( nns.is_created )
  {
    support.push_back(_num_nodes);
    create_klut_node( nns.support, nns.tt );
    add_tt_to_hash( nns.tt );
    X.push_back( new_node );
    std::cout <<  "created f(A[" << nns.indeces[0] << "],A[" << nns.indeces[1] << "])=f(" << nns.support[0] << "," << nns.support[1] << ")=" << nns.tt << std::endl;
    
    if ( nns.rc_del )
    {
      _cnt.Frc++;
      X.erase(X.begin()+std::max(nns.indeces[0], nns.indeces[1]) );
      X.erase(X.begin()+std::min(nns.indeces[0], nns.indeces[1]) );
      support.erase(support.begin()+std::max(nns.indeces[0], nns.indeces[1]));
      support.erase(support.begin()+std::min(nns.indeces[0], nns.indeces[1]));
    }
    else
    {
      _cnt.Fo++;
    }
  }
// END ################################
  return nns.is_created;
  }

    bool try_bottom_decomposition( std::vector<uint64_t>& support, dbitset_vector& X, dbitset_vector const& Y, double Imax )
    {
      new_nodes_storage nns;
      nns.is_created = false;

      std::vector<uint64_t> original_support = support;

      std::vector<uint64_t> indeces2;
      std::vector<uint64_t> support2;

      double Isupp, Ifnew, Ifr, Ifc, Ifrc;  
      dbitset_vector Xtmp;
      for( uint64_t r = 0; r < (X.size()-1) ; ++r )
      {
        for( uint64_t c = r+1; c < X.size() ; ++c )
        {
          indeces2 = {r,c};
          support2 = {original_support[r], original_support[c]};

          std::string Sr, Sc;
          uint64_t Sr_64t = original_support[r];
          uint64_t Sc_64t = original_support[c];
          Sr = std::to_string( Sr_64t );
          Sc = std::to_string( Sc_64t );
          std::string support_key = Sr + " " + Sc;

          Xtmp = {X[r], X[c]};
          auto tt = create_function( Xtmp, Y, indeces2 );

          if( (_Icoll.Frc).find( support_key ) == (_Icoll.Frc).end() )
          {
            Isupp = MI( {X[r], X[c]}, Y, support2 ); 
            Ifnew = MI( { X[X.size()-1] }, Y ); // support[k] -> k
            Ifr = MI( {  X[X.size()-1], X[r] }, Y );
            Ifc = MI( {  X[X.size()-1], X[c] }, Y );
            Ifrc = MI( {  X[X.size()-1], X[r], X[c] }, Y );
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
 

          if( ( Isupp == Ifnew ) && (Ifrc == Ifnew) && ( Ifr == Ifnew ) && ( Ifc == Ifnew) )
          {
            _cnt.Frc++;
            support.push_back(_num_nodes);
            create_klut_node( support2, tt );
            add_tt_to_hash( tt );
            X.push_back(Xtmp[Xtmp.size()-1]);
            X.erase(X.begin()+std::max(indeces2[0], indeces2[1]) );
            X.erase(X.begin()+std::min(indeces2[0], indeces2[1]) );
            support.erase(support.begin()+std::max(indeces2[0], indeces2[1]));
            support.erase(support.begin()+std::min(indeces2[0], indeces2[1]));
            return true;
          }
        }
      }
      return false;
    }

    bool try_bottom_decomposition_S(  std::vector<uint64_t>& support, dbitset_vector& X, dbitset_vector const& Y, double Imax,
                                      std::vector<double>& Ivect, std::vector<uint64_t>& IDXvect )
    {
      quicksort_by_attribute( IDXvect, Ivect,  0, Ivect.size()-1 );
      new_nodes_storage nns;
      nns.is_created = false;

      std::vector<uint64_t> original_support = support;

      std::vector<uint64_t> indeces2;
      std::vector<uint64_t> support2;

      double Isupp, Ifnew, Ifr, Ifc, Ifrc;  
      dbitset_vector Xtmp;
      //for( uint64_t r = 0; r < (X.size()-1) ; ++r )
      //{
        for( uint64_t i=0;  i< Ivect.size()-1 ; ++i )
        {
          auto r = IDXvect[i];
          auto c = IDXvect[i+1];
          indeces2 = {r,c};
          support2 = {original_support[r], original_support[c]};

          std::string Sr, Sc;
          uint64_t Sr_64t = original_support[r];
          uint64_t Sc_64t = original_support[c];
          Sr = std::to_string( Sr_64t );
          Sc = std::to_string( Sc_64t );
          std::string support_key = Sr + " " + Sc;

          Xtmp = {X[r], X[c]};
          auto tt = create_function( Xtmp, Y, indeces2 );

          if( (_Icoll.Frc).find( support_key ) == (_Icoll.Frc).end() )
          {
            Isupp = MI( {X[r], X[c]}, Y, support2 ); 
            Ifnew = MI( { X[X.size()-1] }, Y ); // support[k] -> k
            Ifr = MI( {  X[X.size()-1], X[r] }, Y );
            Ifc = MI( {  X[X.size()-1], X[c] }, Y );
            Ifrc = MI( {  X[X.size()-1], X[r], X[c] }, Y );
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
 

          if( ( Isupp == Ifnew ) && (Ifrc == Ifnew) && ( Ifr == Ifnew ) && ( Ifc == Ifnew) )
          {
            _cnt.Frc++;
            support.push_back(_num_nodes);
            create_klut_node( support2, tt );
            add_tt_to_hash( tt );
            X.push_back(Xtmp[Xtmp.size()-1]);
            X.erase(X.begin()+std::max(indeces2[0], indeces2[1]) );
            X.erase(X.begin()+std::min(indeces2[0], indeces2[1]) );
            support.erase(support.begin()+std::max(indeces2[0], indeces2[1]));
            support.erase(support.begin()+std::min(indeces2[0], indeces2[1]));
            return true;
          }
        }
      //}
      return false;
    }

    bool c2try_bottom_decomposition( std::vector<uint64_t>& support, 
                                    dbitset_vector& X, dbitset_vector& Y )
    {
      if(support.size()<3)
        return false;
      
      for( uint64_t j{1}; j < support.size(); ++j )
      {
        for( uint64_t i = 0; i < j; ++i )
        {

          auto XY0 = compute_cofactor( X, Y, i, 0 );
          auto XY00 = compute_cofactor( XY0.first, XY0.second, j-1, 0 );
          auto XY01 = compute_cofactor( XY0.first, XY0.second, j-1, 1 );

          auto XY1 = compute_cofactor( X, Y, i, 1 );
          auto XY10 = compute_cofactor( XY1.first, XY1.second, j-1, 0 );
          auto XY11 = compute_cofactor( XY1.first, XY1.second, j-1, 1 );

          if( (XY00.first.size() == 0)||(XY01.first.size() == 0)||(XY10.first.size() == 0)||(XY11.first.size() == 0) )
            return false;
          if( (XY00.first[0].size() == 1)||(XY01.first[0].size() == 1)||(XY10.first[0].size() == 1)||(XY11.first[0].size() == 1) )
            return false;
          std::vector<uint64_t> weq01 = {};
          std::vector<uint64_t> weq02 = {};
          std::vector<uint64_t> weq03 = {};
          std::vector<uint64_t> weq12 = {};
          std::vector<uint64_t> weq13 = {};
          std::vector<uint64_t> weq23 = {};

          uint32_t min_corr = 0;
          bool eq01 = is_F1_F0( XY00, XY01, weq01, min_corr ); // tells you where in XY01 there are repetitions
          bool eq02 = is_F1_F0( XY00, XY10, weq02, min_corr );
          bool eq03 = is_F1_F0( XY00, XY11, weq03, min_corr );
          bool eq12 = is_F1_F0( XY01, XY10, weq12, min_corr );
          bool eq13 = is_F1_F0( XY01, XY11, weq13, min_corr );
          bool eq23 = is_F1_F0( XY10, XY11, weq23, min_corr );

          auto num_pairs =  static_cast<uint32_t>(eq01)+
                            static_cast<uint32_t>(eq02)+
                            static_cast<uint32_t>(eq03)+
                            static_cast<uint32_t>(eq12)+
                            static_cast<uint32_t>(eq13)+
                            static_cast<uint32_t>(eq23);

          if( (num_pairs != 2) && (num_pairs != 3) )
            return false;
          if( eq12 && eq13 && eq23 ) // F00 is different
          {
            auto fxy = _klut.create_or( _itos.storage[support[i]], _itos.storage[support[j]] );
            _cntOR++;
            _itos.insert( _num_nodes, fxy );
            support.push_back(_num_nodes);
            _num_nodes++;
            support.erase(support.begin()+j);
            support.erase(support.begin()+i);
            X.push_back( X[i]|X[j] );
            X.erase(X.begin()+j);
            X.erase(X.begin()+i);
            return true;
          }
          else if( eq02 && eq03 && eq23 ) // F01 is different
          {
            auto fxy = _klut.create_lt( _itos.storage[support[i]], _itos.storage[support[j]] );
            _itos.insert( _num_nodes, fxy );
            support.push_back(_num_nodes);
            _num_nodes++;
            support.erase(support.begin()+j);
            support.erase(support.begin()+i);
            _cntLT++;
            X.push_back( (~X[i]) & X[j] );
            X.erase(X.begin()+j);
            X.erase(X.begin()+i);     
            return true;
          }
          else if( eq01 && eq03 && eq13 ) // F10 is different
          {
            auto fxy = _klut.create_le( _itos.storage[support[i]], _itos.storage[support[j]] );
            _cntLE++;
            _itos.insert( _num_nodes, fxy );
            support.push_back(_num_nodes);
            _num_nodes++;
            support.erase(support.begin()+j);
            support.erase(support.begin()+i);
            X.push_back( (~X[i]) | X[j] );
            X.erase(X.begin()+j);
            X.erase(X.begin()+i); 
            return true;
          }
          else if( eq01 && eq02 && eq12 ) // F11 is different
          {
            auto fxy = _klut.create_and( _itos.storage[support[i]], _itos.storage[support[j]] );
            _cntAND++;
            _itos.insert( _num_nodes, fxy );
            support.push_back(_num_nodes);
            _num_nodes++;
            support.erase(support.begin()+j);
            support.erase(support.begin()+i);
            X.push_back( X[i] & X[j] );
            X.erase(X.begin()+j);
            X.erase(X.begin()+i); 
            return true;
          }
          else if( eq03 && eq12 )
          {

            auto fxy = _klut.create_xor( _itos.storage[support[i]], _itos.storage[support[j]] );
            _cntXOR++;
            _itos.insert( _num_nodes, fxy );
            support.push_back(_num_nodes);
            _num_nodes++;
            support.erase(support.begin()+j);
            support.erase(support.begin()+i);
            X.push_back( X[i] ^ X[j] );
            X.erase(X.begin()+j);
            X.erase(X.begin()+i); 
            return true; 
          }
          else
          {
            return false;
          }
        }
      }
    }
    #pragma endregion

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

    #pragma region it-decomposition
    uint64_t idsd_step( std::vector<uint64_t> support, dbitset_vector& X, dbitset_vector& Y )
    {
      if( X.size() == 0 )
        return _klut.get_constant( false );
      assert( (support.size() == X.size()) );
  
      assert( (X[0].size() == Y[0].size()) );
      if( X[0].size() == 0 )
        return _klut.get_constant( false );
      
      if( Y[0].count() == 0 ) // tautology
        return _klut.get_constant( false );
      else if( Y[0].count() == Y[0].size() ) // contradiction
        return _klut.get_constant( true ); 
      
      if( support.size() <= _max_sup )
      {
        _cnt.Chj++;
        return create_klut_node( support, create_function( X, Y ) );
      }

      double Imax = 0;
      uint64_t idx_max = 0;
      double Inew;
      std::vector<double> Ivect;
      std::vector<uint64_t> IDXvect;
      if(_informed)
      {
        for( size_t i{0u}; i<support.size(); ++i )
        {
          Inew = MI( {X[i]}, Y, support );
          IDXvect.push_back(i);
          Ivect.push_back(Inew);
          bool is_new_better = ( Inew > Imax );
          Imax = is_new_better ? Inew : Imax;
          idx_max = is_new_better ? i : idx_max;
        }
      }
      auto XY0 = compute_cofactor( X, Y, idx_max, 0 );
      auto XY1 = compute_cofactor( X, Y, idx_max, 1 );

      std::vector<uint64_t> reduced_support = support;
      reduced_support.erase( reduced_support.begin() + idx_max );
      if(_top_decompose)
      {
        if( (XY0.second.size() != 0) && (XY0.second[0].count() == XY0.second[0].size()) ) // F0 = 1
        {
          _cnt.F0T++;
          _Icoll.clear();
          auto F1 = idsd_step( reduced_support, XY1.first, XY1.second );
          return _klut.create_le( _itos.storage[support[idx_max]], F1 );
        }
        else if( (XY0.second.size() == 0) || (XY0.second[0].count() == 0) ) // F0 = 0
        {
          _cnt.F0C++;
          _Icoll.clear();
          auto F1 = idsd_step( reduced_support, XY1.first, XY1.second );
          return _klut.create_and( _itos.storage[support[idx_max]], F1 );
        }
        else if( (XY1.second.size() != 0) && (XY1.second[0].count() == XY1.second[0].size()) ) // F1 = 1
        {
          _cnt.F1T++;
          _Icoll.clear();
          auto F0 = idsd_step( reduced_support, XY0.first, XY0.second );
          return _klut.create_or( _itos.storage[support[idx_max]], F0 );
        }
        else if( (XY1.second.size() == 0) || (XY1.second[0].count() == 0) ) // F1 = 0
        {
          _cnt.F1C++;
          _Icoll.clear();
          auto F0 = idsd_step( reduced_support, XY0.first, XY0.second );
          return _klut.create_lt( _itos.storage[support[idx_max]], F0 );
        }

        double P1 = (double)XY1.second[0].count()/XY1.second[0].size();
        double dP1 = 5*sqrt(P1*(1-P1))/(sqrt(XY1.second[0].size()));
        double P0 = (double)XY0.second[0].count()/XY0.second[0].size();
        double dP0 = 5*sqrt(P0*(1-P0))/(sqrt(XY0.second[0].size()));
        bool is_xor0 = ( ( P0-dP0 <= 0.5 ) && ( P0+dP0 >= 0.5 ) );
        bool is_xor1 = ( ( P1-dP1 <= 0.5 ) && ( P1+dP1 >= 0.5 ) );
        if ( is_xor1 && is_xor0 )
        double P1 = (double)Y[0].count()/Y[0].size();
            //3.89
        //double dP = 3.89*sqrt(P1*(1-P1))/(sqrt(Y[0].size()) );
        //if( ( P1-dP <= 0.5 ) && ( P1+dP >= 0.5 ) )
        //{
        //std::cout <<"Xs "<< X[0].size() << " XY0s" << XY0.first[0].size()<< " XY1s" << XY1.first[0].size()<<std::endl;
        if(  is_F1_not_F0( XY0, XY1 ) )
        {
          _cnt.XOR++;
          _Icoll.clear();
          auto pi_sig = _itos.storage[support[idx_max]];
          remove_column_and_invert( X, Y, idx_max ); // checked correct
          auto f0bar = idsd_step( reduced_support, X, Y );
          return _klut.create_xor( pi_sig , f0bar );
        }
      }
      if(_bottom_decompose)
      {
        if( try_bottom_decomposition_S( support, X, Y, Imax, Ivect, IDXvect ) ) //try_ME_step
          return idsd_step( support, X, Y );
      }

      _Icoll.clear();
      auto F0 = idsd_step( reduced_support, XY0.first, XY0.second );
      _Icoll.clear();
      auto F1 = idsd_step( reduced_support, XY1.first, XY1.second );

      auto f0 = _klut.create_and( _klut.create_not(_itos.storage[support[idx_max]]), F0 );
      auto f1 = _klut.create_and(_itos.storage[support[idx_max]], F1 );

      return _klut.create_or( f1, f0 );
    }
    
    bool simulate_input( dbitset const& input_pattern )
    {

      std::vector<bool> inpt_v;
      for( uint64_t k{0u}; k<input_pattern.size();++k )
      {
        inpt_v.push_back( ( ( input_pattern[k] == 1 ) ? true : false ) );
      }

      return simulate<bool>( _aig, default_simulator<bool>( inpt_v ) )[0];
    }

    double compute_accuracy( dbitset_vector const& X, dbitset_vector const& Y )
    {
      double acc = 0;
      double delta_acc;
      for( uint64_t k {0u}; k < X[0].size(); ++k )
      {
        dbitset ipattern;
        for ( uint64_t j {0u}; j < X.size(); ++j )
          ipattern.push_back( X[j][k] );
        
        delta_acc = ( ( simulate_input( ipattern ) == Y[0][k] ) ? (double)1.0/X[0].size() : 0.0 );
        acc += delta_acc;
      }
      return acc;
    }

    void ME(  dbitset_vector const& Xtrain, dbitset_vector const& Ytrain,
              dbitset_vector const& Xtest , dbitset_vector const& Ytest ,
              dbitset_vector const& Xvalid, dbitset_vector const& Yvalid )
    {
      double duration;
      std::clock_t start;
      start = std::clock();
      
      auto nodes = _nodes;
      auto outputs = _outputs;
      std::vector<uint64_t> support;
      for( size_t i {0}; i<nodes.size(); ++i )
        support.push_back(i);

      std::cout << "perform decomposition " << std::endl;
      _klut.create_po( idsd_step( support, nodes, outputs ) );
      std::cout << "convert klut to aig " << std::endl;
      _aig = convert_klut_to_graph<aig_network>( _klut );
      std::cout << "compute train accuracy " << std::endl;
      _cnt.train_acc = compute_accuracy( Xtrain, Ytrain );
      std::cout << "compute test accuracy " << std::endl;
      _cnt.test_acc  = compute_accuracy( Xtest, Ytest );
      std::cout << "compute valid accuracy " << std::endl;
      _cnt.valid_acc = compute_accuracy( Xvalid, Yvalid );
      duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
      _duration = duration;
      print_features();

    }
    #pragma endregion

    #pragma region decomposition

    void combine_covers( std::pair<dbitset_vector, dbitset_vector>& XY, std::vector<std::pair<dbitset_vector, dbitset_vector>> XYs )
    {
      if( XYs.size() != 0 )
      {
      for( size_t i {0}; i < XYs.size(); ++i ) // all covers one by one
      {
        if( XYs[i].first.size() != 0)
        {
          for( size_t k{0}; k < XYs[i].first.size(); ++k ) // all nodes in each cover
          {
            for( size_t j{0}; j < XYs[i].first[k].size(); ++j ) // all elements in the cover
            {
              XY.first[k].push_back( XYs[i].first[k][j] );
            }
          }
          for( size_t j{0}; j < XYs[i].second[0].size(); ++j )
          {
            XY.second[0].push_back( XYs[i].second[0][j] );
          }
        }
      }
      }

    }

    bool erase_from_cover( std::pair<dbitset_vector, dbitset_vector>& XY, std::vector<uint64_t> where )
    {
      if( where.size()==0 )
        return true;
        
      bool is_first = true;
      bool is_new = false;
      std::pair<dbitset_vector, dbitset_vector> XYnew;
      XYnew.first = {};
      XYnew.second = {};
      for( uint64_t i{0}; i < XY.first[0].size(); ++i )
      {
        if( std::find (where.begin(), where.end(), i) == where.end() )
        {
          if( is_first )
          {
            dbitset B;
            B.push_back( XY.second[0][i] );
            XYnew.second.push_back( B );
            for( auto j{0}; j<XY.first.size(); ++j )
            {
              dbitset A;
              A.push_back(XY.first[j][i]);
              XYnew.first.push_back( A );
            }
            is_first = false;
          }
          else
          {
            XYnew.second[0].push_back( XY.second[0][i]) ;
            for( auto j{0}; j<XY.first.size(); ++j )
            {
              XYnew.first[j].push_back( XY.first[j][i] );
            }
          }
        }
      }
      if( !is_first ) // at least one removed
        XY = XYnew;

      return ( !is_first );
      
    }

    bool ctry_bottom_decomposition( std::vector<uint64_t>& support, 
                                    dbitset_vector& X, dbitset_vector& Y, uint64_t& new_signal )
    {
      if(support.size()<3)
        return false;
      for( uint64_t j{1}; j < support.size(); ++j )
      {
        for( uint64_t i = 0; i < j; ++i )
        {
          auto XY0 = compute_cofactor( X, Y, i, 0 );

          auto XY00 = compute_cofactor( XY0.first, XY0.second, j-1, 0 );
          auto XY01 = compute_cofactor( XY0.first, XY0.second, j-1, 1 );
        //  std::cout << "b1" << std::endl;

          auto XY1 = compute_cofactor( X, Y, i, 1 );
      //    std::cout << "b2" << std::endl;

          auto XY10 = compute_cofactor( XY1.first, XY1.second, j-1, 0 );
          auto XY11 = compute_cofactor( XY1.first, XY1.second, j-1, 1 );
        //  std::cout << "b3" << std::endl;

          if( (XY00.first.size() == 0)||(XY01.first.size() == 0)||(XY10.first.size() == 0)||(XY11.first.size() == 0) )
            return false;
          if( (XY00.first[0].size() == 1)||(XY01.first[0].size() == 1)||(XY10.first[0].size() == 1)||(XY11.first[0].size() == 1) )
            return false;
          //std::cout << "b4" << std::endl;
          
          std::vector<uint64_t> weq01 = {};
          std::vector<uint64_t> weq02 = {};
          std::vector<uint64_t> weq03 = {};
          std::vector<uint64_t> weq12 = {};
          std::vector<uint64_t> weq13 = {};
          std::vector<uint64_t> weq23 = {};

          uint32_t min_corr = 0;
          bool eq01 = is_F1_F0( XY00, XY01, weq01, min_corr ); // tells you where in XY01 there are repetitions
          bool eq02 = is_F1_F0( XY00, XY10, weq02, min_corr );
          bool eq03 = is_F1_F0( XY00, XY11, weq03, min_corr );
          bool eq12 = is_F1_F0( XY01, XY10, weq12, min_corr );
          bool eq13 = is_F1_F0( XY01, XY11, weq13, min_corr );
          bool eq23 = is_F1_F0( XY10, XY11, weq23, min_corr );

         // std::cout << "b5" << std::endl;


          auto num_pairs =  static_cast<uint32_t>(eq01)+
                            static_cast<uint32_t>(eq02)+
                            static_cast<uint32_t>(eq03)+
                            static_cast<uint32_t>(eq12)+
                            static_cast<uint32_t>(eq13)+
                            static_cast<uint32_t>(eq23);
          //std::cout << "eq01=" << eq01 << " eq02=" << eq02 <<  " eq03=" << eq03 <<  " eq12=" << eq12 <<  " eq13=" << eq13 <<  " eq23=" << eq23 << std::endl; 
          //std::cout << "a1" << std::endl;

          if( (num_pairs != 2) && (num_pairs != 3) )
            return false;
          //std::cout << "num_pairs " << num_pairs << std::endl;

          if( eq12 && eq13 && eq23 ) // F00 is different
          {
          //std::cout << "a2" << std::endl;

            auto fxy = _klut.create_or( _itos.storage[support[i]], _itos.storage[support[j]] );
            _itos.insert( _num_nodes, fxy );
            _num_nodes++;
            support.erase(support.begin()+j);
            support.erase(support.begin()+i);
            
            std::vector<std::pair<dbitset_vector, dbitset_vector>> XYs;// = { XY01, XY10 };
            if( erase_from_cover( XY10, weq12 ))
            {
              XYs.push_back( XY10 );
            }
            if( erase_from_cover( XY11, weq23 ) )
            {
              XYs.push_back( XY11 );
            }
            combine_covers( XY01, XYs );
            auto F11 = cdsd_step( support, XY01.first, XY01.second );
            auto F00 = cdsd_step( support, XY00.first, XY00.second );
            _cntOR++;
            auto f1 = _klut.create_and( fxy, F11 );
            auto f0 = _klut.create_lt( fxy, F00 );
            new_signal = _klut.create_or( f1, f0 );
            return true;
          }
          else if( eq02 && eq03 && eq23 ) // F01 is different
          {
          //std::cout << "a3" << std::endl;

            auto fxy = _klut.create_lt( _itos.storage[support[i]], _itos.storage[support[j]] );
            _itos.insert( _num_nodes, fxy );
            _num_nodes++;
            support.erase(support.begin()+j);
            support.erase(support.begin()+i);
            _cntLT++;
            std::vector<std::pair<dbitset_vector, dbitset_vector>> XYs;// = { XY01, XY10 };
            if( erase_from_cover( XY10, weq02 ) )
              XYs.push_back( XY10 );
            
            for(auto k {0}; k < weq03.size(); ++k )
              weq23.push_back( weq03[k] );
            if( erase_from_cover( XY11, weq23 ) )
              XYs.push_back( XY11 );

            combine_covers( XY00, XYs );

            auto F10 = cdsd_step( support, XY00.first, XY00.second );
            auto F01 = cdsd_step( support, XY01.first, XY01.second );
           
            auto f1 = _klut.create_and( fxy, F01 );
            auto f0 = _klut.create_lt( fxy, F10 );
            new_signal = _klut.create_or( f1, f0 );

            return true;
          }
          else if( eq01 && eq03 && eq13 ) // F10 is different
          {
//          std::cout << "a4" << std::endl;

            auto fxy = _klut.create_le( _itos.storage[support[i]], _itos.storage[support[j]] );
            _cntLE++;
            _itos.insert( _num_nodes, fxy );
            _num_nodes++;
            support.erase(support.begin()+j);
            support.erase(support.begin()+i);

            std::vector<std::pair<dbitset_vector, dbitset_vector>> XYs;// = { XY01, XY10 };
            if( erase_from_cover( XY01, weq01 ) )
              XYs.push_back( XY01 );
            
            for(auto k {0}; k < weq03.size(); ++k )
              weq13.push_back( weq03[k] );
            if( erase_from_cover( XY11, weq13 ) )
              XYs.push_back( XY11 );
            combine_covers( XY00, XYs );

            auto F01 = cdsd_step( support, XY00.first, XY00.second );
            auto F10 = cdsd_step( support, XY10.first, XY10.second );
           
            auto f1 = _klut.create_and( fxy, F01 );
            auto f0 = _klut.create_lt( fxy, F10 );
            new_signal = _klut.create_or( f1, f0 );
            return true;
          }
          else if( eq01 && eq02 && eq12 ) // F11 is different
          {
  //        std::cout << "a5" << std::endl;

            auto fxy = _klut.create_and( _itos.storage[support[i]], _itos.storage[support[j]] );
            _cntAND++;
            _itos.insert( _num_nodes, fxy );
            _num_nodes++;
            support.erase(support.begin()+j);
            support.erase(support.begin()+i);

            std::vector<std::pair<dbitset_vector, dbitset_vector>> XYs;// = { XY01, XY10 };
            if( erase_from_cover( XY01, weq01 ) )
              XYs.push_back( XY01 );
            
            for(auto k {0}; k < weq02.size(); ++k )
              weq12.push_back( weq02[k] );
            if( erase_from_cover( XY10, weq12 ) )
              XYs.push_back( XY10 );
            combine_covers( XY00, XYs );

            auto F00 = cdsd_step( support, XY00.first, XY00.second );
            auto F11 = cdsd_step( support, XY11.first, XY11.second );
            auto f1 = _klut.create_and( fxy, F11 );
            auto f0 = _klut.create_lt( fxy, F00 );
            new_signal = _klut.create_or( f1, f0 );

            return true;
          }
          else if( eq03 && eq12 )
          {
    //      std::cout << "a6" << std::endl;

            auto fxy = _klut.create_xor( _itos.storage[support[i]], _itos.storage[support[j]] );
            _cntXOR++;
            _itos.insert( _num_nodes, fxy );
            _num_nodes++;
            support.erase(support.begin()+j);
            support.erase(support.begin()+i);

            //std::vector<std::pair<dbitset_vector, dbitset_vector>> XYs = {XY11};
            if( erase_from_cover( XY11, weq03 ) )
              combine_covers( XY00, {XY11} );
            if( erase_from_cover( XY10, weq12 ) )
              combine_covers( XY01, {XY10} );
            
            auto F00 = cdsd_step( support, XY00.first, XY00.second );
            auto F01 = cdsd_step( support, XY01.first, XY01.second );
           
            auto f1 = _klut.create_and( fxy, F01 );
            auto f0 = _klut.create_lt( fxy, F00 );
            new_signal = _klut.create_or( f1, f0 );

            return true; 
          }
          else
          {
            return false;
          }
        }
      }
    }


    uint64_t cdsd_step( std::vector<uint64_t> support, dbitset_vector& X, dbitset_vector& Y )
    {
      //std::cout << 1 << std::endl;
      if( X.size() == 0 )
        return _klut.get_constant( false );
      assert( (support.size() == X.size()) );
      //std::cout << 2 << std::endl;
  
      assert( (X[0].size() == Y[0].size()) );
      if( X[0].size() == 0 )
        return _klut.get_constant( false );
      //std::cout << 3 << std::endl;
      
      if( Y[0].count() == 0 ) // tautology
        return _klut.get_constant( false );
      else if( Y[0].count() == Y[0].size() ) // contradiction
        return _klut.get_constant( true ); 
      //std::cout << 3 << std::endl;
      
      if( support.size() <= 1 )
      {
        _cnt.Chj++;
        return create_klut_node( support, create_function( X, Y ) );
      }

      //std::cout << 4 << std::endl;

      std::vector<uint64_t> reduced_support;
      double Imax = 0;
      uint64_t idx_max = 0;
      double Inew;
      std::vector<double> Ivect;
      std::vector<uint64_t> IDXvect;
      for( size_t i{0u}; i<support.size(); ++i )
      {
        Inew = MI( {X[i]}, Y, support );
        IDXvect.push_back(i);
        Ivect.push_back(Inew);
        bool is_new_better = ( Inew > Imax );
        Imax = is_new_better ? Inew : Imax;
        idx_max = is_new_better ? i : idx_max;
      } 
      //std::cout << 5 << std::endl;

      //for( size_t i{0}; i < support.size(); ++i )
      //{
        reduced_support = support;
        auto XY0 = compute_cofactor( X, Y, idx_max, 0 );
        auto XY1 = compute_cofactor( X, Y, idx_max, 1 );

        reduced_support.erase( reduced_support.begin() + idx_max );
      //std::cout << 6 << std::endl;

        if( (XY0.second.size() != 0) && (XY0.second[0].count() == XY0.second[0].size()) ) // F0 = 1
        {
          _cnt.F0T++;
          _Icoll.clear();
          auto F1 = cdsd_step( reduced_support, XY1.first, XY1.second );
          return _klut.create_le( _itos.storage[support[idx_max]], F1 );
        }
        else if( (XY0.second.size() == 0) || (XY0.second[0].count() == 0) ) // F0 = 0
        {
          _cnt.F0C++;
          _Icoll.clear();
          auto F1 = cdsd_step( reduced_support, XY1.first, XY1.second );
          return _klut.create_and( _itos.storage[support[idx_max]], F1 );
        }
        else if( (XY1.second.size() != 0) && (XY1.second[0].count() == XY1.second[0].size()) ) // F1 = 1
        {
          _cnt.F1T++;
          _Icoll.clear();
          auto F0 = cdsd_step( reduced_support, XY0.first, XY0.second );
          return _klut.create_or( _itos.storage[support[idx_max]], F0 );
        }
        else if( (XY1.second.size() == 0) || (XY1.second[0].count() == 0) ) // F1 = 0
        {
          _cnt.F1C++;
          _Icoll.clear();
          auto F0 = cdsd_step( reduced_support, XY0.first, XY0.second );
          return _klut.create_lt( _itos.storage[support[idx_max]], F0 );
        }
        if( is_F1_not_F0( XY0, XY1, 0 ) )
        {
          _cnt.XOR++;
          _Icoll.clear();
          auto pi_sig = _itos.storage[support[idx_max]];
          remove_column_and_invert( X, Y, idx_max ); // checked correct
          auto f0bar = cdsd_step( reduced_support, X, Y );
          return _klut.create_xor( pi_sig , f0bar );
        }
      //}
      uint64_t new_signal;
      //std::cout << "A1" << std::endl;
      //if( ctry_bottom_decomposition( support, X, Y, new_signal ) )//try_ME_step
      //  return new_signal;
      //std::cout << 7 << std::endl;

      if( c2try_bottom_decomposition( support, X, Y ) )//try_ME_step
        return cdsd_step( support, X, Y );
      //std::cout << 8 << std::endl;


      //auto XY0 = compute_cofactor( X, Y, idx_max, 0 );
      //auto XY1 = compute_cofactor( X, Y, idx_max, 1 );
      _Icoll.clear();

      reduced_support = support;
      reduced_support.erase(reduced_support.begin()+idx_max);
      auto F0 = cdsd_step( reduced_support, XY0.first, XY0.second );
      _Icoll.clear();
      auto F1 = cdsd_step( reduced_support, XY1.first, XY1.second );

      auto f0 = _klut.create_and( _klut.create_not(_itos.storage[support[idx_max]]), F0 );
      auto f1 = _klut.create_and(_itos.storage[support[idx_max]], F1 );

      return _klut.create_or( f1, f0 );
    }

    void cdsd(  dbitset_vector const& Xtrain, dbitset_vector const& Ytrain,
              dbitset_vector const& Xtest , dbitset_vector const& Ytest ,
              dbitset_vector const& Xvalid, dbitset_vector const& Yvalid )
    {
      double duration;
      std::clock_t start;
      start = std::clock();
      
      auto nodes = _nodes;
      auto outputs = _outputs;
      std::vector<uint64_t> support;
      for( size_t i {0}; i<nodes.size(); ++i )
        support.push_back(i);

      std::cout << "perform decomposition " << std::endl;
      _klut.create_po( cdsd_step( support, nodes, outputs ) );
      std::cout << "convert klut to aig " << std::endl;
      _aig = convert_klut_to_graph<aig_network>( _klut );
      std::cout << "compute train accuracy " << std::endl;
      _cnt.train_acc = compute_accuracy( Xtrain, Ytrain );
      std::cout << "compute test accuracy " << std::endl;
      _cnt.test_acc  = compute_accuracy( Xtest, Ytest );
      std::cout << "compute valid accuracy " << std::endl;
      _cnt.valid_acc = compute_accuracy( Xvalid, Yvalid );
      duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
      _duration = duration;
      print_features();

    }
    #pragma endregion

    public:
    dbitset_vector _nodes;   /* storage element: value of the output at each example */
    dbitset_vector _outputs; /* storage element: value of the output at each example */
    uint64_t _num_nodes;
    index_to_signal _itos;  
    klut_network _klut;
    aig_network _aig;
    uint64_t _max_sup;
    calls_counter _cnt;
    std::unordered_map<std::string,double> _HTX;
    std::unordered_map<std::string,double> _HTY;
    std::unordered_map<std::string,double> _HTXY;
    std::unordered_map<std::string, uint64_t> _tt_counter;
    bool _has_file = false;
    std::string _pathTOfile;
    std::string _ID = "";
    std::string _IDs = "";
    double _duration;
    Istorage _Icoll;
    uint32_t _cntOR = 0;
    uint32_t _cntLT = 0;
    uint32_t _cntLE = 0;
    uint32_t _cntAND = 0;
    uint32_t _cntXOR = 0;
    bool _top_decompose = false;
    bool _bottom_decompose = false;
    bool _dontknows = false;
    bool _informed = false;

  };

}