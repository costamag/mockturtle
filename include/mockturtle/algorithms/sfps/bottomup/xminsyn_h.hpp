/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
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

#include "../../../traits.hpp"

#include <fmt/format.h>
#include <kitty/constructors.hpp>
#include <kitty/decomposition.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operators.hpp>
#include <kitty/print.hpp>

namespace mockturtle
{

/*! \brief Parameters for dsd_decomposition */
struct xminsyn_h_params
{
  /*! \brief Apply XOR decomposition. */
};

namespace detail
{

template<class Ntk, class TT>
class xminsyn_h_impl
{
public:
  xminsyn_h_impl( Ntk& ntk, TT const& func, std::vector<signal<Ntk>> const& children, xminsyn_h_params const& ps )
      : _ntk( ntk ),
        remainder( func ),
        pis( children ),
        _ps( ps )
  {
    mask = ~remainder.construct();
    for ( uint32_t i = 0u; i < func.num_vars(); ++i )
    {
      if ( kitty::has_var( func, i ) )
      {
        support.push_back( i );
	      TT ipattern = func.construct();
	      kitty::create_nth_var( ipattern, i );
	      X.push_back( ipattern );
      }
    }
  }

private:
  enum class symmetry_types
  {
    NONE,
    ES_,
    NES_,
    MS_,
    SVS0X_,
    SVS1X_,
    SVSX0_,
    SVSX1_,
    CSVS00_,
    CSVS01_,
    CSVS10_,
    CSVS11_
  };

  std::string type2string( symmetry_types type )
  {
    switch( type )
    {
      case symmetry_types::ES_: return "ES";
      case symmetry_types::NES_: return "NES";
      case symmetry_types::MS_: return "MS";
      case symmetry_types::SVS0X_: return "SVS";
      case symmetry_types::SVS1X_: return "SVS";
      case symmetry_types::SVSX0_: return "SVS";
      case symmetry_types::SVSX1_: return "SVS";
      case symmetry_types::CSVS00_: return "CSVS";
      case symmetry_types::CSVS01_: return "CSVS";
      case symmetry_types::CSVS10_: return "CSVS";      
      case symmetry_types::CSVS11_: return "CSVS";
      default: return "NONE";
    }
  }


  struct symmetry_info_t{
    symmetry_types type { symmetry_types::NONE };
    uint32_t var1;
    uint32_t var2;
    uint32_t id{0}; 
    TT func; 
    TT mask;
    };

  TT cofactorA( TT fn, uint32_t A, uint32_t i, uint32_t j )
  {
    switch ( A )
    {
      case 0u: return cofactor0(cofactor0( fn, j), i );
      case 1u: return cofactor0(cofactor1( fn, j), i );
      case 2u: return cofactor1(cofactor0( fn, j), i );
      case 3u: return cofactor1(cofactor1( fn, j), i );
    }
  }

  TT cube_generator( uint32_t cube, uint32_t i, uint32_t j )
  {
    switch ( cube )
    {
    case 0: 
      return ~X[i] & ~X[j];
      break;
    case 1: 
      return ~X[i] & X[j];
      break;
    case 2: 
      return X[i] & ~X[j];
      break;
    case 3: 
      return X[i] & X[j];
      break;
    }
  }

  symmetry_info_t simple_remapping( uint32_t from, uint32_t to, uint32_t i, uint32_t j, symmetry_types type )
  {
    symmetry_info_t res;
    res.type = type;
    res.var1 = i;
    res.var2 = j; 
    res.id = 1u*( from > to );

    TT A = cube_generator( from, i, j );
    TT B = cube_generator( to, i, j  );

    res.mask = ( ~A & mask ) | ( B & cofactorA( mask, from, i, j ) );
    res.func = ( ~B & remainder ) | ( B & cofactorA( remainder & mask, from, i, j ) );
    
    return res;
  }



  symmetry_info_t multiform_remapping( uint32_t from1, uint32_t i, uint32_t j, symmetry_types type )
  {
    symmetry_info_t res;
    res.type = type;
    res.var1 = i;
    res.var2 = j; 
    res.id = from1;
    
    TT A = cube_generator( from1, i, j  );
    TT C = cube_generator( uint32_t( 3 - from1 ), i, j  );

    uint32_t from2 = ( from1 == 0u )*2u + ( from1 == 2 )*3u + ( from1 == 3u )*1u;

    TT B = cube_generator( from2, i, j  );
    TT D = cube_generator( uint32_t( 3 - from2 ), i, j  );

    res.mask = ( ~B & ~A & mask ) | ( C & cofactorA( mask, from1, i, j ) ) | ( D & cofactorA( mask, from2, i, j ) );
    res.func = ( ~( C | D ) & remainder ) | ( C & cofactorA( mask & remainder, from1, i, j ) ) | ( D & cofactorA( mask & remainder, from2, i, j ) );
    
    return res;
  }

  symmetry_info_t compatible_remapping( uint32_t from1, uint32_t from2, uint32_t to, uint32_t i, uint32_t j, symmetry_types type )
  {
    symmetry_info_t res;
    res.type = type;
    res.var1 = i;
    res.var2 = j; 
    res.id = 1u*( from1 > from2 );
    uint32_t excluded = uint32_t( 6u - from1 - from2 - to );

    TT A = cube_generator( from1, i, j );
    TT B = cube_generator( from2, i, j );
    TT C = cube_generator( to, i, j );
    TT D = cube_generator( excluded, i, j );

    res.mask = ( ~A & ~B & mask ) | ( C & cofactorA( mask, from1, i, j ) ) | ( C & cofactorA( mask, from2, i, j ) ) ;
    res.func = ( ~C & remainder ) | ( C & cofactorA( remainder & mask, from1, i, j ) ) | ( C & cofactorA( remainder & mask, from1, i, j ) );
    
    return res;
  }

  std::vector<symmetry_info_t> check_symmetry_type( TT tt, uint32_t i, uint32_t j )
  {
    std::vector<symmetry_info_t> res;
    symmetry_info_t sym_info;
    sym_info.var1 = i;
    sym_info.var2 = j;
    const auto tt0 = cofactor0( tt, i );
    const auto tt1 = cofactor1( tt, i );

    const auto tt00 = cofactor0( tt0, j );
    const auto tt01 = cofactor1( tt0, j );
    const auto tt10 = cofactor0( tt1, j );
    const auto tt11 = cofactor1( tt1, j );

    const auto eq01 = equal( mask&tt00, mask&tt01 );
    const auto eq02 = equal( mask&tt00, mask&tt10 );
    const auto eq03 = equal( mask&tt00, mask&tt11 );
    const auto eq12 = equal( mask&tt01, mask&tt10 );
    const auto eq13 = equal( mask&tt01, mask&tt11 );
    const auto eq23 = equal( mask&tt10, mask&tt11 );

    const auto num_pairs =
      static_cast<uint32_t>( eq01 ) +
      static_cast<uint32_t>( eq02 ) +
      static_cast<uint32_t>( eq03 ) +
      static_cast<uint32_t>( eq12 ) +
      static_cast<uint32_t>( eq13 ) +
      static_cast<uint32_t>( eq23 );
    
    
    if ( num_pairs == 0  )
    {
      sym_info.type = symmetry_types::NONE;
      res.push_back( sym_info );
      return res;
    }
    TT A, B, C, D;
    if ( eq12 ) // F01 = F10 NES
    {
      res.push_back( simple_remapping( 1, 2, i, j, symmetry_types::NES_ ) );
      res.push_back( simple_remapping( 2, 1, i, j, symmetry_types::NES_ ) );
    }
    if ( eq03 ) // F00 = F11 ES
    {
      res.push_back( simple_remapping( 0, 3, i, j, symmetry_types::ES_ ) );
      res.push_back( simple_remapping( 3, 0, i, j, symmetry_types::ES_ ) );
    }
    if ( eq02 ) // F10=F00
    {
      res.push_back( simple_remapping( 0, 2, i, j, symmetry_types::SVSX0_ ) );
      res.push_back( simple_remapping( 2, 0, i, j, symmetry_types::SVSX0_ ) );    
    }
    if ( eq13 ) // F11=F01
    {
      res.push_back( simple_remapping( 1, 3, i, j, symmetry_types::SVSX1_ ) );
      res.push_back( simple_remapping( 3, 1, i, j, symmetry_types::SVSX1_ ) ); 
    }
    if ( eq23 ) // F10=F11
    {
      res.push_back( simple_remapping( 2, 3, i, j, symmetry_types::SVS1X_ ) );
      res.push_back( simple_remapping( 3, 2, i, j, symmetry_types::SVS1X_ ) );
    }
    if ( eq01 ) // F00=F01
    {
      res.push_back( simple_remapping( 0, 1, i, j, symmetry_types::SVS0X_ ) );
      res.push_back( simple_remapping( 1, 0, i, j, symmetry_types::SVS0X_ ) );
    }
    if ( eq12 && eq03 ) // F01=F10 and F00=F11
    {
      res.push_back( multiform_remapping( 0, i, j, symmetry_types::MS_ ) );
      res.push_back( multiform_remapping( 1, i, j, symmetry_types::MS_ ) );      
      res.push_back( multiform_remapping( 2, i, j, symmetry_types::MS_ ) );      
      res.push_back( multiform_remapping( 3, i, j, symmetry_types::MS_ ) );      
    }
    if( eq02 && eq01 && eq12 )
    {
      res.push_back( compatible_remapping( 0, 1, 2, i, j, symmetry_types::CSVS00_ ) );
      res.push_back( compatible_remapping( 2, 0, 1, i, j, symmetry_types::CSVS00_ ) );
    }
    if( eq02 && eq23 && eq03 )
    {
      res.push_back( compatible_remapping( 3, 2, 0, i, j, symmetry_types::CSVS10_ ) );
      res.push_back( compatible_remapping( 2, 0, 3, i, j, symmetry_types::CSVS10_ ) );
    }
    if( eq13 && eq01 && eq03 )
    {
      res.push_back( compatible_remapping( 0, 1, 3, i, j, symmetry_types::CSVS01_ ) );
      res.push_back( compatible_remapping( 1, 3, 0, i, j, symmetry_types::CSVS01_ ) );
    }
    if( eq13 && eq23 && eq12 )
    {
      res.push_back( compatible_remapping( 1, 3, 2, i, j, symmetry_types::CSVS11_ ) );
      res.push_back( compatible_remapping( 3, 2, 1, i, j, symmetry_types::CSVS11_ ) );
    }
    return res;
  }

  void remap( symmetry_info_t symmetry )
  {

    switch( symmetry.type )
    {
      case symmetry_types::NES_:
        signal<Ntk> sig_i, sig_j;
        if( symmetry.id == 0u )
        {
          sig_i = ntk_.create_or( pis[symmetry.var1], pis[symmetry.var2] );
          sig_j = ntk_.create_and( pis[symmetry.var1], pis[symmetry.var2] );
        }
        else
        {
          signal<Ntk> sig_i = ntk_.create_and( pis[symmetry.var1], pis[symmetry.var2] );
          signal<Ntk> sig_j = ntk_.create_or( pis[symmetry.var1], pis[symmetry.var2] );
        }
        pis[symmetry.var1] = sig_i;
        pis[symmetry.var2] = sig_j;
        break; 

      case symmetry_types::ES_:
        signal<Ntk> sig_i, sig_j;
        if( symmetry.id == 0u )
        {
          sig_i = ntk_.create_nand( !pis[symmetry.var1], pis[symmetry.var2] );
          sig_j = ntk_.create_nand( pis[symmetry.var1], !pis[symmetry.var2] );
        }
        else
        {
          signal<Ntk> sig_i = ntk_.create_and( pis[symmetry.var1], !pis[symmetry.var2] );
          signal<Ntk> sig_j = ntk_.create_and( !pis[symmetry.var1], pis[symmetry.var2] );
        }
        pis[symmetry.var1] = sig_i;
        pis[symmetry.var2] = sig_j;
        break;      
    }
  } 

public:
  //signal<Ntk> run()
  void run()
  {
    while( support.size() > 1 )
    {
      auto km_tt = kitty::karnaugh_map( remainder );
      km_tt.print(mask);

      /* detect symmetries */
      std::vector<symmetry_info_t> symmetries;
      uint32_t k = 0;
      uint32_t choice;
      std::string order;

      for ( auto j = 1u; j < support.size(); ++j )
      {
        for ( auto i = 0u; i < j; ++i )
        {
	        std::vector<symmetry_info_t> sym_v = check_symmetry_type( remainder, support[i], support[j] );

	        if( !(sym_v.size()==1u && sym_v[0].type== symmetry_types::NONE ))
	        {
	          for( auto s : sym_v )
	          {
              std::cout << k++ << " " << type2string(s.type) << " " << s.var1 << " " << s.var2 << " |DC|=" << kitty::count_ones( ~s.mask ) << std::endl;
              symmetries.push_back( s );
	          }
	        }
        }
      }
      std::cout << "Choose the symmetry to exploit: " << std::endl;
      std::cin >> choice;
      std::cout << "Remapping " << " " << type2string(symmetries[choice].type) << " " << symmetries[choice].var1 << " " << symmetries[choice].var2 << std::endl;
      remainder = symmetries[choice].func;
      mask = symmetries[choice].mask;
      remap( symmetries[choice] );
    }
  }

private:
  Ntk& _ntk;
  TT remainder;
  TT mask;
  std::vector<uint32_t> support;
  std::vector<signal<Ntk>> pis;
  std::vector<TT> X;
  xminsyn_h_params const& _ps;
};

} // namespace detail

/*! \brief crossing minimization synthesis engine
 *
 * This function heuristically minimizes the number of crossings by using a 
 * symmetry-based synthesis heuristic.
 */
template<class Ntk, class TT>
void xminsyn_h( Ntk& ntk, TT const& func, std::vector<signal<Ntk>> const& children, xminsyn_h_params const& ps = {} )
{
  detail::xminsyn_h_impl impl( ntk, func, children, ps );
  impl.run();
}

} // namespace mockturtle
