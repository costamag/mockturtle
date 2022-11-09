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
 * THE SOFTWARE IS PROVid_ordED "AS IS", WITHOUT WARRANTY OF ANY KIND,
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
  bool try_top_decomposition{false};
};

namespace detail
{

template<class Ntk, class TT>
class xminsyn_h_impl
{
public:
  xminsyn_h_impl( Ntk& ntk, TT const& func, std::vector<signal<Ntk>> const& children, xminsyn_h_params const& ps )
      : ntk_( ntk ),
        remainder( func ),
        target( func ),
        pis( children ),
        _ps( ps )
  {
    mask = ~remainder.construct();
    for ( uint32_t i = 0u; i < func.num_vars(); ++i )
    {
      //if ( kitty::has_var( func, i ) )
      //{
      support.push_back( i );
	    TT ipattern = func.construct();
	    kitty::create_nth_var( ipattern, i );
	    X.push_back( ipattern );
      //}
    }
    initialize_gate_library();
  }

private:

  #pragma region decomposition 

  enum class decomposition_types
  {
    NONE,
    OR_,
    AND_,
    LT_,
    LE_,
    XOR_,
    OR2_,
    AND2_
  };

  struct decomposition_info_t{
    decomposition_types type { decomposition_types::NONE };
    uint32_t i;
    uint32_t j;
    uint32_t id{0}; 
    TT func; 
    TT mask;
    };

  void print_decomposition( decomposition_info_t s )
  {
    std::string remap_str = "";
    std::string i_str = std::to_string(s.i);
    std::string j_str = std::to_string(s.j);
    i_str.resize(2, ' ');
    j_str.resize(2, ' ');

    switch( s.type )
    {
      case decomposition_types::OR_:
        remap_str += "    X " + i_str + " ";
        remap_str.resize(20, ' ');
        remap_str = "[ " + i_str + " OR F0" + "]:" + remap_str;
        std::cout << remap_str ;
        break;
      case decomposition_types::LE_:
        remap_str += "    X " + i_str + " ";
        remap_str.resize(20, ' ');
        remap_str = "[ " + i_str + " LE F1" + "]:" + remap_str;
        std::cout << remap_str ;
        break;
      case decomposition_types::LT_:
        remap_str += "    X " + i_str + " ";
        remap_str.resize(20, ' ');
        remap_str = "[ " + i_str + " LT F0" + "]:" + remap_str;
        std::cout << remap_str ;
        break;
      case decomposition_types::AND_:
        remap_str += "    X " + i_str + " ";
        remap_str.resize(20, ' ');
        remap_str = "[ " + i_str + "AND F1" + "]:" + remap_str;
        std::cout << remap_str ;
        break;
      case decomposition_types::XOR_:
        remap_str += "    X " + i_str + " ";
        remap_str.resize(20, ' ');
        remap_str = "[ " + i_str + "XOR F0" + "]:" + remap_str;
        std::cout << remap_str ;
        break;
      case decomposition_types::OR2_:
        if( s.id == 0u )
        {
          remap_str += "   00 => *  ";
          remap_str.resize(20, ' ');
          remap_str = "[!" + j_str + "!" + i_str + " OR R " + "]:" + remap_str;
          std::cout << remap_str ;
          break;
        }
        else if( s.id == 1u )
        {
          remap_str += "   01 => *  ";
          remap_str.resize(20, ' ');
          remap_str = "[!" + j_str + " " + i_str + " OR R " + "]:" + remap_str;
          std::cout << remap_str ;
          break;
        }
        else if( s.id == 2u )
        {
          remap_str += "   10 => *  ";
          remap_str.resize(20, ' ');
          remap_str = "[ " + j_str + "!" + i_str + " OR R " + "]:" + remap_str;
          std::cout << remap_str ;
          break;
        }
        else if( s.id == 3u )
        {
          remap_str += "   11 => *  ";
          remap_str.resize(20, ' ');
          remap_str = "[ " + j_str + " " + i_str + " OR R " + "]:" + remap_str;
          std::cout << remap_str ;
          break;
        }
      case decomposition_types::AND2_:
        if( s.id == 0u )
        {
          remap_str += "   00 => *  ";
          remap_str.resize(20, ' ');
          remap_str = "[!" + j_str + "!" + i_str + " AND R " + "]:" + remap_str;
          std::cout << remap_str ;
          break;
        }
        else if( s.id == 1u )
        {
          remap_str += "   01 => *  ";
          remap_str.resize(20, ' ');
          remap_str = "[!" + j_str + " " + i_str + " AND R " + "]:" + remap_str;
          std::cout << remap_str ;
          break;
        }
        else if( s.id == 2u )
        {
          remap_str += "   10 => *  ";
          remap_str.resize(20, ' ');
          remap_str = "[ " + j_str + "!" + i_str + " AND R " + "]:" + remap_str;
          std::cout << remap_str ;
          break;
        }
        else if( s.id == 3u )
        {
          remap_str += "   11 => *  ";
          remap_str.resize(20, ' ');
          remap_str = "[ " + j_str + " " + i_str + " AND R " + "]:" + remap_str;
          std::cout << remap_str ;
          break;
        }
      default: 
        std::cout <<  "NONE";
        break;
    }
  }

  std::vector<decomposition_info_t> check_decomposition_type( TT tt, uint32_t i, uint32_t j )
  {
    assert( i < j );
    std::vector<decomposition_info_t> res;
    decomposition_info_t dec_info;
    dec_info.i = i;
    dec_info.j = j;

    uint32_t xi = support[i];
    uint32_t xj = support[j];

    const auto tt0 = cofactor0( tt, xj );
    const auto tt1 = cofactor1( tt, xj );

    const auto tt00 = cofactor0( tt0, xi );
    const auto tt01 = cofactor1( tt0, xi );
    const auto tt10 = cofactor0( tt1, xi );
    const auto tt11 = cofactor1( tt1, xi );

    const auto mk0 = cofactor0( mask, xj );
    const auto mk1 = cofactor1( mask, xj );

    const auto mk00 = cofactor0( mk0, xi );
    const auto mk01 = cofactor1( mk0, xi );
    const auto mk10 = cofactor0( mk1, xi );
    const auto mk11 = cofactor1( mk1, xi );

    const auto eq00_0 = is_const0( mask & tt00 );
    const auto eq01_0 = is_const0( mask & tt01 );
    const auto eq10_0 = is_const0( mask & tt10 );
    const auto eq11_0 = is_const0( mask & tt11 );
    const auto eq00_1 = is_const0( mask & ~tt00 );
    const auto eq01_1 = is_const0( mask & ~tt01 );
    const auto eq10_1 = is_const0( mask & ~tt10 );
    const auto eq11_1 = is_const0( mask & ~tt11 );

    const auto num_pairs =
      static_cast<uint32_t>( eq00_0 ) +
      static_cast<uint32_t>( eq01_0 ) +
      static_cast<uint32_t>( eq10_0 ) +
      static_cast<uint32_t>( eq11_0 ) +
      static_cast<uint32_t>( eq00_1 ) +
      static_cast<uint32_t>( eq01_1 ) +
      static_cast<uint32_t>( eq10_1 ) +
      static_cast<uint32_t>( eq11_1 );
    

    if ( num_pairs == 0  )
    {
      dec_info.type = decomposition_types::NONE;
      res.push_back( dec_info );
      return res;
    }
    if ( eq00_0 ) // F01 = F10 NES
      res.push_back( top2_remapping( 0u, i, j, decomposition_types::AND2_ ) );
    if ( eq01_0 ) // F00 = F11 ES
      res.push_back( top2_remapping( 1u, i, j, decomposition_types::AND2_ ) );
    if ( eq10_0 ) // F00 = F11 ES
      res.push_back( top2_remapping( 2u, i, j, decomposition_types::AND2_ ) );
    if ( eq11_0 ) // F00 = F11 ES
      res.push_back( top2_remapping( 3u, i, j, decomposition_types::AND2_ ) );
    if ( eq00_1 ) // F01 = F10 NES
      res.push_back( top2_remapping( 0u, i, j, decomposition_types::OR2_ ) );
    if ( eq01_1 ) // F00 = F11 ES
      res.push_back( top2_remapping( 1u, i, j, decomposition_types::OR2_ ) );
    if ( eq10_1 ) // F00 = F11 ES
      res.push_back( top2_remapping( 2u, i, j, decomposition_types::OR2_ ) );
    if ( eq11_1 ) // F00 = F11 ES
      res.push_back( top2_remapping( 3u, i, j, decomposition_types::OR2_ ) );

    return res;
  }

  std::vector<decomposition_info_t> check_decomposition_type( TT tt, uint32_t i )
  {
    std::vector<decomposition_info_t> res;
    decomposition_info_t dec_info;
    dec_info.i = i;
    dec_info.j = i;

    uint32_t xi = support[i];

    const auto tt0 = cofactor0( tt, xi );
    const auto tt1 = cofactor1( tt, xi );

    const auto eq0_0 = is_const0( mask & tt0 );
    const auto eq0_1 = is_const0( mask & ~tt0 );
    const auto eq1_0 = is_const0( mask & tt1 );
    const auto eq1_1 = is_const0( mask & ~tt1 );
    const auto eq10_ = equal( mask & tt1 , mask & ~tt0 );

    const auto num_pairs =
      static_cast<uint32_t>( eq0_0 ) +
      static_cast<uint32_t>( eq0_1 ) +
      static_cast<uint32_t>( eq1_0 ) +
      static_cast<uint32_t>( eq1_1 ) +
      static_cast<uint32_t>( eq10_ );
    

    if ( num_pairs == 0  )
    {
      dec_info.type = decomposition_types::NONE;
      res.push_back( dec_info );
      return res;
    }
    if ( eq0_0 ) // F01 = F10 NES
      res.push_back( top1_remapping( i, decomposition_types::AND_ ) );
    if ( eq0_1 ) // F00 = F11 ES
      res.push_back( top1_remapping( i, decomposition_types::LE_ ) );
    if ( eq1_0 ) // F00 = F11 ES
      res.push_back( top1_remapping( i, decomposition_types::LT_ ) );
    if ( eq1_1 ) // F00 = F11 ES
      res.push_back( top1_remapping( i, decomposition_types::OR_ ) );
    if ( eq10_ ) // F01 = F10 NES
      res.push_back( top1_remapping( i, decomposition_types::XOR_ ) );

    return res;
  }

  decomposition_info_t top2_remapping( uint32_t id, uint32_t i, uint32_t j, decomposition_types type )
  {
    assert( i < j );
    decomposition_info_t res;
    res.type = type;
    res.i = i;
    res.j = j; 
    res.id = id;
    
    uint32_t xi = support[i];
    uint32_t xj = support[j];

    TT A = cube_generator( id, xi, xj );

    if( type == decomposition_types::OR2_ )
      res.mask = ~A & mask ;
    else
      res.mask = A & mask ;

    res.func = remainder ;   

    return res;
  }

  decomposition_info_t top1_remapping( uint32_t i, decomposition_types type )
  {
    decomposition_info_t res;
    
    res.type = type;
    res.i = i;
    res.j = i; 
    
    uint32_t xi = support[i];

    TT A = X[xi];
    TT tt0 = cofactor0( remainder, xi );
    TT tt1 = cofactor1( remainder, xi );

    if( type == decomposition_types::AND_ )
    {
      res.func = tt1;
      res.mask = A & mask ;
    }
    else if( type == decomposition_types::OR_ )
    {
      res.func = tt0;
      res.mask = ~A & mask ;
    }
    else if( type == decomposition_types::LT_ )
    {
      res.func = tt0;
      res.mask = ~A & mask ;
    }
    else if( type == decomposition_types::LE_ )
    {
      res.func = tt1;
      res.mask = A & mask ;
    }

    return res;
  }

  signal<Ntk> dec_remap( decomposition_info_t decomposition )
  {
    uint32_t xi = support[decomposition.i];
    uint32_t xj = support[decomposition.j];

    switch( decomposition.type )
    {

      case decomposition_types::AND_:
        support.erase( support.begin() + decomposition.i );
        return ntk_.create_and( pis[xi], run() );
        break; 
      case decomposition_types::OR_:
        support.erase( support.begin() + decomposition.i );
        return ntk_.create_or( pis[xi], run() );
        break; 
      case decomposition_types::LT_:
        support.erase( support.begin() + decomposition.i );
        return ntk_.create_and( ntk_.create_not(pis[xi]), run() );
        break; 
      case decomposition_types::LE_:
        support.erase( support.begin() + decomposition.i );
        return ntk_.create_or( ntk_.create_not(pis[xi]), ntk_.create_not(run()) );
        break; 
      case decomposition_types::XOR_:
        support.erase( support.begin() + decomposition.i );
        return ntk_.create_xor( pis[xi], ntk_.create_not(run()) );
        break;
      case decomposition_types::OR2_:
        if( decomposition.id == 0u )
          return ntk_.create_or( ntk_.create_and(ntk_.create_not(pis[xi]), ntk_.create_not(pis[xj]) ), run()) ;
        else if( decomposition.id == 1u )
          return ntk_.create_or( ntk_.create_and(ntk_.create_not(pis[xj]), pis[xi]) , run() ) ;
        else if( decomposition.id == 2u )
          return ntk_.create_or( ntk_.create_and(pis[xj], ntk_.create_not(pis[xi]) ) , run() ) ;
        else if( decomposition.id == 3u )
          return ntk_.create_or( ntk_.create_and(pis[xj], pis[xi] ) , run() ) ;
        break; 
      case decomposition_types::AND2_:
        if( decomposition.id == 0u )
          return ntk_.create_and( ntk_.create_and(ntk_.create_not(pis[xi]), ntk_.create_not(pis[xj]) ), run()) ;
        else if( decomposition.id == 1u )
          return ntk_.create_and( ntk_.create_and(ntk_.create_not(pis[xj]), pis[xi]) , run() ) ;
        else if( decomposition.id == 2u )
          return ntk_.create_and( ntk_.create_and(pis[xj], ntk_.create_not(pis[xi]) ) , run() ) ;
        else if( decomposition.id == 3u )
          return ntk_.create_and( ntk_.create_and(pis[xj], pis[xi] ) , run() ) ;
        break; 

    }
  }
  
  #pragma endregion decomposition


  #pragma region symmetries

  enum class symmetry_types
  {
    NONE,
    ES_,
    NES_,
    MS_,
    SVS_,
    CSVS_
  };

  struct symmetry_info_t{
    symmetry_types type { symmetry_types::NONE };
    uint32_t i;
    uint32_t j;
    uint32_t id_ord{0}; 
    uint32_t id_sym{0}; 
    TT func; 
    TT mask;
    };

  void print_symmetry( symmetry_info_t s )
  {
    std::string remap_str = "";
    std::string i_str = std::to_string(s.i);
    std::string j_str = std::to_string(s.j);
    i_str.resize(2, ' ');
    j_str.resize(2, ' ');



    switch( s.type )
    {
      case symmetry_types::ES_:
        remap_str += ( s.id_ord == 0u ? "    00 -> 11 " : "    11 -> 00 " );
        remap_str.resize(20, ' ');
        remap_str = " ES [ " + i_str + ", " + j_str + "]:" + remap_str;
        std::cout << remap_str ;
        break;
      case symmetry_types::NES_: 
        remap_str += ( s.id_ord == 0u ? "    01 -> 10 " : "    10 -> 01 " );
        remap_str.resize(20, ' ');
        remap_str = "NES [ " + i_str + ", " + j_str + "]:" + remap_str;
        std::cout << remap_str ;
        break;
      case symmetry_types::MS_: 
        if( s.id_ord == 0u )
        {
          remap_str += "    00|01 -> 11|10 ";
          remap_str.resize(20, ' ');
          remap_str = " MS [ " + i_str + ", " + j_str + "]:" + remap_str;
          std::cout << remap_str ;
        }
        else if( s.id_ord == 1u )
        {
          remap_str += "    01|11 -> 10|00 ";
          remap_str.resize(20, ' ');
          remap_str = " MS [ " + i_str + ", " + j_str + "]:" + remap_str;
          std::cout << remap_str ;
        }
        else if( s.id_ord == 2u )
        {
          remap_str += "    10|00 -> 01|11 ";
          remap_str.resize(20, ' ');
          remap_str = " MS [ " + i_str + ", " + j_str + "]:" + remap_str;
          std::cout << remap_str ;
        }
        else if( s.id_ord == 3u )
        {
          remap_str += "    11|10 -> 00|01 ";
          remap_str.resize(20, ' ');
          remap_str = " MS [ " + i_str + ", " + j_str + "]:" + remap_str;
          std::cout << remap_str ;
        }        
        break;
      case symmetry_types::SVS_: 
        if( s.id_sym == 0u )
        {
          remap_str += ( s.id_ord == 0u ? "    00 -> 10 " : "    10 -> 00 " );
          remap_str.resize(20, ' ');
          remap_str = "[ SVS " + j_str + "]!" + i_str + " :" + remap_str;
          std::cout << remap_str ;
        }
        else if( s.id_sym == 1u )
        {
          remap_str += ( s.id_ord == 0u ? "    01 -> 11 " : "    11 -> 01 " );
          remap_str.resize(20, ' ');
          remap_str = "[ SVS " + j_str + "] " + i_str + " :" + remap_str;
          std::cout << remap_str ;
        }
        else if( s.id_sym == 2u )
        {
          remap_str += ( s.id_ord == 0u ? "    00 -> 01 " : "    01 -> 10 " );
          remap_str.resize(20, ' ');
          remap_str = "[ SVS " + i_str + "] " + j_str + " :" + remap_str;
          std::cout << remap_str ;
        }
        else if( s.id_sym == 3u )
        {
          remap_str += ( s.id_ord == 0u ? "    10 -> 11 " : "    11 -> 10 " );
          remap_str.resize(20, ' ');
          remap_str = "[ SVS " + i_str + "] " + j_str + " :" + remap_str;
          std::cout << remap_str ;
        }
        break;

      case symmetry_types::CSVS_: 
        if( s.id_sym == 0u )
        {
          remap_str += ( s.id_ord == 0u ? "    00|01 -> 10 " : "    10|00 -> 01 " );
          remap_str.resize(20, ' ');
          remap_str = "CSVS[!" + i_str + ",!" + j_str + "]:" + remap_str;
          std::cout << remap_str ;
        }
        else if( s.id_sym == 1u )
        {
          remap_str += ( s.id_ord == 0u ? "    00|01 -> 11 " : "    11|01 -> 00 " );
          remap_str.resize(20, ' ');
          remap_str = "CSVS[!" + i_str + ", " + j_str + "]:" + remap_str;
          std::cout << remap_str ;
        }
        else if( s.id_sym == 2u )
        {
          remap_str += ( s.id_ord == 0u ? "    00|10 -> 11 " : "    11|10 -> 00 " );
          remap_str.resize(20, ' ');
          remap_str = "CSVS[ " + i_str + ",!" + j_str + "]:" + remap_str;
          std::cout << remap_str ;
        }
        else if( s.id_sym == 3u )
        {
          remap_str += ( s.id_ord == 0u ? "    01|11 -> 10 " : "    11|10 -> 01 " );
          remap_str.resize(20, ' ');
          remap_str = "CSVS[ " + i_str + ", " + j_str + "]:" + remap_str;
          std::cout << remap_str ;
        }
        break;
      default: 
        std::cout <<  "NONE";
        break;
    }
  }
  
  symmetry_info_t simple_remapping( uint32_t from, uint32_t to, uint32_t i, uint32_t j, symmetry_types type, uint32_t id_symmetry = 0u )
  {
    assert( i < j );
    symmetry_info_t res;
    res.type = type;
    res.i = i;
    res.j = j; 
    res.id_ord = 1u*( from > to );
    res.id_sym = id_symmetry;
    
    uint32_t xi = support[i];
    uint32_t xj = support[j];

    TT A = cube_generator( from, xi, xj );
    TT B = cube_generator( to, xi, xj );
    TT ttA = cofactorG( remainder, from, xi, xj );
    TT ttB = cofactorG( remainder, to, xi, xj );
    TT mkA = cofactorG( mask, from, xi, xj );
    TT mkB = cofactorG( mask, to, xi, xj );

    res.mask = ( ~A & mask ) | ( B & mkA );

    TT TA = A & remainder; /* remapping for A */
    TT TB = B & ( ( mkB & remainder ) | ( mkA & ttA ) ); /* remapping for B */
    TT TR = ( ~A & ~B & remainder );

    res.func = TA | TB | TR ;
    
    return res;
  }

  symmetry_info_t multiform_remapping( uint32_t from1, uint32_t i, uint32_t j, symmetry_types type )
  {
    assert( i < j );
    symmetry_info_t res;
    res.type = type;
    res.i = i;
    res.j = j; 
    res.id_ord = from1;

    uint32_t xi = support[i];
    uint32_t xj = support[j];
    
    TT A = cube_generator( from1, xi, xj );
    uint32_t to1 = uint32_t( 3 - from1 );
    TT C = cube_generator( to1, xi, xj );

    uint32_t from2 = ( from1 == 0u )*1u + ( from1 == 1u )*3u + ( from1 == 3u )*2u;

    TT B = cube_generator( from2, xi, xj  );
    uint32_t to2 = uint32_t( 3 - from2 );
    TT D = cube_generator( to2, xi, xj  );

    TT ttA = cofactorG( remainder, from1, xi, xj );
    TT ttB = cofactorG( remainder, from2, xi, xj );
    TT ttC = cofactorG( remainder, to1, xi, xj );
    TT ttD = cofactorG( remainder, to2, xi, xj );

    TT mkA = cofactorG( mask, from1, xi, xj );
    TT mkB = cofactorG( mask, from2, xi, xj );
    TT mkC = cofactorG( mask, to1, xi, xj );
    TT mkD = cofactorG( mask, to2, xi, xj );

    res.mask = ( ~B & ~A & mask ) | ( ( C & mkA ) | ( D & mkB ) );
    
    TT preserved = (~A&~B&~C&~D)&remainder;

    TT modifiedA = A & remainder;
    TT modifiedB = B & remainder;

    TT modifiedC = C & ( ( mkA & ~mkC & ttA ) | ( mkC & remainder ) );
    TT modifiedD = D & ( ( mkB & ~mkD & ttB ) | ( mkD & remainder ) );

    res.func = preserved | modifiedA | modifiedC | modifiedB | modifiedD;
    
    return res;
  }

  symmetry_info_t compatible_remapping( uint32_t from1, uint32_t from2, uint32_t to, uint32_t i, uint32_t j, symmetry_types type, uint32_t id_symmetry = 0 )
  {
    assert( i < j );
    symmetry_info_t res;
    res.type = type;
    res.i = i;
    res.j = j; 
    res.id_ord = 1u*( from1 > from2 );
    res.id_sym = id_symmetry;

    uint32_t xi = support[i];
    uint32_t xj = support[j];

    uint32_t excluded = uint32_t( 6u - from1 - from2 - to );

    TT A = cube_generator( from1, xi, xj );
    TT B = cube_generator( from2, xi, xj );
    TT C = cube_generator( to, xi, xj );
    TT D = cube_generator( excluded, xi, xj );

    TT ttA = cofactorG( remainder, from1, xi, xj );
    TT ttB = cofactorG( remainder, from2, xi, xj );
    TT ttC = cofactorG( remainder, to, xi, xj );

    TT mkA = cofactorG( mask, from1, xi, xj );
    TT mkB = cofactorG( mask, from2, xi, xj );
    TT mkC = cofactorG( mask, to, xi, xj );

    res.mask = ( ~B & ~A & mask ) | ( C & ( mkA | mkB ) );
    
    TT TA = A & remainder;
    TT TB = B & remainder;
    TT TC = C & ( ( mkA & ttA ) | ( mkB & ttB ) | ( mkC & remainder ) );
    TT TR = ~A & ~B & ~C & remainder;

    res.func = TA | TB | TC | TR ;
    
    return res;
  }

  std::vector<symmetry_info_t> check_symmetry_type( TT tt, uint32_t i, uint32_t j )
  {
    assert( i < j );
    std::vector<symmetry_info_t> res;
    symmetry_info_t sym_info;
    sym_info.i = i;
    sym_info.j = j;

    uint32_t xi = support[i];
    uint32_t xj = support[j];

    const auto tt0 = cofactor0( tt, xj );
    const auto tt1 = cofactor1( tt, xj );

    const auto tt00 = cofactor0( tt0, xi );
    const auto tt01 = cofactor1( tt0, xi );
    const auto tt10 = cofactor0( tt1, xi );
    const auto tt11 = cofactor1( tt1, xi );

    const auto mk0 = cofactor0( mask, xj );
    const auto mk1 = cofactor1( mask, xj );

    const auto mk00 = cofactor0( mk0, xi );
    const auto mk01 = cofactor1( mk0, xi );
    const auto mk10 = cofactor0( mk1, xi );
    const auto mk11 = cofactor1( mk1, xi );

    const auto eq01 = equal( mk00&mk01&tt00, mk00&mk01&tt01 );
    const auto eq02 = equal( mk00&mk10&tt00, mk00&mk10&tt10 );
    const auto eq03 = equal( mk00&mk11&tt00, mk00&mk11&tt11 );
    const auto eq12 = equal( mk10&mk01&tt01, mk10&mk01&tt10 );
    const auto eq13 = equal( mk01&mk11&tt01, mk01&mk11&tt11 );
    const auto eq23 = equal( mk10&mk11&tt10, mk10&mk11&tt11 );

    const auto num_pairs =
      static_cast<uint32_t>( eq01 ) +
      static_cast<uint32_t>( eq02 ) +
      static_cast<uint32_t>( eq03 ) +
      static_cast<uint32_t>( eq12 ) +
      static_cast<uint32_t>( eq13 ) +
      static_cast<uint32_t>( eq23 );
    
    symmetry_info_t candid_ordate;

    if ( num_pairs == 0  )
    {
      sym_info.type = symmetry_types::NONE;
      res.push_back( sym_info );
      return res;
    }
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
    if ( eq02 ) // F00=F10
    {
      res.push_back( simple_remapping( 0, 2, i, j, symmetry_types::SVS_, 0u ) ); // SVS0X_
      res.push_back( simple_remapping( 2, 0, i, j, symmetry_types::SVS_, 0u ) );
    }
    if ( eq13 ) // F01=F11
    {
      res.push_back( simple_remapping( 1, 3, i, j, symmetry_types::SVS_, 1u ) ); // SVS1X_
      res.push_back( simple_remapping( 3, 1, i, j, symmetry_types::SVS_, 1u ) );
    }
    if ( eq01 ) // F01=F00
    {
      res.push_back( simple_remapping( 0, 1, i, j, symmetry_types::SVS_, 2u ) ); // SVSX0
      res.push_back( simple_remapping( 1, 0, i, j, symmetry_types::SVS_, 2u ) );    
    }
    if ( eq23 ) // F11=F10
    {
      res.push_back( simple_remapping( 2, 3, i, j, symmetry_types::SVS_, 3u ) ); // SVSX1_
      res.push_back( simple_remapping( 3, 2, i, j, symmetry_types::SVS_, 3u ) ); 
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
      res.push_back( compatible_remapping( 0, 1, 2, i, j, symmetry_types::CSVS_, 0u ) ); // CSVS00_
      res.push_back( compatible_remapping( 2, 0, 1, i, j, symmetry_types::CSVS_, 0u ) );
    }
    if( eq13 && eq01 && eq03 )
    {
      res.push_back( compatible_remapping( 0, 1, 3, i, j, symmetry_types::CSVS_, 1u ) ); // CSVS10_ old ij, new ji
      res.push_back( compatible_remapping( 3, 1, 0, i, j, symmetry_types::CSVS_, 1u ) );
    }
    if( eq02 && eq23 && eq03 )
    {
      res.push_back( compatible_remapping( 0, 2, 3, i, j, symmetry_types::CSVS_, 2u ) ); // CSVS01_ old ij, new ji
      res.push_back( compatible_remapping( 3, 2, 0, i, j, symmetry_types::CSVS_, 2u ) );
    }
    if( eq13 && eq23 && eq12 )
    {
      res.push_back( compatible_remapping( 1, 3, 2, i, j, symmetry_types::CSVS_, 3u ) ); // CSVS11_
      res.push_back( compatible_remapping( 3, 2, 1, i, j, symmetry_types::CSVS_, 3u ) );
    }
    return res;
  }

  void remap( symmetry_info_t symmetry )
  {
    uint32_t xi = support[symmetry.i];
    uint32_t xj = support[symmetry.j];

    switch( symmetry.type )
    {
      signal<Ntk> new_fi, new_fj;

      case symmetry_types::NES_:
        if( symmetry.id_ord == 0u )
        {
          new_fi = ntk_.create_and( pis[xi], pis[xj] );
          new_fj = ntk_.create_or( pis[xi], pis[xj] );
        }
        else if( symmetry.id_ord == 1u )
        {
          new_fi = ntk_.create_or( pis[xi], pis[xj] );
          new_fj = ntk_.create_and( pis[xi], pis[xj] );
        }
        else
          std::cerr << "id_ord not valid_ord" << std::endl;

        pis[xi] = new_fi;
        pis[xj] = new_fj;
        break; 

      case symmetry_types::ES_:
        if( symmetry.id_ord == 0u )
        {
          new_fj = ntk_.create_not( ntk_.create_and( pis[xi], ntk_.create_not( pis[xj] ) ) );
          new_fi = ntk_.create_not( ntk_.create_and( ntk_.create_not( pis[xi] ), pis[xj] ) );
        }
        else if( symmetry.id_ord == 1u )
        {
          new_fj = ntk_.create_and( ntk_.create_not( pis[xi] ), pis[xj] );
          new_fi = ntk_.create_and( pis[xi], ntk_.create_not( pis[xj] ) );
        }
        else
          std::cerr << "id_ord not valid_ord" << std::endl;

        pis[xi] = new_fi;
        pis[xj] = new_fj;
        break; 

      case symmetry_types::SVS_:
        if( symmetry.id_sym == 0u ) // SVS0X_ { SVS Xj } Xi'
        {
          new_fi = pis[xi];
          if( symmetry.id_ord == 0u )
            new_fj = ntk_.create_or( ntk_.create_not( pis[xi] ), pis[xj] );  // ( Xi & Xj' )'
          else
            new_fj = ntk_.create_and( pis[xi], pis[xj] ); // ( Xi & Xj )
        }
        else if( symmetry.id_sym == 1u ) // SVS1X  { SVS Xj } Xi
        {
          new_fi = pis[xi];
          if( symmetry.id_ord == 0u )
            new_fj = ntk_.create_or( pis[xi], pis[xj] ); 
          else
            new_fj = ntk_.create_and( ntk_.create_not( pis[xi] ), pis[xj] );
        }
        else if( symmetry.id_sym == 2u ) // SVSX0_  { SVS Xi } Xj'
        {
          if( symmetry.id_ord == 0u )
            new_fi = ntk_.create_or( pis[xi], ntk_.create_not( pis[xj] ) );
          else
            new_fi = ntk_.create_and( pis[xi], pis[xj] );
          new_fj = pis[xj];
        }
        else if( symmetry.id_sym == 3u ) //SVSX1_  { SVS Xi } Xj
        {
          if( symmetry.id_ord == 0u )
            new_fi = ntk_.create_or( pis[xi], pis[xj] );
          else
            new_fi = ntk_.create_and( pis[xi], ntk_.create_not( pis[xj] ) );
          new_fj = pis[xj];
        }
        else
          std::cerr << "wrong symmetry identifier for SVS" ;
        
        pis[xi] = new_fi;
        pis[xj] = new_fj;
        
        break;

      case symmetry_types::MS_:
        if( symmetry.id_ord == 0u )
        {
          pis[xi] = ntk_.create_not( ntk_.create_xor( pis[xi], pis[xj] ) );
          support.erase( support.begin() + symmetry.j );
        }
        else if( symmetry.id_ord == 1u )
        {
          pis[xj] = ntk_.create_xor( pis[xi], pis[xj] );
          support.erase( support.begin() + symmetry.i );
        }        
        else if( symmetry.id_ord == 2u )
        {
          pis[xj] = ntk_.create_not( ntk_.create_xor( pis[xi], pis[xj] ) );
          support.erase( support.begin() + symmetry.i );
        }        
        else if( symmetry.id_ord == 3u )
        {
          pis[xi] = ntk_.create_xor( pis[xi], pis[xj] );
          support.erase( support.begin() + symmetry.j );
        }  
        break; 

      case symmetry_types::CSVS_:
      {
        if( symmetry.id_sym == 0u )
        {
          if( symmetry.id_ord == 0u )
          {
            pis[xi] = ntk_.create_and( pis[xi], pis[xj] );
            support.erase( support.begin() + symmetry.j );
          }
          else
          {
            pis[xj] = ntk_.create_and( pis[xi], pis[xj] );
            support.erase( support.begin() + symmetry.i );
          }
        }
        else if( symmetry.id_sym == 1u )
        {
          if( symmetry.id_ord == 0u )
          {
            pis[xi] = ntk_.create_or( pis[xi], ntk_.create_not( pis[xj] ) );
            support.erase( support.begin() + symmetry.j );
          }
          else
          {
            pis[xj] = ntk_.create_and( ntk_.create_not( pis[xi] ), pis[xj] );
            support.erase( support.begin() + symmetry.i );
          }
        }
        else if( symmetry.id_sym == 2u )
        {
          if( symmetry.id_ord == 0u )
          {
            pis[xj] = ntk_.create_or( ntk_.create_not( pis[xi] ), pis[xj] );
            support.erase( support.begin() + symmetry.i );
          }
          else
          {
            pis[xi] = ntk_.create_and( pis[xi], ntk_.create_not( pis[xj] ) );
            support.erase( support.begin() + symmetry.j );
          }
        }
        else if( symmetry.id_sym == 3u )
        {
          if( symmetry.id_ord == 0u )
          {
            pis[xj] =  ntk_.create_or( pis[xi], pis[xj] );
            support.erase( support.begin() + symmetry.i );
          }
          else
          {
            pis[xi] =  ntk_.create_or( pis[xi], pis[xj] );
            support.erase( support.begin() + symmetry.j );
          }
        }
        else
          std::cerr << "wrong symmetry identifier for CSVS" ;
      }
      break; 
    }
  }

  #pragma end region symmetries

  /*
   * Compute the general cofactor with respect to the cube G
   *            ji 
   *   G = 0 -> 00 : F( Xi=0, Xj=0 )
   *   G = 1 -> 01 : F( Xi=0, Xj=1 )
   *   G = 2 -> 10 : F( Xi=1, Xj=0 )
   *   G = 3 -> 11 : F( Xi=1, Xj=1 )
   * NB: notice the swap of (i,j), assuming that i < j in the given label
  */
  TT cofactorG( TT fn, uint32_t G, uint32_t i, uint32_t j )
  {
    switch ( G )
    {
      case 0u: return cofactor0(cofactor0( fn, j), i ); /* F00 */
      case 1u: return cofactor1(cofactor0( fn, j), i ); /* F01 */
      case 2u: return cofactor0(cofactor1( fn, j), i ); /* F10 */
      case 3u: return cofactor1(cofactor1( fn, j), i ); /* F11 */
    }
  }

  /*
   * Compute the truth table associated to the cube cube for variables Xi Xj
   *               ji 
   *   cube = 0 -> 00 : Xi' & Xj'
   *   cube = 1 -> 01 : Xi  & Xj'
   *   cube = 2 -> 10 : Xi' & Xj
   *   cube = 3 -> 11 : Xi  & Xj
   * NB: notice the swap of (i,j), assuming that i < j in the given label
  */
  TT cube_generator( uint32_t cube, uint32_t i, uint32_t j )
  {
    switch ( cube )
    {
    case 0u: 
      return ~X[j] & ~X[i];
      break;
    case 1u: 
      return ~X[j] & X[i];
      break;
    case 2u: 
      return X[j] & ~X[i];
      break;
    case 3u: 
      return X[j] & X[i];
      break;
    }
  }


  #pragma region cost evaluation
  enum class gate_names
  {
    OR_,
    CRO_,
    INV_,
    SPL_,
    BUF_,
    AND_,
    XOR_
  };

  struct gate_with_cost{
    
    gate_with_cost( )
    {}

    gate_with_cost( gate_names type, uint32_t area, uint32_t delay ) :
    name( name ), area( area ), delay( delay )
    {}

    gate_names name;
    uint32_t area;
    uint32_t delay;
  };
  
  void initialize_gate_library()
  {
    inv_ = gate_with_cost( gate_names::INV_, 1u, 1u );
    buf_ = gate_with_cost( gate_names::BUF_, 1u, 1u );
    cro_ = gate_with_cost( gate_names::CRO_, 1u, 1u );
    spl_ = gate_with_cost( gate_names::SPL_, 1u, 1u );
    and_ = gate_with_cost( gate_names::AND_, 1u, 1u );
    xor_ = gate_with_cost( gate_names::XOR_, 1u, 1u );
    or_ = gate_with_cost( gate_names::OR_, 1u, 1u );
  }

  uint32_t cost_remapping_cell( symmetry_info_t symmetry )
  {
    uint32_t d1, d2, d3, d4, dM, nbuf;
    switch ( symmetry.type )
    {
      case symmetry_types::NES_:
        d1 = spl_.delay + and_.delay;
        d2 = spl_.delay + cro_.delay + or_.delay;
        d3 = spl_.delay + cro_.delay + and_.delay;
        d4 = spl_.delay + or_.delay;
        dM = std::max( { d1, d2, d3, d4 } );
        nbuf = 4*dM - d1 - d2 - d3 - d4;
        return and_.area + or_.area + cro_.area + 2*spl_.area + nbuf * buf_.area;
      case symmetry_types::ES_:
        if( symmetry.id_ord == 0u )
        {
          d1 = spl_.delay + or_.delay;
          d2 = spl_.delay + cro_.delay + or_.delay + inv_.delay;
          d3 = spl_.delay + cro_.delay + or_.delay + inv_.delay;
          d4 = spl_.delay + or_.delay;
          dM = std::max( { d1, d2, d3, d4 } );
          nbuf = 4*dM - d1 - d2 - d3 - d4;  
          return 2*or_.area + cro_.area + 2*spl_.area + 2*inv_.area + nbuf * buf_.area;
        }
        else
        {
          d1 = spl_.delay + and_.delay;
          d2 = spl_.delay + cro_.delay + or_.delay;
          d3 = spl_.delay + cro_.delay + and_.delay;
          d4 = spl_.delay + or_.delay;
          dM = std::max( { d1, d2, d3, d4 } );
          nbuf = 4*dM - d1 - d2 - d3 - d4;
          return 2*and_.area + cro_.area + 2*spl_.area + 2*inv_.area + nbuf * buf_.area;
        }
      case symmetry_types::SVS_:
        if( symmetry.id_sym == 0u )
        {
          if( symmetry.id_ord == 0u )
          {
            d1 = spl_.delay;
            d2 = spl_.delay + inv_.delay + or_.delay;
            d3 = or_.delay;
            dM = std::max( { d1, d2, d3 } );
            nbuf = 3*dM - d1 - d2 - d3 ;
            return or_.area + spl_.area + inv_.area + nbuf * buf_.area;            
          }
          else
          {
            d1 = spl_.delay;
            d2 = spl_.delay + and_.delay;
            d3 = and_.delay;
            dM = std::max( { d1, d2, d3 } );
            nbuf = 3*dM - d1 - d2 - d3 ;
            return or_.area + spl_.area + inv_.area + nbuf * buf_.area;           
          }
        }
        else if( symmetry.id_sym == 1u )
        {
          if( symmetry.id_ord == 0u )
          {
            d1 = spl_.delay;
            d2 = spl_.delay + or_.delay;
            d3 = or_.delay;
            dM = std::max( { d1, d2, d3 } );
            nbuf = 3*dM - d1 - d2 - d3 ;
            return or_.area + spl_.area + nbuf * buf_.area;            
          }
          else
          {
            d1 = spl_.delay;
            d2 = spl_.delay + and_.delay + inv_.delay;
            d3 = and_.delay;
            dM = std::max( { d1, d2, d3 } );
            nbuf = 3*dM - d1 - d2 - d3 ;
            return and_.area + spl_.area + inv_.area + nbuf * buf_.area;           
          }
        }
        else if( symmetry.id_sym == 2u )
        {
          if( symmetry.id_ord == 0u )
          {
            d1 = or_.delay;
            d2 = spl_.delay + or_.delay + inv_.delay;
            d3 = spl_.delay;
            dM = std::max( { d1, d2, d3 } );
            nbuf = 3*dM - d1 - d2 - d3 ;
            return or_.area + spl_.area + inv_.area + nbuf * buf_.area;            
          }
          else
          {
            d1 = and_.delay;
            d2 = spl_.delay + and_.delay;
            d3 = spl_.delay;
            dM = std::max( { d1, d2, d3 } );
            nbuf = 3*dM - d1 - d2 - d3 ;
            return and_.area + spl_.area + nbuf * buf_.area;           
          }
        }
        else
        {
          if( symmetry.id_ord == 0u )
          {
            d1 = or_.delay;
            d2 = spl_.delay + or_.delay;
            d3 = spl_.delay;
            dM = std::max( { d1, d2, d3 } );
            nbuf = 3*dM - d1 - d2 - d3 ;
            return or_.area + spl_.area + nbuf * buf_.area;            
          }
          else
          {
            d1 = and_.delay;
            d2 = spl_.delay + and_.delay + inv_.delay;
            d3 = spl_.delay;
            dM = std::max( { d1, d2, d3 } );
            nbuf = 3*dM - d1 - d2 - d3 ;
            return and_.area + spl_.area + inv_.area + nbuf * buf_.area;           
          }
        }
      case symmetry_types::MS_:
        if( symmetry.id_ord == 0u )
          return xor_.area + inv_.area ;           
        else if( symmetry.id_ord == 1u )
          return xor_.area ;
        else if( symmetry.id_ord == 2u )
          return xor_.area + inv_.area ; 
        else
          return xor_.area ; 

      case symmetry_types::CSVS_:
        if( symmetry.id_sym == 0u )
          return and_.area;
        else if( ( symmetry.id_sym == 1u ) || ( symmetry.id_sym == 2u  ) )
        {
          nbuf = inv_.delay;
          if( symmetry.id_ord == 0u )
            return or_.area + inv_.area + nbuf*buf_.area;
          else 
            return and_.area + inv_.area + nbuf*buf_.area;
        }   
        else
          return or_.area ; 

    }
  }
  #pragma end region cost evaluation

  #pragma region erase redundants
  void erase_redundant()
  {    
    uint32_t init = support.size()-1;
    for( uint32_t i = init; i > 0; i-- )
    {
      uint32_t xi = support[i];
      TT mk0 = cofactor0( mask, xi );
      TT mk1 = cofactor1( mask, xi );
      TT tt0 = cofactor0( remainder, xi );
      TT tt1 = cofactor1( remainder, xi );

      if( equal( mk0&mk1&tt1, mk0&mk1&tt0 ) )
      {
        std::cout << "erase " << xi << std::endl;
        support.erase( support.begin() + i );
      }
    }
  }
  #pragma endregion erase redundant

public:
  //signal<Ntk> run()
  signal<Ntk> run()
  {
    kitty::print_binary( remainder );
    std::cout << std::endl;
    bool game_over {false};
    symmetry_types last_type;

    auto km_tt = kitty::karnaugh_map( remainder );
    km_tt.print(mask);
    uint32_t COST = 0;
    
    while( !game_over && ( support.size() > 1 ) )
    {

      if( is_const0( remainder&mask ) )
        return ntk_.get_constant( false );
      else if( is_const0( ~remainder&mask ) )
        return ntk_.get_constant( true );

      for( uint32_t i{0u}; i < support.size(); ++i )
      {
        if( equal( remainder&mask, X[support[i]]&mask ) )
          return pis[support[i]];
        if( equal( remainder&mask, ~X[support[i]]&mask ) )
          return ntk_.create_not( pis[support[i]] );
      }

      /* detect symmetries */
      std::vector<symmetry_info_t> symmetries;
      std::vector<decomposition_info_t> decompositions1;
      uint32_t k = 0;
      uint32_t choice;
      std::string order;

      std::vector<decomposition_info_t> dec1_v;

      for ( auto j = 0u; j < support.size(); ++j )
      {
        dec1_v = check_decomposition_type( remainder, j );
        if( !(dec1_v.size()==1u && dec1_v[0].type== decomposition_types::NONE ))
        {
          for( auto s : dec1_v )
          {
            std::string ref = std::to_string( k++ );
            ref.resize(3, ' ');

            std::string dc_string = " |DC|= " + std::to_string( kitty::count_ones( ~s.mask ) );
            dc_string.resize(11, ' ');

            std::string cost = " |G|= " + std::to_string( COST + 0 );//cost_remapping_cell( s ) );

            std::cout << ref << " "; 
            print_decomposition(s) ; 
            std::cout << dc_string << cost <<  std::endl;

            decompositions1.push_back( s );
          }
        }
      }

      for ( auto j = 1u; j < support.size(); ++j )
      {
        for ( auto i = 0u; i < j; ++i )
        {
	        std::vector<decomposition_info_t> dec2_v = check_decomposition_type( remainder, i, j );

	        if( !(dec2_v.size()==1u && dec2_v[0].type== decomposition_types::NONE ))
	        {
	          for( auto s : dec2_v )
	          {
              std::string ref = std::to_string( k++ );
              ref.resize(3, ' ');

              std::string dc_string = " |DC|= " + std::to_string( kitty::count_ones( ~s.mask ) );
              dc_string.resize(11, ' ');

              std::string cost = " |G|= " + std::to_string( COST + 2 );//cost_remapping_cell( s ) );

              std::cout << ref << " "; 
              print_decomposition(s) ; 
              std::cout << dc_string << cost <<  std::endl;

              decompositions1.push_back( s );
	          }
	        }
        }
      }
      std::cout << "remappings" << std::endl;
      k=0;
      for ( auto j = 1u; j < support.size(); ++j )
      {
        for ( auto i = 0u; i < j; ++i )
        {
	        std::vector<symmetry_info_t> sym_v = check_symmetry_type( remainder, i, j );

	        if( !(sym_v.size()==1u && sym_v[0].type== symmetry_types::NONE ))
	        {
	          for( auto s : sym_v )
	          {
              std::string ref = std::to_string( k++ );
              ref.resize(3, ' ');

              std::string dc_string = " |DC|= " + std::to_string( kitty::count_ones( ~s.mask ) );
              dc_string.resize(11, ' ');

              std::string cost = " |G|= " + std::to_string( COST + cost_remapping_cell( s ) );

              std::cout << ref << " "; 
              print_symmetry(s) ; 
              std::cout << dc_string << cost <<  std::endl;

              symmetries.push_back( s );
	          }
	        }
        }
      }
      game_over = ( ( symmetries.size() < 1u ) && ( decompositions1.size() < 1u ) ); 

      if( !game_over )
      {
        bool decompose{false};
        if( decompositions1.size() > 0 )
        {
          std::string schoice;
          std::cout << "Decompose? [y/n] " << std::endl;
          std::cin >> schoice; 
          decompose = (schoice == "y");       
        }
        if( decompose )
        {
          std::cout << "Choose the decompo to exploit: " << std::endl;
          std::cin >> choice;
          std::cout << "Remapping " << " "; print_decomposition( decompositions1[choice] ); 
          std::cout << std::endl;
          COST += 1;//cost_remapping_cell( symmetries[choice] );   

          remainder = decompositions1[choice].func;
          kitty::print_binary( remainder );
          std::cout << std::endl;
          mask = decompositions1[choice].mask;

          return dec_remap( decompositions1[choice] );

        }
        else
        {
          std::cout << "Choose the remapping to exploit: " << std::endl;
          std::cin >> choice;
          std::cout << "Remapping " << " "; print_symmetry( symmetries[choice] ); 
          std::cout << std::endl;
          COST += cost_remapping_cell( symmetries[choice] );

          remainder = symmetries[choice].func;
          kitty::print_binary( remainder );
          std::cout << std::endl;
          mask = symmetries[choice].mask;
          remap( symmetries[choice] );
          erase_redundant();
          for( auto x : support )
            std::cout << pis[x].index << " ";
          std::cout << std::endl;
          last_type = symmetries[choice].type;
          auto km_tt = kitty::karnaugh_map( remainder );
          km_tt.print(mask);
        }

      }
      else
      {
        std::cout << "GAME OVER" << std::endl;
      }
    }
    
    return is_const0( remainder & X[support[0]] & mask ) ? !pis[support[0]] : pis[support[0]];
  }

private:
  Ntk& ntk_;
  TT remainder;
  TT target;
  TT mask;
  std::vector<uint32_t> support;
  std::vector<signal<Ntk>> pis;
  std::vector<TT> X;
  xminsyn_h_params const& _ps;

  gate_with_cost inv_;
  gate_with_cost buf_;
  gate_with_cost spl_;
  gate_with_cost cro_;
  gate_with_cost xor_;
  gate_with_cost and_;
  gate_with_cost or_;
};

} // namespace detail

/*! \brief crossing minimization synthesis engine
 *
 * This function heuristically minimizes the number of crossings by using a 
 * symmetry-based synthesis heuristic.
 */
template<class Ntk, class TT>
signal<Ntk> xminsyn_h( Ntk& ntk, TT const& func, std::vector<signal<Ntk>> const& children, xminsyn_h_params const& ps = {} )
{
  detail::xminsyn_h_impl impl( ntk, func, children, ps );
  return impl.run();
}

} // namespace mockturtle
