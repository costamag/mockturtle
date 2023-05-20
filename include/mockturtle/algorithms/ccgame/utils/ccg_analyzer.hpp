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
  \file ccg_analyzer.hpp
  \brief cut analyzer for the ccgame

  \author AnxRea Costamagna
*/
#pragma once

#include "ccg_node.hpp"
#include "ccg_cut.hpp"
#include <kitty/partial_truth_table.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/bit_operations.hpp>
#include <kitty/operators.hpp>
#include <stdio.h>
#include <stack>

namespace mockturtle
{

namespace ccgame
{

using TT = kitty::partial_truth_table;
using DTT = kitty::dynamic_truth_table;

DTT cube_generator( uint32_t cube, DTT Xr, DTT Xl )
{
  switch ( cube )
  {
    case 0: return ~Xl & ~Xr; break;
    case 1: return ~Xl &  Xr; break;
    case 2: return  Xl & ~Xr; break;
    case 3: return  Xl &  Xr; break;
  }
}

DTT cofactorG( DTT fn, uint32_t G, uint32_t idL, uint32_t idR )
{
  switch ( G )
  {
    case 0u: return cofactor0(cofactor0( fn, idL), idR ); /* F00 */
    case 1u: return cofactor1(cofactor0( fn, idL), idR ); /* F01 */
    case 2u: return cofactor0(cofactor1( fn, idL), idR ); /* F10 */
    case 3u: return cofactor1(cofactor1( fn, idL), idR ); /* F11 */
  }
}

/*! \brief Gate type in the ccgame namespace. Convention Xl=1100, Xr=1010
 */
struct symmetry_t
{
  uint8_t type;
  int idL;
  int idR;
  DTT tt;
  DTT mk;
  int rwd;

  symmetry_t(){}
  symmetry_t( uint8_t type, int idL, int idR ): type(type), idL(idL), idR(idR){}

  void remapping_equations( std::vector<DTT> X, DTT func, DTT mask )
  {
    uint32_t idA = ((type & 0xC0) >> 6u ) & 3u;
    uint32_t idC = ((type & 0x30) >> 4u ) & 3u;
    uint32_t idB = ((type & 0x0C) >> 2u ) & 3u;
    uint32_t idD = type & 3u;

    DTT A = cube_generator( idA, X[idR], X[idL] );
    DTT B = cube_generator( idB, X[idR], X[idL] );
    DTT C = cube_generator( idC, X[idR], X[idL] );
    DTT D = cube_generator( idD, X[idR], X[idL] );

    DTT ttA = cofactorG( func, idA, idL, idR );
    DTT ttB = cofactorG( func, idB, idL, idR );
    DTT ttC = cofactorG( func, idC, idL, idR );
    DTT ttD = cofactorG( func, idD, idL, idR );

    DTT mkA = cofactorG( mask, idA, idL, idR );
    DTT mkB = cofactorG( mask, idB, idL, idR );
    DTT mkC = cofactorG( mask, idC, idL, idR );
    DTT mkD = cofactorG( mask, idD, idL, idR );

    if( idA == idB && idC == idD ) //simple remapping
    {
      mk = ( mask & ~A ) | ( C & mkA );
      rwd = kitty::count_zeros( mk );
      DTT TA = A & func; /* remapping for A */
      DTT TC = C & ( ( mkC & func ) | ( mkA & ttA ) ); /* remapping for C */
      DTT TR = ( ~A & ~C & func );
      tt =  TA | TC | TR;
    }
    else if( idA == idB ) // compatible remapping
    {
      mk = ( ~B & ~A & mask ) | ( C & ( mkA | mkB ) ) ;
      rwd = kitty::count_zeros( mk ); //@ltry
      DTT TA = A & func;
      DTT TB = B & func;
      DTT TC = C & ( ( mkA & ttA ) | ( mkB & ttB ) | ( mkC & func ) );
      DTT TR = ~A & ~B & ~C & func;
      tt = TA | TB | TC | TR;
    }
    else // multiform remapping
    {
      mk = ( ~B & ~A & mask ) | ( ( C & mkA ) | ( D & mkB ) );
      rwd = kitty::count_zeros( mk ); // @ltry
      
      DTT preserved = (~A&~B&~C&~D)&func;
      DTT modifiedA = A & func;
      DTT modifiedB = B & func;
      DTT modifiedC = C & ( ( mkA & ~mkC & ttA ) | ( mkC & func ) );
      DTT modifiedD = D & ( ( mkB & ~mkD & ttB ) | ( mkD & func ) );
      tt = preserved | modifiedA | modifiedB | modifiedC | modifiedD;
    }
  }
};



class analyzer_t
{

public:

  analyzer_t();
  ~analyzer_t();
  /* exhaustive analyzers */
  cut_t enumerate_divs( cut_t );
  std::vector<symmetry_t> find_symmetries( std::vector<DTT>, DTT, DTT, std::vector<int> );
  /* symmetry-based */
  void print_symmetries( std::vector<symmetry_t> );
};

#pragma region constructors
analyzer_t::analyzer_t(){}
analyzer_t::~analyzer_t(){}
#pragma enxRegion

#pragma region exhaustive analysis
/*! \brief combines the gates in the last cut to propose all possible nodes*/
cut_t analyzer_t::enumerate_divs( cut_t cut )
{
    cut_t divs;

    for( int i{0}; i < cut.nodes.size() ; ++i )
    {
        node_t d = cut.nodes[i];
        divs.add_node( d.tt, PRJL, d.id, d.id );
    }

    node_t xL;
    node_t xR;

    for( int iR{0}; iR < cut.nodes.size()-1 ; ++iR )
    {
        for( int iL{iR+1}; iL < cut.nodes.size(); ++iL ) // iL > iR
        {
            xL = cut.nodes[iL];
            xR = cut.nodes[iR];
            divs.add_node(  xL.tt &  xR.tt, AI11, xL.id, xR.id );
            divs.add_node(  xL.tt & ~xR.tt, AI10, xL.id, xR.id );
            divs.add_node( ~xL.tt &  xR.tt, AI01, xL.id, xR.id );
            divs.add_node( ~xL.tt & ~xR.tt, AI00, xL.id, xR.id );
            divs.add_node(  xL.tt ^  xR.tt, EXOR, xL.id, xR.id );
        }
    }
    return divs;
} 
#pragma enxRegion exhaustive analysis

#pragma region symmetry analysis

/*! \brief symmetry analysis of variable iR with all the ones on the right 
std::vector<cut_t> analyzer_t::symmetry_analysis( cut_t cut, int idEnd )
{
  std::vector<cut_t> res;
  /* two variables symmetry check 
  for( int iL{0}; iL < cut.size(); ++iL )
  {
  printf("iR = %d\n", iL);

    node_t xL = cut.nodes[iL];
    printf( "L el %d anc %d is remapped ? %d\n", iL, xL.remapped_pi(), xL.is_remapped() );
    if( !xL.is_remapped() ) continue;

    DTT ttR0  = kitty::cofactor0( cut.tt, xL.remapped_pi() );
    DTT ttR1  = kitty::cofactor1( cut.tt, xL.remapped_pi() );
    DTT mkR0  = kitty::cofactor0( cut.mk, xL.remapped_pi() );
    DTT mkR1  = kitty::cofactor1( cut.mk, xL.remapped_pi() );
    /* todo: add top-decomposition check 

    /* symmetry check 
    for( int iR{iL+1}; iR < cut.size(); ++iR )
    {
      node_t xR = cut.nodes[iR];
      printf( "R el %d anc %d is remapped ? %d\n", iR, xR.remapped_pi(), xR.is_remapped() );
      if( !xL.is_remapped() ) continue;
      const auto tt00 = kitty::cofactor0( ttR0, xL.remapped_pi() );
      const auto tt01 = kitty::cofactor0( ttR1, xL.remapped_pi() );
      const auto tt10 = kitty::cofactor1( ttR0, xL.remapped_pi() );
      const auto tt11 = kitty::cofactor1( ttR1, xL.remapped_pi() );
      const auto mk00 = kitty::cofactor0( mkR0, xL.remapped_pi() );
      const auto mk01 = kitty::cofactor0( mkR1, xL.remapped_pi() );
      const auto mk10 = kitty::cofactor1( mkR0, xL.remapped_pi() );
      const auto mk11 = kitty::cofactor1( mkR1, xL.remapped_pi() );

      const auto eq01 = kitty::equal( mk00&mk01&tt00, mk00&mk01&tt01 );
      const auto eq02 = kitty::equal( mk00&mk10&tt00, mk00&mk10&tt10 );
      const auto eq03 = kitty::equal( mk00&mk11&tt00, mk00&mk11&tt11 );
      const auto eq12 = kitty::equal( mk10&mk01&tt01, mk10&mk01&tt10 );
      const auto eq13 = kitty::equal( mk01&mk11&tt01, mk01&mk11&tt11 );
      const auto eq23 = kitty::equal( mk10&mk11&tt10, mk10&mk11&tt11 );

      const auto num_pairs =  static_cast<uint32_t>( eq01 ) + static_cast<uint32_t>( eq02 ) + static_cast<uint32_t>( eq03 ) + 
                              static_cast<uint32_t>( eq12 ) + static_cast<uint32_t>( eq13 ) + static_cast<uint32_t>( eq23 );

      if( num_pairs > 0 )   
      {     
        if ( eq12 ) // F01 = F10 NES
        {
          res.push_back( remap( cut, xL, xR, 0x66 ) ); // 01 -> 10
          res.push_back( remap( cut, xL, xR, 0x99 ) ); // 10 -> 01
        }
        if ( eq03 ) // F00 = F11 ES
        {
          res.push_back( remap( cut, xL, xR, 0x33 ) ) ; // 00 -> 11
          res.push_back( remap( cut, xL, xR, 0xCC ) ) ; // 11 -> 00
        }
        if ( eq01 ) // F01=F00
        {
          res.push_back( remap( cut, xL, xR, 0x11 ) ); // 1:00->01
          res.push_back( remap( cut, xL, xR, 0x44 ) ); // 4:01->00
        }
        if ( eq02 ) // F00=F10
        {
          res.push_back( remap( cut, xL, xR, 0x22 ) ); // 2:00->10
          res.push_back( remap( cut, xL, xR, 0x88 ) ); // 8:10->00
        }
        if ( eq13 ) // F01=F11
        {
          res.push_back( remap( cut, xL, xR, 0x77 ) ); // 7:01->11
          res.push_back( remap( cut, xL, xR, 0xDD ) ); // D:11->01
        }
        if ( eq23 ) // F11=F10
        {
          res.push_back( remap( cut, xL, xR, 0xBB ) ); // B:10->11
          res.push_back( remap( cut, xL, xR, 0xEE ) ); // E:11->10
        }
        if ( eq12 && eq03 ) // F01=F10 and F00=F11
        {
          res.push_back( remap( cut, xL, xR, 0x36 ) ); // 3:00->11 6:01->10 
          res.push_back( remap( cut, xL, xR, 0x6C ) ); // 6:01->10 C:11->00      
          res.push_back( remap( cut, xL, xR, 0x9C ) ); // 9:10->01 C:11->00     
          res.push_back( remap( cut, xL, xR, 0x39 ) ); // 3:00->11 9:10->01
        }
        if( eq02 && eq01 && eq12 )
        {
          res.push_back( remap( cut, xL, xR, 0x19 ) ); // 1:00->01 9:10->01
          res.push_back( remap( cut, xL, xR, 0x26 ) ); // 2:00->10 6:01->10
        }
        if( eq13 && eq01 && eq03 )
        {
          res.push_back( remap( cut, xL, xR, 0x37 ) ); // 3:00->11 7:01->11
          res.push_back( remap( cut, xL, xR, 0x4C ) ); // 4:01->00 C:11->00
        }
        if( eq02 && eq23 && eq03 )
        {
          res.push_back( remap( cut, xL, xR, 0x8C ) ); // 8:10->00 C:11->00
          res.push_back( remap( cut, xL, xR, 0x3B ) ); // 3:00->11 B:10->11
        }
        if( eq13 && eq23 && eq12 )
        {
          res.push_back( remap( cut, xL, xR, 0x6E ) ); // 6:01->10 E:11->10
          res.push_back( remap( cut, xL, xR, 0x9D ) ); // 9:10->01 D:11->01
        }
      }
    }
  }
  return res;
}

/*! \brief symmetry analysis of variable iR with all the ones on the right */
std::vector<symmetry_t> analyzer_t::find_symmetries( std::vector<DTT> xs, DTT tt, DTT mk, std::vector<int> ids )
{
  std::vector<symmetry_t> res;
  int idL, idR;
  for( int iR{0}; iR < ids.size(); ++iR )
  {
    idR = ids[iR];
    if( idR < 0 ) continue; 
    DTT tt0  = kitty::cofactor0( tt, idR );
    DTT tt1  = kitty::cofactor1( tt, idR );
    DTT mk0  = kitty::cofactor0( mk, idR );
    DTT mk1  = kitty::cofactor1( mk, idR );
    /* todo: add top-decomposition check */

    /* symmetry check */
    for( int iL{iR+1}; iL < ids.size(); ++iL )
    {
      idL = ids[iL];
      if( idL < 0 ) continue; 

      assert( idL > idR );

      const auto tt00 = kitty::cofactor0( tt0, idL );
      const auto tt01 = kitty::cofactor0( tt1, idL );
      const auto tt10 = kitty::cofactor1( tt0, idL );
      const auto tt11 = kitty::cofactor1( tt1, idL );
      const auto mk00 = kitty::cofactor0( mk0, idL );
      const auto mk01 = kitty::cofactor0( mk1, idL );
      const auto mk10 = kitty::cofactor1( mk0, idL );
      const auto mk11 = kitty::cofactor1( mk1, idL );

      const auto eq01 = kitty::equal( mk00&mk01&tt00, mk00&mk01&tt01 );
      const auto eq02 = kitty::equal( mk00&mk10&tt00, mk00&mk10&tt10 );
      const auto eq03 = kitty::equal( mk00&mk11&tt00, mk00&mk11&tt11 );
      const auto eq12 = kitty::equal( mk10&mk01&tt01, mk10&mk01&tt10 );
      const auto eq13 = kitty::equal( mk01&mk11&tt01, mk01&mk11&tt11 );
      const auto eq23 = kitty::equal( mk10&mk11&tt10, mk10&mk11&tt11 );

      const auto num_pairs =  static_cast<uint32_t>( eq01 ) + static_cast<uint32_t>( eq02 ) + static_cast<uint32_t>( eq03 ) + 
                              static_cast<uint32_t>( eq12 ) + static_cast<uint32_t>( eq13 ) + static_cast<uint32_t>( eq23 );

      if( num_pairs > 0 )   
      {     
        if ( eq12 ) // F01 = F10 NES
        {
          symmetry_t s66 {0x66, idL, idR};
          s66.remapping_equations( xs, tt, mk );
          res.push_back( s66 );
          symmetry_t s99 {0x99, idL, idR};
          s99.remapping_equations( xs, tt, mk );
          res.push_back( s99 );
        }
        if ( eq03 ) // F00 = F11 ES
        { 
          symmetry_t s33 {0x33, idL, idR}; // 00 -> 11
          s33.remapping_equations( xs, tt, mk );
          res.push_back( s33 );
          symmetry_t sCC {0xCC, idL, idR}; // 11 -> 00
          sCC.remapping_equations( xs, tt, mk );
          res.push_back( sCC );
        }
        if ( eq01 ) // F01=F00
        {
          symmetry_t s11 {0x11, idL, idR}; // 1:00->01
          s11.remapping_equations( xs, tt, mk );
          res.push_back( s11 );
          symmetry_t s44 {0x44, idL, idR}; // 4:01->00
          s44.remapping_equations( xs, tt, mk );
          res.push_back( s44 );
        }
        if ( eq02 ) // F00=F10
        {
          symmetry_t s22 {0x22, idL, idR}; // 2:00->10
          s22.remapping_equations( xs, tt, mk );
          res.push_back( s22 );
          symmetry_t s88 {0x88, idL, idR}; // 8:10->00
          s88.remapping_equations( xs, tt, mk );
          res.push_back( s88 );
        }
        if ( eq13 ) // F01=F11
        {
          symmetry_t s77 {0x77, idL, idR}; // 7:01->11
          s77.remapping_equations( xs, tt, mk );
          res.push_back( s77 );
          symmetry_t sDD {0xDD, idL, idR}; // D:11->01
          sDD.remapping_equations( xs, tt, mk );
          res.push_back( sDD );
        }
        if ( eq23 ) // F11=F10
        {
          symmetry_t sBB {0xBB, idL, idR}; // B:10->11
          sBB.remapping_equations( xs, tt, mk );
          res.push_back( sBB );
          symmetry_t sEE {0xEE, idL, idR}; // E:11->10
          sEE.remapping_equations( xs, tt, mk );
          res.push_back( sEE );

        }
        if ( eq12 && eq03 ) // F01=F10 and F00=F11
        {
          res.emplace_back( 0x36, idL, idR ); // 3:00->11 6:01->10 
          res.emplace_back( 0x6C, idL, idR ); // 6:01->10 C:11->00      
          res.emplace_back( 0x9C, idL, idR ); // 9:10->01 C:11->00     
          res.emplace_back( 0x39, idL, idR ); // 3:00->11 9:10->01
        }
        if( eq02 && eq01 && eq12 )
        {
          res.emplace_back( 0x19, idL, idR ); // 1:00->01 9:10->01
          res.emplace_back( 0x26, idL, idR ); // 2:00->10 6:01->10
        }
        if( eq13 && eq01 && eq03 )
        {
          res.emplace_back( 0x37, idL, idR ); // 3:00->11 7:01->11
          res.emplace_back( 0x4C, idL, idR ); // 4:01->00 C:11->00
        }
        if( eq02 && eq23 && eq03 )
        {
          res.emplace_back( 0x8C, idL, idR ); // 8:10->00 C:11->00
          res.emplace_back( 0x3B, idL, idR ); // 3:00->11 B:10->11
        }
        if( eq13 && eq23 && eq12 )
        {
          res.emplace_back( 0x6E, idL, idR ); // 6:01->10 E:11->10
          res.emplace_back( 0x9D, idL, idR ); // 9:10->01 D:11->01
        }
      }
    }
  }  
  return res;
}



/*
cut_t ccg_analyzer::simple_remap( cut_t cut, node_t xL, node_t xR, uint8_t irem )

cut_t remap( cut_t cut, node_t xL, node_t xR, uint8_t irem )
{
  cut_t res;
  res.set_id( cut.get_id() + 1 );
  for( uint32_t i{0}; i < cut.size(); ++i )
  {
    node_t nd = cut.nodes[i];
    res.add_node( nd.tt, PRJL, nd.id, nd.id );
  }
  uint32_t iL = xL.get_loc_id();
  uint32_t iR = xR.get_loc_id();
  res.nodes[iR].idPi = xR.idPi;
  res.nodes[iL].idPi = xL.idPi;
  res.nodes[iL].idR = res.nodes[iR].idR;
  res.nodes[iR].idL = res.nodes[iL].idL;

  switch ( irem )
  {
    case 0x33:
    {
      res.nodes[iL].gate = OI01;
      res.nodes[iR].gate = OI10;
      res.nodes[iL].tt = ~( ~xL.tt &  xR.tt );
      res.nodes[iR].tt = ~(  xL.tt & ~xR.tt );
      break;
    }
    case 0xCC:
    {
      res.nodes[iL].gate = AI10;
      res.nodes[iR].gate = AI01;
      res.nodes[iL].tt =  xL.tt & ~xR.tt;
      res.nodes[iR].tt = ~xL.tt &  xR.tt;
      break;
    }
    case 0x66:
    {
      res.nodes[iL].gate = OI00;
      res.nodes[iR].gate = AI11;
      res.nodes[iL].tt =  xL.tt | xR.tt;
      res.nodes[iR].tt =  xL.tt & xR.tt;
      break;
    }
    case 0x99:
    {
      res.nodes[iL].gate = AI11;
      res.nodes[iR].gate = OI00;
      res.nodes[iL].tt =  xL.tt & xR.tt;
      res.nodes[iR].tt =  xL.tt | xR.tt;
      break;
    }
    case 0x44:
    {
      res.nodes[iL].gate = PRJL;
      res.nodes[iR].gate = AI11;
      res.nodes[iL].tt = xL.tt;
      res.nodes[iR].tt = xL.tt & xR.tt;
      break;
    }
    case 0x11:
    {
      res.nodes[iL].gate = PRJL;
      res.nodes[iR].gate = OI10;
      res.nodes[iL].tt = xL.tt;
      res.nodes[iR].tt = ~(xL.tt & ~xR.tt);
      break;
    }
    case 0x77:
    {
      res.nodes[iL].gate = OI00;
      res.nodes[iR].gate = PRJR;
      res.nodes[iL].tt = xL.tt | xR.tt;
      res.nodes[iR].tt = xR.tt;
      break;
    }
    case 0xDD:
    {
      res.nodes[iL].gate = AI10;
      res.nodes[iR].gate = PRJR;
      res.nodes[iL].tt = xL.tt & ~xR.tt;
      res.nodes[iR].tt = xR.tt;
      break;
    }
    case 0x88:
    {
      res.nodes[iL].gate = AI11;
      res.nodes[iR].gate = PRJR;
      res.nodes[iL].tt = xL.tt & xR.tt;
      res.nodes[iR].tt = xR.tt;
      break;
    }
    case 0x22:
    {
      res.nodes[iL].gate = OI01;
      res.nodes[iR].gate = PRJR;
      res.nodes[iL].tt = ~( ~xL.tt & xR.tt );
      res.nodes[iR].tt = xR.tt;
      break;
    }
    case 0xBB:
    {
      res.nodes[iL].gate = PRJL;
      res.nodes[iR].gate = OI00;
      res.nodes[iL].tt = xL.tt;
      res.nodes[iR].tt = xL.tt | xR.tt;
      break;
    }
    case 0xEE:
    {
      res.nodes[iL].gate = PRJL;
      res.nodes[iR].gate = AI01;
      res.nodes[iL].tt = xL.tt;
      res.nodes[iR].tt = ~xL.tt & xR.tt;
      break;
    }
    default:
      break;
  }

  
  uint32_t A = ((irem & 0xC0) >> 6u ) & 3u;
  uint32_t C = ((irem & 0x30) >> 4u ) & 3u;
  uint32_t B = ((irem & 0x0C) >> 2u ) & 3u;
  uint32_t D = irem & 3u;

  if( A == B && C == D ) //simple remapping
  {

  }
  else if( A == B ) // compatible remapping
  {

  }
  else // multiform remapping
  {

  }
}


*/

/*cut_t ccg_analyzer::simple_multiform( cut_t cut, node_t xL, node_t xR, uint8_t irem )

cut_t remap( cut_t cut, node_t xL, node_t xR, uint8_t irem )
{
  cut_t res;
  res.set_id( cut.get_id() + 1 );
  for( uint32_t i{0}; i < cut.size(); ++i )
  {
    node_t nd = cut.nodes[i];
    res.add_node( nd.tt, PRJL, nd.id, nd.id );
  }
  uint32_t iL = xL.get_loc_id();
  uint32_t iR = xR.get_loc_id();
  res.nodes[iR].idPi = xR.idPi;
  res.nodes[iL].idPi = xL.idPi;
  res.nodes[iL].idR = res.nodes[iR].idR;
  res.nodes[iR].idL = res.nodes[iL].idL;

  switch ( irem )
  {
    case 0x36:
    {
      res.nodes[iR].gate = XNOR;
      res.nodes[iR].tt = ~xL.tt ^ xR.tt;
      res.nodes.erase( res.nodes.begin() + iL );
      break;
    }
    case 0x6C:
    {
      res.nodes[iL].gate = EXOR;
      res.nodes[iL].tt = xL.tt ^ xR.tt;
      res.nodes.erase( res.nodes.begin() + iR );
      break;
    }
    case 0x9C:
    {
      res.nodes[iR].gate = EXOR;
      res.nodes[iR].tt = xL.tt ^ xR.tt;
      res.nodes.erase( res.nodes.begin() + iL );
      break;
    }
    case 0x39:
    {
      res.nodes[iL].gate = XNOR;
      res.nodes[iL].tt = ~xL.tt ^ xR.tt;
      res.nodes.erase( res.nodes.begin() + iR );
      break;
    }
    case 0x19:
    {
      res.nodes[iL].gate = AI11;
      res.nodes[iL].tt = xL.tt & xR.tt;
      res.nodes.erase( res.nodes.begin() + iR );
      break;
    }
    case 0x26:
    {
      res.nodes[iR].gate = AI11;
      res.nodes[iR].tt = xL.tt & xR.tt;
      res.nodes.erase( res.nodes.begin() + iL );
      break;
    }
    case 0x37:
    {
      res.nodes[iR].gate = OI10;
      res.nodes[iR].tt = ~( xL.tt & ~xR.tt );
      res.nodes.erase( res.nodes.begin() + iL );
      break;
    }
    case 0x4C:
    {
      res.nodes[iR].gate = AI10;
      res.nodes[iR].tt =  xL.tt & ~xR.tt ;
      res.nodes.erase( res.nodes.begin() + iR );
      break;
    }
    case 0x8C:
    {
      res.nodes[iR].gate = AI01;
      res.nodes[iR].tt =  ~xL.tt & xR.tt ;
      res.nodes.erase( res.nodes.begin() + iL );
      break;
    }
    case 0x3B:
    {
      res.nodes[iL].gate = OI01;
      res.nodes[iL].tt =  ~(~xL.tt & xR.tt) ;
      res.nodes.erase( res.nodes.begin() + iR );
      break;
    }
    case 0x6E:
    {
      res.nodes[iL].gate = OI00;
      res.nodes[iL].tt =  xL.tt | xR.tt ;
      res.nodes.erase( res.nodes.begin() + iR );
      break;
    }
    case 0x9D:
    {
      res.nodes[iR].gate = OI00;
      res.nodes[iR].tt =  xL.tt | xR.tt ;
      res.nodes.erase( res.nodes.begin() + iL );
      break;
    }
    default:
      break;
  }

  
  uint32_t A = ((irem & 0xC0) >> 6u ) & 3u;
  uint32_t C = ((irem & 0x30) >> 4u ) & 3u;
  uint32_t B = ((irem & 0x0C) >> 2u ) & 3u;
  uint32_t D = irem & 3u;

  if( A == B && C == D ) //simple remapping
  {

  }
  else if( A == B ) // compatible remapping
  {

  }
  else // multiform remapping
  {

  }
}


*/
/*cut_t ccg_analyzer::simple_compatible( cut_t cut, node_t xL, node_t xR, uint8_t irem )
{
  cut_t res;
  res.set_id( cut.get_id() + 1 );
  for( uint32_t i{0}; i < cut.size(); ++i )
  {
    node_t nd = cut.nodes[i];
    res.add_node( nd.tt, PRJL, nd.id, nd.id );
  }
  uint32_t iL = xL.get_loc_id();
  uint32_t iR = xR.get_loc_id();
  res.nodes[iR].idPi = xR.idPi;
  res.nodes[iL].idPi = xL.idPi;
  res.nodes[iL].idR = res.nodes[iR].idR;
  res.nodes[iR].idL = res.nodes[iL].idL;

  switch ( irem )
  {
    case 0x19:
    {
      res.nodes[iL].gate = AI11;
      res.nodes[iL].tt = xL.tt & xR.tt;
      res.nodes.erase( res.nodes.begin() + iR );
      break;
    }
    case 0x26:
    {
      res.nodes[iR].gate = AI11;
      res.nodes[iR].tt = xL.tt & xR.tt;
      res.nodes.erase( res.nodes.begin() + iL );
      break;
    }
    case 0x37:
    {
      res.nodes[iR].gate = OI10;
      res.nodes[iR].tt = ~( xL.tt & ~xR.tt );
      res.nodes.erase( res.nodes.begin() + iL );
      break;
    }
    case 0x4C:
    {
      res.nodes[iR].gate = AI10;
      res.nodes[iR].tt =  xL.tt & ~xR.tt ;
      res.nodes.erase( res.nodes.begin() + iR );
      break;
    }
    case 0x8C:
    {
      res.nodes[iR].gate = AI01;
      res.nodes[iR].tt =  ~xL.tt & xR.tt ;
      res.nodes.erase( res.nodes.begin() + iL );
      break;
    }
    case 0x3B:
    {
      res.nodes[iL].gate = OI01;
      res.nodes[iL].tt =  ~(~xL.tt & xR.tt) ;
      res.nodes.erase( res.nodes.begin() + iR );
      break;
    }
    case 0x6E:
    {
      res.nodes[iL].gate = OI00;
      res.nodes[iL].tt =  xL.tt | xR.tt ;
      res.nodes.erase( res.nodes.begin() + iR );
      break;
    }
    case 0x9D:
    {
      res.nodes[iR].gate = OI00;
      res.nodes[iR].tt =  xL.tt | xR.tt ;
      res.nodes.erase( res.nodes.begin() + iL );
      break;
    }
    default:
      break;
  }

  
  uint32_t A = ((irem & 0xC0) >> 6u ) & 3u;
  uint32_t C = ((irem & 0x30) >> 4u ) & 3u;
  uint32_t B = ((irem & 0x0C) >> 2u ) & 3u;
  uint32_t D = irem & 3u;

  if( A == B && C == D ) //simple remapping
  {

  }
  else if( A == B ) // compatible remapping
  {

  }
  else // multiform remapping
  {

  }
}*/

cut_t remap( cut_t cut, node_t xL, node_t xR, uint8_t irem )
{
  cut_t res;
  res.set_id( cut.get_id() + 1 );
  for( uint32_t i{0}; i < cut.size(); ++i )
    res.add_node( cut.nodes[i].tt, PRJL, cut.nodes[i].id, cut.nodes[i].id );
  uint32_t A = ((irem & 0xC0) >> 6u ) & 3u;
  uint32_t C = ((irem & 0x30) >> 4u ) & 3u;
  uint32_t B = ((irem & 0x0C) >> 2u ) & 3u;
  uint32_t D = irem & 3u;

  if( A == B && C == D ) //simple remapping
  {
    printf("simple remapping\n");
  }
  else if( A == B ) // compatible remapping
  {
    printf("compatible remapping\n");
  }
  else // multiform remapping
  {
    printf("multiform remapping\n");
  }
  return res;
}

/*! \brief remap */
/*cut_t remap( cut_t cut, node_t xL, node_t xR, uint8_t irem )
{
  cut_t res;
  res.set_id( cut.get_id() + 1 );
  for( uint32_t i{0}; i < cut.size(); ++i )
  {
    node_t nd = cut.nodes[i];
    res.add_node( nd.tt, PRJL, nd.id, nd.id );
  }
  uint32_t iL = xL.get_loc_id();
  uint32_t iR = xR.get_loc_id();
  res.nodes[iR].idPi = xR.idPi;
  res.nodes[iL].idPi = xL.idPi;
  res.nodes[iL].idR = res.nodes[iR].idR;
  res.nodes[iR].idL = res.nodes[iL].idL;

  switch ( irem )
  {
    case 0x33:
    {
      res.nodes[iL].gate = OI01;
      res.nodes[iR].gate = OI10;
      res.nodes[iL].tt = ~( ~xL.tt &  xR.tt );
      res.nodes[iR].tt = ~(  xL.tt & ~xR.tt );
      break;
    }
    case 0xCC:
    {
      res.nodes[iL].gate = AI10;
      res.nodes[iR].gate = AI01;
      res.nodes[iL].tt =  xL.tt & ~xR.tt;
      res.nodes[iR].tt = ~xL.tt &  xR.tt;
      break;
    }
    case 0x66:
    {
      res.nodes[iL].gate = OI00;
      res.nodes[iR].gate = AI11;
      res.nodes[iL].tt =  xL.tt | xR.tt;
      res.nodes[iR].tt =  xL.tt & xR.tt;
      break;
    }
    case 0x99:
    {
      res.nodes[iL].gate = AI11;
      res.nodes[iR].gate = OI00;
      res.nodes[iL].tt =  xL.tt & xR.tt;
      res.nodes[iR].tt =  xL.tt | xR.tt;
      break;
    }
    case 0x44:
    {
      res.nodes[iL].gate = PRJL;
      res.nodes[iR].gate = AI11;
      res.nodes[iL].tt = xL.tt;
      res.nodes[iR].tt = xL.tt & xR.tt;
      break;
    }
    case 0x11:
    {
      res.nodes[iL].gate = PRJL;
      res.nodes[iR].gate = OI10;
      res.nodes[iL].tt = xL.tt;
      res.nodes[iR].tt = ~(xL.tt & ~xR.tt);
      break;
    }
    case 0x77:
    {
      res.nodes[iL].gate = OI00;
      res.nodes[iR].gate = PRJR;
      res.nodes[iL].tt = xL.tt | xR.tt;
      res.nodes[iR].tt = xR.tt;
      break;
    }
    case 0xDD:
    {
      res.nodes[iL].gate = AI10;
      res.nodes[iR].gate = PRJR;
      res.nodes[iL].tt = xL.tt & ~xR.tt;
      res.nodes[iR].tt = xR.tt;
      break;
    }
    case 0x88:
    {
      res.nodes[iL].gate = AI11;
      res.nodes[iR].gate = PRJR;
      res.nodes[iL].tt = xL.tt & xR.tt;
      res.nodes[iR].tt = xR.tt;
      break;
    }
    case 0x22:
    {
      res.nodes[iL].gate = OI01;
      res.nodes[iR].gate = PRJR;
      res.nodes[iL].tt = ~( ~xL.tt & xR.tt );
      res.nodes[iR].tt = xR.tt;
      break;
    }
    case 0xBB:
    {
      res.nodes[iL].gate = PRJL;
      res.nodes[iR].gate = OI00;
      res.nodes[iL].tt = xL.tt;
      res.nodes[iR].tt = xL.tt | xR.tt;
      break;
    }
    case 0xEE:
    {
      res.nodes[iL].gate = PRJL;
      res.nodes[iR].gate = AI01;
      res.nodes[iL].tt = xL.tt;
      res.nodes[iR].tt = ~xL.tt & xR.tt;
      break;
    }
    case 0x36:
    {
      res.nodes[iR].gate = XNOR;
      res.nodes[iR].tt = ~xL.tt ^ xR.tt;
      res.nodes.erase( res.nodes.begin() + iL );
      break;
    }
    case 0x6C:
    {
      res.nodes[iL].gate = EXOR;
      res.nodes[iL].tt = xL.tt ^ xR.tt;
      res.nodes.erase( res.nodes.begin() + iR );
      break;
    }
    case 0x9C:
    {
      res.nodes[iR].gate = EXOR;
      res.nodes[iR].tt = xL.tt ^ xR.tt;
      res.nodes.erase( res.nodes.begin() + iL );
      break;
    }
    case 0x39:
    {
      res.nodes[iL].gate = XNOR;
      res.nodes[iL].tt = ~xL.tt ^ xR.tt;
      res.nodes.erase( res.nodes.begin() + iR );
      break;
    }
    case 0x19:
    {
      res.nodes[iL].gate = AI11;
      res.nodes[iL].tt = xL.tt & xR.tt;
      res.nodes.erase( res.nodes.begin() + iR );
      break;
    }
    case 0x26:
    {
      res.nodes[iR].gate = AI11;
      res.nodes[iR].tt = xL.tt & xR.tt;
      res.nodes.erase( res.nodes.begin() + iL );
      break;
    }
    case 0x37:
    {
      res.nodes[iR].gate = OI10;
      res.nodes[iR].tt = ~( xL.tt & ~xR.tt );
      res.nodes.erase( res.nodes.begin() + iL );
      break;
    }
    case 0x4C:
    {
      res.nodes[iR].gate = AI10;
      res.nodes[iR].tt =  xL.tt & ~xR.tt ;
      res.nodes.erase( res.nodes.begin() + iR );
      break;
    }
    case 0x8C:
    {
      res.nodes[iR].gate = AI01;
      res.nodes[iR].tt =  ~xL.tt & xR.tt ;
      res.nodes.erase( res.nodes.begin() + iL );
      break;
    }
    case 0x3B:
    {
      res.nodes[iL].gate = OI01;
      res.nodes[iL].tt =  ~(~xL.tt & xR.tt) ;
      res.nodes.erase( res.nodes.begin() + iR );
      break;
    }
    case 0x6E:
    {
      res.nodes[iL].gate = OI00;
      res.nodes[iL].tt =  xL.tt | xR.tt ;
      res.nodes.erase( res.nodes.begin() + iR );
      break;
    }
    case 0x9D:
    {
      res.nodes[iR].gate = OI00;
      res.nodes[iR].tt =  xL.tt | xR.tt ;
      res.nodes.erase( res.nodes.begin() + iL );
      break;
    }
    default:
      break;
  }

  
  uint32_t A = ((irem & 0xC0) >> 6u ) & 3u;
  uint32_t C = ((irem & 0x30) >> 4u ) & 3u;
  uint32_t B = ((irem & 0x0C) >> 2u ) & 3u;
  uint32_t D = irem & 3u;

  if( A == B && C == D ) //simple remapping
  {

  }
  else if( A == B ) // compatible remapping
  {

  }
  else // multiform remapping
  {

  }
}

*/
#pragma enxRegion symmetry analysis

#pragma region visulize
void analyzer_t::print_symmetries( std::vector<symmetry_t> sym )
{
  for( auto x : sym )
  {
    switch (x.type)
    {
      case 0x33: printf("l = %2d r = %2d :  ES{ l, r } : l <- nand( l', r )  r <- nand( l , r') : 0x33 : 00->11        : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0xCC: printf("l = %2d r = %2d :  ES{ l, r } : l <- nand( l , r')  r <- nand( l', r ) : 0xCC : 11->00        : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0x66: printf("l = %2d r = %2d : NES{ l, r } : l <-   or( l , r )  r <-  and( l , r ) : 0x66 : 01->10        : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0x99: printf("l = %2d r = %2d : NES{ l, r } : l <-  and( l , r )  r <-   or( l , r ) : 0x99 : 10->01        : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0x44: printf("l = %2d r = %2d : { SVS r }l' : l <- l              r <-  and( l , r ) : 0x44 : 01->00        : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0x11: printf("l = %2d r = %2d : { SVS r }l' : l <- l              r <- nand( l , r') : 0x11 : 00->01        : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0x77: printf("l = %2d r = %2d : { SVS l }r  : l <-   or( l , r )  r <- r             : 0x77 : 01->11        : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0xDD: printf("l = %2d r = %2d : { SVS l }r  : l <-  and( l , r')  r <- r             : 0xDD : 11->01        : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0x88: printf("l = %2d r = %2d : { SVS l }r' : l <-  and( l , r )  r <- r             : 0x88 : 10->00        : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0x22: printf("l = %2d r = %2d : { SVS l }r' : l <- nand( l', r )  r <- r             : 0x22 : 00->10        : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0xBB: printf("l = %2d r = %2d : { SVS r }l  : l <- l              r <-   or( l , r ) : 0xBB : 10->11        : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0xEE: printf("l = %2d r = %2d : { SVS r }l  : l <- l              r <-  and( l', r ) : 0xEE : 11->10        : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0x36: printf("l = %2d r = %2d :  MS{ l, r } : l <- ]              r <- xnor( l , r ) : 0x36 : 00->11 01->10 : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0x6C: printf("l = %2d r = %2d :  MS{ l, r } : l <-  xor( l , r )  r <- ]             : 0x6C : 01->10 11->00 : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0x9C: printf("l = %2d r = %2d :  MS{ l, r } : l <- ]              r <-  xor( l , r ) : 0x9C : 11->00 10->01 : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0x39: printf("l = %2d r = %2d :  MS{ l, r } : l <- xnor( l , r )  r <- ]             : 0x39 : 10->01 00->11 : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0x19: printf("l = %2d r = %2d :CSVS{ l, r } : l <-  and( l , r )  r <- ]             : 0x19 : 00,10->01     : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0x26: printf("l = %2d r = %2d :CSVS{ l, r } : l <- ]              r <-  and( l , r ) : 0x26 : 00,01->10     : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0x37: printf("l = %2d r = %2d :CSVS{ l, r } : l <- ]              r <- nand( l , r') : 0x37 : 00,01->11     : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0x4C: printf("l = %2d r = %2d :CSVS{ l, r } : l <-  and( l , r')  r <- ]             : 0x4C : 01,11->00     : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0x8C: printf("l = %2d r = %2d :CSVS{ l, r } : l <- ]              r <-  and( l', r ) : 0x8C : 10,11->00     : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0x3B: printf("l = %2d r = %2d :CSVS{ l, r } : l <- nand( l', r )  r <- ]             : 0x3B : 00,10->11     : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0x6E: printf("l = %2d r = %2d :CSVS{ l, r } : l <-   or( l , r )  r <- ]             : 0x6E : 01,11->10     : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0x9D: printf("l = %2d r = %2d :CSVS{ l, r } : l <- ]              r <-   or( l , r ) : 0x9D : 10,11->01     : %2d\n", x.idL, x.idR, x.rwd ); break;
    default:
      break;
    }
  }
}
#pragma endregion visualize

} // namespace ccgame

} // namespace mockturtle