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

using TT = kitty::dynamic_truth_table;

TT cube_generator( uint32_t cube, TT Xr, TT Xl )
{
  switch ( cube )
  {
    case 0: return ~Xl & ~Xr; break;
    case 1: return ~Xl &  Xr; break;
    case 2: return  Xl & ~Xr; break;
    case 3: return  Xl &  Xr; break;
  }
}

TT cofactorG( TT fn, uint32_t G, uint32_t idL, uint32_t idR )
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
  TT tt;
  TT mk;
  int rwd;

  symmetry_t(){}
  symmetry_t( uint8_t type, int idL, int idR ): type(type), idL(idL), idR(idR){}

  void remapping_equations( std::vector<TT> * pXs, TT * pTt, TT * pMk )
  {
    uint32_t idA = ((type & 0xC0) >> 6u ) & 3u;
    uint32_t idC = ((type & 0x30) >> 4u ) & 3u;
    uint32_t idB = ((type & 0x0C) >> 2u ) & 3u;
    uint32_t idD = type & 3u;

    TT A = cube_generator( idA, (*pXs)[idR], (*pXs)[idL] );
    TT B = cube_generator( idB, (*pXs)[idR], (*pXs)[idL] );
    TT C = cube_generator( idC, (*pXs)[idR], (*pXs)[idL] );
    TT D = cube_generator( idD, (*pXs)[idR], (*pXs)[idL] );

    TT ttA = cofactorG( (*pTt), idA, idL, idR );
    TT ttB = cofactorG( (*pTt), idB, idL, idR );
    TT ttC = cofactorG( (*pTt), idC, idL, idR );
    TT ttD = cofactorG( (*pTt), idD, idL, idR );

    TT mkA = cofactorG( (*pMk), idA, idL, idR );
    TT mkB = cofactorG( (*pMk), idB, idL, idR );
    TT mkC = cofactorG( (*pMk), idC, idL, idR );
    TT mkD = cofactorG( (*pMk), idD, idL, idR );

    if( idA == idB && idC == idD ) //simple remapping
    {
      mk = ( *pMk & ~A ) | ( C & mkA );
      rwd = kitty::count_zeros( mk );
      TT TA = A & *pTt; /* remapping for A */
      TT TC = C & ( ( mkC & *pTt ) | ( mkA & ttA ) ); /* remapping for C */
      TT TR = ( ~A & ~C & *pTt );
      tt =  TA | TC | TR;
    }
    else if( idA == idB ) // compatible remapping
    {
      mk = ( ~B & ~A & *pMk ) | ( C & ( mkA | mkB ) ) ;
      rwd = kitty::count_zeros( mk ); //@ltry
      TT TA = A & *pTt;
      TT TB = B & *pTt;
      TT TC = C & ( ( mkA & ttA ) | ( mkB & ttB ) | ( mkC & *pTt ) );
      TT TR = ~A & ~B & ~C & *pTt;
      tt = TA | TB | TC | TR;
    }
    else // multiform remapping
    {
      mk = ( ~B & ~A & *pMk ) | ( ( C & mkA ) | ( D & mkB ) );
      rwd = kitty::count_zeros( mk ); // @ltry
      
      TT preserved = (~A&~B&~C&~D)&*pTt;
      TT modifiedA = A & *pTt;
      TT modifiedB = B & *pTt;
      TT modifiedC = C & ( ( mkA & ~mkC & ttA ) | ( mkC & *pTt ) );
      TT modifiedD = D & ( ( mkB & ~mkD & ttB ) | ( mkD & *pTt ) );
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
  std::vector<symmetry_t> find_symmetries( std::vector<TT> *, TT *, TT *, std::vector<int> * );
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

/*! \brief symmetry analysis of variable iR with all the ones on the right */
std::vector<symmetry_t> analyzer_t::find_symmetries( std::vector<TT> * pXs, TT * pTt, TT * pMk, std::vector<int> * pIds )
{
  std::vector<symmetry_t> res;
  int idL, idR;
  for( int iR{0}; iR < pIds->size(); ++iR )
  {
    idR = (*pIds)[iR];
    if( idR < 0 ) continue; 
    TT tt0  = kitty::cofactor0( *pTt, idR );
    TT tt1  = kitty::cofactor1( *pTt, idR );
    TT mk0  = kitty::cofactor0( *pMk, idR );
    TT mk1  = kitty::cofactor1( *pMk, idR );
    /* todo: add top-decomposition check */

    /* symmetry check */
    for( int iL{iR+1}; iL < pIds->size(); ++iL )
    {
      idL = (*pIds)[iL];
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
          s66.remapping_equations( pXs, pTt, pMk );
          res.push_back( s66 );
          symmetry_t s99 {0x99, idL, idR};
          s99.remapping_equations( pXs, pTt, pMk );
          res.push_back( s99 );
        }
        if ( eq03 ) // F00 = F11 ES
        { 
          symmetry_t s33 {0x33, idL, idR}; // 00 -> 11
          s33.remapping_equations( pXs, pTt, pMk );
          res.push_back( s33 );
          symmetry_t sCC {0xCC, idL, idR}; // 11 -> 00
          sCC.remapping_equations( pXs, pTt, pMk );
          res.push_back( sCC );
        }
        if ( eq01 ) // F01=F00
        {
          symmetry_t s11 {0x11, idL, idR}; // 1:00->01
          s11.remapping_equations( pXs, pTt, pMk );
          res.push_back( s11 );
          symmetry_t s44 {0x44, idL, idR}; // 4:01->00
          s44.remapping_equations( pXs, pTt, pMk );
          res.push_back( s44 );
        }
        if ( eq02 ) // F00=F10
        {
          symmetry_t s22 {0x22, idL, idR}; // 2:00->10
          s22.remapping_equations( pXs, pTt, pMk );
          res.push_back( s22 );
          symmetry_t s88 {0x88, idL, idR}; // 8:10->00
          s88.remapping_equations( pXs, pTt, pMk );
          res.push_back( s88 );
        }
        if ( eq13 ) // F01=F11
        {
          symmetry_t s77 {0x77, idL, idR}; // 7:01->11
          s77.remapping_equations( pXs, pTt, pMk );
          res.push_back( s77 );
          symmetry_t sDD {0xDD, idL, idR}; // D:11->01
          sDD.remapping_equations( pXs, pTt, pMk );
          res.push_back( sDD );
        }
        if ( eq23 ) // F11=F10
        {
          symmetry_t sBB {0xBB, idL, idR}; // B:10->11
          sBB.remapping_equations( pXs, pTt, pMk );
          res.push_back( sBB );
          symmetry_t sEE {0xEE, idL, idR}; // E:11->10
          sEE.remapping_equations( pXs, pTt, pMk );
          res.push_back( sEE );
        }
        if ( eq12 && eq03 ) // F01=F10 and F00=F11
        {
          symmetry_t s36 {0x36, idL, idR}; // 3:00->11 6:01->10
          s36.remapping_equations( pXs, pTt, pMk );
          res.push_back( s36 );
          symmetry_t s6C {0x6C, idL, idR}; // 6:01->10 C:11->00
          s6C.remapping_equations( pXs, pTt, pMk );
          res.push_back( s6C );
          symmetry_t s9C {0x9C, idL, idR}; // 9:10->01 C:11->00
          s9C.remapping_equations( pXs, pTt, pMk );
          res.push_back( s9C );
          symmetry_t s39 {0x39, idL, idR}; // 3:00->11 9:10->01
          s39.remapping_equations( pXs, pTt, pMk );
          res.push_back( s39 );
        }
        if( eq02 && eq01 && eq12 )
        {
          symmetry_t s19 {0x19, idL, idR}; // 1:00->01 9:10->01
          s19.remapping_equations( pXs, pTt, pMk );
          res.push_back( s19 );
          symmetry_t s26 {0x26, idL, idR}; // 2:00->10 6:01->10
          s26.remapping_equations( pXs, pTt, pMk );
          res.push_back( s26 );
        }
        if( eq13 && eq01 && eq03 )
        {
          symmetry_t s37 {0x37, idL, idR}; // 3:00->11 7:01->11
          s37.remapping_equations( pXs, pTt, pMk );
          res.push_back( s37 );
          symmetry_t s4C {0x4C, idL, idR}; // 4:01->00 C:11->00
          s4C.remapping_equations( pXs, pTt, pMk );
          res.push_back( s4C );
        }
        if( eq02 && eq23 && eq03 )
        {
          symmetry_t s8C {0x8C, idL, idR}; // 8:10->00 C:11->00
          s8C.remapping_equations( pXs, pTt, pMk );
          res.push_back( s8C );
          symmetry_t s3B {0x3B, idL, idR}; // 3:00->11 B:10->11
          s3B.remapping_equations( pXs, pTt, pMk );
          res.push_back( s3B );
        }
        if( eq13 && eq23 && eq12 )
        {
          symmetry_t s6E {0x6E, idL, idR}; // 6:01->10 E:11->10
          s6E.remapping_equations( pXs, pTt, pMk );
          res.push_back( s6E );
          symmetry_t s9D {0x9D, idL, idR}; // 9:10->01 D:11->01
          s9D.remapping_equations( pXs, pTt, pMk );
          res.push_back( s9D );
        }
      }
    }
  }  
  return res;
}

#pragma enxRegion symmetry analysis

#pragma region visulize
void analyzer_t::print_symmetries( std::vector<symmetry_t> sym )
{
  for( auto x : sym )
  {
    switch (x.type)
    {
      case 0x33: printf("l = %2d r = %2d :  ES{ l, r } : l <- nand( l', r )  r <- nand( l , r') : 0x33 : 00->11        : %2d\n", x.idL, x.idR, x.rwd ); break;
      case 0xCC: printf("l = %2d r = %2d :  ES{ l, r } : l <-  and( l , r')  r <-  and( l', r ) : 0xCC : 11->00        : %2d\n", x.idL, x.idR, x.rwd ); break;
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