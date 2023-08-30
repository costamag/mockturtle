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
  \file sym_synthesis.hpp
  \brief Symmetry-based synthesis

  \author Andrea Costamagna
*/
#pragma once

#include <kitty/partial_truth_table.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/bit_operations.hpp>
#include <kitty/decomposition.hpp>
#include <kitty/constructors.hpp>
#include <kitty/operators.hpp>
#include <kitty/print.hpp>
#include <fmt/format.h>
#include <stdio.h>
#include <stack>

#include "../../networks/aig.hpp"
#include "../../networks/xag.hpp"
#include "../../utils/stopwatch.hpp"


namespace mockturtle
{

namespace techaware
{

void dprintf( std::string s )
{
  printf( "%s\n", s.c_str() );
}

using TT = kitty::dynamic_truth_table;
uint32_t UNK32 = 0x0FFFFFFF;

template<class Ntk>
uint32_t delay_xor()
{
  return std::is_same<Ntk, xag_network>::value ? 1 : 2;
}

/*! \brief Gate type in the techaware namespace. Convention Xl=1100, Xr=1010 */
enum gate_t : uint8_t
{
  PIS_ = 0xF0,
  CNTR = 0x0, // 0000
  PA00 = 0x1, // 0001
  PA01 = 0x2, // 0010
  CMPL = 0x3, // 0011
  PA10 = 0x4, // 0100
  CMPR = 0X5, // 0101
  EXOR = 0X6, // 0110
  IA11 = 0X7, // 0111
  PA11 = 0X8, // 1000
  XNOR = 0X9, // 1001
  PRJR = 0XA, // 1010
  IA10 = 0XB, // 1011
  PRJL = 0XC, // 1100
  IA01 = 0XD, // 1101
  IA00 = 0XE, // 1110
  TAUT = 0XF, // 1111
  TARG = 0xFF
};

/*! \brief Gate type in the techaware namespace. Convention Xl=1100, Xr=1010 */
enum dec_t : uint8_t
{
  AND = 0x8,
  OR_ = 0xE,
  LE_ = 0xB,
  LT_ = 0x2,
  XOR = 0x6,
  NUL = 0x0
};


#pragma region SYMMETRIES

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

struct symmetry_t
{
  uint8_t type;
  int piL;
  int piR;
  TT tTt;
  TT tMk;
  int reward;

  symmetry_t(){}
  symmetry_t( uint8_t type, int idL, int idR ): type(type), piL(idL), piR(idR){}

  void remapping_equations( std::vector<TT> * pXs, TT * pTt, TT * pMk )
  {
    uint32_t idA = ((type & 0xC0) >> 6u ) & 3u;
    uint32_t idC = ((type & 0x30) >> 4u ) & 3u;
    uint32_t idB = ((type & 0x0C) >> 2u ) & 3u;
    uint32_t idD = type & 3u;

    TT A = cube_generator( idA, (*pXs)[piR], (*pXs)[piL] );
    TT B = cube_generator( idB, (*pXs)[piR], (*pXs)[piL] );
    TT C = cube_generator( idC, (*pXs)[piR], (*pXs)[piL] );
    TT D = cube_generator( idD, (*pXs)[piR], (*pXs)[piL] );

    TT ttA = cofactorG( (*pTt), idA, piL, piR );
    TT ttB = cofactorG( (*pTt), idB, piL, piR );
    TT ttC = cofactorG( (*pTt), idC, piL, piR );
    TT ttD = cofactorG( (*pTt), idD, piL, piR );

    TT mkA = cofactorG( (*pMk), idA, piL, piR );
    TT mkB = cofactorG( (*pMk), idB, piL, piR );
    TT mkC = cofactorG( (*pMk), idC, piL, piR );
    TT mkD = cofactorG( (*pMk), idD, piL, piR );

    if( idA == idB && idC == idD ) //simple remapping
    {
      tMk = ( *pMk & ~A ) | ( C & mkA );
      reward = kitty::count_zeros( tMk );
      TT TA = A & *pTt; /* remapping for A */
      TT TC = C & ( ( mkC & *pTt ) | ( mkA & ttA ) ); /* remapping for C */
      TT TR = ( ~A & ~C & *pTt );
      tTt =  TA | TC | TR;
    }
    else if( idA == idB ) // compatible remapping
    {
      tMk = ( ~B & ~A & *pMk ) | ( C & ( mkA | mkB ) ) ;
      reward = kitty::count_zeros( tMk ); //@ltry
      TT TA = A & *pTt;
      TT TB = B & *pTt;
      TT TC = C & ( ( mkA & ttA ) | ( mkB & ttB ) | ( mkC & *pTt ) );
      TT TR = ~A & ~B & ~C & *pTt;
      tTt = TA | TB | TC | TR;
    }
    else // multiform remapping
    {
      tMk = ( ~B & ~A & *pMk ) | ( ( C & mkA ) | ( D & mkB ) );
      reward = kitty::count_zeros( tMk ); // @ltry
      
      TT preserved = (~A&~B&~C&~D)&*pTt;
      TT modifiedA = A & *pTt;
      TT modifiedB = B & *pTt;
      TT modifiedC = C & ( ( mkA & ~mkC & ttA ) | ( mkC & *pTt ) );
      TT modifiedD = D & ( ( mkB & ~mkD & ttB ) | ( mkD & *pTt ) );
      tTt = preserved | modifiedA | modifiedB | modifiedC | modifiedD;
    }
  }
};


struct decomposition_t
{
  uint8_t type;
  int pi;
  TT tTt;
  TT tMk;

  decomposition_t(){}
  ~decomposition_t(){}

};



void print_symmetries( std::vector<symmetry_t> sym )
{
  for( auto x : sym )
  {
    switch (x.type)
    {
      case 0x33: printf("l = %2d r = %2d :  ES{ l, r } : l <- nand( l', r )  r <- nand( l , r') : 0x33 : 00->11        : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0xCC: printf("l = %2d r = %2d :  ES{ l, r } : l <-  and( l , r')  r <-  and( l', r ) : 0xCC : 11->00        : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0x66: printf("l = %2d r = %2d : NES{ l, r } : l <-   or( l , r )  r <-  and( l , r ) : 0x66 : 01->10        : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0x99: printf("l = %2d r = %2d : NES{ l, r } : l <-  and( l , r )  r <-   or( l , r ) : 0x99 : 10->01        : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0x44: printf("l = %2d r = %2d : { SVS r }l' : l <- l              r <-  and( l , r ) : 0x44 : 01->00        : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0x11: printf("l = %2d r = %2d : { SVS r }l' : l <- l              r <- nand( l , r') : 0x11 : 00->01        : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0x77: printf("l = %2d r = %2d : { SVS l }r  : l <-   or( l , r )  r <- r             : 0x77 : 01->11        : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0xDD: printf("l = %2d r = %2d : { SVS l }r  : l <-  and( l , r')  r <- r             : 0xDD : 11->01        : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0x88: printf("l = %2d r = %2d : { SVS l }r' : l <-  and( l , r )  r <- r             : 0x88 : 10->00        : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0x22: printf("l = %2d r = %2d : { SVS l }r' : l <- nand( l', r )  r <- r             : 0x22 : 00->10        : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0xBB: printf("l = %2d r = %2d : { SVS r }l  : l <- l              r <-   or( l , r ) : 0xBB : 10->11        : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0xEE: printf("l = %2d r = %2d : { SVS r }l  : l <- l              r <-  and( l', r ) : 0xEE : 11->10        : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0x36: printf("l = %2d r = %2d :  MS{ l, r } : l <- ]              r <- xnor( l , r ) : 0x36 : 00->11 01->10 : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0x6C: printf("l = %2d r = %2d :  MS{ l, r } : l <-  xor( l , r )  r <- ]             : 0x6C : 01->10 11->00 : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0x9C: printf("l = %2d r = %2d :  MS{ l, r } : l <- ]              r <-  xor( l , r ) : 0x9C : 11->00 10->01 : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0x39: printf("l = %2d r = %2d :  MS{ l, r } : l <- xnor( l , r )  r <- ]             : 0x39 : 10->01 00->11 : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0x19: printf("l = %2d r = %2d :CSVS{ l, r } : l <-  and( l , r )  r <- ]             : 0x19 : 00,10->01     : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0x26: printf("l = %2d r = %2d :CSVS{ l, r } : l <- ]              r <-  and( l , r ) : 0x26 : 00,01->10     : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0x37: printf("l = %2d r = %2d :CSVS{ l, r } : l <- ]              r <- nand( l , r') : 0x37 : 00,01->11     : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0x4C: printf("l = %2d r = %2d :CSVS{ l, r } : l <-  and( l , r')  r <- ]             : 0x4C : 01,11->00     : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0x8C: printf("l = %2d r = %2d :CSVS{ l, r } : l <- ]              r <-  and( l', r ) : 0x8C : 10,11->00     : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0x3B: printf("l = %2d r = %2d :CSVS{ l, r } : l <- nand( l', r )  r <- ]             : 0x3B : 00,10->11     : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0x6E: printf("l = %2d r = %2d :CSVS{ l, r } : l <-   or( l , r )  r <- ]             : 0x6E : 01,11->10     : %2d\n", x.piL, x.piR, x.reward ); break;
      case 0x9D: printf("l = %2d r = %2d :CSVS{ l, r } : l <- ]              r <-   or( l , r ) : 0x9D : 10,11->01     : %2d\n", x.piL, x.piR, x.reward ); break;
    default:
      break;
    }
  }
}

#pragma endregion SYMMETRIES

struct node_t
{
  /*! \brief simulation pattern */
  TT sTt;
  /*! \brief simulation mask */
  TT sMk;
  /*! \brief [8  bits: gate type] */
  gate_t gate;
  /*! \brief [16 bits level id][16 bits node id] */
  uint32_t id;
  /*! \brief [32 bits: left-fanin identifier] */
  uint32_t idL;
  /*! \brief [32 bits: right-fanin identifier] */
  uint32_t idR;
  /*! \brief [1 bit NOT remapped 32 bits remapped pi] */
  uint32_t idPi{0x80000000};
  /*! \brief delay */
  int level;

  node_t(){}
  ~node_t(){}

  node_t( gate_t gate, uint32_t level, TT sim_tt, TT sim_mk, uint32_t cut_id, uint32_t ref_id ) : gate(gate), level(level), sTt(sim_tt), sMk(sim_mk)
  {
    id = ( cut_id << 16u ) | ref_id;
  }

  gate_t   get_gate_type()  { return gate; }
  uint32_t get_this_ref_id(){ return id & 0x0000FFFF;  }
  uint32_t get_this_cut_id(){ return 0x0000FFFF & (( id & 0xFFFF0000 ) >> 16u);  }
  uint32_t get_linp_ref_id(){ return idL & 0x0000FFFF;  }
  uint32_t get_linp_cut_id(){ return 0x0000FFFF & (( idL & 0xFFFF0000 ) >> 16u);  }
  uint32_t get_rinp_ref_id(){ return idR & 0x0000FFFF;  }
  uint32_t get_rinp_cut_id(){ return 0x0000FFFF & (( idR & 0xFFFF0000 ) >> 16u);  }
};

struct cut_t
{
public:
  uint32_t id{0};
  std::vector<node_t> nodes;
  uint32_t nNodes{0};
  TT tTt;
  TT tMk;
  uint32_t tCut;
  uint32_t tRef;
  std::vector<uint32_t> pi_to_node;
  uint32_t delayCost{0};

  cut_t(){};
  cut_t( uint32_t id ) : id(id) {};
  ~cut_t(){};

  void update_cut_id( uint32_t idCutNew )
  {
    id = idCutNew;
    for( uint32_t iNd{0}; iNd<nodes.size(); ++iNd )
      nodes[iNd].id = (nodes[iNd].id & 0x0000FFFF) | ( (idCutNew << 16u ) & 0xFFFF0000 );
  }

  uint32_t add_node( gate_t gate, uint32_t level, TT sTt, TT sMk )
  { 
    nodes.emplace_back( gate, level, sTt, sMk, id, nNodes++ );
    if( level > delayCost ) delayCost = level;
    return nodes.back().get_this_ref_id();
  }

  uint32_t add_node( gate_t gate, uint32_t level, TT sTt, TT sMk, uint32_t cut_id, uint32_t idPi, uint32_t idL, uint32_t idR )
  { 
    node_t nd( gate, level, sTt, sMk, id, nNodes++ );
    nd.idL = idL;
    nd.idR = idR;
    nd.idPi = idPi;
    
    if( pi_to_node.size() <= idPi )
      for( uint32_t i{pi_to_node.size()}; i<idPi+1; ++i )
        pi_to_node.push_back(UNK32);
    pi_to_node[idPi]=nd.get_this_ref_id();

    if( level > delayCost ) delayCost = level;
    nodes.push_back( nd );
    return nd.id;
  }

  template<class Ntk>
  void add_node_symL( cut_t * pPrevCut, symmetry_t * pSym )
  {
    uint32_t iL = pPrevCut->pi_to_node[pSym->piL];
    uint32_t iR = pPrevCut->pi_to_node[pSym->piR];
    node_t xL = pPrevCut->nodes[iL];
    node_t xR = pPrevCut->nodes[iR];

    switch ( pSym->type )
    {
      case 0x33: add_node( IA01, std::max(xL.level, xR.level)+1, ~( ~xL.sTt &  xR.sTt ), xL.sMk | xR.sMk, id, pSym->piL, xL.id, xR.id ); break;  //nand( l', r )
      case 0xCC: add_node( PA10, std::max(xL.level, xR.level)+1,  (  xL.sTt & ~xR.sTt ), xL.sMk | xR.sMk, id, pSym->piL, xL.id, xR.id ); break;  // and( l , r')
      case 0x66: add_node( IA00, std::max(xL.level, xR.level)+1, ~( ~xL.sTt & ~xR.sTt ), xL.sMk | xR.sMk, id, pSym->piL, xL.id, xR.id ); break;  //or( l , r )
      case 0x99: add_node( PA11, std::max(xL.level, xR.level)+1,  (  xL.sTt &  xR.sTt ), xL.sMk | xR.sMk, id, pSym->piL, xL.id, xR.id ); break;  //and( l , r )
      case 0x44: add_node( PRJL, xL.level, xL.sTt, xL.sMk | xR.sMk, id, pSym->piL, xL.id, xL.id ); break;  // l            
      case 0x11: add_node( PRJL, xL.level, xL.sTt, xL.sMk | xR.sMk, id, pSym->piL, xL.id, xL.id ); break;  // l            
      case 0x77: add_node( IA00, std::max(xL.level, xR.level)+1, ~( ~xL.sTt & ~xR.sTt ), xL.sMk | xR.sMk, id, pSym->piL, xL.id, xR.id ); break;  //   or( l , r )
      case 0xDD: add_node( PA10, std::max(xL.level, xR.level)+1,  (  xL.sTt & ~xR.sTt ), xL.sMk | xR.sMk, id, pSym->piL, xL.id, xR.id ); break;  //  and( l , r')
      case 0x88: add_node( PA11, std::max(xL.level, xR.level)+1,  (  xL.sTt &  xR.sTt ), xL.sMk | xR.sMk, id, pSym->piL, xL.id, xR.id ); break;  //  and( l , r )
      case 0x22: add_node( IA01, std::max(xL.level, xR.level)+1, ~( ~xL.sTt &  xR.sTt ), xL.sMk | xR.sMk, id, pSym->piL, xL.id, xR.id ); break;  // nand( l', r )
      case 0xBB: add_node( PRJL, xL.level, xL.sTt, xL.sMk | xR.sMk, id, pSym->piL, xL.id, xL.id ); break;  // l            
      case 0xEE: add_node( PRJL, xL.level, xL.sTt, xL.sMk | xR.sMk, id, pSym->piL, xL.id, xL.id ); break;  // l            
      case 0x36: break;  // ]            
      case 0x6C: add_node( EXOR, std::max(xL.level, xR.level)+delay_xor<Ntk>(), xL.sTt ^ xR.sTt, xL.sMk | xR.sMk, id, pSym->piL, xL.id, xR.id ); break;  //  xor( l , r )
      case 0x9C: break;  // ]            
      case 0x39: add_node( XNOR, std::max(xL.level, xR.level)+delay_xor<Ntk>(), ~(  xL.sTt ^  xR.sTt ), xL.sMk | xR.sMk, id, pSym->piL, xL.id, xR.id ); break;  // xnor( l , r )
      case 0x19: add_node( PA11, std::max(xL.level, xR.level)+1, (  xL.sTt &  xR.sTt ), xL.sMk | xR.sMk, id, pSym->piL, xL.id, xR.id ); break;  //  and( l , r )
      case 0x26: break;  // ]            
      case 0x37: break;  // ]            
      case 0x4C: add_node( PA10, std::max(xL.level, xR.level)+1, (  xL.sTt & ~xR.sTt ), xL.sMk | xR.sMk, id, pSym->piL, xL.id, xR.id ); break;  //  and( l , r')
      case 0x8C: break;  // ]            
      case 0x3B: add_node( PA01, std::max(xL.level, xR.level)+1,  ( ~xL.sTt &  xR.sTt ), xL.sMk | xR.sMk, id, pSym->piL, xL.id, xR.id ); break;  // nand( l', r )
      case 0x6E: add_node( IA00, std::max(xL.level, xR.level)+1, ~( ~xL.sTt & ~xR.sTt ), xL.sMk | xR.sMk, id, pSym->piL, xL.id, xR.id ); break;  //   or( l , r )
      case 0x9D: break;  // ]            
    }
  }

  template<class Ntk>
  void add_node_symR( cut_t * pPrevCut, symmetry_t * pSym )
  {
    uint32_t iL = pPrevCut->pi_to_node[pSym->piL];
    uint32_t iR = pPrevCut->pi_to_node[pSym->piR];
    node_t xL = pPrevCut->nodes[iL];
    node_t xR = pPrevCut->nodes[iR];

    switch ( pSym->type )
    {
      case 0x33: add_node( IA10, std::max(xL.level, xR.level)+1, ~(  xL.sTt & ~xR.sTt ), xL.sMk | xR.sMk, id, pSym->piR, xL.id, xR.id ); break;// nand( l , r')
      case 0xCC: add_node( PA01, std::max(xL.level, xR.level)+1,  ( ~xL.sTt &  xR.sTt ), xL.sMk | xR.sMk, id, pSym->piR, xL.id, xR.id ); break;//  and( l', r )
      case 0x66: add_node( PA11, std::max(xL.level, xR.level)+1,  (  xL.sTt &  xR.sTt ), xL.sMk | xR.sMk, id, pSym->piR, xL.id, xR.id ); break;//  and( l , r )
      case 0x99: add_node( IA00, std::max(xL.level, xR.level)+1, ~( ~xL.sTt & ~xR.sTt ), xL.sMk | xR.sMk, id, pSym->piR, xL.id, xR.id ); break;//   or( l , r )
      case 0x44: add_node( PA11, std::max(xL.level, xR.level)+1,  (  xL.sTt &  xR.sTt ), xL.sMk | xR.sMk, id, pSym->piR, xL.id, xR.id ); break;//  and( l , r )
      case 0x11: add_node( IA10, std::max(xL.level, xR.level)+1, ~(  xL.sTt & ~xR.sTt ), xL.sMk | xR.sMk, id, pSym->piR, xL.id, xR.id ); break;// nand( l , r')
      case 0x77: add_node( PRJR, xR.level, xR.sTt, xL.sMk | xR.sMk, id, pSym->piR, xR.id, xR.id ); break;// r            
      case 0xDD: add_node( PRJR, xR.level, xR.sTt, xL.sMk | xR.sMk, id, pSym->piR, xR.id, xR.id ); break;// r            
      case 0x88: add_node( PRJR, xR.level, xR.sTt, xL.sMk | xR.sMk, id, pSym->piR, xR.id, xR.id ); break;// r            
      case 0x22: add_node( PRJR, xR.level, xR.sTt, xL.sMk | xR.sMk, id, pSym->piR, xR.id, xR.id ); break;// r            
      case 0xBB: add_node( IA00, std::max(xL.level, xR.level)+1, ~( ~xL.sTt & ~xR.sTt ), xL.sMk | xR.sMk, id, pSym->piR, xL.id, xR.id ); break;//   or( l , r )
      case 0xEE: add_node( PA01, std::max(xL.level, xR.level)+1,  ( ~xL.sTt &  xR.sTt ), xL.sMk | xR.sMk, id, pSym->piR, xL.id, xR.id ); break;//  and( l', r )
      case 0x36: add_node( XNOR, std::max(xL.level,xR.level)+delay_xor<Ntk>(), ~( xL.sTt ^  xR.sTt ), xL.sMk | xR.sMk, id, pSym->piR, xL.id, xR.id ); break;// xnor( l , r )
      case 0x6C: break;// ]            
      case 0x9C: add_node( EXOR, std::max(xL.level, xR.level)+delay_xor<Ntk>(), ( xL.sTt ^  xR.sTt ), xL.sMk | xR.sMk, id, pSym->piR, xL.id, xR.id ); break;//  xor( l , r )
      case 0x39: break;// ]            
      case 0x19: break;// ]            
      case 0x26: add_node( PA11, std::max(xL.level, xR.level)+1,  (  xL.sTt &  xR.sTt ), xL.sMk | xR.sMk, id, pSym->piR, xL.id, xR.id ); break;//  and( l , r )
      case 0x37: add_node( IA10, std::max(xL.level, xR.level)+1, ~(  xL.sTt & ~xR.sTt ), xL.sMk | xR.sMk, id, pSym->piR, xL.id, xR.id ); break;// nand( l , r')
      case 0x4C: break;// ]            
      case 0x8C: add_node( PA01, std::max(xL.level, xR.level)+1, ( ~xL.sTt &  xR.sTt ), xL.sMk | xR.sMk, id, pSym->piR, xL.id, xR.id ); break;//  and( l', r )
      case 0x3B: break;// ]            
      case 0x6E: break;// ]            
      case 0x9D: add_node( IA00, std::max(xL.level, xR.level)+1, ~( ~xL.sTt & ~xR.sTt ), xL.sMk | xR.sMk, id, pSym->piR, xL.id, xR.id ); break;//   or( l , r )
    }

  }

  void set_target( TT func, TT mask, uint32_t idCutTrg, uint32_t idRefTrg )
  {
    tTt = func;
    tMk = mask;
    tCut = idCutTrg;
    tRef = idRefTrg;
  }

  void erase_node_from_pi( uint32_t idPi )
  {
    uint32_t idNd = pi_to_node[idPi];
    nodes.erase( nodes.begin() + idNd );
    pi_to_node[idPi] = UNK32;
  }

  void fill_pi_to_node()
  {
    uint32_t idPi;
    for( uint32_t iNd{0}; iNd<nodes.size(); ++iNd )
    {
      idPi = nodes[iNd].idPi;
      if( pi_to_node.size() <= idPi )
        for( uint32_t i{pi_to_node.size()}; i<idPi+1; ++i )
          pi_to_node.push_back(UNK32);
      pi_to_node[idPi]=iNd;
    }
  }

  std::vector<symmetry_t> symmetry_analysis( std::vector<TT> * pX )
  {
    std::vector<symmetry_t> res;
    int piL, piR;
    for( int iR{0}; iR < nNodes-1; ++iR )
    {
      piR = nodes[iR].idPi;
      TT tt0  = kitty::cofactor0( tTt, piR );
      TT tt1  = kitty::cofactor1( tTt, piR );
      TT mk0  = kitty::cofactor0( tMk, piR );
      TT mk1  = kitty::cofactor1( tMk, piR );

      /* symmetry check */
      for( int iL{iR+1}; iL < nNodes; ++iL )
      {
        piL = nodes[iL].idPi;
        assert( piL > piR );

        const auto tt00 = kitty::cofactor0( tt0, piL );
        const auto tt01 = kitty::cofactor0( tt1, piL );
        const auto tt10 = kitty::cofactor1( tt0, piL );
        const auto tt11 = kitty::cofactor1( tt1, piL );
        const auto mk00 = kitty::cofactor0( mk0, piL );
        const auto mk01 = kitty::cofactor0( mk1, piL );
        const auto mk10 = kitty::cofactor1( mk0, piL );
        const auto mk11 = kitty::cofactor1( mk1, piL );

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
            symmetry_t s66 {0x66, piL, piR};
            s66.remapping_equations( pX, &tTt, &tMk );
            res.push_back( s66 );
            symmetry_t s99 {0x99, piL, piR};
            s99.remapping_equations( pX, &tTt, &tMk );
            res.push_back( s99 );
          }
          if ( eq03 ) // F00 = F11 ES
          { 
            symmetry_t s33 {0x33, piL, piR}; // 00 -> 11
            s33.remapping_equations( pX, &tTt, &tMk );
            res.push_back( s33 );
            symmetry_t sCC {0xCC, piL, piR}; // 11 -> 00
            sCC.remapping_equations( pX, &tTt, &tMk );
            res.push_back( sCC );
          }
          if ( eq01 ) // F01=F00
          {
            symmetry_t s11 {0x11, piL, piR}; // 1:00->01
            s11.remapping_equations( pX, &tTt, &tMk );
            res.push_back( s11 );
            symmetry_t s44 {0x44, piL, piR}; // 4:01->00
            s44.remapping_equations( pX, &tTt, &tMk );
            res.push_back( s44 );
          }
          if ( eq02 ) // F00=F10
          {
            symmetry_t s22 {0x22, piL, piR}; // 2:00->10
            s22.remapping_equations( pX, &tTt, &tMk );
            res.push_back( s22 );
            symmetry_t s88 {0x88, piL, piR}; // 8:10->00
            s88.remapping_equations( pX, &tTt, &tMk );
            res.push_back( s88 );
          }
          if ( eq13 ) // F01=F11
          {
            symmetry_t s77 {0x77, piL, piR}; // 7:01->11
            s77.remapping_equations( pX, &tTt, &tMk );
            res.push_back( s77 );
            symmetry_t sDD {0xDD, piL, piR}; // D:11->01
            sDD.remapping_equations( pX, &tTt, &tMk );
            res.push_back( sDD );
          }
          if ( eq23 ) // F11=F10
          {
            symmetry_t sBB {0xBB, piL, piR}; // B:10->11
            sBB.remapping_equations( pX, &tTt, &tMk );
            res.push_back( sBB );
            symmetry_t sEE {0xEE, piL, piR}; // E:11->10
            sEE.remapping_equations( pX, &tTt, &tMk );
            res.push_back( sEE );
          }
          if ( eq12 && eq03 ) // F01=F10 and F00=F11
          {
            symmetry_t s36 {0x36, piL, piR}; // 3:00->11 6:01->10
            s36.remapping_equations( pX, &tTt, &tMk );
            res.push_back( s36 );
            symmetry_t s6C {0x6C, piL, piR}; // 6:01->10 C:11->00
            s6C.remapping_equations( pX, &tTt, &tMk );
            res.push_back( s6C );
            symmetry_t s9C {0x9C, piL, piR}; // 9:10->01 C:11->00
            s9C.remapping_equations( pX, &tTt, &tMk );
            res.push_back( s9C );
            symmetry_t s39 {0x39, piL, piR}; // 3:00->11 9:10->01
            s39.remapping_equations( pX, &tTt, &tMk );
            res.push_back( s39 );
          }
          if( eq02 && eq01 && eq12 )
          {
            symmetry_t s19 {0x19, piL, piR}; // 1:00->01 9:10->01
            s19.remapping_equations( pX, &tTt, &tMk );
            res.push_back( s19 );
            symmetry_t s26 {0x26, piL, piR}; // 2:00->10 6:01->10
            s26.remapping_equations( pX, &tTt, &tMk );
            res.push_back( s26 );
          }
          if( eq13 && eq01 && eq03 )
          {
            symmetry_t s37 {0x37, piL, piR}; // 3:00->11 7:01->11
            s37.remapping_equations( pX, &tTt, &tMk );
            res.push_back( s37 );
            symmetry_t s4C {0x4C, piL, piR}; // 4:01->00 C:11->00
            s4C.remapping_equations( pX, &tTt, &tMk );
            res.push_back( s4C );
          }
          if( eq02 && eq23 && eq03 )
          {
            symmetry_t s8C {0x8C, piL, piR}; // 8:10->00 C:11->00
            s8C.remapping_equations( pX, &tTt, &tMk );
            res.push_back( s8C );
            symmetry_t s3B {0x3B, piL, piR}; // 3:00->11 B:10->11
            s3B.remapping_equations( pX, &tTt, &tMk );
            res.push_back( s3B );
          }
          if( eq13 && eq23 && eq12 )
          {
            symmetry_t s6E {0x6E, piL, piR}; // 6:01->10 E:11->10
            s6E.remapping_equations( pX, &tTt, &tMk );
            res.push_back( s6E );
            symmetry_t s9D {0x9D, piL, piR}; // 9:10->01 D:11->01
            s9D.remapping_equations( pX, &tTt, &tMk );
            res.push_back( s9D );
          }
        }
      }
    }  
    return res;
  }


  decomposition_t decomposition_analysis( std::vector<TT> * pX )
  {
    decomposition_t res;
    int pi;

    uint32_t levelWorst{0};
    std::vector<uint32_t> vIdsNodesCritical;
    for( uint32_t iNd{0}; iNd < nNodes-1; ++iNd )
    {
      if( nodes[iNd].level > levelWorst )
      {
        levelWorst = nodes[iNd].level;
        vIdsNodesCritical = {iNd};
      }
      else if( nodes[iNd].level == levelWorst )
        vIdsNodesCritical.push_back(iNd);
    }

    for( auto iNd : vIdsNodesCritical )
    {
      pi = nodes[iNd].idPi;
      TT tt0  = kitty::cofactor0( tTt, pi );
      TT tt1  = kitty::cofactor1( tTt, pi );
      TT mk0  = kitty::cofactor0( tMk, pi );
      TT mk1  = kitty::cofactor1( tMk, pi );

      if( kitty::is_const0( tt0 & mk0 ) ) // f0 = 0
      {
        res.pi = pi;
        res.tMk = mk1;
        res.tTt = tt1;
        res.type = 0x8;
        return res;
      }
      else if( kitty::is_const0( tt1 & mk1 ) ) // f1 = 0
      {
        res.pi = pi;
        res.tMk = mk0;
        res.tTt = tt0;
        res.type = 0x2;
        return res;
      }
      else if( kitty::equal( tt0 & mk0, mk0 ) ) // f0 = 1
      {
        res.pi = pi;
        res.tMk = mk1;
        res.tTt = tt1;
        res.type = 0xB;
        return res;
      } 
      else if( kitty::equal( tt1 & mk1, mk1 ) ) // f1 = 1
      {
        res.pi = pi;
        res.tMk = mk0;
        res.tTt = tt0;
        res.type = 0xE;
        return res;
      } 
      else if( kitty::equal( tt1 & mk1 & mk0, ~tt0 & mk1 & mk0 ) ) // f1 = f0'
      {
        res.pi = pi;
        res.tMk = mk0 | mk1;
        res.tTt = ( mk0 & tt0 ) | ( mk1 & ~tt1 );
        res.type = 0x6;
        return res;
      }  
    }  
    res.type = 0x0;
    return res;
  }



  void print()
  {
    for( int j{0}; j < nodes.size(); ++j )
    {
        node_t node = nodes[j];
        uint32_t x  = node.get_this_ref_id();
        uint32_t xL = node.get_linp_ref_id();
        uint32_t xR = node.get_rinp_ref_id();
        uint32_t c  = node.get_this_cut_id();
        uint32_t cL = node.get_linp_cut_id();
        uint32_t cR = node.get_rinp_cut_id();
        std::string sLevel = node.level == UNK32 ? "?" : std::to_string(node.level);
        switch ( node.gate )
        {
          case gate_t::PIS_ : { printf("[ PI %d.%2d @ %s ]", c, x, sLevel.c_str() ); break; }
          case gate_t::CNTR : { printf("[00 %d @ %s ]", x, sLevel.c_str() ); break; }
          case gate_t::PA00 : { printf("[%d.%d=and( %d.%2d', %d.%2d' ) @ %s ]", c, x, cL, xL, cR, xR, sLevel.c_str()  ); break; }
          case gate_t::PA01 : { printf("[%d.%d=and( %d.%2d', %d.%2d  ) @ %s ]", c, x, cL, xL, cR, xR, sLevel.c_str()  ); break; }
          case gate_t::CMPL : { printf("[%d.%d=not(    %d.%2d     ) @ %s ]", c, x, cL, xL, sLevel.c_str() ); break; }
          case gate_t::PA10 : { printf("[%d.%d=and( %d.%2d , %d.%2d' ) @ %s ]", c, x, cL, xL, cR, xR, sLevel.c_str() ); break; }
          case gate_t::CMPR : { printf("[%d.%d=not(    %d.%2d     ) @ %s ]", c, x, cR, xR, sLevel.c_str() ); break; }
          case gate_t::EXOR : { printf("[%d.%d=xor( %d.%2d , %d.%2d  ) @ %s ]", c, x, cL, xL, cR, xR, sLevel.c_str() ); break; }
          case gate_t::IA11 : { printf("[%d.%d= or( %d.%2d', %d.%2d' ) @ %s ]", c, x, cL, xL, cR, xR, sLevel.c_str() ); break; }
          case gate_t::PA11 : { printf("[%d.%d=and( %d.%2d , %d.%2d  ) @ %s ]", c, x, cL, xL, cR, xR, sLevel.c_str() ); break; }
          case gate_t::XNOR : { printf("[%d.%d=xor( %d.%2d', %d.%2d' ) @ %s ]", c, x, cL, xL, cR, xR, sLevel.c_str() ); break; }
          case gate_t::PRJR : { printf("[%d.%d=buf(    %d.%2d     ) @ %s ]", c, x, cR, xR, sLevel.c_str() ); break; }
          case gate_t::IA10 : { printf("[%d.%d= or( %d.%2d', %d.%2d  ) @ %s ]", c, x, cL, xL, cR, xR, sLevel.c_str() ); break; }
          case gate_t::PRJL : { printf("[%d.%d=buf(    %d.%2d     ) @ %s ]", c, x, cL, xL, sLevel.c_str() ); break; }
          case gate_t::IA01 : { printf("[%d.%d= or( %d.%2d , %d.%2d' ) @ %s ]", c, x, cL, xL, cR, xR, sLevel.c_str() ); break; }
          case gate_t::IA00 : { printf("[%d.%d= or( %d.%2d , %d.%2d  ) @ %s ]", c, x, cL, xL, cR, xR, sLevel.c_str() ); break; }
          case gate_t::TAUT : { printf("[11 %d.%2d @ %s ]", c, x, sLevel.c_str() ); break; }
          case gate_t::TARG : { printf("[ PO %d.%2d @ %s]", c,  x, sLevel.c_str() ); break; }
          default:  break;
        }
    }
    printf("\n");
  }
};


template<class Ntk>
struct net_t
{
public:
  std::vector<cut_t> cuts;
  uint32_t nCuts{0u};
  std::vector<uint32_t> vCutsEdge;
  uint32_t idCutPo{0};
  uint32_t idRefPo{0};
  bool error{false};
  std::vector<TT> X;

public:
  net_t( TT const&, std::vector<uint32_t> const& );
  ~net_t();

  void add_remapping_cut( uint32_t, cut_t );
  void add_decomposition_cut( uint32_t, decomposition_t );
  cut_t candidate_cut_from_symmetry( cut_t *, symmetry_t * );
  void print();
};

#pragma region NET_t

template<class Ntk>
net_t<Ntk>::net_t( TT const& func, std::vector<uint32_t> const& levels )
{
  uint32_t nVars = func.num_vars();
  TT mk = ~func.construct();

  /* output cut */
  cut_t ocut(nCuts++);
  uint32_t iOut = ocut.add_node( gate_t::TARG, UNK32, func, mk );
  cuts.push_back( ocut );
  idCutPo = ocut.id;
  idRefPo = iOut;

  /* input cut */
  cut_t icut(nCuts++);
  uint32_t iNd;
  for( uint32_t i{0}; i<nVars; ++i )
  {
    X.emplace_back(nVars);
    kitty::create_nth_var( X[i], i );
    iNd = icut.add_node( gate_t::PIS_, levels[i], X[i], mk );
    icut.nodes[iNd].idPi = i;
  }
  icut.set_target( func, mk, ocut.id, iOut );
  icut.fill_pi_to_node();

  cuts.push_back( icut );
  vCutsEdge.push_back( icut.id );
}

template<class Ntk>
net_t<Ntk>::~net_t()
{
}

template<class Ntk>
void net_t<Ntk>::add_remapping_cut( uint32_t idxCutEdge, cut_t cutRemap )
{
  cutRemap.update_cut_id( nCuts++ );
  cuts.push_back( cutRemap );
  vCutsEdge.erase( vCutsEdge.begin() + idxCutEdge );
  vCutsEdge.push_back( cutRemap.id );
}

template<class Ntk>
void net_t<Ntk>::add_decomposition_cut( uint32_t idxCutEdge, decomposition_t dec )
{
  /* create a cut containing the remainder*/
  uint32_t idCutEdge = vCutsEdge[idxCutEdge];
  cut_t * pCutPrev = &cuts[idCutEdge];
  node_t * pNdDiv = &(pCutPrev->nodes[pCutPrev->pi_to_node[dec.pi]]);
  node_t * pTrgPrev = &cuts[pCutPrev->tCut].nodes[pCutPrev->tRef];

  cut_t tCut(nCuts++);

  uint32_t iOut;
  TT sTtDiv = pCutPrev->nodes[pCutPrev->pi_to_node[dec.pi]].sTt;
  switch (dec.type)
  {
  case 0x8: // top and
    iOut = tCut.add_node( gate_t::TARG, UNK32, pTrgPrev->sTt, pTrgPrev->sMk & ~sTtDiv );
    pTrgPrev->gate = PA11;
    pTrgPrev->level = pNdDiv->level + 1;
    break;
  case 0xE: // top or
    iOut = tCut.add_node( gate_t::TARG, UNK32, pTrgPrev->sTt, pTrgPrev->sMk & sTtDiv );
    pTrgPrev->gate = IA00;
    pTrgPrev->level = pNdDiv->level + 1;
    break;
  case 0xB: // top le
    iOut = tCut.add_node( gate_t::TARG, UNK32, pTrgPrev->sTt, pTrgPrev->sMk & ~sTtDiv );
    pTrgPrev->gate = IA10;
    pTrgPrev->level = pNdDiv->level + 1;
    break;
  case 0x2:
    iOut = tCut.add_node( gate_t::TARG, UNK32, pTrgPrev->sTt, pTrgPrev->sMk & sTtDiv );
    pTrgPrev->gate = PA01;
    pTrgPrev->level = pNdDiv->level + 1;
    break;
  case 0x6:
    iOut = tCut.add_node( gate_t::TARG, UNK32, sTtDiv ^ pTrgPrev->sTt , pTrgPrev->sMk );
    pTrgPrev->gate = EXOR;
    pTrgPrev->level = pNdDiv->level + delay_xor<Ntk>();
    break;
  default:
    break;
  }
  cuts.push_back( tCut );

  cut_t rCut = *pCutPrev;
  rCut.update_cut_id( nCuts++ );
  rCut.erase_node_from_pi( dec.pi );
  rCut.tCut = tCut.id;
  rCut.tRef = iOut;
  rCut.tMk = dec.tMk;
  rCut.tTt = dec.tTt;
  cuts.push_back( rCut );

  vCutsEdge.erase( vCutsEdge.begin() + idxCutEdge );
  vCutsEdge.push_back( rCut.id );
}


template<class Ntk>
cut_t net_t<Ntk>::candidate_cut_from_symmetry( cut_t * pCut, symmetry_t * pSym )
{
  cut_t newCut( nCuts );

  for( uint32_t i{0}; i < pCut->nodes.size(); ++i )
  {
    uint32_t idPi    = pCut->nodes[i].idPi;

    if( idPi == pSym->piL ) 
    {
      newCut.add_node_symL<Ntk>( pCut, pSym );
    }
    else if( idPi == pSym->piR ) 
    {
      newCut.add_node_symR<Ntk>( pCut, pSym );
    }
    else
    {
      newCut.add_node( gate_t::PRJL, pCut->nodes[i].level, pCut->nodes[i].sTt, pCut->nodes[i].sMk, pCut->id, pCut->nodes[i].idPi, pCut->nodes[i].idL, pCut->nodes[i].idR );
    }
  }
  newCut.tTt = pSym->tTt;
  newCut.tMk = pSym->tMk;
  newCut.tCut = pCut->tCut;
  newCut.tRef = pCut->tRef;

  return newCut;
}

template<class Ntk>
void net_t<Ntk>::print()
{
  for( int i{0}; i < cuts.size(); ++i )
  {
    printf(" CUT %d\n", i );
    cuts[i].print();
  }
  if( vCutsEdge.size() > 0 )
  {
    printf("active cuts:\n");
    for( auto c : vCutsEdge )
      printf( "%d->[%d %d] ", c, cuts[c].tCut, cuts[c].tRef );
  }
  printf("\n");
}


#pragma endregion SYMMETRIES_t


template<class Ntk>
class sym_synthesis
{
public:
    net_t<Ntk> net;

public:
    sym_synthesis( TT const&, std::vector<uint32_t> const& );
    ~sym_synthesis();

    bool try_functionality_matching();
    bool try_top_decomposition_on_critical();
    bool try_symmetry_remapping();
    void run();
};

#pragma region SYNTHESIS_t

template<class Ntk>
sym_synthesis<Ntk>::sym_synthesis( TT const& func, std::vector<uint32_t> const& levels ) : net(func, levels)
{
  run();
}

template<class Ntk>
sym_synthesis<Ntk>::~sym_synthesis()
{
}

template<class Ntk>
bool sym_synthesis<Ntk>::try_functionality_matching()
{
  for( uint32_t idxCutEdge{0}; idxCutEdge<net.vCutsEdge.size(); ++idxCutEdge )
  {
    uint32_t idCutEdge = net.vCutsEdge[idxCutEdge];
    cut_t * pCut = &net.cuts[idCutEdge];
    node_t * pTrg = &(net.cuts[pCut->tCut].nodes[pCut->tRef]);
        
    for( uint32_t iNd{0}; iNd<pCut->nodes.size(); ++iNd )
    {
      if( kitty::equal( pCut->nodes[iNd].sTt & pTrg->sMk, pTrg->sTt & pTrg->sMk ) ) // equal matching
      {
        pTrg->gate = gate_t::PRJL;
        pTrg->idL = pCut->nodes[iNd].id;
        pTrg->level = pCut->nodes[iNd].level;
        net.vCutsEdge.erase( net.vCutsEdge.begin() + idxCutEdge );
        return true;
      }
      else if ( kitty::equal( pCut->nodes[iNd].sTt & pTrg->sMk, ~pTrg->sTt & pTrg->sMk ) ) // equal matching
      {
        pTrg->gate = gate_t::CMPL;
        pTrg->idL = pCut->nodes[iNd].id;
        pTrg->level = pCut->nodes[iNd].level;
        net.vCutsEdge.erase( net.vCutsEdge.begin() + idxCutEdge );
        return true;
      }
    }
  }
  return false;
}

template<class Ntk>
bool sym_synthesis<Ntk>::try_symmetry_remapping()
{
  uint32_t idxCutEdge{net.vCutsEdge.size()-1};
  uint32_t idCutEdge{net.vCutsEdge.back()};
  /* perform symmetry analysis */
  cut_t * pCut = &net.cuts[idCutEdge];
  std::vector<symmetry_t> candidates = pCut->symmetry_analysis( &net.X );

  if( candidates.size() == 0 )  return false;

  /* select the symmetry */
  cut_t cutBest;
  uint32_t delayBest {0xFFFFFFFF};
  uint32_t rewardBest {0x00000000};

  for( int iCand{0}; iCand<candidates.size(); ++iCand )
  {
    cut_t candidateCut = net.candidate_cut_from_symmetry( pCut, &candidates[iCand] );
    if( candidates[iCand].reward > rewardBest || ( candidates[iCand].reward == rewardBest && candidateCut.delayCost < delayBest ) )
    {
      rewardBest = candidates[iCand].reward;
      delayBest = candidateCut.delayCost;
      cutBest = candidateCut;
    }
  }
  net.add_remapping_cut( idxCutEdge, cutBest );

  return true;
}

template<class Ntk>
bool sym_synthesis<Ntk>::try_top_decomposition_on_critical()
{
  uint32_t idxCutEdge{net.vCutsEdge.size()-1};
  uint32_t idCutEdge{net.vCutsEdge.back()};
  cut_t * pCut = &net.cuts[idCutEdge];
  decomposition_t dec = pCut->decomposition_analysis( &net.X );
  if( dec.type == 0x0 )  return false;
  net.add_decomposition_cut( idxCutEdge, dec );
  return true;
}

template<class Ntk>
void sym_synthesis<Ntk>::run()
{
  while( net.vCutsEdge.size() > 0 && !net.error )
  {
    if( try_functionality_matching() ) continue;
    //if( try_top_decomposition_on_critical() ) continue;
    if( try_symmetry_remapping() )  continue;
    net.error = true;
    break;
  }
  if( net.error )
    printf("ERROR: network not found\n");
  net.print();
}

#pragma endregion SYNTHESIS_t

} // namespace techaware

} // namespace mockturtle