/* mocktufuncle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
 *
 * Pemaskission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to pemaskit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this pemaskission notice shall be
 * included in all copies or substantial pofuncions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PAfuncICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TOfunc OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file DecSims.hpp
  \brief data structure for storing the Sims

  \author Andrea Costamagna
*/
#pragma once

#include <stdio.h>
#include <stack>
#include "DecNet.hpp"
#include <kitty/operators.hpp>
#include <kitty/print.hpp>


namespace mockturtle
{

enum class DecAct_t
{
  ERASE_,
  /* top decompositions */
  D1_AND_,
  D1_OR_,
  D1_XOR_,
  D1_LT_,
  D1_LE_,
  /* top decompositions with target */
  T_AND_,
  T_OR_,
  T_LT_,
  T_LE_,
  T_XOR_,
  T1_AND_,
  T1_OR_,
  T1_XOR_,
  T1_LT_,
  T1_LE_,
  /* closures */
  CL_AND_,
  CL_NAND_,
  CL_OR_,
  CL_NOR_,
  CL_XOR_,
  CL_XNOR_,
  CL_LT_,
  CL_LE_,
  XCL_AND_,
  XCL_NAND_,
  XCL_OR_,
  XCL_NOR_,
  XCL_XOR_,
  XCL_XNOR_,
  XCL_LT_,
  XCL_LE_,
  /* shannon decomposition */
  SD_,
  /* symmetry-remappings */
  NES_,
  ES_,
  MS_,
  SVS_,
  CSVS_,
  STR_,
  /* termination */
  BUF_,
  INV_
};

template<class TT>
struct action_t
{
  DecAct_t type;
  std::vector<int> sigs; // signals characterizing the move TD [x,T]
  std::vector<TT> mask;                 // setmaskapped mask
  std::vector<TT> func;                 // setmaskapped truth
  int reward;
  int cost;
  int id_ord;
  int id_sym;
};

template<class TT, class Ntk>
class DecAnalyzer
{
private:
  DecNet<TT, Ntk> * pNet;
  std::vector<signal_t> *pX;
  std::vector<signal_t> *pY;
  std::vector<std::vector<int>>      *pvS;
  std::vector<std::vector<signal_t>> *pvD;
  /* only inputs */
  std::vector<action_t<TT>> set_topdec; // top decs when only input patterns
  std::vector<action_t<TT>> set_remap;

  std::vector<action_t<TT>> set_remove; // nodes that can be removed 
  std::vector<action_t<TT>> set_closure;
  /* divisors */
  std::vector<action_t<TT>> set_trgclosure;
  std::vector<action_t<TT>> set_2stepsdec; //

  int cost_xor  = std::is_same<Ntk, xag_network>::value ? 1 : 3;//8;
  int cost_xnor = std::is_same<Ntk, xag_network>::value ? 1 : 3;//7;
  int cost_and  = 1;//1;
  int cost_lt   = 1;//2;
  int cost_le   = 1;//3;
  int cost_or   = 1;//4;

public:
  DecAnalyzer(DecNet<TT, Ntk> *, std::vector<signal_t> *, std::vector<signal_t> *, std::vector<std::vector<int>> *, std::vector<std::vector<signal_t>> * );
  ~DecAnalyzer();
  /* remapping */
  action_t<TT> simple_remapping( int, int, int, int, int, DecAct_t, int );
  action_t<TT> multiform_remapping( int, int, int, int, DecAct_t );
  action_t<TT> compatible_remapping( int, int, int, int, int, int, DecAct_t, int );
  action_t<TT> spectral_translation( int, int, int, DecAct_t, int );

  /* analyze */
  void check_tarclosure();
  void check_divclosure();
  void check_2stepsdec();
  void check1();
  void check2();
  void check_tdecore();
  void check_str();
  void check_topdec( std::vector<signal_t> );

  bool check_sat( std::vector<signal_t>, int );
  bool check_satX( std::vector<signal_t>, int );

  std::vector<action_t<TT>> get_topdec();
  std::vector<action_t<TT>> get_remove();
  std::vector<action_t<TT>> get_remap();
  std::vector<action_t<TT>> get_divclosure();
  std::vector<action_t<TT>> get_closure();
  std::vector<action_t<TT>> get_trgclosure();
  std::vector<action_t<TT>> get_2stepsdec();

  int CountNodes( DecAct_t );

  TT cube_generator( int, TT, TT );

  void print_actions(std::vector<action_t<TT>>);
};

#pragma region propefuncies
template<class TT, class Ntk> 
DecAnalyzer<TT, Ntk>::DecAnalyzer( DecNet<TT, Ntk> * pNet, std::vector<signal_t> (*pX), std::vector<signal_t> ((*pY)), std::vector<std::vector<int>> (*pvS), std::vector<std::vector<signal_t>> (*pvD)):
pNet(pNet),
pX(pX),
pY(pY),
pvS(pvS),
pvD(pvD)
{}
template<class TT, class Ntk> DecAnalyzer<TT, Ntk>::~DecAnalyzer(){}

template<class TT, class Ntk> std::vector<action_t<TT>> DecAnalyzer<TT, Ntk>::get_topdec(){ return set_topdec; }
template<class TT, class Ntk> std::vector<action_t<TT>> DecAnalyzer<TT, Ntk>::get_remove(){ return set_remove; }
template<class TT, class Ntk> std::vector<action_t<TT>> DecAnalyzer<TT, Ntk>::get_remap(){ return set_remap; }
template<class TT, class Ntk> std::vector<action_t<TT>> DecAnalyzer<TT, Ntk>::get_divclosure(){ return set_closure; }
template<class TT, class Ntk> std::vector<action_t<TT>> DecAnalyzer<TT, Ntk>::get_closure(){ return set_closure; }
template<class TT, class Ntk> std::vector<action_t<TT>> DecAnalyzer<TT, Ntk>::get_trgclosure(){ return set_trgclosure; }
template<class TT, class Ntk> std::vector<action_t<TT>> DecAnalyzer<TT, Ntk>::get_2stepsdec(){ return set_2stepsdec; }
#pragma endregion

#pragma region termination
template<class TT, class Ntk>
void DecAnalyzer<TT, Ntk>::check_divclosure()
{
  for( int iTrg{0}; iTrg < pY->size(); ++iTrg )
  {
    for( int iS{0}; iS < (*pvS)[iTrg].size(); ++iS )
    {
      int iDiv = (*pvS)[iTrg][iS];
      TT * pFT = pNet->getFuncP( (*pY)[iTrg] );
      TT * pFM = pNet->getMaskP( (*pY)[iTrg] );
      TT * pDT = pNet->getFuncP( (*pX)[iDiv] );
      TT * pDM = pNet->getMaskP( (*pX)[iDiv] );

      bool isFM_C_DM = kitty::is_const0( (*pFM)&(~*pDM) );
      
      if(isFM_C_DM)
      {
        if( kitty::equal( *pFT & *pFM, *pDT & *pFM ) )
        {
          action_t<TT> act = { DecAct_t::BUF_, {iTrg, iS}, {*pFM}, {*pFT}, kitty::count_zeros( *pFM & ~(*pFM) ) };
          set_closure.push_back( act );
        }
        else if( kitty::equal( (~*pFT) & *pFM, *pDT & *pFM ) )
        {
          action_t<TT> act = { DecAct_t::INV_, {iTrg, iS}, {*pFM}, {*pFT}, kitty::count_zeros( *pFM & ~(*pFM) ) };
          set_closure.push_back( act );
        }
      }
    }
  }
}

template<class TT, class Ntk>
void DecAnalyzer<TT, Ntk>::check_tarclosure()
{
  if( pY->size() > 1 )
  {
    for( int iTrg2{1}; iTrg2 < pY->size(); ++iTrg2 )
    {
      for( int iTrg1{0}; iTrg1 < iTrg2; ++iTrg1 )
      {
        TT * pF1T = pNet->getFuncP( ((*pY))[iTrg1] );
        TT * pF1M = pNet->getMaskP( ((*pY))[iTrg1] );
        TT * pF2T = pNet->getFuncP( ((*pY))[iTrg2] );
        TT * pF2M = pNet->getMaskP( ((*pY))[iTrg2] );

        TT   s1C2 = ~(*pF2M) | *pF1M; // M1 >= M2, M1 C M2 M1 contained in M1 : can only F2 x--> F1
        TT   s2C1 = ~(*pF1M) | *pF2M; // M2 >= M1, M2 C M1 M2 contained in M2 : can only F1 x--> F2

        if( kitty::count_zeros( s1C2 ) == 0 ) // masks are implied M1 >= M2
        {
          if( kitty::equal( *pF1T & *pF2M, *pF2T & *pF2M ) )
          {
            action_t<TT> act = { DecAct_t::BUF_, {iTrg1, iTrg2}, {*pF1M}, {*pF1T}, kitty::count_zeros( *pF1M & ~(*pF1M) ) };
            set_trgclosure.push_back( act );
          }
          else if( kitty::equal( (~*pF1T) & *pF2M, *pF2T & *pF2M ) )
          {
            action_t<TT> act = { DecAct_t::INV_, {iTrg1, iTrg2}, {*pF1M}, {*pF1T}, kitty::count_zeros( *pF1M & ~(*pF1M) ) };
            set_trgclosure.push_back( act );
          }
        }
        if( kitty::count_zeros( s2C1 ) == 0 ) // masks are implied M2 >= M1
        {
          if( kitty::equal( *pF1T & *pF1M, *pF2T & *pF1M ) )
          {
            action_t<TT> act = { DecAct_t::BUF_, {iTrg2, iTrg1}, {*pF2M}, {*pF2T}, kitty::count_zeros( *pF2M & ~(*pF2M) ) };
            set_trgclosure.push_back( act );
          }
          else if( kitty::equal( (~*pF1T) & *pF1M, *pF2T & *pF1M ) )
          {
            action_t<TT> act = { DecAct_t::INV_, {iTrg2, iTrg1}, {*pF2M}, {*pF2T}, kitty::count_zeros( *pF2M & ~(*pF2M) ) };
            set_trgclosure.push_back( act );
          }
        }
      }
    }
  }
}
#pragma endregion

template<class TT, class Ntk> 
int DecAnalyzer<TT, Ntk>::CountNodes( DecAct_t type)
{
  switch (type)
  {
  case DecAct_t::NES_ :
    return 2;
    break;
  case DecAct_t::ES_ :
    return 2;
    break;
  case DecAct_t::SVS_:
    return 1;
    break;
  case DecAct_t::CSVS_:
    return 1;
    break;
  case DecAct_t::MS_:
    return 1;
    break;
  default:
    return -1;
    break;
  }
}

#pragma region decomposability
template<class TT, class Ntk>
void DecAnalyzer<TT, Ntk>::check1()
{
  for( int iTrg{0}; iTrg < pY->size(); ++iTrg )
  {
    for( int iS{0}; iS < (*pvS)[iTrg].size(); ++iS )
    {
      int iDiv = (*pvS)[iTrg][iS];
      TT * pFT = pNet->getFuncP( ((*pY))[iTrg] );
      TT * pFM = pNet->getMaskP( ((*pY))[iTrg] );
      TT * pDT = pNet->getFuncP( (*pX)[iDiv] );
      TT * pDM = pNet->getMaskP( (*pX)[iDiv] );

      TT tt0  = kitty::cofactor0( * pFT, (*pvS)[iS] );
      TT tt1  = kitty::cofactor1( * pFT, (*pvS)[iS] );
      TT mk0  = kitty::cofactor0( * pFM, (*pvS)[iS] );
      TT mk1  = kitty::cofactor1( * pFM, (*pvS)[iS] );

      bool eq0_0 = is_const0( *pFM & ~*pDT & *pFT );
      bool eq0_1 = is_const0( *pFM & ~*pDT & ~*pFT );
      bool eq1_0 = is_const0( *pFM & *pDT & *pFT );
      bool eq1_1 = is_const0( *pFM & *pDT & ~*pFT );
      bool eq10_ = equal( mk0 & mk1 & ~tt1 , mk0 & mk1 & tt0 );
      bool mask = equal( mk0 & mk1 & tt1 , mk0 & mk1 & tt0 );

      if( mask )
      {
        action_t<TT> act = { DecAct_t::ERASE_, {iTrg, iS}, {*pFM}, {*pFT}, kitty::count_ones(~*pFM ) };
        set_remove.push_back( act );
      }
      else
      {
        if ( eq0_0 ) // F0 = 0  =>  F = Xi & F1
        {
          action_t<TT> act = { DecAct_t::D1_AND_, {iTrg, iS}, {*pDT & *pFM}, {*pFT}, kitty::count_ones(~(*pDT & *pFM)) };
          set_topdec.push_back( act );
        }
        if ( eq0_1 ) // F0 = 1  =>  F = Xi' + F1
        {
          action_t<TT> act = { DecAct_t::D1_LE_, {iTrg, iS}, {*pDT & *pFM}, {*pFT}, kitty::count_ones(~(*pDT & *pFM)) };
          set_topdec.push_back( act );
        }
        if ( eq1_0 ) // F1 = 0  =>  F = Xi' & F0
        {
          action_t<TT> act = { DecAct_t::D1_LT_, {iTrg, iS}, {~*pDT & *pFM}, {*pFT}, kitty::count_ones(~(~*pDT & *pFM)) };
          set_topdec.push_back( act );
        }
        if ( eq1_1 ) // F1 = 1  =>  F = Xi + F1
        {
          action_t<TT> act = { DecAct_t::D1_OR_, {iTrg, iS}, {~*pDT & *pFM}, {*pFT}, kitty::count_ones(~(~*pDT & *pFM)) };
          set_topdec.push_back( act );
        }
        if ( eq10_ )// F1 = F0' =>  F = Xi ^ F0
        {
          action_t<TT> act = { DecAct_t::D1_XOR_, {iTrg, iS}, {*pFM}, {( *pFT & *pFM & ~*pDT )|( ~*pFT & *pFM & *pDT )}, kitty::count_ones(*pFM) };
          set_topdec.push_back( act );      
        }
      }
    }
  }
}

template<class TT, class Ntk>
void DecAnalyzer<TT, Ntk>::check_2stepsdec()
{
  assert( (*pvS).size() == (*pY).size() ); // require that the set of signals is given 

  for( int iTrg{0}; iTrg < pY->size(); ++iTrg )
  {
    for( int iS{0}; iS < (*pvS)[iTrg].size(); ++iS )
    {
      int iDiv = (*pvS)[iTrg][iS];
      TT * pFT = pNet->getFuncP( (*pY)[iTrg] );
      TT * pFM = pNet->getMaskP( (*pY)[iTrg] );
      TT * pDT = pNet->getFuncP( (*pvD)[iTrg][iDiv] );
      TT * pDM = pNet->getMaskP( (*pvD)[iTrg][iDiv] );

      bool mask = kitty::is_const0((*pDM)&(*pFM)&(*pDT)) || 
                  kitty::is_const0((*pDM)&(*pFM)&(~*pFT));
      if( mask )
      {
        action_t<TT> act = { DecAct_t::ERASE_, {iTrg, iS}, {*pFM}, {*pFT}, kitty::count_ones(~*pFM ) };
        set_remove.push_back( act );
      }
      else
      {
        TT T1 =  *pDT & *pFT;
        TT T0 = ~*pDT & *pFT;
        TT M1 =  *pDT & *pDM & *pFM;
        TT M0 = ~*pDT & *pDM & *pFM;
        int reward = std::max( kitty::count_ones( (*pDT ^ *pFT)&(*pFM) ), kitty::count_ones(((~*pDT) ^ *pFT)&(*pFM)) );
          action_t<TT> act = { DecAct_t::SD_, {iTrg, iS}, 
                              {M0, M1}, {T0, T1}, 
                              reward };
          set_2stepsdec.push_back( act );
      }
    }
  }
}
#pragma endregion

#pragma region remappinng

/*
   * Compute the general cofactor with respect to the cube G
   *            ji 
   *   G = 0 -> 00 : F( Xi=0, Xj=0 )
   *   G = 1 -> 01 : F( Xi=0, Xj=1 )
   *   G = 2 -> 10 : F( Xi=1, Xj=0 )
   *   G = 3 -> 11 : F( Xi=1, Xj=1 )
   * NB: notice the swap of (i,j), assuming that i < j in the given label
  */
 template<class TT>
  TT cofactorG( TT  fn, uint32_t G, uint32_t i, uint32_t j )
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
template<class TT, class Ntk>
  TT DecAnalyzer<TT, Ntk>::cube_generator( int cube, TT Xi, TT Xj )
  {
    switch ( cube )
    {
    case 0: 
      return ~Xj & ~Xi;
      break;
    case 1: 
      return ~Xj &  Xi;
      break;
    case 2: 
      return  Xj & ~Xi;
      break;
    case 3: 
      return  Xj & Xi;
      break;
    }
  }

template<class TT, class Ntk>
  action_t<TT> DecAnalyzer<TT, Ntk>::simple_remapping( int from, int to, int i, int j, int iTrg, DecAct_t type, int id_symmetry )
  {
    assert( i < j );
    action_t<TT> res;
    res.type = type;
    res.sigs = {};
    res.sigs.push_back(iTrg);
    res.sigs.push_back(i);
    res.sigs.push_back(j); 

    res.id_ord = 1u*( from > to );
    res.id_sym = id_symmetry;
    
    uint32_t xi = (*pvS)[iTrg][i];
    uint32_t xj = (*pvS)[iTrg][j];

    TT DTi = *pNet->getFuncP( (*pX)[xi] );
    TT DTj = *pNet->getFuncP( (*pX)[xj] );
    //TT FT  = *pNet->getFuncP( (*pY)[iTrg] );
    TT Ms  = *pNet->getMaskP( (*pY)[iTrg] );
    TT FT = pNet->getFuncOSY();
    TT FM = pNet->getMaskCDC();

    TT A = cube_generator( from, DTi, DTj );
    TT B = cube_generator( to, DTi, DTj );
    TT ttA = cofactorG( FT, from, xi, xj );
    TT ttB = cofactorG( FT, to, xi, xj );
    TT mkA = cofactorG( FM, from, xi, xj );
    TT mkB = cofactorG( FM, to, xi, xj );

    res.mask.push_back( ( FM & ~A ) | ( B & mkA ) );//( ~A & *pFM ) | ( B ));//& mkA ));  //
    res.reward = kitty::count_zeros( res.mask[iTrg]&Ms );

    TT TA = A & FT; /* remapping for A */
    TT TB = B & ( ( mkB & FT ) | ( mkA & ttA ) ); /* remapping for B */
    TT TR = ( ~A & ~B & FT );

    res.func.push_back( TA | TB | TR );//FT

    switch (type)
    {
    case DecAct_t::NES_:
        res.cost = cost_and + cost_or;
        break;
    case DecAct_t::ES_:
        res.cost = res.id_ord == 0 ? 2*cost_le : 2*cost_and;
    case DecAct_t::SVS_:
        if( res.id_sym == 0 )
          res.cost = res.id_ord == 0 ? 2*cost_le : 2*cost_lt;
        else if( res.id_sym == 1 )
          res.cost = res.id_ord == 0 ? cost_or : cost_lt;
        else if( res.id_sym == 2 )
          res.cost = res.id_ord == 0 ? cost_le : cost_and;
        else if( res.id_sym == 3 )
          res.cost = res.id_ord == 0 ? cost_or : cost_lt;
    default:
      break;
    }
    
    return res;
  }

template<class TT, class Ntk>
action_t<TT> DecAnalyzer<TT, Ntk>::multiform_remapping( int from1, int i, int j, int iTrg, DecAct_t type )
  {
    assert( i < j );
    action_t<TT> res;
    res.type = type;
    res.sigs.push_back( iTrg );
    res.sigs.push_back(i); 
    res.sigs.push_back(j); 
    res.id_ord = from1;

    uint32_t xi = (*pvS)[iTrg][i];
    uint32_t xj = (*pvS)[iTrg][j];
    
    TT DTi = *pNet->getFuncP( (*pX)[xi] );
    TT DTj = *pNet->getFuncP( (*pX)[xj] );
    //TT FT  = *pNet->getFuncP( (*pY)[iTrg] );
    TT Ms  = *pNet->getMaskP( (*pY)[iTrg] );
    TT FT = pNet->getFuncOSY();
    TT FM = pNet->getMaskCDC();

    TT A = cube_generator( from1, DTi, DTj );
    uint32_t to1 = uint32_t( 3 - from1 );
    TT C = cube_generator( to1, DTi, DTj );

    uint32_t from2 = ( from1 == 0u )*1u + ( from1 == 1u )*3u + ( from1 == 3u )*2u;

    TT B = cube_generator( from2, DTi, DTj  );
    uint32_t to2 = uint32_t( 3 - from2 );
    TT D = cube_generator( to2, DTi, DTj  );

    TT ttA = cofactorG( FT, from1, xi, xj );
    TT ttB = cofactorG( FT, from2, xi, xj );
    TT ttC = cofactorG( FT, to1, xi, xj );
    TT ttD = cofactorG( FT, to2, xi, xj );

    TT mkA = cofactorG( FM, from1, xi, xj );
    TT mkB = cofactorG( FM, from2, xi, xj );
    TT mkC = cofactorG( FM, to1, xi, xj );
    TT mkD = cofactorG( FM, to2, xi, xj );

    res.mask.push_back(( ~B & ~A & FM ) | ( ( C & mkA ) | ( D & mkB ) ));
    res.reward = kitty::count_zeros( res.mask[0]&Ms ); // @ltry
    
    TT preserved = (~A&~B&~C&~D)&(FT);

    TT modifiedA = A & FT;
    TT modifiedB = B & FT;

    TT modifiedC = C & ( ( mkA & ~mkC & ttA ) | ( mkC & FT ) );
    TT modifiedD = D & ( ( mkB & ~mkD & ttB ) | ( mkD & FT ) );

    res.func.push_back( preserved | modifiedA | modifiedC | modifiedB | modifiedD);//
    
    if( res.id_ord == 0 || res.id_ord == 2 )
      res.cost = cost_xor;
    else
      res.cost = cost_xnor;

    return res;
  }

template<class TT, class Ntk>
action_t<TT> DecAnalyzer<TT, Ntk>::compatible_remapping( int from1, int from2, int to, int i, int j, int iTrg, DecAct_t type, int id_symmetry )
  {
    assert( i < j );
    action_t<TT> res;
    res.type = type;
    res.sigs.push_back(iTrg);
    res.sigs.push_back( i );
    res.sigs.push_back( j );
    res.id_ord = 1u*( from1 > from2 );
    res.id_sym = id_symmetry;

    uint32_t xi = (*pvS)[iTrg][i];
    uint32_t xj = (*pvS)[iTrg][j];

    TT DTi = *pNet->getFuncP( (*pX)[xi] );
    TT DTj = *pNet->getFuncP( (*pX)[xj] );
    //TT FT  = *pNet->getFuncP( (*pY)[iTrg] );
    TT Ms = *pNet->getMaskP( (*pY)[iTrg] );
    TT FT = pNet->getFuncOSY();
    TT FM = pNet->getMaskCDC();

    uint32_t excluded = uint32_t( 6u - from1 - from2 - to );

    TT A = cube_generator( from1, DTi, DTj );
    TT B = cube_generator( from2, DTi, DTj );
    TT C = cube_generator( to, DTi, DTj );
    TT D = cube_generator( excluded, DTi, DTj );

    TT ttA = cofactorG( FT, from1, xi, xj );
    TT ttB = cofactorG( FT, from2, xi, xj );
    TT ttC = cofactorG( FT, to, xi, xj );

    TT mkA = cofactorG( FM, from1, xi, xj );
    TT mkB = cofactorG( FM, from2, xi, xj );
    TT mkC = cofactorG( FM, to, xi, xj );

    res.mask.push_back( ( ~B & ~A & FM ) | ( C & ( mkA | mkB ) ) );
    
    res.reward = kitty::count_zeros( res.mask[iTrg]&Ms ); //@ltry

    TT TA = A & FT;
    TT TB = B & FT;
    TT TC = C & ( ( mkA & ttA ) | ( mkB & ttB ) | ( mkC & FT ) );
    TT TR = ~A & ~B & ~C & FT;

    res.func.push_back( TA | TB | TC | TR );//FT );/
    
    if( res.id_sym == 0 )
      res.cost = cost_and;
    else if( res.id_sym == 1 )
      res.cost = res.id_ord == 0 ? cost_le : cost_lt ;
    else if( res.id_sym == 2 )
      res.cost = res.id_ord == 0 ? cost_le : cost_lt ;
    else if( res.id_sym == 3 )
      res.cost = cost_or;

    return res;
  }

template<class TT, class Ntk>
void DecAnalyzer<TT, Ntk>::check2()
{
  for( int iTrg{0}; iTrg<pY->size(); ++iTrg )
  {
    for ( auto j = 1u; j < (*pvS)[iTrg].size(); ++j )
    {
      for ( auto i = 0u; i < j; ++i )
      {
        uint32_t xi = (*pvS)[iTrg][i];
        uint32_t xj = (*pvS)[iTrg][j];

        //TT FT = * pNet->getFuncP( (*pY)[iTrg] );
        //TT FM = * pNet->getMaskP( (*pY)[iTrg] );
        TT FT = pNet->getFuncOSY();
        TT FM = pNet->getMaskCDC();

        TT tti0  = kitty::cofactor0( FT, xi );
        TT tti1  = kitty::cofactor1( FT, xi );
        TT mki0  = kitty::cofactor0( FM, xi );
        TT mki1  = kitty::cofactor1( FM, xi );

        TT ttj0  = kitty::cofactor0( FT, xj );
        TT ttj1  = kitty::cofactor1( FT, xj );
        TT mkj0  = kitty::cofactor0( FM, xj );
        TT mkj1  = kitty::cofactor1( FM, xj );

        const auto tt00 = cofactor0( ttj0, xi );
        const auto tt01 = cofactor1( ttj0, xi );
        const auto tt10 = cofactor0( ttj1, xi );
        const auto tt11 = cofactor1( ttj1, xi );
        const auto mk00 = cofactor0( mkj0, xi );
        const auto mk01 = cofactor1( mkj0, xi );
        const auto mk10 = cofactor0( mkj1, xi );
        const auto mk11 = cofactor1( mkj1, xi );

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

        //if(  )
        if( (num_pairs != 0) )//&& !equal( mkj0 & mkj1 & ttj1 , mkj0 & mkj1 & ttj0 ) )   
        {     
          //if( !equal( mki0 & mki1 & tti1 , mki0 & mki1 & tti0 ) )
          {
            if ( eq12 ) // F01 = F10 NES
            {
              set_remap.push_back( simple_remapping( 1, 2, i, j, iTrg, DecAct_t::NES_, 0 ) );
              set_remap.push_back( simple_remapping( 2, 1, i, j, iTrg, DecAct_t::NES_, 0 ) );
            }
            if ( eq03 ) // F00 = F11 ES
            {
              set_remap.push_back( simple_remapping( 0, 3, i, j, iTrg, DecAct_t::ES_, 0 ) ) ;
              set_remap.push_back( simple_remapping( 3, 0, i, j, iTrg, DecAct_t::ES_, 0 ) ) ;
            }
            if ( eq02 ) // F00=F10
            {
              set_remap.push_back( simple_remapping( 0, 2, i, j, iTrg, DecAct_t::SVS_, 0u ) ); // SVS0X_
              set_remap.push_back( simple_remapping( 2, 0, i, j, iTrg, DecAct_t::SVS_, 0u ) );
            }
            if ( eq13 ) // F01=F11
            {
              set_remap.push_back(  simple_remapping( 1, 3, i, j, iTrg, DecAct_t::SVS_, 1 ) ); // SVS1X_
              set_remap.push_back(  simple_remapping( 3, 1, i, j, iTrg, DecAct_t::SVS_, 1 ) );
            }
            if ( eq01 ) // F01=F00
            {
              set_remap.push_back(  simple_remapping( 0, 1, i, j, iTrg, DecAct_t::SVS_, 2 ) ); // SVSX0
              set_remap.push_back(  simple_remapping( 1, 0, i, j, iTrg, DecAct_t::SVS_, 2 ) );    
            }
            if ( eq23 ) // F11=F10
            {
              set_remap.push_back(  simple_remapping( 2, 3, i, j, iTrg, DecAct_t::SVS_, 3 ) ); // SVSX1_
              set_remap.push_back(  simple_remapping( 3, 2, i, j, iTrg, DecAct_t::SVS_, 3 ) );
            }
            if ( eq12 && eq03 ) // F01=F10 and F00=F11
            {
              set_remap.push_back(  multiform_remapping( 0, i, j, iTrg, DecAct_t::MS_ ) );
              set_remap.push_back(  multiform_remapping( 1, i, j, iTrg, DecAct_t::MS_ ) );      
              set_remap.push_back(  multiform_remapping( 2, i, j, iTrg, DecAct_t::MS_ ) );      
              set_remap.push_back(  multiform_remapping( 3, i, j, iTrg, DecAct_t::MS_ ) );
            }
            if( eq02 && eq01 && eq12 )
            {
              set_remap.push_back(  compatible_remapping( 0, 1, 2, i, j, iTrg, DecAct_t::CSVS_, 0 ) ); // CSVS00_
              set_remap.push_back(  compatible_remapping( 2, 0, 1, i, j, iTrg, DecAct_t::CSVS_, 0 ) );
            }
            if( eq13 && eq01 && eq03 )
            {
              set_remap.push_back(  compatible_remapping( 0, 1, 3, i, j, iTrg, DecAct_t::CSVS_, 1 ) ); // CSVS10_ old ij, new ji
              set_remap.push_back(  compatible_remapping( 3, 1, 0, i, j, iTrg, DecAct_t::CSVS_, 1 ) );
            }
            if( eq02 && eq23 && eq03 )
            {
              set_remap.push_back(  compatible_remapping( 0, 2, 3, i, j, iTrg, DecAct_t::CSVS_, 2 ) ); // CSVS01_ old ij, new ji
              set_remap.push_back(  compatible_remapping( 3, 2, 0, i, j, iTrg, DecAct_t::CSVS_, 2 ) );
            }
            if( eq13 && eq23 && eq12 )
            {
              set_remap.push_back(  compatible_remapping( 1, 3, 2, i, j, iTrg, DecAct_t::CSVS_, 3 ) ); // CSVS11_
              set_remap.push_back(  compatible_remapping( 3, 2, 1, i, j, iTrg, DecAct_t::CSVS_, 3 ) );
            }
          }
        }
      }
    }
  }
}

template<class TT, class Ntk>
void DecAnalyzer<TT, Ntk>::check_tdecore()
{
  for( int iTrg{0}; iTrg<pY->size(); ++iTrg )
  {
    for ( auto i = 0u; i < (*pvS)[iTrg].size(); ++i )
    {
      uint32_t xi = (*pvS)[iTrg][i];

      TT FT = pNet->getFuncOSY();
      TT FM = pNet->getMaskCDC();

      TT tt0  = kitty::cofactor0( FT, xi );
      TT tt1  = kitty::cofactor1( FT, xi );
      TT mk0  = kitty::cofactor0( FM, xi );
      TT mk1  = kitty::cofactor1( FM, xi );

      if ( kitty::is_const0( tt0 & mk0 ) ) // F0=0 => AND
      {
        action_t<TT> act = { DecAct_t::T1_AND_, {iTrg, i}, { mk1 }, { tt1 }, kitty::count_ones( mk1 ) };
        set_remap.push_back( act );
      }
      if ( kitty::is_const0( tt1 & mk1 ) ) // F1=0 => LT
      {
        action_t<TT> act = { DecAct_t::T1_LT_, {iTrg, i}, { mk0 }, { tt0 }, kitty::count_ones( mk0 ) };
        set_remap.push_back( act );
      }
      if ( kitty::equal( tt0 & mk0, mk0 ) ) // F1=1 => x+R
      {
        action_t<TT> act = { DecAct_t::T1_OR_, {iTrg, i}, { ~mk1 }, { tt0 }, kitty::count_ones( ~mk1 ) };
        set_remap.push_back( act );
      }
      if ( kitty::equal( tt1 & mk1, mk1 ) ) // F1=1 => x+R
      {
        action_t<TT> act = { DecAct_t::T1_LE_, {iTrg, i}, { ~mk0 }, { tt1 }, kitty::count_ones( ~mk0 ) };
        set_remap.push_back( act );
      }
      if ( kitty::equal( tt1 & mk1 & mk0, ~tt0 & mk1 & mk0 ) ) // F1=1 => x+R
      {
        action_t<TT> act = { DecAct_t::T1_XOR_, {iTrg, i}, { mk0 | mk1 }, { (mk1&tt1) | (mk0&(~tt0)) }, kitty::count_ones( FM ) };
        set_remap.push_back( act );
      }

    }
  }
}


template<class TT, class Ntk>
  action_t<TT> DecAnalyzer<TT, Ntk>::spectral_translation( int i, int j, int iTrg, DecAct_t type, int id_symmetry )
  {
//printf("@@a\n");
    action_t<TT> res;
    res.type = type;
    res.sigs = {};
    res.sigs.push_back(iTrg);
    res.sigs.push_back(i);
    res.sigs.push_back(j); 
//printf("@@a2\n");
    uint32_t xi = (*pvS)[iTrg][i];
    uint32_t xj = (*pvS)[iTrg][j];
//printf("@@a3\n");
    //TT FT  = *pNet->getFuncP( (*pY)[iTrg] );
    TT Ms   = *pNet->getMaskP( (*pY)[iTrg] );
    TT FT   = pNet->getFuncOSY();
    TT FM   = pNet->getMaskCDC();
    TT DTi  = *pNet->getFuncP( (*pX)[xi] );
    TT DTj  = *pNet->getFuncP( (*pX)[xj] );
    TT tti0 = kitty::cofactor0( FT, xi );
    TT tti1 = kitty::cofactor1( FT, xi );
    TT ttj0 = kitty::cofactor0( FT, xj );
    TT ttj1 = kitty::cofactor1( FT, xj );
    TT mki0 = kitty::cofactor0( FM, xi );
    TT mki1 = kitty::cofactor1( FM, xi );
    TT mkj0 = kitty::cofactor0( FM, xj );
    TT mkj1 = kitty::cofactor1( FM, xj );
    res.reward = kitty::count_zeros( FM&Ms );
    res.func.push_back( (~DTj & FT) | (DTj&((DTi&~FT)|(~DTi&FT))) );//FT
    res.mask.push_back( (~DTj & FM) | (DTj&((DTi&~FM)|(~DTi&FM))) );//FT
    return res;
  }

template<class TT, class Ntk>
void DecAnalyzer<TT, Ntk>::check_str()
{
//printf("@@c\n");
  for( int iTrg{0}; iTrg<pY->size(); ++iTrg )
  {
    for ( auto j = 1u; j < (*pvS)[iTrg].size(); ++j )
    {
      for ( auto i = 0u; i < j; ++i )
      {
        set_remap.push_back( spectral_translation( i, j, iTrg, DecAct_t::STR_, 0 ) );
        set_remap.push_back( spectral_translation( j, i, iTrg, DecAct_t::STR_, 0 ) );
      }
    }
  }
//printf("@@d\n");

}


template<class TT, class Ntk>
void DecAnalyzer<TT, Ntk>::check_topdec( std::vector<signal_t> boolean_chain )
{
  assert( pY->size() == 1 );
  int r = boolean_chain.size();

  TT Tf = * pNet->getFuncP( (*pY)[0] );
  TT Mf = * pNet->getMaskP( (*pY)[0] ); 
  TT Mo = pNet->getMaskCDC(); 

  for( int i{r-1}; i>=0; --i )
  {
    signal_t xi = boolean_chain[i];
    TT Tdi = * pNet->getFuncP( xi );
    if( kitty::is_const0( (~Tdi) & Tf & Mf ) ) // F0 = 0 -> top and
    {
      TT M1 = Mf&Tdi;
      TT M2 = Mo&Mf&Tdi;
      action_t<TT> act = { DecAct_t::T1_AND_, {0, i }, { M1 }, { Tf }, kitty::count_zeros(M2), 1 };
      set_remap.push_back(act);
    }
    if( kitty::is_const0( Tdi & Tf & Mf ) ) // F1 = 0 -> top lt
    {
      TT M1 = Mf&(~Tdi);
      TT M2 = Mo&Mf&(~Tdi);
      action_t<TT> act = { DecAct_t::T1_LT_, {0, i }, { M1 }, {Tf}, kitty::count_zeros(M2), 1 };
      set_remap.push_back(act);
    }
    if( kitty::equal( Tdi & Mf, Tdi & Tf & Mf ) ) // F = 0 -> top or
    {      
      TT M1 = Mf&(~Tdi);
      TT M2 = Mo&Mf&(~Tdi);
      action_t<TT> act = { DecAct_t::T1_OR_, {0, i }, { M1 }, {Tf}, kitty::count_zeros(M2), 1 };
      set_remap.push_back(act);
    }    
    if( kitty::equal( (~Tdi) & Mf, (~Tdi) & Tf & Mf ) ) // F = 0 -> top or
    {
      TT M1 = Mf&(Tdi);
      TT M2 = Mo&Mf&(Tdi);
      action_t<TT> act = { DecAct_t::T1_LE_, {0, i }, { M1 }, {Tf}, kitty::count_zeros(M2), 1 };
      set_remap.push_back(act);
    }
  }
}


template<class TT, class Ntk>
bool DecAnalyzer<TT, Ntk>::check_sat( std::vector<signal_t> boolean_chain, int newBound )
{
  TT Tf = * pNet->getFuncP( (*pY)[0] );
  TT Mf = * pNet->getMaskP( (*pY)[0] );
  int r = boolean_chain.size();
  int nbits = kitty::count_zeros( Tf & ~Tf );

  for( int i{r-1}; i>=newBound; --i )
  {
    signal_t xi = boolean_chain[i];
    TT Tdi = * pNet->getFuncP( xi );
    if( kitty::equal( Mf&Tdi, Mf&Tf ) )
    {
      action_t<TT> act = { DecAct_t::BUF_, {0, i }, { Mf }, { Tf }, nbits, 0 };
      //printf("COST %d BEC BUF\n", act.cost);
      set_closure.push_back( act );
    }
    else if( kitty::equal( Mf&~Tdi, Mf&Tf ) )
    {
      action_t<TT> act = { DecAct_t::INV_, {0, i }, { Mf }, { Tf }, nbits, 0 };
      set_closure.push_back( act );
      //printf("COST %d BEC INV\n", act.cost);

    }

    // second loop: the previous variables
    if( i > 0 )
    {
      for( int j{i-1}; j>=0; --j )
      {
        signal_t xj = boolean_chain[j];
        TT Tdj = * pNet->getFuncP( xj );
        if( kitty::equal( Mf&( Tdi & Tdj ), Mf&Tf ) ) // AND
        {
          action_t<TT> act = { DecAct_t::CL_AND_, {0, i, j }, { Mf }, { Tf }, nbits, 1 };
          //printf("COST %d bec AND\n", act.cost );
          set_closure.push_back( act );
        }
        else if( kitty::equal( Mf&( ~(Tdi & Tdj) ), Mf&Tf ) ) // NAND
        {
          action_t<TT> act = { DecAct_t::CL_NAND_, {0, i, j }, { Mf }, { Tf }, nbits, 1 };
          //printf("COST %d bec NAND\n", act.cost );
          set_closure.push_back( act );
        }
        else if( kitty::equal( Mf&( Tdi | Tdj ), Mf&Tf ) ) // NAND
        {
          action_t<TT> act = { DecAct_t::CL_OR_, {0, i, j }, { Mf }, { Tf }, nbits, 1 };
          //printf("COST %d bec OR\n", act.cost );
          set_closure.push_back( act );
        }
        else if( kitty::equal( Mf&( ~(Tdi | Tdj) ), Mf&Tf ) ) // NAND
        {
          action_t<TT> act = { DecAct_t::CL_NOR_, {0, i, j }, { Mf }, { Tf }, nbits, 1 };
          //printf("COST %d bec NOR\n", act.cost );
          set_closure.push_back( act );
        }
        else if( kitty::equal( Mf&( ~Tdi & Tdj ), Mf&Tf ) ) // NAND
        {
          action_t<TT> act = { DecAct_t::CL_LT_, {0, i, j }, { Mf }, { Tf }, nbits, 1 };
          //printf("COST %d bec LT\n", act.cost );
          set_closure.push_back( act );
        }
        else if( kitty::equal( Mf&( Tdi& ~Tdj), Mf&Tf ) )// NAND
        {
          action_t<TT> act = { DecAct_t::CL_LT_, {0, j, i }, { Mf }, { Tf }, nbits, 1 };
          //printf("COST %d bec LT\n", act.cost );
          set_closure.push_back( act );
        }
        else if( kitty::equal( Mf&( ~Tdi | Tdj ), Mf&Tf ) ) // NAND
        {
          action_t<TT> act = { DecAct_t::CL_LE_, {0, i, j }, { Mf }, { Tf }, nbits, 1 };
          //printf("COST %d bec LE\n", act.cost );
          set_closure.push_back( act );
        }
        else if( kitty::equal( Mf&( Tdi | ~Tdj), Mf&Tf ) )// NAND
        {
          action_t<TT> act = { DecAct_t::CL_LE_, {0, j, i }, { Mf }, { Tf }, nbits, 1 };
          //printf("COST %d bec LE\n", act.cost );
          set_closure.push_back( act );
        }
        else if( kitty::equal( Mf&( Tdi ^ Tdj ), Mf&Tf ) ) // NAND
        {
          action_t<TT> act = { DecAct_t::CL_XOR_, {0, i, j }, { Mf }, { Tf }, nbits, std::is_same<Ntk, xag_network>::value ? 1 : 3 };
          //printf("COST %d bec NAND\n", act.cost );
          set_closure.push_back( act );
        }
        else if( kitty::equal( Mf&( Tdi ^ ~Tdj), Mf&Tf ) )// NAND
        {
          action_t<TT> act = { DecAct_t::CL_XNOR_, {0, i, j }, { Mf }, { Tf }, nbits, std::is_same<Ntk, xag_network>::value ? 1 : 3 };
          //printf("COST %d bec NAND\n", act.cost );
          set_closure.push_back( act );
        }
      }
    }

  }
  return false;

}


template<class TT, class Ntk>
bool DecAnalyzer<TT, Ntk>::check_satX( std::vector<signal_t> boolean_chain, int newBound )
{
  int COST_XOR = std::is_same<Ntk, xag_network>::value ? 1 : 3;
  TT Tf = * pNet->getFuncP( (*pY)[0] );
  TT Mf = * pNet->getMaskP( (*pY)[0] );
  int nbits = kitty::count_zeros( Mf& ~Mf );
  int r = boolean_chain.size();
  for( int d{0}; d<r; ++d )
  {
    signal_t xd = boolean_chain[d];
    TT Td = * pNet->getFuncP( xd );

    for( int i{r-1}; i>=newBound; --i )
    {
      signal_t xi = boolean_chain[i];
      TT Tdi = * pNet->getFuncP( xi );

      for( int j{i-1}; j>=0; --j )
      {
        signal_t xj = boolean_chain[j];
        TT Tdj = * pNet->getFuncP( xi );

        if( kitty::equal( Mf&( Tdi & Tdj ), Mf&(Td^Tf) ) ) // AND
        {
          action_t<TT> act = { DecAct_t::XCL_AND_, {0, i, j, d }, { Mf }, { Tf }, nbits, 1+COST_XOR };
          //printf("COST %d bec AND\n", act.cost );
          set_closure.push_back( act );
        }
        else if( kitty::equal( Mf&( ~(Tdi & Tdj) ), Mf&(Td^Tf)  ) ) // NAND
        {
          action_t<TT> act = { DecAct_t::XCL_NAND_, {0, i, j, d }, { Mf }, { Tf }, nbits, 1+COST_XOR };
          //printf("COST %d bec NAND\n", act.cost );
          set_closure.push_back( act );
        }
        else if( kitty::equal( Mf&( Tdi | Tdj ), Mf&(Td^Tf)  ) ) // NAND
        {
          action_t<TT> act = { DecAct_t::XCL_OR_, {0, i, j, d }, { Mf }, { Tf }, nbits, 1+COST_XOR };
          //printf("COST %d bec OR\n", act.cost );
          set_closure.push_back( act );
        }
        else if( kitty::equal( Mf&( ~(Tdi | Tdj) ), Mf&(Td^Tf)  ) ) // NAND
        {
          action_t<TT> act = { DecAct_t::XCL_NOR_, {0, i, j, d }, { Mf }, { Tf }, nbits, 1+COST_XOR };
          //printf("COST %d bec NOR\n", act.cost );
          set_closure.push_back( act );
        }
        else if( kitty::equal( Mf&( ~Tdi & Tdj ), Mf&(Td^Tf)  ) ) // NAND
        {
          action_t<TT> act = { DecAct_t::XCL_LT_, {0, i, j, d }, { Mf }, { Tf }, nbits, 1+COST_XOR };
          //printf("COST %d bec LT\n", act.cost );
          set_closure.push_back( act );
        }
        else if( kitty::equal( Mf&( Tdi& ~Tdj), Mf&(Td^Tf)  ) )// NAND
        {
          action_t<TT> act = { DecAct_t::XCL_LT_, {0, j, i }, { Mf }, { Tf }, nbits, 1+COST_XOR };
          //printf("COST %d bec LT\n", act.cost );
          set_closure.push_back( act );
        }
        else if( kitty::equal( Mf&( ~Tdi | Tdj ), Mf&(Td^Tf)  ) ) // NAND
        {
          action_t<TT> act = { DecAct_t::XCL_LE_, {0, i, j, d }, { Mf }, { Tf }, nbits, 1+COST_XOR };
          //printf("COST %d bec LE\n", act.cost );
          set_closure.push_back( act );
        }
        else if( kitty::equal( Mf&( Tdi | ~Tdj), Mf&(Td^Tf)  ) )// NAND
        {
          action_t<TT> act = { DecAct_t::XCL_LE_, {0, j, i }, { Mf }, { Tf }, nbits, 1+COST_XOR };
          //printf("COST %d bec LE\n", act.cost );
          set_closure.push_back( act );
        }
        else if( kitty::equal( Mf&( Tdi ^ Tdj ), Mf&(Td^Tf)  ) ) // NAND
        {
          action_t<TT> act = { DecAct_t::XCL_XOR_, {0, i, j, d }, { Mf }, { Tf }, nbits, 2*COST_XOR };
          //printf("COST %d bec NAND\n", act.cost );
          set_closure.push_back( act );
        }
        else if( kitty::equal( Mf&( Tdi ^ ~Tdj), Mf&(Td^Tf)  ) )// NAND
        {
          action_t<TT> act = { DecAct_t::XCL_XNOR_, {0, i, j, d }, { Mf }, { Tf }, nbits, 2*COST_XOR };
          //printf("COST %d bec NAND\n", act.cost );
          set_closure.push_back( act );
        }
      }
    }
  }
  return false;

}



#pragma endregion

#pragma region visualize
template<class TT, class Ntk>
void DecAnalyzer<TT, Ntk>::print_actions( std::vector<action_t<TT>> actions )
{
  printf("==========================================================================================\n");
  int cMV = 0;
  for( auto act : actions )
  {
    switch ( act.type )
    {
    case DecAct_t::ERASE_:
      printf( "%3d | targ(%2d ):    %d --x \n", cMV++, act.sigs[0], act.sigs[1] );
      break;
    case DecAct_t::D1_AND_:
      printf( "%3d | targ(%2d ): %4d  and R : %d\n", cMV++, act.sigs[0], act.sigs[1], act.reward );
      break;
    case DecAct_t::D1_OR_:
      printf( "%3d | targ(%2d ): %4d  or  R : %d\n", cMV++, act.sigs[0], act.sigs[1], act.reward );
      break;
    case DecAct_t::D1_LT_:
      printf( "%3d | targ(%2d ): %4d' and R : %d\n", cMV++, act.sigs[0], act.sigs[1], act.reward );
      break;
    case DecAct_t::D1_LE_:
      printf( "%3d | targ(%2d ): %4d' or  R : %d\n", cMV++, act.sigs[0], act.sigs[1], act.reward );
      break;
    case DecAct_t::D1_XOR_:
      printf( "%3d | targ(%2d ): %4d  xor R : %d\n", cMV++, act.sigs[0], act.sigs[1], act.reward );
      break;
    case DecAct_t::T1_AND_:
      printf( "%3d | targ(%2d ): %4d  and R : %d\n", cMV++, act.sigs[0], act.sigs[1], act.reward );
      break;
    case DecAct_t::T1_OR_:
      printf( "%3d | targ(%2d ): %4d  or  R : %d\n", cMV++, act.sigs[0], act.sigs[1], act.reward );
      break;
    case DecAct_t::T1_LT_:
      printf( "%3d | targ(%2d ): %4d' and R : %d\n", cMV++, act.sigs[0], act.sigs[1], act.reward );
      break;
    case DecAct_t::T1_LE_:
      printf( "%3d | targ(%2d ): %4d' or  R : %d\n", cMV++, act.sigs[0], act.sigs[1], act.reward );
      break;
    case DecAct_t::T1_XOR_:
      printf( "%3d | targ(%2d ): %4d  xor R : %d\n", cMV++, act.sigs[0], act.sigs[1], act.reward );
      break;
    case DecAct_t::NES_:
      printf( "%3d | targ(%2d ):    NES{%2d, %2d }  :  %4d    %s\n",cMV++, act.sigs[0], act.sigs[1], act.sigs[2], act.reward, 
              ( act.id_ord == 0u ? "{ j ; i } -> {  or(  j,  i ); and(  j,  i ) }" : "{ j ; i } -> { and(  j,  i );  or(  j,  i ) }" ) );
      break;
    case DecAct_t::ES_:
      printf( "%3d | targ(%2d ):     ES{%2d, %2d }  :  %4d    %s\n", cMV++, act.sigs[0], act.sigs[1], act.sigs[2], act.reward, 
              ( act.id_ord == 0u ? "{ j ; i } -> {  or(  j, ~i );  or( ~j,  i ) }" : "{ j ; i } -> { and(  j, ~i ); and( ~j,  i ) }" ) );
      break;
    case DecAct_t::SVS_:
      if( act.id_sym == 0u )
        {
          printf( "%3d | targ(%2d ):    { SVS %2d }%2d' :  %4d    %s\n", cMV++, act.sigs[0], act.sigs[2], act.sigs[1], act.reward, 
                  ( act.id_ord == 0u ? "{ j ; i } -> {  or(  j, ~i );        i      }" : "{ j ; i } -> { and(  j,  i );        i      }" ) );
        }
        else if( act.id_sym == 1u )
        {
          printf( "%3d | targ(%2d ):    { SVS %2d }%2d  :  %4d    %s\n", cMV++, act.sigs[0], act.sigs[2], act.sigs[1], act.reward, 
                  ( act.id_ord == 0u ? "{ j ; i } -> {  or( ~j, ~i );        i      }" : "{ j ; i } -> { and(  j, ~i );        i      }" ) );
        }
        else if( act.id_sym == 2u )
        {
          printf( "%3d | targ(%2d ):    { SVS %2d }%2d  :  %4d    %s\n", cMV++, act.sigs[0], act.sigs[1], act.sigs[2], act.reward, 
                  ( act.id_ord == 0u ? "{ j ; i } -> {        j     ;  or( ~j,  i ) }" : "{ j ; i } -> {        j     ; and(  j,  i ) }" ) );
        }
        else if( act.id_sym == 3u )
        {
          printf( "%3d | targ(%2d ):    { SVS %2d }%2d  :  %4d    %s\n", cMV++, act.sigs[0], act.sigs[1], act.sigs[2], act.reward, 
                  ( act.id_ord == 0u ? "{ j ; i } -> {        j     ;  or(  j,  i ) }" : "{ j ; i } -> {        j     ; and( ~j,  i ) }" ) );
        }
        break;
    case DecAct_t::MS_:
      if( act.id_ord == 0u )
      {
        printf( "%3d | targ(%2d ):     MS{%2d, %2d }  :  %4d    %s\n", cMV++, act.sigs[0], act.sigs[1], act.sigs[2], act.reward, 
                "{ j ; i } -> {              ; xor( ~j,  i ) }" );
      }
      else if( act.id_ord == 1u )
      {
        printf( "%3d | targ(%2d ):     MS{%2d, %2d }  :  %4d    %s\n", cMV++, act.sigs[0], act.sigs[1], act.sigs[2], act.reward, 
                "{ j ; i } -> { xor(  j,  i );               }" );
      }
      else if( act.id_ord == 2u )
      {
        printf( "%3d | targ(%2d ):     MS{%2d, %2d }  :  %4d    %s\n", cMV++, act.sigs[0], act.sigs[1], act.sigs[2], act.reward, 
                "{ j ; i } -> { xor( ~j,  i );               }" );
      }
      else if( act.id_ord == 3u )
      {
        printf( "%3d | targ(%2d ):     MS{%2d, %2d }  :  %4d    %s\n", cMV++, act.sigs[0], act.sigs[1], act.sigs[2], act.reward, 
                "{ j ; i } -> {              ; xor( ~j,  i ) }" );
      }        
      break;
    case DecAct_t::CSVS_: 
      if( act.id_sym == 0 )
      {
          printf( "%3d | targ(%2d ):    CSVS{ %1d', %2d'}:  %4d    %s\n", cMV++, act.sigs[0], act.sigs[1], act.sigs[2], act.reward, 
                  ( act.id_ord == 0u ? "{ j ; i } -> {              ; and(  j,  i ) }" : "{ j ; i } -> { and( ~j,  i ) ;              }" ) );
      }
      else if( act.id_sym == 1 )
      {
         printf( "%3d | targ(%2d ):    CSVS{ %1d', %2d }:  %4d    %s\n", cMV++, act.sigs[0], act.sigs[1], act.sigs[2], act.reward, 
                  ( act.id_ord == 0u ? "{ j ; i } -> {              ;  or( ~j,  i ) }" : "{ j ; i } -> { and(  j, ~i ) ;              }" ) );
      }
      else if( act.id_sym == 2 )
      {
        printf( "%3d | targ(%2d ):    CSVS{ %2d, %1d' }:  %4d    %s\n", cMV++, act.sigs[0], act.sigs[1], act.sigs[2], act.reward, 
                  ( act.id_ord == 0u ? "{ j ; i } -> {  or( j, ~i ) ;              }" : "{ j ; i } -> {              ; and( ~j,  i ) }" ) );
      }
      else if( act.id_sym == 3 )
      {
        printf( "%3d | targ(%2d ):    CSVS{ %2d, %2d }:  %4d    %s\n", cMV++, act.sigs[0], act.sigs[1], act.sigs[2], act.reward, 
                  ( act.id_ord == 0u ? "{ j ; i } -> {  or( j,  i ) ;              }" : "{ j ; i } -> {              ;  or(  j,  i ) }" ) );
      }
      break;
    case DecAct_t::SD_:
      printf( "%3d | targ(%2d ):   f=%1d'f0+ %1d'f1:  %4d    %s\n", cMV++, act.sigs[0], act.sigs[1],act.sigs[1], act.reward, 
                  ( "ite( i, f1, f0 )" ) );

    default:
      break;
    }
  }
  printf("==========================================================================================\n");

}

#pragma endregion visualize


} // namespace check_