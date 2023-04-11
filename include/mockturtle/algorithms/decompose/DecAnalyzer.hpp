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
  /* symmetry-remappings */
  NES_,
  ES_,
  MS_,
  SVS_,
  CSVS_,
  /* termination */
  BUF_,
  INV_
};

template<class TT>
struct action_t
{
  DecAct_t type;
  std::vector<int> sigs; // signals characterizing the move TD [x,T]
  TT func;                 // setmaskapped truth
  TT mask;                 // setmaskapped mask
  int reward;
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
  std::vector<int>      *pV;
  std::vector<action_t<TT>> set_topdec;
  std::vector<action_t<TT>> set_remove;
  std::vector<action_t<TT>> set_remap;
  std::vector<action_t<TT>> set_closure;

public:
  DecAnalyzer(DecNet<TT, Ntk> *, std::vector<signal_t> *, std::vector<signal_t> *, std::vector<int> *);
  ~DecAnalyzer();
  /* remapping */
  action_t<TT> simple_remapping( int, int, int, int, int, DecAct_t, int );
  action_t<TT> multiform_remapping( int, int, int, int, DecAct_t );
  action_t<TT> compatible_remapping( int, int, int, int, int, int, DecAct_t, int );

  /* analyze */
  void check0();
  void check1();
  void check2();
  std::vector<action_t<TT>> get_topdec();
  std::vector<action_t<TT>> get_remove();
  std::vector<action_t<TT>> get_remap();
  std::vector<action_t<TT>> get_closure();
  void print_actions(std::vector<action_t<TT>>);
};

#pragma region propefuncies
template<class TT, class Ntk> 
DecAnalyzer<TT, Ntk>::DecAnalyzer( DecNet<TT, Ntk> * pNet, std::vector<signal_t> (*pX), std::vector<signal_t> ((*pY)), std::vector<int> (*pV)):
pNet(pNet),
pX(pX),
pY(pY),
pV(pV)
{}
template<class TT, class Ntk> DecAnalyzer<TT, Ntk>::~DecAnalyzer(){}

template<class TT, class Ntk> std::vector<action_t<TT>> DecAnalyzer<TT, Ntk>::get_topdec(){ return set_topdec; }
template<class TT, class Ntk> std::vector<action_t<TT>> DecAnalyzer<TT, Ntk>::get_remove(){ return set_remove; }
template<class TT, class Ntk> std::vector<action_t<TT>> DecAnalyzer<TT, Ntk>::get_remap(){ return set_remap; }
template<class TT, class Ntk> std::vector<action_t<TT>> DecAnalyzer<TT, Ntk>::get_closure(){ return set_closure; }
#pragma endregion

#pragma region termination
template<class TT, class Ntk>
void DecAnalyzer<TT, Ntk>::check0()
{
  for( int iDiv{0}; iDiv < pX->size(); ++iDiv )
  {
    for( int iTrg{0}; iTrg < pY->size(); ++iTrg )
    {
      TT * pFT = pNet->getFuncP( ((*pY))[iTrg] );
      TT * pFM = pNet->getMaskP( ((*pY))[iTrg] );
      TT * pDT = pNet->getFuncP( (*pX)[iDiv] );
      TT * pDM = pNet->getMaskP( (*pX)[iDiv] );
      
      if( kitty::equal( *pFT & *pFM, *pDT & *pDM ) )
      {
        action_t<TT> act = { DecAct_t::BUF_, {iTrg, iDiv}, *pFM, *pFT, kitty::count_zeros( *pFM & ~(*pFM) ) };
        set_closure.push_back( act );
      }
      else if( kitty::equal( (~*pFT) & *pFM, *pDT & *pFM ) )
      {
        action_t<TT> act = { DecAct_t::INV_, {iTrg, iDiv}, *pFM, *pFT, kitty::count_zeros( *pFM & ~(*pFM) ) };
        set_closure.push_back( act );
      }
    }
  }
}
#pragma endregion

#pragma region decomposability
template<class TT, class Ntk>
void DecAnalyzer<TT, Ntk>::check1()
{
  for( int iDiv{0}; iDiv < pX->size(); ++iDiv )
  {
    for( int iTrg{0}; iTrg < pY->size(); ++iTrg )
    {
      TT * pFT = pNet->getFuncP( ((*pY))[iTrg] );
      TT * pFM = pNet->getMaskP( ((*pY))[iTrg] );
      TT * pDT = pNet->getFuncP( (*pX)[iDiv] );
      TT * pDM = pNet->getMaskP( (*pX)[iDiv] );

      TT tt0  = kitty::cofactor0( * pFT, (*pV)[iDiv] );
      TT tt1  = kitty::cofactor1( * pFT, (*pV)[iDiv] );
      TT mk0  = kitty::cofactor0( * pFM, (*pV)[iDiv] );
      TT mk1  = kitty::cofactor1( * pFM, (*pV)[iDiv] );

      bool eq0_0 = is_const0( *pFM & ~*pDT & *pFT );
      bool eq0_1 = is_const0( *pFM & ~*pDT & ~*pFT );
      bool eq1_0 = is_const0( *pFM & *pDT & *pFT );
      bool eq1_1 = is_const0( *pFM & *pDT & ~*pFT );
      bool eq10_ = equal( mk0 & mk1 & ~tt1 , mk0 & mk1 & tt0 );
      bool mask = equal( mk0 & mk1 & tt1 , mk0 & mk1 & tt0 );

      if( mask )
      {
        action_t<TT> act = { DecAct_t::ERASE_, {iTrg, iDiv}, *pFM, *pFT, kitty::count_ones(~*pFM ) };
        set_remove.push_back( act );
      }
      else
      {
        if ( eq0_0 ) // F0 = 0  =>  F = Xi & F1
        {
          action_t<TT> act = { DecAct_t::D1_AND_, {iTrg, iDiv}, *pDT & *pFM, *pFT, kitty::count_ones(~(*pDT & *pFM)) };
          set_topdec.push_back( act );
        }
        if ( eq0_1 ) // F0 = 1  =>  F = Xi' + F1
        {
          action_t<TT> act = { DecAct_t::D1_LE_, {iTrg, iDiv}, *pDT & *pFM, *pFT, kitty::count_ones(~(*pDT & *pFM)) };
          set_topdec.push_back( act );
        }
        if ( eq1_0 ) // F1 = 0  =>  F = Xi' & F0
        {
          action_t<TT> act = { DecAct_t::D1_LT_, {iTrg, iDiv}, ~*pDT & *pFM, *pFT, kitty::count_ones(~(~*pDT & *pFM)) };
          set_topdec.push_back( act );
        }
        if ( eq1_1 ) // F1 = 1  =>  F = Xi + F1
        {
          action_t<TT> act = { DecAct_t::D1_OR_, {iTrg, iDiv}, ~*pDT & *pFM, *pFT, kitty::count_ones(~(~*pDT & *pFM)) };
          set_topdec.push_back( act );
        }
        if ( eq10_ )// F1 = F0' =>  F = Xi ^ F0
        {
          action_t<TT> act = { DecAct_t::D1_XOR_, {iTrg, iDiv}, *pFM, ( *pFT & *pFM & ~*pDT )|( ~*pFT & *pFM & *pDT ), kitty::count_ones(*pFM) };
          set_topdec.push_back( act );      
        }
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
  TT cofactorG( TT * fn, uint32_t G, uint32_t i, uint32_t j )
  {
    switch ( G )
    {
      case 0u: return cofactor0(cofactor0( *fn, j), i ); /* F00 */
      case 1u: return cofactor1(cofactor0( *fn, j), i ); /* F01 */
      case 2u: return cofactor0(cofactor1( *fn, j), i ); /* F10 */
      case 3u: return cofactor1(cofactor1( *fn, j), i ); /* F11 */
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
 template<class TT>
  TT cube_generator( int cube, TT * pDTi, TT * pDTj )
  {
    switch ( cube )
    {
    case 0: 
      return ~(*pDTj) & ~(*pDTi);
      break;
    case 1: 
      return ~(*pDTj) & *pDTi;
      break;
    case 2: 
      return *pDTj & ~(*pDTi);
      break;
    case 3: 
      return *pDTj & *pDTi;
      break;
    }
  }

template<class TT, class Ntk>
  action_t<TT> DecAnalyzer<TT, Ntk>::simple_remapping( int from, int to, int i, int j, int iTrg, DecAct_t type, int id_symmetry )
  {
    assert( i < j );
    action_t<TT> res;
    res.type = type;

    res.sigs.push_back(iTrg);
    res.sigs.push_back(i);
    res.sigs.push_back(j); 

    res.id_ord = 1u*( from > to );
    res.id_sym = id_symmetry;
    
    int xi = (*pV)[i];
    int xj = (*pV)[j];

    TT * pDTi = pNet->getFuncP( (*pX)[i] );
    TT * pDTj = pNet->getFuncP( (*pX)[j] );
    TT * pFT  = pNet->getFuncP( (*pY)[iTrg] );
    TT * pFM  = pNet->getMaskP( (*pY)[iTrg] );

    TT A = cube_generator( from, pDTi, pDTj );
    TT B = cube_generator( to, pDTi, pDTj );
    TT ttA = cofactorG( pFT, from, xi, xj );
    TT ttB = cofactorG( pFT, to, xi, xj );
    TT mkA = cofactorG( pFM, from, xi, xj );
    TT mkB = cofactorG( pFM, to, xi, xj );

    res.mask = ( ~A & *pFM ) | ( B & mkA );
    res.reward = kitty::count_zeros( res.mask );

    TT TA = A & *pFT; /* remapping for A */
    TT TB = B & ( ( mkB & *pFT ) | ( mkA & ttA ) ); /* remapping for B */
    TT TR = ( ~A & ~B & *pFT );

    res.func = TA | TB | TR ;
    
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

    uint32_t xi = (*pV)[i];
    uint32_t xj = (*pV)[j];
    
    TT * pDTi = pNet->getFuncP( (*pX)[i] );
    TT * pDTj = pNet->getFuncP( (*pX)[j] );
    TT * pFT  = pNet->getFuncP( (*pY)[iTrg] );
    TT * pFM  = pNet->getMaskP( (*pY)[iTrg] );

    TT A = cube_generator( from1, pDTi, pDTj );
    uint32_t to1 = uint32_t( 3 - from1 );
    TT C = cube_generator( to1, pDTi, pDTj );

    uint32_t from2 = ( from1 == 0u )*1u + ( from1 == 1u )*3u + ( from1 == 3u )*2u;

    TT B = cube_generator( from2, pDTi, pDTj  );
    uint32_t to2 = uint32_t( 3 - from2 );
    TT D = cube_generator( to2, pDTi, pDTj  );

    TT ttA = cofactorG( pFT, from1, xi, xj );
    TT ttB = cofactorG( pFT, from2, xi, xj );
    TT ttC = cofactorG( pFT, to1, xi, xj );
    TT ttD = cofactorG( pFT, to2, xi, xj );

    TT mkA = cofactorG( pFM, from1, xi, xj );
    TT mkB = cofactorG( pFM, from2, xi, xj );
    TT mkC = cofactorG( pFM, to1, xi, xj );
    TT mkD = cofactorG( pFM, to2, xi, xj );

    res.mask = ( ~B & ~A & *pFM ) | ( ( C & mkA ) | ( D & mkB ) );
    res.reward = kitty::count_zeros( res.mask );
    
    TT preserved = (~A&~B&~C&~D)&(*pFT);

    TT modifiedA = A & *pFT;
    TT modifiedB = B & *pFT;

    TT modifiedC = C & ( ( mkA & ~mkC & ttA ) | ( mkC & *pFT ) );
    TT modifiedD = D & ( ( mkB & ~mkD & ttB ) | ( mkD & *pFT ) );

    res.func = preserved | modifiedA | modifiedC | modifiedB | modifiedD;
    
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

    uint32_t xi = (*pV)[i];
    uint32_t xj = (*pV)[j];

    TT * pDTi = pNet->getFuncP( (*pX)[i] );
    TT * pDTj = pNet->getFuncP( (*pX)[j] );
    TT * pFT  = pNet->getFuncP( (*pY)[iTrg] );
    TT * pFM  = pNet->getMaskP( (*pY)[iTrg] );

    uint32_t excluded = uint32_t( 6u - from1 - from2 - to );

    TT A = cube_generator( from1, pDTi, pDTj );
    TT B = cube_generator( from2, pDTi, pDTj );
    TT C = cube_generator( to, pDTi, pDTj );
    TT D = cube_generator( excluded, pDTi, pDTj );

    TT ttA = cofactorG( pFT, from1, xi, xj );
    TT ttB = cofactorG( pFT, from2, xi, xj );
    TT ttC = cofactorG( pFT, to, xi, xj );

    TT mkA = cofactorG( pFM, from1, xi, xj );
    TT mkB = cofactorG( pFM, from2, xi, xj );
    TT mkC = cofactorG( pFM, to, xi, xj );

    res.mask = ( ~B & ~A & *pFM ) | ( C & ( mkA | mkB ) );
    
    res.reward = kitty::count_zeros( res.mask );

    TT TA = A & *pFT;
    TT TB = B & *pFT;
    TT TC = C & ( ( mkA & ttA ) | ( mkB & ttB ) | ( mkC & *pFT ) );
    TT TR = ~A & ~B & ~C & *pFT;

    res.func = TA | TB | TC | TR ;
    
    return res;
  }

template<class TT, class Ntk>
void DecAnalyzer<TT, Ntk>::check2()
{

  for ( auto j = 1u; j < pV->size(); ++j )
  {
    for ( auto i = 0u; i < j; ++i )
    {
      for( int iTrg{0}; iTrg<pY->size(); ++iTrg )
      {
        TT * pFT = pNet->getFuncP( (*pY)[iTrg] );
        TT * pFM = pNet->getMaskP( (*pY)[iTrg] );
        
        TT tti0  = kitty::cofactor0( * pFT, (*pV)[i] );
        TT tti1  = kitty::cofactor1( * pFT, (*pV)[i] );
        TT mki0  = kitty::cofactor0( * pFM, (*pV)[i] );
        TT mki1  = kitty::cofactor1( * pFM, (*pV)[i] );
        TT ttj0  = kitty::cofactor0( * pFT, (*pV)[j] );
        TT ttj1  = kitty::cofactor1( * pFT, (*pV)[j] );
        TT mkj0  = kitty::cofactor0( * pFM, (*pV)[j] );
        TT mkj1  = kitty::cofactor1( * pFM, (*pV)[j] );
        if( !equal( mkj0 & mkj1 & ttj1 , mkj0 & mkj1 & ttj0 ) )   
        {     
          if( !equal( mki0 & mki1 & tti1 , mki0 & mki1 & tti0 ) )
          {
            const auto tt00 = cofactor0( ttj0, (*pV)[i] );
            const auto tt01 = cofactor1( ttj0, (*pV)[i] );
            const auto tt10 = cofactor0( ttj1, (*pV)[i] );
            const auto tt11 = cofactor1( ttj1, (*pV)[i] );
            const auto mk00 = cofactor0( mkj0, (*pV)[i] );
            const auto mk01 = cofactor1( mkj0, (*pV)[i] );
            const auto mk10 = cofactor0( mkj1, (*pV)[i] );
            const auto mk11 = cofactor1( mkj1, (*pV)[i] );

            const auto eq01 = equal( mk00&mk01&tt00, mk00&mk01&tt01 );
            const auto eq02 = equal( mk00&mk10&tt00, mk00&mk10&tt10 );
            const auto eq03 = equal( mk00&mk11&tt00, mk00&mk11&tt11 );
            const auto eq12 = equal( mk10&mk01&tt01, mk10&mk01&tt10 );
            const auto eq13 = equal( mk01&mk11&tt01, mk01&mk11&tt11 );
            const auto eq23 = equal( mk10&mk11&tt10, mk10&mk11&tt11 );

            if ( eq12 ) // F01 = F10 NES
            {
              set_remap.push_back( simple_remapping( 1, 2, i, j, iTrg, DecAct_t::NES_, 0 ) );
              set_remap.push_back( simple_remapping( 2, 1, i, j, iTrg, DecAct_t::NES_, 0 ) );
            }
            if ( eq03 ) // F00 = F11 ES
            {
              set_remap.push_back( simple_remapping( 0, 3, i, j, iTrg, DecAct_t::ES_, 0 ) );
              set_remap.push_back( simple_remapping( 3, 0, i, j, iTrg, DecAct_t::ES_, 0 ) );
            }
            if ( eq02 ) // F00=F10
            {
              set_remap.push_back( simple_remapping( 0, 2, i, j, iTrg, DecAct_t::SVS_, 0u ) ); // SVS0X_
              set_remap.push_back( simple_remapping( 2, 0, i, j, iTrg, DecAct_t::SVS_, 0u ) );
            }
            if ( eq13 ) // F01=F11
            {
              set_remap.push_back( simple_remapping( 1, 3, i, j, iTrg, DecAct_t::SVS_, 1 ) ); // SVS1X_
              set_remap.push_back( simple_remapping( 3, 1, i, j, iTrg, DecAct_t::SVS_, 1 ) );
            }
            if ( eq01 ) // F01=F00
            {
              set_remap.push_back( simple_remapping( 0, 1, i, j, iTrg, DecAct_t::SVS_, 2 ) ); // SVSX0
              set_remap.push_back( simple_remapping( 1, 0, i, j, iTrg, DecAct_t::SVS_, 2 ) );    
            }
            if ( eq23 ) // F11=F10
            {
              set_remap.push_back( simple_remapping( 2, 3, i, j, iTrg, DecAct_t::SVS_, 3 ) ); // SVSX1_
              set_remap.push_back( simple_remapping( 3, 2, i, j, iTrg, DecAct_t::SVS_, 3 ) ); 
            }
            if ( eq12 && eq03 ) // F01=F10 and F00=F11
            {
              set_remap.push_back( multiform_remapping( 0, i, j, iTrg, DecAct_t::MS_ ) );
              set_remap.push_back( multiform_remapping( 1, i, j, iTrg, DecAct_t::MS_ ) );      
              set_remap.push_back( multiform_remapping( 2, i, j, iTrg, DecAct_t::MS_ ) );      
              set_remap.push_back( multiform_remapping( 3, i, j, iTrg, DecAct_t::MS_ ) );      
            }
            if( eq02 && eq01 && eq12 )
            {
              set_remap.push_back( compatible_remapping( 0, 1, 2, i, j, iTrg, DecAct_t::CSVS_, 0 ) ); // CSVS00_
              set_remap.push_back( compatible_remapping( 2, 0, 1, i, j, iTrg, DecAct_t::CSVS_, 0 ) );
            }
            if( eq13 && eq01 && eq03 )
            {
              set_remap.push_back( compatible_remapping( 0, 1, 3, i, j, iTrg, DecAct_t::CSVS_, 1 ) ); // CSVS10_ old ij, new ji
              set_remap.push_back( compatible_remapping( 3, 1, 0, i, j, iTrg, DecAct_t::CSVS_, 1 ) );
            }
            if( eq02 && eq23 && eq03 )
            {
              set_remap.push_back( compatible_remapping( 0, 2, 3, i, j, iTrg, DecAct_t::CSVS_, 2 ) ); // CSVS01_ old ij, new ji
              set_remap.push_back( compatible_remapping( 3, 2, 0, i, j, iTrg, DecAct_t::CSVS_, 2 ) );
            }
            if( eq13 && eq23 && eq12 )
            {
              set_remap.push_back( compatible_remapping( 1, 3, 2, i, j, iTrg, DecAct_t::CSVS_, 3 ) ); // CSVS11_
              set_remap.push_back( compatible_remapping( 3, 2, 1, i, j, iTrg, DecAct_t::CSVS_, 3 ) );
            }
          }
        }
      }
    }
  }
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

    default:
      break;
    }
  }
  printf("==========================================================================================\n");

}

#pragma endregion visualize


} // namespace mocktufuncle