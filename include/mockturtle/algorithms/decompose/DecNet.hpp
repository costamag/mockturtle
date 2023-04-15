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
  \file ChSimNet.hpp
  \brief data structure for combining simulations and nodes

  \author Andrea Costamagna
*/
#pragma once

#include <stdio.h>
#include <stack>
#include "DecSims.hpp"
#include "DecNodes.hpp"
#include "../../networks/detail/foreach.hpp"
#include "../../networks/aig.hpp"
#include "../../networks/xag.hpp"
#include <kitty/print.hpp>

namespace mockturtle
{

struct signal_t
{
  sim_t sim ;
  node_t node ;
  signal_t( uint32_t sim, uint32_t node ): sim(sim), node(node) {} // user-defined default constructor

};

template<class TT, class Ntk>
class DecNet
{
private:
  /* nodes & sims */
  DecNodes<Ntk>         nodes;
  DecSims<TT>           sims;
  /* basic envelope interface */
  std::vector<signal_t> vPIs;
  std::vector<signal_t> vPOs;
  std::vector<signal_t> vTargs;
  int                   nIns;
  int                   nOut;
  TT                    FuncOSY;
  TT                    MaskOSY;

public:
  DecNet();
  ~DecNet();

public:
  DecNet( const DecNet& );              // Declare copy constructor.
  /* modify */
public:
  signal_t create_target( const TT&, const TT& );
  void     close_target( signal_t, signal_t, int );
  signal_t create_PI( const TT& );
  signal_t create_PO( signal_t );
  signal_t create_xor( signal_t, signal_t );
  signal_t create_and( signal_t, signal_t );
  signal_t create_or( signal_t, signal_t );
  signal_t create_lt( signal_t, signal_t );
  signal_t create_le( signal_t, signal_t );
  signal_t create_ge( signal_t, signal_t );
  signal_t create_gt( signal_t, signal_t );
  signal_t create_nand( signal_t, signal_t );
  signal_t create_nor( signal_t, signal_t );
  signal_t create_xnor( signal_t, signal_t );
  signal_t create_not( signal_t );
  signal_t create_buf( signal_t );
  void     init( const std::vector<TT>&, const std::vector<TT>& );
  void change_sim_info( signal_t, TT, TT );
  void print_net();
  void print_net_rec(signal_t);

public:
  /* iterate */
  template<typename Fn> void foreach_po( Fn&& );
  template<typename Fn> void foreach_pi( Fn&& );
  template<typename Fn> void foreach_fanin( signal_t, Fn&& );

  /* read */
public:
  int numPOs();
  int numPIs();
  int numTargets();
  sim_t NodeToSim( node_t );
  signal_t NodeToSig( node_t );
  std::vector<node_t> * SimToNodes( sim_t );
  TT * getFuncP( signal_t );
  TT * getMaskP( signal_t );
  TT * getTargetFuncP( sim_t );
  TT * getTargetMaskP( sim_t );
  signal<Ntk> getNtkSig( signal_t );
  DecFunc_t getFnType( signal_t );
  std::vector<signal_t> getTargets();
  std::vector<signal_t> getPIs();

  void setOSY( TT, TT );
  TT   getFuncOSY();
  TT   getMaskOSY();

  /* properties */
public:
  bool isPI( signal_t );
  int isSynt( signal_t );
  void setSig( node_t, signal<Ntk> );

};

#pragma region constructors
template<class TT, class Ntk> DecNet<TT, Ntk>::DecNet(){}
template<class TT, class Ntk> DecNet<TT, Ntk>::~DecNet(){}
template<class TT, class Ntk> 
DecNet<TT, Ntk>::DecNet( const DecNet<TT, Ntk>& other )
{
  nodes = other.nodes;
  sims  = other.sims;
  vPIs  = other.vPIs;
  vPOs  = other.vPOs;
  vTargs = other.vTargs;
}
#pragma endregion

#pragma region copy

#pragma endregion

#pragma region explore
template<class TT, class Ntk> 
template<typename Fn> 
void DecNet<TT, Ntk>::foreach_po( Fn&& fn ){ detail::foreach_element( vPOs.begin(), vPOs.end(), fn );  }

template<class TT, class Ntk>
template<typename Fn>
void DecNet<TT, Ntk>::foreach_pi( Fn&& fn ){ detail::foreach_element( vPIs.begin(), vPIs.end(), fn );  }

template<class TT, class Ntk>
template<typename Fn>
void DecNet<TT, Ntk>::foreach_fanin( signal_t sig, Fn&& fn )
{ 
  std::vector<node_t> * pFanins = nodes.getFanInsP( sig.node );
  detail::foreach_element( pFanins->begin(), pFanins->end(), fn ); 
}
#pragma endregion

#pragma region properties
template<class TT, class Ntk> bool DecNet<TT,Ntk>::isPI( signal_t sig ){ return nodes.isPI( sig.node ); }
template<class TT, class Ntk> int DecNet<TT,Ntk>::isSynt( signal_t sig ){ return nodes.isSynt( sig.node ); }
template<class TT, class Ntk> void DecNet<TT,Ntk>::setSig( node_t node, signal<Ntk> ntk_sig ){   nodes.setSig( node, ntk_sig ); }
#pragma endregion

#pragma region read
template<class TT, class Ntk>  int DecNet<TT, Ntk>::numPOs() { return vPOs.size(); }
template<class TT, class Ntk>  int DecNet<TT, Ntk>::numPIs() { return vPIs.size(); }
template<class TT, class Ntk>  int DecNet<TT, Ntk>::numTargets() { return vTargs.size(); }
template<class TT, class Ntk>  sim_t DecNet<TT, Ntk>::NodeToSim( node_t node ) { return nodes.getSim( node ); }
template<class TT, class Ntk>  signal_t DecNet<TT, Ntk>::NodeToSig( node_t node ) { sim_t sim = nodes.getSim( node ); return signal_t{ sim, node }; }
template<class TT, class Ntk>  std::vector<node_t> * DecNet<TT, Ntk>::SimToNodes( sim_t sim ) { return sims.getNodesP( sim ); }
template<class TT, class Ntk>  TT * DecNet<TT, Ntk>::getFuncP( signal_t sig ){ return sims.getFuncP( sig.sim ); };
template<class TT, class Ntk>  TT * DecNet<TT, Ntk>::getMaskP( signal_t sig ){ return sims.getMaskP( sig.sim );  };
template<class TT, class Ntk>  TT * DecNet<TT, Ntk>::getTargetFuncP( sim_t sim ){ return sims.getFuncP( sim ); };
template<class TT, class Ntk>  TT * DecNet<TT, Ntk>::getTargetMaskP( sim_t sim ){ return sims.getMaskP( sim ); };
template<class TT, class Ntk>  signal<Ntk> DecNet<TT,Ntk>::getNtkSig( signal_t sig ){ return nodes.getNtkSig( sig.node ); };
template<class TT, class Ntk>  DecFunc_t DecNet<TT,Ntk>::getFnType( signal_t sig ){ return nodes.getFunc( sig.node ); };

#pragma endregion

#pragma region modify
template<class TT, class Ntk>
signal_t DecNet<TT, Ntk>::create_target( const TT& func, const TT& mask )
{
  sim_t  sim = sims.addSim( func, mask );
  node_t node = nodes.addHungNode( sim );
  signal_t sig { sim, node };
  return sig;
}

template<class TT, class Ntk>
void DecNet<TT, Ntk>::close_target( signal_t sTrg, signal_t sDiv, int isInv )
{
  if( isInv )
    nodes.attachHunging( sDiv.node, sDiv.sim, sTrg.node, DecFunc_t::NOT_ );
  else
    nodes.attachHunging( sDiv.node, sDiv.sim, sTrg.node, DecFunc_t::BUF_ );
}

template<class TT, class Ntk>
signal_t DecNet<TT, Ntk>::create_PI( const TT& func )
{
  sim_t  simPI = sims.addSim( func, func | ~func );
  node_t nodePI = nodes.addNode( {}, simPI, DecFunc_t::PI_ );
  signal_t sig { simPI, nodePI };
  vPIs.push_back( sig );
  return sig;
}

template<class TT, class Ntk>
signal_t DecNet<TT, Ntk>::create_PO( signal_t sig )
{
  vPOs.push_back( sig );
  return sig;
}

template<class TT, class Ntk>
signal_t DecNet<TT, Ntk>::create_not( signal_t a )
{
  sim_t  sim = sims.addSim( ~*getFuncP(a), *getMaskP(a) );
  node_t node = nodes.addNode( { a.node}, sim, DecFunc_t::NOT_ );
  signal_t sig { sim, node };
  return sig;
}

template<class TT, class Ntk>
signal_t DecNet<TT, Ntk>::create_buf( signal_t a )
{
  sim_t  sim = sims.addSim( *getFuncP(a), *getMaskP(a) );
  node_t node = nodes.addNode( { a.node}, sim, DecFunc_t::BUF_ );
  signal_t sig { sim, node };
  return sig;
}

template<class TT, class Ntk>
signal_t DecNet<TT, Ntk>::create_xor( signal_t a, signal_t b )
{
  sim_t  sim = sims.addSim( *getFuncP(a) ^ *getFuncP(b), *getMaskP(a) & *getMaskP(b) );
  node_t node = nodes.addNode( { a.node, b.node }, sim, DecFunc_t::XOR_ );
  signal_t sig { sim, node };
  return sig;
}

template<class TT, class Ntk>
signal_t DecNet<TT, Ntk>::create_and( signal_t a, signal_t b )
{
  sim_t  sim = sims.addSim( *getFuncP(a) & *getFuncP(b), *getMaskP(a) & *getMaskP(b) );
  node_t node = nodes.addNode( { a.node, b.node }, sim, DecFunc_t::AND_ );
  signal_t sig { sim, node };
  return sig;
}

template<class TT, class Ntk>
signal_t DecNet<TT, Ntk>::create_or( signal_t a, signal_t b )
{
  sim_t  sim = sims.addSim( *getFuncP(a) | *getFuncP(b), *getMaskP(a) & *getMaskP(b) );
  node_t node = nodes.addNode( { a.node, b.node }, sim, DecFunc_t::OR_ );
  signal_t sig { sim, node };
  return sig;
}

template<class TT, class Ntk>
signal_t DecNet<TT, Ntk>::create_lt( signal_t a, signal_t b )
{
  sim_t  sim = sims.addSim( ~*getFuncP(a) & *getFuncP(b), *getMaskP(a) & *getMaskP(b) );
  node_t node = nodes.addNode( { a.node, b.node }, sim, DecFunc_t::LT_ );
  signal_t sig { sim, node };
  return sig;
}

template<class TT, class Ntk>
signal_t DecNet<TT, Ntk>::create_gt( signal_t a, signal_t b )
{
  sim_t  sim = sims.addSim( *getFuncP(a) & ~*getFuncP(b), *getMaskP(a) & *getMaskP(b) );
  node_t node = nodes.addNode( { a.node, b.node }, sim, DecFunc_t::GT_ );
  signal_t sig { sim, node };
  return sig;
}

template<class TT, class Ntk>
signal_t DecNet<TT, Ntk>::create_le( signal_t a, signal_t b )
{
  sim_t  sim = sims.addSim( ~*getFuncP(a) | *getFuncP(b), *getMaskP(a) & *getMaskP(b) );
  node_t node = nodes.addNode( { a.node, b.node }, sim, DecFunc_t::LE_ );
  signal_t sig { sim, node };
  return sig;
}

template<class TT, class Ntk>
signal_t DecNet<TT, Ntk>::create_ge( signal_t a, signal_t b )
{
  sim_t  sim = sims.addSim( *getFuncP(a) | ~*getFuncP(b), *getMaskP(a) & *getMaskP(b) );
  node_t node = nodes.addNode( { a.node, b.node }, sim, DecFunc_t::GE_ );
  signal_t sig { sim, node };
  return sig;
}

template<class TT, class Ntk>
signal_t DecNet<TT, Ntk>::create_nand( signal_t a, signal_t b )
{
  sim_t  sim = sims.addSim( ~(*getFuncP(a) & *getFuncP(b)), *getMaskP(a) & *getMaskP(b) );
  node_t node = nodes.addNode( { a.node, b.node }, sim, DecFunc_t::NAND_ );
  signal_t sig { sim, node };
  return sig;
}

template<class TT, class Ntk>
signal_t DecNet<TT, Ntk>::create_nor( signal_t a, signal_t b )
{
  sim_t  sim = sims.addSim( ~(*getFuncP(a) | *getFuncP(b)), *getMaskP(a) & *getMaskP(b) );
  node_t node = nodes.addNode( { a.node, b.node }, sim, DecFunc_t::NOR_ );
  signal_t sig { sim, node };
  return sig;
}

template<class TT, class Ntk>
signal_t DecNet<TT, Ntk>::create_xnor( signal_t a, signal_t b )
{
  sim_t  sim = sims.addSim( ~(*getFuncP(a) ^ *getFuncP(b)), *getMaskP(a) & *getMaskP(b) );
  node_t node = nodes.addNode( { a.node, b.node }, sim, DecFunc_t::XNOR_ );
  signal_t sig { sim, node };
  return sig;
}

template<class TT, class Ntk> std::vector<signal_t> DecNet<TT, Ntk>::getTargets(){ return vTargs; }
template<class TT, class Ntk> std::vector<signal_t> DecNet<TT, Ntk>::getPIs(){ return vPIs; }


template<class TT, class Ntk>
void DecNet<TT, Ntk>::setOSY( TT func, TT mask )
{
  FuncOSY = func;
  MaskOSY = mask;
}

template<class TT, class Ntk>
TT DecNet<TT, Ntk>::getFuncOSY( )
{
  return FuncOSY;
}

template<class TT, class Ntk>
TT DecNet<TT, Ntk>::getMaskOSY( )
{
  return MaskOSY;
}

template<class TT, class Ntk>
void DecNet<TT, Ntk>::init( const std::vector<TT>& vTruths, const std::vector<TT>& vMasks )
{
  assert( vTruths.size() == vMasks.size() );
  assert( vTruths[0].num_vars() == vMasks[0].num_vars() );
  nIns = vTruths[0].num_vars();
  nOut = vTruths.size();
  /* create normal divisors */
	TT  x = vTruths[0].construct();
  for( int i{0}; i < nIns; ++i )
  {
	  kitty::create_nth_var( x, i );
    create_PI( x );
  }
  /* create targets */
  for( int i{0}; i<nOut; ++i )
  {
    signal_t target = create_target( vTruths[i], vMasks[i] );
    vTargs.push_back( target );
    create_PO( target );
  }
  /*  */
}

template<class TT, class Ntk>
void DecNet<TT, Ntk>::change_sim_info( signal_t sig, TT func, TT mask )
{
  sims.change_mask( sig.sim, mask );
  sims.change_func( sig.sim, func );
}


#pragma endregion modify

#pragma region print
template<class TT, class Ntk>
void DecNet<TT, Ntk>::print_net_rec( signal_t sig )
{
        foreach_fanin( sig, [&]( node_t x ) 
        {   
            signal_t child = NodeToSig(x);
            print_net_rec( child );
        } );

        switch ( getFnType( sig ) )
        {
            case DecFunc_t::NOT_:
                printf("%d = NOT ", sig.node );
                break;
            case DecFunc_t::BUF_:
                printf("%d = BUF", sig.node );
                break;
            case DecFunc_t::AND_:
                printf("%d = AND ", sig.node );
                break;
            case DecFunc_t::NAND_:
                printf("%d = NAND ", sig.node );
                break;
            case DecFunc_t::OR_:
                printf("%d = OR ", sig.node );
                break;
            case DecFunc_t::NOR_:
                printf("%d = NOR ", sig.node );
                break;
            case DecFunc_t::XOR_:
                printf("%d = XOR ", sig.node );
                break;
            case DecFunc_t::XNOR_:
                printf("%d = XNOR ", sig.node );
                break;
            case DecFunc_t::LT_:
                printf("%d = LT ", sig.node );
                break;
            case DecFunc_t::GE_:
                printf("%d = GE ", sig.node );
                break;
            case DecFunc_t::LE_:
                printf("%d = LE ", sig.node );
                break;
            case DecFunc_t::GT_:
                printf("%d = GT ", sig.node );
                break;
            default:
                break;
        }
        
        foreach_fanin( sig, [&]( node_t x ) 
        {   
          printf(" %d ", x);
        } );
        printf("\n");

}

template<class TT, class Ntk>
void DecNet<TT, Ntk>::print_net()
{
  printf("INPUTS\n");
    foreach_pi( [&]( const auto& x, auto index ) 
    {    
      printf("%d: id %d ", index, x.node );
      kitty::print_binary( *getFuncP(x) );
      printf("\n");
    } );
    foreach_po( [&]( const auto& x, auto index ) 
    {   
        print_net_rec( x );
    } );
    printf("OUTPUTS\n");
    foreach_po( [&]( const auto& x, auto index ) 
    {   
      printf("%d: id %d ", index, x.node );
      kitty::print_binary( *getFuncP(x) );
      printf("\n");
    } );
}
#pragma endregion print

} // namespace mockturtle