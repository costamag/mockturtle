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
  \file DecSims.hpp
  \brief data structure for synthesizing an network

  \author Andrea Costamagna
*/

#pragma once

#include <stdio.h>
#include <kitty/print.hpp>
#include "DecNet.hpp"
#include "DecAnalyzer.hpp"
#include "DecChsToGraph.hpp"
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <random>
#include <ctime>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"


namespace mockturtle
{

template<class TT, class Ntk>
class DecSolver
{
private:
    std::vector<TT>  vTruths;
    std::vector<TT>  vMasks;
    /* solver view */
    // using only the inputs
    std::vector<signal_t>               X ; // ( remapped ) inputs signals |X| = n
    std::vector<std::vector<int>>       vS;
    // using divisors
    std::vector<std::vector<signal_t>>  vD;
    std::vector<signal_t>               boolean_chain;

    std::vector<signal_t> Y; // targets signals             |Y| = m
    TT                    remainder;
    TT                    mask;

    std::random_device rand_dev;

public:
    DecSolver( const std::vector<TT>&, const std::vector<TT>& );
    ~DecSolver();
    /* solve */
    Ntk man_sym_solve();
    Ntk man_sym_solve_rs();
    Ntk aut_sym_solve_rs( int );
    Ntk aut_sym_solve_xor( int );
    Ntk ccg_relax( int );
    Ntk ccg_spectral( int, int );
    Ntk ccg_xor( int );
    Ntk ccgX( int );
    Ntk man_decsym_solve();
    Ntk aut_sym_solve(int);
    Ntk aut_symGT_solve(int);
    Ntk man_rdec_solve();
    Ntk aut_rdec_solve(int);

    action_t<TT> select_uniformly( std::vector<action_t<TT>> * );

    void remap_nes( DecNet<TT, Ntk> *, action_t<TT> );
    void remap_es( DecNet<TT, Ntk> *, action_t<TT> );
    void remap_svs( DecNet<TT, Ntk> *, action_t<TT> );
    void remap_csvs( DecNet<TT, Ntk> *, action_t<TT> );
    void remap_ms( DecNet<TT, Ntk> *, action_t<TT> );
    void remap_and( DecNet<TT, Ntk> *, action_t<TT> );
    void remap_or( DecNet<TT, Ntk> *, action_t<TT> );
    void remap_lt( DecNet<TT, Ntk> *, action_t<TT> );
    void remap_le( DecNet<TT, Ntk> *, action_t<TT> );
    void remap_str( DecNet<TT, Ntk> *, action_t<TT> );
    void remap( DecNet<TT, Ntk> *, action_t<TT> );
    void close( DecNet<TT, Ntk> *, action_t<TT> );
    void decompose( DecNet<TT, Ntk> *, action_t<TT> );
    void close_div( DecNet<TT, Ntk> *, std::vector<action_t<TT>> );
    void close_tar( DecNet<TT, Ntk> *, std::vector<action_t<TT>> );
    void so_terminate( DecNet<TT, Ntk> *, int );

    bool check_sat( DecNet<TT, Ntk> *, int, int ); // WIP
    bool check_satX( DecNet<TT, Ntk> *, int, int ); // WIP
    bool check_2sat( DecNet<TT, Ntk> *, signal_t, signal_t );
    bool check_2satX( DecNet<TT, Ntk> *, signal_t, signal_t, signal_t );

    /* visualize */
    void PrintSpecs();
    void show_state( DecNet<TT, Ntk> *, std::vector<signal_t> * );

};

template<class TT, class Ntk>
DecSolver<TT, Ntk>::DecSolver( const std::vector<TT>& vTruths, const std::vector<TT>& vMasks ) :
vTruths(vTruths),
vMasks( vMasks )
{
}

template<class TT, class Ntk>
DecSolver<TT, Ntk>::~DecSolver()
{
}

#pragma region solve

template<class TT, class Ntk>
Ntk DecSolver<TT, Ntk>::man_sym_solve()
{
  Ntk ntk;
  Ntk ntkBest;
  int nBest = 100000;
  for( int it{0}; it<1; ++it )
  {
    /* init */
    DecNet<TT, Ntk>  net;
    net.init( vTruths, vMasks );
    /* info on the outputs */
    Y = net.getTargets();
    //assert( Y.size()==1 );
    net.setOSY( *net.getFuncP(Y[0]), *net.getMaskP(Y[0]) );
    /* info on the inputs */
    X = net.getPIs();
    std::vector<int> S;
    for( int i{0}; i < X.size(); ++i )  
      S.push_back( i );
    vD={};
    vS = {};
    for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
    {
      vD.push_back( X );
      vS.push_back( S );
    }

    /* solve */
    while( vS[0].size() > 1 )
    {
      show_state( &net, &Y );

      DecAnalyzer<TT, Ntk> checker( &net, &X, &Y, &vS, &vD );

      //checker.check_divclosure();
      //std::vector<action_t<TT>> CS = checker.get_divclosure();
      //if( CS.size() > 0 )
      //  close_div( &net, CS );
      //printf("|CS|=%d\n", CS.size());
      if( vS[0].size() > 1 )
      {
        checker.check2();
        std::vector<action_t<TT>> RP = checker.get_remap();
        checker.print_actions( RP );
        int iEnd = RP.size()-1;
        printf( "Choose the move[%2d-%3d]:", 0, iEnd  );
        int MV;
        std::cin >> MV;
        remap( &net, RP[MV] );
      }
    }
    show_state( &net, &Y );
    so_terminate( &net, 0 );

    net.print_net();
    /* convert result */
    DecChsToGraph<TT, Ntk> conv( net );  
    ntk = conv.convert();
    printf("num gates = %d", ntk.num_gates());
    //if( ntk.num_gates() < nBest )
    //    ntkBest = ntk;
  }
  return ntk;
}


template<class TT, class Ntk>
Ntk DecSolver<TT, Ntk>::man_sym_solve_rs()
{
  Ntk ntk;
  Ntk ntkBest;
  int nBest = 100000;
  for( int it{0}; it<1; ++it )
  {
    /* init */
    DecNet<TT, Ntk>  net;
    net.init( vTruths, vMasks );
    /* info on the outputs */
    Y = net.getTargets();
    assert( Y.size()==1 );
    net.setOSY( *net.getFuncP(Y[0]), *net.getMaskP(Y[0]) );
    /* info on the inputs */
    X = net.getPIs();
    std::vector<int> S;
    for( int i{0}; i < X.size(); ++i )  
    {
      boolean_chain.push_back( X[i] );
      S.push_back( i );
    }
    vD={};
    vS = {};
    for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
    {
      vD.push_back( X );
      vS.push_back( S );
    }

    bool is_sat = false;
    /* solve */
    while( !is_sat )
    {
      show_state( &net, &Y );
      is_sat = check_sat( &net, 0, 0 );

      DecAnalyzer<TT, Ntk> checker( &net, &X, &Y, &vS, &vD );
      if( !is_sat )
      {
        checker.check2();
        std::vector<action_t<TT>> RP = checker.get_remap();
        checker.print_actions( RP );
        int iEnd = RP.size()-1;
        printf( "Choose the move[%2d-%3d]:", 0, iEnd  );
        int MV;
        std::cin >> MV;
        remap( &net, RP[MV] );
      }
    }
    //so_terminate( &net, 0 );

    net.print_net();
    /* convert result */
    DecChsToGraph<TT, Ntk> conv( net );  
    ntk = conv.convert();
    printf("num gates = %d", ntk.num_gates());
    //if( ntk.num_gates() < nBest )
    //    ntkBest = ntk;
  }
  return ntk;
}


template<class TT, class Ntk>
Ntk DecSolver<TT, Ntk>::man_decsym_solve()
{
  Ntk ntk;
  Ntk ntkBest;
  int nBest = 100000;
  for( int it{0}; it<1; ++it )
  {
    /* init */
    DecNet<TT, Ntk>  net;
    net.init( vTruths, vMasks );
    /* info on the outputs */
    Y = net.getTargets();
    assert( Y.size()==1 );
    net.setOSY( *net.getFuncP(Y[0]), *net.getMaskP(Y[0]) );
    /* info on the inputs */
    X = net.getPIs();
    std::vector<int> S;
    for( int i{0}; i < X.size(); ++i )  
      S.push_back( i );
    vD={};
    vS = {};
    for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
    {
      vD.push_back( X );
      vS.push_back( S );
    }

    /* solve */
    while( vS[0].size() > 1 )
    {
      show_state( &net, &Y );

      DecAnalyzer<TT, Ntk> checker( &net, &X, &Y, &vS, &vD );

      //checker.check_divclosure();
      //std::vector<action_t<TT>> CS = checker.get_divclosure();
      //if( CS.size() > 0 )
      //  close_div( &net, CS );
      //printf("|CS|=%d\n", CS.size());
      if( vS[0].size() > 1 )
      {
        checker.check_mixed();
        std::vector<action_t<TT>> RP = checker.get_remap();
        checker.print_actions( RP );
        int iEnd = RP.size()-1;
        printf( "Choose the move[%2d-%3d]:", 0, iEnd  );
        int MV;
        std::cin >> MV;
        remap( &net, RP[MV] );
      }
    }
    so_terminate( &net, 0 );

    net.print_net();
    /* convert result */
    DecChsToGraph<TT, Ntk> conv( net );  
    ntk = conv.convert();
    printf("num gates = %d", ntk.num_gates());
    //if( ntk.num_gates() < nBest )
    //    ntkBest = ntk;
  }
  return ntk;
}

template<class TT, class Ntk>
Ntk DecSolver<TT, Ntk>::man_rdec_solve()
{
  Ntk ntk;
  Ntk ntkBest;
  int nBest = 100000;
  for( int it{0}; it<10; ++it )
  {
    /* init */
    /* init */
    DecNet<TT, Ntk>  net;
    net.init( vTruths, vMasks );
    /* info on the outputs */
    Y = net.getTargets();
    /* info on the inputs */
    X = net.getPIs();
    std::vector<int> S;
    for( int i{0}; i < X.size(); ++i )  
      S.push_back( i );
    vS = {};
    vD={};
    for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
    {
      vD.push_back( X );
      vS.push_back( S );
    }

    /* solve */
    while( Y.size() > 0 )
    {
      DecAnalyzer<TT, Ntk> checker( &net, &X, &Y, &vS, &vD );

      checker.check_tarclosure();
      std::vector<action_t<TT>> CT = checker.get_trgclosure();
      if( CT.size() > 0 )
        close_tar( &net, CT );
      printf("|close tar|=%d\n", CT.size());      

      checker.check_divclosure();
      std::vector<action_t<TT>> CD = checker.get_divclosure();
      if( CD.size() > 0 )
        close_div( &net, CD );
      printf("|close div|=%d\n", CD.size());
      if( Y.size() > 0 )
      {
        checker.check_2stepsdec();
        std::vector<action_t<TT>> SD = checker.get_2stepsdec();
      //  checker.check2();
      //  std::vector<action_t<TT>> RP = checker.get_remap();
        checker.print_actions( SD );
        int iEnd = SD.size()-1;
        printf( "Choose the move[%2d-%3d]:", 0, iEnd  );
        int MV;
        std::cin >> MV;
        decompose( &net, SD[MV] );
      }
    }
    //net.print_net();

    /* convert result */
    DecChsToGraph<TT, Ntk> conv( net );  
    ntk = conv.convert();
    printf("**********************************\n");
    printf("num gates =%d\n", ntk.num_gates());
    printf("**********************************\n");
    //if( ntk.num_gates() < nBest )
    //    ntkBest = ntk;
  }
  return ntk;
}

template<class TT, class Ntk>
action_t<TT> DecSolver<TT, Ntk>::select_uniformly( std::vector<action_t<TT>> * pvA )
{
  std::mt19937       generator(rand_dev());
  std::uniform_int_distribution<int>  distr(0, (*pvA).size()-1);  
  int MOVE = distr(generator);
  //printf( "MOVE %d\n", MOVE );
  return (*pvA)[MOVE];
}

template<class TT, class Ntk>
Ntk DecSolver<TT, Ntk>::aut_sym_solve(int nIters)
{
  bool verbose{false};

  std::mt19937       generator(5);

  double avg_time{0};
  double avg_num_gates{0};
  int max_num_gates{0};

  std::clock_t global_start;
  double global_duration;
  global_start = std::clock();

  Ntk ntk;
  Ntk ntkBest;
  int nBest = 100000;

  int nNodesBest = 100000;
  int nActualIters = 0;
  int initial_nNodes = 0;
  for( int it{0}; it<nIters; ++it )
  {
    int nNodes  = 0;
    bool bBELOW = true;

    std::clock_t start;
    double duration;
    start = std::clock();
    bool CEC;
    
    /* init */
    DecNet<TT, Ntk>  net;
    net.init( vTruths, vMasks );
    /* info on the outputs */
    Y = net.getTargets();
    
    assert( Y.size() == 1 );
    net.setOSY( *net.getFuncP(Y[0]), *net.getMaskP(Y[0]) );

    /* info on the inputs */
    X = net.getPIs();
    std::vector<int> S;
    for( int i{0}; i < X.size(); ++i )  
      S.push_back( i );
    vS = {};
    vD = {};
    for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
    {
      vD.push_back( X );
      vS.push_back( S );
    }

    int best_reward{0};
    int best_cost{1000};
    /* solve */
    while(  (vS[0].size() > 1) ) //bBELOW &&
    {
      DecAnalyzer<TT, Ntk> checker( &net, &X, &Y, &vS, &vD );

      if( vS[0].size() > 1 )
      {
        checker.check2();
        std::vector<action_t<TT>> RP = checker.get_remap();
        best_reward=0;
        best_cost = 1000;
        std::vector<int> moves; 
        for( int iRP{0}; iRP<RP.size(); ++iRP )
        {
          if( (RP[iRP].reward > best_reward) | ( (RP[iRP].reward == best_reward) && ( nNodes+RP[iRP].cost < best_cost ) ) )
          {
            moves = {iRP};
            best_reward = RP[iRP].reward;
            best_cost = nNodes+RP[iRP].cost;
          }
          else if( (RP[iRP].reward == best_reward) && ( nNodes+RP[iRP].cost == best_cost ) )
            moves.push_back( iRP );
        }
        std::uniform_int_distribution<>  distr(0, moves.size()-1);  

        int MOVE = it == 0 ? moves[0] : moves[distr(generator)];//
        MOVE = it == 1 ? moves[moves.size()-1] : MOVE;//

        remap( &net, RP[MOVE] );
        nNodes += RP[MOVE].cost;
        if( nNodes > nNodesBest )
          bBELOW = false; 
      }
    }
    //if( bBELOW )
    {
      nActualIters++;
      so_terminate( &net, 0 );

      if( bBELOW )
        nNodesBest = nNodes;

      /* convert result */
      DecChsToGraph<TT, Ntk> conv( net ); 
  
      ntk = conv.convert();
      Y = net.getTargets();
      default_simulator<kitty::dynamic_truth_table> sim( (*net.getFuncP(Y[0])).num_vars() );
      const auto tt = simulate<kitty::dynamic_truth_table>( ntk, sim )[0];
  if(verbose)
  {
      printf("**********************************\n");
      printf("num gates =%d\n", ntk.num_gates());
  }
      for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
      {
        TT F = *net.getFuncP( Y[iTrg] );
        CEC = kitty::equal( tt, F );
        if( kitty::equal( tt, F ) )
        {
            if(verbose) printf( ANSI_COLOR_GREEN " CEC SUCCESFUL " ANSI_COLOR_RESET "\n" );
        }
        else
        {
          printf( ANSI_COLOR_RED " CEC FAILED " ANSI_COLOR_RESET "\n" );
          assert(0);
        }
      }

      duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

      avg_time += duration;
      avg_num_gates += 1.0*ntk.num_gates();
      if( ntk.num_gates() > max_num_gates )
        max_num_gates = ntk.num_gates();

if(verbose)
{ 
      printf("%5.5f\n", duration );
      printf("**********************************\n");
}
      if( (ntk.num_gates() < nBest) && CEC )
      {
        nBest = ntk.num_gates();
        ntkBest = ntk;
      }
      if( it == 0 )
        initial_nNodes = ntk.num_gates();
    }

    if( !bBELOW )
    {
      if( verbose ) printf( ANSI_COLOR_RED " LARGER THAN BEST " ANSI_COLOR_RESET "\n" );
    }
  }
  global_duration = ( std::clock() - global_start ) / (double) CLOCKS_PER_SEC;

if(verbose)
{ 
  printf("SUMMARY:\n");
  printf("%10s %10s %10s %10s %10s\n %10d %5.5f %d %5.5f %5.5f\n", 
          "min G", "E[G]", "max G", "E[t]", "T\n",
           nBest, avg_num_gates/nActualIters, max_num_gates, avg_time/nActualIters, global_duration );
  printf( "init %d, delta(max)=+%d, delta(min)=-%d\n", nBest, max_num_gates-nBest, initial_nNodes-nBest );
}
  return ntkBest;
}

template<class TT, class Ntk>
Ntk DecSolver<TT, Ntk>::aut_sym_solve_rs( int nIters )
{
  //if( vTruths[0].num_vars() < 10 )
  //  PrintSpecs();

  std::mt19937       generator(5);

  double avg_time{0};
  double avg_num_gates{0};
  int max_num_gates{0};

  Ntk ntk;
  Ntk ntkBest;
  int nBest = 100000;
  for( int it{0}; it<nIters; ++it )
  {

    std::clock_t start;
    double duration;
    start = std::clock();
    bool CEC;

    /* init */
    DecNet<TT, Ntk>  net;
    net.init( vTruths, vMasks );
    /* info on the outputs */
    Y = net.getTargets();
    assert( Y.size()==1 );

    net.setOSY( *net.getFuncP(Y[0]), *net.getMaskP(Y[0]) );

    /* info on the inputs */
    X = net.getPIs();
    std::vector<int> S;

    boolean_chain = {}; // clean the vector
    for( int i{0}; i < X.size(); ++i )  
    {
      boolean_chain.push_back( X[i] );
      S.push_back( i );
    }

    for( int i{1}; i < X.size(); ++i )  
    {
      for( int j{0}; j < i; ++j )  
      {
        boolean_chain.push_back(net.create_xor( X[i], X[j] ));
      }
    }


    vD={};
    vS = {};
    for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
    {
      vD.push_back( X );
      vS.push_back( S );
    }

    int best_reward{0};
    int old_reward{0};
    bool is_sat = false;
    bool is_stuck = false;
    /* solve */
    while( !is_sat && !is_stuck )
    {
//show_state( &net, &Y );
      is_sat = check_sat( &net, 0, 0 );

      DecAnalyzer<TT, Ntk> checker( &net, &X, &Y, &vS, &vD );
      if( !is_sat  )
      {
        checker.check2();
        checker.check_topdec( boolean_chain );
        
        std::vector<action_t<TT>> RP = checker.get_remap();
//checker.print_actions( RP );
        best_reward=0;
        std::vector<int> moves; 
//printf("b \n" );

        for( int iRP{0}; iRP<RP.size(); ++iRP )
        {
          int OLD_RWD = 0;
          if( OLD_RWD == 1 )
          {
            if( RP[iRP].reward >= old_reward )
              moves.push_back( iRP );
          }
          else
          {
            if( RP[iRP].reward > best_reward )
            {
              moves = {iRP};
              best_reward = RP[iRP].reward;
            }
            else if( RP[iRP].reward == best_reward )
              moves.push_back( iRP );
          }
        }
//printf("a \n" );
        if( moves.size() == 0 )
        {
          is_stuck = true;
        }
        else
        {
          std::uniform_int_distribution<>  distr(0, moves.size()-1);  
          int MOVE = moves[distr(generator)];
          old_reward = RP[MOVE].reward;
  //printf("MOVE %d\n", MOVE );
          remap( &net, RP[MOVE] );
        }
      }
    }
//net.print_net();
    /* convert result */
    if( is_stuck )
    {
      printf( ANSI_COLOR_RED " EMPTY SET OF MOVES " ANSI_COLOR_RESET "\n" );
    }
    else
    {
      DecChsToGraph<TT, Ntk> conv( net ); 
  
      ntk = conv.convert();
      Y = net.getTargets();
      default_simulator<kitty::dynamic_truth_table> sim( (*net.getFuncP(Y[0])).num_vars() );
      const auto tt = simulate<kitty::dynamic_truth_table>( ntk, sim )[0];
      printf("**********************************\n");
      printf("num gates =%d\n", ntk.num_gates());
      for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
      {
        TT F = *net.getFuncP( Y[iTrg] );
        CEC = kitty::equal( tt, F );
        if( kitty::equal( tt, F ) )
          printf( ANSI_COLOR_GREEN " CEC SUCCESFUL " ANSI_COLOR_RESET "\n" );
        else
        {
          printf( ANSI_COLOR_RED " CEC FAILED " ANSI_COLOR_RESET "\n" );
          assert(0);
        }
      }

      duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
      avg_time += duration/nIters;
      avg_num_gates += 1.0*ntk.num_gates()/nIters;
      if( ntk.num_gates() > max_num_gates )
        max_num_gates = ntk.num_gates();

      printf("%5.5f\n", duration );
      printf("**********************************\n");
      if( (ntk.num_gates() < nBest) && CEC )
      {
        nBest = ntk.num_gates();
        ntkBest = ntk;
      }
    }
  }
  printf("SUMMARY:\n");
  printf("%10s %10s %10s %10s\n %10d %5.5f %d %5.5f\n", 
          "min G", "E[G]", "max G", "E[t]\n",
           nBest, avg_num_gates, max_num_gates, avg_time );
  return ntkBest;
}

template<class TT, class Ntk>
Ntk DecSolver<TT, Ntk>::aut_sym_solve_xor( int nIters )
{
  printf("xor\n");
  //if( vTruths[0].num_vars() < 10 )
  //  PrintSpecs();

  std::mt19937       generator(rand_dev());

  double avg_time{0};
  double avg_num_gates{0};
  int max_num_gates{0};

  Ntk ntk;
  Ntk ntkBest;
  int nBest = 100000;
  for( int it{0}; it<nIters; ++it )
  {

    std::clock_t start;
    double duration;
    start = std::clock();
    bool CEC;

    /* init */
    DecNet<TT, Ntk>  net;
    net.init( vTruths, vMasks );
    /* info on the outputs */
    Y = net.getTargets();
    assert( Y.size()==1 );

    net.setOSY( *net.getFuncP(Y[0]), *net.getMaskP(Y[0]) );

    /* info on the inputs */
    X = net.getPIs();
    std::vector<int> S;

    boolean_chain = {}; // clean the vector
    for( int i{0}; i < X.size(); ++i )  
    {
      boolean_chain.push_back( X[i] );
      S.push_back( i );
    }

    for( int i{1}; i < X.size(); ++i )  
    {
      for( int j{0}; j < i; ++j )  
      {
        boolean_chain.push_back(net.create_xor( X[i], X[j] ));
      }
    }


    vD={};
    vS = {};
    for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
    {
      vD.push_back( X );
      vS.push_back( S );
    }

    int best_reward{0};
    int old_reward{0};
    bool is_sat = false;
    bool is_stuck = false;
    int nNewNodes = boolean_chain.size();
    /* solve */
    std::uniform_int_distribution<int>  distr2(0, 100);  

    while( !is_sat && !is_stuck )
    {
//show_state( &net, &Y );
      is_sat = check_sat( &net, 0, boolean_chain.size()-nNewNodes );

      DecAnalyzer<TT, Ntk> checker( &net, &X, &Y, &vS, &vD );
      if( !is_sat  )
      {
        checker.check2();
        checker.check_topdec( boolean_chain );
        checker.check_str( );
        
        std::vector<action_t<TT>> RP = checker.get_remap();
//checker.print_actions( RP );
        best_reward=0;
        std::vector<int> moves; 
//printf("b \n" );

        for( int iRP{0}; iRP<RP.size(); ++iRP )
        {
          int OLD_RWD = 1;
          if( OLD_RWD == 1 )
          {
            if( RP[iRP].reward >= old_reward && RP[iRP].type != DecAct_t::STR_ ) // there was =
              moves.push_back( iRP );
          }
          else
          {
            if( RP[iRP].reward > best_reward )
            {
              moves = {iRP};
              best_reward = RP[iRP].reward;
            }
            else if( RP[iRP].reward == best_reward )
              moves.push_back( iRP );
          }
        }
        int TCONS = -1;
        if(TCONS < 100)
        {
          for( int iRP{0}; iRP<RP.size(); ++iRP )
          {
            if( RP[iRP].type == DecAct_t::STR_ )
            {
              if( distr2(generator) > TCONS )
                moves.push_back( iRP );
            }
          }
        }
//printf("a \n" );
        if( moves.size() == 0 )
        {
          is_stuck = true;
        }
        else
        {
          std::uniform_int_distribution<int>  distr(0, moves.size()-1);  
          int MOVE = moves[distr(generator)];
          old_reward = RP[MOVE].reward;
  //printf("MOVE %d\n", MOVE );
          remap( &net, RP[MOVE] );
          nNewNodes = checker.CountNodes(RP[MOVE].type);

          if( nNewNodes < 0 )
            nNewNodes = boolean_chain.size();

        }
      }
    }
//net.print_net();
    /* convert result */
    if( is_stuck )
    {
      printf( ANSI_COLOR_RED " EMPTY SET OF MOVES " ANSI_COLOR_RESET "\n" );
    }
    else
    {
      DecChsToGraph<TT, Ntk> conv( net ); 
  
      ntk = conv.convert();
      Y = net.getTargets();
      default_simulator<kitty::dynamic_truth_table> sim( (*net.getFuncP(Y[0])).num_vars() );
      const auto tt = simulate<kitty::dynamic_truth_table>( ntk, sim )[0];
      printf("**********************************\n");
      printf("num gates =%d\n", ntk.num_gates());
      for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
      {
        TT F = *net.getFuncP( Y[iTrg] );
        CEC = kitty::equal( tt, F );
        if( kitty::equal( tt, F ) )
          printf( ANSI_COLOR_GREEN " CEC SUCCESFUL " ANSI_COLOR_RESET "\n" );
        else
        {
          printf( ANSI_COLOR_RED " CEC FAILED " ANSI_COLOR_RESET "\n" );
          assert(0);
        }
      }

      duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
      avg_time += duration/nIters;
      avg_num_gates += 1.0*ntk.num_gates()/nIters;
      if( ntk.num_gates() > max_num_gates )
        max_num_gates = ntk.num_gates();

      printf("%5.5f\n", duration );
      printf("**********************************\n");
      if( (ntk.num_gates() < nBest) && CEC )
      {
        nBest = ntk.num_gates();
        ntkBest = ntk;
      }
    }
  }
  printf("SUMMARY:\n");
  printf("%10s %10s %10s %10s\n %10d %5.5f %d %5.5f\n", 
          "min G", "E[G]", "max G", "E[t]\n",
           nBest, avg_num_gates, max_num_gates, avg_time );
  return ntkBest;
}



template<class TT, class Ntk>
Ntk DecSolver<TT, Ntk>::ccg_relax( int nIters )
{
  printf("CCG-RELAX\n");
  //if( vTruths[0].num_vars() < 10 )
  //  PrintSpecs();

  std::mt19937       generator(5);

  double avg_time{0};
  double avg_num_gates{0};
  int max_num_gates{0};

  Ntk ntk;
  Ntk ntkBest;
  int nBest = 100000;
  int nInit;

  std::clock_t start_global;
  double duration_global;
  start_global = std::clock();

  //for( int it{0}; it<nIters; ++it )
  int it = 0;
  while( ( std::clock() - start_global ) < nIters*(double) CLOCKS_PER_SEC )
  {
    std::clock_t start;
    double duration;
    start = std::clock();
    bool CEC;

    /* init */
    DecNet<TT, Ntk>  net;
    net.init( vTruths, vMasks );
    /* info on the outputs */
    Y = net.getTargets();
    assert( Y.size()==1 );

    net.setOSY( *net.getFuncP(Y[0]), *net.getMaskP(Y[0]) );

    /* info on the inputs */
    X = net.getPIs();
    std::vector<int> S;

    boolean_chain = {}; // clean the vector
    for( int i{0}; i < X.size(); ++i )  
    {
      boolean_chain.push_back( X[i] );
      S.push_back( i );
    }

    vD={};
    vS = {};
    for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
    {
      vD.push_back( X );
      vS.push_back( S );
    }

    int best_reward{0};
    int old_reward{0};
    bool is_sat = false;
    bool is_stuck = false;
    int nNewNodes = boolean_chain.size();
    /* solve */
    std::uniform_int_distribution<int>  distr2(0, 100);  

    while( !is_sat && !is_stuck )
    {
//show_state( &net, &Y );
      is_sat = check_sat( &net, 0, boolean_chain.size()-nNewNodes );

      DecAnalyzer<TT, Ntk> checker( &net, &X, &Y, &vS, &vD );
      if( !is_sat  )
      {
        checker.check2();
        //checker.check_topdec( boolean_chain ); CCG-DEC
        //checker.check_str( ); CCG-SPEC
        
        std::vector<action_t<TT>> RP = checker.get_remap();
//checker.print_actions( RP );
        best_reward=0;
        std::vector<int> moves; 
//printf("b \n" );

        for( int iRP{0}; iRP<RP.size(); ++iRP )
        {
          int OLD_RWD = 1;
          if( OLD_RWD == 1 )
          {
            if( RP[iRP].reward > old_reward )
              moves.push_back( iRP );
          }
          else
          {
            if( RP[iRP].reward > best_reward )
            {
              moves = {iRP};
              best_reward = RP[iRP].reward;
            }
            else if( RP[iRP].reward == best_reward )
              moves.push_back( iRP );
          }
        }
      

        if( moves.size() == 0 )
        {
          is_stuck = true;
        }
        else
        {
          std::uniform_int_distribution<int>  distr(0, moves.size()-1);  
          int MOVE = moves[distr(generator)];//
          old_reward = RP[MOVE].reward;
          remap( &net, RP[MOVE] );
          nNewNodes = checker.CountNodes(RP[MOVE].type);

          if( nNewNodes < 0 )
            nNewNodes = boolean_chain.size();

        }
      }
    }
    /* convert result */
    if( is_stuck )
    {
      printf( ANSI_COLOR_RED " EMPTY SET OF MOVES " ANSI_COLOR_RESET "\n" );
    }
    else
    {
      DecChsToGraph<TT, Ntk> conv( net ); 
  
      ntk = conv.convert();
      Y = net.getTargets();
      default_simulator<kitty::dynamic_truth_table> sim( (*net.getFuncP(Y[0])).num_vars() );
      const auto tt = simulate<kitty::dynamic_truth_table>( ntk, sim )[0];
      printf("**********************************\n");
      printf("num gates =%d\n", ntk.num_gates());
      
      for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
      {
        TT F = *net.getFuncP( Y[iTrg] );
        CEC = kitty::equal( tt, F );
        if( kitty::equal( tt, F ) )
          printf( ANSI_COLOR_GREEN " CEC SUCCESFUL " ANSI_COLOR_RESET "\n" );
        else
        {
          printf( ANSI_COLOR_RED " CEC FAILED " ANSI_COLOR_RESET "\n" );
          assert(0);
        }
      }

      duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
      avg_time += duration/nIters;
      avg_num_gates += 1.0*ntk.num_gates()/nIters;
      if( ntk.num_gates() > max_num_gates )
        max_num_gates = ntk.num_gates();

      printf("%5.5f\n", duration );
      printf("**********************************\n");
      if( (ntk.num_gates() < nBest) && CEC )
      {
        nBest = ntk.num_gates();
        ntkBest = ntk;
      }
    }
    if( it == 0 )
      nInit = ntk.num_gates();
    it++;
  }
  duration_global = ( std::clock() - start_global ) / (double) CLOCKS_PER_SEC;

  printf("SUMMARY:\n");
  printf("%10s %10s %10s %10s %10s\n %10d %5.5f %d %5.5f %5.5f\n", 
          "min G", "E[G]", "max G", "E[t]\n" , "T[s]\n",
           nBest, avg_num_gates, max_num_gates, avg_time, duration_global );
  
  printf("%10s %10s %10s %10s\n %10d %10d %10d %10d\n", 
          "nIts", "dinit", "+dmin", "+dmax",
          it, nInit-nBest, nBest, max_num_gates-nBest );
  return ntkBest;
}

template<class TT, class Ntk>
Ntk DecSolver<TT, Ntk>::ccg_xor( int nIters )
{
  printf("CCG-XOR\n");
  //if( vTruths[0].num_vars() < 10 )
  //  PrintSpecs();

  std::mt19937       generator(5);

  double avg_time{0};
  double avg_num_gates{0};
  int max_num_gates{0};

  Ntk ntk;
  Ntk ntkBest;
  int nBest = 100000;
  int nInit;

  std::clock_t start_global;
  double duration_global;
  start_global = std::clock();

  //for( int it{0}; it<nIters; ++it )
  int it = 0;
  while((( std::clock() - start_global ) < nIters*(double) CLOCKS_PER_SEC ) && (nBest>6) )
  {
    std::clock_t start;
    double duration;
    start = std::clock();
    bool CEC;

    /* init */
    DecNet<TT, Ntk>  net;
    net.init( vTruths, vMasks );
    /* info on the outputs */
    Y = net.getTargets();
    assert( Y.size()==1 );

    net.setOSY( *net.getFuncP(Y[0]), *net.getMaskP(Y[0]) );

    /* info on the inputs */
    X = net.getPIs();
    std::vector<int> S;

    boolean_chain = {}; // clean the vector
    for( int i{0}; i < X.size(); ++i )  
    {
      boolean_chain.push_back( X[i] );
      S.push_back( i );
    }

    for( int i{1}; i < X.size(); ++i )  
    {
      for( int j{0}; j < i; ++j )  
      {
        boolean_chain.push_back(net.create_xor( X[i], X[j] ));
      }
    }

    vD={};
    vS = {};
    for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
    {
      vD.push_back( X );
      vS.push_back( S );
    }

    int best_reward{0};
    int old_reward{0};
    bool is_sat = false;
    bool is_stuck = false;
    int nNewNodes = boolean_chain.size();
    /* solve */
    std::uniform_int_distribution<int>  distr2(0, 100);  

    while( !is_sat && !is_stuck )
    {
//show_state( &net, &Y );
      is_sat = check_sat( &net, 0, boolean_chain.size()-nNewNodes );

      DecAnalyzer<TT, Ntk> checker( &net, &X, &Y, &vS, &vD );
      if( !is_sat  )
      {
        checker.check2();
        checker.check_topdec( boolean_chain ); //CCG-DEC
        //checker.check_str( ); CCG-SPEC
        
        std::vector<action_t<TT>> RP = checker.get_remap();
//checker.print_actions( RP );
        best_reward=0;
        std::vector<int> moves; 
//printf("b \n" );

        for( int iRP{0}; iRP<RP.size(); ++iRP )
        {
          int OLD_RWD = 1;
          if( OLD_RWD == 1 )
          {
            if( RP[iRP].reward > old_reward )
              moves.push_back( iRP );
          }
          else
          {
            if( RP[iRP].reward > best_reward )
            {
              moves = {iRP};
              best_reward = RP[iRP].reward;
            }
            else if( RP[iRP].reward == best_reward )
              moves.push_back( iRP );
          }
        }
      

        if( moves.size() == 0 )
        {
          is_stuck = true;
        }
        else
        {
          std::uniform_int_distribution<int>  distr(0, moves.size()-1);  
          int MOVE = moves[distr(generator)];//
          old_reward = RP[MOVE].reward;
          remap( &net, RP[MOVE] );
          nNewNodes = checker.CountNodes(RP[MOVE].type);

          if( nNewNodes < 0 )
            nNewNodes = boolean_chain.size();

        }
      }
    }
    /* convert result */
    if( is_stuck )
    {
      printf( ANSI_COLOR_RED " EMPTY SET OF MOVES " ANSI_COLOR_RESET "\n" );
    }
    else
    {
      net.print_net();
      DecChsToGraph<TT, Ntk> conv( net ); 
  
      ntk = conv.convert();
      Y = net.getTargets();
      default_simulator<kitty::dynamic_truth_table> sim( (*net.getFuncP(Y[0])).num_vars() );
      const auto tt = simulate<kitty::dynamic_truth_table>( ntk, sim )[0];
      printf("**********************************\n");
      printf("num gates =%d\n", ntk.num_gates());
      
      for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
      {
        TT F = *net.getFuncP( Y[iTrg] );
        CEC = kitty::equal( tt, F );
        if( kitty::equal( tt, F ) )
          printf( ANSI_COLOR_GREEN " CEC SUCCESFUL " ANSI_COLOR_RESET "\n" );
        else
        {
          printf( ANSI_COLOR_RED " CEC FAILED " ANSI_COLOR_RESET "\n" );
          assert(0);
        }
      }

      duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
      avg_time += duration/nIters;
      avg_num_gates += 1.0*ntk.num_gates()/nIters;
      if( ntk.num_gates() > max_num_gates )
        max_num_gates = ntk.num_gates();

      printf("%5.5f\n", duration );
      printf("**********************************\n");
      if( (ntk.num_gates() < nBest) && CEC )
      {
        nBest = ntk.num_gates();
        ntkBest = ntk;
      }
    }
    if( it == 0 )
      nInit = ntk.num_gates();
    it++;
  }
  duration_global = ( std::clock() - start_global ) / (double) CLOCKS_PER_SEC;

  printf("SUMMARY:\n");
  printf("%10s %10s %10s %10s %10s\n %10d %5.5f %d %5.5f %5.5f\n", 
          "min G", "E[G]", "max G", "E[t]\n" , "T[s]\n",
           nBest, avg_num_gates, max_num_gates, avg_time, duration_global );
  
  printf("%10s %10s %10s %10s\n %10d %10d %10d %10d\n", 
          "nIts", "dinit", "nmin", "+dmax",
          it, nInit-nBest, nBest, max_num_gates-nBest );
  return ntkBest;
}


template<class TT, class Ntk>
Ntk DecSolver<TT, Ntk>::ccgX( int nIters )
{
  printf("CCG-X\n");

  std::mt19937       generator(5);

  double avg_time{0};
  double avg_num_gates{0};
  int max_num_gates{0};

  Ntk ntk;
  Ntk ntkBest;
  int nBest = 100000;
  int nInit;

  std::clock_t start_global;
  double duration_global;
  start_global = std::clock();

  //for
  int it = 0;
  for( int it{0}; it<nIters; ++it )//(( std::clock() - start_global ) < nIters*(double) CLOCKS_PER_SEC ) && (nBest>6) )
  {
    std::clock_t start;
    double duration;
    start = std::clock();
    bool CEC;

    /* init */
    DecNet<TT, Ntk>  net;
    net.init( vTruths, vMasks );
    /* info on the outputs */
    Y = net.getTargets();
    assert( Y.size()==1 );

    net.setOSY( *net.getFuncP(Y[0]), *net.getMaskP(Y[0]) );

    /* info on the inputs */
    X = net.getPIs();
    std::vector<int> S;

    boolean_chain = {}; // clean the vector
    for( int i{0}; i < X.size(); ++i )  
    {
      boolean_chain.push_back( X[i] );
      S.push_back( i );
    }

    for( int i{1}; i < X.size(); ++i )  
    {
      for( int j{0}; j < i; ++j )  
      {
        boolean_chain.push_back(net.create_xor( X[i], X[j] ));
      }
    }

    vD={};
    vS = {};
    for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
    {
      vD.push_back( X );
      vS.push_back( S );
    }

    int best_reward{0};
    int old_reward{0};
    bool is_sat = false;
    bool is_stuck = false;
    int nNewNodes = boolean_chain.size();
    /* solve */
    std::uniform_int_distribution<int>  distr2(0, 100);  
    std::vector<action_t<TT>> CLOSE;
    while( CLOSE.size() == 0 && !is_stuck )
    {
//show_state( &net, &Y );
      DecAnalyzer<TT, Ntk> checker( &net, &X, &Y, &vS, &vD );
      checker.check_sat( boolean_chain, boolean_chain.size()-nNewNodes );
      checker.check_satX( boolean_chain, boolean_chain.size()-nNewNodes );
      CLOSE = checker.get_closure();
      is_sat = CLOSE.size() != 0;
      if( CLOSE.size() == 0  )
      {
        checker.check2();
        checker.check_topdec( boolean_chain ); //CCG-DEC
        //checker.check_str( ); CCG-SPEC
        
        std::vector<action_t<TT>> RP = checker.get_remap();
//checker.print_actions( RP );
        best_reward=0;
        std::vector<int> moves; 
//printf("b \n" );

        for( int iRP{0}; iRP<RP.size(); ++iRP )
        {
          int OLD_RWD = 0;
          if( OLD_RWD == 1 )
          {
            if( RP[iRP].reward > old_reward )
              moves.push_back( iRP );
          }
          else
          {
            if( RP[iRP].reward > best_reward )
            {
              moves = {iRP};
              best_reward = RP[iRP].reward;
            }
            else if( RP[iRP].reward == best_reward )
              moves.push_back( iRP );
          }
        }
      

        if( moves.size() == 0 )
        {
          is_stuck = true;
        }
        else
        {
          std::uniform_int_distribution<int>  distr(0, moves.size()-1);  
          int MOVE = moves[distr(generator)];//
          old_reward = RP[MOVE].reward;
          remap( &net, RP[MOVE] );
          nNewNodes = checker.CountNodes(RP[MOVE].type);

          if( nNewNodes < 0 )
            nNewNodes = boolean_chain.size();

        }
      }
    }
    /* convert result */
    if( is_stuck )
      printf( ANSI_COLOR_RED " EMPTY SET OF MOVES " ANSI_COLOR_RESET "\n" );
    else if(CLOSE.size() == 0)
      printf( ANSI_COLOR_RED " NO CLOSURE " ANSI_COLOR_RESET "\n" );
    else
    {
      Ntk ntk_tmp;
      Ntk ntk_best;
      int nBest_loc = 1000;
      for( int iClosure = 0; iClosure<CLOSE.size(); ++iClosure )
      {
        DecNet<TT, Ntk> net_tmp = net;
        close( &net_tmp, CLOSE[iClosure] );
        DecChsToGraph<TT, Ntk> conv( net_tmp ); 
        ntk_tmp = conv.convert();
        ntk_tmp = cleanup_dangling( ntk_tmp );
        //printf("this is %d\n", ntk_tmp.num_gates());
        if( ntk_tmp.num_gates() < nBest_loc && ntk_tmp.num_gates() != 0 )
        {
          ntk_best = ntk_tmp;
          nBest_loc = ntk_tmp.num_gates();
        }
      }
net.print_net();
      DecChsToGraph<TT, Ntk> conv( net ); 

      ntk=ntk_best;

      Y = net.getTargets();
      default_simulator<kitty::dynamic_truth_table> sim( (*net.getFuncP(Y[0])).num_vars() );
      const auto tt = simulate<kitty::dynamic_truth_table>( ntk, sim )[0];
      printf("**********************************\n");
      printf("num gates =%d\n", ntk.num_gates());
      
      for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
      {
        TT F = *net.getFuncP( Y[iTrg] );
        kitty::print_binary( F ); printf("\n");
        kitty::print_binary( tt ); printf("\n");

        CEC = kitty::equal( tt, F );
        if( kitty::equal( tt, F ) )
          printf( ANSI_COLOR_GREEN " CEC SUCCESFUL " ANSI_COLOR_RESET "\n" );
        else
        {
          printf( ANSI_COLOR_RED " CEC FAILED " ANSI_COLOR_RESET "\n" );
          assert(0);
        }
      }

      duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
      avg_time += duration/nIters;
      avg_num_gates += 1.0*ntk.num_gates()/nIters;
      if( ntk.num_gates() > max_num_gates )
        max_num_gates = ntk.num_gates();

      printf("%5.5f\n", duration );
      printf("**********************************\n");
      if( (ntk.num_gates() < nBest) && CEC )
      {
        nBest = ntk.num_gates();
        ntkBest = ntk;
      }
    }
    if( it == 0 )
      nInit = ntk.num_gates();
    it++;
  }
  duration_global = ( std::clock() - start_global ) / (double) CLOCKS_PER_SEC;

  printf("SUMMARY:\n");
  printf("%10s %10s %10s %10s %10s\n %10d %5.5f %d %5.5f %5.5f\n", 
          "min G", "E[G]", "max G", "E[t]\n" , "T[s]\n",
           nBest, avg_num_gates, max_num_gates, avg_time, duration_global );
  
  printf("%10s %10s %10s %10s\n %10d %10d %10d %10d\n", 
          "nIts", "dinit", "nmin", "+dmax",
          it, nInit-nBest, nBest, max_num_gates-nBest );
  return ntkBest;
}


template<class TT, class Ntk>
Ntk DecSolver<TT, Ntk>::ccg_spectral( int nIters, int PRC )
{
  printf("CCG-SPECTRAL\n");
  //if( vTruths[0].num_vars() < 10 )
  //  PrintSpecs();

  std::mt19937       generator(5);

  double avg_time{0};
  double avg_num_gates{0};
  int max_num_gates{0};

  Ntk ntk;
  Ntk ntkBest;
  int nBest = 100000;

  std::clock_t start_global;
  double duration_global;
  start_global = std::clock();
  int it=0;
  while( ( std::clock() - start_global ) < nIters*(double) CLOCKS_PER_SEC )//int it{0}; it<nIters; ++it )
  {
    
    std::clock_t start;
    double duration;
    start = std::clock();
    bool CEC;

    /* init */
    DecNet<TT, Ntk>  net;
    net.init( vTruths, vMasks );
    /* info on the outputs */
    Y = net.getTargets();
    assert( Y.size()==1 );

    net.setOSY( *net.getFuncP(Y[0]), *net.getMaskP(Y[0]) );

    /* info on the inputs */
    X = net.getPIs();
    std::vector<int> S;

    boolean_chain = {}; // clean the vector
    for( int i{0}; i < X.size(); ++i )  
    {
      boolean_chain.push_back( X[i] );
      S.push_back( i );
    }

    for( int i{1}; i < X.size(); ++i )  
    {
      for( int j{0}; j < i; ++j )  
      {
        boolean_chain.push_back(net.create_xor( X[i], X[j] ));
      }
    }


    vD={};
    vS = {};
    for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
    {
      vD.push_back( X );
      vS.push_back( S );
    }

    int best_reward{0};
    int old_reward{0};
    bool is_sat = false;
    bool is_stuck = false;
    int nNewNodes = boolean_chain.size();
    /* solve */
    std::uniform_int_distribution<int>  distr2(0, 100);  

    while( !is_sat && !is_stuck && (( std::clock() - start ) < 1*(double)CLOCKS_PER_SEC) )
    {
//show_state( &net, &Y );
      is_sat = check_sat( &net, 0, boolean_chain.size()-nNewNodes );

      DecAnalyzer<TT, Ntk> checker( &net, &X, &Y, &vS, &vD );
      if( !is_sat  )
      {
        checker.check2();
        checker.check_topdec( boolean_chain );
        checker.check_str( );
        
        std::vector<action_t<TT>> RP = checker.get_remap();
//checker.print_actions( RP );
        best_reward=0;
        std::vector<int> moves; 
//printf("b \n" );
        double nSTR=0;
        for( int iRP{0}; iRP<RP.size(); ++iRP )
        {
          if( RP[iRP].type == DecAct_t::STR_ )
            nSTR++;

          int OLD_RWD = 1;
          if( OLD_RWD == 1 )
          {
            if( RP[iRP].reward >= old_reward && RP[iRP].type != DecAct_t::STR_ ) // there was =
              moves.push_back( iRP );
          }
          else
          {
            if( RP[iRP].reward > best_reward )
            {
              moves = {iRP};
              best_reward = RP[iRP].reward;
            }
            else if( RP[iRP].reward == best_reward )
              moves.push_back( iRP );
          }
        }

        int TCONS = -1;
        double r = 0.05;
        TCONS = (int)PRC;//*1.*RP.size()/nSTR;
        if(TCONS < 100)
        {
          for( int iRP{0}; iRP<RP.size(); ++iRP )
          {
            if( RP[iRP].type == DecAct_t::STR_ )
            {
              if( distr2(generator) < TCONS )
                moves.push_back( iRP );
            }
          }
        }
//printf("a \n" );
        if( moves.size() == 0 )
        {
          is_stuck = true;
        }
        else
        {
          std::uniform_int_distribution<int>  distr(0, moves.size()-1);  
          int MOVE = moves[distr(generator)];
          old_reward = RP[MOVE].reward;
  //printf("MOVE %d\n", MOVE );
          remap( &net, RP[MOVE] );
          nNewNodes = checker.CountNodes(RP[MOVE].type);

          if( nNewNodes < 0 )
            nNewNodes = boolean_chain.size();

        }
      }
    }
//net.print_net();
    /* convert result */
    if( is_stuck )
      printf( ANSI_COLOR_RED " EMPTY SET OF MOVES " ANSI_COLOR_RESET "\n" );
    else if( ( std::clock() - start ) > 1*(double)CLOCKS_PER_SEC )
      printf( ANSI_COLOR_RED " OVERTIME " ANSI_COLOR_RESET "\n" );
    else
    {
      DecChsToGraph<TT, Ntk> conv( net ); 
  
      ntk = conv.convert();
      Y = net.getTargets();
      default_simulator<kitty::dynamic_truth_table> sim( (*net.getFuncP(Y[0])).num_vars() );
      const auto tt = simulate<kitty::dynamic_truth_table>( ntk, sim )[0];
      printf("**********************************\n");
      printf("num gates =%d\n", ntk.num_gates());
      for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
      {
        TT F = *net.getFuncP( Y[iTrg] );
        CEC = kitty::equal( tt, F );
        if( kitty::equal( tt, F ) )
          printf( ANSI_COLOR_GREEN " CEC SUCCESFUL " ANSI_COLOR_RESET "\n" );
        else
        {
          printf( ANSI_COLOR_RED " CEC FAILED " ANSI_COLOR_RESET "\n" );
          assert(0);
        }
      }

      duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
      avg_time += duration/nIters;
      avg_num_gates += 1.0*ntk.num_gates()/nIters;
      if( ntk.num_gates() > max_num_gates )
        max_num_gates = ntk.num_gates();

      printf("%5.5f\n", duration );
      printf("**********************************\n");
      if( (ntk.num_gates() < nBest) && CEC )
      {
        nBest = ntk.num_gates();
        ntkBest = ntk;
      }
    }
    it++;
  }
  duration_global = ( std::clock() - start_global ) / (double) CLOCKS_PER_SEC;

  printf("SUMMARY:\n");
  printf("%10s %10s %10s %10s %10s\n %10d %5.5f %d %5.5f %5.5f\n", 
          "min G", "E[G]", "max G", "E[t]", "T[s]\n",
           nBest, avg_num_gates, max_num_gates, avg_time, duration_global );
  printf("nIters %d", it);
  return ntkBest;
}



template<class TT, class Ntk>
Ntk DecSolver<TT, Ntk>::aut_symGT_solve(int nIters)
{

  std::mt19937       generator(rand_dev());


  Ntk ntk;
  Ntk ntkBest;
  int nBest = 100000;
  for( int it{0}; it<nIters; ++it )
  {
    std::clock_t start;
    double duration;
    start = std::clock();
    bool CEC;
    
    /* init */
    DecNet<TT, Ntk>  net;
    net.init( vTruths, vMasks );
    /* info on the outputs */
    Y = net.getTargets();
    assert( Y.size() == 1 );
    net.setOSY( *net.getFuncP(Y[0]), *net.getMaskP(Y[0]) );

    /* info on the inputs */
    X = net.getPIs();
    std::vector<int> S;
    for( int i{0}; i < X.size(); ++i )  
      S.push_back( i );
    vS = {};
    vD = {};
    for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
    {
      vD.push_back( X );
      vS.push_back( S );
    }

    int best_reward{0};
    /* solve */
    while( vS[0].size() > 1 )
    {
      DecAnalyzer<TT, Ntk> checker( &net, &X, &Y, &vS, &vD );
      /* check for divisors synthesizing the function */
      //checker.check_divclosure();
      //std::vector<action_t<TT>> CS = checker.get_divclosure();
      //if( CS.size() > 0 )
      //  close_div( &net, CS );

      /* if synthesis is not finished */
      if( vS[0].size() > 1 )
      {
        checker.check2();
        std::vector<action_t<TT>> RP = checker.get_remap();
        /* select uniformly at random among the ones with maximum |CDC| */
        //checker.print_actions( RP );

        best_reward=0;
        std::vector<int> moves; 
        std::vector<int> rewards; 
        for( int iRP{0}; iRP<RP.size(); ++iRP )
        {
          if( RP[iRP].reward > best_reward  )
          {
            moves.push_back( iRP );
            rewards.push_back( RP[iRP].reward );
          }
        }
        std::uniform_int_distribution<int>  distr(0, moves.size()-1);  
        int MOVE = moves[distr(generator)];
        best_reward = rewards[MOVE];
        //printf("MOVE %d\n", MOVE );
        remap( &net, RP[MOVE] );
      }
    }
    so_terminate( &net, 0 );
    /* convert result */
    DecChsToGraph<TT, Ntk> conv( net ); 
 
    ntk = conv.convert();
    net.print_net();
    Y = net.getTargets();
    default_simulator<kitty::dynamic_truth_table> sim( (*net.getFuncP(Y[0])).num_vars() );
    const auto tt = simulate<kitty::dynamic_truth_table>( ntk, sim )[0];
    printf("**********************************\n");
    printf("num gates =%d\n", ntk.num_gates());
    for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
    {
      TT F = vTruths[iTrg];
      
      //printf("\n");
     //kitty::print_binary(F);
      //printf("\n");
      //kitty::print_binary(tt);
      //printf("\n");
      CEC = kitty::equal( tt, F );

      if( kitty::equal( tt, F ) )
        printf( ANSI_COLOR_GREEN " CEC SUCCESFUL " ANSI_COLOR_RESET "\n" );
      else
      {
        printf( ANSI_COLOR_RED " CEC FAILED " ANSI_COLOR_RESET "\n" );
      }
    }

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    printf("%5.5f\n", duration );
    printf("**********************************\n");
    if( (ntk.num_gates() < nBest) && CEC )
        ntkBest = ntk;
  }
  return ntkBest;
}

template<class TT, class Ntk>
Ntk DecSolver<TT, Ntk>::aut_rdec_solve(int nIters)
{
  Ntk ntk;
  Ntk ntkBest;
  int nBest = 100000;
  for( int it{0}; it<nIters; ++it )
  {
    std::clock_t start;
    double duration;
    start = std::clock();
    bool CEC;

    /* init */
    DecNet<TT, Ntk>  net;
    net.init( vTruths, vMasks );
    /* info on the outputs */
    Y = net.getTargets();
    /* info on the inputs */
    X = net.getPIs();
    std::vector<int> S;
    for( int i{0}; i < X.size(); ++i )  
      S.push_back( i );
    vS = {};
    vD={};
    for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
    {
      vD.push_back( X );
      vS.push_back( S );
    }

    /* solve */
    while( Y.size() > 0 )
    {
      DecAnalyzer<TT, Ntk> checker( &net, &X, &Y, &vS, &vD );

      checker.check_tarclosure();
      std::vector<action_t<TT>> CT = checker.get_trgclosure();
      if( CT.size() > 0 )
        close_tar( &net, CT );

      checker.check_divclosure();
      std::vector<action_t<TT>> CD = checker.get_divclosure();
      if( CD.size() > 0 )
        close_div( &net, CD );

      if( Y.size() > 0 )
      {
        checker.check_2stepsdec();
        std::vector<action_t<TT>> SD = checker.get_2stepsdec();
        action_t<TT> act = select_uniformly( &SD );        
        decompose( &net, act );
      }
    }
    //net.print_net();

    Y = net.getTargets();
    /* convert result */
    DecChsToGraph<TT, Ntk> conv( net );  
    ntk = conv.convert();
    default_simulator<kitty::dynamic_truth_table> sim( (*net.getFuncP(Y[0])).num_vars() );
    printf("**********************************\n");
    printf("num gates =%d\n", ntk.num_gates());
    for( int iTrg{0}; iTrg<Y.size(); ++iTrg )
    {
      const auto tt = simulate<kitty::dynamic_truth_table>( ntk, sim )[iTrg];

      //kitty::print_binary(F[iTrg]);
      //printf("\n");
      //kitty::print_binary(tt);
      //printf("\n");

      TT F = *net.getFuncP( Y[iTrg] );
      CEC = kitty::equal(tt,F);
      if( kitty::equal( tt, F ) )
        printf( ANSI_COLOR_GREEN " CEC SUCCESFUL " ANSI_COLOR_RESET "\n"  );
      else
      {
        printf( ANSI_COLOR_RED " CEC FAILED " ANSI_COLOR_RESET "\n" );
        assert(0);
      }
    }
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    printf("%5.2f\n", duration );
    printf("**********************************\n");
    if( ntk.num_gates() < nBest )
    {
      nBest = ntk.num_gates();
        ntkBest = ntk;
    }
  }
  return ntkBest;
}


#pragma endregion solve

#pragma region closure
template<class TT, class Ntk>
void DecSolver<TT, Ntk>::close_div( DecNet<TT, Ntk> * pNet, std::vector<action_t<TT>> actions )
{
  std::vector<int> rm_targs;
  for( int i{0}; i < Y.size(); ++i )
    rm_targs.push_back(0);

  for( auto act : actions )
  {
    if( rm_targs[act.sigs[0]] == 0 )
    {
      int iTrg = act.sigs[0];
      int iS   = act.sigs[1];
      int iDiv = vS[iTrg][iS];

      if( act.type == DecAct_t::BUF_ )
        pNet->close_target( Y[iTrg], vD[iTrg][iDiv], 0 );
      else if( act.type == DecAct_t::INV_ )
        pNet->close_target( Y[iTrg], vD[iTrg][iDiv], 1 );
      rm_targs[act.sigs[0]] = 1;
    }
  }
  for( int i{Y.size()-1}; i >= 0; --i )
  {
    if( rm_targs[i] == 1 )
    {
      Y.erase( Y.begin() + i ); 
      vS.erase( vS.begin() + i );
    } 
  }
}

template<class TT, class Ntk>
void DecSolver<TT, Ntk>::so_terminate( DecNet<TT, Ntk> * pNet, int iTrg )
{
  assert( vS.size() > 0 );
  assert( vS[0].size() == 1 );

  int iDiv = vS[iTrg][0];
  TT DT = * pNet->getFuncP( X[iDiv] );
  TT DM = * pNet->getMaskP( X[iDiv] );
  TT FT = pNet->getFuncOSY( );
  TT FM = pNet->getMaskCDC( );


  if( kitty::equal( DT & FM, FT & FM ) )
    pNet->close_target( Y[iTrg], vD[iTrg][iDiv], 0 );
  else if( kitty::equal( ~DT & FM, FT & FM ) )
    pNet->close_target( Y[iTrg], vD[iTrg][iDiv], 1 );
  else
  {
    printf("BAD TERMINATION\n");
    assert(0);
  }

}

template<class TT, class Ntk>
bool DecSolver<TT, Ntk>::check_2sat( DecNet<TT, Ntk> * pNet, signal_t xi, signal_t xj )
{
  TT Tf = * pNet->getFuncP(Y[0]); // WARNING 0
  TT Mf = * pNet->getMaskP(Y[0]);

  TT Tdi = * pNet->getFuncP( xi );
  TT Tdj = * pNet->getFuncP( xj );

  if( kitty::equal( Mf&( Tdi & Tdj ), Mf&Tf ) ) // AND
  {
    printf( ANSI_COLOR_BLUE " %d AND %d" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    signal_t f0 = pNet->create_and( xi, xj );
    pNet->close_target( Y[0], f0, 0 );
    return true;
  }
  else if( kitty::equal( Mf&( ~(Tdi & Tdj) ), Mf&Tf ) ) // NAND
  {
    printf( ANSI_COLOR_BLUE " %d NAND %d" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    signal_t f0 = pNet->create_nand( xi, xj );
    pNet->close_target( Y[0], f0, 0 );
    return true;
  }
  else if( kitty::equal( Mf&( Tdi | Tdj ), Mf&Tf ) ) // NAND
  {
    printf( ANSI_COLOR_BLUE " %d OR %d" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    signal_t f0 = pNet->create_or( xi, xj );
    pNet->close_target( Y[0], f0, 0 );
    return true;
  }
  else if( kitty::equal( Mf&( ~(Tdi | Tdj) ), Mf&Tf ) ) // NAND
  {
    signal_t f0 = pNet->create_nor( xi, xj );
    printf( ANSI_COLOR_BLUE " %d NOR %d" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    pNet->close_target( Y[0], f0, 0 );
    return true;
  }
  else if( kitty::equal( Mf&( ~Tdi & Tdj ), Mf&Tf ) ) // NAND
  {
    printf( ANSI_COLOR_BLUE " %d LT %d" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    signal_t f0 = pNet->create_lt( xi, xj );
    pNet->close_target( Y[0], f0, 0 );
    return true;
  }
  else if( kitty::equal( Mf&( Tdi& ~Tdj), Mf&Tf ) )// NAND
  {
    printf( ANSI_COLOR_BLUE " %d AND %d" ANSI_COLOR_RESET "\n", xj.node, xi.node );
    signal_t f0 = pNet->create_lt( xj, xi );
    pNet->close_target( Y[0], f0, 0 );
    return true;
  }
  else if( kitty::equal( Mf&( ~Tdi | Tdj ), Mf&Tf ) ) // NAND
  {
    printf( ANSI_COLOR_BLUE " %d LE %d" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    signal_t f0 = pNet->create_le( xi, xj );
    pNet->close_target( Y[0], f0, 0 );
    return true;
  }
  else if( kitty::equal( Mf&( Tdi | ~Tdj), Mf&Tf ) )// NAND
  {
    printf( ANSI_COLOR_BLUE " %d LE %d" ANSI_COLOR_RESET "\n", xj.node, xi.node );
    signal_t f0 = pNet->create_le( xj, xi );
    pNet->close_target( Y[0], f0, 0 );
    return true;
  }
  else if( kitty::equal( Mf&( Tdi ^ Tdj ), Mf&Tf ) ) // NAND
  {
    printf( ANSI_COLOR_BLUE " %d XOR %d" ANSI_COLOR_RESET "\n", xi.node, xj.node );

    signal_t f0 = pNet->create_xor( xi, xj );
    pNet->close_target( Y[0], f0, 0 );
    return true;
  }
  else if( kitty::equal( Mf&( Tdi ^ ~Tdj), Mf&Tf ) )// NAND
  {
    printf( ANSI_COLOR_BLUE " %d XNOR %d" ANSI_COLOR_RESET "\n", xi.node, xj.node );

    signal_t f0 = pNet->create_xnor( xi, xj );
    pNet->close_target( Y[0], f0, 0 );
    return true;
  }
  else
  {
    return false;
  }
}

template<class TT, class Ntk>
bool DecSolver<TT, Ntk>::check_2satX( DecNet<TT, Ntk> * pNet, signal_t xi, signal_t xj, signal_t xd )
{
  TT Tf = * pNet->getFuncP(Y[0]); // WARNING 0
  TT Mf = * pNet->getMaskP(Y[0]);

  TT Tdi = * pNet->getFuncP( xi );
  TT Tdj = * pNet->getFuncP( xj );
  TT Tdd = * pNet->getFuncP( xd );

  if( kitty::equal( Mf&( Tdi & Tdj ), Mf&(Tf^Tdd) ) ) // AND
  {
    printf( ANSI_COLOR_BLUE "^ %d AND %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    signal_t f0 = pNet->create_and( xi, xj );
    signal_t fx = pNet->create_xor( f0, xd );
    pNet->close_target( Y[0], fx, 0 );
    return true;
  }
  else if( kitty::equal( Mf&( ~(Tdi & Tdj) ), Mf&(Tf^Tdd) ) ) // NAND
  {
    printf( ANSI_COLOR_BLUE "^%d NAND %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    signal_t f0 = pNet->create_nand( xi, xj );
    signal_t fx = pNet->create_xor( f0, xd );
    pNet->close_target( Y[0], fx, 0 );
    return true;
  }
  else if( kitty::equal( Mf&( Tdi | Tdj ), Mf&(Tf^Tdd) ) ) // NAND
  {
    printf( ANSI_COLOR_BLUE "^%d OR %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    signal_t f0 = pNet->create_or( xi, xj );
    signal_t fx = pNet->create_xor( f0, xd );
    pNet->close_target( Y[0], fx, 0 );
    return true;
  }
  else if( kitty::equal( Mf&( ~(Tdi | Tdj) ), Mf&(Tf^Tdd) ) ) // NAND
  {
    signal_t f0 = pNet->create_nor( xi, xj );
    printf( ANSI_COLOR_BLUE "^%d NOR %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    signal_t fx = pNet->create_xor( f0, xd );
    pNet->close_target( Y[0], fx, 0 );
    return true;
  }
  else if( kitty::equal( Mf&( ~Tdi & Tdj ), Mf&(Tf^Tdd) ) ) // NAND
  {
    printf( ANSI_COLOR_BLUE "^%d LT %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    signal_t f0 = pNet->create_lt( xi, xj );
    signal_t fx = pNet->create_xor( f0, xd );
    pNet->close_target( Y[0], fx, 0 );
    return true;
  }
  else if( kitty::equal( Mf&( Tdi& ~Tdj), Mf&(Tf^Tdd) ) )// NAND
  {
    printf( ANSI_COLOR_BLUE "^%d AND %d^" ANSI_COLOR_RESET "\n", xj.node, xi.node );
    signal_t f0 = pNet->create_lt( xj, xi );
    signal_t fx = pNet->create_xor( f0, xd );
    pNet->close_target( Y[0], fx, 0 );
    return true;
  }
  else if( kitty::equal( Mf&( ~Tdi | Tdj ), Mf&(Tf^Tdd) ) ) // NAND
  {
    printf( ANSI_COLOR_BLUE "^%d LE %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    signal_t f0 = pNet->create_le( xi, xj );
    signal_t fx = pNet->create_xor( f0, xd );
    pNet->close_target( Y[0], fx, 0 );
    return true;
  }
  else if( kitty::equal( Mf&( Tdi | ~Tdj), Mf&(Tf^Tdd) ) )// NAND
  {
    printf( ANSI_COLOR_BLUE "^%d LE %d^" ANSI_COLOR_RESET "\n", xj.node, xi.node );
    signal_t f0 = pNet->create_le( xj, xi );
    signal_t fx = pNet->create_xor( f0, xd );
    pNet->close_target( Y[0], fx, 0 );
    return true;
  }
  else if( kitty::equal( Mf&( Tdi ^ Tdj ), Mf&(Tf^Tdd) ) ) // NAND
  {
    printf( ANSI_COLOR_BLUE "^%d XOR %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );

    signal_t f0 = pNet->create_xor( xi, xj );
    signal_t fx = pNet->create_xor( f0, xd );
    pNet->close_target( Y[0], fx, 0 );
    return true;
  }
  else if( kitty::equal( Mf&( Tdi ^ ~Tdj), Mf&(Tf^Tdd) ) )// NAND
  {
    printf( ANSI_COLOR_BLUE "^%d XNOR %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );

    signal_t f0 = pNet->create_xnor( xi, xj );
    signal_t fx = pNet->create_xor( f0, xd );
    pNet->close_target( Y[0], fx, 0 );
    return true;
  }
  else
  {
    return false;
  }

}

template<class TT, class Ntk>
bool DecSolver<TT, Ntk>::check_sat( DecNet<TT, Ntk> * pNet, int iTrg, int newBound )
{
  TT Tf = * pNet->getFuncP(Y[iTrg]);
  TT Mf = * pNet->getMaskP(Y[iTrg]);
  int r = boolean_chain.size();
  for( int i{r-1}; i>=newBound; --i )
  {
    signal_t xi = boolean_chain[i];
    TT Tdi = * pNet->getFuncP( xi );
    //printf("CHECK SAT [%d %d]\n",xi.node, xi.sim);
    //kitty::print_binary( Td );
    //printf("\n");
    if( kitty::equal( Mf&Tdi, Mf&Tf ) )
    {
      pNet->close_target( Y[iTrg], xi, 0 );
      return true;
    }
    else if( kitty::equal( Mf&~Tdi, Mf&Tf ) )
    {
      pNet->close_target( Y[iTrg], xi, 1 );
      return true;
    }

    // second loop: the previous variables
    if( i > 0 )
    {
      for( int j{i-1}; j>=0; --j )
      {
        signal_t xj = boolean_chain[j];
        TT Tdj = * pNet->getFuncP( xj );
        if( check_2sat( pNet, xi, xj ) )
        {
          return true;
        }
      }
    }

  }
  return false;

}


template<class TT, class Ntk>
bool DecSolver<TT, Ntk>::check_satX( DecNet<TT, Ntk> * pNet, int iTrg, int newBound )
{
  TT Tf = * pNet->getFuncP(Y[iTrg]);
  TT Mf = (* pNet->getMaskP(Y[iTrg]));
  int r = boolean_chain.size();
  for( int d{0}; d<r; ++d )
  {
    signal_t xd = boolean_chain[d];
    for( int i{r-1}; i>=newBound; --i )
    {
      signal_t xi = boolean_chain[i];
      for( int j{i-1}; j>=0; --j )
      {
        signal_t xj = boolean_chain[j];
        if( check_2satX( pNet, xi, xj, xd ) )
          return true;
      }
    }
  }
  return false;

}



template<class TT, class Ntk>
void DecSolver<TT, Ntk>::close_tar( DecNet<TT, Ntk> * pNet, std::vector<action_t<TT>> actions )
{
  std::vector<int> rm_targs;
  for( int i{0}; i < Y.size(); ++i )
    rm_targs.push_back(0);

  for( auto act : actions )
  {
    if( rm_targs[act.sigs[1]] == 0 )
    {
      int iTrg1 = act.sigs[0];
      int iTrg2 = act.sigs[1];

      if( act.type == DecAct_t::BUF_ )
        pNet->close_target( Y[iTrg2], Y[iTrg1], 0 );
      else if( act.type == DecAct_t::INV_ )
        pNet->close_target( Y[iTrg2], Y[iTrg1], 1 );
      rm_targs[iTrg2] = 1;
      rm_targs[iTrg1] = 2;
      
    }
  }
  for( int i{Y.size()-1}; i >= 0; --i )
  {
    if( rm_targs[i] == 1 )
    {
      Y.erase( Y.begin() + i );  
      vS.erase( vS.begin() + i );  
    }
  }
}
#pragma endregion closure

#pragma region remap
template <class TT, class Ntk>
void DecSolver<TT, Ntk>::remap_nes(  DecNet<TT, Ntk> * pNet, action_t<TT> act  )
{
  int iTrg    = act.sigs[0];
  int i       = act.sigs[1];
  int j       = act.sigs[2];
  int iS      = vS[iTrg][i];  
  int jS      = vS[iTrg][j];
  if( act.id_ord == 0 )     
  { 
    signal_t ri = pNet->create_and( vD[iTrg][jS], vD[iTrg][iS] );
    signal_t rj = pNet->create_or ( vD[iTrg][jS], vD[iTrg][iS] );  
    vD[iTrg][jS] = rj;
    vD[iTrg][iS] = ri;
  }
  else if( act.id_ord == 1 )
  { 
    signal_t ri = pNet->create_or ( vD[iTrg][jS], vD[iTrg][iS] ); 
    signal_t rj = pNet->create_and( vD[iTrg][jS], vD[iTrg][iS] );  
    vD[iTrg][jS] = rj;
    vD[iTrg][iS] = ri;
  }
  else  
    std::cerr << "id_ord not valid_ord" << std::endl;

  boolean_chain.push_back( vD[iTrg][iS] );
  boolean_chain.push_back( vD[iTrg][jS] );
  
  pNet->setOSY( act.func[0], act.mask[0] );

}

template <class TT, class Ntk>
void DecSolver<TT, Ntk>::remap_es(  DecNet<TT, Ntk> * pNet, action_t<TT> act  )
{
  int iTrg    = act.sigs[0];
  int i       = act.sigs[1];
  int j       = act.sigs[2];
  int iS      = vS[iTrg][i];  
  int jS      = vS[iTrg][j];
  if( act.id_ord == 0 )     
  { 
    signal_t ri = pNet->create_le( vD[iTrg][jS], vD[iTrg][iS] );
    signal_t rj = pNet->create_le( vD[iTrg][iS], vD[iTrg][jS] );   
    vD[iTrg][jS] = rj;
    vD[iTrg][iS] = ri;
  }
  else if( act.id_ord == 1 )
  { 
    signal_t ri = pNet->create_lt( vD[iTrg][jS], vD[iTrg][iS] );  
    signal_t rj = pNet->create_lt( vD[iTrg][iS], vD[iTrg][jS] );  

    vD[iTrg][jS] = rj;
    vD[iTrg][iS] = ri;
  }
  else  
    std::cerr << "id_ord not valid_ord" << std::endl;
  
  boolean_chain.push_back( vD[iTrg][iS] );
  boolean_chain.push_back( vD[iTrg][jS] );

  pNet->setOSY( act.func[0], act.mask[0] );


}

template<class TT, class Ntk>
void DecSolver<TT, Ntk>::remap_svs(  DecNet<TT, Ntk> * pNet, action_t<TT> act  )
{
  int iTrg    = act.sigs[0];
  int i       = act.sigs[1];
  int j       = act.sigs[2];
  int iS      = vS[iTrg][i];  
  int jS      = vS[iTrg][j];
  if( act.id_sym == 0 ) // SVS0X_ { SVS Xj } Xi'
  {
    if( act.id_ord == 0 )
      vD[iTrg][jS] = pNet->create_le( vD[iTrg][iS], vD[iTrg][jS] );  // ( Xi & Xj' )'
    else
      vD[iTrg][jS] = pNet->create_and( vD[iTrg][jS], vD[iTrg][iS] ); // ( Xi & Xj )
    
    boolean_chain.push_back( vD[iTrg][jS] );
  }
  else if( act.id_sym == 1 ) // SVS1X  { SVS Xj } Xi
  {
    if( act.id_ord == 0 )
      vD[iTrg][jS] = pNet->create_or( vD[iTrg][jS], vD[iTrg][iS] ); 
    else
      vD[iTrg][jS] = pNet->create_lt( vD[iTrg][iS], vD[iTrg][jS] );

    boolean_chain.push_back( vD[iTrg][jS] );
  }
  else if( act.id_sym == 2 ) // SVSX0_  { SVS Xi } Xj'
  {
    if( act.id_ord == 0 )
      vD[iTrg][iS] = pNet->create_le( vD[iTrg][jS], vD[iTrg][iS] );
    else
      vD[iTrg][iS] = pNet->create_and( vD[iTrg][jS], vD[iTrg][iS] );

    boolean_chain.push_back( vD[iTrg][iS] );
  }
  else if( act.id_sym == 3 ) //SVSX1_  { SVS Xi } Xj
  {
    if( act.id_ord == 0 )
      vD[iTrg][iS] = pNet->create_or( vD[iTrg][jS], vD[iTrg][iS] );
    else
      vD[iTrg][iS] = pNet->create_lt( vD[iTrg][jS], vD[iTrg][iS] );

    boolean_chain.push_back( vD[iTrg][iS] );
  }
  else
    std::cerr << "wrong symmetry identifier for SVS" ;

  pNet->setOSY( act.func[0], act.mask[0] );

}

template<class TT, class Ntk>
void DecSolver<TT, Ntk>::remap_ms( DecNet<TT, Ntk> * pNet, action_t<TT> act )
{
  int iTrg    = act.sigs[0];
  int i       = act.sigs[1];
  int j       = act.sigs[2];
  int iS      = vS[iTrg][i];  
  int jS      = vS[iTrg][j];
  if( act.id_ord == 0 )
  {
    vD[iTrg][iS] = pNet->create_xnor( vD[iTrg][jS], vD[iTrg][iS] ) ;
    vS[iTrg].erase( vS[iTrg].begin() + j );
    boolean_chain.push_back( vD[iTrg][iS] );
  }
  else if( act.id_ord == 1 )
  {
    vD[iTrg][jS] = pNet->create_xor( vD[iTrg][jS], vD[iTrg][iS] );
    vS[iTrg].erase( vS[iTrg].begin() + i );
    boolean_chain.push_back( vD[iTrg][jS] );
  }        
  else if( act.id_ord == 2 )
  {
    vD[iTrg][jS] = pNet->create_xnor( vD[iTrg][jS], vD[iTrg][iS] );
    vS[iTrg].erase( vS[iTrg].begin() + i );
    boolean_chain.push_back( vD[iTrg][jS] );
  }        
  else if( act.id_ord == 3 )
  {
    vD[iTrg][iS] = pNet->create_xor( vD[iTrg][jS], vD[iTrg][iS] );
    vS[iTrg].erase( vS[iTrg].begin() + j );
    boolean_chain.push_back( vD[iTrg][iS] );
  }  
  
  pNet->setOSY( act.func[0], act.mask[0] );

}

template< class TT, class Ntk >
void DecSolver<TT, Ntk>::remap_csvs( DecNet<TT, Ntk> * pNet, action_t<TT> act )
{
  int iTrg    = act.sigs[0];
  int i       = act.sigs[1];
  int j       = act.sigs[2];
  int iS      = vS[iTrg][i];  
  int jS      = vS[iTrg][j];

  if( act.id_sym == 0 )
  {
    if( act.id_ord == 0 )
    {
      vD[iTrg][iS] = pNet->create_and( vD[iTrg][jS], vD[iTrg][iS] );
      vS[iTrg].erase( vS[iTrg].begin() + j );  
    }
    else
    {
      vD[iTrg][jS] = pNet->create_and( vD[iTrg][jS], vD[iTrg][iS] );
      vS[iTrg].erase( vS[iTrg].begin() + i );
    }
  }
  else if( act.id_sym == 1 )
  {
    if( act.id_ord == 0 )
    {
      vD[iTrg][iS] = pNet->create_le( vD[iTrg][jS], vD[iTrg][iS] );
      vS[iTrg].erase( vS[iTrg].begin() + j );
    }
    else
    {
      vD[iTrg][jS] = pNet->create_lt( vD[iTrg][iS], vD[iTrg][jS] );
      vS[iTrg].erase( vS[iTrg].begin() + i );
    }
  }
  else if( act.id_sym == 2 )
  {
    if( act.id_ord == 0 )
    {
      vD[iTrg][jS] = pNet->create_le( vD[iTrg][iS], vD[iTrg][jS] );
      vS[iTrg].erase( vS[iTrg].begin() + i );
    }
    else
    {
      vD[iTrg][iS] = pNet->create_lt( vD[iTrg][jS], vD[iTrg][iS] );
      vS[iTrg].erase( vS[iTrg].begin() + j );
    }
  }
  else if( act.id_sym == 3 )
  {
    if( act.id_ord == 0 )
    {
      vD[iTrg][jS] =  pNet->create_or( vD[iTrg][jS], vD[iTrg][iS] );
      vS[iTrg].erase( vS[iTrg].begin() + i );
    }
    else
    {
      vD[iTrg][iS] =  pNet->create_or( vD[iTrg][jS], vD[iTrg][iS] );
      vS[iTrg].erase( vS[iTrg].begin() + j );
    }
  }
  else
    std::cerr << "wrong symmetry identifier for CSVS" ;

  boolean_chain.push_back( vD[iTrg][iS] );
  boolean_chain.push_back( vD[iTrg][jS] );

  pNet->setOSY( act.func[0], act.mask[0] );

}

template<class TT, class Ntk>
void DecSolver<TT, Ntk>::remap_and( DecNet<TT, Ntk> * pNet, action_t<TT> act )
{
  printf( ANSI_COLOR_MAGENTA "TOP AND" ANSI_COLOR_RESET "\n" );
  int iTrg = act.sigs[0];
  assert( iTrg == 0 );
  int iDiv = act.sigs[1];
  signal_t xi = boolean_chain[iDiv];
  signal_t targ1 = pNet->create_target( act.func[0], act.mask[0] );
  signal_t targ_source = pNet->create_and( xi, targ1 );
  pNet->close_target( Y[iTrg], targ_source, 0 );
  Y[iTrg] = targ1;
  boolean_chain.erase( boolean_chain.begin() + iDiv );

}

template<class TT, class Ntk>
void DecSolver<TT, Ntk>::remap_or( DecNet<TT, Ntk> * pNet, action_t<TT> act )
{
  printf( ANSI_COLOR_MAGENTA "TOP OR" ANSI_COLOR_RESET "\n" );
  //TT Mf = pNet->getMaskCDC();

  int iTrg = act.sigs[0];
  assert( iTrg == 0 );
  int iDiv = act.sigs[1];
  signal_t xi = boolean_chain[iDiv];
  signal_t targ1 = pNet->create_target( act.func[0], act.mask[0] );
  signal_t targ_source = pNet->create_or( xi, targ1 );
  pNet->close_target( Y[iTrg], targ_source, 0 );
  Y[iTrg] = targ1;
  boolean_chain.erase( boolean_chain.begin() + iDiv );
}

template<class TT, class Ntk>
void DecSolver<TT, Ntk>::remap_lt( DecNet<TT, Ntk> * pNet, action_t<TT> act )
{
  printf( ANSI_COLOR_MAGENTA "TOP LT" ANSI_COLOR_RESET "\n" );

  int iTrg = act.sigs[0];
  assert( iTrg == 0 );
  int iDiv = act.sigs[1];
  signal_t xi = boolean_chain[iDiv];
  signal_t targ1 = pNet->create_target( act.func[0], act.mask[0] );
  signal_t targ_source = pNet->create_lt( xi, targ1 );
  pNet->close_target( Y[iTrg], targ_source, 0 );
  Y[iTrg] = targ1;
  boolean_chain.erase( boolean_chain.begin() + iDiv );

  //TT Mf = pNet->getMaskCDC();

}

template<class TT, class Ntk>
void DecSolver<TT, Ntk>::remap_le( DecNet<TT, Ntk> * pNet, action_t<TT> act )
{
  printf( ANSI_COLOR_MAGENTA "TOP LE" ANSI_COLOR_RESET "\n" );

  int iTrg = act.sigs[0];
  assert( iTrg == 0 );
  int iDiv = act.sigs[1];
  signal_t xi = boolean_chain[iDiv];
  signal_t targ1 = pNet->create_target( act.func[0], act.mask[0] );
  signal_t targ_source = pNet->create_le( xi, targ1 );
  pNet->close_target( Y[iTrg], targ_source, 0 );
  Y[iTrg] = targ1;
  boolean_chain.erase( boolean_chain.begin() + iDiv );

  //TT Mf = pNet->getMaskCDC();

}

template <class TT, class Ntk>
void DecSolver<TT, Ntk>::remap_str(  DecNet<TT, Ntk> * pNet, action_t<TT> act  )
{
//printf("@@e\n");
  printf( ANSI_COLOR_YELLOW "STR" ANSI_COLOR_RESET "\n" );

  int iTrg    = act.sigs[0];
  int i       = act.sigs[1];
  int j       = act.sigs[2];
  int iS      = vS[iTrg][i];  
  int jS      = vS[iTrg][j];
  
  vD[iTrg][iS] = pNet->create_xor( vD[iTrg][jS], vD[iTrg][iS] );

  boolean_chain.push_back( vD[iTrg][iS] );

  pNet->setOSY( act.func[0], act.mask[0] );

//printf("@@f\n");

}

template<class TT, class Ntk>
 void DecSolver<TT, Ntk>::remap( DecNet<TT, Ntk> * pNet, action_t<TT> act )
  {
    int iTrg    = act.sigs[0];

    switch( act.type )
    {
      case DecAct_t::NES_:
        remap_nes( pNet, act );
        break; 
      case DecAct_t::ES_:
        remap_es( pNet, act );
        break; 
      case DecAct_t::SVS_:
        remap_svs( pNet, act );
        break;
      case DecAct_t::MS_:
        remap_ms( pNet, act );
        break; 
      case DecAct_t::CSVS_:
        remap_csvs( pNet, act );
        break;
      case DecAct_t::T1_AND_:
        remap_and( pNet, act );
        break; 
      case DecAct_t::T1_OR_:
        remap_or( pNet, act );
        break; 
      case DecAct_t::T1_LT_:
        remap_lt( pNet, act );
        break; 
      case DecAct_t::T1_LE_:
        remap_le( pNet, act );
        break;
      case DecAct_t::STR_:
        remap_str( pNet, act );
        break; 
    }
    //pNet->change_sim_info( Y[iTrg], act.func[0], act.mask[0] );
  }

template<class TT, class Ntk>
void DecSolver<TT, Ntk>::close( DecNet<TT, Ntk> * pNet, action_t<TT> act )
{
  int iTrg    = act.sigs[0];
  int i       = act.sigs[1];
  int j       = act.sigs[2];
  int d       = act.sigs.size() > 3 ? act.sigs[3] : 0;
  signal_t xi = boolean_chain[i];
  signal_t xj = boolean_chain[j];
  signal_t xd = boolean_chain[d];

  TT Tf = * pNet->getFuncP(Y[0]); // WARNING 0
  TT Mf = * pNet->getMaskP(Y[0]);

  TT Tdi = * pNet->getFuncP( xi );
  TT Tdj = * pNet->getFuncP( xj );
  TT Tdd = * pNet->getFuncP( xd );

  signal_t f0 {0,0};
  signal_t fx {0,0};

  switch ( act.type )
  {
  case DecAct_t::XCL_AND_:
    //printf( ANSI_COLOR_BLUE "^ %d AND %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    f0 = pNet->create_and( xi, xj );
    fx = pNet->create_xor( f0, xd );
    pNet->close_target( Y[0], fx, 0 );
    break;
  case DecAct_t::XCL_NAND_:
    //printf( ANSI_COLOR_BLUE "^%d NAND %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    f0 = pNet->create_nand( xi, xj );
    fx = pNet->create_xor( f0, xd );
    pNet->close_target( Y[0], fx, 0 );
    break;
  case DecAct_t::XCL_OR_:
    //printf( ANSI_COLOR_BLUE "^%d OR %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    f0 = pNet->create_or( xi, xj );
    fx = pNet->create_xor( f0, xd );
    pNet->close_target( Y[0], fx, 0 );
    break;
  case DecAct_t::XCL_NOR_:
    f0 = pNet->create_nor( xi, xj );
    //printf( ANSI_COLOR_BLUE "^%d NOR %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    fx = pNet->create_xor( f0, xd );
    pNet->close_target( Y[0], fx, 0 );
    break;
  case DecAct_t::XCL_LT_:
    //printf( ANSI_COLOR_BLUE "^%d LT %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    f0 = pNet->create_lt( xi, xj );
    fx = pNet->create_xor( f0, xd );
    pNet->close_target( Y[0], fx, 0 );
    break;
  case DecAct_t::XCL_LE_:
    //printf( ANSI_COLOR_BLUE "^%d LE %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    f0 = pNet->create_le( xi, xj );
    fx = pNet->create_xor( f0, xd );
    pNet->close_target( Y[0], fx, 0 );
    break;
  case DecAct_t::XCL_XOR_:
    //printf( ANSI_COLOR_BLUE "^%d XOR %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    f0 = pNet->create_xor( xi, xj );
    fx = pNet->create_xor( f0, xd );
    pNet->close_target( Y[0], fx, 0 );
    break;
  case DecAct_t::XCL_XNOR_:
    //printf( ANSI_COLOR_BLUE "^%d XNOR %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    f0 = pNet->create_xnor( xi, xj );
    fx = pNet->create_xor( f0, xd );
    pNet->close_target( Y[0], fx, 0 );
    break;
  case DecAct_t::CL_AND_:
    //printf( ANSI_COLOR_BLUE "^ %d AND %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    f0 = pNet->create_and( xi, xj );
    pNet->close_target( Y[0], f0, 0 );
    break;
  case DecAct_t::CL_NAND_:
    //printf( ANSI_COLOR_BLUE "^%d NAND %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    f0 = pNet->create_nand( xi, xj );
    pNet->close_target( Y[0], f0, 0 );
    break;
  case DecAct_t::CL_OR_:
    //printf( ANSI_COLOR_BLUE "^%d OR %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    f0 = pNet->create_or( xi, xj );
    pNet->close_target( Y[0], f0, 0 );
    break;
  case DecAct_t::CL_NOR_:
    f0 = pNet->create_nor( xi, xj );
    //printf( ANSI_COLOR_BLUE "^%d NOR %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    pNet->close_target( Y[0], f0, 0 );
    break;
  case DecAct_t::CL_LT_:
    //printf( ANSI_COLOR_BLUE "^%d LT %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    f0 = pNet->create_lt( xi, xj );
    pNet->close_target( Y[0], f0, 0 );
    break;
  case DecAct_t::CL_LE_:
    //printf( ANSI_COLOR_BLUE "^%d LE %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    f0 = pNet->create_le( xi, xj );
    pNet->close_target( Y[0], f0, 0 );
    break;
  case DecAct_t::CL_XOR_:
    //printf( ANSI_COLOR_BLUE "^%d XOR %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    f0 = pNet->create_xor( xi, xj );
    pNet->close_target( Y[0], f0, 0 );
    break;
  case DecAct_t::CL_XNOR_:
    //printf( ANSI_COLOR_BLUE "^%d XNOR %d^" ANSI_COLOR_RESET "\n", xi.node, xj.node );
    f0 = pNet->create_xnor( xi, xj );
    pNet->close_target( Y[0], f0, 0 );
    break;
  default:
    break;
  }

}

#pragma endregion remap

#pragma region decompose
template<class TT, class Ntk>
 void DecSolver<TT, Ntk>::decompose( DecNet<TT, Ntk> * pNet, action_t<TT> act )
  {
    //pNet->change_sim_info( Y[act.sigs[0]], act.func[0], act.mask[0] );
    TT * pT1 = &act.func[1];
    TT * pT0 = &act.func[0];
    TT * pM1 = &act.mask[1];
    TT * pM0 = &act.mask[0];

    bool isAND = kitty::is_const0(   *pT0  & *pM0 );
    bool isLT  = kitty::is_const0(   *pT1  & *pM1 );
    bool isOR  = kitty::is_const0( ~(*pT1) & *pM1 );
    bool isLE  = kitty::is_const0( ~(*pT0) & *pM0 );

    int iTrg = act.sigs[0];  
    int iS = act.sigs[1];
    int iDiv = vS[iTrg][iS];
    signal_t x = vD[iTrg][iDiv];

    if( isAND )
    {
      signal_t targ1 = pNet->create_target( act.func[1], act.mask[1] );
      signal_t targ_source = pNet->create_and( x, targ1 );
      pNet->close_target( Y[iTrg], targ_source, 0 );
      Y[iTrg] = targ1;
      vS[iTrg].erase( vS[iTrg].begin() + iS );
    }
    else if( isLT )
    {
      signal_t targ0 = pNet->create_target( act.func[0], act.mask[0] );
      signal_t targ_source = pNet->create_lt( x, targ0 );
      pNet->close_target( Y[iTrg], targ_source, 0 );
      Y[iTrg] = targ0;
      vS[iTrg].erase( vS[iTrg].begin() + iS );
    }
    else if( isOR )
    {
      signal_t targ0 = pNet->create_target( act.func[0], act.mask[0] );
      signal_t targ_source = pNet->create_or( x, targ0 );
      pNet->close_target( Y[iTrg], targ_source, 0 );
      Y[iTrg] = targ0;
      vS[iTrg].erase( vS[iTrg].begin() + iS );
    }
    else if( isLE )
    {
      signal_t targ1 = pNet->create_target( act.func[1], act.mask[1] );
      signal_t targ_source = pNet->create_le( x, targ1 );
      pNet->close_target( Y[iTrg], targ_source, 0 );
      Y[iTrg] = targ1;
      vS[iTrg].erase( vS[iTrg].begin() + iS );
    }
    else
    {
      signal_t targ0 = pNet->create_target( act.func[0], act.mask[0] );
      signal_t targ1 = pNet->create_target( act.func[1], act.mask[1] );
      signal_t targ_replacement  = pNet->create_or( targ0, targ1 );
      pNet->close_target( Y[iTrg], targ_replacement, 0 );

      signal_t deeptarg0 = pNet->create_target( act.func[0], act.mask[0] );
      signal_t deeptarg1 = pNet->create_target( act.func[1], act.mask[1] );
      signal_t targ_replacement0 = pNet->create_lt( x, deeptarg0 );
      signal_t targ_replacement1 = pNet->create_and( x, deeptarg1 );
      pNet->close_target( targ0, targ_replacement0, 0 );
      pNet->close_target( targ1, targ_replacement1, 0 );


      Y[iTrg] = deeptarg0;
      vS[iTrg].erase( vS[iTrg].begin() + iS );
      vD.push_back( vD[iTrg] );
      vS.push_back( vS[iTrg] );
      Y.push_back(deeptarg1);
    }

  }
#pragma endregion decompose

#pragma region visualize
template<class TT, class Ntk>
void DecSolver<TT, Ntk>::PrintSpecs()
{
    printf("TRUTHS:\n");
    for( int i{0}; i<vTruths.size(); ++i ){ printf("%d ", i ); kitty::print_binary( vTruths[i] ); printf("\n");}
    printf("MASKS:\n");
    for( int i{0}; i<vMasks.size(); ++i ){ printf("%d ", i ); kitty::print_binary( vMasks[i] ); printf("\n");}
}

template<class TT, class Ntk>
void DecSolver<TT, Ntk>::show_state( DecNet<TT, Ntk> * pNet, std::vector<signal_t> * pY )
{
  for( int i{0}; i<(*pY).size(); ++i )
  {
    printf( ANSI_COLOR_YELLOW " TARGET #%d" ANSI_COLOR_RESET "", i );
    printf( ANSI_COLOR_YELLOW );
    TT F = pNet->getFuncOSY(  );
    TT M = pNet->getMaskCDC(  );
    printf("\n");
    kitty::print_binary( F );
    printf("\n");
    auto km_tt = kitty::karnaugh_map( F );
    km_tt.print( M );
    printf( ANSI_COLOR_RESET );
  }
}
#pragma endregion visualize


} // namespace mockturtle