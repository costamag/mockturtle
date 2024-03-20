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
  \file spfd_utils.hpp
  \brief Utilities for spfd-manipulation

  \author Andrea Costamagna
*/

#pragma once

#include "index_list.hpp"
#include <kitty/constructors.hpp>
#include <kitty/operators.hpp>
#include <kitty/print.hpp>
#include <bit>

namespace mockturtle
{

  template<uint32_t NUM_FANINS, uint32_t MAX_WIDTH>
  struct lut_resynthesis_t
  {
    using truth_table_t = kitty::dynamic_truth_table;
    lut_resynthesis_t()
    {}

    std::vector<uint32_t> find_support( truth_table_t const& tt, truth_table_t const& mk )
    {
      std::vector<uint32_t> supp;
      for( int i{0}; i<sim_target.num_vars(); ++i )
      {
        auto tt1 = kitty::cofactor1( tt, i );
        auto tt0 = kitty::cofactor0( tt, i );
        auto mk1 = kitty::cofactor1( mk, i );
        auto mk0 = kitty::cofactor0( mk, i );
        if( !kitty::equal( tt1&mk0&mk1, tt0&mk0&mk1 ) )
        {
          supp.push_back( i );
        }
      }
      return supp;
    }

    std::tuple<truth_table_t, truth_table_t> extract_lut( std::vector<uint32_t> const& cut, truth_table_t const& tt, truth_table_t const& mk )
    {
      assert( cut.size() <= NUM_FANINS && "Netlist size exceeds maximum fanin size" );
      truth_table_t lut(cut.size());
      auto sim = tt.construct();
      auto tmp = tt.construct();
      std::uniform_real_distribution<> U01(0, 1);

      for( uint32_t m{0}; m < ( 1u << cut.size() ); ++m )
      {
        tmp = tmp | ~tmp;
        for( int i{0}; i<cut.size(); ++i )
        {
          if( ( ( m >> i ) & 0x1 ) == 0x1 )
          {
            tmp &= sims[cut[i]];
          }
          else
          {
            tmp &= ~sims[cut[i]];
          }  
        }
        int n0 = kitty::count_ones( ~tt & mk & tmp );
        int n1 = kitty::count_ones(  tt & mk & tmp );
        if( n0 > n1 )
        {
          kitty::clear_bit( lut, m );
        }
        else if( n0 < n1 )
        {
          kitty::set_bit( lut, m );
          sim |= tmp;
        }
        else
        {
          if( U01(RNG) < 0.5 )
          {
            kitty::clear_bit( lut, m );
          }
          else
          {
            kitty::set_bit( lut, m );
            sim |= tmp;
          }
        }
      }
      return std::make_tuple( lut, sim );
    } 

    struct comb_t
    {
      comb_t( int N, int K ): N(N), K(K)
      {
      }

      uint32_t popcount( uint32_t m )
      {
          m = m - ((m >> 1) & 0x55555555);                 // add pairs of bits
          m = (m & 0x33333333) + ((m >> 2) & 0x33333333);  // quads
          m = (m + (m >> 4)) & 0x0F0F0F0F;                 // groups of 8
          m *= 0x01010101;                                 // horizontal sum of bytes

          return m >> 24;
      }

      std::optional<std::vector<uint32_t>> get()
      {
        std::vector<uint32_t> res;
        for( uint32_t i{state}; i < (1u << N); ++i )
        {
          if( popcount(i) == K )
          {
            for( int j{0}; j<N; ++j )
            {
              if( ( i >> j ) & 0x1 == 0x1 )
              {
                res.push_back(j);
              }
            }
            state++;
            return res;
          }
        }
        return std::nullopt;
      }

      int N;
      int K;
      uint32_t state{0};
    };

    void sort_nlist_by_I( std::vector<uint32_t> & nlist, uint32_t act, truth_table_t const& tt, truth_table_t const& mk )
    { 
      uint32_t best, cost, best_cost;
      manager.init( tt, mk );
      for( int i{0}; i<act; ++i )
      {
        best_cost = std::numeric_limits<uint32_t>::max();

        for( uint32_t j{i}; j<nlist.size(); ++j )
        {
          cost = manager.evaluate( sims[nlist[j]] );
          if( cost < best_cost  )
          {
            best_cost = cost;
            best = j;
          }
        }
        uint32_t tmp = nlist[i];
        nlist[i]=nlist[best];
        nlist[best]=tmp;

        manager.update( sims[nlist[i]] );
      }
    } 

    uint32_t _1decompose( std::vector<uint32_t> const& supp, truth_table_t const& tt, truth_table_t const& mk )
    {
      auto [lut,sim] = extract_lut( supp, tt, mk );

      funcs.push_back( lut );
      supps.push_back( supp );
      sims.push_back( sim );
      return sims.size()-1;
    }

    void clear()
    {
      sims.clear();
      supps.clear();
      funcs.clear();
      manager.reset();
    }

    void reset()
    {
      for( int i{sims.size()-1}; i>=sim_target.num_vars(); --i )
      {
        funcs.erase( funcs.begin() + funcs.size() -1 );
        supps.erase( supps.begin() + supps.size() -1 );
        sims.erase( sims.begin() + sims.size() -1 );
      }
      manager.reset();
    }

    std::optional<uint32_t> _2decompose( std::vector<uint32_t>& supp, truth_table_t const& tt, truth_table_t const& mk )
    {
      if( KILLER > _num_inserts ) return std::nullopt;
      sort_nlist_by_I( supp, supp.size(), tt, mk );

      comb_t combs( supp.size(), NUM_FANINS-1 );
      std::vector<uint32_t> free_supp;
      
      while( true )
      {
        auto comb = combs.get();
        if( !comb )
        {
          return std::nullopt;
        }

        free_supp.clear();
        manager.init( tt, mk );

        for( int i{0}; i<(*comb).size(); ++i )
        {
          free_supp.push_back( supp[(*comb)[i]] );
          manager.update( sims[supp[(*comb)[i]]] );
        }
        
        for( int m{0}; m<(1u<<manager.nMasks); ++m )
        {
          auto [tt_new, mk_new] = manager.extract_reminder();

          auto supp_bound = find_support( tt_new, mk_new );
          if( supp_bound.size() > NUM_FANINS )
            continue;

          uint32_t lit_b;
          if( supp_bound.size() <= NUM_FANINS )
          {
            auto [lut,sim] = extract_lut( supp_bound, tt_new, mk_new );
            funcs.push_back( lut );
            supps.push_back( supp_bound );
            sims.push_back( sim );
            free_supp.push_back( sims.size()-1 );

            auto [lut_f,sim_f] = extract_lut( free_supp, tt, mk );
            funcs.push_back( lut_f );
            supps.push_back( free_supp );
            sims.push_back( sim_f );
            if( kitty::equal(sim_f&mk, mk&tt) )
            {
              return sims.size()-1;
            }
            else
            {
              funcs.erase( funcs.begin() + funcs.size() - 1 );
              supps.erase( supps.begin() + supps.size() - 1 );
              sims.erase( sims.begin() + sims.size() - 1 );
              funcs.erase( funcs.begin() + funcs.size() - 1 );
              supps.erase( supps.begin() + supps.size() - 1 );
              sims.erase( sims.begin() + sims.size() - 1 );
              free_supp.erase( free_supp.begin() + free_supp.size() - 1 );
              continue;
            }
          }
        }
      }
      return std::nullopt;
    }

    std::optional<uint32_t> _Kdecompose( std::vector<uint32_t>& supp, truth_table_t const& tt, truth_table_t const& mk )
    {
      if( KILLER > _num_inserts ) return std::nullopt;
      sort_nlist_by_I( supp, supp.size(), tt, mk );

      std::vector<uint32_t> supp_f;
      for( int i{0}; i<NUM_FANINS-2; ++i )
      {
        supp_f.push_back( supp[i] );
      }

      auto tt1 = kitty::cofactor1( tt, supp[0] );
      auto mk1 = kitty::cofactor1( mk, supp[0] );
      manager.init( tt1, mk1 );
      for( int i{1}; i<supp_f.size(); ++i )
      {
        manager.update( sims[supp_f[i]] );
      }
      auto [tt1_r, mk1_r] = manager.extract_reminder();
      auto res1 = decompose_rec( tt1, mk1 );
      if( !res1 )
      {
        reset();
        return std::nullopt;
      }

      auto tt0 = kitty::cofactor0( tt, supp[0] );
      auto mk0 = kitty::cofactor0( mk, supp[0] );
      manager.init( tt0, mk0 );
      for( int i{1}; i<supp_f.size(); ++i )
      {
        manager.update( sims[supp_f[i]] );
      }
      auto [tt0_r, mk0_r] = manager.extract_reminder();
      auto res0 = decompose_rec( tt0, mk0 );
      if( !res0 )
      {
        reset();
        return std::nullopt;
      }
      supp_f.push_back( *res1 );
      supp_f.push_back( *res0 );

      auto [lut,sim] = extract_lut( supp_f, tt, mk );
      
      funcs.push_back( lut );
      supps.push_back( supp_f );
      sims.push_back( sim );
      return sims.size()-1;
    }

    std::optional<uint32_t> _Tdecompose( std::vector<uint32_t>& supp, truth_table_t const& tt, truth_table_t const& mk )
    {
      if( KILLER > _num_inserts ) return std::nullopt;
      bool upd{true};

      truth_table_t tt_r = tt;
      truth_table_t mk_r = mk;
      std::vector<uint32_t> supp_f;

      while( upd && ( supp_f.size() < NUM_FANINS - 1 ) )
      {
        upd=false;
        for( int i{0}; i<supp.size(); ++i )
        {
          if( std::find( supp_f.begin(), supp_f.end(), supp[i] ) != supp_f.end() )
            continue;

          auto tt0 = kitty::cofactor0( tt_r&mk, i );
          auto tt1 = kitty::cofactor1( tt_r&mk, i );
          auto mk0 = kitty::cofactor0( mk, i );
          auto mk1 = kitty::cofactor1( mk, i );

          if( kitty::is_const0( tt0 ) )
          {
            upd = true;
            tt_r = tt1;
            supp_f.push_back( supp[i] );
            break;
          }
          else if( kitty::is_const0( tt1 ) )
          {
            upd = true;
            tt_r = tt0;
            supp_f.push_back( supp[i] );
            break;
          }
          else if( kitty::equal( tt1&mk1, mk1 ) )
          {
            upd = true;
            tt_r = tt0;
            supp_f.push_back( supp[i] );
            break;
          }
          else if( kitty::equal( tt0&mk0, mk0 ) )
          {
            upd = true;
            tt_r = tt1;
            supp_f.push_back( supp[i] );
            break;
          }
          else if( kitty::equal( (~tt1)&mk0&mk1, tt0&mk0&mk1 ) )
          {
            upd = true;
            tt_r = (tt0 & mk0) | (~tt1&mk1);
            mk_r = mk0 | mk1;
            supp_f.push_back( supp[i] );
            break;
          }
        }
      }

      if( supp_f.size() > 0 )
      {
        auto lit_r = decompose_rec( tt_r, mk_r );
        if( !lit_r )
          return std::nullopt;

        supp_f.push_back(*lit_r);
        auto [lut,sim] = extract_lut( supp_f, tt, mk );
        if( kitty::equal( sim&mk, tt&mk ) )
        {
          funcs.push_back( lut );
          supps.push_back( supp_f );
          sims.push_back( sim );
          return sims.size()-1;
        }
        else
        {
          reset();
          return std::nullopt;
        }
      }

      return std::nullopt;
    }

    std::optional<uint32_t> decompose_rec( truth_table_t const& tt, truth_table_t const& mk )
    {      
      KILLER++;
      if( KILLER > _num_inserts ) return std::nullopt;
      auto supp = find_support( tt, mk );

      /* 0-resynthesis */
      if( supp.size() == 1 )
      {
        return supp[0];
      }
      /* 1-resynthesis */
      if( supp.size() <= NUM_FANINS )
      {
        return _1decompose( supp, tt, mk );
      }

      /* 2-resynthesis */
      if( supp.size() <= 2*NUM_FANINS-1 )
      {
        auto lit2 = _2decompose( supp, tt, mk );

        if( lit2 )
          return *lit2;
      }

      auto litT = _Tdecompose( supp, tt, mk );

      if( litT )
      {
        return *litT;
      }

      auto litK = _Kdecompose( supp, tt, mk ); 

      if( litK )
        return *litK;
        
      return std::nullopt;

    }

    std::optional<uint32_t> decompose( truth_table_t const& tt, truth_table_t const& mk, uint32_t num_inserts )
    {
      KILLER=0;
      _num_inserts = 20;

      /* initialization */
      sims.clear();
      sim_target = tt;

      funcs.clear();
      supps.clear();

      for( int i{0}; i<tt.num_vars(); ++i )
      {
        funcs.emplace_back(1u);
        supps.push_back({i});
        sims.emplace_back(tt.num_vars());
        kitty::create_nth_var( sims[i], i );
      }

      return decompose_rec( tt, mk );

    }

    std::optional<uint32_t> decompose( truth_table_t const& tt, uint32_t num_inserts )
    {
      auto mk = ~tt.construct();
      return decompose( tt, mk, num_inserts );
    }

    uint32_t num_luts()
    {
        return (uint32_t)funcs.size() - sim_target.num_vars();
    }

    template<class Ntk>
    signal<Ntk> add_to_network( Ntk& ntk, std::vector<node<Ntk>> nodes )
    {
      signal<Ntk> res = 0;
      std::vector<signal<Ntk>> children;
      std::vector<signal<Ntk>> supp;
      for( auto s : nodes )
        supp.push_back( ntk.make_signal(s) );

      for( int i{sim_target.num_vars()}; i<sims.size(); ++i )
      {
        children.clear();
        for( auto x : supps[i] )
        {
          children.push_back( supp[x] );
        }
        res = ntk.create_node( children, funcs[i] );
        supp.push_back( res );
      }

      return res;
    }

    void print()
    {
        for( int i{0}; i<sims.size(); ++i )
        {
            printf( "%3d ", i );
            kitty::print_binary( sims[i] );
            printf("\n");
        }

        printf( "    " );
        kitty::print_binary( sim_target );
        printf("\n");
    }

    large_lig_index_list::element_type to_index_list( large_lig_index_list& index_list, std::vector<large_lig_index_list::element_type> lits )
    {
      std::vector<large_lig_index_list::element_type> supp;
      large_lig_index_list::element_type lit;

      for( int i{sim_target.num_vars()}; i<sims.size(); ++i )
      {
        supp.clear();
        for( auto x : supps[i] )
        {
          supp.push_back( lits[x] );
        }
        lit = index_list.add_function( supp, funcs[i] );
        lits.push_back( lit );
      }

      return lit;
    }

    public:
        truth_table_t sim_target;
        spfd_covering_manager_t<truth_table_t, 1u << MAX_WIDTH> manager;

        std::vector<std::vector<uint32_t>> supps;
        std::vector<truth_table_t> funcs;
        std::vector<truth_table_t> sims;
        uint32_t KILLER{0};

        uint32_t _num_inserts;

  };

} // namespace mockturtle


