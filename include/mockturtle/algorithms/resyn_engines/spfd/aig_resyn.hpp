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
  \file aig_resyn.hpp
  \brief Resynthesis by recursive decomposition for AIGs or aigs.
  (based on ABC's implementation in `giaResub.c` by Alan Mishchenko)

  \author Andrea Costamagna
*/

#pragma once

#include "../../../utils/index_list.hpp"
#include "../../../utils/node_map.hpp"
#include "../../../utils/stopwatch.hpp"
#include "../../../utils/tech_library.hpp"
#include "../../node_resynthesis/xag_npn.hpp"

#include <abcresub/abcresub.hpp>
#include <fmt/format.h>
#include <kitty/kitty.hpp>
#include <kitty/constructors.hpp>

#include <algorithm>
#include <optional>
#include <type_traits>
#include <vector>

#include <random>
#include <limits>

namespace mockturtle
{

namespace spfd
{

namespace aig
{

template<class DIV>
bool comparator(const DIV& lhs, const DIV& rhs) 
{
  return lhs.cost < rhs.cost;
}


bool VERBOSE{false};
bool VERBOSE1{false};

template<class TT>
void print_tt_with_dcs( TT tt, TT mk )
{
  for( auto m = tt.num_bits()-1; m >= 0; --m )
  {
    if( kitty::get_bit( mk, m ) == 1 )
    {
      if( kitty::get_bit( tt, m ) == 1 )
      {
        printf("1");
      }
      else
      {
        printf("0");
      }
    }
    else
    {
      printf("*");
    }
  }
  printf("\n");
}

std::mt19937 RNG(5);

template<class TT> TT compute_buff( TT const& tt1, TT const& tt2 ){ return tt1; }
template<class TT> TT compute_pa00( TT const& tt1, TT const& tt2 ){ return ~tt1 & ~tt2; }
template<class TT> TT compute_pa01( TT const& tt1, TT const& tt2 ){ return ~tt1 &  tt2; }
template<class TT> TT compute_pa10( TT const& tt1, TT const& tt2 ){ return  tt1 & ~tt2; }
template<class TT> TT compute_pa11( TT const& tt1, TT const& tt2 ){ return  tt1 &  tt2; }
template<class TT> TT compute_exor( TT const& tt1, TT const& tt2 ){ return  tt1 ^  tt2; }

template<class LIST> uint32_t add_buff_to_list( LIST& list, uint32_t lit1, uint32_t lit2 ){ return lit1; }
template<class LIST> uint32_t add_pa00_to_list( LIST& list, uint32_t lit1, uint32_t lit2 ){ return list.add_and( lit1 ^ 0x1, lit2 ^ 0x1 ); }
template<class LIST> uint32_t add_pa01_to_list( LIST& list, uint32_t lit1, uint32_t lit2 ){ return list.add_and( lit1 ^ 0x1, lit2 ); }
template<class LIST> uint32_t add_pa10_to_list( LIST& list, uint32_t lit1, uint32_t lit2 ){ return list.add_and( lit1, lit2 ^ 0x1 ); }
template<class LIST> uint32_t add_pa11_to_list( LIST& list, uint32_t lit1, uint32_t lit2 ){ return list.add_and( lit1, lit2 ); }
template<class LIST> uint32_t add_exor_to_list( LIST& list, uint32_t lit1, uint32_t lit2 ){ return list.add_xor( lit1, lit2 ); }


struct aig_resyn_static_params
{
  using base_type = aig_resyn_static_params;

  /*! \brief Maximum number of binate divisors to be considered. */
  static constexpr uint32_t max_binates{ 50u };

  /*! \brief Reserved capacity for divisor truth tables (number of divisors). */
  static constexpr uint32_t reserve{ 200u };

  /*! \brief Whether to consider single XOR gates (i.e., using aigs instead of AIGs). */
  static constexpr bool use_xor{ true };

  /*! \brief Whether to copy truth tables. */
  static constexpr bool copy_tts{ false };

  /*! \brief Whether to preserve depth. */
  static constexpr bool preserve_depth{ false };

  /*! \brief Whether the divisors have uniform costs (size and depth, whenever relevant). */
  static constexpr bool uniform_div_cost{ true };

  /*! \brief Size cost of each AND gate. */
  static constexpr uint32_t size_cost_of_and{ 1u };

  /*! \brief Size cost of each XOR gate (only relevant when `use_xor = true`). */
  static constexpr uint32_t size_cost_of_xor{ 1u };

  /*! \brief Depth cost of each AND gate (only relevant when `preserve_depth = true`). */
  static constexpr uint32_t depth_cost_of_and{ 1u };

  /*! \brief Depth cost of each XOR gate (only relevant when `preserve_depth = true` and `use_xor = true`). */
  static constexpr uint32_t depth_cost_of_xor{ 1u };


  /*! \brief Maximum number of support variables. */
  static constexpr uint32_t max_support_size{ 4u };
  static constexpr uint32_t max_num_support_samplings{ 20u };
  static constexpr uint32_t max_resyn_attempts{ 1 };

  static constexpr double beta_support{100};
  static constexpr double beta_synthesis{100};

  static constexpr bool try_boolean_matching{ false };

  static constexpr bool use_greedy_support{ false };
  static constexpr bool use_boltz{ false };
  static constexpr bool use_enum{ false };


  static constexpr bool use_spfd_synthesis{ true };
  static constexpr bool use_decomposition{ false };

  using truth_table_storage_type = void;
  using node_type = void;
};

template<class TT>
struct aig_resyn_static_params_default : public aig_resyn_static_params
{
  using truth_table_storage_type = std::vector<TT>;
  using node_type = uint32_t;
  static constexpr bool use_xor = false;
};

template<class Ntk, uint32_t SUPP_SIZE, uint32_t N_SAMPL, uint32_t N_RESYN, bool IS_BMATCH>
struct aig_resyn_static_params_for_sim_resub : public aig_resyn_static_params
{
  using truth_table_storage_type = incomplete_node_map<kitty::partial_truth_table, Ntk>;
  using node_type = typename Ntk::node;
  static constexpr bool use_xor = false;
  static constexpr uint32_t max_support_size = SUPP_SIZE;
  static constexpr uint32_t max_num_support_samplings = N_SAMPL;
  static constexpr uint32_t max_resyn_attempts = N_RESYN;
  static constexpr uint32_t try_boolean_matching = IS_BMATCH;
  static constexpr uint32_t use_greedy_support = true;

};

struct aig_resyn_stats
{
  /*! \brief Time for finding 0-resub and collecting unate literals. */
  stopwatch<>::duration time_unate{ 0 };

  /*! \brief Time for finding 1-resub. */
  stopwatch<>::duration time_resub1{ 0 };

  /*! \brief Time for finding 2-resub. */
  stopwatch<>::duration time_resub2{ 0 };

  /*! \brief Time for finding 3-resub. */
  stopwatch<>::duration time_resub3{ 0 };

  /*! \brief Time for sorting unate literals and unate pairs. */
  stopwatch<>::duration time_sort{ 0 };

  /*! \brief Time for collecting unate pairs. */
  stopwatch<>::duration time_collect_pairs{ 0 };

  /*! \brief Time for dividing the target and recursive call. */
  stopwatch<>::duration time_divide{ 0 };

  void report() const
  {
    fmt::print( "[i]         <xag_resyn_decompose>\n" );
    fmt::print( "[i]             0-resub      : {:>5.2f} secs\n", to_seconds( time_unate ) );
    fmt::print( "[i]             1-resub      : {:>5.2f} secs\n", to_seconds( time_resub1 ) );
    fmt::print( "[i]             2-resub      : {:>5.2f} secs\n", to_seconds( time_resub2 ) );
    fmt::print( "[i]             3-resub      : {:>5.2f} secs\n", to_seconds( time_resub3 ) );
    fmt::print( "[i]             sort         : {:>5.2f} secs\n", to_seconds( time_sort ) );
    fmt::print( "[i]             collect pairs: {:>5.2f} secs\n", to_seconds( time_collect_pairs ) );
    fmt::print( "[i]             dividing     : {:>5.2f} secs\n", to_seconds( time_divide ) );
    fmt::print( "[i]             sort         : {:>5.2f} secs\n", to_seconds( time_sort ) );
  }
};

/*! \brief Logic resynthesis engine for AIGs or aigs.
 *
 * The algorithm is based on ABC's implementation in `giaResub.c` by Alan Mishchenko.
 *
 * Divisors are classified as positive unate (not overlapping with target offset),
 * negative unate (not overlapping with target onset), or binate (overlapping with
 * both onset and offset). Furthermore, pairs of binate divisors are combined with
 * an AND operation and considering all possible input polarities and again classified
 * as positive unate, negative unate or binate. Simple solutions of zero cost
 * (one unate divisor), one node (two unate divisors), two nodes (one unate divisor +
 * one unate pair), and three nodes (two unate pairs) are exhaustively examined.
 * When no simple solutions can be found, the algorithm heuristically chooses an unate
 * divisor or an unate pair to divide the target function with and recursively calls
 * itself to decompose the remainder function.
   \verbatim embed:rst

   Example

   .. code-block:: c++

      using TT = kitty::static_truth_table<6>;
      const std::vector<aig_network::node> divisors = ...;
      const node_map<TT, aig_network> tts = ...;
      const TT target = ..., care = ...;
      aig_resyn_stats st;
      aig_resyn<TT, node_map<TT, aig_network>, false, false, aig_network::node> resyn( st );
      auto result = resyn( target, care, divisors.begin(), divisors.end(), tts );
   \endverbatim
 */
template<class TT, class database_t, class static_params = aig_resyn_static_params_default<TT>>
class aig_resyn
{
public:
  using stats = aig_resyn_stats;
  using index_list_t = large_xag_index_list;
  using truth_table_t = TT;
  using truth_table4_t = kitty::static_truth_table<4u>;
  using truth_tableK_t = kitty::static_truth_table<static_params::max_support_size>;
  using divisor_id_t = uint32_t;

private:
  struct unate_lit
  {
    unate_lit( uint32_t l )
        : lit( l )
    {}

    bool operator==( unate_lit const& other ) const
    {
      return lit == other.lit;
    }

    uint32_t lit;
    uint32_t score{ 0 };
  };

  struct fanin_pair
  {
    fanin_pair( uint32_t l1, uint32_t l2 )
        : lit1( l1 < l2 ? l1 : l2 ), lit2( l1 < l2 ? l2 : l1 )
    {}

    fanin_pair( uint32_t l1, uint32_t l2, bool is_xor )
        : lit1( l1 > l2 ? l1 : l2 ), lit2( l1 > l2 ? l2 : l1 )
    {
      (void)is_xor;
    }

    bool operator==( fanin_pair const& other ) const
    {
      return lit1 == other.lit1 && lit2 == other.lit2;
    }

    uint32_t lit1, lit2;
    uint32_t score{ 0 };
  };

private:

  struct divisort_t;
  struct candidate_t;
  template<class LTT, uint32_t CAP> struct spfd_manager_t;
  struct functional_library_t;


  struct gate_t
  {
    gate_t(){};
    gate_t( uint32_t id, truth_tableK_t (*pF)( truth_tableK_t const&, truth_tableK_t const& ), uint32_t (*pG)( index_list_t&, uint32_t, uint32_t ) ) : id(id), pF(pF), pG(pG){}

    truth_tableK_t compute( truth_tableK_t const& tt1, truth_tableK_t const& tt2 ){ return pF( tt1, tt2 ); };
    truth_tableK_t compute( truth_tableK_t const& tt1 ){ return pF( tt1, tt1 ); };

    uint32_t add_to_list( index_list_t& list, uint32_t lit1, uint32_t lit2 ){ return pG( list, lit1, lit2 ); }
    uint32_t add_to_list( index_list_t& list, uint32_t lit1 ){ return pG( list, lit1, lit1 ); }

    bool is_buffer()
    {
      return id == 0x0;
    }

    bool is_pa00()
    {
      return id == 0x1;
    }

    bool is_pa01()
    {
      return id == 0x2;
    }

    bool is_pa10()
    {
      return id == 0x4;
    }

    bool is_pa11()
    {
      return id == 0x8;
    }

    bool is_exor()
    {
      return id == 0x6;
    }

    uint32_t id;
    truth_tableK_t (*pF)( truth_tableK_t const&, truth_tableK_t const& );
    uint32_t (*pG)( index_list_t&, uint32_t, uint32_t );
  };

  struct functional_library_t
  {
    functional_library_t()
    {
      gates1[0] = gate_t{ 0x0, &compute_buff<truth_tableK_t>, &add_buff_to_list<index_list_t> };
      gates2[0] = gate_t{ 0x1, &compute_pa00<truth_tableK_t>, &add_pa00_to_list<index_list_t> };
      gates2[1] = gate_t{ 0x2, &compute_pa01<truth_tableK_t>, &add_pa01_to_list<index_list_t> };
      gates2[2] = gate_t{ 0x4, &compute_pa10<truth_tableK_t>, &add_pa10_to_list<index_list_t> };
      gates2[3] = gate_t{ 0x8, &compute_pa11<truth_tableK_t>, &add_pa11_to_list<index_list_t> };
    }

    std::array<gate_t, 1u> gates1;
    std::array<gate_t, 4u> gates2;
  };

  struct scored_divisor_t
  {
    scored_divisor_t( uint32_t d, double c )
        : div( d ), cost( c )
    {}

    bool operator==( scored_divisor_t const& other ) const
    {
      return div == other.div;
    }

    bool operator<( const scored_divisor_t & other ) const
    {
      return cost < other.cost;
    }

    uint32_t div;
    double cost;
  };

  struct scored_divisors_t
  {
    scored_divisors_t(){};

    void emplace_back( uint32_t div, double cost )
    {
      divs.emplace_back( div, cost );
    } 

    void sort()
    {
      std::sort( divs.begin(), divs.end() );
    }

    void print()
    {
      for( auto div : divs )
      {
        printf("(%d,%f) ", div.div, div.cost );
      }
      printf("\n");
    }

    std::vector<scored_divisor_t> divs;
  };

  struct divisor_t
  {
    divisor_t( truth_tableK_t func, uint32_t lit ) : func(func), lit(lit) {}
    divisor_t( truth_tableK_t func ) : func( func ) {}
    divisor_t(){}

    truth_tableK_t func;
    uint32_t lit;
  };

  struct divisors_t
  {
    divisors_t(){}

    void emplace_back( truth_tableK_t func, uint32_t lit )
    {
      divs.emplace_back( func, lit );
    }

    uint32_t size()
    {
      return divs.size();
    }

    inline divisor_t const& operator[]( uint32_t idx ) const
    {
      return divs[idx];
    }

    inline truth_tableK_t const& get_div( uint32_t idx ) const
    {
      return divs[idx].func;
    } 

    void set_support( std::vector<uint32_t> const& supp, std::array<truth_tableK_t,static_params::max_support_size> const& funcs )
    {
      divs.clear();

      for( auto i{0}; i<supp.size(); ++i )
      {
        divs.emplace_back( funcs[i], supp[i] << 1u );
      }
    }

    void set_target( truth_tableK_t const& func, truth_tableK_t const& care )
    {
      spfd.init( func, care );
    }

    void clear()
    {
      divs.clear();
      spfd.reset();
    }

    bool update( index_list_t & list, functional_library_t const& functional_library, uint32_t max_num_gates )
    {
      uint32_t num_buffers{0};
      std::vector<divisor_t> new_divs;

      std::vector<candidate_t> candidates;
      uint32_t cand_id{0};
      for( uint32_t v1{0}; v1 < divs.size(); ++v1 )
      {
        for( auto gate : functional_library.gates1 )
        {
          candidates.emplace_back( cand_id++, gate, divs[v1] );
        }

        for( uint32_t v2 = v1 + 1; v2 < divs.size(); ++v2 )
        {
          for( auto gate : functional_library.gates2 )
          {
            //if( gate.is_exor() && ( ( list.num_gates() + 3 > max_num_gates ) || ( divs.size() != 2 ) ) ) continue;
            candidates.emplace_back( cand_id++, gate, divs[v1], divs[v2] );
          }
        }
      }

      double cost;
      double min_cost = std::numeric_limits<double>::max();
      double max_cost = std::numeric_limits<double>::min();
      std::set<uint32_t> set_used;
      std::uniform_real_distribution<> U01(0,1);

      spfd.reset();

      while( !spfd.is_covered() && new_divs.size() < static_params::max_support_size )
      {
        for( auto & cand : candidates )
        {
          cost = spfd.evaluate( cand.compute() );
          cand.cost = cost;
          if( cost < min_cost && set_used.find(cand.id) == set_used.end() ) min_cost = cost;
          if( cost > max_cost && set_used.find(cand.id) == set_used.end() ) max_cost = cost;
        }

        double Z{0};
        bool copy_previous;
        for( auto & cand : candidates )
        {
          copy_previous = set_used.find(cand.id) != set_used.end();
          copy_previous |= (cand.gate.is_buffer() && ( num_buffers >= divs.size()-1 ));
          Z = cand.update_cost( Z, min_cost, max_cost, copy_previous );
        }

        double rnd = U01(RNG);
        bool is_updated{false};

        for( auto & cand : candidates )
        {
          if( rnd*Z <= cand.cost )
          {
            set_used.insert( cand.id );
            if( cand.gate.is_buffer() ) num_buffers++;

            truth_tableK_t tt = cand.compute();
            new_divs.emplace_back( tt, cand.add_to_list( list ) );
            spfd.update( tt );
            is_updated = true;
            break;
          }
        }
        if( !is_updated )
        {
          return false;
        }
      }
      if( spfd.is_covered() )
      {
        divs = new_divs;
        return true;
      }

      return false;
    }

    bool update_f( index_list_t & list, functional_library_t const& functional_library, uint32_t max_num_gates )
    {
      uint32_t num_buffers{0};
      std::vector<divisor_t> new_divs;

      std::vector<candidate_t> candidates;
      std::vector<uint32_t> best_candidates;;
      uint32_t cand_id{0};
      for( uint32_t v1{0}; v1 < divs.size(); ++v1 )
      {
        for( auto gate : functional_library.gates1 )
        {
          candidates.emplace_back( cand_id++, gate, divs[v1] );
        }

        for( uint32_t v2 = v1 + 1; v2 < divs.size(); ++v2 )
        {
          for( auto gate : functional_library.gates2 )
          {
            candidates.emplace_back( cand_id++, gate, divs[v1], divs[v2] );
          }
        }
      }

      double cost;
      std::set<uint32_t> set_used;

      spfd.reset();

      while( !spfd.is_covered() && new_divs.size() < static_params::max_support_size )
      {
        double best_cost = std::numeric_limits<double>::max();
        for( auto & cand : candidates )
        {
          cost = spfd.evaluate( cand.compute() );
          cand.cost = cost;
          if( !(cand.gate.is_buffer() && ( num_buffers >= divs.size()-1 )) )
          {
            if( cost < best_cost )
            {
              best_candidates.clear();
              best_cost = cost;
              best_candidates = { cand.id };
            }
            else if( cost == best_cost )
            {
              best_candidates.push_back( cand.id );
            }
          }
        }
        if( best_candidates.size() == 0 ) return false;
        std::uniform_int_distribution<> U01( 0, best_candidates.size() - 1 );

        uint32_t rnd = U01(RNG);

        set_used.insert( best_candidates[rnd] );
        if( candidates[best_candidates[rnd]].gate.is_buffer() ) num_buffers++;

        truth_tableK_t tt = candidates[best_candidates[rnd]].compute();
        new_divs.emplace_back( tt, candidates[best_candidates[rnd]].add_to_list( list ) );
        spfd.update( tt );
      }
      if( spfd.is_covered() )
      {
        divs = new_divs;
        return true;
      }

      return false;
    }

    std::vector<divisor_t> divs;
    spfd_manager_t<truth_tableK_t, 1<<static_params::max_support_size> spfd;

  };

  template<class LTT, uint32_t CAP>
  struct spfd_manager_t
  {
    spfd_manager_t(){}
    
    void init( LTT const& target, LTT const& careset )
    {
      care = careset;
      func[1] =  target & careset;
      func[0] = ~target & careset;
      reset();
    }

    void reset()
    {
      masks[0] = care;
      nMasks = 1;
      nEdges = kitty::count_ones( func[1] ) * kitty::count_ones( func[0] );  
      killed[0] = nEdges > 0 ? false : true;
      nKills = nEdges > 0 ? 0 : 1;
    }

    bool update( LTT const& tt )
    {
      nEdges = 0;
      for( uint32_t iMask{0}; iMask < nMasks; ++iMask )
      {
        if( killed[iMask] )
        {
          killed[nMasks+iMask] = true;
          nKills++;
        }
        else
        {
          masks[nMasks+iMask] = masks[iMask] & tt;
          masks[iMask] &= ~tt;

          if( ( kitty::count_ones( masks[nMasks+iMask] & func[1] ) == 0 ) || ( kitty::count_ones( masks[nMasks+iMask] & func[0] ) == 0 ) )
          {
            killed[nMasks+iMask] = true;
            nKills++;
          }
          else
          {
            killed[nMasks+iMask] = false;
            nEdges += kitty::count_ones( func[1] & masks[nMasks+iMask] ) * kitty::count_ones( func[0] & masks[nMasks+iMask] );  
          }

          if( kitty::count_ones( masks[iMask] & func[1] ) == 0 || kitty::count_ones( masks[iMask] & func[0] ) == 0 )
          {
            killed[iMask] = true;
            nKills++;
          }
          else
          {
            killed[iMask] = false;
            nEdges += kitty::count_ones( func[1] & masks[iMask] ) * kitty::count_ones( func[0] & masks[iMask] );  
          }
        }
      }
      nMasks = nMasks * 2;
      return true;
    }

    double evaluate( LTT const& tt )
    {
      double res=0;
      for( auto iMask{0}; iMask<nMasks; ++iMask )
      {
        if( !killed[iMask] )
        {
          res+= kitty::count_ones( func[1] & masks[iMask] &  tt ) * kitty::count_ones( func[0] & masks[iMask] & tt )/nEdges;  
          res+= kitty::count_ones( func[1] & masks[iMask] & ~tt ) * kitty::count_ones( func[0] & masks[iMask] & ~tt )/nEdges;  
        }
      }
      return res;
    } 

    double evaluate( LTT const& tt1, LTT const& tt2 )
    {
      double res=0;
      for( auto iMask{0}; iMask<nMasks; ++iMask )
      {
        if( !killed[iMask] )
        {
          res+= kitty::count_ones( func[1] & masks[iMask] &  tt1 & tt2 ) * kitty::count_ones( func[0] & masks[iMask] &  tt1 & tt2 )/nEdges;  
          res+= kitty::count_ones( func[1] & masks[iMask] & ~tt1 & tt2 ) * kitty::count_ones( func[0] & masks[iMask] & ~tt1 & tt2 )/nEdges;  
          res+= kitty::count_ones( func[1] & masks[iMask] & ~tt1 &~tt2 ) * kitty::count_ones( func[0] & masks[iMask] & ~tt1 &~tt2 )/nEdges;  
          res+= kitty::count_ones( func[1] & masks[iMask] &  tt1 &~tt2 ) * kitty::count_ones( func[0] & masks[iMask] &  tt1 &~tt2 )/nEdges;  
        }
      }
      return res;
    } 

    bool is_covered()
    {
      return nMasks <= nKills;
    }

    bool is_saturated()
    {
      return nMasks >= CAP;
    }

    std::array<LTT, CAP> masks;
    std::array<bool, CAP> killed;
    uint32_t nMasks;
    uint32_t nKills;
    double nEdges;
    LTT care;
    std::array<LTT, 2u> func;
  };

  template<class LTT, uint32_t CAP>
  struct u_spfd_manager_t
  {
    u_spfd_manager_t(){}
    
    void init( LTT const& target, LTT const& careset )
    {
      care = careset;
      func[1] =  target & careset;
      func[0] = ~target & careset;
      reset();
    }

    void reset()
    {
      masks[0] = care;
      nMasks = 1;
      nEdges = kitty::count_ones( func[1] ) * kitty::count_ones( func[0] );  
      killed[0] = nEdges > 0 ? false : true;
      nKills = nEdges > 0 ? 0 : 1;
    }

    bool update( LTT const& tt )
    {
      nEdges = 0;
      for( uint32_t iMask{0}; iMask < nMasks; ++iMask )
      {
        if( killed[iMask] )
        {
          killed[nMasks+iMask] = true;
          nKills++;
        }
        else
        {
          masks[nMasks+iMask] = masks[iMask] & tt;
          masks[iMask] &= ~tt;

          if( ( kitty::count_ones( masks[nMasks+iMask] & func[1] ) == 0 ) || ( kitty::count_ones( masks[nMasks+iMask] & func[0] ) == 0 ) )
          {
            killed[nMasks+iMask] = true;
            nKills++;
          }
          else
          {
            killed[nMasks+iMask] = false;
            nEdges += kitty::count_ones( func[1] & masks[nMasks+iMask] ) * kitty::count_ones( func[0] & masks[nMasks+iMask] );  
          }

          if( kitty::count_ones( masks[iMask] & func[1] ) == 0 || kitty::count_ones( masks[iMask] & func[0] ) == 0 )
          {
            killed[iMask] = true;
            nKills++;
          }
          else
          {
            killed[iMask] = false;
            nEdges += kitty::count_ones( func[1] & masks[iMask] ) * kitty::count_ones( func[0] & masks[iMask] );  
          }
        }
      }
      nMasks = nMasks * 2;
      return true;
    }

    uint32_t evaluate( LTT const& tt )
    {
      uint32_t res=0;
      for( auto iMask{0}; iMask<nMasks; ++iMask )
      {
        if( !killed[iMask] )
        {
          res+= kitty::count_ones( func[1] & masks[iMask] &  tt ) * kitty::count_ones( func[0] & masks[iMask] & tt );  
          res+= kitty::count_ones( func[1] & masks[iMask] & ~tt ) * kitty::count_ones( func[0] & masks[iMask] & ~tt );  
        }
      }
      return res;
    } 

    uint32_t evaluate( LTT const& tt1, LTT const& tt2 )
    {
      uint32_t res=0;
      for( auto iMask{0}; iMask<nMasks; ++iMask )
      {
        if( !killed[iMask] )
        {
          res+= kitty::count_ones( func[1] & masks[iMask] &  tt1 & tt2 ) * kitty::count_ones( func[0] & masks[iMask] & tt1 & tt2 );  
          res+= kitty::count_ones( func[1] & masks[iMask] & ~tt1 & tt2 ) * kitty::count_ones( func[0] & masks[iMask] & ~tt1 & tt2 );  
          res+= kitty::count_ones( func[1] & masks[iMask] & ~tt1 &~tt2 ) * kitty::count_ones( func[0] & masks[iMask] & ~tt1 & ~tt2 );  
          res+= kitty::count_ones( func[1] & masks[iMask] & tt1 & ~tt2 ) * kitty::count_ones( func[0] & masks[iMask] & tt1 &~tt2 );  
        }
      }
      return res;
    } 

    bool is_covered()
    {
      return nMasks <= nKills;
    }

    bool is_saturated()
    {
      return nMasks >= CAP;
    }

    std::array<LTT, CAP> masks;
    std::array<bool, CAP> killed;
    uint32_t nMasks;
    uint32_t nKills;
    uint32_t nEdges;
    LTT care;
    std::array<LTT, 2u> func;
  };

  struct candidate_t
  {
    candidate_t(){}
    candidate_t( uint32_t id, gate_t gate, divisor_t const& div1 ) : id(id), gate(gate), div1(div1), div2(div1){}
    candidate_t( uint32_t id, gate_t gate, divisor_t const& div1, divisor_t const& div2 ) : id(id), gate(gate), div1(div1), div2(div2){}

    uint32_t add_to_list( index_list_t & list ){ return gate.add_to_list( list, div1.lit, div2.lit ); }

    truth_tableK_t compute(){ return gate.compute( div1.func, div2.func ); }

    double update_cost( double const& cost_previous, double const& min_cost, double const& max_cost, bool copy_previous )
    {
      if( copy_previous )
      {
        cost = cost_previous;
      }
      else
      {
        cost = cost_previous + exp( -static_params::beta_synthesis*( cost - min_cost )/( max_cost - min_cost ) );
      }
      return cost;
    }

    uint32_t id;
    gate_t gate;
    double cost;
    divisor_t const& div1;
    divisor_t const& div2;
  };

public:
  explicit aig_resyn( database_t database, stats& st ) noexcept
      : st( st ), _database(database)
  {
    static_assert( std::is_same_v<typename static_params::base_type, aig_resyn_static_params>, "Invalid static_params type" );
    static_assert( !( static_params::uniform_div_cost && static_params::preserve_depth ), "If depth is to be preserved, divisor depth cost must be provided (usually not uniform)" );
    divisors.reserve( static_params::reserve );
    //lits.reserve( static_params::reserve );
    _scored_divs.reserve( static_params::reserve );
  }

  /*! \brief Perform aig resynthesis.
   *
   * `tts[*begin]` must be of type `TT`.
   * Moreover, if `static_params::copy_tts = false`, `*begin` must be of type `static_params::node_type`.
   *
   * \param target Truth table of the target function.
   * \param care Truth table of the care set.
   * \param begin Begin iterator to divisor nodes.
   * \param end End iterator to divisor nodes.
   * \param tts A data structure (e.g. std::vector<TT>) that stores the truth tables of the divisor functions.
   * \param max_size Maximum number of nodes allowed in the dependency circuit.
   */
  template<class iterator_type,
           bool enabled = static_params::uniform_div_cost && !static_params::preserve_depth, typename = std::enable_if_t<enabled>>
  std::optional<index_list_t> operator()( TT const& target, TT const& care, iterator_type begin, iterator_type end, typename static_params::truth_table_storage_type const& tts, uint32_t max_size = std::numeric_limits<uint32_t>::max() )
  {
    static_assert( static_params::copy_tts || std::is_same_v<typename std::iterator_traits<iterator_type>::value_type, typename static_params::node_type>, "iterator_type does not dereference to static_params::node_type" );
    _scored_divs.clear();

    _past_supports_count.clear();

    ptts = &tts;
    on_off_sets[0] = ~target & care;
    on_off_sets[1] = target & care;

    _gSPFD.init( target, care );
    _uSPFD.init( target, care );

    divisors.resize( 1 ); /* clear previous data and reserve 1 dummy node for constant */
    while ( begin != end )
    {
      if constexpr ( static_params::copy_tts )
      {
        divisors.emplace_back( ( *ptts )[*begin] );
      }
      else
      {
        divisors.emplace_back( *begin );
      }
      ++begin;
    }

    double cost;
    for( uint32_t iDiv{1}; iDiv<divisors.size(); ++iDiv )
    {
      _scored_divs.emplace_back( iDiv, _gSPFD.evaluate( get_div( iDiv ) ) );  
    }

    //_scored_divs.print();
    std::sort( _scored_divs.begin(), _scored_divs.end() );
    //_scored_divs.print();

    for( uint32_t i{0}; i<4; ++i )
      kitty::create_nth_var( _xs4[i], i );

    for( uint32_t i{0}; i<static_params::max_support_size; ++i )
      kitty::create_nth_var( _xsK[i], i );

    return compute_function( max_size );
  }

  template<class iterator_type, class Fn,
           bool enabled = !static_params::uniform_div_cost && !static_params::preserve_depth, typename = std::enable_if_t<enabled>>
  std::optional<index_list_t> operator()( TT const& target, TT const& care, iterator_type begin, iterator_type end, typename static_params::truth_table_storage_type const& tts, Fn&& size_cost, uint32_t max_size = std::numeric_limits<uint32_t>::max() )
  {}

  template<class iterator_type, class Fn,
           bool enabled = !static_params::uniform_div_cost && static_params::preserve_depth, typename = std::enable_if_t<enabled>>
  std::optional<index_list_t> operator()( TT const& target, TT const& care, iterator_type begin, iterator_type end, typename static_params::truth_table_storage_type const& tts, Fn&& size_cost, Fn&& depth_cost, uint32_t max_size = std::numeric_limits<uint32_t>::max(), uint32_t max_depth = std::numeric_limits<uint32_t>::max() )
  {}

private:
  std::optional<index_list_t> compute_function( uint32_t num_inserts )
  {
    index_list.clear();
    index_list.add_inputs( divisors.size() - 1 );
    auto const lit = compute_function_rec( num_inserts );
    if ( lit )
    {
      assert( index_list.num_gates() <= num_inserts );
      index_list.add_output( *lit );
    //  std::cout << to_index_list_string( index_list ) << std::endl; 

      return index_list;
    }
    return std::nullopt;
  }

  std::optional<uint32_t> compute_function_rec( uint32_t num_inserts )
  {

    pos_unate_lits.clear();
    neg_unate_lits.clear();
    binate_divs.clear();
    pos_unate_pairs.clear();
    neg_unate_pairs.clear();

    /* try 0-resub and collect unate literals */
    auto const res0 = call_with_stopwatch( st.time_unate, [&]() {
      return find_one_unate();
    } );
    if ( res0 )
    {
      return *res0;
    }
    if ( num_inserts == 0u )
    {
      return std::nullopt;
    }

    if( static_params::use_decomposition )
    {
      printf("DECOMPOSE!!!\n");
      /* sort unate literals and try 1-resub */
      call_with_stopwatch( st.time_sort, [&]() {
        sort_unate_lits( pos_unate_lits, 1 );
        sort_unate_lits( neg_unate_lits, 0 );
      } );
      auto const res1or = call_with_stopwatch( st.time_resub1, [&]() {
        return find_div_div( pos_unate_lits, 1 );
      } );
      if ( res1or )
      {
        return *res1or;
      }
      auto const res1and = call_with_stopwatch( st.time_resub1, [&]() {
        return find_div_div( neg_unate_lits, 0 );
      } );
      if ( res1and )
      {
        return *res1and;
      }

      if ( binate_divs.size() > static_params::max_binates )
      {
        binate_divs.resize( static_params::max_binates );
      }

      if constexpr ( static_params::use_xor )
      {
        /* collect XOR-type unate pairs and try 1-resub with XOR */
        auto const res1xor = find_xor();
        if ( res1xor )
        {
          return *res1xor;
        }
      }

      if ( num_inserts > 1u )
      {

        /* collect AND-type unate pairs and sort (both types), then try 2- and 3-resub */
        call_with_stopwatch( st.time_collect_pairs, [&]() {
          collect_unate_pairs();
        } );
        call_with_stopwatch( st.time_sort, [&]() {
          sort_unate_pairs( pos_unate_pairs, 1 );
          sort_unate_pairs( neg_unate_pairs, 0 );
        } );
        auto const res2or = call_with_stopwatch( st.time_resub2, [&]() {
          return find_div_pair( pos_unate_lits, pos_unate_pairs, 1 );
        } );
        if ( res2or )
        {
          return *res2or;
        }
        auto const res2and = call_with_stopwatch( st.time_resub2, [&]() {
          return find_div_pair( neg_unate_lits, neg_unate_pairs, 0 );
        } );
        if ( res2and )
        {
          return *res2and;
        }

        if ( num_inserts >= 3u )
        {
          auto const res3or = call_with_stopwatch( st.time_resub3, [&]() {
            return find_pair_pair( pos_unate_pairs, 1 );
          } );
          if ( res3or )
          {
            return *res3or;
          }
          auto const res3and = call_with_stopwatch( st.time_resub3, [&]() {
            return find_pair_pair( neg_unate_pairs, 0 );
          } );
          if ( res3and )
          {
            return *res3and;
          }
        }

        /* choose something to divide and recursive call on the remainder */
        /* Note: dividing = AND the on-set (if using positive unate) or the off-set (if using negative unate)
                            with the *negation* of the divisor/pair (subtracting) */
        uint32_t on_off_div, on_off_pair;
        uint32_t score_div = 0, score_pair = 0;

        call_with_stopwatch( st.time_divide, [&]() {
          if ( pos_unate_lits.size() > 0 )
          {
            on_off_div = 1; /* use pos_lit */
            score_div = pos_unate_lits[0].score;
            if ( neg_unate_lits.size() > 0 && neg_unate_lits[0].score > pos_unate_lits[0].score )
            {
              on_off_div = 0; /* use neg_lit */
              score_div = neg_unate_lits[0].score;
            }
          }
          else if ( neg_unate_lits.size() > 0 )
          {
            on_off_div = 0; /* use neg_lit */
            score_div = neg_unate_lits[0].score;
          }

          if ( num_inserts > 3u )
          {
            if ( pos_unate_pairs.size() > 0 )
            {
              on_off_pair = 1; /* use pos_pair */
              score_pair = pos_unate_pairs[0].score;
              if ( neg_unate_pairs.size() > 0 && neg_unate_pairs[0].score > pos_unate_pairs[0].score )
              {
                on_off_pair = 0; /* use neg_pair */
                score_pair = neg_unate_pairs[0].score;
              }
            }
            else if ( neg_unate_pairs.size() > 0 )
            {
              on_off_pair = 0; /* use neg_pair */
              score_pair = neg_unate_pairs[0].score;
            }
          }
        } );

        if ( score_div > score_pair / 2 ) /* divide with a divisor */
        {
          /* if using pos_lit (on_off_div = 1), modify on-set and use an OR gate on top;
            if using neg_lit (on_off_div = 0), modify off-set and use an AND gate on top
          */
          uint32_t const lit = on_off_div ? pos_unate_lits[0].lit : neg_unate_lits[0].lit;
          call_with_stopwatch( st.time_divide, [&]() {
            on_off_sets[on_off_div] &= lit & 0x1 ? get_div( lit >> 1 ) : ~get_div( lit >> 1 );
          } );

          auto const res_remain_div = compute_function_rec( num_inserts - 1 );
          if ( res_remain_div )
          {
            auto const new_lit = index_list.add_and( ( lit ^ 0x1 ), *res_remain_div ^ on_off_div );
            return new_lit + on_off_div;
          }
        }
        else if ( score_pair > 0 ) /* divide with a pair */
        {
          fanin_pair const pair = on_off_pair ? pos_unate_pairs[0] : neg_unate_pairs[0];
          call_with_stopwatch( st.time_divide, [&]() {
            if constexpr ( static_params::use_xor )
            {
              if ( pair.lit1 > pair.lit2 ) /* XOR pair: ~(lit1 ^ lit2) = ~lit1 ^ lit2 */
              {
                on_off_sets[on_off_pair] &= ( pair.lit1 & 0x1 ? get_div( pair.lit1 >> 1 ) : ~get_div( pair.lit1 >> 1 ) ) ^ ( pair.lit2 & 0x1 ? ~get_div( pair.lit2 >> 1 ) : get_div( pair.lit2 >> 1 ) );
              }
              else /* AND pair: ~(lit1 & lit2) = ~lit1 | ~lit2 */
              {
                on_off_sets[on_off_pair] &= ( pair.lit1 & 0x1 ? get_div( pair.lit1 >> 1 ) : ~get_div( pair.lit1 >> 1 ) ) | ( pair.lit2 & 0x1 ? get_div( pair.lit2 >> 1 ) : ~get_div( pair.lit2 >> 1 ) );
              }
            }
            else /* AND pair: ~(lit1 & lit2) = ~lit1 | ~lit2 */
            {
              on_off_sets[on_off_pair] &= ( pair.lit1 & 0x1 ? get_div( pair.lit1 >> 1 ) : ~get_div( pair.lit1 >> 1 ) ) | ( pair.lit2 & 0x1 ? get_div( pair.lit2 >> 1 ) : ~get_div( pair.lit2 >> 1 ) );
            }
          } );

          auto const res_remain_pair = compute_function_rec( num_inserts - 2 );
          if ( res_remain_pair )
          {
            uint32_t new_lit1;
            if constexpr ( static_params::use_xor )
            {
              new_lit1 = ( pair.lit1 > pair.lit2 ) ? index_list.add_xor( pair.lit1, pair.lit2 ) : index_list.add_and( pair.lit1, pair.lit2 );
            }
            else
            {
              new_lit1 = index_list.add_and( pair.lit1, pair.lit2 );
            }
            auto const new_lit2 = index_list.add_and( new_lit1 ^ 0x1, *res_remain_pair ^ on_off_pair );
            return new_lit2 + on_off_pair;
          }
        }
      }
    }

    if( static_params::use_spfd_synthesis )
    {
      index_list_t copy_id_list = index_list;
      _past_supports.clear();

      if( static_params::use_greedy_support )
      {
        for( uint32_t i{0}; i<static_params::max_num_support_samplings; ++i )
        {
          RNG.seed(i);
          auto supp = find_support_greedy( i );
          if( supp )
          {
            _support = *supp;
            auto syn = find_resynthesis( *supp, num_inserts );
            if( syn )
            {
              return syn;
            }
          }
        }
      }
      else if( static_params::use_enum )
      {
        auto suppG = find_support_greedy( 0 );
        if( suppG )
        {
          _support = *suppG;
          auto synG = find_resynthesis( *suppG, num_inserts );
          if( synG )
          {
            return synG;
          }
        }


        if( divisors.size() < 4 ) return std::nullopt;

        std::array<uint32_t, 4u> ref3 = {0,1,2,2};
        std::vector<uint32_t> supp3 {0, 0, 0,0};

        for( uint32_t i{0}; i<static_params::max_num_support_samplings; ++i )
        {
          RNG.seed(i);

          while( find_next_support4( ref3, supp3 ) )
          {
//            printf("Se) ");
//            for( auto x : supp3 )
//            {
//              printf("%d ", x );
//            }
//            printf("\n");
            auto synE = find_resynthesis( supp3, num_inserts );
            if( synE )
            {
              return synE;
            }
          }
        }
        index_list = copy_id_list;
      }
      else if( static_params::use_boltz )
      {
        for( uint32_t i{0}; i<static_params::max_num_support_samplings; ++i )
        {
          RNG.seed(i);
          if( i == 0 || _support.size() == 0)
          {
            auto supp = find_support_greedy( i );

            if( supp )
            {
//              printf("Se) ");
//              for( auto x : *supp )
//              {
//                printf("%d ", x );
//              }
//              printf("\n");
              _support = *supp;
              auto syn = find_resynthesis( *supp, num_inserts );
              if( syn )
              {
                return syn;
              }
            }
          }
          else
          {

            std::vector<uint32_t> partial_support = _support;
            // randomly erase an element an element
            int nerase=2;
            std::vector<uint32_t> erased;
            while( partial_support.size() > 0 && (partial_support.size()+nerase > _support.size()) )
            {
              std::uniform_int_distribution<> distrib(0, partial_support.size()-1);
              int idx = distrib(RNG);
              erased.push_back( partial_support[idx] );
              partial_support.erase( partial_support.begin() + idx );
            }

            auto supp = find_support_boltz( {}, i, {} );

            if( supp )
            {
    //          printf("Se) ");
    //          for( auto x : *supp )
    //          {
    //            printf("%d ", x );
    //          }
    //          printf("\n");
              _support = *supp;
              auto syn = find_resynthesis( *supp, num_inserts );
              if( syn )
              {
                return syn;
              }
            }
          }
        }
      }
    }




    return std::nullopt;
  }

  
  /* See if there is a Boolean cut and a candidate resynthesis function
     1. Perform iterative greedy support selection
     2. Extract the local functionality
     3. Resynthesize using Boolean matching with don't cares
   */
  std::optional<uint32_t> find_resynthesis( std::vector<uint32_t> const& supp, uint32_t max_num_gates )
  {
    index_list_t index_list_copy = index_list;

    if( supp.size() == 0 || supp.size() > static_params::max_support_size ) return std::nullopt;
    if( static_params::try_boolean_matching )
    {
      if( supp.size() > 4u )
      {
        for( int iTry{0}; iTry<static_params::max_resyn_attempts; ++iTry )
        {
          index_list_copy=index_list;
          auto [funcK, careK] = extract_functionalityK_from_signatures( supp );
          if( find_spfd_remapping( supp, funcK, careK, max_num_gates ) )
          {
            auto [lits4, func4, care4] = extract_functionality4_from_Kdivs( funcK, careK );
            const auto res = find_boolean_matching( lits4, func4, care4, max_num_gates ); 
            if( res && index_list.num_gates() <= max_num_gates )
            {
              return *res;
            }
            else
              index_list = index_list_copy;
          }
        }
      }
      else
      {
        auto [func4, care4] = extract_functionality4_from_signatures( supp );
        std::array<uint32_t, 4> lits = compute_literals( supp );
        const auto res =  find_boolean_matching( lits, func4, care4, max_num_gates ); 
        if( res && index_list.num_gates() <= max_num_gates )
        {
          return res;
        }
        else
          index_list = index_list_copy;
      }
      return std::nullopt;
    }
    else       
    {
      if( supp.size() == 0u )  return std::nullopt;
      auto [funcK, careK] = extract_functionalityK_from_signatures( supp );
      const auto res = find_spfd_resynthesis( supp, funcK, careK, max_num_gates );
      if( res && index_list.num_gates() <= max_num_gates )
      {
        return *res;
      }
      else
        index_list = index_list_copy;
    }
    
    index_list = index_list_copy;
    return std::nullopt;
  }

  bool find_spfd_remapping( std::vector<uint32_t> supp, truth_tableK_t const& funcK, truth_tableK_t const& careK, uint32_t max_num_gates )
  {
    _divsK.clear();
    _divsK.set_target( funcK, careK );
    _divsK.set_support( supp, _xsK );

    while( _divsK.size() > 4 && index_list.num_gates() <= max_num_gates )
    {
      if( !_divsK.update( index_list, _functional_library, max_num_gates ) )  return false;
    }
    return _divsK.size() <= 4 ;
    
  }

  std::tuple<std::array<uint32_t, 4>, truth_table4_t, truth_table4_t> extract_functionality4_from_Kdivs( truth_tableK_t const& funcK, truth_tableK_t const& careK )
  {
    if( _divsK.size() > 4 ) printf("[w] divisors size exceeds the limit \n");
    std::array<uint32_t, 4> lits {0};

    for( int i{0}; i<_divsK.size(); ++i )
      lits[i] = _divsK[i].lit;
    
    truth_table4_t func4;
    truth_table4_t care4;
    truth_table4_t temp4;
    truth_tableK_t temp = _divsK[0].func.construct();

    for( uint32_t m{0u}; m < 16u; ++m )
    {
      if( m < ( 1 << _divsK.size() ) )
      {
        temp = temp | ~temp;
        temp4 = temp4 | ~temp4;

        for( uint32_t l{0u}; l < _divsK.size(); ++l )
        {
          if( ( m >> l ) & 0x1 == 0x1 )
          {
            temp &= _divsK[l].func;
            temp4 &= _xs4[l];
          }
          else
          {
            temp &= ~_divsK[l].func;
            temp4 &= ~_xs4[l];
          }
        }

        if( kitty::count_ones( temp & careK ) > 0 ) // care value
        {
          care4 |= temp4;

          if( kitty::count_ones( temp & funcK ) > 0 )
          {
            func4 |= temp4;
          }
        }
      }
      else
        kitty::clear_bit( care4, m );
    }

    return std::make_tuple( lits, func4, care4 );
  }

  std::array<uint32_t, 4> compute_literals( std::vector<uint32_t> supp )
  {
    std::array<uint32_t, 4> lits {0};
    for( int i{0}; i<supp.size(); ++i )
    {
      lits[i] = supp[i] << 1u;
    }
    return lits;
  }


  /* See if there is a constant or divisor covering all on-set bits or all off-set bits.
     1. Check constant-resub
     2. Collect unate literals
     3. Find 0-resub (both positive unate and negative unate) and collect binate (neither pos nor neg unate) divisors
   */
  std::optional<uint32_t> find_one_unate()
  {
    num_bits[0] = kitty::count_ones( on_off_sets[0] ); /* off-set */
    num_bits[1] = kitty::count_ones( on_off_sets[1] ); /* on-set */
    if ( num_bits[0] == 0 )
    {
      return 1;
    }
    if ( num_bits[1] == 0 )
    {
      return 0;
    }

    for ( auto v = 1u; v < divisors.size(); ++v )
    {
      bool unateness[4] = { false, false, false, false };
      /* check intersection with off-set */
      if ( kitty::intersection_is_empty<TT, 1, 1>( get_div( v ), on_off_sets[0] ) )
      {
        pos_unate_lits.emplace_back( v << 1 );
        unateness[0] = true;
      }
      else if ( kitty::intersection_is_empty<TT, 0, 1>( get_div( v ), on_off_sets[0] ) )
      {
        pos_unate_lits.emplace_back( v << 1 | 0x1 );
        unateness[1] = true;
      }

      /* check intersection with on-set */
      if ( kitty::intersection_is_empty<TT, 1, 1>( get_div( v ), on_off_sets[1] ) )
      {
        neg_unate_lits.emplace_back( v << 1 );
        unateness[2] = true;
      }
      else if ( kitty::intersection_is_empty<TT, 0, 1>( get_div( v ), on_off_sets[1] ) )
      {
        neg_unate_lits.emplace_back( v << 1 | 0x1 );
        unateness[3] = true;
      }

      /* 0-resub */
      if ( unateness[0] && unateness[3] )
      {
        return ( v << 1 );
      }
      if ( unateness[1] && unateness[2] )
      {
        return ( v << 1 ) + 1;
      }
      /* useless unate literal */
      if ( ( unateness[0] && unateness[2] ) || ( unateness[1] && unateness[3] ) )
      {
        pos_unate_lits.pop_back();
        neg_unate_lits.pop_back();
      }
      /* binate divisor */
      else if ( !unateness[0] && !unateness[1] && !unateness[2] && !unateness[3] )
      {
        binate_divs.emplace_back( v );
      }
    }
    return std::nullopt;
  }

  /* Sort the unate literals by the number of minterms in the intersection.
     - For `pos_unate_lits`, `on_off` = 1, sort by intersection with on-set;
     - For `neg_unate_lits`, `on_off` = 0, sort by intersection with off-set
   */
  void sort_unate_lits( std::vector<unate_lit>& unate_lits, uint32_t on_off )
  {
    for ( auto& l : unate_lits )
    {
      l.score = kitty::count_ones( ( l.lit & 0x1 ? ~get_div( l.lit >> 1 ) : get_div( l.lit >> 1 ) ) & on_off_sets[on_off] );
    }
    std::sort( unate_lits.begin(), unate_lits.end(), [&]( unate_lit const& l1, unate_lit const& l2 ) {
      return l1.score > l2.score; // descending order
    } );
  }

  void sort_unate_pairs( std::vector<fanin_pair>& unate_pairs, uint32_t on_off )
  {
    for ( auto& p : unate_pairs )
    {
      if constexpr ( static_params::use_xor )
      {
        p.score = ( p.lit1 > p.lit2 ) ? kitty::count_ones( ( ( p.lit1 & 0x1 ? ~get_div( p.lit1 >> 1 ) : get_div( p.lit1 >> 1 ) ) ^ ( p.lit2 & 0x1 ? ~get_div( p.lit2 >> 1 ) : get_div( p.lit2 >> 1 ) ) ) & on_off_sets[on_off] )
                                      : kitty::count_ones( ( p.lit1 & 0x1 ? ~get_div( p.lit1 >> 1 ) : get_div( p.lit1 >> 1 ) ) & ( p.lit2 & 0x1 ? ~get_div( p.lit2 >> 1 ) : get_div( p.lit2 >> 1 ) ) & on_off_sets[on_off] );
      }
      else
      {
        p.score = kitty::count_ones( ( p.lit1 & 0x1 ? ~get_div( p.lit1 >> 1 ) : get_div( p.lit1 >> 1 ) ) & ( p.lit2 & 0x1 ? ~get_div( p.lit2 >> 1 ) : get_div( p.lit2 >> 1 ) ) & on_off_sets[on_off] );
      }
    }
    std::sort( unate_pairs.begin(), unate_pairs.end(), [&]( fanin_pair const& p1, fanin_pair const& p2 ) {
      return p1.score > p2.score; // descending order
    } );
  }

  /* See if there are two unate divisors covering all on-set bits or all off-set bits.
     - For `pos_unate_lits`, `on_off` = 1, try covering all on-set bits by combining two with an OR gate;
     - For `neg_unate_lits`, `on_off` = 0, try covering all off-set bits by combining two with an AND gate
   */
  std::optional<uint32_t> find_div_div( std::vector<unate_lit>& unate_lits, uint32_t on_off )
  {
    for ( auto i = 0u; i < unate_lits.size(); ++i )
    {
      uint32_t const& lit1 = unate_lits[i].lit;
      if ( unate_lits[i].score * 2 < num_bits[on_off] )
      {
        break;
      }
      for ( auto j = i + 1; j < unate_lits.size(); ++j )
      {
        uint32_t const& lit2 = unate_lits[j].lit;
        if ( unate_lits[i].score + unate_lits[j].score < num_bits[on_off] )
        {
          break;
        }
        auto const ntt1 = lit1 & 0x1 ? get_div( lit1 >> 1 ) : ~get_div( lit1 >> 1 );
        auto const ntt2 = lit2 & 0x1 ? get_div( lit2 >> 1 ) : ~get_div( lit2 >> 1 );
        if ( kitty::intersection_is_empty( ntt1, ntt2, on_off_sets[on_off] ) )
        {
          auto const new_lit = index_list.add_and( ( lit1 ^ 0x1 ), ( lit2 ^ 0x1 ) );
          return new_lit + on_off;
        }
      }
    }
    return std::nullopt;
  }

  std::optional<uint32_t> find_div_pair( std::vector<unate_lit>& unate_lits, std::vector<fanin_pair>& unate_pairs, uint32_t on_off )
  {
    for ( auto i = 0u; i < unate_lits.size(); ++i )
    {
      uint32_t const& lit1 = unate_lits[i].lit;
      for ( auto j = 0u; j < unate_pairs.size(); ++j )
      {
        fanin_pair const& pair2 = unate_pairs[j];
        if ( unate_lits[i].score + pair2.score < num_bits[on_off] )
        {
          break;
        }
        auto const ntt1 = lit1 & 0x1 ? get_div( lit1 >> 1 ) : ~get_div( lit1 >> 1 );
        TT ntt2;
        if constexpr ( static_params::use_xor )
        {
          if ( pair2.lit1 > pair2.lit2 )
          {
            ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) ) ^ ( pair2.lit2 & 0x1 ? ~get_div( pair2.lit2 >> 1 ) : get_div( pair2.lit2 >> 1 ) );
          }
          else
          {
            ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) ) | ( pair2.lit2 & 0x1 ? get_div( pair2.lit2 >> 1 ) : ~get_div( pair2.lit2 >> 1 ) );
          }
        }
        else
        {
          ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) ) | ( pair2.lit2 & 0x1 ? get_div( pair2.lit2 >> 1 ) : ~get_div( pair2.lit2 >> 1 ) );
        }

        if ( kitty::intersection_is_empty( ntt1, ntt2, on_off_sets[on_off] ) )
        {
          uint32_t new_lit1;
          if constexpr ( static_params::use_xor )
          {
            if ( pair2.lit1 > pair2.lit2 )
            {
              new_lit1 = index_list.add_xor( pair2.lit1, pair2.lit2 );
            }
            else
            {
              new_lit1 = index_list.add_and( pair2.lit1, pair2.lit2 );
            }
          }
          else
          {
            new_lit1 = index_list.add_and( pair2.lit1, pair2.lit2 );
          }
          auto const new_lit2 = index_list.add_and( ( lit1 ^ 0x1 ), new_lit1 ^ 0x1 );
          return new_lit2 + on_off;
        }
      }
    }
    return std::nullopt;
  }

  std::optional<uint32_t> find_pair_pair( std::vector<fanin_pair>& unate_pairs, uint32_t on_off )
  {
    for ( auto i = 0u; i < unate_pairs.size(); ++i )
    {
      fanin_pair const& pair1 = unate_pairs[i];
      if ( pair1.score * 2 < num_bits[on_off] )
      {
        break;
      }
      for ( auto j = i + 1; j < unate_pairs.size(); ++j )
      {
        fanin_pair const& pair2 = unate_pairs[j];
        if ( pair1.score + pair2.score < num_bits[on_off] )
        {
          break;
        }
        TT ntt1, ntt2;
        if constexpr ( static_params::use_xor )
        {
          if ( pair1.lit1 > pair1.lit2 )
          {
            ntt1 = ( pair1.lit1 & 0x1 ? get_div( pair1.lit1 >> 1 ) : ~get_div( pair1.lit1 >> 1 ) ) ^ ( pair1.lit2 & 0x1 ? ~get_div( pair1.lit2 >> 1 ) : get_div( pair1.lit2 >> 1 ) );
          }
          else
          {
            ntt1 = ( pair1.lit1 & 0x1 ? get_div( pair1.lit1 >> 1 ) : ~get_div( pair1.lit1 >> 1 ) ) | ( pair1.lit2 & 0x1 ? get_div( pair1.lit2 >> 1 ) : ~get_div( pair1.lit2 >> 1 ) );
          }
          if ( pair2.lit1 > pair2.lit2 )
          {
            ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) ) ^ ( pair2.lit2 & 0x1 ? ~get_div( pair2.lit2 >> 1 ) : get_div( pair2.lit2 >> 1 ) );
          }
          else
          {
            ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) ) | ( pair2.lit2 & 0x1 ? get_div( pair2.lit2 >> 1 ) : ~get_div( pair2.lit2 >> 1 ) );
          }
        }
        else
        {
          ntt1 = ( pair1.lit1 & 0x1 ? get_div( pair1.lit1 >> 1 ) : ~get_div( pair1.lit1 >> 1 ) ) | ( pair1.lit2 & 0x1 ? get_div( pair1.lit2 >> 1 ) : ~get_div( pair1.lit2 >> 1 ) );
          ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) ) | ( pair2.lit2 & 0x1 ? get_div( pair2.lit2 >> 1 ) : ~get_div( pair2.lit2 >> 1 ) );
        }

        if ( kitty::intersection_is_empty( ntt1, ntt2, on_off_sets[on_off] ) )
        {
          uint32_t fanin_lit1, fanin_lit2;
          if constexpr ( static_params::use_xor )
          {
            if ( pair1.lit1 > pair1.lit2 )
            {
              fanin_lit1 = index_list.add_xor( pair1.lit1, pair1.lit2 );
            }
            else
            {
              fanin_lit1 = index_list.add_and( pair1.lit1, pair1.lit2 );
            }
            if ( pair2.lit1 > pair2.lit2 )
            {
              fanin_lit2 = index_list.add_xor( pair2.lit1, pair2.lit2 );
            }
            else
            {
              fanin_lit2 = index_list.add_and( pair2.lit1, pair2.lit2 );
            }
          }
          else
          {
            fanin_lit1 = index_list.add_and( pair1.lit1, pair1.lit2 );
            fanin_lit2 = index_list.add_and( pair2.lit1, pair2.lit2 );
          }
          uint32_t const output_lit = index_list.add_and( fanin_lit1 ^ 0x1, fanin_lit2 ^ 0x1 );
          return output_lit + on_off;
        }
      }
    }
    return std::nullopt;
  }

  std::optional<uint32_t> find_xor()
  {
    /* collect XOR-type pairs (d1 ^ d2) & off = 0 or ~(d1 ^ d2) & on = 0, selecting d1, d2 from binate_divs */
    for ( auto i = 0u; i < binate_divs.size(); ++i )
    {
      for ( auto j = i + 1; j < binate_divs.size(); ++j )
      {
        auto const tt_xor = get_div( binate_divs[i] ) ^ get_div( binate_divs[j] );
        bool unateness[4] = { false, false, false, false };
        /* check intersection with off-set; additionally check intersection with on-set is not empty (otherwise it's useless) */
        if ( kitty::intersection_is_empty<TT, 1, 1>( tt_xor, on_off_sets[0] ) && !kitty::intersection_is_empty<TT, 1, 1>( tt_xor, on_off_sets[1] ) )
        {
          pos_unate_pairs.emplace_back( binate_divs[i] << 1, binate_divs[j] << 1, true );
          unateness[0] = true;
        }
        if ( kitty::intersection_is_empty<TT, 0, 1>( tt_xor, on_off_sets[0] ) && !kitty::intersection_is_empty<TT, 0, 1>( tt_xor, on_off_sets[1] ) )
        {
          pos_unate_pairs.emplace_back( ( binate_divs[i] << 1 ) + 1, binate_divs[j] << 1, true );
          unateness[1] = true;
        }

        /* check intersection with on-set; additionally check intersection with off-set is not empty (otherwise it's useless) */
        if ( kitty::intersection_is_empty<TT, 1, 1>( tt_xor, on_off_sets[1] ) && !kitty::intersection_is_empty<TT, 1, 1>( tt_xor, on_off_sets[0] ) )
        {
          neg_unate_pairs.emplace_back( binate_divs[i] << 1, binate_divs[j] << 1, true );
          unateness[2] = true;
        }
        if ( kitty::intersection_is_empty<TT, 0, 1>( tt_xor, on_off_sets[1] ) && !kitty::intersection_is_empty<TT, 0, 1>( tt_xor, on_off_sets[0] ) )
        {
          neg_unate_pairs.emplace_back( ( binate_divs[i] << 1 ) + 1, binate_divs[j] << 1, true );
          unateness[3] = true;
        }

        if ( unateness[0] && unateness[2] )
        {
          return index_list.add_xor( ( binate_divs[i] << 1 ), ( binate_divs[j] << 1 ) );
        }
        if ( unateness[1] && unateness[3] )
        {
          return index_list.add_xor( ( binate_divs[i] << 1 ) + 1, ( binate_divs[j] << 1 ) );
        }
      }
    }

    return std::nullopt;
  }

  /* collect AND-type pairs (d1 & d2) & off = 0 or ~(d1 & d2) & on = 0, selecting d1, d2 from binate_divs */
  void collect_unate_pairs()
  {
    for ( auto i = 0u; i < binate_divs.size(); ++i )
    {
      for ( auto j = i + 1; j < binate_divs.size(); ++j )
      {
        collect_unate_pairs_detail<1, 1>( binate_divs[i], binate_divs[j] );
        collect_unate_pairs_detail<0, 1>( binate_divs[i], binate_divs[j] );
        collect_unate_pairs_detail<1, 0>( binate_divs[i], binate_divs[j] );
        collect_unate_pairs_detail<0, 0>( binate_divs[i], binate_divs[j] );
      }
    }
  }

  template<bool pol1, bool pol2>
  void collect_unate_pairs_detail( uint32_t div1, uint32_t div2 )
  {
    /* check intersection with off-set; additionally check intersection with on-set is not empty (otherwise it's useless) */
    if ( kitty::intersection_is_empty<TT, pol1, pol2>( get_div( div1 ), get_div( div2 ), on_off_sets[0] ) && !kitty::intersection_is_empty<TT, pol1, pol2>( get_div( div1 ), get_div( div2 ), on_off_sets[1] ) )
    {
      pos_unate_pairs.emplace_back( ( div1 << 1 ) + (uint32_t)( !pol1 ), ( div2 << 1 ) + (uint32_t)( !pol2 ) );
    }
    /* check intersection with on-set; additionally check intersection with off-set is not empty (otherwise it's useless) */
    else if ( kitty::intersection_is_empty<TT, pol1, pol2>( get_div( div1 ), get_div( div2 ), on_off_sets[1] ) && !kitty::intersection_is_empty<TT, pol1, pol2>( get_div( div1 ), get_div( div2 ), on_off_sets[0] ) )
    {
      neg_unate_pairs.emplace_back( ( div1 << 1 ) + (uint32_t)( !pol1 ), ( div2 << 1 ) + (uint32_t)( !pol2 ) );
    }
  }

  #pragma region support_sampling

  std::optional<std::vector<uint32_t>> find_support( uint32_t iTry )
  {
    if( static_params::use_greedy_support )
    {
      auto supp = find_support_greedy( iTry );
      if( supp )
      {
        return *supp;
      }
      return std::nullopt;
    }
    else
    {
      auto supp = find_support_greedy( iTry );
      if( supp )
      {
        return *supp;
      }

      auto supps = enumerate_boolean_cuts();
      if( supps.size() > 0 )
      {
        return supps[0];
      }
      return std::nullopt;
    }

    return std::nullopt;
  }

  bool find_next_support4( uint32_t & ref0, uint32_t & ref1, uint32_t & ref2, uint32_t & ref3, std::vector<uint32_t> & supp4 )
  {
    if( ref3 + 1 < _scored_divs.size() )
    {
      ref3++;
    }
    else if( ref2 + 2 < _scored_divs.size() )
    {
      ref2++;
      ref3=ref2+1;
    }
    else if( ref1 + 3 < _scored_divs.size() )
    {
      ref1++;
      ref2=ref1+1;
      ref3=ref2+1;
    }
    else if( ref0 + 4 < _scored_divs.size() )
    {
      ref0++;
      ref1=ref0+1;
      ref2=ref1+1;
      ref3=ref2+1;
    }
    else
    {
      return false;
    }

    if( ref0 > _scored_divs.size() || ref1 > _scored_divs.size() || ref2 > _scored_divs.size() || ref3 > _scored_divs.size() ) return false;
    _gSPFD.reset();

    std::array<TT,2u>  masks0;
    std::array<TT,4u>  masks1;
    std::array<TT,8u>  masks2;
    std::array<TT,16u> masks3;
    std::array<bool,16> is_killed;
    uint32_t nKills0;
    uint32_t nKills1;
    uint32_t nKills2;
    uint32_t nKills3;

    bool is_valid;
    for( int i0{ref0}; i0<_scored_divs.size(); ++i0 )
    {
      nKills0 = 0;
      masks0[0] = _gSPFD.care & get_div( _scored_divs[i0].div );
      masks0[1] = _gSPFD.care & ~get_div( _scored_divs[i0].div );

      if( kitty::is_const0(masks0[0]) ) { is_killed[0] = true; nKills0++; } else { is_killed[0] = false; }
      if( kitty::is_const0(masks0[1]) ) { is_killed[1] = true; nKills0++; } else { is_killed[1] = false; }
      if( nKills0 == 2u ) continue;
      //if( _scored_divs[i0].cost > _gSPFD.nEdges/2 ) return false;

      for( int i1{ref1}; i1<_scored_divs.size(); ++i1 )
      {
        //if( _scored_divs[i0].cost + _scored_divs[i1].cost > _gSPFD.nEdges/2 ) return false; 

        for( int k1{0}; k1<2; k1++ )
        {
          masks1[k1] = masks0[k1] &  get_div( _scored_divs[i1].div );
          masks1[k1+2] = masks0[k1] & ~get_div( _scored_divs[i1].div );
          if( is_killed[k1] )
          { 
            is_killed[k1+2] = true; 
            nKills1+=2; 
          } 
          else 
          { 
            if( kitty::is_const0(masks1[k1]) ) { is_killed[k1] = true; nKills1++; } else { is_killed[k1] = false; }
            if( kitty::is_const0(masks1[k1+2]) ) { is_killed[k1+2] = true; nKills1++; } else { is_killed[k1+2] = false; }
          }
        }
        if( nKills1 == 4u ) continue;

        for( int i2{ref2}; i2<_scored_divs.size(); ++i2 )
        {            
          for( int k2{0}; k2<4; k2++ )
          {
            masks2[k2] = masks1[k2] &  get_div( _scored_divs[i2].div );
            masks2[k2+4] = masks1[k2] & ~get_div( _scored_divs[i2].div );
            if( is_killed[k2] )
            { 
              is_killed[k2+4] = true; 
              nKills2+=2; 
            } 
            else 
            { 
              if( kitty::is_const0(masks2[k2]) ) { is_killed[k2] = true; nKills2++; } else { is_killed[k2] = false; }
              if( kitty::is_const0(masks2[k2+4]) ) { is_killed[k2+4] = true; nKills2++; } else { is_killed[k2+4] = false; }
            }
          }
          if( nKills2 == 8u ) continue;

          for( int i3{ref3}; i3<_scored_divs.size(); ++i3 )
          {

            for( int k3{0}; k3<8; k3++ )
            {
              masks3[k3] = masks2[k3] &  get_div( _scored_divs[i3].div );
              masks3[k3+8] = masks2[k3] & ~get_div( _scored_divs[i3].div );
              if( is_killed[k3] )
              { 
                is_killed[k3+8] = true; 
                nKills3+=2; 
              } 
              else 
              { 
                if( kitty::is_const0(masks3[k3]) ) { is_killed[k3] = true; nKills3++; } else { is_killed[k3] = false; }
                if( kitty::is_const0(masks3[k3+8]) ) { is_killed[k3+2] = true; nKills3++; } else { is_killed[k3+8] = false; }
              }
            }
            if( nKills3 == 16u ) continue;

            if( _scored_divs[i3].cost + _scored_divs[i2].cost + _scored_divs[i1].cost + _scored_divs[i0].cost > _gSPFD.nEdges ) return false;

            is_valid = true;
            for( uint32_t m{0}; m<16; ++m )
            {
              is_valid = kitty::is_const0( masks3[m] & _gSPFD.func[1] ) || kitty::equal( masks3[m] & _gSPFD.func[1], masks3[m] );
              if( !is_valid ) break;
            }
            if( is_valid )
            {
              ref0 = i0;
              ref1 = i1;
              ref2 = i2;
              ref3 = i3;
              supp4[0] = _scored_divs[i0].div;
              supp4[1] = _scored_divs[i1].div;
              supp4[2] = _scored_divs[i2].div;
              supp4[3] = _scored_divs[i3].div;
              return true;
            }
          }
        } 
      } 
    }
    return false;
  }

  bool find_next_support4( std::array<uint32_t, 4> & ref4, std::vector<uint32_t> & supp4 )
  {
    uint64_t key;

    if( ref4[3] + 1 < _scored_divs.size() )
    {
      ref4[3]++;
    }
    else if( ref4[2] + 2 < _scored_divs.size() )
    {
      ref4[2]++;
      ref4[3]=ref4[1]+1;
    }
    else if( ref4[1] + 3 < _scored_divs.size() )
    {
      ref4[1]++;
      ref4[2]=ref4[1]+1;
      ref4[3]=ref4[2]+1;
    }
    else if( ref4[0] + 4 < _scored_divs.size() )
    {
      ref4[0]++;
      ref4[1]=ref4[0]+1;
      ref4[2]=ref4[1]+1;
      ref4[3]=ref4[2]+1;
    }
    else
    {
      return false;
    }

    if( ref4[0] > _scored_divs.size() || ref4[1] > _scored_divs.size() || ref4[2] > _scored_divs.size() || ref4[3] > _scored_divs.size() ) return false;
    _gSPFD.reset();

    std::array<TT,2u>  masks0;
    std::array<TT,4u>  masks1;
    std::array<TT,8u>  masks2;
    std::array<TT,16u>  masks3;
    std::array<bool,16> is_killed;
    uint32_t nKills0;
    uint32_t nKills1;
    uint32_t nKills2;
    uint32_t nKills3;

    bool is_valid;
    for( int i0{ref4[0]}; i0<_scored_divs.size(); ++i0 )
    {
      if( _scored_divs[i0].cost > 0.65 ) break;

      uint64_t key0=i0;

      nKills0 = 0;
      masks0[0] = _gSPFD.care & get_div( _scored_divs[i0].div );
      masks0[1] = _gSPFD.care & ~get_div( _scored_divs[i0].div );

      if( kitty::is_const0(masks0[0]) ) { is_killed[0] = true; nKills0++; } else { is_killed[0] = false; }
      if( kitty::is_const0(masks0[1]) ) { is_killed[1] = true; nKills0++; } else { is_killed[1] = false; }
      if( nKills0 == 2u ) continue;

      for( int i1{ref4[1]}; i1<_scored_divs.size(); ++i1 )
      {
        if( _scored_divs[i1].cost > 0.65 ) break;

        //if( _scored_divs[i0].cost + _scored_divs[i1].cost > _gSPFD.nEdges/2 ) return false; 
        uint64_t key1 = i1;
        for( int k1{0}; k1<2; k1++ )
        {
          masks1[k1] = masks0[k1] &  get_div( _scored_divs[i1].div );
          masks1[k1+2] = masks0[k1] & ~get_div( _scored_divs[i1].div );
          if( is_killed[k1] )
          { 
            is_killed[k1+2] = true; 
            nKills1+=2; 
          } 
          else 
          { 
            if( kitty::is_const0(masks1[k1]) ) { is_killed[k1] = true; nKills1++; } else { is_killed[k1] = false; }
            if( kitty::is_const0(masks1[k1+2]) ) { is_killed[k1+2] = true; nKills1++; } else { is_killed[k1+2] = false; }
          }
        }
        if( nKills1 == 4u ) continue;

        for( int i2{ref4[2]}; i2<_scored_divs.size(); ++i2 )
        {    
          if( _scored_divs[i2].cost > 0.65 )  break ;

          for( int k2{0}; k2<4; k2++ )
          {
            masks2[k2] = masks1[k2] &  get_div( _scored_divs[i2].div );
            masks2[k2+4] = masks1[k2] & ~get_div( _scored_divs[i2].div );
            if( is_killed[k2] )
            { 
              is_killed[k2+4] = true; 
              nKills2+=2; 
            } 
            else 
            { 
              if( kitty::is_const0(masks2[k2]) ) { is_killed[k2] = true; nKills2++; } else { is_killed[k2] = false; }
              if( kitty::is_const0(masks2[k2+4]) ) { is_killed[k2+4] = true; nKills2++; } else { is_killed[k2+4] = false; }
            }
          }
          if( nKills2 == 8u ) continue;

          for( int i3{ref4[3]}; i3<_scored_divs.size(); ++i3 )
          {    
            if( _scored_divs[i3].cost > 0.65 )  break ;

            uint64_t key3 = i3;
            key = ( key3 << 60u ) | ( key0 << 40u ) | ( key1 << 20u ) | ( key0 );
            uint32_t count;
            if( auto search = _past_supports_count.find( key ); search != _past_supports_count.end() )
            {
              count = search->second;
              _past_supports_count[key]++;
            }
            if( _past_supports_count[key] > 16 )
              continue;

            for( int k2{0}; k2<8; k2++ )
            {
              masks3[k2] = masks2[k2] &  get_div( _scored_divs[i3].div );
              masks3[k2+8] = masks2[k2] & ~get_div( _scored_divs[i3].div );
              if( is_killed[k2] )
              { 
                is_killed[k2+8] = true; 
                nKills3+=2; 
              } 
              else 
              { 
                if( kitty::is_const0(masks3[k2]) ) { is_killed[k2] = true; nKills3++; } else { is_killed[k2] = false; }
                if( kitty::is_const0(masks3[k2+8]) ) { is_killed[k2+8] = true; nKills3++; } else { is_killed[k2+8] = false; }
              }
            }
            if( nKills3 == 16u ) continue;

            if( _scored_divs[i2].cost + _scored_divs[i1].cost + _scored_divs[i0].cost + _scored_divs[i3].cost > _gSPFD.nEdges ) return false;

            is_valid = true;
            for( uint32_t m{0}; m<16; ++m )
            {
              is_valid = kitty::is_const0( masks3[m] & _gSPFD.func[1] ) || kitty::equal( masks3[m] & _gSPFD.func[1], masks3[m] );
              if( !is_valid ) break;
            }
            if( is_valid )
            {
              ref4[0] = i0;
              ref4[1] = i1;
              ref4[2] = i2;
              ref4[3] = i3;
              supp4[0] = _scored_divs[i0].div;
              supp4[1] = _scored_divs[i1].div;
              supp4[2] = _scored_divs[i2].div;
              supp4[3] = _scored_divs[i3].div;
              return true;
            }
          }
        } 
      } 
    }
    return false;
  }


  bool find_next_support3( std::array<uint32_t, 3> & ref3, std::vector<uint32_t> & supp3 )
  {
    uint64_t key;

    if( ref3[2] + 1 < _scored_divs.size() )
    {
      ref3[2]++;
    }
    else if( ref3[1] + 2 < _scored_divs.size() )
    {
      ref3[1]++;
      ref3[2]=ref3[1]+1;
    }
    else if( ref3[0] + 3 < _scored_divs.size() )
    {
      ref3[0]++;
      ref3[1]=ref3[0]+1;
      ref3[2]=ref3[1]+1;
    }
    else
    {
      return false;
    }

    if( ref3[0] > _scored_divs.size() || ref3[1] > _scored_divs.size() || ref3[2] > _scored_divs.size() ) return false;
    _gSPFD.reset();

    std::array<TT,2u>  masks0;
    std::array<TT,4u>  masks1;
    std::array<TT,8u>  masks2;
    std::array<bool,8> is_killed;
    uint32_t nKills0;
    uint32_t nKills1;
    uint32_t nKills2;

    bool is_valid;
    for( int i0{ref3[0]}; i0<_scored_divs.size(); ++i0 )
    {
      if( _scored_divs[i0].cost > 0.65 ) break;

      uint64_t key0=i0;

      nKills0 = 0;
      masks0[0] = _gSPFD.care & get_div( _scored_divs[i0].div );
      masks0[1] = _gSPFD.care & ~get_div( _scored_divs[i0].div );

      if( kitty::is_const0(masks0[0]) ) { is_killed[0] = true; nKills0++; } else { is_killed[0] = false; }
      if( kitty::is_const0(masks0[1]) ) { is_killed[1] = true; nKills0++; } else { is_killed[1] = false; }
      if( nKills0 == 2u ) continue;

      for( int i1{ref3[1]}; i1<_scored_divs.size(); ++i1 )
      {
        if( _scored_divs[i1].cost > 0.65 ) break;

        //if( _scored_divs[i0].cost + _scored_divs[i1].cost > _gSPFD.nEdges/2 ) return false; 
        uint64_t key1 = i1;
        for( int k1{0}; k1<2; k1++ )
        {
          masks1[k1] = masks0[k1] &  get_div( _scored_divs[i1].div );
          masks1[k1+2] = masks0[k1] & ~get_div( _scored_divs[i1].div );
          if( is_killed[k1] )
          { 
            is_killed[k1+2] = true; 
            nKills1+=2; 
          } 
          else 
          { 
            if( kitty::is_const0(masks1[k1]) ) { is_killed[k1] = true; nKills1++; } else { is_killed[k1] = false; }
            if( kitty::is_const0(masks1[k1+2]) ) { is_killed[k1+2] = true; nKills1++; } else { is_killed[k1+2] = false; }
          }
        }
        if( nKills1 == 4u ) continue;

        for( int i2{ref3[2]}; i2<_scored_divs.size(); ++i2 )
        {    
          if( _scored_divs[i2].cost > 0.65 )  break ;

          uint64_t key2 = i2;
          key = ( key0 << 40u ) | ( key1 << 20u ) | ( key0 );
          uint32_t count;
          if( auto search = _past_supports_count.find( key ); search != _past_supports_count.end() )
          {
            count = search->second;
            _past_supports_count[key]++;
          }
          if( _past_supports_count[key] > 16 )
            continue;

          for( int k2{0}; k2<4; k2++ )
          {
            masks2[k2] = masks1[k2] &  get_div( _scored_divs[i2].div );
            masks2[k2+4] = masks1[k2] & ~get_div( _scored_divs[i2].div );
            if( is_killed[k2] )
            { 
              is_killed[k2+4] = true; 
              nKills2+=2; 
            } 
            else 
            { 
              if( kitty::is_const0(masks2[k2]) ) { is_killed[k2] = true; nKills2++; } else { is_killed[k2] = false; }
              if( kitty::is_const0(masks2[k2+4]) ) { is_killed[k2+4] = true; nKills2++; } else { is_killed[k2+4] = false; }
            }
          }
          if( nKills2 == 8u ) continue;

          if( _scored_divs[i2].cost + _scored_divs[i1].cost + _scored_divs[i0].cost > _gSPFD.nEdges ) return false;

          is_valid = true;
          for( uint32_t m{0}; m<8; ++m )
          {
            is_valid = kitty::is_const0( masks2[m] & _gSPFD.func[1] ) || kitty::equal( masks2[m] & _gSPFD.func[1], masks2[m] );
            if( !is_valid ) break;
          }
          if( is_valid )
          {
            ref3[0] = i0;
            ref3[1] = i1;
            ref3[2] = i2;
            supp3[0] = _scored_divs[i0].div;
            supp3[1] = _scored_divs[i1].div;
            supp3[2] = _scored_divs[i2].div;
            return true;
          }
        } 
      } 
    }
    return false;
  }

  std::vector<std::vector<uint32_t>> enumerate_boolean_cuts()
  {
    _gSPFD.reset();
    std::vector<std::vector<uint32_t>> boolean_cuts;
    std::array<TT,8u> masks;
    bool is_valid;
    for( int i0{1}; i0<_scored_divs.size(); ++i0 )
    {
      for( int i1{1}; i1<_scored_divs.size(); ++i1 )
      {
        for( int i2{1}; i2<_scored_divs.size(); ++i2 )
        {
          if( _scored_divs[i2].cost + _scored_divs[i1].cost + _scored_divs[i2].cost < _gSPFD.nEdges ) break;
          is_valid = true;
          for( uint32_t m{0}; m<8; ++m )
          {
            masks[m] = _gSPFD.care;
            masks[m] = ( m & 0x1 == 0x0 ) ? ~get_div( i0 ) : get_div( i0 );
            masks[m] = ( (m>>1) & 0x1 == 0x0 ) ? ~get_div( i1 ) : get_div( i1 );
            masks[m] = ( (m>>2) & 0x1 == 0x0 ) ? ~get_div( i2 ) : get_div( i2 );

            is_valid = kitty::is_const0( masks[m] & _gSPFD.func[1] ) || kitty::equal( masks[m] & _gSPFD.func[1], masks[m] );
            if( !is_valid ) break;
          }
          if( is_valid )
          {
            //printf("found %d %d %d \n", i0, i1, i2 );
            boolean_cuts.emplace_back( std::vector<uint32_t>{i0, i1, i2} );
          }
        } 
      } 
    }
    return boolean_cuts;
  }

  std::optional<std::vector<uint32_t>> find_support_greedy( uint32_t iteration )
  {
    _level_1_divisors.clear();
    uint32_t cost, best_cost;
    std::vector<uint32_t> best_candidates;
    std::vector<uint32_t> supp;
    
    _uSPFD.reset();

    /* add recomputation of the support */

    while( !_uSPFD.is_covered() )
    {
      best_cost = std::numeric_limits<double>::max();
      if( _uSPFD.is_saturated() ) break;
      for( uint32_t iCnd{1}; iCnd<divisors.size(); ++iCnd )
      {
        cost = _uSPFD.evaluate( get_div( iCnd ) );
        if( cost < best_cost  )
        {
          best_cost = cost;
          best_candidates = {iCnd};
        }
        else if( cost == best_cost )
        {
          best_candidates.push_back( iCnd );
        }
      }
      if( best_candidates.size() == 0 ) break;

      std::uniform_int_distribution<> distrib(0, best_candidates.size()-1);
      int idx = distrib(RNG);
      supp.push_back( best_candidates[idx] );
      _level_1_divisors.insert( best_candidates[idx] );
      _uSPFD.update( get_div( best_candidates[idx] ) );
    }
    if( _uSPFD.is_covered() )
    {
      std::sort( supp.begin(), supp.end() );
      if( _past_supports.insert(supp).second )
      {
        return supp;
      }
    }
    return std::nullopt;
  }


  std::optional<std::vector<uint32_t>> find_support_boltz( std::vector<uint32_t> const& partial_support, int iteration , std::vector<uint32_t> erased={} )
  {
    double BETA = 1u << 11u;
    _gSPFD.reset();
    double min_cost;
    double max_cost;
    std::vector<double> costs;
    std::vector<uint32_t> supp;
    for( auto div : partial_support )
    {
      if( _gSPFD.is_covered() ) break;
      supp.push_back( div );
      _gSPFD.update( get_div( div ) );
    }
    /* add recomputation of the support */

    int iter{0};
    while( !_gSPFD.is_covered() )
    {
      BETA = 1 << ( static_params::max_support_size + 2 - iter );
      iter++;
      if( _gSPFD.is_saturated() ) break;
      costs.clear();
      costs.push_back(0);
      min_cost = std::numeric_limits<double>::max();
      max_cost = std::numeric_limits<double>::min();

      for( uint32_t iDiv{1}; iDiv<divisors.size(); ++iDiv )
      {
        costs.push_back(_gSPFD.evaluate( get_div( iDiv ) ));
        if( costs[iDiv] < min_cost ) min_cost = costs[iDiv];
        if( costs[iDiv] > max_cost ) max_cost = costs[iDiv];
      }

      for( uint32_t i{1}; i<costs.size(); ++i )
      {
        costs[i] = exp( -BETA*(costs[i]-min_cost)/(max_cost-min_cost) );
      }
      
      for( auto div : supp )
        costs[div] = 0;
      for( auto er : erased )
        costs[er] = 0;

      for( uint32_t i{1}; i<costs.size(); ++i )
      {
        costs[i] += costs[i-1];
      }

      std::uniform_real_distribution<> distrib(0, 1);
      double rnd = distrib(RNG);

      bool found{false};
      for( uint32_t i{1}; i<costs.size(); ++i )
      {
        if( rnd*costs.back() <= costs[i] )
        {
          supp.push_back(i);
          _gSPFD.update( get_div(i) );
          found = true;
          break;
        }
      }

      if( !found ) return std::nullopt;
    }
    if( _gSPFD.is_covered() )
    {
      std::sort( supp.begin(), supp.end() );
      if( _past_supports.insert(supp).second )
      {
        _support = supp;
        return supp;
      }
    }

    return std::nullopt;
  }


  std::optional<std::vector<uint32_t>> find_support_greedy_with_offset( uint32_t offset )
  {
    offset = offset % _scored_divs.size();

    uint32_t cost, best_cost;
    std::vector<uint32_t> best_candidates;
    std::vector<uint32_t> supp;
    
    _uSPFD.reset();
    //_scored_divs.print();

    supp = {_scored_divs[offset].div};
    _uSPFD.update( get_div( _scored_divs[offset].div ) );

    /* add recomputation of the support */

    while( !_uSPFD.is_covered() )
    {
      best_cost = std::numeric_limits<double>::max();
      if( _uSPFD.is_saturated() ) break;

      for( uint32_t iDiv{offset+1}; iDiv < _scored_divs.size(); ++iDiv )
      {
        uint32_t div = _scored_divs[iDiv].div;
        cost = _uSPFD.evaluate( get_div( div ) );
        if( cost < best_cost  )
        {
          best_cost = cost;
          best_candidates = {div};
        }
        else if( cost == best_cost )
        {
          best_candidates.push_back( div );
        }
      }
      if( best_candidates.size() == 0 ) break;

      std::uniform_int_distribution<> distrib(0, best_candidates.size()-1);
      int idx = distrib(RNG);
      supp.push_back( best_candidates[idx] );
      _uSPFD.update( get_div( best_candidates[idx] ) );
    }

    while( !_uSPFD.is_covered() )
    {
      best_cost = std::numeric_limits<double>::max();
      if( _uSPFD.is_saturated() ) break;

      for( uint32_t iDiv{0}; iDiv < offset; ++iDiv )
      {
        uint32_t div = _scored_divs[iDiv].div;
        cost = _uSPFD.evaluate( get_div( div ) );
        if( cost < best_cost  )
        {
          best_cost = cost;
          best_candidates = {div};
        }
        else if( cost == best_cost )
        {
          best_candidates.push_back( div );
        }
      }
      if( best_candidates.size() == 0 ) break;

      std::uniform_int_distribution<> distrib(0, best_candidates.size()-1);
      int idx = distrib(RNG);
      supp.push_back( best_candidates[idx] );
      _uSPFD.update( get_div( best_candidates[idx] ) );
    }

    if( _uSPFD.is_covered() )
    {
      std::sort( supp.begin(), supp.end() );
      if( _past_supports.insert(supp).second )
      {
        return supp;
      }
    }
    return std::nullopt;
  }

  std::optional<std::vector<uint32_t>> find_support_greedy_with_offset_s( uint32_t offset )
  {

    double BETA = 100;
    offset = offset % _scored_divs.size();
    double cost;
    std::vector<uint32_t> supp;
    std::vector<double> costs;
    double max_cost = std::numeric_limits<double>::max();
    double min_cost = std::numeric_limits<double>::min();
    
    _gSPFD.reset();
    //_scored_divs.print();

    supp = {_scored_divs[offset].div};
    _gSPFD.update( get_div( _scored_divs[offset].div ) );

    /* add recomputation of the support */
    bool found{false};

    while( !_gSPFD.is_covered() && supp.size() < static_params::max_support_size)
    {
      costs.clear();
      for( uint32_t i{0}; i<divisors.size(); ++i )
      {
        costs.push_back(0);
      }

      if( _gSPFD.is_saturated() ) break;

      for( uint32_t iDiv{offset+1}; iDiv < _scored_divs.size(); ++iDiv )
      {
        uint32_t div = _scored_divs[iDiv].div;
        cost = _gSPFD.evaluate( get_div( div ) );
        costs[div] = cost ;
        if( cost < min_cost ) min_cost = cost;
        if( cost > max_cost ) max_cost = cost;
      }

      for( uint32_t iDiv{0}; iDiv < _scored_divs.size(); ++iDiv )
      {
        uint32_t div = _scored_divs[iDiv].div;
        if( ( div > offset ) && (std::find( supp.begin(), supp.end(), div ) == supp.end()) )
        {
          //costs[div] = exp(-BETA*((costs[div]-min_cost)/(_gSPFD.nEdges-min_cost)));
          costs[div] = exp(-BETA*(costs[div]/(_gSPFD.nEdges)));
        }
        else
        {
          costs[div] = 0;
        }
      }

      for( uint32_t iDiv{1}; iDiv < divisors.size(); ++iDiv )
      {
        costs[iDiv]+=costs[iDiv-1];
      }

      std::uniform_real_distribution<> distrib(0, 1);
      double rnd = distrib(RNG);

      found=false;
      if( costs.back()>0 )
      {
        for( uint32_t i{1}; i<costs.size(); ++i )
        {
          if( rnd*costs.back() <= costs[i] )
          {
            supp.push_back(i);
            _gSPFD.update( get_div( i ) );
            found = true;
            break;
          }
        }
        
  //      for( auto x : supp)
  //        printf("%d ", x );
  //      printf(" : %f\n", _gSPFD.nEdges );
      }
      if( !found ) break;
    }

    while( found && !_gSPFD.is_covered() && supp.size() < static_params::max_support_size )
    {
      costs.clear();
      for( uint32_t i{0}; i<divisors.size(); ++i )
      {
        costs.push_back(0);
      }

      if( _gSPFD.is_saturated() ) break;

      for( uint32_t iDiv{offset+1}; iDiv < _scored_divs.size(); ++iDiv )
      {
        uint32_t div = _scored_divs[iDiv].div;
        cost = _gSPFD.evaluate( get_div( div ) );
        costs[div] = cost ;
        if( cost < min_cost ) min_cost = cost;
        if( cost > max_cost ) max_cost = cost;
      }

      for( uint32_t iDiv{0}; iDiv < _scored_divs.size(); ++iDiv )
      {
        uint32_t div = _scored_divs[iDiv].div;
        if( ( div < offset ) && (std::find( supp.begin(), supp.end(), div ) == supp.end()) )
        {
          //costs[div] = exp(-BETA*((costs[div]-min_cost)/(_gSPFD.nEdges-min_cost) ));
          costs[div] = exp(-BETA*(costs[div]/(_gSPFD.nEdges)));
        }
        else
        {
          costs[div] = 0;
        }
      }

      for( uint32_t iDiv{1}; iDiv < divisors.size(); ++iDiv )
      {
        costs[iDiv]+=costs[iDiv-1];
      }

      std::uniform_real_distribution<> distrib(0, 1);
      double rnd = distrib(RNG);

      found=false;
      if( costs.back()>0 )
      {
        for( uint32_t i{1}; i<costs.size(); ++i )
        {
          if( rnd*costs.back() <= costs[i] )
          {
            supp.push_back(i);
            _gSPFD.update( get_div( i ) );
            found = true;
            break;
          }
        }
        
  //      for( auto x : supp)
  //        printf("%d ", x );
  //      printf(" : %f\n", _gSPFD.nEdges );
      }
      if( !found ) break;
    }

    if( _gSPFD.is_covered() )
    {
      std::sort( supp.begin(), supp.end() );
      if( _past_supports.insert(supp).second )
      {
        return supp;
      }
    }
    return std::nullopt;
  }


  std::optional<std::vector<uint32_t>> find_support_stats( uint32_t offset )
  {

    double BETA = 10000;
    offset = offset % _scored_divs.size();
    double cost;
    double mean_cost;
    std::vector<uint32_t> supp;
    std::vector<double> costs;
    double max_cost = std::numeric_limits<double>::max();
    double min_cost = std::numeric_limits<double>::min();
    
    _gSPFD.reset();
    //_scored_divs.print();

    /* add recomputation of the support */
    bool found{false};

    while( !_gSPFD.is_covered() && supp.size() < static_params::max_support_size)
    {
      costs.clear();
      for( uint32_t i{0}; i<divisors.size(); ++i )
      {
        costs.push_back(0);
      }

      if( _gSPFD.is_saturated() ) break;

      for( uint32_t iDiv{1}; iDiv < divisors.size(); ++iDiv )
      {
        cost = _gSPFD.evaluate( get_div( iDiv ) );
        costs[iDiv] = cost ;
        if( cost < min_cost ) min_cost = cost;
        if( cost > max_cost ) max_cost = cost;

      }

      for( uint32_t iDiv{1}; iDiv < _scored_divs.size(); ++iDiv )
      {
        if( std::find( supp.begin(), supp.end(), iDiv ) == supp.end())
        {
          //costs[div] = exp(-BETA*((costs[div]-min_cost)/(_gSPFD.nEdges-min_cost)));
          costs[iDiv] = exp(-BETA*(((costs[iDiv]-min_cost)*(costs[iDiv]-min_cost)/((max_cost-min_cost)*(max_cost-min_cost)))));
        }
        else
        {
          costs[iDiv] = 0;
        }
      }
      costs[0]=0;
      for( uint32_t iDiv{1}; iDiv < divisors.size(); ++iDiv )
      {
        costs[iDiv]+=costs[iDiv-1];
      }

      std::uniform_real_distribution<> distrib(0, 1);
      double rnd = distrib(RNG);

      found=false;
      if( costs.back()>0 )
      {
        for( uint32_t i{1}; i<costs.size(); ++i )
        {
          if( rnd*costs.back() <= costs[i] )
          {
            supp.push_back(i);
            _gSPFD.update( get_div( i ) );
            found = true;
            break;
          }
        }
        
  //      for( auto x : supp)
  //        printf("%d ", x );
  //      printf(" : %f\n", _gSPFD.nEdges );
      }
      if( !found ) break;
    }


    if( _gSPFD.is_covered() )
    {
      std::sort( supp.begin(), supp.end() );
      if( _past_supports.insert(supp).second )
      {
        return supp;
      }
    }
    return std::nullopt;
  }



  scored_divisors_t find_divisors_subset( double center, double threshold )
  {
    double distance;
    scored_divisors_t new_divs;

    for( uint32_t i{0}; i<_scored_divs.size(); ++i )
    {
      distance = abs( _scored_divs[i].cost - center );
      if( distance < threshold )
        new_divs.emplace_back( _scored_divs[i].div, distance );
    }
    new_divs.sort();
    return new_divs;
  }

  std::optional<std::vector<uint32_t>> find_support_greedy_centered( uint32_t iteration, double center, double threshold )
  {
    uint32_t cost, best_cost;
    std::vector<uint32_t> best1, best2;
    std::vector<uint32_t> supp;
    _gSPFD.reset();
    //_scored_divs.print();
    /* find a good starting point */
    scored_divisors_t subset = find_divisors_subset( center, threshold );
    if( subset.divs.size() < 2 ) return std::nullopt;
    //subset.print();
    bool found_one{false};
    /* find a pair to use */
    best_cost = std::numeric_limits<double>::max();
    for( uint32_t i1{0}; i1<subset.divs.size(); ++i1 )
    {
      if( _level_1_divisors.find( subset.divs[i1].div ) != _level_1_divisors.end() ) continue;
      for( uint32_t i2{i1+1}; i2<subset.divs.size(); ++i2 )
      {
        if( _level_1_divisors.find( subset.divs[i2].div ) != _level_1_divisors.end() ) continue;

        double cost = _gSPFD.evaluate( get_div( subset.divs[i1].div ), get_div( subset.divs[i2].div ) );
        if( cost < best_cost  )
        {
          best_cost = cost;
          best1 = {i1};
          best2 = {i2};
          found_one=true;
        }
        else if( cost == best_cost )
        {
          best1.push_back( i1 );
          best2.push_back( i2 );
        }
      }
    }
    //printf("COST(0)=%f\n", _gSPFD.nEdges  );
    if( !found_one )
    {
      return std::nullopt;
    }

    std::uniform_int_distribution<> distrib2(0, best1.size()-1);
    int idx2 = distrib2(RNG);
    supp.push_back( subset.divs[best1[idx2]].div );
    _gSPFD.update( get_div( subset.divs[best1[idx2]].div ) );
    supp.push_back( subset.divs[best2[idx2]].div );
    _gSPFD.update( get_div( subset.divs[best2[idx2]].div ) );
    /* add recomputation of the support */
    //printf("chosen %d %d from %d choices\n", subset.divs[best1[idx2]].div, subset.divs[best2[idx2]].div, subset.divs.size() );

    while( !_gSPFD.is_covered() )
    {
      best1.clear();
      best_cost = std::numeric_limits<double>::max();
      if( _gSPFD.is_saturated() ) break;
      for( uint32_t iCnd{1}; iCnd<divisors.size(); ++iCnd )
      {
        if( _level_1_divisors.find( iCnd ) != _level_1_divisors.end() ) continue;

        cost = _gSPFD.evaluate( get_div( iCnd ) );
        if( cost < best_cost )
        {
          best_cost = cost;
          best1 = {iCnd};
        }
        else if( cost == best_cost && (std::find( supp.begin(), supp.end(), iCnd ) == supp.end()) )
        {
          best1.push_back( iCnd );
        }
      }
      if( best1.size() == 0 ) break;

      std::uniform_int_distribution<> distrib(0, best1.size()-1);
      int idx = distrib(RNG);
      supp.push_back( best1[idx] );
      _gSPFD.update( get_div( best1[idx] ) );
    }
    if( _gSPFD.is_covered() )
    {
      std::sort( supp.begin(), supp.end() );
      if( _past_supports.insert(supp).second )
      {
        return supp;
      }
    }
    return std::nullopt;
  }

  #pragma endregion support_sampling

  #pragma region function_extraction

  std::tuple<truth_table4_t,truth_table4_t> extract_functionality4_from_signatures( std::vector<uint32_t> const& supp )
  {
    truth_table4_t func4;
    truth_table4_t care4;
    truth_table4_t temp4;
    truth_table_t  temp = _gSPFD.care.construct();

    for( uint32_t m{0u}; m < 16u; ++m )
    {
      if( m < ( 1 << supp.size() ) )
      {
        temp = temp | ~temp;
        temp4 = temp4 | ~temp4;

        for( uint32_t l{0u}; l < supp.size(); ++l )
        {
          if( ( m >> l ) & 0x1 == 0x1 )
          {
            temp &= get_div(supp[l]);
            temp4 &= _xs4[l];
          }
          else
          {
            temp &= ~get_div(supp[l]);
            temp4 &= ~_xs4[l];
          }
        }

        if( kitty::count_ones( temp & _gSPFD.care ) > 0 ) // care value
        {
            care4 |= temp4;
            //kitty::set_bit( care4, m );


          if( kitty::count_ones( temp & _gSPFD.func[1] ) > 0 )
          {
            func4 |= temp4;
            //kitty::set_bit( func4, m );
          }
        }
      }
      else
        kitty::clear_bit( care4, m );
    }

    return { func4, care4 };
  }

  std::tuple<truth_tableK_t,truth_tableK_t> extract_functionalityK_from_signatures( std::vector<uint32_t> const& supp )
  {
    truth_tableK_t funcK;
    truth_tableK_t careK;
    truth_tableK_t tempK;
    truth_table_t  temp = _gSPFD.care.construct();

    for( uint32_t m{0u}; m < (1u << static_params::max_support_size); ++m )
    {
      if( m < ( 1 << supp.size() ) )
      {
        temp = temp | ~temp;
        tempK = tempK | ~tempK;

        for( uint32_t l{0u}; l < supp.size(); ++l )
        {
          if( ( m >> l ) & 0x1 == 0x1 )
          {
            temp &= get_div(supp[l]);
            tempK &= _xsK[l];
          }
          else
          {
            temp &= ~get_div(supp[l]);
            tempK &= ~_xsK[l];
          }
        }

        if( kitty::count_ones( temp & _gSPFD.care ) > 0 ) // care value
        {
          careK |= tempK;

          if( kitty::count_ones( temp & _gSPFD.func[1] ) > 0 )
          {
            funcK |= tempK;
          }
        }
      }
      else
        kitty::clear_bit( careK, m );
    }

    return { funcK, careK };
  }

  #pragma endregion function_extraction

  #pragma region boolean_matching_resynthesis


  std::optional<uint32_t> find_boolean_matching( std::array<uint32_t, 4> lits, truth_table4_t const& func4, truth_table4_t const& care4, uint32_t max_num_gates )
  {
    if(VERBOSE)
    {
      printf("TT(0):"); print_tt_with_dcs(func4, care4);
    }

    auto [func_npn, neg, perm] = exact_npn_canonization( func4 );
    if(VERBOSE)
    {
      printf("neg  = ");
      for( auto iBit=3; iBit>=0; iBit-- )
        printf("%d", (neg >> iBit) & 0x1);

      printf(" | perm  =");
      for( auto iBit=0; iBit<4; iBit++ )
        printf("%d ", perm[iBit]);
      printf("\n");
      
      for( auto iBit=0; iBit<4; iBit++ )
      {
        if( neg >> iBit & 0x1 == 0x1 )
          printf("%2d : ~X[%d] <= X[%d]  <<  X[%d] <= P[%d]\n", lits[iBit] ^ 0x1, iBit, iBit, perm[iBit], iBit );
        else
          printf("%2d :  X[%d] <= X[%d]  <<  X[%d] <= P[%d]\n", lits[iBit], iBit, iBit, perm[iBit], iBit );
      }
    }
    
    auto const care_npn = apply_npn_transformation( care4, neg & ~( 1 << 4u ), perm );
    if(VERBOSE)
    {
      printf("npn(TT)"); print_tt_with_dcs(func_npn, care_npn);
    }

    auto const structures = _database.get_supergates( func_npn, ~care_npn, neg, perm );
    if( structures == nullptr )
    {
      return std::nullopt;
    }
    if(VERBOSE)
    {
      printf("neg* = ");
      for( auto iBit=3; iBit>=0; iBit-- )
        printf("%d", (neg >> iBit) & 0x1);

      printf(" | perm* =");
      for( auto iBit=0; iBit<4; iBit++ )
        printf("%d ", perm[iBit]);
      printf("\n");
    }
    bool phase = ( neg >> 4 == 1 ) ? true : false;

    for( auto i{0}; i<lits.size(); ++i )
    {
      if( ( neg >> i ) & 0x1 == 0x1 )
        lits[i] ^= 0x1;
    }

    std::array<uint32_t, 4> leaves {0};

    for( auto i{0}; i<4; ++i )
    {
      leaves[i] = lits[perm[i]];
    }

    auto & db = _database.get_database();

    std::unordered_map<uint64_t, uint32_t> existing_nodes;

    index_list_t index_list_copy = index_list;

    auto res = create_index_list( db.get_node( structures->at(0).root ), leaves, existing_nodes );

    if(VERBOSE)
    {
      printf(" || --> [%d <?= %d]\n", index_list.num_gates(), max_num_gates );
    }
    if( index_list.num_gates() <= max_num_gates )
    {
      return phase != db.is_complemented( structures->at(0).root ) ? (res ^ 0x1 ) : res;//create_index
    }
    //else
    //  index_list = index_list_copy;

    return std::nullopt;
  }

  uint32_t create_index_list( node<aig_network> const& n, std::array<uint32_t, 4u> const& leaves, std::unordered_map<uint64_t, uint32_t> & existing_nodes )
  {
    return create_index_list_rec( n, leaves, existing_nodes );
  }

  uint32_t create_index_list_rec( node<aig_network> const& n, std::array<uint32_t, 4u> const& leaves, std::unordered_map<uint64_t, uint32_t> & existing_nodes )
  {
    auto& db = _database.get_database();
 
    std::array<uint32_t, 2u> node_data;

    int i{0};
    db.foreach_fanin( n, [&]( auto const& f, auto i ) 
    {
      node<aig_network> g = db.get_node( f );
      if( db.is_pi( g ) )
      {
        node_data[i] = db.is_complemented( f ) ? leaves[f.index-1] ^ 0x1 : leaves[f.index-1];
        return;
      }
      else if( db.is_and( g ) )
      {
        auto res = create_index_list_rec( g, leaves, existing_nodes );

        if( res )
        {
          node_data[i] = db.is_complemented( f ) ? (res) ^ 0x1 : res;
        }
      }
    } );

//    if( node_data[0] == 1 && node_data[1] > 1 )
//    {
//      return node_data[1];
//    }
//    else if( node_data[1] == 1 && node_data[0] > 1 )
//    {
//      return node_data[0];
//    }
//    else if( node_data[1] == 0 || node_data[0] == 0 )
//    {
//      return 0;
//    }
    //else 
    if( db.is_and( n ) )
    {
      uint32_t new_lit;
      uint64_t key0 = node_data[0];
      uint64_t key1 = node_data[1];

      if( key0 < key1 )
      {
        key0 = ( key0 << 32u ) | key1;
      }
      else
      {
        key0 = key0 | ( key1 << 32 );
      }

      if( auto search = existing_nodes.find( key0 ); search != existing_nodes.end() )
      {
        new_lit = search->second;
        if(VERBOSE)
        {
          printf("%d=and(%d,%d)* ", new_lit, node_data[0], node_data[1]);
        }
      }
      else
      {
        new_lit = index_list.add_and( node_data[0], node_data[1] );
        if(VERBOSE)
        {
          printf("%d=and(%d,%d) ", new_lit, node_data[0], node_data[1]);
        }
        existing_nodes[key0] = new_lit;
      }
      return new_lit;
    }
    else
    {
      return 0;
    }
  }

  #pragma endregion boolean_matching_resynthesis

  #pragma region spfd_resynthesis
  std::optional<uint32_t> find_spfd_resynthesis( std::vector<uint32_t> const& supp, truth_tableK_t const& funcK, truth_tableK_t const& careK, uint32_t max_num_gates )
  {
    index_list_t index_list_copy = index_list;
    uint32_t max_num_gates_copy = max_num_gates;
    _divsK.set_target( funcK, careK );

    for( auto iTry{0}; iTry<static_params::max_resyn_attempts; ++iTry )
    {
      index_list = index_list_copy;
      max_num_gates = max_num_gates_copy;

      _divsK.set_support( supp, _xsK );
      while( _divsK.size() > 1 && index_list.num_gates() <= max_num_gates )
      {
        if( !_divsK.update( index_list, _functional_library, max_num_gates ) )  break;
      }
      if( _divsK.spfd.is_covered() && _divsK.size() == 1 )
      {
        if( kitty::equal( _divsK.get_div(0) & _divsK.spfd.care, _divsK.spfd.func[1] ) )
        {
          return _divsK[0].lit;
        }
        else if( kitty::equal( _divsK.get_div(0) & _divsK.spfd.care, _divsK.spfd.func[0] ) )
        {
          return _divsK[0].lit ^ 0x1;
        }
        else
        {
          printf( "[w]: one divisor not matching\n" );
        }
      }
    }
    return std::nullopt;
  }
  #pragma endregion spfd_resynthesis

  inline TT const& get_div( uint32_t idx ) const
  {
    if constexpr ( static_params::copy_tts )
    {
      return divisors[idx];
    }
    else
    {
      return ( *ptts )[divisors[idx]];
    }
  }

private:
  std::array<TT, 2> on_off_sets;
  std::array<uint32_t, 2> num_bits; /* number of bits in on-set and off-set */
  uint32_t num_edges;
  
  const typename static_params::truth_table_storage_type* ptts;
  std::vector<std::conditional_t<static_params::copy_tts, TT, typename static_params::node_type>> divisors;
 
  index_list_t index_list;
  //std::vector<scored_lit> lits;

  spfd_manager_t<truth_table_t, 1<<static_params::max_support_size> _gSPFD;
  u_spfd_manager_t<truth_table_t, 1<<static_params::max_support_size> _uSPFD;
  std::array<truth_table4_t,4> _xs4;
  std::array<truth_tableK_t,static_params::max_support_size> _xsK;
  std::set<std::vector<uint32_t>> _past_supports;
  std::unordered_map<uint64_t, uint32_t> _past_supports_count;
  std::vector<uint32_t> _support;
  std::set<uint32_t> _level_1_divisors;
  divisors_t _divsK;
  std::vector<scored_divisor_t> _scored_divs;

  /* positive unate: not overlapping with off-set
     negative unate: not overlapping with on-set */
  std::vector<unate_lit> pos_unate_lits, neg_unate_lits;
  std::vector<uint32_t> binate_divs;
  std::vector<fanin_pair> pos_unate_pairs, neg_unate_pairs;

  functional_library_t _functional_library;

  database_t _database;


  stats& st;
}; /* aig_resyn */ 

} /* namespace aig */

} /* namespace spfd */

} /* namespace mockturtle */


//./experiments/bmatch_resub_aig
//| benchmark | size | gates(SOA) | gates(SPFD) | time(SOA) | time(SPFD) | eq(SOA) | eq(SPFD) |
//|       c17 |    6 |          6 |           6 |      0.00 |       0.16 |    true |     true |
//|      c432 |  208 |        167 |         169 |      0.00 |       0.13 |    true |     true |
//|      c499 |  398 |        392 |         398 |      0.00 |       0.14 |    true |     true |
//|      c880 |  325 |        306 |         320 |      0.00 |       0.14 |    true |     true |
//|     c1355 |  502 |        456 |         502 |      0.00 |       0.16 |    true |     true |
//|     c1908 |  341 |        287 |         302 |      0.00 |       0.14 |    true |     true |
//|     c2670 |  716 |        567 |         634 |      0.01 |       0.16 |    true |     true |
//|     c3540 | 1024 |        841 |         904 |      0.02 |       0.19 |    true |     true |
//|     c5315 | 1776 |       1380 |        1616 |      0.03 |       0.27 |    true |     true |
//|     c6288 | 2337 |       1887 |        2334 |      0.05 |       0.43 |    true |     true |
//|     c7552 | 1469 |       1381 |        1435 |      0.01 |       0.20 |    true |     true |
//
//| benchmark | size | gates(SOA) | gates(SPFD) | time(SOA) | time(SPFD) | eq(SOA) | eq(SPFD) |
//|       c17 |    6 |          6 |           6 |      0.00 |       0.13 |    true |     true |
//|      c432 |  208 |        167 |         169 |      0.00 |       0.13 |    true |     true |
//|      c499 |  398 |        392 |         394 |      0.00 |       0.14 |    true |     true |
//|      c880 |  325 |        306 |         315 |      0.00 |       0.14 |    true |     true |
//|     c1355 |  502 |        456 |         462 |      0.00 |       0.14 |    true |     true |
//|     c1908 |  341 |        287 |         292 |      0.00 |       0.14 |    true |     true |
//|     c2670 |  716 |        567 |         586 |      0.01 |       0.15 |    true |     true |
//|     c3540 | 1024 |        841 |         858 |      0.02 |       0.17 |    true |     true |
//|     c5315 | 1776 |       1380 |        1470 |      0.03 |       0.19 |    true |     true |
//|     c6288 | 2337 |       1887 |        2313 |      0.05 |       0.16 |    true |     true |
//|     c7552 | 1469 |       1381 |        1423 |      0.01 |       0.17 |    true |     true |

