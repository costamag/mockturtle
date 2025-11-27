/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2023  EPFL
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
  \file mapped_list.hpp
  \brief Index list for mapped networks

  \author Andrea Costamagna
*/

#pragma once

#include <vector>

#include <kitty/hash.hpp>
#include <kitty/constructors.hpp>

#include <parallel_hashmap/phmap.h>

#include "mapped_list.hpp"
#include "../list_simulator.hpp"
#include "../mapped_list_simulator.hpp"
#include "../networks/bindings.hpp"

namespace mockturtle
{
  enum choice_order : uint32_t
  {
    area_delay_power = 0,
    area_power_delay = 1,
    power_delay_area = 2,
    power_area_delay = 3,
    delay_area_power = 4,
    delay_power_area = 5
  };

  template<uint32_t num_vars, uint32_t capacity, choice_order order>
  class homo_entry
  {
    using TT = kitty::static_truth_table<num_vars>;
    static constexpr double max_double = std::numeric_limits<double>::max();

  public:
    homo_entry( TT const& tt ):
    tt( tt ), symm( num_vars )
    {
      compute_symmetries();
    }

    homo_entry(){}

    ~homo_entry() {}

    homo_entry( const homo_entry& other )
      : tt( other.tt ), data( other.data ), symm( other.symm ), var_to_symm( other.var_to_symm )
        {}

    homo_entry( homo_entry&& other ) noexcept
      : tt( other.tt ), data( other.data ), symm( other.symm ), var_to_symm( other.var_to_symm )
        {}

    homo_entry& operator=( const homo_entry& other ) {
      if (this != &other) {
        data = other.data;
        tt = other.tt;
        symm = other.symm;
        var_to_symm = other.var_to_symm;
      }
      return *this;
    }

    homo_entry& operator=( homo_entry&& other ) noexcept {
        if (this != &other) {
          data = other.data;
          tt = other.tt;
          symm = other.symm;
          var_to_symm = other.var_to_symm;
        }
        return *this;
    }

    void compute_symmetries()
    {
      symm.clear();
      for( uint32_t i{0}; i < num_vars - 1; ++i )
      {
        var_to_symm[i] = -1;
      }

      std::vector<uint32_t> symmetry;
      uint32_t done {0u};
      for( uint32_t i{0}; i < num_vars - 1; ++i )
      {
        if( ( ( done >> i ) & 0x1 ) == 0x1 )
          continue;
        if ( !kitty::has_var( tt, i ) )
        {
          continue;
        }
        symmetry = { i };
        TT const F1 = kitty::cofactor1( tt, i );
        TT const F0 = kitty::cofactor0( tt, i );

        for( uint32_t j{ i + 1 }; j < num_vars; ++j )
        {
          if ( !kitty::has_var( tt, j ) )
          {
            continue;
          }
          TT const tt01 = kitty::cofactor0( F1, j );
          TT const tt10 = kitty::cofactor1( F0, j );
          if( kitty::equal( tt01, tt10 ) )
          {
            symmetry.push_back( j );
            done |= ( 1u << j );
          }
        }
        if( symmetry.size() > 1 )
        {
          int const index = symm.size();

          for( auto v : symmetry )
          {
            var_to_symm[v] = index;
          }

          symm.push_back( symmetry );

        }
      }
    }

    /*! \brief insert a mapped list in a database entry after Time-canonization
     * It assumes that the list has the same functionality of the entry where
     * we are adding it, and that the symmetries in the function have already been 
     * extracted in the vector symm.
     * This method performs T-canonization. The variables in the same
     * symmetry group are sorted in decreasing order of propagation delay.
     */
    void insert( mapped_list const& list )
    {
      //std::cout << "try inserting" << std::endl;
      //list.print();
      /* perform t-canonization of the list */
      mapped_list mlist( list );
      t_canonize( mlist );
      /* extract the metrics */
      double const & area = mlist.get_area();
      double const& power = mlist.get_power();
      auto const & pin_delays = mlist.get_pin_delays();

      double worst_power = max_double;
      double worst_area = max_double;
      double worst_delay = max_double;
      int worst = -1;
      double average_delay;
      /* check for dominance */
      for( int i = 0; i < data.size(); ++i )
      {
        double const & area_curr = data[i].get_area();
        double const& power_curr = data[i].get_power();
        auto const & pin_delays_curr = data[i].get_pin_delays();

        bool same_pin_delays = true;
        for( int p = 0; p < pin_delays.size(); ++p )
        {
          same_pin_delays &= std::abs( pin_delays[p] - pin_delays_curr[p] ) <= 0.01;
        }

        bool const is_same = ( std::abs(area - area_curr) <= 0.01 ) && ( std::abs(power_curr - power) <= 0.01 ) && same_pin_delays;
        if( is_same )
        {
          return;
        }

        int win_curr = 0;
        int win_cand = 0;
        average_delay = 0;
        for( int j = 0; j < num_vars; ++j )
        {
          if( pin_delays_curr[j] > pin_delays[j] )
          {
            win_cand += 1;
          }
          else if( pin_delays_curr[j] < pin_delays[j] )
          {
            win_curr += 1;
          }
          average_delay += ( pin_delays_curr[j] - pin_delays[j] )/num_vars;
        }
        bool area_dominance = area_curr > area;
        bool power_dominance = power_curr > power;
        bool delay_dominance = win_curr == 0 && win_cand > 0;
        bool area_dominated = area_curr < area;
        bool power_dominated = power_curr < power;
        bool delay_dominated = win_curr > 0 && win_cand == 0;
        if( area_dominance && power_dominance && delay_dominance )
        {

          data[i] = mlist;
          return;
        }
        else if( area_dominated && power_dominated && delay_dominated )
        {
          return;
        }
        /* none of the above: accumulate information for replacement */
        if( data.size() < capacity )
        {
          bool cond1 {false};
          bool cond2 {false};
          bool cond3 {false};
          /* which one sucks the most? */
          if( order == area_delay_power )
          {
            cond1 = area_curr > worst_area;
            cond2 = ( area_curr == worst_area ) && ( average_delay > worst_delay );
            cond3 = ( area_curr == worst_area ) && ( average_delay == worst_delay ) && ( power_curr > worst_power );
          }
          else if( order == area_power_delay )
          {
            cond1 = area_curr > worst_area;
            cond2 = ( area_curr == worst_area ) && ( power_curr > worst_power );
            cond3 = ( area_curr == worst_area ) && ( power_curr == worst_power ) && ( average_delay > worst_delay );
          }
          else if( order == delay_area_power )
          {
            cond1 = ( average_delay > worst_delay );
            cond2 = ( average_delay == worst_delay ) && ( area_curr > worst_area );
            cond3 = ( average_delay == worst_delay ) && ( area_curr == worst_area ) && ( power_curr > worst_power );
          }
          else if( order == delay_power_area )
          {
            cond1 = ( average_delay > worst_delay );
            cond2 = ( average_delay == worst_delay ) && ( power_curr > worst_power );
            cond3 = ( average_delay == worst_delay ) && ( power_curr == worst_power ) && ( area_curr > worst_area );
          }
          else if( order == power_delay_area )
          {
            cond1 = ( power_curr > worst_power );
            cond2 = ( power_curr == worst_power ) && ( average_delay > worst_delay );
            cond3 = ( power_curr == worst_power ) && ( average_delay == worst_delay ) && ( area_curr > worst_area );
          }
          else if( order == power_area_delay )
          {
            cond1 = ( power_curr > worst_power );
            cond2 = ( power_curr == worst_power ) && ( area_curr > worst_area );
            cond3 = ( power_curr == worst_power ) && ( area_curr == worst_area ) && ( average_delay > worst_delay );
          }
          if( cond1 || cond2 || cond3 )
          {
            worst_area = area_curr;
            worst_delay = average_delay;
            worst_power = power_curr;
            worst = i;
          }
        }
      }
      if( data.size() < capacity  )
      {
        data.push_back( mlist );
      }
      else
      {
        /* use the heuristic to find the best replacement */
        if( worst >= 0 )
        {
          data[worst] = mlist;
        }
      }
    }

  struct timed_index
  {
    uint32_t index;
    double delay;
  };

  /*! \brief perform t-canonization of a mapped list for a database entry
   * The variables in the same symmetry group are sorted from the slowest to the fastest
   */
  void t_canonize( mapped_list & mlist )
  {
    /* extract the propagation time */
    auto const pin_delays = mlist.get_pin_delays();
    std::array<timed_index, num_vars> inputs;
    for( int i = 0; i < num_vars; ++i )
    {
      inputs[i].index = i;
      inputs[i].delay = pin_delays[i];
    }
    /* consider all the variables */
    bool swapped;
    for( auto const& vars : symm )
    {
      int const n = vars.size();
      /* bubble-sort variables in a symmetry group */
      for( int v1 = 0; v1 < n - 1; ++v1 )
      {
        swapped = false;
        for( int v2 = n - 1; v2 >= v1 + 1; --v2 )
        {
          uint32_t const var1 = vars[v2-1];
          uint32_t const var2 = vars[v2];
          if( inputs[var1].delay < inputs[var2].delay )
          {
            std::swap( inputs[var1], inputs[var2] );
            swapped = true;
          }
        }
      }
    }
    /* now all the variables in a symmetry group are sorted from slowest to fastest */
    std::vector<uint8_t> perm( num_vars );
    for( int i=0; i<num_vars; ++i )
      perm[i] = inputs[i].index;
    mlist.permute_inputs( perm );
  }

  void t_canonization( mapped_list const& mlist, TT const& tt )
  {
    std::array<uint32_t, num_vars> symm;
    for( uint32_t i{0}; i < num_vars; ++i )
    {
      symm[i] = ( 1u << i );
      for( uint32_t j{ 0 }; j < num_vars; ++j )
      {
        if( i == j )
          continue;
        TT const tt01 = kitty::cofactor0( kitty::cofactor1( tt, i ), j );
        TT const tt10 = kitty::cofactor1( kitty::cofactor0( tt, i ), j );
        if( kitty::equal( tt01, tt10 ) )
        {
          symm[i] |= ( 1u << j );
          symm[j] |= ( 1u << i );
        }
      }
    }
    /* now I have all the symmetries */
    auto pin_delays = mlist.get_pin_delays();
    std::vector<uint32_t> res;
    for( uint32_t p = 0; p < num_vars; ++p )
    {
      res.push_back( p );
    }
    for( uint32_t i{0}; i < num_vars; ++i )
    {
      std::vector<std::pair<uint32_t, double>> unsort_group;
      uint32_t simmetry = symm[i];
      for( uint32_t j{ 0 }; j < num_vars; ++j )
      {
        if( ( ( simmetry >> j ) & 0x1 ) == 0x1 )
          unsort_group.push_back( std::make_pair( j, pin_delays[j] ) );
      }
      std::vector<std::pair<uint32_t, double>> sorted_group = unsort_group;
      std::stable_sort( sorted_group.begin(), sorted_group.end(), [&]( auto a, auto b ){
        return a.second > b.second;
      } );

      for( int k = 0; k < sorted_group.size(); ++k )
      {
        res[unsort_group[k].first] = sorted_group[k].first;
      }

      for( uint32_t j{ i }; j < num_vars; ++j )
      {
        symm[j] &= ~simmetry;
      }
    }
    return std::make_tuple( res, symm );
  }


    std::vector<std::vector<uint32_t>> get_symm() const
    {
      return symm;
    }

    size_t size() const
    {
      return data.size();
    }

    /*! \brief iterator to apply a lambda to all the nodes */
    template<typename Fn>
    void foreach_list( Fn&& fn ) const
    {
      for ( int i = 0; i < size(); i++ )
      {
        fn( data[i], i );
      }
    }

  mapped_list const& operator[]( uint32_t index ) const
  {
    mapped_list const & list = data[index];
    return list;
  }

  private:
    std::vector<mapped_list> data;
    TT const tt;
    std::vector<std::vector<uint32_t>> symm;
    std::array<int, num_vars> var_to_symm;
  };


  template<uint32_t num_vars, uint32_t capacity, choice_order order>
  class homo_storage
  {
  public:
    using TT = kitty::static_truth_table<num_vars>;
    using tt_hash = kitty::hash<TT>;
    using entry_t = homo_entry<num_vars, capacity, order>;

  /*! \brief add an entry to the database
   * assumes that the truth table and the list have the same functionality.
   * the function may be P-canonized or not.
   */
  void add_entry( TT const& tt, mapped_list const& list )
  {
    auto match = find( tt );
    if ( match == end() )
    {
      uint32_t index = data.size();
      tt_to_index[tt] = index;
      data.emplace_back( tt );
      data[index].insert( list );
    }
    else // TODO: add the reward-based check
    {
      uint32_t index = match->second;
      data[index].insert( list );
    }
  }

  uint32_t get_index( TT const& tt )
  {
    uint32_t index = tt_to_index[tt];
    return index;
  }

  entry_t const& at( uint32_t index ) const
  {
    entry_t const& entry = data[index];
    return entry;
  }

  auto const get_symm( uint32_t index ) const
  {
    entry_t const& entry = data[index];
    return entry.get_symm();
  }

  entry_t const& operator[]( uint32_t index ) const
  {
    entry_t const & entry = data[index];
    return entry;
  }

  entry_t & operator[]( TT const& tt )
  {
    uint32_t const index = get_index( tt );
    return this->at(index);
  }

  auto find( TT const& tt )
  {
    return tt_to_index.find( tt );
  }

  auto end()
  {
    return tt_to_index.end();
  }

  /*! \brief iterator to apply a lambda to all the nodes */
  template<typename Fn>
  void foreach_entry( Fn&& fn ) const
  {
    for ( int i = 0; i < data.size(); i++ )
    {
      fn( data[i], i );
    }
  }

    /*! \brief iterator to apply a lambda to all the nodes */
    template<typename Fn>
    void foreach_list( uint32_t index, Fn&& fn ) const
    {
      data[index].foreach_list( fn );
    }

  private:
    std::vector<entry_t> data;
    phmap::flat_hash_map<TT, uint32_t, tt_hash> tt_to_index;
  };

} // namespace mockturtle