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

#if 0
  template<uint32_t num_vars, uint32_t capacity>
  class easy_entry
  {
    using TT = kitty::static_truth_table<num_vars>;
    static constexpr double max_double = std::numeric_limits<double>::max();

  public:
    easy_entry( TT const& tt ):
    tt( tt ), symm( num_vars )
    {
      compute_symmetries();
    }

    easy_entry(){}

    ~easy_entry() {}

    easy_entry( const easy_entry& other )
      : tt( other.tt ), data( other.data ), symm( other.symm ), var_to_symm( other.var_to_symm )
        {}

    easy_entry( easy_entry&& other ) noexcept
      : tt( other.tt ), data( other.data ), symm( other.symm ), var_to_symm( other.var_to_symm )
        {}

    easy_entry& operator=( const easy_entry& other ) {
      if (this != &other) {
        data = other.data;
        tt = other.tt;
      }
      return *this;
    }

    easy_entry& operator=( easy_entry&& other ) noexcept {
        if (this != &other) {
          data = other.data;
          tt = other.tt;
        }
        return *this;
    }

    easy_entry& operator=( easy_entry&& other ) noexcept {
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
      /* perform t-canonization of the list */
      mapped_list mlist( list );

      if( data.size() < capacity  )
      {
        data.push_back( mlist );
      }
    }

  struct timed_index
  {
    uint32_t index;
    double delay;
  };

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
#endif


  template<uint32_t num_vars, uint32_t capacity>
  class easy_entry
  {
    using TT = kitty::static_truth_table<num_vars>;
    static constexpr double max_double = std::numeric_limits<double>::max();

  public:
    easy_entry( TT const& tt ):
    tt( tt ), symm( num_vars )
    {
      compute_symmetries(); //TODO remove this comment. It must be computed
    }

    easy_entry(){}

    ~easy_entry() {}

    easy_entry( const easy_entry& other )
      : tt( other.tt ), data( other.data ), symm( other.symm ), var_to_symm( other.var_to_symm )
        {}

    easy_entry( easy_entry&& other ) noexcept
      : tt( other.tt ), data( other.data ), symm( other.symm ), var_to_symm( other.var_to_symm )
        {}

    easy_entry& operator=( const easy_entry& other ) {
      if (this != &other) {
        data = other.data;
        tt = other.tt;
        symm = other.symm;
        var_to_symm = other.var_to_symm;
      }
      return *this;
    }

    easy_entry& operator=( easy_entry&& other ) noexcept {
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
      /* perform t-canonization of the list */
      mapped_list mlist( list );
      /* extract the metrics */
      for ( auto elist : data )
      {
        if ( elist == mlist )
          return;
      }
      if( data.size() < capacity  )
      {
        data.push_back( mlist );
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



  template<uint32_t num_vars, uint32_t capacity>
  class easy_storage
  {
  public:
    using TT = kitty::static_truth_table<num_vars>;
    using tt_hash = kitty::hash<TT>;
    using entry_t = easy_entry<num_vars, capacity>;

  /*! \brief add an entry to the database
   * assumes that the truth table and the list have the same functionality.
   * the function may be P-canonized or not.
   */
  void add_entry( TT const& tt, mapped_list const& list )
  {
    //kitty::print_binary( tt );
    //std::cout << " <- adding " << std::endl;
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