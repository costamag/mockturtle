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
  \file mapped_database.hpp
  \brief Manager for databases of mapped networks

  \author Andrea Costamagna
*/
#pragma once

#include "../../../io/write_verilog.hpp"
#include "../../index_lists/lists/mapped/bound_list.hpp"
#include "../../symm_utils.hpp"

namespace mockturtle
{

/*! \brief Database of mapped networks
 *
 * \tparam NtkDb Network type of the stored database
 */
template<typename NtkDb, uint32_t MaxNumVars = 6u>
class mapped_database
{
  static constexpr bound::design_type_t design_t = NtkDb::design_t;
  using signal_t = typename NtkDb::signal;
  using library_t = bound::augmented_library<design_t>;
  using node_index_t = typename NtkDb::node_index_t;
  using truth_table_t = kitty::static_truth_table<MaxNumVars>;
  using list_simulator_t = list_simulator<bound_list<design_t>, truth_table_t>;

private:
  struct database_entry_t
  {

    bool operator<( database_entry_t const& other ) const
    {
      bool dominates = true;
      dominates &= area < other.area;
      dominates &= switches < other.switches;
      bool one_strict = false;
      for ( int i = 0; i < delays.size(); ++i )
      {
        one_strict |= delays[i] < other.delays[i];
        dominates &= delays[i] <= other.delays[i];
      }
      return dominates & one_strict;
    }

    bool operator>=( database_entry_t const& other ) const
    {
      bool dominated = true;
      dominated &= area >= other.area;
      dominated &= switches >= other.switches;
      for ( int i = 0; i < delays.size(); ++i )
      {
        dominated &= delays[i] >= other.delays[i];
      }
      return dominated;
    }

    /*! \brief Area of the sub-network */
    double area;
    /*! \brief Zero-delay switching */
    uint32_t switches;
    /*! \brief Longest path from each pin to the outputs */
    std::vector<double> delays;
    /*! \brief Node implementing the functionality */
    node_index_t index;
  };

  struct database_row_t
  {
    size_t size() const
    {
      return entries.size();
    }

    void push_back( database_entry_t const& entry )
    {
      entries.push_back( entry );
    }

    database_entry_t& operator[]( uint32_t i )
    {
      return entries[i];
    }

    symmetries_t symm;
    truth_table_t repr;
    std::vector<database_entry_t> entries;
  };

  struct funmap_t
  {
    funmap_t() = default;
    funmap_t( funmap_t&& ) = default;
    funmap_t( funmap_t const& ) = default;

    funmap_t& operator=( funmap_t const& ) = default;
    funmap_t& operator=( funmap_t&& ) = default;

    funmap_t( permutation_t const& perm, uint64_t const& row )
        : perm( perm ), row( row )
    {}

    permutation_t perm;
    uint64_t row;
  };

public:
  mapped_database( library_t& lib )
      : lib_( lib ),
        ntk_( lib ),
        simulator_( lib )
  {
    for ( auto i = 0u; i < MaxNumVars; ++i )
      pis_.push_back( ntk_.create_pi() );

    for ( int i = 0; i < MaxNumVars; ++i )
    {
      kitty::create_nth_var( proj_funcs_[i], i );
    }
    /* extract the list's functionality */
    for ( auto i = 0u; i < MaxNumVars; ++i )
    {
      sims_ptrs_.push_back( &proj_funcs_[i] );
    }
  }

#pragma region Saving

  void commit( std::string const& file )
  {
    write_verilog( ntk_, file );
  }

  void commit( std::ostream& os )
  {
    write_verilog( ntk_, os );
  }

#pragma endregion

#pragma region Getters
  /*! \brief Get the number of rows in the database */
  uint64_t num_rows() const
  {
    return database_.size();
  }

  /*! \brief Get the number of sub-networks stored */
  uint64_t size() const
  {
    return ntk_.num_pos();
  }
#pragma endregion

#pragma region Insert
  uint64_t memoize_func( truth_table_t const& tt )
  {
    uint64_t row;
    const auto it0 = func_to_map_.find( tt );
    if ( it0 != func_to_map_.end() )
    {
      row = ( it0->second ).row;
    }
    else
    {
      auto [repr, _, perm] = kitty::exact_p_canonization( tt );
      const auto it1 = repr_to_row_.find( repr );
      /* insert in the database */
      if ( it1 != repr_to_row_.end() ) // the representative doesn't have any implementation yet
      {
        row = it1->second;
      }
      else
      {
        row = database_.size();
        database_.emplace_back();
        database_.back().symm = symmetries_t( repr );
        database_.back().repr = repr;
        repr_to_row_[repr] = row;
      }

      func_to_map_[tt] = funmap_t{ permutation_t( perm ), row };
    }
    return row;
  }

  /*! \brief Insert a mapped list into the database
   *
   * This method inserts the list into the database the way it is provided.
   * Any check for whether the list should be inserted or not should be handled
   * at a higher level of abstraction.
   */
  bool add( bound_list<design_t> list )
  {
    assert( list.num_pis() == MaxNumVars );
    simulator_( list, sims_ptrs_ );
    truth_table_t const tt = simulator_.get_simulation( list, sims_ptrs_, list.po_at( 0 ) );

    /* perform P-canonization on the list */
    uint64_t row = memoize_func( tt );

    perm_canonize( list, func_to_map_[tt].perm );
    time_canonize( list, lib_, database_[row].symm );

    bool is_inserted = add( list, row );

    simulator_( list, sims_ptrs_ );
    truth_table_t const tt2 = simulator_.get_simulation( list, sims_ptrs_, list.po_at( 0 ) );
    assert( kitty::equal( tt2, database_[row].repr ) );
    /* create the new entry */

    return is_inserted;
  }

private:
  bool add( bound_list<design_t>& list, uint64_t row )
  {
    database_entry_t entry;
    entry.area = list.get_area( lib_ );
    entry.switches = simulator_.get_switches( list );
    entry.delays = get_longest_paths( list, lib_ );

    for ( auto i = 0u; i < database_[row].size(); ++i )
    {
      if ( entry >= database_[row][i] )
      {
        return false;
      }
      else if ( entry < database_[row][i] )
      {
        auto const f = insert( ntk_, pis_, list );
        ntk_.substitute_node( database_[row][i].index, f );
        entry.index = ntk_.get_node( f );
        database_[row][i] = entry;
        return true;
      }
    }
    // TODO: add capacity
    auto const f = insert( ntk_, pis_, list );
    if ( ntk_.is_po( f ) )
      return false; // do not re-insert POs in the database
    ntk_.create_po( f );
    entry.index = ntk_.get_node( f );
    database_[row].push_back( entry );
    return true;
  }
#pragma endregion

private:
  /*! \brief Map a truth table to a storage of nodes and input permutations */
  std::vector<database_row_t> database_;

  /*! \brief Map a completely specified truth table to a database row */
  phmap::flat_hash_map<truth_table_t, funmap_t, kitty::hash<truth_table_t>> func_to_map_;
  phmap::flat_hash_map<truth_table_t, uint64_t, kitty::hash<truth_table_t>> repr_to_row_;

  /*! \brief Database represented as a network */
  NtkDb ntk_;
  std::vector<signal_t> pis_;

  /*! \brief Technology library */
  library_t lib_;

  /*! \brief Simulation engine for mapped lists */
  list_simulator_t simulator_;

  /*! \brief Vector of simulationpatterns */
  std::array<truth_table_t, MaxNumVars> proj_funcs_;
  std::vector<truth_table_t const*> sims_ptrs_;
};

} /* namespace mockturtle */