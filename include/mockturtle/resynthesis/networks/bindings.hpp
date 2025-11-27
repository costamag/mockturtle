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
  \file bind_library.hpp
  \brief Implements methods for handling and evaluating a library of standerd cells.

  This engine can be used for efficient Boolean evaluaton of the gates in a
  standard cell library. Each gate is represented as an AIG index list for
  efficient evaluation. Additionally, in the presence of multiple-output cells,
  this engine identifies which gates belong to a multiple-output and allows handling
  this information.

  \author Andrea Costamagna
*/

#pragma once

#include <map>
#include <vector>

#include "../aig_syn.hpp"
#include "../../io/genlib_reader.hpp"
#include "../../utils/index_list.hpp"
#include "../../utils/standard_cell.hpp"


namespace mockturtle
{

struct gate_with_list // : gate
{
  gate_with_list( gate genlib_gate, large_xag_index_list const& index_list, int cell_group )
    : list( index_list ), group( cell_group ), pin_delays( genlib_gate.num_vars ), pins( genlib_gate.pins )
  {
    name = genlib_gate.name;
    id = genlib_gate.id;
    name = genlib_gate.name;
    expression = genlib_gate.expression;
    num_vars = genlib_gate.num_vars;
    function = genlib_gate.function;
    area = genlib_gate.area;
    output_name = genlib_gate.output_name;

    /* compute the pin delays */
    for( int i = 0; i < num_vars; ++i )
    {
      double const rise_time = genlib_gate.pins[i].rise_block_delay;
      double const fall_time = genlib_gate.pins[i].fall_block_delay;
      double const pin_delay = std::max( rise_time, fall_time );
      pin_delays[i] = pin_delay;
    }

  }

  unsigned int id;
  std::string name;
  std::string expression;
  uint32_t num_vars;
  kitty::dynamic_truth_table function;
  double area;
  std::vector<pin> pins;
  std::string output_name;
  large_xag_index_list list;
  int group;
  std::vector<double> pin_delays;

};

struct multiple_output_group
{
  multiple_output_group( unsigned int id, std::vector<unsigned int> ids )
    :  id( id ), ids( ids )
  {
  }

  /* unique identifier of the multiple-output group */
  unsigned int id;
  /* genlib identifiers of the gates in the group */
  std::vector<unsigned int> ids;
};

class bindings
{
public:
  /* truth table used to express the gate functionality */
  using truth_table_t = kitty::dynamic_truth_table;
  using index_list_t = large_xag_index_list;

  bindings()
    : ggates(),
      cells(),
      name_to_cell(),
      groups(),
      gates()
  {}


  bindings( std::vector<gate> const& genlib )
    : ggates( genlib ),
      cells( get_standard_cells( genlib ) ),
      name_to_cell(),
      groups( identify_groups( genlib, cells ) ),
      gates( construct_gates( genlib, groups ) )
  {
  }

  bindings( bindings const& other )
    : bindings( other.get_genlib_gates() )
  {
  }

  size_t size()
  {
    return gates.size();
  }

  void initialize( std::vector<gate> const& genlib )
  {
    if( gates.size() == genlib.size() )
      return;
    cells = get_standard_cells( genlib ),
    groups = identify_groups( genlib, cells );
    gates = construct_gates( genlib, groups );
  }

  gate_with_list const operator[]( unsigned int id ) const
  {
    return gates[id];
  }

  std::vector<gate_with_list> get_gates() const
  {
    return gates;
  }

  std::vector<gate> get_genlib_gates() const
  {
    return ggates;
  }

  index_list_t get_list( unsigned int id ) const
  {
    return gates[id].list;
  }

  std::string get_name( unsigned int id ) const
  {
    return gates[id].name;
  }

  double get_area( unsigned int id ) const
  {
    return gates[id].area;
  }

  std::string get_expression( unsigned int id ) const
  {
    return gates[id].expression;
  }

  int get_group( unsigned int id ) const
  {
    return gates[id].group;
  }

  int get_pin_delay( unsigned int id, int i ) const
  {
    double const fall = gates[id].pins[i].fall_block_delay;
    double const rise = gates[id].pins[i].rise_block_delay;
    return std::max( rise, fall );
  }

  uint32_t get_group_id( std::string const& name )
  {
    uint32_t cell_id = name_to_cell[name];
    uint32_t gate_id = cells[cell_id].gates[0].id;
    return gates[gate_id].group;
  }

  std::vector<unsigned int> get_ids( std::string const& name )
  {
    uint32_t cell_id = name_to_cell[name];
    auto const& gates = cells[cell_id].gates;
    std::vector<unsigned int> ids( gates.size() );
    std::transform( gates.begin(), gates.end(), ids.begin(), [&]( auto g ){ return g.id; } );
    return ids;
  }

  truth_table_t get_function( unsigned int id ) const
  {
    return gates[id].function;
  }

  truth_table_t get_function( std::string const& name )
  {
    auto id = get_ids( name )[0];
    return gates[id].function;
  }

  gate_with_list get_binding( unsigned int id ) const
  {
    return gates[id];
  }

  bool is_multiple_output( unsigned int id ) const
  {
    return get_group( id ) >= 0;
  }

private:
  /*! \brief identify multiple-output groups
   * 
   * This method takes as input the gates in a genlib library and
   * groups identifiers belonging to the same multiple-output cell
   * into an unique group.
   * 
  */
  std::vector<multiple_output_group> identify_groups( std::vector<gate> const& genlib, std::vector<standard_cell> const& cells )
  {
    /* extract the standard cell-description */
    std::vector<int> gates_groups( genlib.size(), -1 );
    std::vector<multiple_output_group> groups;

    for ( auto cell : cells )
    {
      name_to_cell[cell.name] = cell.id;
      if ( cell.gates.size() > 1 )
      {
        std::vector<unsigned int> ids;
        for( auto g : cell.gates )
        {
          gates_groups[g.id] = groups.size();
          ids.push_back( g.id );
        }
        groups.emplace_back( groups.size(), ids );
      }
    }

    return groups;
  }

  /*! \brief synthesize the gates as index lists
   * 
   * This method synthesizes the index lists for the gates in a library
   * and assigns to them the indices of the groups to which they belong 
  */
  std::vector<gate_with_list> construct_gates( std::vector<gate> const& genlib, std::vector<multiple_output_group> const& groups )
  {
    /* mark each cell as single-output or as part of a multiple-output group */
    std::vector<int> cell_groups( genlib.size(), -1 );

    for( auto const& group : groups )
    {
      for( auto const& id : group.ids )
      {
        cell_groups[id] = group.id;
      }
    }
    /* synthesize each gate and label it based on the multiple-output group */
    std::vector<gate_with_list> gates;
    gates.reserve( genlib.size() );
    mockturtle::aig_syn_decompose<truth_table_t> syn_engine;
    for( auto g : genlib )
    {
//std::cout << g.name << std::endl;
      index_list_t const list = syn_engine( g.function );
//std::cout << to_index_list_string(list) << std::endl;
      gates.emplace_back( g, list, cell_groups[g.id] );
    }

    return gates;
  }

public:
  double get_max_area()
  {
    double res = 0;
    for ( auto g : ggates )
    {
      res = std::max( res, g.area );
    }
    return res;
  }

private:

  std::vector<gate> ggates;
  /* standard cells */
  std::vector<standard_cell> cells;
  /* maps a gate name to the standad cell index */
  std::map<std::string, uint32_t> name_to_cell;
  /* multiple-output groups */
  std::vector<multiple_output_group> groups;
  /* collection of standard cells */
  std::vector<gate_with_list> gates;
};

} // namespace mockturtle