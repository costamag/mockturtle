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

#include "../../../../io/genlib_reader.hpp"
#include "../../../../io/write_verilog.hpp"
#include "../../../../networks/block.hpp"
#include "../../../../networks/mapped/bound_network.hpp"
#include "../../../../symm_utils.hpp"
#include "../../../../utils/index_lists/lists/mapped/bound_list.hpp"
#include "../../../../views/cell_view.hpp"
#include "../../../emap.hpp"
#include "../../../klut_to_graph.hpp"
#include "../../../rewrite.hpp"

namespace mockturtle
{

namespace bound
{

/*! \brief Engine to initialize the database with simple structures
 *
 * \tparam DesignType CELL_BASED or ARRAY_BASED.
 * \tparam MaxNumVars Maximum number of input variables. Works best for 6 or lower
 * \tparam MaxNumOuts Maximum number of outputs in a cell. For now only 1 output is supported
 *
 * \verbatim embed:rst

   Example

   .. code-block:: c++

      database_generator gen( gates );
      gen.area_oriented_generation( "asap7_database" );
   \endverbatim
 */
template<design_type_t DesignType = design_type_t::CELL_BASED, uint32_t MaxNumVars = 6, uint32_t MaxNumOuts = 2>
class database_generator
{
  using library_t = augmented_library<DesignType>;
  using network_t = bound_network<DesignType, MaxNumOuts>;
  using database_t = mapped_database<network_t, MaxNumVars>;
  using Ntk_t = cell_view<block_network>;
  using Signal_t = typename cell_view<block_network>::signal;
  using Tt_t = kitty::dynamic_truth_table;
  using tt_set = std::unordered_set<Tt_t, kitty::hash<Tt_t>>;

public:
  database_generator( std::vector<gate> const& gates )
      : library_( gates ), gates_( gates ), db_( library_ )
  {
    init();
  }

  void area_oriented_generation( std::string const& output_file )
  {
    // Preprocess the AIG
    aig_preprocessing();

    // Perform technology mapping
    auto mapped = map_to_block_network();

    // Create a database from the mapped network
    create_database_from_mapped( mapped, output_file + ".v" );
  }

private:
  void init()
  {
    // Store a representative for each P class
    auto const classes = load_p_representatives();

    // Collect the classes into a single kLUT network
    auto const klut = classes_to_klut( classes );

    // Transform the kLUT network into an AIG
    aig_ = convert_klut_to_graph<aig_network, klut_network>( klut );
  }

  /*! \brief Load the P representatives in a truth-table set
   */
  tt_set load_p_representatives()
  {
    tt_set classes;
    classes.reserve( 222 );

    Tt_t tt( 4 );
    do
    {
      auto const res = kitty::exact_p_canonization( tt );
      classes.emplace( std::move( std::get<0>( res ) ) );
      kitty::next_inplace( tt );
    } while ( !kitty::is_const0( tt ) );

    return classes;
  }

  /*! \brief Convert a truth-table set into a kLUT network where each function is a PO
   */
  klut_network classes_to_klut( tt_set const& classes )
  {
    klut_network klut;

    std::vector<typename klut_network::signal> pis( 4 );
    for ( auto& pi : pis )
      pi = klut.create_pi();

    for ( auto const& entry : classes )
    {
      auto const f = klut.create_node( pis, entry );
      klut.create_po( f );
    }

    return klut;
  }

  /*! \brief Area-oriented AIG minimization
   */
  void aig_preprocessing()
  {
    const xag_npn_resynthesis<aig_network> resyn;
    exact_library_params eps;
    eps.np_classification = false;
    const exact_library<aig_network> exact_lib{ resyn, eps };
    rewrite_params ps;
    ps.preserve_depth = true;

    int iteration = 0;
    while ( iteration++ < 10 )
    {
      const auto size_before = aig_.num_gates();
      rewrite( aig_, exact_lib, ps );
      const auto size_after = aig_.num_gates();
      if ( size_before <= size_after )
        break;
    }
  }

  /*! \brief Area-oriented technology-mapping
   */
  Ntk_t map_to_block_network()
  {
    tech_library_params tps;
    tps.ignore_symmetries = false;
    tps.verbose = true;

    tech_library<9> tech_lib( gates_, tps );

    emap_params mps;
    mps.matching_mode = emap_params::hybrid;
    mps.area_oriented_mapping = true;
    mps.map_multioutput = false;
    mps.relax_required = 0;

    emap_stats mst;

    Ntk_t const ntk = emap<9>( aig_, tech_lib, mps, &mst );
    return ntk;
  }

  void create_database_from_mapped( Ntk_t& ntk, std::string const& output_file )
  {
    std::vector<Signal_t> pis( ntk.num_pis() );
    ntk.foreach_pi( [&]( Signal_t const& f, uint32_t const& i ) {
      pis[i] = f;
    } );
    bound::augmented_library<bound::design_type_t::CELL_BASED> lib( gates_ );
    int i = 0;
    ntk.foreach_po( [&]( auto f ) {
      db_.add( ntk, pis, f );
    } );

    db_.commit( output_file );
  }

private:
  std::vector<gate> gates_;
  library_t library_;
  aig_network aig_;
  database_t db_;
};

} /* namespace bound */

} /* namespace mockturtle */