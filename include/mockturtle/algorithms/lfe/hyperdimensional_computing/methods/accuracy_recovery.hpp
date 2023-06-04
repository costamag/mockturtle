/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2021  EPFL
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
  \file accuracy_recovery.hpp
  \brief Contains the accuracy recovery methods.
  \author Andrea Costamagna
*/

#pragma once

#include "../../chatterjee_method.hpp"
#include "../../sim_create_nodes.hpp"
#include "../../sim_decomposition_fast.hpp"
#include "../../sim_decomposition_fastS.hpp"
#include "../../sim_decomposition_xor.hpp"
#include "../../forest_decomposition.hpp"
#include "../../forest_decompositionx2.hpp"
#include "../../dc_decomposition_fastS.hpp"
#include "selectors.hpp"
#include <kitty/print.hpp>
#include <kitty/properties.hpp>

namespace mockturtle
{
namespace hdc
{
namespace detail
{
enum class arecovery_method
{
  none,
  SD,
  DK_X,
  DK_SD,
  DK_TSD,
  DK_XTSD,
  DK_XTSDS,
  DK_TSDS,
  DC_TSD,
  DC_XTSD,
  DC_IXTSD,
  DK_DSD,
  DK_RDSD,
  forestS,
  forestSx2,
  xforestS,
  xforestSx2,
  DK_DSDS
};

class arecovery_params
{
  public:
    uint32_t output{0};
    bool verbose{true};
    uint32_t max_sup{2};

    uint32_t num_trees{3};
    uint32_t nImpurity{0};
};


/*! \brief 
 * Statistics based decomposition: tries top decomposition and performs shannon decomposition in case of failure
 */
template<class Ntk>
signal<Ntk> best_node( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
{

  std::vector<kitty::partial_truth_table*> X;
  X.push_back( &ntk.sim_patterns[0].pat );

  double Imax = -std::numeric_limits<double>::max();
  double Inew ;
  signal<Ntk> bsig;
  
  for( auto sim : ntk.sim_patterns )
  {
    X[0] = &sim.pat ;
    Inew = kitty::mutual_information(X, &ntk.targets[ps.output]);
    if( Inew > Imax )
    {
      bsig = sim.sig;
      Imax = Inew;
    }
  }
  double accuracy = 100*(double)kitty::count_ones( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[bsig]].pat) )/ntk.targets[ps.output].num_bits();

  if( ps.verbose )
  {
    std::cout << "[o " << ps.output << "] : " << accuracy << "%" <<std::endl;
    kitty::print_binary( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[bsig]].pat) );
  }
  return bsig;
}

/*! \brief 
 * Statistics based decomposition: tries top decomposition and performs shannon decomposition in case of failure
 */
template<class Ntk>
signal<Ntk> SD( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
{
  sim_decomposition_fastS_params decps;
  decps.verbose = ps.verbose;
  decps.max_sup = ps.max_sup;
  decps.is_informed = false;
  decps.is_size_aware = false;
  decps.try_top_decomposition = false;
  decps.try_bottom_decomposition = false;
  decps.use_correlation = false;

  std::vector<kitty::partial_truth_table> examples;
  for( auto sim : ntk.sim_patterns )
    examples.push_back( sim.pat );

  signal<Ntk> osignal = sim_decomposition_fastS( ntk, examples, ntk.targets[ps.output], decps, false );
  double accuracy = 100*(double)kitty::count_ones( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) )/ntk.targets[ps.output].num_bits();

  if( ps.verbose )
  {
    std::cout << "[o " << ps.output << "] : " << accuracy << "%" <<std::endl;
    kitty::print_binary( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) );
  }

  return osignal;
}

/*! \brief 
 * Statistics based decomposition: tries top decomposition and performs shannon decomposition in case of failure
 */
template<class Ntk>
signal<Ntk> DK_SD( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
{
  sim_decomposition_fastS_params decps;
  decps.verbose = ps.verbose;
  decps.max_sup = ps.max_sup;
  decps.is_informed = true;
  decps.is_size_aware = false;
  decps.try_top_decomposition = false;
  decps.try_bottom_decomposition = false;
  decps.use_correlation = false;

  std::vector<kitty::partial_truth_table> examples;
  for( auto sim : ntk.sim_patterns )
    examples.push_back( sim.pat );

  signal<Ntk> osignal = sim_decomposition_fastS( ntk, examples, ntk.targets[ps.output], decps, false );
  double accuracy = 100*(double)kitty::count_ones( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) )/ntk.targets[ps.output].num_bits();

  if( ps.verbose )
  {
    std::cout << "[o " << ps.output << "] : " << accuracy << "%" <<std::endl;
    kitty::print_binary( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) );
  }

  return osignal;
}

/*! \brief 
 * Statistics based decomposition: tries top decomposition and performs shannon decomposition in case of failure
 */
template<class Ntk>
signal<Ntk> DK_TSD( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
{
  sim_decomposition_fastS_params decps;
  decps.verbose = ps.verbose;
  decps.max_sup = ps.max_sup;
  decps.is_informed = true;
  decps.is_size_aware = false;
  decps.try_top_decomposition = true;
  decps.try_bottom_decomposition = false;
  decps.use_correlation = false;
  decps.try_xor = false;

  std::vector<kitty::partial_truth_table> examples;
  for( auto sim : ntk.sim_patterns )
    examples.push_back( sim.pat );

  signal<Ntk> osignal = sim_decomposition_fastS( ntk, examples, ntk.targets[ps.output], decps, false );
  double accuracy = 100*(double)kitty::count_ones( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) )/ntk.targets[ps.output].num_bits();

  if( ps.verbose )
  {
    std::cout << "[o " << ps.output << "] : " << accuracy << "%" <<std::endl;
    kitty::print_binary( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) );
  }

  return osignal;
}

/*! \brief 
 * Statistics based decomposition: tries top decomposition and performs shannon decomposition in case of failure
 */
template<class Ntk>
signal<Ntk> DK_XTSD( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
{
  sim_decomposition_fastS_params decps;
  decps.verbose = ps.verbose;
  decps.max_sup = ps.max_sup;
  decps.is_informed = true;
  decps.is_size_aware = false;
  decps.try_top_decomposition = true;
  decps.try_bottom_decomposition = false;
  decps.use_correlation = false;
  decps.try_xor = true;
  decps.nImpurity = ps.nImpurity;

  std::vector<kitty::partial_truth_table> examples;
  for( auto sim : ntk.sim_patterns )
    examples.push_back( sim.pat );

  signal<Ntk> osignal = sim_decomposition_fastS( ntk, examples, ntk.targets[ps.output], decps, false );
  double accuracy = 100*(double)kitty::count_ones( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) )/ntk.targets[ps.output].num_bits();

  if( ps.verbose )
  {
    std::cout << "[o " << ps.output << "] : " << accuracy << "%" <<std::endl;
    kitty::print_binary( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) );
  }

  return osignal;
}

/*! \brief 
 * Statistics based decomposition: tries top decomposition and performs shannon decomposition in case of failure
 */
template<class Ntk>
signal<Ntk> DK_XTSDS( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
{
  sim_decomposition_fastS_params decps;
  decps.verbose = ps.verbose;
  decps.max_sup = ps.max_sup;
  decps.is_informed = true;
  decps.is_size_aware = true;
  decps.try_top_decomposition = true;
  decps.try_bottom_decomposition = false;
  decps.use_correlation = false;
  decps.try_xor = true;

  std::vector<kitty::partial_truth_table> examples;
  for( auto sim : ntk.sim_patterns )
    examples.push_back( sim.pat );

  signal<Ntk> osignal = sim_decomposition_fastS( ntk, examples, ntk.targets[ps.output], decps, false );
  double accuracy = 100*(double)kitty::count_ones( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) )/ntk.targets[ps.output].num_bits();

  if( ps.verbose )
  {
    std::cout << "[o " << ps.output << "] : " << accuracy << "%" <<std::endl;
    kitty::print_binary( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) );
  }

  return osignal;
}

/*! \brief 
 * Statistics based decomposition: tries top decomposition and performs shannon decomposition in case of failure
 */
template<class Ntk>
signal<Ntk> DK_TSDS( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
{
  sim_decomposition_fastS_params decps;
  decps.verbose = ps.verbose;
  decps.max_sup = ps.max_sup;
  decps.is_informed = true;
  decps.is_size_aware = true;
  decps.try_top_decomposition = true;
  decps.try_bottom_decomposition = false;
  decps.use_correlation = false;
  decps.try_xor = false;

  std::vector<kitty::partial_truth_table> examples;
  for( auto sim : ntk.sim_patterns )
    examples.push_back( sim.pat );

  signal<Ntk> osignal = sim_decomposition_fastS( ntk, examples, ntk.targets[ps.output], decps, false );
  double accuracy = 100*(double)kitty::count_ones( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) )/ntk.targets[ps.output].num_bits();

  if( ps.verbose )
  {
    std::cout << "[o " << ps.output << "] : " << accuracy << "%" <<std::endl;
    kitty::print_binary( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) );
  }

  return osignal;
}

/*! \brief 
 * Statistics based decomposition: tries top decomposition and performs shannon decomposition in case of failure
 */
template<class Ntk>
signal<Ntk> DK_DSD( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
{
  sim_decomposition_fastS_params decps;
  decps.verbose = ps.verbose;
  decps.max_sup = ps.max_sup;
  decps.is_informed = true;
  decps.is_size_aware = false;
  decps.try_top_decomposition = true;
  decps.try_bottom_decomposition = true;
  decps.use_correlation = false;
  decps.try_xor = true;
  decps.is_relaxed = false;

  std::vector<kitty::partial_truth_table> examples;
  for( auto sim : ntk.sim_patterns )
    examples.push_back( sim.pat );

  signal<Ntk> osignal = sim_decomposition_fastS( ntk, examples, ntk.targets[ps.output], decps, false );
  double accuracy = 100*(double)kitty::count_ones( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) )/ntk.targets[ps.output].num_bits();

  if( ps.verbose )
  {
    std::cout << "[o " << ps.output << "] : " << accuracy << "%" <<std::endl;
    kitty::print_binary( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) );
  }

  return osignal;
}

/*! \brief 
 * Statistics based decomposition: tries top decomposition and performs shannon decomposition in case of failure
 */
template<class Ntk>
signal<Ntk> DK_RDSD( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
{
  sim_decomposition_fastS_params decps;
  decps.verbose = ps.verbose;
  decps.max_sup = ps.max_sup;
  decps.is_informed = true;
  decps.is_size_aware = false;
  decps.try_top_decomposition = true;
  decps.try_bottom_decomposition = true;
  decps.use_correlation = false;
  decps.try_xor = true;
  decps.is_relaxed = true;

  std::vector<kitty::partial_truth_table> examples;
  for( auto sim : ntk.sim_patterns )
    examples.push_back( sim.pat );

  signal<Ntk> osignal = sim_decomposition_fastS( ntk, examples, ntk.targets[ps.output], decps, false );
  double accuracy = 100*(double)kitty::count_ones( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) )/ntk.targets[ps.output].num_bits();

  if( ps.verbose )
  {
    std::cout << "[o " << ps.output << "] : " << accuracy << "%" <<std::endl;
    kitty::print_binary( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) );
  }

  return osignal;
}

/*! \brief 
 * Statistics based decomposition: tries top decomposition and performs shannon decomposition in case of failure
 */
template<class Ntk>
signal<Ntk> DK_X( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
{
  sim_decomposition_xor_params decps;
  decps.verbose = ps.verbose;
  decps.max_sup = ps.max_sup;

  std::vector<kitty::partial_truth_table> examples;
  for( auto sim : ntk.sim_patterns )
    examples.push_back( sim.pat );

  signal<Ntk> osignal = sim_decomposition_xor( ntk, examples, ntk.targets[ps.output], decps, false );
  double accuracy = 100*(double)kitty::count_ones( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) )/ntk.targets[ps.output].num_bits();

  if( ps.verbose )
  {
    std::cout << "[o " << ps.output << "] : " << accuracy << "%" <<std::endl;
    kitty::print_binary( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) );
  }

  return osignal;
}

/*! \brief 
 * Statistics based decomposition: tries top decomposition and performs shannon decomposition in case of failure
 */
template<class Ntk>
signal<Ntk> DC_TSD( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
{
  dc_decomposition_fastS_params decps;
  decps.verbose = ps.verbose;
  decps.max_sup = ps.max_sup;
  decps.is_size_aware = false;
  decps.try_top_decomposition = true;
  decps.use_correlation = false;
  decps.try_xor = false;

  std::vector<kitty::partial_truth_table> examples;
  for( auto sim : ntk.sim_patterns )
    examples.push_back( sim.pat );

  signal<Ntk> osignal = dc_decomposition_fastS( ntk, examples, ntk.targets[ps.output], decps, false );
  double accuracy = 100*(double)kitty::count_ones( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) )/ntk.targets[ps.output].num_bits();

  if( ps.verbose )
  {
    std::cout << "[o " << ps.output << "] : " << accuracy << "%" <<std::endl;
    kitty::print_binary( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) );
  }

  return osignal;
}

/*! \brief 
 * Statistics based decomposition: tries top decomposition and performs shannon decomposition in case of failure
 */
template<class Ntk>
signal<Ntk> DC_IXTSD( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
{
  sim_decomposition_fastS_params decps;
  decps.verbose = ps.verbose;
  decps.max_sup = ps.max_sup;
  decps.is_informed = true;
  decps.is_size_aware = false;
  decps.try_top_decomposition = true;
  decps.try_bottom_decomposition = false;
  decps.use_correlation = false;
  decps.try_xor = true;
  decps.is_dc = true;

  std::vector<kitty::partial_truth_table> examples;
  for( auto sim : ntk.sim_patterns )
    examples.push_back( sim.pat );

  signal<Ntk> osignal = sim_decomposition_fastS( ntk, examples, ntk.targets[ps.output], decps, false );
  double accuracy = 100*(double)kitty::count_ones( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) )/ntk.targets[ps.output].num_bits();

  if( ps.verbose )
  {
    std::cout << "[o " << ps.output << "] : " << accuracy << "%" <<std::endl;
    kitty::print_binary( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) );
  }

  return osignal;
}

/*! \brief 
 * Statistics based decomposition: tries top decomposition and performs shannon decomposition in case of failure
 */
template<class Ntk>
signal<Ntk> DC_XTSD( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
{
  dc_decomposition_fastS_params decps;
  decps.verbose = ps.verbose;
  decps.max_sup = ps.max_sup;
  decps.is_size_aware = false;
  decps.try_top_decomposition = true;
  decps.use_correlation = false;
  decps.try_xor = true;

  std::vector<kitty::partial_truth_table> examples;
  for( auto sim : ntk.sim_patterns )
    examples.push_back( sim.pat );

  signal<Ntk> osignal = dc_decomposition_fastS( ntk, examples, ntk.targets[ps.output], decps, false );
  double accuracy = 100*(double)kitty::count_ones( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) )/ntk.targets[ps.output].num_bits();

  if( ps.verbose )
  {
    std::cout << "[o " << ps.output << "] : " << accuracy << "%" <<std::endl;
    kitty::print_binary( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) );
  }

  return osignal;
}

/*! \brief 
 * Statistics based decomposition: tries top decomposition and performs shannon decomposition in case of failure
 */
/*template<class Ntk>
signal<Ntk> DK_TSDS( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
{
  sim_decomposition_fastS_params decps;
  decps.verbose = ps.verbose;
  decps.max_sup = ps.max_sup;
  decps.is_informed = true;
  decps.is_size_aware = true;
  decps.try_top_decomposition = true;
  decps.try_bottom_decomposition = false;
  decps.use_correlation = false;

  std::vector<kitty::partial_truth_table> examples;
  for( auto sim : ntk.sim_patterns )
    examples.push_back( sim.pat );

  signal<Ntk> osignal = sim_decomposition_fastS( ntk, examples, ntk.targets[ps.output], decps, false );
  double accuracy = 100*(double)kitty::count_ones( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) )/ntk.targets[ps.output].num_bits();

  if( ps.verbose )
  {
    std::cout << "[o " << ps.output << "] : " << accuracy << "%" <<std::endl;
    kitty::print_binary( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) );
  }

  return osignal;
}*/

/*! \brief 
 * Statistics based decomposition: tries top decomposition and performs shannon decomposition in case of failure*/
template<class Ntk>
signal<Ntk> DK_DSDS( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
{
  sim_decomposition_fastS_params decps;
  decps.verbose = ps.verbose;
  decps.max_sup = ps.max_sup;
  decps.is_informed = true;
  decps.is_size_aware = false;
  decps.try_top_decomposition = true;
  decps.try_bottom_decomposition = true;
  decps.use_correlation = false;

  std::vector<kitty::partial_truth_table> examples;
  for( auto sim : ntk.sim_patterns )
    examples.push_back( sim.pat );

  signal<Ntk> osignal = sim_decomposition_fastS( ntk, examples, ntk.targets[ps.output], decps, false );
  double accuracy = 100*(double)kitty::count_ones( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) )/ntk.targets[ps.output].num_bits();

  if( ps.verbose )
  {
    std::cout << "[o " << ps.output << "] : " << accuracy << "%" <<std::endl;
    kitty::print_binary( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) );
  }

  return osignal;
}

/*! \brief 
 * Statistics based decomposition: tries top decomposition and performs shannon decomposition in case of failure
 */
template<class Ntk>
signal<Ntk> forestS( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
{
  forest_decomposition_params decps;
  decps.verbose = ps.verbose;
  decps.max_sup = ps.max_sup;
  decps.is_informed = true;
  decps.is_size_aware = true;
  decps.try_top_decomposition = true;
  decps.try_bottom_decomposition = false;
  decps.use_correlation = false;
  decps.try_xor = false;
  decps.num_trees = ps.num_trees;

  std::vector<kitty::partial_truth_table> examples;
  for( auto sim : ntk.sim_patterns )
    examples.push_back( sim.pat );

  signal<Ntk> osignal = forest_decomposition( ntk, examples, ntk.targets[ps.output], decps, false );
  double accuracy = 100*(double)kitty::count_ones( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) )/ntk.targets[ps.output].num_bits();

  if( ps.verbose )
  {
    std::cout << "[o " << ps.output << "] : " << accuracy << "%" <<std::endl;
    kitty::print_binary( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) );
  }

  return osignal;
}

/*! \brief 
 * Statistics based decomposition: tries top decomposition and performs shannon decomposition in case of failure
 */
template<class Ntk>
signal<Ntk> forestSx2( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
{
  forest_decompositionx2_params decps;
  decps.verbose = ps.verbose;
  decps.max_sup = ps.max_sup;
  decps.is_informed = true;
  decps.is_size_aware = true;
  decps.try_top_decomposition = true;
  decps.try_bottom_decomposition = false;
  decps.use_correlation = false;
  decps.try_xor = false;
  decps.num_trees = ps.num_trees;

  std::vector<kitty::partial_truth_table> examples;
  for( auto sim : ntk.sim_patterns )
    examples.push_back( sim.pat );

  signal<Ntk> osignal = forest_decompositionx2( ntk, examples, ntk.targets[ps.output], decps, false );
  double accuracy = 100*(double)kitty::count_ones( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) )/ntk.targets[ps.output].num_bits();

  if( ps.verbose )
  {
    std::cout << "[o " << ps.output << "] : " << accuracy << "%" <<std::endl;
    kitty::print_binary( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) );
  }

  return osignal;
}



/*! \brief 
 * Statistics based decomposition: tries top decomposition and performs shannon decomposition in case of failure
 */
template<class Ntk>
signal<Ntk> xforestS( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
{
  forest_decomposition_params decps;
  decps.verbose = ps.verbose;
  decps.max_sup = ps.max_sup;
  decps.is_informed = true;
  decps.is_size_aware = true;
  decps.try_top_decomposition = true;
  decps.try_bottom_decomposition = false;
  decps.use_correlation = false;
  decps.try_xor = true;
  decps.num_trees = ps.num_trees;

  std::vector<kitty::partial_truth_table> examples;
  for( auto sim : ntk.sim_patterns )
    examples.push_back( sim.pat );

  signal<Ntk> osignal = forest_decomposition( ntk, examples, ntk.targets[ps.output], decps, false );
  double accuracy = 100*(double)kitty::count_ones( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) )/ntk.targets[ps.output].num_bits();

  if( ps.verbose )
  {
    std::cout << "[o " << ps.output << "] : " << accuracy << "%" <<std::endl;
    kitty::print_binary( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) );
  }

  return osignal;
}

/*! \brief 
 * Statistics based decomposition: tries top decomposition and performs shannon decomposition in case of failure
 */
template<class Ntk>
signal<Ntk> xforestSx2( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
{
  forest_decompositionx2_params decps;
  decps.verbose = ps.verbose;
  decps.max_sup = ps.max_sup;
  decps.is_informed = true;
  decps.is_size_aware = true;
  decps.try_top_decomposition = true;
  decps.try_bottom_decomposition = false;
  decps.use_correlation = false;
  decps.try_xor = true;
  decps.num_trees = ps.num_trees;

  std::vector<kitty::partial_truth_table> examples;
  for( auto sim : ntk.sim_patterns )
    examples.push_back( sim.pat );

  signal<Ntk> osignal = forest_decompositionx2( ntk, examples, ntk.targets[ps.output], decps, false );
  double accuracy = 100*(double)kitty::count_ones( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) )/ntk.targets[ps.output].num_bits();

  if( ps.verbose )
  {
    std::cout << "[o " << ps.output << "] : " << accuracy << "%" <<std::endl;
    kitty::print_binary( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) );
  }

  return osignal;
}

} // namespace detail

template<class Ntk>
signal<Ntk> recover_accuracy( simulation_view<Ntk>& ntk, detail::arecovery_method const& arecovery_m, detail::arecovery_params const& arecovery_ps )
{ 
    signal<Ntk> osignal;
    switch(arecovery_m) {
      case detail::arecovery_method::none:
        osignal = detail::best_node( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::SD:
        osignal = detail::SD( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::DK_SD:
        osignal = detail::DK_SD( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::DK_TSD:
        osignal = detail::DK_TSD( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::DK_XTSD:
        osignal = detail::DK_XTSD( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::DK_TSDS:
        osignal = detail::DK_XTSDS( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::DK_XTSDS:
        osignal = detail::DK_XTSDS( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::DC_TSD:
        osignal = detail::DC_TSD( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::DC_XTSD:
        osignal = detail::DC_XTSD( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::DC_IXTSD:
        osignal = detail::DC_IXTSD( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::DK_DSD:
        osignal = detail::DK_DSD( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::DK_RDSD:
        osignal = detail::DK_RDSD( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::DK_X:
        osignal = detail::DK_X( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::forestS:
        osignal = detail::forestS( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::forestSx2:
        osignal = detail::forestSx2( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::xforestS:
        osignal = detail::xforestS( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::xforestSx2:
        osignal = detail::xforestSx2( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::DK_DSDS:
        osignal = detail::DK_DSDS( ntk, arecovery_ps );
        break;
    }

    return osignal; 
}
} /* namespace hyperdimensional computing */
} // namespace mockturtle
