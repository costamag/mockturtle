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
  sdec,
  isdec,
  itsdec,
  ixtsdec,
  ixtsdecS,
  dcsdec,
  dcxsdec,
  itsdecS,
  itdsdec,
  forestS,
  forestSx2,
  idsdS
};

class arecovery_params
{
  public:
    uint32_t output{0};
    bool verbose{true};
    uint32_t max_sup{2};

    uint32_t num_trees{3};
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
signal<Ntk> sdec( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
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
signal<Ntk> isdec( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
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
signal<Ntk> itsdec( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
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
signal<Ntk> ixtsdec( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
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
signal<Ntk> ixtsdecS( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
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
signal<Ntk> idsdS( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
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
signal<Ntk> dcsdec( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
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
signal<Ntk> dcxsdec( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
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
template<class Ntk>
signal<Ntk> itsdecS( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
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
}

/*! \brief 
 * Statistics based decomposition: tries top decomposition and performs shannon decomposition in case of failure
 */
template<class Ntk>
signal<Ntk> itdsdec( simulation_view<Ntk>& ntk, detail::arecovery_params const& ps )
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
      case detail::arecovery_method::sdec:
        osignal = detail::sdec( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::isdec:
        osignal = detail::isdec( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::itsdec:
        osignal = detail::itsdec( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::ixtsdec:
        osignal = detail::ixtsdec( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::ixtsdecS:
        osignal = detail::ixtsdecS( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::dcsdec:
        osignal = detail::dcsdec( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::dcxsdec:
        osignal = detail::dcxsdec( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::itsdecS:
        osignal = detail::itsdec( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::itdsdec:
        osignal = detail::itdsdec( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::forestS:
        osignal = detail::forestS( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::forestSx2:
        osignal = detail::forestSx2( ntk, arecovery_ps );
        break;
      case detail::arecovery_method::idsdS:
        osignal = detail::idsdS( ntk, arecovery_ps );
        break;
    }

    return osignal; 
}
} /* namespace hyperdimensional computing */
} // namespace mockturtle
