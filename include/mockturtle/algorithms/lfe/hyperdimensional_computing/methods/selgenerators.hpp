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
  \file selgenerators.hpp
  \brief Methods alllowing to select variables and to generate the nodes.
  \author Andrea Costamagna
*/

#pragma once

#include "../../chatterjee_method.hpp"
#include "../../muesli.hpp"
#include "../../sim_muesli.hpp"
#include "../../sim_create_nodes.hpp"
#include "../../sim_decomposition_fast.hpp"
#include "selectors.hpp"
#include <kitty/print.hpp>
#include <kitty/properties.hpp>

namespace mockturtle
{
namespace hdc
{
namespace detail
{
enum class selcreation_method
{
  muesli,
  sim_muesli
};

class selcreation_params
{
  public:
    uint32_t output{0};
    bool recover_accuracy{false};
    bool verbose{false};
    bool re_initialize{false};
    uint32_t max_act{5};
};

template<class Ntk>
signal<Ntk> muesli( simulation_view<Ntk>& ntk, selcreation_params const& ps )
{
  muesli_params muesli_ps;
  muesli_ps.verbose = ps.verbose;
  muesli_ps.try_accuracy_recovery = ps.recover_accuracy;
  muesli_ps.re_initialize = ps.re_initialize;
  muesli_ps.max_act = ps.max_act;
  std::vector<kitty::partial_truth_table> examples;
  for( auto sim : ntk.sim_patterns )
    examples.push_back( sim.pat );

  signal<Ntk> osignal = muesli( ntk, examples, ntk.targets[ps.output], muesli_ps );
  if( ps.verbose )
  {
    double accuracy = 100*(double)kitty::count_ones( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) )/ntk.targets[ps.output].num_bits();
    std::cout << "[o " << ps.output << "] : " << accuracy << "%" <<std::endl;
  }
  return osignal;
}

template<class Ntk>
signal<Ntk> sim_muesli( simulation_view<Ntk>& ntk, selcreation_params const& ps )
{
  sim_muesli_params muesli_ps;
  muesli_ps.verbose = ps.verbose;
  muesli_ps.try_accuracy_recovery = ps.recover_accuracy;
  muesli_ps.re_initialize = ps.re_initialize;
  std::vector<kitty::partial_truth_table> examples;
  for( auto sim : ntk.sim_patterns )
    examples.push_back( sim.pat );

  signal<Ntk> osignal = sim_muesli( ntk, examples, ntk.targets[ps.output], muesli_ps );
  if( ps.verbose )
  {
    double accuracy = 100*(double)kitty::count_ones( ~(ntk.targets[ps.output]^ntk.sim_patterns[ntk.nodes_to_patterns[osignal]].pat) )/ntk.targets[ps.output].num_bits();
    std::cout << "[o " << ps.output << "] : " << accuracy << "%" <<std::endl;
  }
  return osignal;
}

} // namespace detail

/*! \brief Select varibles and create nodes
 *
 */
template<class Ntk>
signal<Ntk> selcreate_nodes( simulation_view<Ntk>& ntk, detail::selcreation_method const& selcreation_m, detail::selcreation_params const& selcreation_ps )
{ 
  signal<Ntk> osignal;
  switch(selcreation_m) {
    case detail::selcreation_method::muesli:
      osignal = detail::muesli( ntk, selcreation_ps );
      break;
    case detail::selcreation_method::sim_muesli:
      osignal = detail::sim_muesli( ntk, selcreation_ps );
      break;
  }
  return osignal; 
}
} /* namespace hyperdimensional computing */
} // namespace mockturtle
