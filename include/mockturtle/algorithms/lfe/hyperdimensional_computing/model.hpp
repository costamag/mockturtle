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
  \file model.hpp
  \brief Data structure to treat logic networks as a machine learning model
  \author Andrea Costamagna
*/

#pragma once

#include "../../../networks/klut.hpp"
#include "../simulation_view.hpp"
#include "../../../utils/node_map.hpp"
#include "../sim_patterns.hpp"
#include "methods/selectors.hpp"
#include "methods/generators.hpp"
#include "methods/selgenerators.hpp"
#include "methods/accuracy_recovery.hpp"
#include <kitty/partial_truth_table.hpp>
#include <kitty/print.hpp>

namespace mockturtle
{
namespace hdc
{
template<typename Ntk>
class model
{
public:
  using storage = typename Ntk::storage;
  using node    = typename Ntk::node;
  using signal  = typename Ntk::signal;
  using TT = kitty::partial_truth_table;

public:
model( simulation_view<Ntk>& ntk, std::vector<TT> & X, std::vector<TT> & Y )
: ntk_(ntk)
{
  ntk_.initialize_network( X );
  ntk_.targets = Y;
  ntk_.layer_pointer = 1;
}
#pragma end region initialization

#pragma region reposition pointer
void reposition_pointer( uint32_t new_layer_pointer = 1 )
{
  assert( new_layer_pointer > 0 );
  ntk_.layer_pointer = new_layer_pointer;
}
#pragma endregion reposition pointer

#pragma region add a layer
void add( hdc::detail::selection_method const& selection_m, detail::selection_params const& selection_ps,
          detail::creation_method const& creation_m, detail::creation_params const& creation_ps )
{
  std::vector<std::vector<signal>> divisors = select_variables( ntk_, selection_m, selection_ps );
  create_nodes( ntk_, divisors, creation_m, creation_ps );
}

signal add( detail::selcreation_method const& selcreation_m, detail::selcreation_params const& selcreation_ps )
{
  return selcreate_nodes( ntk_, selcreation_m, selcreation_ps );
}
#pragma endregion add layer

#pragma region accuracy recovery
signal accuracy_recovery(detail::arecovery_method const& arecovery_m, detail::arecovery_params const& arecovery_ps)
{
  return recover_accuracy( ntk_, arecovery_m, arecovery_ps );
}
#pragma endregion accuracy recovery

#pragma region print summary
void print_summary()
{
    std::cout << "============== SUMMARY: =============" << std::endl;

  for( auto i {0}; i < ntk_.summary.size(); ++i )
  {
    std::cout << "=============== LAYER" << i << " ==============" << std::endl;
    std::cout << ntk_.summary[i] << std::endl;
  }
    std::cout << "=====================================" << std::endl;

}
#pragma endregion print summary

public:
  simulation_view<Ntk>& ntk_;
  
};

} /* namespace hyperdimensional computing */
} /* namespace mockturtle */
