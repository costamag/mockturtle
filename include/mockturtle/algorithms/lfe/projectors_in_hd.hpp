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
  \file iwls2022.hpp
  \brief Interfaces contest.cpp with the different initial synthesis methods.

  \author Andrea Costamagna
*/

#pragma once

#include <cstdint>
#include <vector>

#include "../../traits.hpp"
#include "../simulation.hpp"
#include "simulation_view.hpp"
#include "../../views/depth_view.hpp"
#include "sim_patterns.hpp"
#include "muesli.hpp"
#include "hyperdimensional_computing/methods/accuracy_recovery.hpp"
#include "hyperdimensional_computing/methods/generators.hpp"
#include "hyperdimensional_computing/methods/selectors.hpp"
#include "hyperdimensional_computing/methods/selgenerators.hpp"
#include "hyperdimensional_computing/model.hpp"

#include <fmt/format.h>
#include <kitty/constructors.hpp>
#include <kitty/decomposition.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/partial_truth_table.hpp>
#include <kitty/statistics.hpp>
#include <kitty/operators.hpp>
#include <kitty/print.hpp>

#include "chatterjee_method.hpp"


namespace mockturtle
{
namespace hdc
{

klut_network project_in_hd( std::vector<kitty::partial_truth_table> examples, std::vector<kitty::partial_truth_table> targets, int const& topology )
{

  klut_network oklut;
  simulation_view oklut_sim{ oklut };
  model M( oklut_sim, examples, targets );
  std::vector<signal<klut_network>> osignals;
  hdc::detail::creation_params creation_ps;
  hdc::detail::selection_params selection_ps;
  creation_ps.verbose=false;
  selection_ps.verbose=false;

  switch(topology)
  {
    case 0:
    {
      hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::sdec;
      hdc::detail::arecovery_params arecovery_ps;
      arecovery_ps.verbose = false;

      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        arecovery_ps.output=y;
        osignals.push_back(M.accuracy_recovery(arecovery_m,arecovery_ps));
      }
      break;
    }
    case 1:
    {
      hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::isdec;
      hdc::detail::arecovery_params arecovery_ps;
      arecovery_ps.verbose = false;

      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        arecovery_ps.output=y;
        osignals.push_back(M.accuracy_recovery(arecovery_m,arecovery_ps));
      }
      break;
    }
    case 2:
    {
      hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::itsdec;
      hdc::detail::arecovery_params arecovery_ps;
      arecovery_ps.verbose = false;

      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        arecovery_ps.output=y;
        osignals.push_back(M.accuracy_recovery(arecovery_m,arecovery_ps));
      }
      break;
    }
    case 3:
    {
      hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::ixtsdec;
      hdc::detail::arecovery_params arecovery_ps;
      arecovery_ps.verbose = false;

      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        arecovery_ps.output=y;
        osignals.push_back(M.accuracy_recovery(arecovery_m,arecovery_ps));
      }
      break;
    }
    case 4:
    {
      hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::dcsdec;
      hdc::detail::arecovery_params arecovery_ps;
      arecovery_ps.verbose = false;

      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        arecovery_ps.output=y;
        osignals.push_back(M.accuracy_recovery(arecovery_m,arecovery_ps));
      }
      break;
    }
    case 5:
    {
      hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::dcxsdec;
      hdc::detail::arecovery_params arecovery_ps;
      arecovery_ps.verbose = false;

      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        arecovery_ps.output=y;
        osignals.push_back(M.accuracy_recovery(arecovery_m,arecovery_ps));
      }
      break;
    }
    case 6:
    {
      hdc::detail::selcreation_method selcreation_m = hdc::detail::selcreation_method::muesli;
      hdc::detail::selcreation_params selcreation_ps;
      selcreation_ps.re_initialize = true;
      
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        selcreation_ps.output=y;
        M.add( selcreation_m, selcreation_ps );
      }

      hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::none;
      hdc::detail::arecovery_params arecovery_ps;
      arecovery_ps.verbose = false;
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        arecovery_ps.output=y;
        osignals.push_back( M.accuracy_recovery(arecovery_m, arecovery_ps) );
      }
    break;
    }
    case 7:
    {
      hdc::detail::selcreation_method selcreation_m = hdc::detail::selcreation_method::muesli;
      hdc::detail::selcreation_params selcreation_ps;
      selcreation_ps.re_initialize = false;
      
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        selcreation_ps.output=y;
        M.add( selcreation_m, selcreation_ps );
      }

      hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::ixtsdec;
      hdc::detail::arecovery_params arecovery_ps;
      arecovery_ps.verbose = false;
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        arecovery_ps.output=y;
        osignals.push_back( M.accuracy_recovery(arecovery_m, arecovery_ps) );
      }
    break;
    }
    case 8:
    {
      hdc::detail::selcreation_method selcreation_m = hdc::detail::selcreation_method::sim_muesli;
      hdc::detail::selcreation_params selcreation_ps;
      selcreation_ps.re_initialize = false;
      
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        selcreation_ps.output=y;
        M.add( selcreation_m, selcreation_ps );
      }

      hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::ixtsdec;
      hdc::detail::arecovery_params arecovery_ps;
      arecovery_ps.verbose = false;
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        arecovery_ps.output=y;
        osignals.push_back( M.accuracy_recovery(arecovery_m, arecovery_ps) );
      }
    break;
    }
    case 9:
    {
      hdc::detail::selection_method selection_m = hdc::detail::selection_method::depth_selector;
      hdc::detail::selection_params selection_ps;
      selection_ps.max_new_supports = 1024;
      selection_ps.max_selection_attempts = 10000;
      selection_ps.support_size = 2;
      selection_ps.max_search_depth = std::numeric_limits<uint32_t>::max();
      hdc::detail::creation_method creation_m = hdc::detail::creation_method::fgenerator1;
      hdc::detail::creation_params creation_ps;
      creation_ps.max_nodes_total = 1024;
      creation_ps.max_nodes_support = 1;

      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        creation_ps.output=y;
        M.add( selection_m, selection_ps, creation_m, creation_ps );
      }

      hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::ixtsdec;
      hdc::detail::arecovery_params arecovery_ps;
      arecovery_ps.verbose = false;
      std::cout << M.ntk_.num_gates() << " ";
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        arecovery_ps.output=y;
        osignals.push_back( M.accuracy_recovery(arecovery_m, arecovery_ps) );
      }
      std::cout << M.ntk_.num_gates() << std::endl;

    break;
    }
    case 10:
    {
      hdc::detail::selection_method selection_m = hdc::detail::selection_method::depth_selector;
      hdc::detail::selection_params selection_ps;
      selection_ps.max_new_supports = 1024;
      selection_ps.max_selection_attempts = 10000;
      selection_ps.support_size = 2;
      selection_ps.max_search_depth = std::numeric_limits<uint32_t>::max();
      hdc::detail::creation_method creation_m = hdc::detail::creation_method::ifgenerator1;
      hdc::detail::creation_params creation_ps;
      creation_ps.max_nodes_total  = 1024;


      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        creation_ps.output=y;
        M.add( selection_m, selection_ps, creation_m, creation_ps );
      }
      std::cout << M.ntk_.num_gates() << " ";
      hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::ixtsdec;
      hdc::detail::arecovery_params arecovery_ps;
      arecovery_ps.verbose = false;
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        arecovery_ps.output=y;
        osignals.push_back( M.accuracy_recovery(arecovery_m, arecovery_ps) );
      }
      std::cout << M.ntk_.num_gates() << std::endl;

    break;
    }
    /*case 6:
    {
      hdc::detail::selcreation_method selcreation_m = hdc::detail::selcreation_method::sim_muesli;
      hdc::detail::selcreation_params selcreation_ps;
      selcreation_ps.re_initialize = false;
      
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        selcreation_ps.output=y;
        osignals.push_back(M.add( selcreation_m, selcreation_ps ));
      }

    break;
    }
    case 7:
    {
      hdc::detail::selcreation_method selcreation_m = hdc::detail::selcreation_method::sim_muesli;
      hdc::detail::selcreation_params selcreation_ps;
      selcreation_ps.re_initialize = false;
      
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        selcreation_ps.output=y;
        M.add( selcreation_m, selcreation_ps );
      }


      hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::itsdec;
      hdc::detail::arecovery_params arecovery_ps;
      arecovery_ps.verbose = false;
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        arecovery_ps.output=y;
        osignals.push_back( M.accuracy_recovery(arecovery_m, arecovery_ps) );
      }
    break;
    }
    case 8:
    {
      hdc::detail::selection_method selection_m = hdc::detail::selection_method::urandom;
      hdc::detail::selection_params selection_ps;
      selection_ps.max_new_nodes = 500;
      selection_ps.support_size = 4;
      selection_ps.max_attempts = 1000;
      
      hdc::detail::creation_method creation_m = hdc::detail::creation_method::chatterjee1;
      hdc::detail::creation_params creation_ps;
      
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        creation_ps.output=y;
        creation_ps.max_new_nodes = 500;
        
        for( uint32_t l = 0; l < 5; ++l )
        {
          M.add( selection_m, selection_ps, creation_m, creation_ps );
        }
      }

      hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::itsdec;
      hdc::detail::arecovery_params arecovery_ps;
      arecovery_ps.verbose = false;
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        arecovery_ps.output=y;
        osignals.push_back( M.accuracy_recovery(arecovery_m, arecovery_ps) );
      }
    break;
    }
    case 9:
    {
      hdc::detail::selection_method selection_m = hdc::detail::selection_method::urandom;
      hdc::detail::selection_params selection_ps;
      selection_ps.max_new_nodes = 500;
      selection_ps.support_size = 4;
      selection_ps.max_attempts = 1000;
      
      hdc::detail::creation_method creation_m = hdc::detail::creation_method::chatterjee1;
      hdc::detail::creation_params creation_ps;
      
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        creation_ps.output=y;
        creation_ps.max_new_nodes = 500;
        
        for( uint32_t l = 0; l < 5; ++l )
        {
          M.add( selection_m, selection_ps, creation_m, creation_ps );
        }
      }

      hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::itsdecS;
      hdc::detail::arecovery_params arecovery_ps;
      arecovery_ps.verbose = false;
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        arecovery_ps.output=y;
        osignals.push_back( M.accuracy_recovery(arecovery_m, arecovery_ps) );
      }
    break;
    }
    case 10:
    {
      hdc::detail::selection_method selection_m = hdc::detail::selection_method::urandom;
      hdc::detail::selection_params selection_ps;
      selection_ps.max_new_nodes = 1028;
      selection_ps.support_size = 2;
      selection_ps.max_attempts = 10000;
      
      hdc::detail::creation_method creation_m = hdc::detail::creation_method::chatterjee1;
      hdc::detail::creation_params creation_ps;
      
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        creation_ps.output=y;
        creation_ps.max_new_nodes = 1028;
        
        for( uint32_t l = 0; l < 5; ++l )
        {
          M.add( selection_m, selection_ps, creation_m, creation_ps );
        }
      }
      std::cout << M.ntk_.num_gates() << std::endl;

      hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::itsdecS;
      hdc::detail::arecovery_params arecovery_ps;
      arecovery_ps.verbose = false;
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        arecovery_ps.output=y;
        osignals.push_back( M.accuracy_recovery(arecovery_m, arecovery_ps) );
      }
    break;
    }
    case 11:
    {
      hdc::detail::selection_method selection_m = hdc::detail::selection_method::urandom;
      hdc::detail::selection_params selection_ps;
      selection_ps.max_new_nodes = 1024;
      selection_ps.max_attempts = 10000;
      creation_ps.max_new_nodes = 1024;

      selection_ps.support_size = 2;

      
      hdc::detail::creation_method creation_m = hdc::detail::creation_method::assembler1;
      hdc::detail::creation_params creation_ps;
      
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        creation_ps.output=y;
        //for( uint32_t l = 0; l < 10; ++l )
        //{
          M.add( selection_m, selection_ps, creation_m, creation_ps );
        //}
      }
      std::cout << M.ntk_.num_gates() << std::endl;

      hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::itsdecS;
      hdc::detail::arecovery_params arecovery_ps;
      arecovery_ps.verbose = false;
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        arecovery_ps.output=y;
        osignals.push_back( M.accuracy_recovery(arecovery_m, arecovery_ps) );
      }
    break;
    }
    case 12:
    {
      hdc::detail::selection_method selection_m = hdc::detail::selection_method::urandom;
      hdc::detail::selection_params selection_ps;
      selection_ps.max_new_nodes = 1028;
      selection_ps.support_size = 2;
      selection_ps.max_attempts = 10000;
      
      hdc::detail::creation_method creation_m = hdc::detail::creation_method::chatterjee1;
      hdc::detail::creation_params creation_ps;
      
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        creation_ps.output=y;
        creation_ps.max_new_nodes = 1028;
        
        for( uint32_t l = 0; l < 10; ++l )
        {
          M.add( selection_m, selection_ps, creation_m, creation_ps );
        }
      }
      std::cout << M.ntk_.num_gates() << std::endl;

      hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::itsdecS;
      hdc::detail::arecovery_params arecovery_ps;
      arecovery_ps.verbose = false;
      arecovery_ps.max_sup = 8;

      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        arecovery_ps.output=y;
        osignals.push_back( M.accuracy_recovery(arecovery_m, arecovery_ps) );
      }
    break;
    }
    case 13:
    {
      hdc::detail::selection_method selection_m = hdc::detail::selection_method::urandom;
      hdc::detail::selection_params selection_ps;
      selection_ps.max_new_nodes = 1028;
      selection_ps.support_size = 2;
      selection_ps.max_attempts = 10000;
      
      hdc::detail::creation_method creation_m = hdc::detail::creation_method::chatterjee1;
      hdc::detail::creation_params creation_ps;
      
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        creation_ps.output=y;
        creation_ps.max_new_nodes = 1028;
        
        for( uint32_t l = 0; l < 10; ++l )
        {
          M.add( selection_m, selection_ps, creation_m, creation_ps );
        }
      }
      std::cout << M.ntk_.num_gates() << std::endl;

      hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::itsdecS;
      hdc::detail::arecovery_params arecovery_ps;
      arecovery_ps.verbose = false;
      arecovery_ps.max_sup = 15;

      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        arecovery_ps.output=y;
        osignals.push_back( M.accuracy_recovery(arecovery_m, arecovery_ps) );
      }
    break;
    }
    case 14:
    {
      hdc::detail::selection_method selection_m = hdc::detail::selection_method::urandom;
      hdc::detail::selection_params selection_ps;
      selection_ps.max_new_nodes = 1024;
      selection_ps.support_size = 2;
      selection_ps.max_attempts = 10000;
      
      hdc::detail::creation_method creation_m = hdc::detail::creation_method::assembler1;
      hdc::detail::creation_params creation_ps;
      
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        creation_ps.output=y;
        creation_ps.max_new_nodes = 1024;
        M.add( selection_m, selection_ps, creation_m, creation_ps );
      }
      std::cout << M.ntk_.num_gates() << std::endl;

      hdc::detail::arecovery_method arecovery_m = hdc::detail::arecovery_method::itsdecS;
      hdc::detail::arecovery_params arecovery_ps;
      arecovery_ps.verbose = false;
      for( uint32_t y = 0; y < targets.size(); ++y )
      {
        arecovery_ps.output=y;
        osignals.push_back( M.accuracy_recovery(arecovery_m, arecovery_ps) );
      }
    break;
    }*/
  }

  for( size_t i = 0; i < osignals.size(); ++i )
  {
    oklut_sim.create_po(osignals[i]);
  }

  return oklut_sim;
}
}
}