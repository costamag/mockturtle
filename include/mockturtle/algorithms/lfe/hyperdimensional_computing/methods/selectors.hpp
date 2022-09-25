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
  \file selectors.hpp
  \brief Select variables with which to build the new nodes.
  \author Andrea Costamagna
*/

#pragma once
#include <algorithm>
#include <random>
#include <array>
#include <vector>
#include <set>

#include <kitty/statistics.hpp>

namespace mockturtle
{
namespace hdc
{
namespace detail
{
enum class selection_method
{
  depth_selector,
  layer_selector,
  similarity_selector
};

class selection_params
{
  public:
    uint32_t max_new_supports{1};
    uint32_t max_selection_attempts{std::numeric_limits<uint32_t>::max()};
    uint32_t support_size{2};
    // used by depth selector
    uint32_t max_search_depth{std::numeric_limits<uint32_t>::max()};
    // used by layer selector
    uint32_t layer;
    // used by multilayer selector
    uint32_t min_layer;
    uint32_t max_layer;

    bool verbose{false};
};

template<class Ntk>
inline std::vector<std::vector<signal<Ntk>>> depth_selector( simulation_view<Ntk>& ntk, detail::selection_params const& ps )
{
  std::vector<std::vector<signal<Ntk>>> child_signals;
  uint32_t layer_pointer = ntk.layer_to_signals.size();
  uint32_t max_depth = std::min<uint32_t>( layer_pointer, ps.max_search_depth );
  uint32_t min_layer = layer_pointer-max_depth;
  uint32_t max_layer = layer_pointer-1;
  
  uint32_t layer;
  signal<Ntk> child;

  uint32_t num_supp = 0;
  uint32_t num_attempts = 0;
  std::set<std::vector<signal<Ntk>>> supports_set;
  srand(ntk.seed);
  while( ( num_supp < ps.max_new_supports ) && ( num_attempts < ps.max_selection_attempts ) )
  {
    std::vector<signal<Ntk>> child_i;

    for( uint32_t s{0}; s < ps.support_size; ++s )
    {
      layer = rand()%(max_layer-min_layer+1)+min_layer;

      child = rand()%(ntk.layer_to_signals[layer].size());

      if( ( std::find(child_i.begin(), child_i.end(), ntk.layer_to_signals[layer][child]) == child_i.end()) )
        child_i.push_back( ntk.layer_to_signals[layer][child] );
    }
    std::sort(child_i.begin(), child_i.end());

    if( ( child_i.size() > 1 ) && ( supports_set.find(child_i) == supports_set.end()) )
    {
      child_signals.push_back( child_i );
      supports_set.insert( child_i );
      num_supp++;
    }

    num_attempts++;
  }

  std::vector<std::vector<signal<Ntk>>> res;
  res = child_signals;
  ntk.seed*=res.size();
  return res;   
}

template<class Ntk>
inline std::vector<std::vector<signal<Ntk>>> layer_selector( simulation_view<Ntk>& ntk, detail::selection_params const& ps )
{
  assert( ps.layer < ntk.layer_to_signals.size() );

  std::vector<std::vector<signal<Ntk>>> child_signals;
  
  signal<Ntk> child;

  uint32_t num_supp = 0;
  uint32_t num_attempts = 0;
  std::set<std::vector<signal<Ntk>>> supports_set;
  srand(ntk.seed);
  while( ( num_supp < ps.max_new_supports ) && ( num_attempts < ps.max_selection_attempts ) )
  {
    std::vector<signal<Ntk>> child_i;

    for( uint32_t s{0}; s < ps.support_size; ++s )
    {

      child = rand()%(ntk.layer_to_signals[ps.layer].size());

      if( ( std::find(child_i.begin(), child_i.end(), ntk.layer_to_signals[ps.layer][child]) == child_i.end()) )
        child_i.push_back( ntk.layer_to_signals[ps.layer][child] );
    }
    std::sort(child_i.begin(), child_i.end());

    if( ( child_i.size() > 1 ) && ( supports_set.find(child_i) == supports_set.end()) )
    {
      child_signals.push_back( child_i );
      supports_set.insert( child_i );
      num_supp++;
    }

    num_attempts++;
  }

  std::vector<std::vector<signal<Ntk>>> res;
  res = child_signals;
  ntk.seed*=res.size();
  return res;   
}

template<class Ntk>
inline std::vector<std::vector<signal<Ntk>>> similarity_selector( simulation_view<Ntk>& ntk, detail::selection_params const& ps )
{
  std::vector<std::vector<signal<Ntk>>> child_signals;
  uint32_t layer_pointer = ntk.layer_to_signals.size();
  uint32_t max_depth = std::min<uint32_t>( layer_pointer, ps.max_search_depth );
  uint32_t min_layer = layer_pointer-max_depth;
  uint32_t max_layer = layer_pointer-1;
  std::unordered_map<double,std::vector<uint32_t>> i_to_sigs;
  
  double inew;
  uint32_t iref;
  std::vector<uint32_t> sorted_indeces;
  std::vector<double> sorted_mis;
  for( uint32_t i = 0; i < ntk.sim_patterns.size(); ++i )
  {
    inew = kitty::mutual_information( std::vector{&ntk.sim_patterns[i].pat}, &ntk.targets[0] );
    if( sorted_mis.size() == 0 || inew >= sorted_mis[sorted_mis.size()-1] )
    {
      sorted_mis.push_back( inew );
      sorted_indeces.push_back( i );
    }
    else
    {
      for( uint32_t j = 0; j < sorted_indeces.size(); ++j )
      {
        if( inew <= sorted_mis[j] )
        {
          sorted_mis.insert( sorted_mis.begin()+j, inew );
          sorted_indeces.insert( sorted_indeces.begin()+j, i );
          break;
        }
      }
    }
  }

  std::vector<std::vector<signal<Ntk>>> res;
  for( uint32_t i=0; i < sorted_mis.size()-1; ++i )
  {
    std::cout << sorted_mis[i] << " ";
    res.push_back( std::vector{ ntk.sim_patterns[sorted_indeces[i]].sig,ntk.sim_patterns[sorted_indeces[i+1]].sig } );
  }
  std::cout << std::endl;


  ntk.seed*=res.size();
  return res;   
}

/*template<class Ntk>
inline std::vector<std::vector<signal<Ntk>>> urandom_selector( simulation_view<Ntk>& ntk, detail::selection_params const& ps )
{
  std::vector<std::vector<signal<Ntk>>> child_signals;

  uint32_t max_depth = std::min<uint32_t>( ntk.layer_pointer, ps.max_search_depth );
  uint32_t min_layer = ntk.layer_pointer-max_depth;
  uint32_t max_layer = ntk.layer_pointer-1;
  uint32_t min_nodes = 0;
  uint32_t max_nodes;

  uint32_t layer;
  signal<Ntk> child;

  uint32_t n = 0;
  uint32_t c = 0;
  std::set<std::vector<signal<Ntk>>> supports_set;
  srand(ntk.seed);
  
  while( ( n < ps.max_new_nodes ) && ( c < ps.max_attempts ) )
  {
    std::vector<signal<Ntk>> child_i;

    for( uint32_t s{0}; s < ps.support_size; ++s )
    {
      layer = (rand()%(max_layer-min_layer+1))+min_layer;
      max_nodes = ntk.layer_to_signals[layer].size()-1;
      child = (rand()%(max_nodes-min_nodes+1))+min_nodes;

      if( ( std::find(child_i.begin(), child_i.end(), ntk.layer_to_signals[layer][child]) == child_i.end()) )
      {
        child_i.push_back( ntk.layer_to_signals[layer][child] );
      }
    }
    std::sort(child_i.begin(), child_i.end());

    if( ( child_i.size() > 1 ) && ( supports_set.find(child_i)==supports_set.end()) )
    {
      child_signals.push_back( child_i );
      supports_set.insert( child_i );
      n++;
    }
    c++;
  }
  std::vector<std::vector<signal<Ntk>>> res;
  res = child_signals;
  ntk.seed+=res.size();
  return res;   
}*/

} // namespace detail

/*! \brief select variables
 *
 */
template<class Ntk>
std::vector<std::vector<signal<Ntk>>> select_variables( simulation_view<Ntk>& ntk, detail::selection_method const& selection_m, detail::selection_params const& selection_ps )
{ 
    std::vector<std::vector<signal<Ntk>>> divisors;
    switch(selection_m) {
      case detail::selection_method::depth_selector:
        divisors = detail::depth_selector( ntk, selection_ps );
        break;
      case detail::selection_method::layer_selector:
        divisors = detail::layer_selector( ntk, selection_ps );
        break;
      case detail::selection_method::similarity_selector:
        divisors = detail::similarity_selector( ntk, selection_ps );
        break;

    }
    return divisors;
}
} // namespace hyperdimensional computing
} // namespace mockturtle
