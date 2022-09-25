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
  \file generators.hpp
  \brief Methods to create new nodes given the support variables.
  \author Andrea Costamagna
*/

#pragma once

#include "../../create_candidates.hpp"
#include "../../sim_create_nodes.hpp"
#include "../../chatterjee_method.hpp"
#include "selectors.hpp"
#include <kitty/print.hpp>
#include <kitty/properties.hpp>
#include <kitty/statistics.hpp>

namespace mockturtle
{
namespace hdc
{
namespace detail
{
enum class creation_method
{
  fgenerator1,
  ifgenerator1,
  xorgen,
  andgen,
  majgen,
  orthogonal_creator,
  chatterjee1,
  random
};

class creation_params
{
  public:
    uint32_t output{0};
    uint32_t max_nodes_total {std::numeric_limits<uint32_t>::max()};
    uint32_t max_nodes_support {std::numeric_limits<uint32_t>::max()};
    bool verbose{false};
};

template<class Ntk>
void fgenerator1( simulation_view<Ntk>& ntk, std::vector<std::vector<signal<Ntk>>> & supports, 
                           detail::creation_params const& ps )
{
  kitty::partial_truth_table * Y = &ntk.targets[ps.output];

  std::vector<signal<Ntk>> new_signals;
  uint32_t nodes_added_support;
  uint32_t nodes_added_total{0};

  for( uint32_t i=0; i < supports.size(); ++i )
  {

    if( nodes_added_total >= ps.max_nodes_total  ) break;
    if( supports[i].size()>1 )
    {

      std::vector<kitty::partial_truth_table*> X;
      for( auto d : supports[i] )
        X.push_back( &( ntk.sim_patterns[ntk.nodes_to_patterns[ntk.get_node(d)]].pat ) );

      create_candidates_result<kitty::partial_truth_table> Fset = create_candidates_method( X, Y );

      nodes_added_support = 0;
      for( uint32_t j{0}; j< Fset.dtt_v.size(); ++j )
      {
        auto f = std::make_pair( supports[i], Fset.tt_v[j] );
        if( ntk.available_nodes.find( f ) == ntk.available_nodes.end() )
        {
          ntk.create_node( supports[i], Fset.dtt_v[j] );
          ntk.available_nodes.insert( f );
          nodes_added_support++;
          nodes_added_total++;
          if( nodes_added_support >= ps.max_nodes_support ) break;
        }
      }
    }
  }
}

template<typename Ntk>
struct candidate_type{
  std::vector<signal<Ntk>> support;
  kitty::dynamic_truth_table dtt;
  std::string tt;
  double mi;
};

template<class Ntk>
void ifgenerator1( simulation_view<Ntk>& ntk, std::vector<std::vector<signal<Ntk>>> & supports, 
                           detail::creation_params const& ps )
{
  kitty::partial_truth_table * Y = &ntk.targets[ps.output];
  std::vector<candidate_type<Ntk>> candidates;
  candidate_type<Ntk> candidate;

  std::vector<signal<Ntk>> new_signals;
  uint32_t nodes_added_total{0};

  for( uint32_t i=0; i < supports.size(); ++i )
  {
    if( nodes_added_total >= ps.max_nodes_total  ) break;
    if( supports[i].size()>1 )
    {
      std::vector<kitty::partial_truth_table*> X;
      for( auto d : supports[i] )
        X.push_back( &( ntk.sim_patterns[ntk.nodes_to_patterns[ntk.get_node(d)]].pat ) );

      create_candidates_result<kitty::partial_truth_table> Fset = create_candidates_method( X, Y );

      for( uint32_t j{0}; j< Fset.dtt_v.size(); ++j )
      {
        auto f = std::make_pair( supports[i], Fset.tt_v[j] );
        if( ntk.available_nodes.find( f ) == ntk.available_nodes.end() )
        {
          candidate.support = supports[i];
          candidate.dtt = Fset.dtt_v[j];
          candidate.tt = Fset.tt_v[j];
          candidate.mi = kitty::mutual_information( X, Y );
          if( candidates.size() == 0 || ( candidates[candidates.size()-1].mi >= candidate.mi ) )
            candidates.push_back(candidate);
          else
          {
            for( uint32_t j{0}; j < candidates.size(); ++j )
            {
              if( candidate.mi > candidates[j].mi )
              {
                candidates.insert( candidates.begin() + j, candidate );
                break;
              }
            }
          }
        }
      }
    }
  }

  if( candidates.size() > ps.max_nodes_total )
  {
    for( uint32_t i{0}; candidates.size()-ps.max_nodes_total; ++i )
      candidates.erase( candidates.begin() + ps.max_nodes_total );
  }

  for( auto cand : candidates )
  {
    ntk.create_node( cand.support, cand.dtt );
    auto f = std::make_pair( cand.support, cand.tt );
    ntk.available_nodes.insert( f );
  }
}

template<class Ntk>
void andgen( simulation_view<Ntk>& ntk, std::vector<std::vector<signal<Ntk>>> & supports, 
                           detail::creation_params const& ps )
{
  for( uint32_t i=0; i < supports.size(); ++i )
    ntk.create_and( supports[i][0], supports[i][1] );
}


template<class Ntk>
void xorgen( simulation_view<Ntk>& ntk, std::vector<std::vector<signal<Ntk>>> & supports, 
                           detail::creation_params const& ps )
{
  for( uint32_t i=0; i < supports.size(); ++i )
    ntk.create_xor( supports[i][0], supports[i][1] );
}



template<class Ntk>
void majgen( simulation_view<Ntk>& ntk, std::vector<std::vector<signal<Ntk>>> & supports, 
                           detail::creation_params const& ps )
{
  for( uint32_t i=0; i < supports.size(); ++i )
  {
    if( supports[i].size() == 3 )
      ntk.create_maj( supports[i][0], supports[i][1], supports[i][2] );
  }
}

template<class Ntk>
void orthogonal_creator( simulation_view<Ntk>& ntk, std::vector<std::vector<signal<Ntk>>> & supports, 
                           detail::creation_params const& ps )
{
  for( uint32_t i=0; i < supports.size(); ++i )
  {
    ntk.create_xor( supports[i][0], supports[i][1] );
    ntk.create_and( supports[i][0], supports[i][1] );
    ntk.create_and( !supports[i][0], supports[i][1] );
    ntk.create_and( supports[i][0], !supports[i][1] );
    ntk.create_and( !supports[i][0], !supports[i][1] );
  }
}

template<class Ntk>
void chatterjee1( simulation_view<Ntk>& ntk, std::vector<std::vector<signal<Ntk>>> & divisors, 
                           detail::creation_params const& ps )
{
  kitty::partial_truth_table * Y = &ntk.targets[ps.output];
  std::vector<std::pair<std::vector<signal<Ntk>>,std::string>> candidates_s;
  std::vector<std::pair< std::vector<signal<Ntk>>,kitty::dynamic_truth_table>> candidates;

  std::vector<signal<Ntk>> new_signals;
  for( uint32_t i=0; i < divisors.size(); ++i )
  {
    if( divisors[i].size()>1 )
    {
      std::vector<kitty::partial_truth_table*> X;
      for( auto d : divisors[i] )
        X.push_back( &( ntk.sim_patterns[ntk.nodes_to_patterns[ntk.get_node(d)]].pat ) );

      chj_result F = chatterjee_method( X, Y, ntk.seed );

      auto f = std::make_pair( divisors[i], F.tt );
      if( ntk.available_nodes.find( f ) == ntk.available_nodes.end() )
      {
        auto candidate = std::make_pair( divisors[i], F.dtt );
        auto candidate_s = std::make_pair( divisors[i], F.tt );
        candidates.push_back( candidate );
        candidates_s.push_back( candidate_s );
      }
    }
  }
  srand(ntk.seed++);
  if( candidates.size() > ps.max_nodes_total )
  {
    size_t num_candidates = candidates.size() - ps.max_nodes_total;
    for( uint32_t j{0}; j<num_candidates; ++j )
    {
      auto y = rand()%(candidates.size()+1);
      candidates.erase( candidates.begin() + y );
      candidates_s.erase( candidates_s.begin() + y );
    }
  }
  for( uint32_t j=0; j<candidates.size(); ++j )
  {
    auto fnew = ntk.create_node( candidates[j].first, candidates[j].second );
    new_signals.push_back( fnew );
    ntk.available_nodes.insert( candidates_s[j] );
    candidates_s[j].second = candidates_s[j].second + " -> " + std::to_string(fnew) + "~" + std::to_string(ps.output); 
    if( ps.verbose )
      std::cout << candidates_s[j].second << std::endl;
  }
  
}


template<class Ntk>
void random( simulation_view<Ntk>& ntk, std::vector<std::vector<signal<Ntk>>> & divisors, 
                           detail::creation_params const& ps )
{

  std::vector<signal<Ntk>> new_signals;
  for( uint32_t i=0; i < divisors.size(); ++i )
  {
    if( divisors[i].size()>1 )
    {
      uint32_t nvars = divisors[i].size();
      kitty::dynamic_truth_table ttn(nvars);
      kitty::create_random( ttn, ntk.seed++ );
      //kitty::print_binary(ttn);
      //std::cout << std::endl;
      ntk.create_node( divisors[i], ttn );
    }

  }

  
}

} // namespace detail

/*! \brief create nodes from a given support
 *
 */
template<class Ntk>
void create_nodes( simulation_view<Ntk>& ntk, std::vector<std::vector<signal<Ntk>>> & supports, detail::creation_method const& creation_m, detail::creation_params const& creation_ps )
{ 
    switch(creation_m) {
      case detail::creation_method::fgenerator1:
        detail::fgenerator1( ntk, supports, creation_ps );
        break;
      case detail::creation_method::ifgenerator1:
        detail::ifgenerator1( ntk, supports, creation_ps );
        break;
      case detail::creation_method::andgen:
        detail::andgen( ntk, supports, creation_ps );
        break;
      case detail::creation_method::xorgen:
        detail::xorgen( ntk, supports, creation_ps );
        break;
      case detail::creation_method::majgen:
        detail::majgen( ntk, supports, creation_ps );
        break;
      case detail::creation_method::orthogonal_creator:
        detail::orthogonal_creator( ntk, supports, creation_ps );
        break;
      case detail::creation_method::chatterjee1:
        detail::chatterjee1( ntk, supports, creation_ps );
        break;
      case detail::creation_method::random:
        detail::random( ntk, supports, creation_ps );
        break;
    }
      
}
} /* namespace hyperdimensional computing */
} // namespace mockturtle
