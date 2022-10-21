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
  \file simulation_view.hpp
  \brief Implements depth, level, simulation and fanin-size for a network

  \author Andrea Costamagna
*/

#pragma once

#include "../networks/events.hpp"
#include "../traits.hpp"
#include "../utils/cost_functions.hpp"
#include "../utils/node_map.hpp"

#include <kitty/print.hpp>

#include <cstdint>
#include <vector>
#include <set>

namespace mockturtle
{

struct simulation_view_params
{
  /*! \brief Take complemented edges into account for depth computation. */
  bool count_complements{ false };

  /*! \brief Whether PIs have costs. */
  bool pi_cost{ false };
};

/*! \brief Implements `depth`, `level`, `simulation` and `fanin-size` for networks.
 *
 * 
 * This view originates from `depth_view` and is designed for simulation-guided synthesis.
 * Each node maps to two objects  
 * - simulation: A truth table whose *i*-th entry is the value of the node function at the *i*-th simulation pattern.
 * - faninsize: An integer measuring the number of gates in the transitive fanin of the node.
 * 
 * This view implements the network interface methods `level`, `depth`, `simulation`, and `faninsize`. 
 * All methods are computed at construction and can be recomputed by calling the `update_<method>.
 * It also automatically updates network features when creating nodes or
 * creating a PI, or creating a PO on a simulation_view, however, it does not update the information,
 * when modifying or deleting nodes, neither will the critical paths be
 * recalculated (due to efficiency reasons). In order to recalculate levels,
 * depth, and critical paths, one can call `update_levels` instead.
 * To update the simulation patterns call `update_simulations` and to update the faninsizes call `update_faninsizes`.
 * To change the input simulation patterns use the `set_input_simulations` method.
 * `simulation_view` supports any network type and any truth table for the simulation.  
 *
 * **Required network functions:**
 * - `size`
 * - `get_node`
 * - `visited`
 * - `set_visited`
 * - `foreach_fanin`
 * - `foreach_po`
 *
 * Example
 *
   \verbatim embed:rst

   .. code-block:: c++

      // define simulation type and the network type
      using TT = kitty::static_truth_table<2u>;
      using Ntk = aig_network;

      // create network
      Ntk aig = ...;
      
      // create a simulation view on the network
      simulation_view<Ntk,TT> saig{aig};

   \endverbatim
 */
template<class Ntk, class TT, class NodeCostFn = unit_cost<Ntk>, bool has_depth_interface = has_depth_v<Ntk>&& has_level_v<Ntk>&& has_update_levels_v<Ntk> && has_simulation_v<Ntk>&& has_update_simulations_v<Ntk>&& has_faninsize_v<Ntk>&& has_update_faninsizes_v<Ntk>>
class simulation_view
{
};

template<class Ntk, class TT, class NodeCostFn>
class simulation_view<Ntk, TT, NodeCostFn, true> : public Ntk
{
public:
  simulation_view( Ntk const& ntk, simulation_view_params const& ps = {} ) : Ntk( ntk )
  {
    (void)ps;
  }
};

template<class Ntk, class TT, class NodeCostFn>
class simulation_view<Ntk, TT, NodeCostFn, false> : public Ntk
{
public:
  using storage = typename Ntk::storage;
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

  explicit simulation_view( NodeCostFn const& cost_fn = {}, simulation_view_params const& ps = {} )
      : Ntk(), _ps( ps ), _levels( *this ), _crit_path( *this ), _faninsizes( *this ), _simulations( *this ), _input_simulations( *this ), _tfis( *this ), _cost_fn( cost_fn )
  {
    static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
    static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
    static_assert( has_compute_v<Ntk,TT>, "Ntk does not implement the compute method" );
    static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
    static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement the is_complemented method" );
    static_assert( has_visited_v<Ntk>, "Ntk does not implement the visited method" );
    static_assert( has_set_visited_v<Ntk>, "Ntk does not implement the set_visited method" );
    static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
    static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin method" );

    add_event = Ntk::events().register_add_event( [this]( auto const& n ) { on_add( n ); } );
  }

  /*! \brief Standard constructor requiring explicit `simulation_type` declaration
   *
   * \param ntk Base network
   */
  explicit simulation_view( Ntk const& ntk, NodeCostFn const& cost_fn = {}, simulation_view_params const& ps = {} )
      : Ntk( ntk ), _ps( ps ), _levels( ntk ), _crit_path( ntk ), _faninsizes( *this ), _simulations( ntk ), _input_simulations( ntk ), _tfis( ntk ), _cost_fn( cost_fn )
  {
    static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
    static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
    static_assert( has_compute_v<Ntk,TT>, "Ntk does not implement the compute method" );
    static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
    static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement the is_complemented method" );
    static_assert( has_visited_v<Ntk>, "Ntk does not implement the visited method" );
    static_assert( has_set_visited_v<Ntk>, "Ntk does not implement the set_visited method" );
    static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
    static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin method" );

    update_levels();
    update_faninsizes();
    
    add_event = Ntk::events().register_add_event( [this]( auto const& n ) { on_add( n ); } );
  }

  /*! \brief Standard constructor.
   *
   * \param ntk Base network
   * \param input_simulations vector containing one simulation per each input
   */
  explicit simulation_view( Ntk const& ntk, std::vector<TT> input_simulations, NodeCostFn const& cost_fn = {}, simulation_view_params const& ps = {} )
      : Ntk( ntk ), _ps( ps ), _levels( ntk ), _crit_path( ntk ), _faninsizes( *this ), _simulations( ntk ), _input_simulations( ntk ), _tfis( ntk ), _cost_fn( cost_fn )
  {
    static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
    static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
    static_assert( has_compute_v<Ntk,TT>, "Ntk does not implement the compute method" );
    static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
    static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement the is_complemented method" );
    static_assert( has_visited_v<Ntk>, "Ntk does not implement the visited method" );
    static_assert( has_set_visited_v<Ntk>, "Ntk does not implement the set_visited method" );
    static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
    static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin method" );

    update_levels();
    update_faninsizes();
    set_input_simulations( input_simulations );
    update_simulations();
    
    add_event = Ntk::events().register_add_event( [this]( auto const& n ) { on_add( n ); } );
  }

  /*! \brief Copy constructor. */
  explicit simulation_view( simulation_view<Ntk, TT, NodeCostFn, false> const& other )
      : Ntk( other ), _ps( other._ps ), _levels( other._levels ), _crit_path( other._crit_path ), _faninsizes( *this ), _depth( other._depth ), _simulations( other._simulations ), _input_simulations( other._input_simulations ), _tfis( other._tfis ), _cost_fn( other._cost_fn )
  {
    add_event = Ntk::events().register_add_event( [this]( auto const& n ) { on_add( n ); } );
  }

  simulation_view<Ntk, TT, NodeCostFn, false>& operator=( simulation_view<Ntk, TT, NodeCostFn, false> const& other )
  {
    /* delete the event of this network */
    Ntk::events().release_add_event( add_event );

    /* update the base class */
    this->_storage = other._storage;
    this->_events = other._events;

    /* copy */
    _ps = other._ps;
    _levels = other._levels;
    _faninsizes = other._faninsizes;
    _faninsizes = other._tfis;
    _simulations = other._simulations;
    _input_simulations = other._input_simulations;
    _crit_path = other._crit_path;
    _depth = other._depth;
    _cost_fn = other._cost_fn;

    /* register new event in the other network */
    add_event = Ntk::events().register_add_event( [this]( auto const& n ) { on_add( n ); } );

    return *this;
  }

  ~simulation_view()
  {
    Ntk::events().release_add_event( add_event );
  }

  void set_input_simulations( std::vector<TT> input_simulations )
  {
    const0 = input_simulations[0].construct();
    _simulations[get_constant(false)]=const0;
    _input_simulations[get_constant(false)]=const0;

    uint32_t i = 0;
    this->foreach_pi( [&]( auto const& n ) {
      _input_simulations[n] = input_simulations[i];
      _simulations[n] = input_simulations[i++];
    } );
  }

  uint32_t depth() const
  {
    return _depth;
  }

  uint32_t level( node const& n ) const
  {
    return _levels[n];
  }

  uint32_t faninsize( node const& n ) const
  {
    return _faninsizes[n];
  }

  TT simulation( node const& n ) const
  {
    return _simulations[n];
  }

  bool is_on_critical_path( node const& n ) const
  {
    return _crit_path[n];
  }

  void set_level( node const& n, uint32_t level )
  {
    _levels[n] = level;
  }

  void set_faninsize( node const& n, uint32_t faninsize )
  {
    _faninsizes[n] = faninsize;
  }

  void set_simulation( node const& n, TT simulation )
  {
    _simulations[n] = simulation;
  }

  void set_depth( uint32_t level )
  {
    _depth = level;
  }

  void update_levels()
  {
    _levels.reset( 0 );
    _crit_path.reset( false );

    this->incr_trav_id();
    compute_levels();
  }

  void update_faninsizes()
  {
    _faninsizes.reset( 0 );

    this->incr_trav_id();
    compute_faninsizes();
  }

  void resize_levels()
  {
    _levels.resize();
  }

  void resize_faninsizes()
  {
    _faninsizes.resize();
  }

  void update_simulations()
  {
    _simulations.reset();

    this->incr_trav_id();
    compute_simulations();
  }

  void create_po( signal const& f )
  {
    Ntk::create_po( f );
    _depth = std::max( _depth, _levels[f] );
  }

  signal create_pi( TT& input_simulation )
  {
    signal pi_sig = Ntk::create_pi();
    
    _levels.resize();
    _levels[ get_node( pi_sig ) ] = 0;

    _faninsizes.resize();
    _faninsizes[ get_node( pi_sig ) ] = 0;

    _input_simulations.resize();
    _input_simulations[ get_node( pi_sig ) ] = input_simulation;

    _simulations.resize();
    _simulations[ get_node( pi_sig ) ] = input_simulation;

    return pi_sig;
  }

  signal create_pi()
  {
    signal pi_sig = Ntk::create_pi();

    _levels.resize();
    _levels[ get_node( pi_sig ) ] = 0;

    _faninsizes.resize();
    _faninsizes[ get_node( pi_sig ) ] = 0;

    _simulations.resize();
    _input_simulations.resize();

    return pi_sig;
  }

  node get_node( signal f )
  {
    return Ntk::get_node( f );
  }

  signal get_constant( bool value = false ) const
  {
    return Ntk::get_constant( value );
  }

private:
  uint32_t compute_levels( node const& n )
  {
    if ( this->visited( n ) == this->trav_id() )
    {
      return _levels[n];
    }
    this->set_visited( n, this->trav_id() );

    if ( this->is_constant( n ) )
    {
      return _levels[n] = 0;
    }
    if ( this->is_pi( n ) )
    {
      assert( !_ps.pi_cost || _cost_fn( *this, n ) >= 1 );
      return _levels[n] = _ps.pi_cost ? _cost_fn( *this, n ) - 1 : 0;
    }

    uint32_t level{ 0 };
    this->foreach_fanin( n, [&]( auto const& f ) {
      auto clevel = compute_levels( this->get_node( f ) );
      if ( _ps.count_complements && this->is_complemented( f ) )
      {
        clevel++;
      }
      level = std::max( level, clevel );
    } );

    _levels[n] = level + _cost_fn( *this, n );
    return _levels[n];
  }

  void compute_levels()
  {
    _depth = 0;
    this->foreach_po( [&]( auto const& f ) {
      auto clevel = compute_levels( this->get_node( f ) );
      if ( _ps.count_complements && this->is_complemented( f ) )
      {
        clevel++;
      }
      _depth = std::max( _depth, clevel );
    } );

    this->foreach_po( [&]( auto const& f ) {
      const auto n = this->get_node( f );
      if ( _levels[n] == _depth )
      {
        set_critical_path( n );
      }
    } );
  }

  void set_critical_path( node const& n )
  {
    _crit_path[n] = true;
    if ( !this->is_constant( n ) && !( _ps.pi_cost && this->is_pi( n ) ) )
    {
      const auto lvl = _levels[n];
      this->foreach_fanin( n, [&]( auto const& f ) {
        const auto cn = this->get_node( f );
        auto offset = _cost_fn( *this, n );
        if ( _ps.count_complements && this->is_complemented( f ) )
        {
          offset++;
        }
        if ( _levels[cn] + offset == lvl && !_crit_path[cn] )
        {
          set_critical_path( cn );
        }
      } );
    }
  }

  TT compute_simulations( node const& n )
  {
    if ( this->visited( n ) == this->trav_id() )
    {
      return _simulations[n];
    }
    if ( this->is_constant( n ) )
    {
      _simulations[n] = const0;
      return _simulations[n];
    }
    if ( this->is_pi( n ) )
    {
      _simulations[n] = _input_simulations[n];
      return _input_simulations[n];
    }

    TT simulation;
    std::vector<TT> children_simulations;
    this->foreach_fanin( n, [&]( auto const& f ) {
      children_simulations.push_back( compute_simulations( get_node(f) ) );
    } );
    
    _simulations[n] = Ntk::compute( n, children_simulations.begin(), children_simulations.end() );

    this->set_visited( n, this->trav_id() );

    return _simulations[n];
  }

  void compute_simulations()
  {
    _simulations[get_constant(false)]=_input_simulations[get_constant(false)];
    
    this->foreach_pi( [&]( auto const& n ) {
      _simulations[n] = _input_simulations[n];

      this->set_visited( n, this->trav_id() );
    } );

    this->foreach_po( [&]( auto const& f ) {
      compute_simulations( this->get_node( f ) );
    } );
  }

  uint32_t compute_faninsizes( node const& n )
  {
    if ( this->visited( n ) == this->trav_id() )
    {
      return _faninsizes[n];
    }
    if ( this->is_constant( n ) )
    {
      _faninsizes[n] = 0;
      return _faninsizes[n];
    }
    if ( this->is_pi( n ) )
    {
      return 0;
    }

    std::set<node> tfi;
    std::set<node> tfi_child;

    this->foreach_fanin( n, [&]( auto const& f ) {
      compute_faninsizes( get_node(f) );
      tfi_child = _tfis[ get_node(f) ];
      tfi.merge( tfi_child );
      tfi.merge( std::set<node>{get_node(f)} );
    } );
    
    this->set_visited( n, this->trav_id() );
    _tfis[n] = tfi;
    _faninsizes[n] = tfi.size();

    return _faninsizes[n];
  }

  void compute_faninsizes()
  {
    this->foreach_pi( [&]( auto const& n ) {
      _faninsizes[n] = 0;
      this->set_visited( n, this->trav_id() );
    } );

    this->foreach_po( [&]( auto const& f ) {
      compute_faninsizes( this->get_node( f ) );
    } );
  }

  void on_add( node const& n )
  {
    _levels.resize();

    uint32_t level{ 0 };
    this->foreach_fanin( n, [&]( auto const& f ) {
      auto clevel = _levels[f];
      if ( _ps.count_complements && this->is_complemented( f ) )
      {
        clevel++;
      }
      level = std::max( level, clevel );
    } );

    _levels[n] = level + _cost_fn( *this, n );

    _simulations.resize();

    TT simulation;
    std::vector<TT> children_simulations;
    this->foreach_fanin( n, [&]( auto const& f ) {
      children_simulations.push_back( _simulations[f] );
    } );
    
    _simulations[n] = Ntk::compute( n, children_simulations.begin(), children_simulations.end() );

    _tfis.resize();
    _faninsizes.resize();

    std::set<node> tfi;
    std::set<node> tfi_child;    
    
    this->foreach_fanin( n, [&]( auto const& f ) {
      tfi_child = _tfis[ get_node(f) ];
      tfi.merge( tfi_child );
      tfi.merge( std::set<node>{get_node(f)} );
    } );
    
    _tfis[n] = tfi;
    _faninsizes[n] = tfi.size();

  }

  simulation_view_params _ps;
  node_map<uint32_t, Ntk> _levels;
  node_map<uint32_t, Ntk> _crit_path;
  node_map<uint32_t, Ntk> _faninsizes;
  node_map<TT, Ntk> _simulations;
  node_map<TT, Ntk> _input_simulations;
  node_map<std::set<node>, Ntk> _tfis;
  uint32_t _depth{};
  NodeCostFn _cost_fn;
  TT const0;

  std::shared_ptr<typename network_events<Ntk>::add_event_type> add_event;
};

template<class T, class TT>
simulation_view( T const&, std::vector<TT> ) -> simulation_view<T, TT>;

template<class T, class TT, class NodeCostFn = unit_cost<T>>
simulation_view( T const&, std::vector<TT>, NodeCostFn const&, simulation_view_params const& ) -> simulation_view<T, TT, NodeCostFn>;

} // namespace mockturtle