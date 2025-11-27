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
  \file mapped_list.hpp
  \brief Index list for mapped networks

  \author Andrea Costamagna
*/

#pragma once

#include "../networks/bound.hpp"

#include <vector>

namespace mockturtle
{
  /*! \brief mapped index list
   * Area is always maintained up-to-date.
   */
  class mapped_list
  {
public:
    /* space reserved for constants */
    static constexpr uint32_t offset = 2;
    using element_type = uint32_t;
    static constexpr element_type const0 = 0;
    static constexpr element_type const1 = 1;

private:
    /*! \brief node of a mapped index list */
    struct node
    {

      node(){}

      ~node(){}

      node( const node& other )
      {
        id = other.id;
        for( int i = 0; i < other.children.size(); ++i )
        {
          children.push_back( other.children[i] );
        }
      }

      node( node&& other ) noexcept
      {
        id = other.id;
        for( int i = 0; i < other.children.size(); ++i )
        {
          children.push_back( other.children[i] );
        }
      }

      node& operator=( const node& other )
      {
        if ( this != &other )
        {
          children = other.children;
          id = other.id;
        }
        return *this;
      }

      node& operator=( node&& other ) noexcept
      {
        if (this != &other)
        {
          children = other.children;
          id = other.id;
        }
        return *this;
      }

      bool operator==( node const& other )
      {
        bool is_same = true;
        is_same &= ( children == other.children );
        is_same &= ( id == other.id );
        return is_same;
      }

      /* preferred construction from children and gate identifier */
      node( std::vector<element_type>const& children, unsigned int id )
        :  children( children ), id( id )
      {}

      /* cell identifier */
      unsigned int id;
      /* indices of the fanins */
      std::vector<element_type> children;
    };

public:
    mapped_list()
      : nodes(),
        outputs(),
        pin_delays(),
        num_vars( 0 ), 
        area( 0 ),
        power( 0 )
    {
    }

    mapped_list( uint32_t num_vars )
      : nodes(),
        outputs(),
        pin_delays( num_vars ),
        num_vars( num_vars ), 
        area( 0 ),
        power( 0 )
    {
      nodes.reserve( 10 );
    }

    mapped_list( uint32_t num_vars, size_t size )
      : nodes(),
        outputs(),
        pin_delays(),
        num_vars( num_vars ), 
        area( 0 ),
        power( 0 )
    {
      nodes.reserve( size );
    }

    mapped_list( mapped_list const& other )
      : outputs( other.outputs ),
        num_vars( other.num_vars ),
        area( other.area ),
        power( other.power )
    {

      for( int i = 0; i < other.nodes.size(); ++i )
      {
        nodes.push_back( other.nodes[i] );
      }
      for( int i = 0; i < other.pin_delays.size(); ++i )
      {
        pin_delays.push_back( other.pin_delays[i] );
      }
    }

    ~mapped_list(){}

    mapped_list( mapped_list&& other ) noexcept
      : outputs( other.outputs ),
        nodes( other.nodes ),
        num_vars( other.num_vars ),
        pin_delays( other.pin_delays ),
        area( other.area ),
        power( other.power )
    {
    }

      mapped_list& operator=( const mapped_list& other ) {
          if (this != &other) {
            nodes = other.nodes;
            outputs = other.outputs;
            pin_delays = other.pin_delays;
            num_vars = other.num_vars;
            area = other.area;
            power = other.power;
          }
          return *this;
      }

      bool operator==( const mapped_list& other )
      {
        bool is_same = true;
        if ( nodes.size() != other.nodes.size() )
          return false;
          for ( int i = 0; i < nodes.size(); ++i )
            is_same &= ( nodes[i] == other.nodes[i] );
        is_same &= ( outputs == other.outputs );
        is_same &= ( pin_delays == other.pin_delays );
        return is_same;
      }

      mapped_list& operator=( mapped_list&& other ) noexcept {
          if (this != &other) {
            nodes = other.nodes;
            outputs = other.outputs;
            pin_delays = other.pin_delays;
            num_vars = other.num_vars;
            area = other.area;
            power = other.power;
          }
          return *this;
      }

//    template<class Ntk, class Bindings_t>
//    mapped_list( Ntk & ntk, Bindings_t const& lib, std::vector<typename Ntk::node> const& cut, typename Ntk::signal f )
//      : nodes(),
//        outputs(),
//        pin_delays(),
//        num_vars( 0 ), 
//        area( 0 ),
//        power( 0 )
//    {
//      /* mark the cut nodes as visited */
//      ntk.incr_trav_id();
//      set_pis( cut );
////    for( auto const& n : cut )
////    {
////      ntk.set_visited( n, ntk.trav_id() );
////      auto fi = create_pi();
////      ntk.set_value( n, fi );
////    }
//      /* recursively construct */
//      construct_rec( ntk, f );
//      create_po( size()-1 );
//      /* compute the longest path of each input to output */
//      compute_pin_delays( lib );
//    }

    template<class Ntk, class Bindings_t, auto L>
    mapped_list( Ntk & ntk, Bindings_t const& lib, std::array<typename Ntk::node, L> const& cut, typename Ntk::signal f )
      : outputs(),
        pin_delays(),
        num_vars( 0 ), 
        area( 0 ),
        power( 0 )
    {
      nodes.reserve( 100 );
      /* mark the cut nodes as visited */
      ntk.incr_trav_id();
      /* The PIs are assigned in topological order */
      for( auto const& n : cut )
      {
        if( n < ntk.size() )
        {
          ntk.set_visited( n, ntk.trav_id() );
          ntk.set_value( n, -1 );
        }
      }
//      assign_pis( ntk, f );
      for( int i = 0; i < L; ++i ) // alternative to assign_pis
      {
        if( cut[i] < ntk.size() )
        {
          auto fi = create_pi();
          ntk.set_value( cut[i], fi );
        }
      }


      while( num_pis() < L )
      {
        create_pi();
      }
      

      ntk.incr_trav_id();
      for( auto const& n : cut )
      {
        if( n < ntk.size() )
        {
          ntk.set_visited( n, ntk.trav_id() );
        }
      }
//      for( auto const& n : cut )
//      {
//        auto fi = create_pi();
//        if( n < ntk.size() )
//        {
//          ntk.set_visited( n, ntk.trav_id() );
//          ntk.set_value( n, fi );
//        }
//      }
      /* recursively construct */
      construct_rec( ntk, f );
      create_po( ntk.value( ntk.get_node( f ) ) );
      /* compute the longest path of each input to output */
      compute_pin_delays( lib );
    }


    template<class Ntk, class Bindings_t>
    mapped_list( Ntk & ntk, Bindings_t const& lib, std::vector<typename Ntk::node> const& cut, typename Ntk::signal f )
      : outputs(),
        pin_delays(),
        num_vars( 0 ), 
        area( 0 ),
        power( 0 )
    {
      uint32_t L = cut.size();
      nodes.reserve( 100 );
      /* mark the cut nodes as visited */
      ntk.incr_trav_id();
      /* The PIs are assigned in topological order */
      for( auto const& n : cut )
      {
        if( n < ntk.size() )
        {
          ntk.set_visited( n, ntk.trav_id() );
          ntk.set_value( n, -1 );
        }
      }
//      assign_pis( ntk, f );
      for( int i = 0; i < L; ++i ) // alternative to assign_pis
      {
        auto fi = create_pi();
        ntk.set_value( cut[i], fi );
      }


      while( num_pis() < L )
      {
        create_pi();
      }
      

      ntk.incr_trav_id();
      for( auto const& n : cut )
      {
        if( n < ntk.size() )
        {
          ntk.set_visited( n, ntk.trav_id() );
        }
      }
//      for( auto const& n : cut )
//      {
//        auto fi = create_pi();
//        if( n < ntk.size() )
//        {
//          ntk.set_visited( n, ntk.trav_id() );
//          ntk.set_value( n, fi );
//        }
//      }
      /* recursively construct */
      construct_rec( ntk, f );
      create_po( ntk.value( ntk.get_node( f ) ) );
      /* compute the longest path of each input to output */
      compute_pin_delays( lib );
    }

private:

    template<class Ntk>
    double get_area( Ntk const& ntk, typename Ntk::node n )
    {
      if constexpr ( std::is_same_v<typename Ntk::base_type, bound_network> )
      {
        return ntk.get_area( n );
      }
      else
      {
        auto const cell = ntk.get_cell( n );
        return cell.gates[0].area;
      }
    }
    
    template<class Ntk>
    unsigned int get_binding_id( Ntk const& ntk, typename Ntk::node n )
    {
      if constexpr ( std::is_same_v<typename Ntk::base_type, bound_network> )
      {
        return ntk.get_binding_id( n );
      }
      else
      {
        auto const cell = ntk.get_cell( n );
        return cell.gates[0].id;
      }
    }

    /*! \brief recursively construct the list 
     * This function assumes that some nodes in the TFI of f
     * are marked as visited ( visited flag = trav_id ) and
     * their value field is set to the element type in the index list
    */
    template<class Ntk>
    void construct_rec( Ntk & ntk, typename Ntk::signal f )
    {
      auto const n = ntk.get_node( f );
      /* termination conditions */
      if( ntk.visited( n ) == ntk.trav_id() )
      {
        return;
      }
      else if( ntk.is_pi( n ) ) // not marked as visited: illegal cone
      {
        std::cerr << "[e] list construction from ntk failed" << std::endl;
        return;
      }
      /* node construction */
      std::vector<element_type> children;
      
      area += get_area( ntk, n );
      auto const id = get_binding_id( ntk, n );
      bool is_valid = true;
      ntk.foreach_fanin( n, [&]( auto fi, auto ii ){
        auto const ni = ntk.get_node( fi );
        construct_rec( ntk, fi );
        children.emplace_back( ntk.value( ni ) );
      } );
      auto v = create_node( children, id );
      ntk.set_value( n, v );
      ntk.set_visited( n, ntk.trav_id() );
    }

    template<class Ntk>
    void assign_pis( Ntk & ntk, typename Ntk::signal f )
    {
      auto const n = ntk.get_node( f );
      /* termination conditions */
      if( ntk.visited( n ) == ntk.trav_id() )
      {
        if( ntk.value( n ) == -1 )
        {
          auto fi = create_pi();
          ntk.set_value( n, fi );
        }
        return;
      }
      
      ntk.foreach_fanin( n, [&]( auto fi, auto ii ){
        auto const ni = ntk.get_node( fi );
        assign_pis( ntk, fi );
      } );
    }

public:
    template<class Bindings_t>
    void compute_pin_delays( Bindings_t const& lib )
    {
      std::vector<double> delays( offset + num_pis() + num_gates(), 0 );
      uint32_t index = offset + num_pis() + num_gates() - 1;
      foreach_gate_rev( [&]( auto const& children, unsigned int id ){
        for( int i = 0; i < children.size(); ++i )
        {
          double const fall = lib[id].pins[i].fall_block_delay;
          double const rise = lib[id].pins[i].rise_block_delay;
          double const delta = std::max( rise, fall );
          delays[children[i]] = std::max( delays[children[i]], delays[index] + delta );
        }
        index--;
      } );
      pin_delays.clear();
      pin_delays.resize( num_vars, 0 );
      int i = 0;
      
      for( int i = 0; i < num_vars; ++i )
      {
        element_type lit = i + offset;
        pin_delays[i] = delays[lit];
      }
    }

  public:
    element_type create_pi()
    {
      if( nodes.size() > 0 )
      {
        std::cerr << "[e] illegal PI construction in the presence of nodes" << std::endl;
        return std::numeric_limits<element_type>::max();
      }
      const element_type f{ num_vars + offset };
      num_vars++;
      return f;
    }

    element_type pi_at( uint32_t index ) const
    {
      const element_type f{ index + offset };
      return f;
    }

    element_type po_at( uint32_t index ) const
    {
      assert( index < outputs.size() );
      return outputs[index];
    }

    void create_po( element_type const v )
    {
      outputs.push_back( v );
    }

    /*! \brief create a node in the list 
     * the element_type is an unique identifier to reference an element of the list
    */
    element_type create_node( std::vector<element_type> const& children, unsigned int id )
    {
      const element_type f = nodes.size() + num_vars + offset;
      nodes.emplace_back( children, id );
      return f;
    }

    bool is_constant( element_type f ) const
    {
      return f < 2;
    }

    uint32_t get_index( uint32_t lit )
    {
      return lit - offset;
    }

    uint32_t get_lit( uint32_t index ) const
    {
      return index + offset;
    }


    element_type get_constant( bool value )
    {
      return value ? 1 : 0;
    }

    bool is_input( element_type f ) const
    {
      return ( ( f >= 2 ) && f < ( 2 + num_vars ) );
    }

    uint32_t node_to_index( element_type f ) const
    {
      return f - num_vars - offset;
    }

    element_type index_to_node( auto index ) const
    {
      return index + num_vars + offset;
    }

    uint32_t num_gates() const
    {
      return nodes.size();
    }

    uint32_t size() const
    {
      return nodes.size() + num_pis() + offset;
    }

    uint32_t num_pis() const
    {
      return num_vars;
    }

    uint32_t num_pos() const
    {
      return outputs.size();
    }

    /* perm are the input permutation received by the P-canonization */
    void permute_inputs( std::vector<uint8_t> const& perm )
    {
      std::vector<uint32_t> permutation( perm.size() );
      for( int i=0; i < perm.size(); ++i )
      {
        permutation[perm[i]] = i;
      }

      if( perm.size() != num_vars )
      {
        std::cout << perm.size() << "!=" << num_vars << std::endl;
        std::cerr << "[e] permutation vector does not match the number of variables" << std::endl;
      }
      int inode = 0;
      foreach_gate( [&]( auto const children, auto const id ){
        for( int ichild = 0; ichild < children.size(); ++ichild )
        {
          auto const child = children[ichild];
          if( ( child < ( num_vars + offset ) ) && ( child >= offset )  )
          {
            nodes[inode].children[ichild] = permutation[child-offset] + offset;
          }
        }
        inode++;
      } );

      auto const old_delays = pin_delays;
      for( int i = 0; i < num_vars; ++i )
      {
        pin_delays[permutation[i]] = old_delays[i];
      }
    }


  struct timed_index
  {
    timed_index( uint32_t index, double delay ) : index(index), delay(delay){}
    timed_index() : index(0), delay(0){}
    ~timed_index(){};

    uint32_t index;
    double delay;
  };

    template<class Times, class uint32_t>
    void match_times( Times const& times, std::vector<std::vector<uint32_t>> const& symm )
    {
      if( times.size() != num_vars )
      {
        std::cerr << "[e] times have wrong size" << std::endl;
      }
      /* extract the propagation time */
      std::vector<timed_index> inputs( num_vars );
      for( int i = 0; i < num_vars; ++i )
      {
        timed_index const tindex{ i, times[i] };
        inputs[i] = tindex;
      }
      /* consider all the variables */
      bool swapped;
      for( auto const& vars : symm )
      {
        int const n = vars.size();
        /* bubble-sort variables in a symmetry group */
        for( int v1 = 0; v1 < n - 1; ++v1 )
        {
          swapped = false;
          for( int v2 = n - 1; v2 >= v1 + 1; --v2 )
          {
            uint32_t const var1 = vars[v2-1];
            uint32_t const var2 = vars[v2];
            if( inputs[var1].delay > inputs[var2].delay )
            {
              std::swap( inputs[var1], inputs[var2] );
              swapped = true;
            }
          }
        }
      }
      /* now all the variables in a symmetry group are sorted from slowest to fastest */
      std::vector<uint8_t> perm( num_vars );
      for( int i=0; i<num_vars; ++i )
        perm[i] = inputs[i].index;
      permute_inputs( perm );
    }


    /*! \brief iterator to apply a lambda to all the nodes */
    template<typename Fn>
    void foreach_pi( Fn&& fn ) const
    {
      for ( uint64_t i = 0; i < num_vars; ++i )
      {
        fn( static_cast<element_type>( offset + i ) );
      }
    }

    /*! \brief iterator to apply a lambda to all the nodes */
    template<typename Fn>
    void foreach_gate( Fn&& fn ) const
    {
      for ( uint64_t i = 0; i < nodes.size(); i ++ )
      {
        fn( nodes[i].children, nodes[i].id );
      }
    }

    /*! \brief iterator to apply a lambda to all the outputs */
    template<typename Fn>
    void foreach_po( Fn&& fn ) const
    {
      for ( auto o : outputs )
      {
        fn( o );
      }
    }

    /*! \brief iterator to apply a lambda to all the nodes */
    template<typename Fn>
    void foreach_gate_rev( Fn&& fn ) const
    {
      for ( int i = nodes.size()-1; i >= 0 ; --i )
      {
        fn( nodes[i].children, nodes[i].id );
      }
    }

    void print() const
    {
      std::cout << "L: [";
      int i = num_pis() + offset;
      for( auto t : pin_delays )
      {
        std::cout << t << " ";
      }
      std::cout << "] ";
      foreach_gate( [&]( auto children, auto id ){
        std::cout << i++ << "= " <<"{";
        for ( auto const& x : children )
        {
          std::cout << x << " ";
        }
        std::cout << "}." << id << " ";
      } );
      std::cout << " pos: ";
      foreach_po( [&]( element_type const& v ){
        std::cout << v << " ";
      } );
      std::cout << std::endl;
    }

    template<class Ntk, class Leaves_t>
    typename Ntk::signal insert( Ntk& ntk, Leaves_t const& leaves ) const
    {
      using signal = typename Ntk::signal;
      using node = typename Ntk::node;

      //if( leaves.size() != num_pis() )
      //{
      //  std::cerr << "[e] cut-size inconsistent with the list's PIs" << std::endl;
      //}
      std::vector<signal> sigs( num_pis() + num_gates() + offset );
      sigs[0] = ntk.get_constant( false );
      sigs[1] = ntk.get_constant( true );
      for( int i = 0; i < num_pis(); ++i )
      {
        sigs[offset+i] = ntk.make_signal( leaves[i] );
      }
      int i = num_pis() + offset;
      foreach_gate( [&]( auto children, auto id ){
        std::vector<signal> fanin;
        for( auto c : children )
        {
          fanin.push_back( sigs[c] );
        }
        auto const g = ntk.get_binding_from_id( id );
        sigs[i] = ntk.create_node( fanin, g.name )[0];
        assert( (children[0] != po_at(0)) );
        i++;
      } );

      return sigs[po_at(0)];
    }

    template<class Bindings_t>
    double get_area( Bindings_t const& binds ) const
    {
      double area = 0;
      foreach_gate( [&]( auto children, auto id ){
        area += binds[id].area;
      } );
      return area;
    }

    double get_area() const
    {
      /* check that area is consistent with num of gates */
      return area;
    }


    template<class Bindings_t>
    double get_power( Bindings_t const& binds ) const
    {
      double power = 0;
      foreach_gate( [&]( auto children, auto id ){
        auto const tt = binds[id].function;
        double const n0 = kitty::count_zeros( tt );
        double const nb = tt.num_bits();
        double const n1 = nb - n0;
        power += 2*n1*n0/( nb*nb );
      } );
      return power;
    }

    std::vector<node> get_nodes() const
    {
      return nodes;
    }

    std::vector<element_type> get_outputs() const
    {
      return outputs;
    }

    element_type get_num_vars() const
    {
      return num_vars;
    }

    std::vector<double> get_pin_delays() const
    {
      return pin_delays;
    }

    double get_power() const
    {
      return power;
    }

    void set_power( double cost )
    {
      power = cost;
    }

    void set_area( double value )
    {
      area = value;
    }

    /* store the nodes  */
    std::vector<node> nodes;
    std::vector<element_type> outputs;
    std::vector<double> pin_delays;
    element_type num_vars;
    double area;
    double power;
  };


  template<class Ntk, class Leaves>
  typename Ntk::signal insert( Ntk & ntk, mapped_list const& list, Leaves const& leaves )
  {
    auto const& lib = ntk.get_bindings();
    std::vector<typename Ntk::signal> fs;
    fs.push_back( ntk.get_constant( false ) );
    fs.push_back( ntk.get_constant( true ) );
    int i = 0;
    list.foreach_pi( [&]( auto lit ){
      fs.push_back( leaves[i++] );
    } );
    std::vector<typename Ntk::signal> fanins;
    typename Ntk::signal fo;
    list.foreach_gate( [&]( auto children, auto id ){
      fanins.clear();
      for( auto child : children )
      {
        fanins.push_back( fs[child] );
      }

      fo = ntk.create_node( fanins, lib[id].name )[0];
      fs.push_back( fo );
    } );
    return fo;
  }
} // namespace mockturtle