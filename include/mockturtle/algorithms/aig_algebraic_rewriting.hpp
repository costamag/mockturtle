/*!
  \file aig_algebraic_rewriting.hpp
  \brief AIG algebraric rewriting

  EPFL CS-472 2021 Final Project Option 1
*/

#pragma once

#include "../networks/aig.hpp"
#include "../views/depth_view.hpp"
#include "../views/topo_view.hpp"

namespace mockturtle
{

namespace detail
{

template<class Ntk>
class aig_algebraic_rewriting_impl
{
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  aig_algebraic_rewriting_impl( Ntk& ntk )
    : ntk( ntk )
  {
    static_assert( has_level_v<Ntk>, "Ntk does not implement depth interface." );
  }

  void run()
  {
    bool cont{true}; /* continue trying */
    while ( cont )
    {
      cont = false; /* break the loop if no updates can be made */
      ntk.foreach_gate( [&]( node n ){
        if ( try_algebraic_rules( n ) )
        {
          ntk.update_levels();
          cont = true;
        }
      });
    }
  }

private:
  /* Try various algebraic rules on node n. Return true if the network is updated. */
  bool try_algebraic_rules( node n )
  {
    if ( try_associativity( n ) )
      return true;
    if ( try_distributivity( n ) )
      return true;
    if ( try_3_levels_3_nodes( n ) )
      return true;
    /* TODO: add more rules here... */

    return false;
  }

  /* Try the associativity rule on node n. Return true if the network is updated. */
  bool try_associativity( node n )
  {

    signal s1, s2, s3, s4; /* signals present in the sub-network under analysis */

    if ( order_children( n, s1, s2 ) && ( ntk.level( n ) > 1 ) ) /* if the children are on different levels sL = s1, sH = s2 */
    {

      node n1 = ntk.get_node( s2 );
      bool has_n1_two_levels = order_children( n1, s3, s4 ); /* order sL = s3 and sH = s4. Save if they are on different levels */
      
      signals_comparator cmp14( s1, s4 ); /* to know if s1 = s4 or s1 = s4' */ 

      if ( ntk.is_complemented( s2 ) ) /* phi2(x) = x' : f(s1,s3,s4)=(s1 (s3' + s4')) */
      {
        if( cmp14.abs == true ) /* if |s1| = |s4| */
        { 
          /* swap s3 and s4. So that s3 is always the one compared to s1 */
          swap_elements( s3, s4 );
        }
        signals_comparator cmp13( s1, s3 );
        if( cmp13.abs == true ) /* |s1|=|s3| */
        {
          if( cmp13.sign == true ) /* s1=s3 */
          {
            /* case f(s1,s3,s4)->(s1 s4')  */
            auto f_new = ntk.create_and( s1, !s4 );
            ntk.substitute_node( n, f_new );
            return true;
          }
          else /* s1=s3' */
          {
            /* case f(s1,s3,s4)->(s1)  */
            ntk.substitute_node( n, s1 );
            return true;
          }
          return false;
        }
        return false;
      }
      else /* phi2(x) = x : f(s1,s3,s4)=(s1(s3 s4)) */
      {
        if( cmp14.abs == true ) // put in s3 the one equal to s1 in absolute value
        {
          swap_elements( s3, s4 );
        }

        signals_comparator cmp13( s1, s3 );
        if( cmp13.abs ) /* |s1|= |s3| */
        {
          if( cmp13.sign ) /* s1 = s3 */
          {
            /* case f(s1,s3,s4)->(s1 s4)  */
            auto f_new = ntk.create_and( s4, s1 );
            ntk.substitute_node( n, f_new );
            return true;
          }
          else /* s1 = s3' */
          {
            /* case f(s1,s3,s4)-> 0  */
            auto f_new = ntk.get_constant( 0 );
            ntk.substitute_node( n, f_new );
            return true;
          }
          return false;
        }
        /* if the following piece of code is executed it means that s4 and s3 have not been swapped and s4 is 
         * higher than s3. In this case it is interesting to apply the algebraic rules on the critical path,
         * to reduce the depth. If node n is critical the critical signal must be s4 since its level is the 
         * highest in the network.
         */
        if ( ntk.is_on_critical_path( n ) && has_n1_two_levels && is_higher( s4, s1) )
        {
          /* f -> s4(s1 s3) */
          auto f1 = ntk.create_and( s1, s3 );
          auto f_new = ntk.create_and( s4, f1 );
          ntk.substitute_node( n, f_new );
          return true;
        }
        return false;
      }
      return false;
    }
    return false;
  }

  /* Try the distributivity rule on node n. Return true if the network is updated. */
  bool try_distributivity( node n )
  {
    signal s1, s2, s3, s4, s5, s6; /* signals to be used for the network manipulation */
    /* load the children of n and obtain the main  condition for the applicability of distributivity */
    bool n_has_two_levels = order_children( n, s1, s2 ); 
    
    node n1 = ntk.get_node( s1 );
    node n2 = ntk.get_node( s2 );
    
    /* first exit conditions */
    if ( ntk.level( n ) < 2 || ntk.is_pi( n1 ) || ntk.is_pi( n2 ) || n_has_two_levels )
      return false;

    order_children( n1, s3, s4 ); /* load the children of n1 */
    order_children( n2, s5, s6 ); /* load the children of n2 */
  

    if ( ( ntk.is_complemented( s1 ) == false ) && ( ntk.is_complemented( s2 ) == false ) ) /* phi1(s)=phi2(s)=s */
    { // S1
      signals_comparator cmp46( s4, s6 );
      signals_comparator cmp36( s3, s6 );

      /* if there is a child of n1 equal to a child of n2 they are labelled as s4 and s5 respectively */
      if( cmp46.abs || cmp36.abs )
        swap_elements( s5, s6 );
      if( cmp36.abs ) 
        swap_elements( s3, s4 );
      
      signals_comparator cmp45( s4, s5 );

      if ( cmp45.abs ) /* S1: |s5| = |s4| */
      {
        if ( cmp45.sign ) /* S1a: s5 = s4 */ 
        {
          /* f -> phi0( s4(s3 s6) )*/
          auto f1 = ntk.create_and( s3, s6 );
          auto f0 = ntk.create_and( s4, f1 );
          ntk.substitute_node( n, f0 );
          ntk.update_levels();
          node n_new = ntk.get_node( f0 ); /* the newly derived formula is a 2 level 2 nodes */
          return try_associativity( n_new ); /* return true */
        }
        else /* S1b: s5 = s4' */
        {
          /* f -> phi0( 0 ) */
          auto f0 = ntk.get_constant( false );
          ntk.substitute_node( n, f0 );
          return true;
        }
      }
      return false;
    }
    else if ( ( ntk.is_complemented( s1 ) == true ) && ( ntk.is_complemented( s2 ) == true ) ) /* phi1(s)=phi2(s)=s */
    {
      /* f = phi0(s3's5'+s3's6'+s4's5'+s4's6') */
      signals_comparator cmp46( s4, s6 );
      signals_comparator cmp36_init( s3, s6 );

      /* if present, place the children of n1 and n2 having the same index in s4 and s5 */
      if( cmp46.abs || cmp36_init.abs )
        swap_elements( s5, s6 );
      if( cmp36_init.abs )
        swap_elements( s3, s4 );
      
      signals_comparator cmp45( s4, s5 );
      signals_comparator cmp36( s3, s6 );


      if ( cmp45.abs ) /* |s5| = |s4| */
      {
        if ( cmp36.abs && ( cmp45.sign != cmp36.sign ) )
        {
          if ( cmp45.sign ) // S2a
          { /* s4=s5 s3=s6': f-> phi0(s4') */
            ntk.substitute_node( n, !s4 );
            return true;
          }
          else // S2b
          { /* s4=s5' s3=s6: f-> phi0(s3') */
            ntk.substitute_node( n, !s3 );
            return true;
          }
        }

        if ( cmp45.sign ) // S2c: s5 = s4 ( tested as Ditributivity )
        {
          /* f -> phi0'(s4(s3's6')') */
          auto f1 = ntk.create_and( !s3, !s6 );
          auto f0 = ntk.create_and( s4, !f1 );

          ntk.substitute_node( n, !f0 );
          ntk.update_levels();
          node n_new = ntk.get_node(f0);
          return try_associativity( n_new );
        }
      }
      return false;
    }
    else
    {
      if ( ( ntk.is_complemented( s1 ) == true ) && ( ntk.is_complemented( s2 ) == false ) )
      { // map the problem to phi1(s) = phi2(s)' = s 
        swap_elements( s1, s2 );
        swap_elements( n1, n2 );
        swap_elements( s3, s5 );
        swap_elements( s4, s6 );
      }
      /* f = phi0(s3s4s5'+s3s4s6') */
      signals_comparator cmp36( s3, s6 );
      signals_comparator cmp46( s4, s6 );

      /* if present, place the children of n1 and n2 having the same index in s4 and s5 */
      if( cmp36.abs || cmp46.abs )
      {
        swap_elements( s5, s6 );
      }
      signals_comparator cmp35( s3, s5 );
      if ( cmp35.abs )
        swap_elements( s3, s4 );

      signals_comparator cmp45( s4, s5 );

      if( cmp45.abs ) /* |s4|=|s5| */
      {
        if ( cmp45.sign ) /* s4=s5 */
        { /* f -> phi0(s6'(s3s4)) */
          auto f1 = ntk.create_and( s3, s4 );
          auto f0 = ntk.create_and( f1, !s6 );
          ntk.substitute_node( n, f0 );
          ntk.update_levels();
          node n_new = ntk.get_node(f0);
          return try_associativity( n_new );
        }
        else /* s4=s5' */
        { /* f -> phi0(s3s4) */
          auto f0 = ntk.create_and( s3, s4 );
          ntk.substitute_node( n, f0 );
          return true;
        }
      }
    }
    return false;
  }
  
  /* try to simplify the three level networks not reduced by the previous methods*/
  bool try_3_levels_3_nodes( node n )
  {
    signal s1, s2, s3, s4, s5, s6;
    order_children( n, s1, s2 );
    node n1 = ntk.get_node( s1 );
    node n2 = ntk.get_node( s2 );

    /* first filtering condition identifying the cases not included in the previous methods */
    if ( ntk.level( n ) < 3 || ntk.level( n2 ) < 2 || ntk.is_pi( n2 ) )
      return false;

    order_children( n2, s3, s4 );
    node n4 = ntk.get_node( s4 );
    auto lev1 = ntk.level( ntk.get_node( s1 ) );
    auto lev3 = ntk.level( ntk.get_node( s3 ) );
    auto lev4 = ntk.level( n4 );

    /* second filtering condition identifying the cases not included in the previous methods */
    if ( ntk.level( n4 ) < 1 || ntk.is_pi( n4 ) || ( lev4 <= lev3 ) || ( lev4 <= lev1 ) )
      return false;
    order_children( n4, s5, s6 );

    /* place in s5 the signal related to s1, if present */
    signals_comparator cmp16( s1, s6 );
    if( cmp16.abs )
      swap_elements( s5, s6 );

    signals_comparator cmp15( s1, s5 );
    
    if ( ntk.is_complemented( s2 ) ) // S1 phi2(x) = x' 
    {
      if( ntk.is_complemented( s4 ) ) // S1a phi4(x) = x' tested
      { /* f = phi0(s1s3'+s1s5s6) */
        if( cmp15.abs )
        {
          if ( cmp15.sign ) // S1a1 s1 = s5
          { /* f -> phi0(s1(s6's3)') */
            auto f1 = ntk.create_or( s6, !s3 );
            auto f0 = ntk.create_and( s1, f1 );

            ntk.substitute_node( n, f0 );
            ntk.update_levels();
            node n_new = ntk.get_node( f0 );
            try_associativity( n_new );
            return true;
          }
          else  //S1a2 s1 = s5'
          { /* f -> phi0((s1s3')) */
            auto f0 = ntk.create_and( s1, !s3 );
            ntk.substitute_node( n, f0 );
            return true;
          }
        }
        /* S1a3: f -> phi0'((s1s3')'(s1(s5s6))') */
        auto f0L = ntk.create_and( s1, !s3 );
        auto f1R = ntk.create_and( s1, s6 );
        auto f0R = ntk.create_and( f1R, s5 );
        auto f0 = ntk.create_and( !f0L, !f0R );
        ntk.substitute_node( n, !f0 );
        ntk.update_levels();
        node nR_new = ntk.get_node( f0R );
        return try_associativity( nR_new );
      }
      else // S1b: phi4(s) = s
      { /* f = phi0(s1s3'+s1s5'+s1s6') */
        if( cmp15.abs ) /* |s1|=|s5| */
        {
          if( cmp15.sign ) // S1b1: s1 = s5
          { /* f -> phi0(s1(s3s6)') */
            auto f1 = ntk.create_and( s3, s6 );
            auto f0 = ntk.create_and( s1, !f1 );
            ntk.substitute_node( n, f0 );
            ntk.update_levels();
            node n_new = ntk.get_node(f0);
            try_associativity( n_new );
            return true;
          }
          else // S1b2 : s1 = s5'
          { /* f -> phi0(s1) */
            ntk.substitute_node( n, s1 );
            return true;
          }
        }
      }
    }
    else // S2 phi2(s) = s
    {
      if( ntk.is_complemented( s4 ) ) // S2a: phi4(s)=s'
      { /* f = phi0(s1s3s5'+s1s3s6') */
        if ( cmp15.abs ) /* |s1|=|s5| */
        {
          if( cmp15.sign ) // S2a: s1=s5
          { /* f -> phi0(s1(s3s6')) */
            auto f1 = ntk.create_and( s3, !s6 );
            auto f0 = ntk.create_and( f1, s1 );
            ntk.substitute_node( n, f0 );
            ntk.update_levels();
            node n_new = ntk.get_node(f0);
            try_associativity( n_new );
            return true;
          }
          else // S2b : 
          { /* f -> phi0(s1s3) */
            auto f0 = ntk.create_and( s1, s3 );
            ntk.substitute_node( n, f0 );
            return true;
          }
        }
      }
      return false;
    }
    return false;
  }

#pragma region helper_functions
private:
  class signals_comparator
  {
    public:
      bool eq;    // true if the signals are equal
      bool abs;   // true if the signals have the same index
      bool sign;  // true if the signals are both complemented or not complemented

    signals_comparator( signal s1, signal s2 )
    {
      eq = ( s1 == s2 );
      abs = ( s1.index == s2.index );
      sign = ( s1.complement == s2.complement );
    }

  };

/* returns true if signal s1 originates from a node that is on a higher 
 * level with respect to the node from which originates s2
*/
  bool is_higher( signal s1, signal s2 ) 
  {
    return ( ntk.level( ntk.get_node( s1 ) ) > ntk.level( ntk.get_node( s2 ) ) );
  }

/* places the children signals of n in sL and sH, so that the levels of the corresponding nodes
 * satisfy level(nH) >= level(nL). If level(nH) > level(nL) returns true. 
*/
  bool order_children( node n, signal& sL, signal& sH )
  {

    bool has_two_levels = false;
    std::vector<signal> signals;
    std::vector<node> nodes;


    ntk.foreach_fanin( n, [&]( signal const& fi ){ 
      auto const& child = ntk.get_node( fi );
      nodes.emplace_back( child );
      signals.emplace_back( fi );
    });
    if ( ntk.level( nodes[0] ) > ntk.level( nodes[1] ) )
    {
      sL = signals[1];
      sH = signals[0];
      has_two_levels = true;
    }
    else if ( ntk.level( nodes[0] ) < ntk.level( nodes[1] ) )
    {
      sL = signals[0];
      sH = signals[1];
      has_two_levels = true;
    }
    else
    {
      sL = signals[0];
      sH = signals[1];
    }

    return has_two_levels;
  }
 
 /* swap two any elements */
  template<class T>
  void swap_elements( T& e1, T& e2 )
  {
    T e1_tmp = e1;
    e1 = e2;
    e2 = e1_tmp;
  }

#pragma endregion helper_functions
private:
  Ntk& ntk;
};

} // namespace detail

/* Entry point for users to call */
template<class Ntk>
void aig_algebraic_rewriting( Ntk& ntk )
{
  static_assert( std::is_same_v<typename Ntk::base_type, aig_network>, "Ntk is not an AIG" );

  depth_view dntk{ntk};
  detail::aig_algebraic_rewriting_impl p( dntk );
  p.run();
}

} /* namespace mockturtle */