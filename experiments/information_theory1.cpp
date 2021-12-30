#include <iostream>
//#include <catch.hpp>

#include <sstream>
#include <string>
#include <vector>
#include <bitset>
#include <cmath>
#include <random>
#include <boost/dynamic_bitset.hpp>

#include <mockturtle/networks/klut.hpp>

using namespace mockturtle;
using namespace kitty;

/* ##################### start: PLA Ntk #######################*/
struct index_to_signal
{
  index_to_signal()
  {
    storage.reserve( 10000u );
  }
  void insert( uint64_t pla_index, uint64_t klut_signal )
  {
    storage[pla_index] = klut_signal;
  }

  std::unordered_map<uint64_t, uint64_t> storage;
};

  class pla_network
  {
    #pragma region Types and constructors
    public:
      pla_network( std::vector<uint64_t> input_nodes, std::vector<uint64_t> output_nodes, uint32_t Nin )
      : _nodes( input_nodes ),
        _outputs( output_nodes ),
        _num_data( input_nodes.size()),
        _num_nodes( Nin )
        {
          _init();
        }

    protected:
      inline void _init()
      {

        for ( uint32_t i {1u}; i < _num_nodes ; ++i )
        {
          auto pi = klut.create_pi();
          _itos.insert( i, pi );
        }
        _act = 0;

      }
    #pragma endregion
    
    #pragma region visual
    public:
      void print_pla()
      {
        for ( uint32_t i {0u}; i < _num_data; ++i )
          {
            boost::dynamic_bitset<> BS ( _num_nodes, _nodes.at(i) );
            std::cout << BS << std::endl ;
          }
      }

      void print_probabilities( std::vector<float> probabilities )
      {
        uint32_t num_vars = probabilities.size();
        for ( uint32_t mask {0u}; mask < num_vars; ++mask )
        {
          boost::dynamic_bitset<> BS ( (uint32_t)log2(num_vars), mask );
          std::cout << "|P(" << BS << ") = " << probabilities.at(mask) << std::endl ;
        }
      }
    #pragma endregion
      
    #pragma region basic_function
    std::vector<float> Pr( std::vector<uint64_t> indeces_nodes, std::vector<uint64_t> indeces_outputs )
    {

      uint32_t size_P_space = std::pow( 2, indeces_nodes.size() + indeces_outputs.size() );
      std::vector<float> probabilities;
      float_t proba;
      double eq_flag_nodes, eq_flag_outputs; 

      for ( uint64_t xin {0u}; xin < size_P_space; ++xin )
      {
        uint64_t mask_nodes {0u};
        uint64_t X_nodes {0u};
        uint32_t jeff;
        for ( uint32_t j {0u}; j < indeces_nodes.size(); ++j )
        {
          jeff = indeces_outputs.size()+j;
          mask_nodes |= ( 1u << indeces_nodes.at(j) );
          X_nodes |= ( ( ( ( 1u << jeff ) & xin ) >> jeff ) << indeces_nodes.at(j) );
        }
        
        uint64_t mask_outputs {0u};
        uint64_t X_outputs {0u};
        for ( uint32_t j {0u}; j < indeces_outputs.size(); ++j )
        {
          mask_outputs |= ( 1u << indeces_outputs.at(j) );
          X_outputs |= ( ( ( ( 1u << j ) & xin ) >> j ) << indeces_outputs.at(j) );
        }
        proba = 0;
        //std::cout << xin << "-th pattern\n\n";
        for ( uint32_t i {0u}; i < _num_data; ++i )
        {
          //std::cout << "pat" << ( mask & _nodes.at(i) ) << "\n" ;
          //std::cout << "X=" << X << "\n";
          if ( ( indeces_nodes.size() != 0 ) && ( indeces_outputs.size() != 0 ) )
          {
            //std::cout << "both\n";
            eq_flag_nodes = ( X_nodes == ( mask_nodes & _nodes.at(i) ) ) ? 1 : 0;
            eq_flag_outputs = ( X_outputs == ( mask_outputs & _outputs.at(i) ) ) ? 1 : 0;
          }
          else if ( indeces_nodes.size() == 0 )
          {
            //std::cout << "out only\n";
            eq_flag_nodes = 1;
            eq_flag_outputs = ( X_outputs == ( mask_outputs & _outputs.at(i) ) ) ? 1 : 0;
          }
          else if ( indeces_outputs.size() == 0 )
          {
            //std::cout << "nodes only\n";
            eq_flag_nodes = ( X_nodes == ( mask_nodes & _nodes.at(i) ) ) ? 1 : 0;
            eq_flag_outputs = 1;
          }
          proba += eq_flag_outputs*eq_flag_nodes/_num_data;

          //std::cout << "p: " << proba << std::endl;
        }
        probabilities.push_back( proba );
      }
      return probabilities;

    }

    float H( std::vector<uint64_t> indeces_nodes, std::vector<uint64_t> indeces_outputs )
    {
      auto proba = Pr( indeces_nodes, indeces_outputs );
      uint32_t size_P_space = proba.size();
      float entropy { 0 };
      float deltaH { 0 };

      for ( uint64_t xin {0u}; xin < size_P_space; ++xin )
      {
        deltaH = ( proba[xin] == 0 ) ? 0 : -1*proba[xin]*log2( proba[xin] );
        entropy += deltaH;
      } 
      return entropy;
    }

    float MI ( std::vector<uint64_t> Xindeces, std::vector<uint64_t> Yindeces )
    {
      auto Hx = H( Xindeces, {} );
      auto Hy = H( {}, Yindeces );
      auto Hxy = H( Xindeces, Yindeces ); 
      //std::cout << "H(X)=" << Hx << std::endl;
      //std::cout << "H(Y)=" << Hy << std::endl;
      //std::cout << "H(X, Y)=" << Hxy << std::endl;

      return ( Hx + Hy - Hxy );
    }

    #pragma endregion

    #pragma region new_node
    void fill_active_list( )
    {

      float mi_loc;
      float mi_max = 0;
      uint32_t idx;
      /* first active variable */
      for( uint32_t i {0u}; i < _num_nodes; ++i )
      {
        mi_loc = MI( {i}, {0} );
        uint32_t idx; 
        //std::cout << "mi_loc " << mi_loc << std::endl;
        //std::cout << "mi_max " << mi_max << std::endl;
        if ( mi_loc >= mi_max )
        {
          mi_max = mi_loc;
          idx = i;
        }
        //std::cout << "idx: " << idx << std::endl;
        _active_list = {idx};
      }

      /*std::cout << "AL pre :\n";
      for( uint32_t k {0u}; k < _active_list.size(); ++k )
      {
        std::cout << _active_list[k] _active_list<< " ";
      }*/
      /* next active variables */
      std::vector<uint64_t> inv_indeces;

      for ( uint32_t i {1u}; i < _num_nodes; ++i )
      {
        mi_max = 0;
        inv_indeces = _active_list;
        inv_indeces.emplace_back(0);
        for ( uint32_t j {0u}; j < _num_nodes; ++j )
        { 
          //std::cout << "in " << j << std::endl;

          if ( std::find( _active_list.begin(), _active_list.end(), j ) != _active_list.end() )
            continue;
          //std::cout << "not cont " << j << std::endl;
          
          inv_indeces.at(i) = j;

          /*std::cout << "inv indeces:\n";
          for( uint32_t k {0u}; k <= i; ++k )
          {
            std::cout << inv_indeces.at(k) << " ";
          }
          */
          mi_loc = MI( inv_indeces, {0} );
          if ( mi_loc >= mi_max )
          {
            mi_max = mi_loc;
            idx = j;
          }
        }
      _active_list.push_back(idx);

      }
    }

    std::string create_fn( std::vector<uint32_t> support )
    {

      /*std::cout << "AL:\n";
      for( uint32_t k {0u}; k < nact; ++k )
      {
        std::cout << _active_list[k] << " ";
      }
      std::cout << "\n\n";*/
      uint32_t nin_node = support.size();
      std::cout << "supp size = " << nin_node << std::endl;
      uint32_t domain_size = pow( 2, nin_node );
      uint32_t Ci0, Ci1;
      uint64_t mask {0u};
      uint64_t X {0u};
      std::default_random_engine generator;
      std::bernoulli_distribution distribution(0.5);
      std::string tt_str;
      
      uint64_t mask0 = ~( 1u << _num_nodes );
      for ( uint32_t j {0u}; j < _num_data; ++j )
      {
        _nodes.at(j) &= mask0; 
      }

      for ( uint32_t xin {0u}; xin < domain_size; ++xin )
      {
        Ci0 = 0;
        Ci1 = 0;
        mask = 0u;
        X = 0u;

        for ( uint32_t j {0u}; j < support.size(); ++j )
        {
          //std::cout << "m: " << mask << std::endl;
          mask |= ( 1u << support.at(j) );
          //std::cout << "m: " << mask << std::endl;

          X |= ( ( ( ( 1u << j ) & xin ) >> j ) << support.at(j) );
        }

        for ( uint32_t j {0u}; j < _num_data; ++j )
        {
          if ( X == ( mask & _nodes.at(j) ) )
            ( ( _outputs.at(j) & 1u ) == 1u ) ? Ci1++ : Ci0++;
        }
        uint64_t new_val {0u};
        if( Ci1 > Ci0 )
        {
          new_val = 1u << _num_nodes;
          tt_str = "1" + tt_str;
        }
        else if( Ci1 == Ci0 )
        {
          if (distribution(generator))
          {
            new_val = 1u << _num_nodes;
            tt_str = "1" + tt_str;
          }
          else
          {
            tt_str = "0" + tt_str;
          }
        }
        else
        {
          tt_str = "0" + tt_str;
        }

        for ( uint32_t j {0u}; j < _num_data; ++j )
        {
          if ( X == ( mask & _nodes.at(j) ) )
            _nodes.at(j) |= new_val;
        }
      }

      return tt_str;

    }

    void create_klut_node( std::vector<uint32_t> support, std::string tt_str )
    {
      auto supp_size = support.size();
      dynamic_truth_table tt( supp_size );
      create_from_binary_string( tt, tt_str );
      std::vector<uint64_t> klut_signals;
      for ( uint32_t i {0u}; i < supp_size; ++i )
        {
          klut_signals.push_back( _itos.storage[support[i]] );
        }
      auto f0 = klut.create_node( klut_signals, tt );
      _itos.insert( _num_nodes ,f0 );
      _num_nodes++;
    }

    void muesli( uint32_t nact )
    {
      std::string tt_str;

      fill_active_list( );
      std::cout << "\n";

      for(uint32_t k{0u}; k<_active_list.size(); ++k)
      {
        std::cout << _active_list[k] << " ";
      }
      std::cout << "\n";
      
      std::vector<uint32_t> support;
      for( uint32_t act_idx {_act}; act_idx < ( _num_nodes - nact ) ; ++act_idx )
      {
        support = {};
        std::cout << "AL:\n";
        for ( uint32_t k{0u}; k < nact; ++k )
        {
          support.push_back( _active_list.at( act_idx + k ));
          std::cout << _active_list.at( act_idx + k ) << std::endl;
        }
        auto mi_old = MI( {support.at( 0 )}, {0} );

        tt_str = create_fn( support );
        std::cout << tt_str << std::endl;

        auto mi_new = MI( {_num_nodes}, {0} );
        std::cout << "mi_new " << mi_new << std::endl;
        std::cout << "mi_old " << mi_old << std::endl;

        if ( mi_new > mi_old )
        {
          create_klut_node( support, tt_str );
          break;
        }
        else
        {
          _act++;
        }
      }
    }
    #pragma endregion


    public:
      std::vector<uint64_t> _nodes; /* storage element: value of the output at each example */
      std::vector<uint64_t> _outputs; /* storage element: value of the output at each example */
      uint32_t _num_data; /* number of examples */
      uint32_t _num_nodes;
      klut_network klut;
      std::vector<uint64_t> _active_list;
      index_to_signal _itos;
      uint32_t _act;

  };
/* #####################  end: PLA Ntk #######################*/


int main()
{

  std::vector<uint64_t> input_nodes;
  for ( uint32_t i {0u}; i < 32; ++i )
  {
    input_nodes.emplace_back(i);
  }

  std::vector<uint64_t> output_nodes;
  output_nodes = {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1}; /* ab + cde */
  //output_nodes = {0,0,0,1,0,0,0,1,0,0,0,1,1,1,1,1};

  pla_network pla( input_nodes, output_nodes, 5 );
  pla.print_pla();
  pla.muesli(2);
  pla.print_pla();
  pla.muesli(2);
  pla.print_pla();
  pla.muesli(2);
  pla.print_pla();
  pla._act = 0;
  pla.muesli(2);

  std::cout << "|P(f):\n";
  auto probs = pla.Pr( {0}, {} );
  pla.print_probabilities( probs );
  probs = pla.Pr( {}, {0} );
  pla.print_probabilities( probs );
  probs = pla.Pr( {0}, {0} );
  pla.print_probabilities( probs );
  /*std::cout << "####### pla created: some stats ############\n";
  std::cout << "|P(f):\n";
  auto probs = pla.Pr( {}, {0} );
  pla.print_probabilities( probs );
  std::cout << "\n|P(a):\n";
  probs = pla.Pr( {4}, {} );
  pla.print_probabilities( probs );
  std::cout << "\n|P(e):\n";
  probs = pla.Pr( {0}, {} );
  pla.print_probabilities( probs );
  std::cout << "\n|P(a,b):\n";
  probs = pla.Pr( {2,3}, {} );
  pla.print_probabilities( probs );
  std::cout << "\n####### Look for first function ############\n";

  std::cout << "MI(a,f)=" << pla.MI( {4}, {0} ) << std::endl;
  std::cout << "MI(b,f)=" << pla.MI( {3}, {0} ) << std::endl;
  std::cout << "MI(c,f)=" << pla.MI( {2}, {0} ) << std::endl;
  std::cout << "MI(d,f)=" << pla.MI( {1}, {0} ) << std::endl;
  std::cout << "MI(e,f)=" << pla.MI( {0}, {0} ) << std::endl;

  std::cout << "This means that either a or b are picked as first choice. Suppose a is picked. The second variable is chosen to maximize MI(a||x;f)\n";
  std::cout << "MI(a||b,f)=" << pla.MI( {4,3}, {0} ) << std::endl;
  std::cout << "MI(a||c,f)=" << pla.MI( {4,2}, {0} ) << std::endl;
  std::cout << "MI(a||d,f)=" << pla.MI( {4,1}, {0} ) << std::endl;
  std::cout << "MI(a||e,f)=" << pla.MI( {4,0}, {0} ) << std::endl;
  std::cout << "\nclearly the next chosen is b. \n";
  std::cout << "This whole thing is handled by the method fill_active list. Indeed, the active list is \n";
  pla.fill_active_list(2);
  for( uint32_t i {0u}; i < pla._active_list.size(); ++i )
  {
    std::cout << pla._active_list.at(i) << std::endl;
  }
  std::cout << "where 0->e, 1->d, 2->c, 3->b, 4->a\n"; 
  pla.create_node(2);
  std::cout << "new node created. Located in 5:\n";
  pla.print_pla();
  std::cout << "If you now want to create a new node, first find the most relevant node:\n";
  std::cout << "\nMI(g5(a,b),f)=" << pla.MI( {5}, {0} ) << std::endl;
  std::cout << "MI(a,f)=" << pla.MI( {4}, {0} ) << std::endl;
  std::cout << "MI(b,f)=" << pla.MI( {3}, {0} ) << std::endl;
  std::cout << "MI(c,f)=" << pla.MI( {2}, {0} ) << std::endl;
  std::cout << "MI(d,f)=" << pla.MI( {1}, {0} ) << std::endl;
  std::cout << "MI(e,f)=" << pla.MI( {0}, {0} ) << std::endl;
  std::cout << "\n=> the first element of the active list is g5. Go for the second:";
  std::cout << "\nMI(g5(a,b)||a,f)=" << pla.MI( {5,4}, {0} ) << std::endl;
  std::cout << "MI(g5(a,b)||b,f)=" << pla.MI( {5,3}, {0} ) << std::endl;
  std::cout << "MI(g5(a,b)||c,f)=" << pla.MI( {5,2}, {0} ) << std::endl;
  std::cout << "MI(g5(a,b)||d,f)=" << pla.MI( {5,1}, {0} ) << std::endl;
  std::cout << "MI(g5(a,b)||e,f)=" << pla.MI( {5,0}, {0} ) << std::endl;
  std::cout << " the second is c, go for the 3-rd";
  std::cout << "\nMI(g5(a,b)||c||a,f)=" << pla.MI( {5,2,4}, {0} ) << std::endl;
  std::cout << "MI(g5(a,b)||c||b,f)=" << pla.MI( {5,2,3}, {0} ) << std::endl;
  std::cout << "MI(g5(a,b)||c||d,f)=" << pla.MI( {5,2,1}, {0} ) << std::endl;
  std::cout << "MI(g5(a,b)||c||e,f)=" << pla.MI( {5,2,0}, {0} ) << std::endl;
  std::cout << "\n OK:NOT CLEAR: what if the active list is to be filled maximizing the info so far?\n";
  std::cout << "\nMI(g5(a,b)||a,f)=" << pla.MI( {5,4}, {0} ) << std::endl;
  std::cout << "MI(g5(a,b)||b,f)=" << pla.MI( {5,3}, {0} ) << std::endl;
  std::cout << "MI(g5(a,b)||c,f)=" << pla.MI( {5,2}, {0} ) << std::endl;
  std::cout << "MI(g5(a,b)||d,f)=" << pla.MI( {5,1}, {0} ) << std::endl;
  std::cout << "MI(g5(a,b)||e,f)=" << pla.MI( {5,0}, {0} ) << std::endl;
  std::cout << "Pick c\n";
  std::cout << "\nMI(g5(a,b)||c||a,f)=" << pla.MI( {5,2,4}, {0} ) << std::endl;
  std::cout << "MI(g5(a,b)||c||b,f)=" << pla.MI( {5,2,3}, {0} ) << std::endl;
  std::cout << "MI(g5(a,b)||c||d,f)=" << pla.MI( {5,2,1}, {0} ) << std::endl;
  std::cout << "MI(g5(a,b)||c||e,f)=" << pla.MI( {5,2,0}, {0} ) << std::endl;
  std::cout << "If we use only c?\n";
  std::cout << "\nMI(c||a,f)=" << pla.MI( {2,4}, {0} ) << std::endl;
  std::cout << "MI(c||b,f)=" << pla.MI( {2,3}, {0} ) << std::endl;
  std::cout << "MI(c||d,f)=" << pla.MI( {2,1}, {0} ) << std::endl;
  std::cout << "MI(c||e,f)=" << pla.MI( {2,0}, {0} ) << std::endl;
  std::cout << "DONE\n";
  std::cout << "MI(g5||c||d||e,f)=" << pla.MI( {5,2,1,0}, {0} ) << std::endl;

  pla.fill_active_list(2);
  pla.create_node(2);

  for( uint32_t i {0u}; i < pla._active_list.size(); ++i )
  {
    std::cout << pla._active_list.at(i) << std::endl;
  }
  std::cout << "MI(g(g5(a,b),c),f)=" << pla.MI( {6}, {0} ) << std::endl;*/

  return 0;
}