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
      : _inputs( input_nodes ),
        _outputs( output_nodes ),
        _num_data( input_nodes.size()),
        _num_nodes( Nin )
        {
          _init();
        }

    protected:
      inline void _init()
      {
        for ( uint32_t i {0u}; i < _num_data; ++i )
        {
          _nodes.push_back(( _inputs.at(i) << 1u ) | _outputs.at(i)) ;
        }

        for ( uint32_t i {1u}; i < _num_nodes ; ++i )
        {
          auto pi = klut.create_pi();
          _itos.insert( i, pi );

        }

        _num_nodes += 1;
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
    std::vector<float> Pr( std::vector<uint64_t> indeces )
    {

      uint32_t size_P_space = std::pow( 2, indeces.size() );
      std::vector<float> probabilities;
      float_t proba;
      double eq_flag; 

      for ( uint64_t xin {0u}; xin < size_P_space; ++xin )
      {
        uint64_t mask {0u};
        uint64_t X {0u};

        for ( uint32_t j {0u}; j < indeces.size(); ++j )
        {
          //std::cout << "m: " << mask << std::endl;
          mask |= ( 1u << indeces.at(j) );
          //std::cout << "m: " << mask << std::endl;

          X |= ( ( ( ( 1u << j ) & xin ) >> j ) << indeces.at(j) );
        }

        proba = 0;
        //std::cout << xin << "-th pattern\n\n";
        for ( uint32_t i {0u}; i < _num_data; ++i )
        {
          //std::cout << "pat" << ( mask & _nodes.at(i) ) << "\n" ;
          //std::cout << "X=" << X << "\n";

          eq_flag = ( X == ( mask & _nodes.at(i) ) ) ? 1 : 0;
          //std::cout << "ef: " << eq_flag << std::endl;
          proba += eq_flag/_num_data;
          //std::cout << "p: " << proba << std::endl;
        }
        probabilities.push_back( proba );
      }
      return probabilities;

    }

    float H( std::vector<uint64_t> indeces )
    {
      auto proba = Pr( indeces );
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
      auto Hx = H( Xindeces );
      auto Hy = H( Yindeces );
      Xindeces.insert(Xindeces.end(), Yindeces.begin(), Yindeces.end()); 
      auto Hxy = H( Xindeces ); 
      std::cout << "H(X)=" << Hx << std::endl;
      std::cout << "H(Y)=" << Hy << std::endl;
      std::cout << "H(X, Y)=" << Hxy << std::endl;

      return ( Hx + Hy - Hxy );
    }

    #pragma endregion

    #pragma region new_node
    void fill_active_list( uint32_t nact )
    {

      float mi_loc;
      float mi_max = 0;
      uint32_t idx;
      /* first active variable */
      for( uint32_t i {1u}; i < _num_nodes; ++i )
      {
        mi_loc = MI( {i}, {0} );
        uint32_t idx; 
        std::cout << "mi_loc " << mi_loc << std::endl;
        std::cout << "mi_max " << mi_max << std::endl;
        if ( mi_loc >= mi_max )
        {
          mi_max = mi_loc;
          idx = i;
        }
        std::cout << "idx: " << idx << std::endl;
        _active_list = {idx};
      }
      /* next active variables */
      std::vector<uint64_t> inv_indeces;

      for ( uint32_t i {1u}; i < nact; ++i )
      {
        mi_max = 0;
        inv_indeces = _active_list;
        inv_indeces.emplace_back(0);
        for ( uint32_t j {1u}; j < _num_nodes; ++j )
        { 
          if ( std::find( _active_list.begin(), _active_list.end(), j ) != _active_list.end() )
            continue;
          inv_indeces.at(i) = j;
          mi_loc = MI( inv_indeces, {0} );
          if ( mi_loc >= mi_max )
          {
            mi_max = mi_loc;
            idx = i;
          }
        }
      _active_list.push_back(idx);

      }
    }

    void create_node( uint32_t nact )
    {
      fill_active_list( nact );
      
      uint32_t nin_node = _active_list.size();
      uint32_t domain_size = pow( 2, nin_node );
      uint32_t Ci0, Ci1;
      uint64_t mask {0u};
      uint64_t X {0u};
      std::default_random_engine generator;
      std::bernoulli_distribution distribution(0.5);
      std::string tt_str;

      for ( uint32_t xin {0u}; xin < domain_size; ++xin )
      {
        Ci0 = 0;
        Ci1 = 0;
        mask = 0u;
        X = 0u;

        for ( uint32_t j {0u}; j < _active_list.size(); ++j )
        {
          //std::cout << "m: " << mask << std::endl;
          mask |= ( 1u << _active_list.at(j) );
          //std::cout << "m: " << mask << std::endl;

          X |= ( ( ( ( 1u << j ) & xin ) >> j ) << _active_list.at(j) );
        }

        for ( uint32_t j {0u}; j < _num_data; ++j )
        {
          if ( X == ( mask & _nodes.at(j) ) )
            ( ( _nodes.at(j) & 1u ) == 1u ) ? Ci1++ : Ci0++;
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

      std::cout << tt_str << std::endl;
      dynamic_truth_table tt( nact );
      create_from_binary_string( tt, tt_str );
      std::vector<uint64_t> klut_signals;
      for ( uint32_t i {0u}; i < nact; ++i )
      {
        klut_signals.push_back( _itos.storage[_active_list[i]] );
      }
      auto f0 = klut.create_node( klut_signals, tt );
      _itos.insert( _num_nodes ,f0 );
      _num_nodes++;

    }
    #pragma endregion


    public:
      std::vector<uint64_t> _inputs;   /* storage element: value of the node at each example */
      std::vector<uint64_t> _outputs; /* storage element: value of the output at each example */
      uint32_t _num_data; /* number of examples */
      uint32_t _num_nodes;
      std::vector<uint64_t> _nodes; /* storage element: value of the output at each example */
      klut_network klut;
      std::vector<uint64_t> _active_list;
      index_to_signal _itos;

  };
/* #####################  end: PLA Ntk #######################*/


int main()
{

  std::vector<uint64_t> input_nodes;
  for ( uint32_t i {0u}; i < 4; ++i )
  {
    input_nodes.emplace_back(i);
  }

  std::vector<uint64_t> output_nodes;
  output_nodes = {0,0,0,1};

  pla_network pla( input_nodes, output_nodes, 2 );
  pla.print_pla();
  auto probs = pla.Pr( {0, 1} );
  pla.print_probabilities( probs );
  std::cout << probs[0] << std::endl;
  std::cout << probs[1] << std::endl;
  std::cout << probs[2] << std::endl;
  std::cout << probs[3] << std::endl;

  std::cout << "Entropy:" << std::endl;
  std::cout << "H=" << pla.H( {1} ) << std::endl;
  std::cout << "H=" << pla.H( {0} ) << std::endl;
  std::cout << "MI=" << pla.MI( {0}, {1} ) << std::endl;

  std::cout << "Fill active list\n";
  pla.fill_active_list(2);
  for( uint32_t i {0u}; i < pla._active_list.size(); ++i )
  {
    std::cout << pla._active_list.at(i) << std::endl;
  }

  pla.create_node(2);
  pla.print_pla();
  

  std::cout << "DONE\n";

  return 0;
}