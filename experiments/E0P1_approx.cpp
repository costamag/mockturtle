#include "experiments.hpp"
#include <mockturtle/traits.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/resyn_engines/xag_resyn.hpp>
//#include <mockturtle/algorithms/experimental/contest.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/detail/mffc_utils.hpp>
#include <mockturtle/algorithms/lfe/muesli.hpp>
#include <mockturtle/algorithms/lfe/sim_muesli.hpp>
//#include <mockturtle/algorithms/lfe/sim_decomposition.hpp>
//#include <mockturtle/algorithms/lfe/sim_decomposition_dp.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <lorina/aiger.hpp>
#include <iostream>
#include <string>
#include <set>

using namespace mockturtle;
//using namespace mockturtle::experimental;
using namespace experiments;

std::vector<int> combination;

void pretty_print(const std::vector<int>& v) {
  static int count = 0;
  std::cout << "combination no " << (++count) << ": [ ";
  for (int i = 0; i < v.size(); ++i) { std::cout << v[i] << " "; }
  std::cout << "] " << std::endl;
}

std::vector<std::vector<int>> combinations; 

void go(std::vector<uint32_t> Bits, int offset, int k) {
  if (k == 0) {
    combinations.push_back( combination );
    //pretty_print(combination);
    return;
  }
  for (int i = offset; i <= Bits.size() - k; ++i) {
    combination.push_back(Bits[i]);
    go(Bits, i+1, k-1);
    combination.pop_back();
  }
}

double erase_and_verify( std::vector<kitty::partial_truth_table> pats, kitty::partial_truth_table Y, uint32_t num_erase, std::string Func )
{
  std::vector<sim_pattern<klut_network>> Xsp;
  for( auto tt : pats )
  {
    sim_pattern<klut_network> xsp( tt );
    Xsp.push_back(xsp);
  } 
  std::vector<uint32_t> reduced_support = {0,1,2};


  std::vector<uint32_t> erasable_bits;
  for (int i = 0; i < Y.num_bits(); ++i) { erasable_bits.push_back(i); }
  combinations = {};

  go( erasable_bits, 0, num_erase );
  
  uint32_t count = 0;
  for( auto x : combinations )
  {
    kitty::partial_truth_table on_f( Y.num_bits() ), amask1, amask0;
    on_f = Y;
    //kitty::print_binary( on_f ); std::cout << std::endl;
    amask1 = pats[3];
    //kitty::print_binary( amask1 ); std::cout << std::endl;
    amask0 = ~pats[3];
    //kitty::print_binary( amask0 ); std::cout << std::endl;

    for( uint32_t k = 0; k < x.size(); ++k )
    {
      kitty::clear_bit( amask0, x[k] );
      kitty::clear_bit( amask1, x[k] );
      kitty::clear_bit( on_f, x[k] );
    }

    std::reverse(x.begin(), x.end());


    //for( uint32_t i = 0; i < x.size(); ++i )
    //{
      kitty::partial_truth_table ytt = Y;
      std::vector<kitty::partial_truth_table> Xtt = pats;
      for( uint32_t j = 0; j < pats.size(); ++j )
      {
        for( uint32_t k = 0; k < x.size(); ++k )
          Xtt[j].erase_bit_shift( x[k] );
      }
  
      std::vector<kitty::partial_truth_table*> X;
      X.push_back( &Xtt[0] );
  

      for( uint32_t k = 0; k < x.size(); ++k )
        ytt.erase_bit_shift( x[k] );

//std::cout << "inputs:" << std::endl;
      std::vector<uint32_t> max_indeces;
      double Imax = 0;
      double Inew;
      std::cout << num_erase << " ";
      for( uint32_t i = 0; i < pats.size(); ++i )
      {
        X[0] = &Xtt[i];
        Inew = kitty::mutual_information( X, &ytt );
        std::cout << i << " " << Inew;
        if( Inew == Imax )
          max_indeces.push_back( i );
        else if( Inew > Imax )
        {
          Imax = Inew;
          max_indeces = std::vector{ i };
        }
      }
      std::cout << std::endl;
      for( auto bidx : max_indeces )
      {
        //std::cout << bidx << std::endl;
        if( bidx == 3 )
        {
          sim_top_decomposition_fast res = is_top_decomposable_fast( Xsp, reduced_support, on_f, amask1, amask0 );

          if( Func == "and" )
          {
            //kitty::print_binary( amask0 & on_f ); std::cout << std::endl;
            count += 1;//double( res == sim_top_decomposition_fast::and_ );
          }
          else if ( Func == "xor" )
            count += 1;
        }
      }
      //for( uint32_t i = 0 ; i < x.size(); ++i )
      //  std::cout << x[i] << " ";
    //}
  }
  std::cout << num_erase << " " << (double)count/(double)combinations.size() << std::endl;
  return (double)count/(double)combinations.size();
}


template<class Ntk, class TT, class RewritingFn>
std::vector<double> database_lookup( TT const& target, RewritingFn const& rewriting_fn_, std::string const& func )
{
  using signal = typename Ntk::signal;
  std::optional<uint32_t> res;

  // first create a network 
  uint32_t num_pis = target.num_vars();
  Ntk ntk;
  std::vector<signal> pis( num_pis );
  std::generate( std::begin(pis), std::end(pis), [&](){ return ntk.create_pi(); } );
  
  signal osig;
  // run npn resynthesis
  auto const on_signal = [&]( auto const& s ) {
    uint32_t _num_nodes = mockturtle::detail::recursive_ref<Ntk>( ntk, ntk.get_node( s ) );
    mockturtle::detail::recursive_deref<Ntk>( ntk, ntk.get_node( s ) );
    if ( !res || *res > _num_nodes ) res = _num_nodes;
    osig = s;
    return true;
  };
  rewriting_fn_( ntk, target, std::begin(pis), std::end(pis), on_signal );
  
  //std::cout << ntk.num_gates() << " -> ";
  signal x = ntk.create_pi();
  signal f0;
  if( func == "xor" )
    f0 = ntk.create_xor( x, osig );
  else if( func == "and" )
    f0 = ntk.create_and( x, osig );
  else if( func == "or" )
    f0 = ntk.create_or( x, osig );
  else if( func == "lt" )
    f0 = ntk.create_lt( x, osig );
  else if( func == "le" )
    f0 = ntk.create_le( x, osig );
  else
    std::cerr << "error" << std::endl;
  ntk.create_po( f0 );
  
  //std::cout << ntk.num_gates() << std::endl;


  std::vector<kitty::partial_truth_table> pats;
  kitty::partial_truth_table sim_pat( std::pow(2,pis.size()+1) );
  for( uint32_t i = 0; i < pis.size()+1; ++i )
  {
    kitty::create_nth_var( sim_pat, i );
    pats.push_back( sim_pat );
  }

  partial_simulator sim( pats );

  unordered_node_map<kitty::partial_truth_table, Ntk> node_to_value( ntk );
  simulate_nodes( ntk, node_to_value, sim );

/*
  kitty::print_binary(node_to_value[f0]);
  std::cout << std::endl; 
*/
  double fraction_right;
  std::vector<double> fractions;
  for( uint32_t i = 0; i < pats[0].num_bits(); ++i )
  {
    kitty::partial_truth_table Y = node_to_value[f0];
    fraction_right = erase_and_verify( pats, Y, i, func );
    fractions.push_back( fraction_right );
  }
  return fractions;
}


template<uint32_t num_vars>
std::vector<std::vector<double>> test_n_var_function( std::string func, bool only_npn = false )
{
  using TT = kitty::dynamic_truth_table;
  using Ntk = xag_network;
  
  // prepare the database for lookup
  xag_npn_resynthesis<Ntk> resyn;
  xag_npn_resynthesis<Ntk, Ntk, xag_npn_db_kind::aig_complete > resyn_complete;

  std::vector<std::vector<double>> npn_fractions;
  TT target(num_vars);
  std::set<TT> reprs;
  do
  {
    const auto repr = kitty::exact_npn_canonization(target);

    if( only_npn && ( reprs.find( std::get<0>(repr) ) != reprs.end() ) )
    {
      kitty::next_inplace( target );
      continue;
    }
    else
    {
      reprs.insert( std::get<0>(repr) );
      npn_fractions.push_back(database_lookup<Ntk, TT>( target, resyn_complete, func ));
      

      kitty::next_inplace( target );
    }
  } while ( !kitty::is_const0( target ) );
  return npn_fractions;
}

int main()
{

  // for ( auto const& benchmark : epfl_benchmarks() )
  // {
  //   // fmt::print( "[i] processing {}\n", benchmark );

  //   aig_network aig;
  //   auto const result = lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) );
  //   assert( result == lorina::return_code::success ); (void)result;
    
  //   uint32_t tot_mffc = 0;
  //   uint32_t big_mffc = 0;
  //   mockturtle::detail::initialize_values_with_fanout( aig );
  //   aig.foreach_node( [&]( auto n ){
  //     if ( mockturtle::detail::mffc_size( aig, n ) > 3 ) big_mffc++;
  //     tot_mffc++; 
  //   } );

  //   fmt::print( "{:.2f}\n", (double)big_mffc/tot_mffc );

  // }

  auto npn_fractions = test_n_var_function<3>("and", true);
  for( auto f : npn_fractions )
  {
    for( uint32_t i = 0; i < f.size(); ++i )
      std::cout << f[i] << " ";
    std::cout << std::endl;
  }

  return 0;
}