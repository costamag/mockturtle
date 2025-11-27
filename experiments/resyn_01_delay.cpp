/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2024  EPFL
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

#include <iostream>
#include <string>
#include <vector>
#include <chrono>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/genlib.hpp>
#include <mockturtle/algorithms/aig_balancing.hpp>
#include <mockturtle/algorithms/emap.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/rewrite.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/block.hpp>
#include <mockturtle/networks/klut.hpp>

#include <mockturtle/utils/name_utils.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/cell_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/names_view.hpp>

#include <mockturtle/resynthesis/aig_syn.hpp>
#include <mockturtle/resynthesis/algorithms/resynthesize.hpp>
//#include <mockturtle/resynthesis/analyzers/homo_xpx_analyzer.hpp>
#include <mockturtle/resynthesis/power/node_simulator.hpp>
#include <mockturtle/resynthesis/power/analyzers/power_analyzer.hpp>
#include <mockturtle/resynthesis/database/homo_database.hpp>
#include <mockturtle/resynthesis/database/mapped_list.hpp>
#include <mockturtle/resynthesis/list_simulator.hpp>
#include <mockturtle/resynthesis/networks/bindings.hpp>
#include <mockturtle/resynthesis/networks/bound.hpp>
#include <mockturtle/resynthesis/database/homo_storage.hpp>
#include <mockturtle/resynthesis/networks/mapped_to_bound.hpp>
#include <mockturtle/resynthesis/networks/cleanup_bound.hpp>
#include <mockturtle/resynthesis/power/analyzers/convo_analyzer.hpp>

#include <experiments.hpp>

template<class Ntk>
inline void abc_baseline( Ntk const& ntk, std::string const& library_full_path, std::string const& benchmark )
{

  mockturtle::write_verilog_with_binding( ntk, "/tmp/ntk.v" );
  std::string command = fmt::format( "abc -q \"read_genlib {}; read_verilog -m /tmp/ntk.v; ps -p; strash; dch; if -p; strash; dch; if -p; mfs -p; mfs -p; map -ap; ps -p; write_blif /home/costamag/POST_MAPPING/mockturtle/build/glitch/{}_abc.blif; \"", library_full_path, benchmark );

  std::array<char, 128> buffer;
  std::string result;
#if WIN32
  std::unique_ptr<FILE, decltype( &_pclose )> pipe( _popen( command.c_str(), "r" ), _pclose );
#else
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "w" ), pclose );
#endif
  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  return;
}

template<class Ntk>
inline void abc_powermap( Ntk const& ntk, std::string const& library_full_path, std::string const& benchmark )
{

  mockturtle::write_verilog_with_binding( ntk, "/tmp/ntk2.v" );
  std::string command = fmt::format( "abc -q \"read_genlib {}; read_verilog -m /tmp/ntk2.v; ps -p; strash; dch; if -p; strash; dch; if -p; map -p; ps -p; write_blif /home/costamag/POST_MAPPING/mockturtle/build/glitch/{}_abc2.blif; \"", library_full_path, benchmark );

  std::array<char, 128> buffer;
  std::string result;
#if WIN32
  std::unique_ptr<FILE, decltype( &_pclose )> pipe( _popen( command.c_str(), "r" ), _pclose );
#else
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "w" ), pclose );
#endif
  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  return;
}

template<class Ntk>
inline bool abc_power_eval( Ntk const& ntk, std::string const& library_full_path )
{
  mockturtle::write_verilog_with_binding( ntk, "/tmp/ntk.v" );
  std::string command = fmt::format( "abc -q \"read_genlib {}; read_verilog -m /tmp/ntk.v; ps -pg;\"", library_full_path );

  std::array<char, 128> buffer;
  std::string result;
#if WIN32
  std::unique_ptr<FILE, decltype( &_pclose )> pipe( _popen( command.c_str(), "r" ), _pclose );
#else
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
#endif
  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }

  /* search for one line which says "Networks are equivalent" and ignore all other debug output from ABC */
  std::stringstream ss( result );
  std::string line;
  while ( std::getline( ss, line, '\n' ) )
  {
    std::cout << line << std::endl;
    if ( line.size() >= 23u && line.substr( 0u, 23u ) == "Networks are equivalent" )
    {
      return true;
    }
  }

  return false;
}

void build_recursive( auto const& klut, auto& bound, auto f, auto & old_to_new, auto tt_to_id, auto gates )
{
  auto n = klut.get_node( f );
  if( old_to_new.find( n ) != old_to_new.end() )
    return;
  std::vector<mockturtle::bound_network::signal> children;
  klut.foreach_fanin( n, [&]( auto fi, auto ii ){
    auto ni = klut.get_node( fi );
    build_recursive( klut, bound, fi, old_to_new, tt_to_id, gates );
    children.push_back( old_to_new[ni] );
  } );
  auto tt = klut.node_function( n );
  uint32_t id = tt_to_id[tt];
  std::string name = gates[ id ].name;
  old_to_new[n] = bound.create_node( children, name )[0];
  return;
}

mockturtle::bound_network klut_to_bound( mockturtle::klut_network const& klut, std::vector<mockturtle::gate> const& gates )
{
  using namespace mockturtle;
  using TT = kitty::dynamic_truth_table;
  std::unordered_map<TT, uint32_t, kitty::hash<TT>> tt_to_id;
  std::unordered_map<klut_network::node, bound_network::signal> old_to_new;
  int i = 0;
  for( auto g : gates )
    tt_to_id[g.function] = i++;

  mockturtle::bound_network bound( gates );
  klut.foreach_pi( [&]( auto pi ){
    old_to_new[pi] = bound.create_pi();
  } );
  std::cout << bound.num_pis() << std::endl;

  klut.foreach_po( [&]( auto fo, auto io ){
    build_recursive( klut, bound, fo, old_to_new, tt_to_id, gates );
    bound.create_po( old_to_new[fo] );
  } );
  std::cout << bound.num_gates() << std::endl;


  return bound;

}

mockturtle::bound_network read_verilog( std::string const& path, std::vector<mockturtle::gate> const& gates )
{
  using namespace std;
  mockturtle::bound_network res( gates );

  std::ifstream file ( path );
  std::string line;

  enum state_t : uint32_t
  {
    start,
    collect_input,
    collect_output,
    collect_gates,
    convert,
    end
  };


  vector<string> pos;
  std::unordered_map<std::string, mockturtle::bound_network::signal> map;
  if (file.is_open()) 
  {
    state_t state = state_t::start;
    vector<string> signal_names;
    string name;
    bool collecting = true;
    while ( std::getline(file, line) )
    {
      std::cout << line << std::endl;
      istringstream iss( line );
      string s;
      while ( getline( iss, s, ' ' ) )
      {
        s.erase (std::remove(s.begin(), s.end(), ','), s.end());
        if( s == "endmodule" )
          break;
        if( s == "input" )
        {
//std::cout << "change state to input collection" << std::endl;
          state = collect_input;
          continue;
        }
        if( s == "output" )
        {
//std::cout << res.num_pis() << " number of PIs" << std::endl;
//std::cout << "change state to output collection" << std::endl;
          state = collect_output;
          continue;
        }
        if( state == collect_gates )
        {
          //if( collecting )
          {
            if( s[0] == 'g' )
            {
//std::cout << ".>" << s << std::endl;

              while( s[0] != '(' )
                s.erase( s.begin() );
              s.erase( s.begin() );
              s.erase( s.begin() + std::size( s ) - 1 );
            }
            if( s[0] == '.' )
            {
              while( s[0] != '(' )
                s.erase( s.begin() );
              s.erase( s.begin() );
              if( s[std::size(s)-1] == ';' )
              {
                state = convert;
                s.erase( s.begin() + std::size( s ) - 1 ); // - ;
                s.erase( s.begin() + std::size( s ) - 1 ); // - )
              }
              if( s[std::size(s)-1] == ')' )
                s.erase( s.begin() + std::size( s ) - 1 );
              if( s[std::size(s)-1] == ')' )
                s.erase( s.begin() + std::size( s ) - 1 );
//std::cout << "-->" << s << std::endl;
              signal_names.push_back( s );
            }
            if( state == convert )
            {
              //s.erase (std::remove(s.begin(), s.end(), ';'), s.end());
//std::cout << "#names " << signal_names.size() << std::endl;
              std::vector<mockturtle::bound_network::signal> children;
              int i = 0;
              while( i < ( signal_names.size() - 1 ) )
              {
                children.push_back( map[signal_names[i++]] );
              }
              map[signal_names[i]] = res.create_node( children, name )[0];
//std::cout << "out is " << signal_names[i] << "-- ->" << map[signal_names[i]] << std::endl;

            }
          }
        }
        else if( state == collect_input )
        {
          if( s[0] == 'x' )
          {
            if( s[std::size(s)-1] == ';' )
              s.erase( s.begin() + std::size( s ) - 1 );
            std::cout << s << std::endl;
            map[s] = res.create_pi();
          }
        }
        else if( state == collect_output )
        {
          if( s[0] == 'y' )
          {
            if( s[std::size(s)-1] == ';' )
              s.erase( s.begin() + std::size( s ) - 1 );
//std::cout << "output: " << s << std::endl;
            pos.push_back( s );
          }
        }
        for( auto g : gates )
        {
          if( g.name == s )
          {
            signal_names.clear();
            state = collect_gates;
            collecting = true;
//std::cout << "change to gates collection" << std::endl;
            name = s;
          }
        }
      }
    }
    file.close();
  }
  else
  {
    std::cerr << "Unable to open file!" << std::endl;
  }


  for( auto po_name : pos )
    res.create_po( map[po_name] );
  return res;
}

int main()
{
  using namespace experiments;
  using namespace mockturtle;//klut_to_bound

  experiment<std::string, double, double, double, double, double, double, double, double, float, bool> exp(
      "dac_001", "benchmark", "A(bef)", "A(aft)", "D(bef)", "D(aft)", "G(bef)", "G(aft)", "P(bef)", "P(aft)", "runtime", "cec" );


  /* library to map to technology */
  //fmt::print( "[i] processing technology library\n" );
  std::string library = "asap7";
  std::vector<gate> gates;
  std::ifstream in( cell_libraries_path( library ) );

  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  bindings binds( gates );
  /* get the database */
  static constexpr uint32_t num_vars = 4;
  static constexpr uint32_t capacity = 10;
  static constexpr mockturtle::choice_order order = power_area_delay;
  using database_t = mockturtle::homo_database<num_vars, capacity, order>;
  /*
  _base_
  _emap_area_
  _lazy_dap_
  */
  std::string path = "/home/costamag/POST_MAPPING/mockturtle/build/" + library + "_lazy_pad_" + std::to_string( num_vars ) + ".db";

  database_t database( binds, path );// = create_resyn_database<num_vars, capacity, order>();
//database.print();
  tech_library_params tps;
  tps.ignore_symmetries = false; // set to true to drastically speed-up mapping with minor delay increase
  tps.verbose = true;
  tech_library<9> tech_lib( gates, tps );

  /* run resynthesis */
  for ( auto const& benchmark : all_benchmarks( iscas ) ) // experiments::voter ))//
  {
    if( benchmark == "hyp" || benchmark == "div" || benchmark == "sqrt" )
      continue;
    fmt::print( "[i] processing {}\n", benchmark );

    names_view<aig_network> aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    mockturtle::xag_npn_resynthesis<aig_network, aig_network, mockturtle::xag_npn_db_kind::aig_complete> resyn;
    mockturtle::exact_library_params eps;
    eps.np_classification = false;
    mockturtle::exact_library<aig_network> exact_lib( resyn, eps );

    mockturtle::rewrite_params rwps;
    rwps.preserve_depth = true;
    std::cout << aig.num_gates() << " " << aig.num_pos() << std::endl;
    mockturtle::rewrite( aig, exact_lib, rwps );
    std::cout << aig.num_gates() << " " << aig.num_pos() << std::endl;
    /* remove structural redundancies */
    aig_balancing_params bps;
    bps.minimize_levels = true;
    bps.fast_mode = true;
    aig_balance( aig, bps );
    const uint32_t size_before = aig.num_gates();
    const uint32_t depth_before = depth_view( aig ).depth();

    emap_params ps;
    ps.matching_mode = emap_params::hybrid;
    ps.area_oriented_mapping = false;
    ps.map_multioutput = false;
    emap_stats st;
    cell_view<block_network> res = emap<9>( aig, tech_lib, ps, &st );

    auto ntk = convert_mapped_to_bound( res, gates );
//bound.print();

    //write_verilog_with_binding( bound, "glitch/" + benchmark + "_bef.v" );
    std::string const library_full_path = "/home/costamag/POST_MAPPING/mockturtle/experiments/cell_libraries/asap7.genlib";
    abc_baseline( ntk, library_full_path, benchmark );

  using namespace mockturtle;
  klut_network klut;
  std::string path = "glitch/" + benchmark + "_abc.blif";
  std::istringstream in( path );
  std::cout << path << std::endl;
  const auto resr = lorina::read_blif( path, blif_reader( klut ) );
  if( resr != lorina::return_code::success )
    std::cout << "failed opening" << std::endl;
  std::cout << "#pis   " << klut.num_pis() << std::endl;
  std::cout << "#gates " << klut.num_gates() << std::endl;

    auto bound = klut_to_bound( klut, gates );
//bound.print();
    //abc_power_eval( bound, library_full_path );

    uint32_t constexpr num_steps = 12;
    uint32_t constexpr cube_size = 12;
    uint32_t constexpr max_leaves = 10;

    using Ntk = bound_network;
    
    using dNtk = depth_view<bound_network>;
    dNtk dbound{bound};

    using node_simulator_t = mockturtle::node_simulator<Ntk, num_steps, cube_size>;

    std::vector<node_simulator_t> pis_sims;
    std::default_random_engine::result_type seed = 1;
    std::array<kitty::static_truth_table<cube_size>, 2> workload;
    int i = 0;
    bound.foreach_pi([&]( auto n ){
      kitty::create_random( workload[0], seed + i );
      kitty::create_random( workload[1], seed + i + 1 );
      i += 2;
      node_simulator_t sim( n, workload , -1, -0.5, 0 );
      pis_sims.push_back( sim );
    });

    double power_simple_1_zero_delay;
    double power_simple_1_dyn_before;
    double power_simple_1_pow_before;
    double power_simple_1_glitching ;
    double delay_bef ;
    double area_bef       ;
    double zero_delay_bef ;
    double dyn_switch_bef ;
    double glitching_bef  ;

    {
      dNtk dbound{bound};
      mockturtle::power_analyzer<Ntk, num_steps, cube_size> power_simple_1( bound );
      power_simple_1( pis_sims );
      power_simple_1_zero_delay = power_simple_1.get_zero_delay();
      power_simple_1_dyn_before = power_simple_1.get_dyn_switch();
      power_simple_1_pow_before = power_simple_1.get_dyn_power();
      power_simple_1_glitching = power_simple_1.get_dyn_switch() - power_simple_1_zero_delay;
      arrival_analyzer<Ntk> arrival1( bound );
      arrival1.compute_arrival_times();
      delay_bef = arrival1.worst_delay();
      area_bef = bound.area();
    }
    std::cout << "ZD(0) " << power_simple_1_zero_delay << std::endl;
    std::cout << "GL(0) " << power_simple_1_glitching << std::endl;
    std::cout << "DS(0) " << power_simple_1_dyn_before << std::endl;
    std::cout << "PW(0) " << power_simple_1_pow_before << std::endl;


    auto beg = std::chrono::high_resolution_clock::now();
    mockturtle::resynthesis_params rps;
    rps.preserve_depth = true;
    rps.dynamic_database = false;
    rps.try_struct = true;
    rps.try_window = true;
    rps.try_simula = true;
uint32_t map_size_before = bound.size();

std::cout << "Optimizing" << std::endl;
    mockturtle::delay_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves>( bound, database, rps );
    auto bound4 = mockturtle::cleanup_bound( bound, gates );
std::cout << "Done optimizing" << std::endl;
//abc_power_eval( bound3, library_full_path );

std::cout << "Optimizing" << std::endl;
    mockturtle::delay_resynthesize<Ntk, database_t, num_vars, num_steps, cube_size, max_leaves>( bound4, database, rps );
    auto bound5 = mockturtle::cleanup_bound( bound4, gates );
std::cout << "Done optimizing" << std::endl;

//abc_power_eval( bound4, library_full_path );

    auto end = std::chrono::high_resolution_clock::now();

    dNtk dbound2{bound5};

    double delay_aft ;
    double area_aft ;
    double power_simple_2_zero_delay ;
    double power_simple_2_dyn_before ;
    double power_simple_2_pow_before ;
    double power_simple_2_glitching ;

    auto result = mockturtle::cleanup_bound( bound5, gates );
   // abc_power_eval( result, library_full_path );

//result.print();
    {
      dNtk dresult{result};
      arrival_analyzer<Ntk> arrival2( result );
      arrival2.compute_arrival_times();
      delay_aft = arrival2.worst_delay();
      area_aft = result.area();

      mockturtle::power_analyzer<Ntk, num_steps, cube_size> power_simple_2( result );
      power_simple_2( pis_sims );
      power_simple_2_zero_delay = power_simple_2.get_zero_delay();
      power_simple_2_pow_before = power_simple_2.get_dyn_power();
      power_simple_2_dyn_before = power_simple_2.get_dyn_switch();
      power_simple_2_glitching = power_simple_2.get_dyn_switch() - power_simple_2_zero_delay;
    }
    std::cout << "ZD(2) " << power_simple_2_zero_delay << std::endl;
    std::cout << "GL(2) " << power_simple_2_glitching << std::endl;
    std::cout << "DS(2) " << power_simple_2_dyn_before << std::endl;
    std::cout << "PW(2) " << power_simple_2_pow_before << std::endl;
std::cout << map_size_before << std::endl;
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - beg);

//result.print();

    const auto cec = benchmark == "hyp" ? true : abc_cec( result, benchmark );

    /* write verilog netlist */
    //write_verilog_with_binding( result, "glitch/" + benchmark + "_aft.v" );
    exp( benchmark, area_bef, area_aft, delay_bef, delay_aft, power_simple_1_glitching, power_simple_2_glitching, power_simple_1_dyn_before, power_simple_2_dyn_before, duration.count()*0.001, cec );
  }

  exp.save();
  exp.table();

  return 0;
}