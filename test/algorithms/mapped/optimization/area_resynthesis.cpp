#include <catch.hpp>

#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>

#include <mockturtle/algorithms/mapped/optimization/resynthesize.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/mapped/bound_network.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>

std::string const test_library = "GATE   and2    1.0 O=a*b;                 PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   or2     1.0 O=a+b;                 PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   xor2    0.5 O=a^b;                 PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   or3     1.0 O=a+b+c;               PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   and3    1.0 O=((a*b)*c);           PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   maj3    1.0 O=(a*b)+(b*c)+(a*c);   PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   fa      1.0 C=a*b+a*c+b*c;         PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   fa      1.0 S=a^b^c;               PIN * INV 1   999 1.0 0.0 1.0 0.0";

struct custom_area_rewire_params : mockturtle::default_resynthesis_params
{
  static constexpr bool try_rewire = true;
  static constexpr bool try_struct = false;
  static constexpr bool try_window = false;
  static constexpr bool try_simula = false;
  static constexpr bool dynamic_database = false;
  static constexpr int32_t odc_levels = 0;
  static constexpr uint32_t max_input_win = 8u;
  static constexpr uint32_t max_cuts_size = 6u;
};

TEST_CASE( "Area resynthesis via rewiring - single-output gate without don't cares", "[area_resynthesis]" )
{
  using Ntk = mockturtle::bound_network<mockturtle::bound::design_type_t::CELL_BASED, 2>;
  using signal = typename Ntk::signal;
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bound::augmented_library<mockturtle::bound::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  using Db = mockturtle::mapped_database<Ntk, MaxNumVars>;
  Db db( lib );

  Ntk ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const e = ntk.create_pi();
  auto const f1 = ntk.create_node( { c, d, e }, 4u );
  auto const f2 = ntk.create_node( { a, b }, 0u );
  auto const f3 = ntk.create_node( { c, d }, 0u );
  auto const f4 = ntk.create_node( { e, f3 }, 0u );
  auto const f5 = ntk.create_node( { f2, f4 }, 0u );

  ntk.create_po( f1 );
  ntk.create_po( f5 );

  using DNtk = mockturtle::depth_view<Ntk>;
  mockturtle::window_manager_stats st;
  DNtk dntk( ntk );
  custom_area_rewire_params ps;
  mockturtle::area_resynthesize<DNtk, Db, custom_area_rewire_params>( dntk, db, ps );
  CHECK( ntk.area() == 3 );
}

TEST_CASE( "Area resynthesis via rewiring - single-output gate with don't cares", "[area_resynthesis]" )
{
  using Ntk = mockturtle::bound_network<mockturtle::bound::design_type_t::CELL_BASED, 2>;
  using signal = typename Ntk::signal;
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bound::augmented_library<mockturtle::bound::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  using Db = mockturtle::mapped_database<Ntk, MaxNumVars>;
  Db db( lib );

  Ntk ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const f1 = ntk.create_node( { a, b }, 0u );
  auto const f2 = ntk.create_node( { c, d }, 1u );
  auto const f3 = ntk.create_node( { c, d }, 0u );
  auto const f4 = ntk.create_node( { c, d }, 2u );
  auto const f5 = ntk.create_node( { f1, f2 }, 0u );
  auto const f6 = ntk.create_node( { f3, f5 }, 1u );
  auto const f7 = ntk.create_node( { f3, f4 }, 0u );

  ntk.create_po( f6 );
  ntk.create_po( f7 );

  using DNtk = mockturtle::depth_view<Ntk>;
  mockturtle::window_manager_stats st;
  DNtk dntk( ntk );
  custom_area_rewire_params ps;
  ps.window_manager_ps.odc_levels = 3;
  mockturtle::area_resynthesize<DNtk, Db, custom_area_rewire_params>( dntk, db, ps );
  CHECK( ntk.area() == 5.5 );
}

TEST_CASE( "Area resynthesis via rewiring - multiple-output gate without don't cares", "[area_resynthesis]" )
{
  using Ntk = mockturtle::bound_network<mockturtle::bound::design_type_t::CELL_BASED, 2>;
  using signal = typename Ntk::signal;
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bound::augmented_library<mockturtle::bound::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  using Db = mockturtle::mapped_database<Ntk, MaxNumVars>;
  Db db( lib );

  Ntk ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const e = ntk.create_pi();
  auto const f1 = ntk.create_node( { c, d, e }, 4u );
  auto const f2 = ntk.create_node( { a, b }, 0u );
  auto const f3 = ntk.create_node( { c, d }, 0u );
  auto const f4 = ntk.create_node( { e, f3 }, 0u );
  auto const f5 = ntk.create_node( { c, f2, f4 }, { 6u, 7u } );

  ntk.create_po( f1 );
  ntk.create_po( { f5.index, 0 } );
  ntk.create_po( { f5.index, 1 } );

  using DNtk = mockturtle::depth_view<Ntk>;
  mockturtle::window_manager_stats st;
  DNtk dntk( ntk );
  custom_area_rewire_params ps;
  mockturtle::area_resynthesize<DNtk, Db, custom_area_rewire_params>( dntk, db, ps );
  CHECK( ntk.area() == 3 );
}

TEST_CASE( "Area resynthesis via rewiring - multiple-output gate with don't cares", "[area_resynthesis]" )
{
  using Ntk = mockturtle::bound_network<mockturtle::bound::design_type_t::CELL_BASED, 2>;
  using signal = typename Ntk::signal;
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  mockturtle::bound::augmented_library<mockturtle::bound::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  using Db = mockturtle::mapped_database<Ntk, MaxNumVars>;
  Db db( lib );

  Ntk ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const e = ntk.create_pi();
  auto const f1 = ntk.create_node( { a, b }, 0u );
  auto const f2 = ntk.create_node( { c, d }, 1u );
  auto const f3 = ntk.create_node( { c, d }, 0u );
  auto const f4 = ntk.create_node( { c, d }, 2u );
  auto const f5 = ntk.create_node( { e, f1, f2 }, { 6u, 7u } );
  auto const f6 = ntk.create_node( { { f5.index, 0 }, e }, 0u );
  auto const f7 = ntk.create_node( { { f5.index, 1 }, f4 }, 0u );
  auto const f8 = ntk.create_node( { f6, f7 }, 0u );
  auto const f9 = ntk.create_node( { f3, f8 }, 1u );

  ntk.create_po( f9 );

  using DNtk = mockturtle::depth_view<Ntk>;
  mockturtle::window_manager_stats st;
  DNtk dntk( ntk );
  custom_area_rewire_params ps;
  mockturtle::area_resynthesize<DNtk, Db, custom_area_rewire_params>( dntk, db, ps );
  CHECK( ntk.area() == 6.0 );
}