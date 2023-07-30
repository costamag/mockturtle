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

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <mockturtle/algorithms/balancing.hpp>
#include <mockturtle/algorithms/balancing/sop_balancing.hpp>
#include <mockturtle/algorithms/balancing/mct_balancing.hpp>
#include <mockturtle/algorithms/collapse_mapped.hpp>
#include <mockturtle/algorithms/equivalence_checking.hpp>
#include <mockturtle/algorithms/lut_mapping.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/algorithms/miter.hpp>
#include <mockturtle/views/mapping_view.hpp>
#include <string>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, double, bool, uint32_t, uint32_t, double, bool> exp( "mct_balancing", "benchmark", "size", "depth", "size 4", "depth 4", "RT 4", "cec 4", "size 6", "depth 6", "RT 6", "cec 6" );

  mct_rebalancing<aig_network> mct_balancing;
  sop_rebalancing<aig_network> sop_balancing;

  for ( auto const& benchmark : epfl_benchmarks( ~experiments::hyp ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    depth_view daig{ aig };

    balancing_params ps;
    balancing_stats st4_mct, st4_sop;

    ps.progress = true;
    ps.only_on_critical_path = true;
    ps.cut_enumeration_ps.cut_size = 4u;

    lut_mapping_stats st_0;
    mapping_view<aig_network, true> mapped_aig_0{ aig };
    lut_mapping<decltype( mapped_aig_0 ), true>( mapped_aig_0, {}, &st_0 );
    const auto klut_0 = *collapse_mapped_network<klut_network>( mapped_aig_0 );
    equivalence_checking_stats ecst_0;
    depth_view dklut_0{ klut_0 };

    const auto aig4_sop = balancing( aig, { sop_balancing }, ps, &st4_sop );
    depth_view daig_sop{ aig4_sop };
    lut_mapping_stats st_sop;
    mapping_view<aig_network, true> mapped_aig_sop{ aig4_sop };
    lut_mapping<decltype( mapped_aig_sop ), true>( mapped_aig_sop, {}, &st_sop );
    const auto klut_sop = *collapse_mapped_network<klut_network>( mapped_aig_sop );
    equivalence_checking_stats ecst_sop;
    auto cec_sop = *equivalence_checking( *miter<klut_network>( aig, klut_sop ), {}, &ecst_sop );
    depth_view dklut_sop{ klut_sop };

    const auto aig4_mct = balancing( aig, { mct_balancing }, ps, &st4_mct );
    depth_view daig_mct{ aig4_mct };
    lut_mapping_stats st_mct;
    
    mapping_view<aig_network, true> mapped_aig_mct{ aig4_mct };
    lut_mapping<decltype( mapped_aig_mct ), true>( mapped_aig_mct, {}, &st_mct );
    const auto klut_mct = *collapse_mapped_network<klut_network>( mapped_aig_mct );
    depth_view dklut_mct{ klut_mct };

    equivalence_checking_stats ecst_mct;
    auto cec_mct = *equivalence_checking( *miter<klut_network>( aig, klut_mct ), {}, &ecst_mct );

    printf( "[%d %d][%d %d][%d %d]\n",
         klut_0.num_gates(), dklut_0.depth(),
         klut_sop.num_gates(), dklut_sop.depth(),
         klut_mct.num_gates(), dklut_mct.depth()
         );

    exp( benchmark,
         klut_0.num_gates(), dklut_0.depth(),
         klut_sop.num_gates(), dklut_sop.depth(),
         to_seconds( st4_sop.time_total ), cec_sop,
         klut_mct.num_gates(), dklut_mct.depth(),
         to_seconds( st4_mct.time_total ), cec_mct );

  }

  exp.save();
  exp.table();

  return 0;
}