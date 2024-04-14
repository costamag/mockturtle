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

#include <iostream>
#include <string>
#include <fstream>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/utils/spfd_utils.hpp>


#include <experiments.hpp>

using namespace experiments;
using namespace mockturtle;


int main()
{
  using namespace experiments;
  using namespace mockturtle;

  using TT = kitty::dynamic_truth_table;
  static constexpr uint32_t K = 7u;
  static constexpr uint32_t k = 6u;

  lut_resynthesis_t<k, K, TT> _resyn;
  
  using tt_hash = kitty::hash<TT>;
  std::unordered_set<TT, tt_hash> classes;
  TT tt(k);
  uint32_t i=0;
  do
  {
    const auto res = kitty::exact_npn_canonization( tt );
    classes.insert( std::get<0>( res ) );
    kitty::next_inplace( tt );
    
  } while ( i++<1000 && !kitty::is_const0( tt ) );

  TT ttb(k);
  TT ttf(k);
  int D{0};
  for( auto const& ttf : classes )
  {
    for( auto const& ttb : classes )
    {

        D++;

        klut_network klut;

        auto x1 = klut.create_pi();
        auto x2 = klut.create_pi();
        auto x3 = klut.create_pi();
        auto x4 = klut.create_pi();
        auto x5 = klut.create_pi();
        auto x6 = klut.create_pi();
        auto x7 = klut.create_pi();

        auto fb = klut.create_node( {x2,x3,x4,x5,x6,x7}, ttb );
        auto ff = klut.create_node( {x1, x2, x3,x4,x5, fb}, ttf );
        klut.create_po( ff );

        default_simulator<kitty::dynamic_truth_table> sim( K );
        const auto tt = simulate<kitty::dynamic_truth_table>( klut, sim )[0];

        auto res = _resyn.decompose( tt, 3u );
    //  kitty::print_binary(ttf);
    //  printf("<-F\n");
    //  kitty::print_binary(ttb);
    //  printf("<-B\n");

        if( !res )
        {
            printf("error\n");
        }
        
        if( !kitty::equal( tt, _resyn.sims[*res] ) && !kitty::equal( ~tt, _resyn.sims[*res] ) )
        {

          _resyn.print();
            printf("\n    ");
            kitty::print_binary(_resyn.sims.back());
            printf("\n    ");
            kitty::print_binary(tt);
            printf("%d\n    ",D);
            printf("\n");
          return 0;

        }
        if( _resyn.num_luts()>2 )
        {
            printf("%d\n",D);

            printf("%d\n",_resyn.num_luts() );
            return 0;
        }
        printf("%d\n",_resyn.num_luts());

        //return 0;

    }
  }


  return 0;
}


