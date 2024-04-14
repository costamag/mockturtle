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

#include <experiments.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/static_truth_table.hpp>
#include <kitty/print.hpp>
#include <mockturtle/utils/spfd_utils.hpp>
#include <chrono>
using namespace std::chrono;

using namespace experiments;
using namespace mockturtle;


int main()
{
  using namespace experiments;
  using namespace mockturtle;

  using TT = kitty::dynamic_truth_table;
  static constexpr uint32_t K = 5u;
  TT tt(K);

  std::string line;
  
  double nSucc=0;
  double nIter=0;
  std::ifstream practical ("../experiments/NPN_practical/NPN_practical/"+std::to_string(K)+".txt");

  auto start = high_resolution_clock::now();
  
  if (practical.is_open())
  {
    while ( std::getline (practical,line) )
    {
      lut_resynthesis_t<4, 11> _resyn;

      kitty::create_from_hex_string( tt, line );
      //kitty::print_binary(tt);
      auto res = _resyn.decompose( tt, 3u );
      if( res )
      {
        if( _resyn.num_luts() <= 2 )
        {
          if( !kitty::equal( tt, _resyn.sims.back() ) )
          {
            _resyn.print();
          }
          else
          {
            nSucc += 1;
          }
        }
        else
        {
          //printf("%d\n", _resyn.num_luts() );
        }
      }
      nIter+=1;

      //std::cout << nIter << " " << nSucc << std::endl;
    }
    practical.close();
  }
  auto stop = high_resolution_clock::now();
  auto duration = (duration_cast<milliseconds>(stop - start)).count()*0.001;
  std::cout << nSucc/nIter << std::endl;
  std::cout << duration << " seconds" << std::endl;


  


  return 0;
}



