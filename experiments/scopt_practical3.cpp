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
#include <kitty/constructors.hpp>
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
  static constexpr uint32_t K = 7u;
  TT tt(K);

  std::string line;
  
  double nSucc=0;
  double nIter=0;

  double nSucc1=0;
  double nSucc2=0;
  double nSucc3=0;
  double nIter1=0;
  double nIterDC=0;

  std::ifstream practical ("../experiments/NPN_practical/NPN_practical/"+std::to_string(K)+".txt");

  auto start = high_resolution_clock::now();
  
  if (practical.is_open())
  {
    while ( std::getline (practical,line) )
    {
      std::cout << nIter << std::endl ;

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
            //return 1;
          }
          else
          {
            nSucc += 1;
          }
        }
        else
        {
          printf("%d\n", _resyn.num_luts() );
          for( int i{0}; i<10; ++i )
          {
            nIterDC++;
            TT mk(K);
            kitty::create_random(mk);
            //kitty::print_binary(mk);printf("\n");

            /* synthesize with dcs */
            {
              auto resDC = _resyn.decompose( tt, mk, 3u );

              if( resDC )
              {
                //printf(" DC.%d ", _resyn.num_luts());
                if( !kitty::equal( tt&mk, _resyn.sims.back()&mk ) )
                {
                  printf("dc mistake\n");
                  //_resyn.print();
                  //return 1;
                }
                else if( _resyn.num_luts() <= 2 )
                {
                  nSucc1 += 1;
                }
              }
            }

            /* synthesize with set to 0 */
            {
              auto resF0 = _resyn.decompose( tt&mk, mk | ~mk, 3u );
              if( resF0 )
              {
                //printf(" F0.%d ", _resyn.num_luts());
                if( !kitty::equal( tt&mk, _resyn.sims.back()&mk ) )
                {
                  printf("f0 mistake\n");

                  //_resyn.print();
                  //return 1;
                }
                else if( _resyn.num_luts() <= 2 )
                {
                  nSucc2 += 1;
                }
              }
            }

            /* synthesize with set to random */
            {
              TT rd(K);
              kitty::create_random(rd);
              auto resRD = _resyn.decompose( (tt&mk) | (rd&(~mk)), mk | ~mk, 3u );
              if( resRD )
              {

                //printf(" F0.%d ", _resyn.num_luts());
                if( !kitty::equal( tt&mk, _resyn.sims.back()&mk ) )
                {
                  printf("rd mistake\n");
                  //_resyn.print();
                  //return 1;
                }
                else if ( _resyn.num_luts() <= 2 )
                {
                  nSucc3 += 1;
                }
              }
            }
          }
          std::cout << " DC " << nSucc1/nIterDC ;
          std::cout << " F0 " << nSucc2/nIterDC ;
          std::cout << " RD " << nSucc3/nIterDC ;
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
  std::cout << " DC " << nSucc1/nIterDC << std::endl;
  std::cout << " F0 " << nSucc2/nIterDC << std::endl;
  std::cout << " RD " << nSucc3/nIterDC << std::endl;
  std::cout << duration << " seconds" << std::endl;
  std::cout << nIterDC << " iters" << std::endl;

  return 0;
}



