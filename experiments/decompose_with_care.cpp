#include <fstream>
#include <iostream>
#include <string>

#include <experiments.hpp>

#include <kitty/dynamic_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <kitty/print.hpp>

#include <mockturtle/utils/spfd_utils.hpp>

using namespace mockturtle;

int main()
{
    std::ifstream input("../experiments/benchmarks/luts_6_4.txt");
    std::string line;
    std::vector<uint32_t> counter = {0,0,0,0,0,0,0,0,0,0,0,0,0};

    int cnt{0};
    kitty::dynamic_truth_table tt(6u);
    kitty::dynamic_truth_table mk(6u);
    uint32_t reward;
    for( std::string line; getline( input, line ); )
    {
        if( cnt == 1000 )
            break;
        if( cnt % 3 == 0 ) // reward
        {
            //std::cout << line << std::endl;
            reward = std::stoi( line );
        }
        else if( cnt % 3 == 1 ) // func
        {
            //std::cout << line << std::endl;
            kitty::create_from_binary_string( tt, line );
        }
        else if( cnt % 3 == 2 ) // care
        {
            //std::cout << line << std::endl;
            kitty::create_from_binary_string( mk, line );
            
            lut_resynthesis_t<4, 10u> resyn;
            auto lit_out = resyn.decompose( tt, mk, 20 );
            if( lit_out && resyn.num_luts() <= reward )
            {
                counter[reward - resyn.num_luts()]++;
            }
        }
        cnt++;

    }
    for( int i{0}; i<counter.size(); ++i )
        printf("[%2d %2d]", i, counter[i]);
    printf("\n");
    return 0;
}