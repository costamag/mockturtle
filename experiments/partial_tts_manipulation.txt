#include <kitty/partial_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <kitty/statistics.hpp>
#include <kitty/print.hpp>
int main()
{
    kitty::partial_truth_table tt( 4 );
    kitty::partial_truth_table x( 4 );
    kitty::partial_truth_table y( 4 );

    std::cout << "create random" << std::endl;
    kitty::create_from_binary_string( tt, "1100" );
    kitty::print_binary(tt);   
    std::cout << std::endl;
    std::cout << "get bit" << std::endl;
    for( auto i=0; i<tt.num_bits(); ++i )
        std::cout << get_bit( tt, i ) << std::endl;
    auto P = probability( tt );
    std::cout << P[0] << std::endl;
    std::cout << P[1] << std::endl;


    return 0;
}