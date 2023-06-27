#include <catch.hpp>

#include <mockturtle/algorithms/mcts/mnist_manager.hpp>
#include <kitty/partial_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <kitty/print.hpp>


using namespace mockturtle;
using namespace mcts;


TEST_CASE( "binary reader [0,1,2,3,4] [5,6,7,8,9]", "[MNIST]" )
{
    using PTT = kitty::partial_truth_table;
    std::vector<PTT> x_train = read_mnist_image_bin( "../experiments/MNIST/train-images.idx3-ubyte" );
    std::vector<PTT> y_train = read_mnist_label_04_59( "../experiments/MNIST/train-labels.idx1-ubyte" );
    std::vector<PTT> x_test = read_mnist_image_bin( "../experiments/MNIST/t10k-images.idx3-ubyte" );
    std::vector<PTT> y_test  = read_mnist_label_04_59( "../experiments/MNIST/t10k-labels.idx1-ubyte" );
    
    PTT x0_train_inv(28*28);
kitty::create_from_binary_string(x0_train_inv, 
"0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000011111111111100000000000011111111111111110000000000011111111111111110000000000001111111111100000000000000000011111110110000000000000000000101110000000000000000000000000111100000000000000000000000011110000000000000000000000000111110000000000000000000000001111110000000000000000000000011111100000000000000000000000111110000000000000000000000000111100000000000000000000011111110000000000000000000111111110000000000000000001111111110000000000000000011111111110000000000000000111111111100000000000000001111111111000000000000000000111111110000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
PTT x0_train;
for( auto i{28*28-1}; i>=0; --i )
    x0_train.add_bit( kitty::get_bit( x0_train_inv, i ) );

PTT x0;
for( int i{0}; i<x_train.size(); ++i )
    x0.add_bit( kitty::get_bit( x_train[i], 0 ) );
CHECK( kitty::equal(x0_train, x0 ));
CHECK( kitty::get_bit(y_train[0],0) == 1 );
CHECK( kitty::get_bit(y_train[0],1) == 0 );
CHECK( kitty::get_bit(y_train[0],2) == 0 );
CHECK( kitty::get_bit(y_train[0],3) == 0 );
CHECK( kitty::get_bit(y_train[0],4) == 1 );
CHECK( kitty::get_bit(y_train[0],5) == 0 );
CHECK( kitty::get_bit(y_train[0],6) == 0 );
CHECK( kitty::get_bit(y_train[0],7) == 0 );
CHECK( kitty::get_bit(y_train[0],8) == 0 );
CHECK( kitty::get_bit(y_train[0],9) == 0 );
CHECK( kitty::get_bit(y_test[0],0) == 1 );
CHECK( kitty::get_bit(y_test[0],1) == 0 );
CHECK( kitty::get_bit(y_test[0],2) == 0 );
CHECK( kitty::get_bit(y_test[0],3) == 0 );
CHECK( kitty::get_bit(y_test[0],4) == 0 );
CHECK( kitty::get_bit(y_test[0],5) == 0 );
CHECK( kitty::get_bit(y_test[0],6) == 0 );
CHECK( kitty::get_bit(y_test[0],7) == 1 );
CHECK( kitty::get_bit(y_test[0],8) == 1 );
CHECK( kitty::get_bit(y_test[0],9) == 1 );

}

TEST_CASE( "binary reader 10", "[MNIST]" )
{
    using PTT = kitty::partial_truth_table;
    std::vector<PTT> x_train = read_mnist_image_bin( "../experiments/MNIST/train-images.idx3-ubyte" );
    std::vector<PTT> y_train = read_mnist_label_10( "../experiments/MNIST/train-labels.idx1-ubyte" );
    std::vector<PTT> x_test = read_mnist_image_bin( "../experiments/MNIST/t10k-images.idx3-ubyte" );
    std::vector<PTT> y_test  = read_mnist_label_10( "../experiments/MNIST/t10k-labels.idx1-ubyte" );
    
/* 5 */
CHECK( kitty::get_bit(y_train[0],0) == 0 );
CHECK( kitty::get_bit(y_train[1],0) == 0 );
CHECK( kitty::get_bit(y_train[2],0) == 1 );
CHECK( kitty::get_bit(y_train[3],0) == 0 );
CHECK( kitty::get_bit(y_train[4],0) == 1 );
/* 0 */
CHECK( kitty::get_bit(y_train[0],1) == 1 );
CHECK( kitty::get_bit(y_train[1],1) == 0 );
CHECK( kitty::get_bit(y_train[2],1) == 0 );
CHECK( kitty::get_bit(y_train[3],1) == 0 );
CHECK( kitty::get_bit(y_train[4],1) == 0 );
/* 4 */
CHECK( kitty::get_bit(y_train[0],2) == 0 );
CHECK( kitty::get_bit(y_train[1],2) == 0 );
CHECK( kitty::get_bit(y_train[2],2) == 1 );
CHECK( kitty::get_bit(y_train[3],2) == 0 );
CHECK( kitty::get_bit(y_train[4],2) == 0 );
/* 1 */
CHECK( kitty::get_bit(y_train[0],3) == 1 );
CHECK( kitty::get_bit(y_train[1],3) == 0 );
CHECK( kitty::get_bit(y_train[2],3) == 0 );
CHECK( kitty::get_bit(y_train[3],3) == 0 );
CHECK( kitty::get_bit(y_train[4],3) == 1 );
/* 9 */
CHECK( kitty::get_bit(y_train[0],4) == 1 );
CHECK( kitty::get_bit(y_train[1],4) == 1 );
CHECK( kitty::get_bit(y_train[2],4) == 0 );
CHECK( kitty::get_bit(y_train[3],4) == 0 );
CHECK( kitty::get_bit(y_train[4],4) == 1 );
/* 2 */
CHECK( kitty::get_bit(y_train[0],5) == 0 );
CHECK( kitty::get_bit(y_train[1],5) == 0 );
CHECK( kitty::get_bit(y_train[2],5) == 0 );
CHECK( kitty::get_bit(y_train[3],5) == 1 );
CHECK( kitty::get_bit(y_train[4],5) == 0 );
/* 3 */
CHECK( kitty::get_bit(y_train[0],7) == 0 );
CHECK( kitty::get_bit(y_train[1],7) == 0 );
CHECK( kitty::get_bit(y_train[2],7) == 0 );
CHECK( kitty::get_bit(y_train[3],7) == 1 );
CHECK( kitty::get_bit(y_train[4],7) == 1 );
/* 7 */
CHECK( kitty::get_bit(y_test[0],0) == 0 );
CHECK( kitty::get_bit(y_test[1],0) == 0 );
CHECK( kitty::get_bit(y_test[2],0) == 1 );
CHECK( kitty::get_bit(y_test[3],0) == 1 );
CHECK( kitty::get_bit(y_test[4],0) == 1 );

}