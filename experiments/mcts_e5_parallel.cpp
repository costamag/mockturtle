#include <kitty/constructors.hpp>
#include <mockturtle/algorithms/mcts/parallelize.hpp>
#include <string>
#include <bit>
#include <bitset>
#include <cstdint>
#include <fmt/format.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <chrono>

using namespace mockturtle;
using namespace mcts;

int main( int argc, char * argv [] )
{
  auto start = std::chrono::high_resolution_clock::now();
  test_parallelism(1);
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);  
  printf ("It took me %lld seconds.\n", duration.count());
  
  start = std::chrono::high_resolution_clock::now();
  test_parallelism(-1);
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);  
  printf ("It took me %lld seconds.\n", duration.count());
  return 0;
}