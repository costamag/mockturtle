#include <thread>
#include <mutex>
#include <algorithm>
#include <set>
#include <random>

#include <iostream>

namespace mockturtle{

namespace mcts{

#pragma region mutex
std::atomic<uint32_t> exp_id{0};
#pragma endregion

void thread_run()
{
    int id = exp_id++;
    while( id < 100 )
    {
        for(int i = 0; i < 10000000; ++i)
        {
        }
        id = exp_id++;
    }
}

void test_parallelism(int nTHREADS)
{
  const auto processor_count = nTHREADS < 0 ? std::thread::hardware_concurrency() : nTHREADS;

  /* update json method */
  exp_id.store( 0 );
  std::vector<std::thread> threads;

  /* generate threads */
  printf( "[i] Running on %d threads\n", processor_count );
  for ( auto i = 0u; i < processor_count; ++i )
  {
    threads.emplace_back( thread_run );//, ps_contest, run_only_one );
  }

  /* wait threads */
  for ( auto i = 0u; i < processor_count; ++i )
  {
    threads[i].join();
  }

  // exp_res.table( "best" );
}

}

}