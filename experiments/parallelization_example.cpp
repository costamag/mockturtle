#include <thread>
#include <mutex>
#include <string>
#include <vector>
#include <iostream>
#include <atomic>


#pragma region mutex
std::atomic<uint32_t> exp_id{0};
std::mutex exp_mutex;
#pragma endregion


void thread_run( std::string const& run_only_one )
{

  uint32_t id = exp_id++;

  while ( id < 100 )
  {
    /* read benchmark */
    std::string benchmark =  "ex" + std::to_string( id );
    if ( run_only_one != "" && benchmark != run_only_one )
    {
      id = exp_id++;
      continue;
    }
    std::cout << "[i] processing " << benchmark << "\n";
    for( uint32_t i{0u}; i<1000000000; ++i )
      {}
    std::cout << "[i] done " << benchmark << std::endl;
    id = exp_id++;
  }
      
    
}


int main( int argc, char* argv[] )
{

  std::string run_only_one = "";

  if ( argc == 2 )
    run_only_one = std::string( argv[1] );

  const auto processor_count = run_only_one != "" ? 1 : std::thread::hardware_concurrency();

  /* starting benchmark id */
  exp_id.store( 0 );

  std::vector<std::thread> threads;

  /* generate threads */
  std::cout << "[i] Running on " << processor_count <<  "threads\n" <<std::endl;
  for ( auto i = 0u; i < processor_count; ++i )
  {
    threads.emplace_back( thread_run, run_only_one );
  }

  /* wait threads */
  for ( auto i = 0u; i < processor_count; ++i )
  {
    threads[i].join();
  }


  return 0;
}
