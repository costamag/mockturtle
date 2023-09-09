#include <catch.hpp>

#include <algorithm>
#include <vector>
#include <ctime>

#include <kitty/constructors.hpp>
#include <mockturtle/algorithms/techaware/sym_synthesis.hpp>

using namespace mockturtle;
using namespace techaware;

TEST_CASE( "Input matching", "[techaware]" )
{
  techaware::TT F(3u);
  kitty::create_nth_var(F, 0);
  std::vector<uint32_t> T {1, 5, 2};

  sym_synthesis<xag_network> synt( F, T );

  std::vector<techaware::TT> xs;
  xag_network xag;
  std::vector<xag_network::signal> signals;
  for( uint32_t i{0}; i<3u; ++i )
  {
    signals.push_back(xag.create_pi());
    xs.emplace_back(3u);
    kitty::create_nth_var(xs[i], i);
  }
  synt.rewrite( &xag, signals );
  printf("%d\n", xag.num_gates());

}

TEST_CASE( "Majority of 3", "[techaware]" )
{
  techaware::TT F(3u);
  kitty::create_majority(F);
  std::vector<uint32_t> T {0,0,0};

  sym_synthesis<xag_network> synt( F, T );

  std::vector<techaware::TT> xs;
  xag_network xag;
  std::vector<xag_network::signal> signals;
  for( uint32_t i{0}; i<3u; ++i )
  {
    signals.push_back(xag.create_pi());
    xs.emplace_back(3u);
    kitty::create_nth_var(xs[i], i);
  }
  synt.rewrite( &xag, signals );
  printf("%d\n", xag.num_gates());
}

TEST_CASE( "topdec", "[techaware]" )
{
  techaware::TT F(3u);
  std::vector<techaware::TT> xs;
  xag_network xag;
  std::vector<xag_network::signal> signals;
  for( uint32_t i{0}; i<3u; ++i )
  {
    signals.push_back(xag.create_pi());
    xs.emplace_back(3u);
    kitty::create_nth_var(xs[i], i);
  }
  F = xs[0] & xs[1];
  std::vector<uint32_t> T {0,0,0};

  sym_synthesis<xag_network> synt( F, T );
  synt.rewrite( &xag, signals );
  printf("%d\n", xag.num_gates());
}

TEST_CASE( "Majority of 4", "[techaware]" )
{
  techaware::TT F(4u);
  kitty::create_majority(F);

    std::clock_t start;
    double duration;

    start = std::clock();

  std::vector<uint32_t> T {0,0,0,0};

  sym_synthesis<xag_network> synt( F, T );

  std::vector<techaware::TT> xs;
  xag_network xag;
  std::vector<xag_network::signal> signals;
  for( uint32_t i{0}; i<4u; ++i )
  {
    signals.push_back(xag.create_pi());
    xs.emplace_back(4u);
    kitty::create_nth_var(xs[i], i);
  }
  synt.rewrite( &xag, signals );
  printf("%d\n", xag.num_gates());

    /* Your algorithm here */

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    std::cout<<"printf: "<< duration <<'\n';
}