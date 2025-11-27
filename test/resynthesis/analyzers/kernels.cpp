#include <catch.hpp>

#include <kitty/kitty.hpp>
#include <math.h>

#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/resynthesis/aig_syn.hpp>
#include <mockturtle/resynthesis/networks/bindings.hpp>
#include <mockturtle/resynthesis/placement/kernels.hpp>

TEST_CASE( "Basic kernel functionalities", "[kernels]" )
{
  int static constexpr N = 7;
  mockturtle::gaussian_kernel_t<N> kernel;
  double sigma = 0.5;
  double delta = 0.8;
  double border = 20;
  kernel.set( delta, sigma, border );
  kernel.print_field();

  mockturtle::vector_t<double> m { 6, 6 };
  for ( int x = -N/2; x <= N/2; ++x )
  {
    for ( int y = -N/2; y < N/2; ++y )
    {
      mockturtle::vector_t<double> q { m[0] + x*delta, m[1] + y*delta };
      mockturtle::vector_t<double> r = q - m;
      double const den = 2 * M_PI * pow( sigma, 2 );
      double const exponent = -r.squared_norm2() / ( 2 * pow( sigma, 2 ) );
      double const measure = pow( delta, 2 );
      double const computed = exp( exponent ) * measure / den;
      double stored = kernel.kernel_field( q, m );
      std::cout << computed << "/" << stored << " ";
      CHECK( ( std::abs( stored - computed ) < 0.01 ) );
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  kernel.print_gradx();
  kernel.print_grady();

}
