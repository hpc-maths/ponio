#include <iostream>
#include <valarray>
#include <numeric>
#include <random>
#include <sstream>

#include <solver/problem.hpp>
#include <solver/solver.hpp>
#include <solver/observer.hpp>
#include <solver/butcher_methods.hpp>

int
main (int argc,char** argv)
{
  std::size_t N = 10;
  if (argc > 1) { N = std::atoi(argv[1]); }

  using namespace observer;
  using state_t = std::valarray<double>;

  std::random_device rd;
  std::mt19937 gen(rd());

  std::normal_distribution<> d{0.,2};

  double dt = 1e-3;

  auto brownian_pb = ode::make_simple_problem(
    [&](double t, state_t const& y)->state_t {
      return {
        d(gen),
        d(gen)
      };
    }
  );

  state_t yini = { 0., 0. };

  for ( auto i=0u ; i<N ; ++i ) {
    std::stringstream filename; filename << "example5/brownian_"<< i << ".dat";
    ode::solve( brownian_pb , ode::butcher::rk_33<>() , yini , {0.,10.} , dt , observer::file_observer(filename.str()) );
  }

  return 0;
}
