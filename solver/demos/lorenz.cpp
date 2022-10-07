#include <iostream>
#include <valarray>
#include <numeric>
#include <vector>

#include <solver/problem.hpp>
#include <solver/solver.hpp>
#include <solver/observer.hpp>
#include <solver/butcher_methods.hpp>

int
main (int,char**)
{
  using namespace observer;
  using state_t = std::valarray<double>;

  double sigma=10. , rho=28., beta=8./3.;
  auto lorenz = [=](double t, state_t const& u) -> state_t
  {
    return {
      sigma*( u[1] - u[0] ),
      rho*u[0] - u[1] - u[0]*u[2],
      u[0]*u[1] - beta*u[2]
    };
  };

  state_t u0 = {1.,1.,1.};
  std::vector<double> tspan = {0.,20.};
  double dt = 0.01;

  ode::solve( lorenz , ode::butcher::rk_nssp_53<>() , u0 , tspan , dt , "example6/lorenz.dat"_fobs );

  return 0;
}
