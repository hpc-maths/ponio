#include <iostream>
#include <tuple>

#include <solver/problem.hpp>
#include <solver/solver.hpp>
#include <solver/observer.hpp>
#include <solver/butcher_methods.hpp>

// solve $\dot{u} = u$ with $u(t=0) = 1$, and $t\in[0,2]$.

int
main (int,char**)
{
  using namespace observer;

  auto identity = []( double t, double u ){ return u; };
  double x0 = 1.0;
  double dt = 0.1;

  ode::solve( identity , ode::butcher::rk_nssp_21<>() , x0 , {0.,1.15,2.0} , dt , "example1/exp.dat"_fobs );

  return 0;
}
