#include <iostream>
#include <tuple>

#include <solver/problem.hpp>
#include <solver/solver.hpp>
#include <solver/observer.hpp>
#include <solver/butcher_methods.hpp>

// First example to solve

int
main (int,char**)
{
  using namespace observer;

  auto pb_exp = ode::make_simple_problem(
          []( double t, double u ){ return u; }
        );

  ode::solve( pb_exp , ode::butcher::rk_nssp_21() , 1.0 , {0.,1.15,2.0} , 0.1 , "exp.dat"_fobs );

  return 0;
}
