#include <iostream>
#include <valarray>
#include <numeric>
#include <numbers>

#include <solver/problem.hpp>
#include <solver/solver.hpp>
#include <solver/observer.hpp>
#include <solver/butcher_methods.hpp>

/*
solve Arenstorf orbit problem :

*/

struct arenstorf_model
{
  using state_t = std::valarray<double>;

  double mu;

  arenstorf_model( double m )
  : mu(m)
  {}

  double
  D1 ( double x , double y )
  {
    return std::pow( (x + mu)*(x + mu) + y*y , 1.5 );
  }

  double
  D2 ( double x , double y )
  {
    return std::pow( (x-1.+mu)*(x-1.+mu) + y*y , 1.5 );
  }

  state_t
  operator () ( double t , state_t const& u )
  {
    double x  = u[0], y  = u[2];
    double xp = u[1], yp = u[3];
    return {
      xp ,
      x + 2.*yp - (1.-mu)*(x+mu)/D1(x,y) - mu*(x-1.+mu)/D2(x,y) ,
      yp ,
      y - 2.*xp - (1.-mu)*y/D1(x,y) - mu*y/D2(x,y)
    };
  }

};

int
main (int,char**)
{
  using namespace observer;
  using state_t = std::valarray<double>;

  double Tf = 17.0652165601579625588917206249;
  double dt = 1e-5;

  double mu = 0.012277471;

  auto arenstorf_pb = ode::make_problem( arenstorf_model(mu) );

  state_t uini = { 0.994 , 0. , 0. , -2.00158510637908252240537862224 };

  ode::solve( arenstorf_pb , ode::butcher::rk54_6m<>(1e-5) , uini , {0.,Tf} , dt , "example4/arenstorf_rk546m.dat"_fobs );
  
  ode::solve( arenstorf_pb , ode::butcher::rk54_7m<>(1e-5) , uini , {0.,Tf} , dt , "example4/arenstorf_rk547m.dat"_fobs );

  ode::solve( arenstorf_pb , ode::butcher::rk54_7s<>(1e-5) , uini , {0.,Tf} , dt , "example4/arenstorf_rk547s.dat"_fobs );

  return 0;
}
