#pragma once

#include <vector>
#include <tuple>
#include <numeric>
#include <cmath>
#include <valarray>

//#define DEBUG

#ifdef DEBUG

#include <sstream>
#include <filesystem>
namespace fs = std::filesystem;
#include <solver/observer.hpp>
#include <string>
#include <string_view>
using namespace std::string_literals;

#endif

#include <solver/solver.hpp>
#include <solver/detail.hpp>
#include <solver/problem.hpp>

template <typename Container>
auto
mayor_method ( Container const& x , Container const& y )
{
  using value_t = typename Container::value_type;

  auto x_mid = std::begin(x) + (std::end(x)-std::begin(x))/2;
  auto y_mid = std::begin(y) + (std::end(y)-std::begin(y))/2;
  
  value_t x1_avg = std::accumulate(x.begin(), x_mid, value_t{0.0})/(x_mid-x.begin());
  value_t x2_avg = std::accumulate(x_mid, x.end(),   value_t{0.0})/(x.end()-x_mid);
  value_t y1_avg = std::accumulate(y.begin(), y_mid, value_t{0.0})/(y_mid-y.begin());
  value_t y2_avg = std::accumulate(y_mid, y.end(),   value_t{0.0})/(y.end()-y_mid);
  
  value_t a = (y2_avg-y1_avg)/(x2_avg-x1_avg);
  value_t b = y1_avg - a*x1_avg;

  return std::make_tuple( a, b );
}

template <typename T=double>
auto
error (T u, T v)
{
  return std::abs(u-v);
}

template <typename T=double>
auto
relative_error (T u, T v)
{
  return std::abs((u-v)/u);
}

template <typename Algorithm_t, typename T=double>
auto
solve_exp (Algorithm_t & algo, T dt , T Tf )
{
  using state_t = T;

  auto pb = [](T t, state_t y)->state_t {
    return y;
  };

  state_t y0 = 1.0;
  std::vector<T> t_span = {0.,Tf};

  #ifdef DEBUG
    std::stringstream ss; ss << "debug_info/" << std::string(Algorithm_t::id) << "/dt_" << dt << ".dat";
    std::string filename = ss.str();
    auto obs = observer::file_observer(filename);
  #else
    auto obs = [](T,state_t,T){};
  #endif
    return ode::solve( pb , algo , y0 , t_span , dt , obs );
}


template <typename Algorithm_t, typename T=double>
T
short_time_check_order (Algorithm_t algo=Algorithm_t())
{
  using state_t = T;

  std::vector<T> errors, dts;

  T Tf = 1.0;

  state_t u_exa = std::exp(Tf); // solution of Brusselator at time 1.0

  #ifdef DEBUG
    fs::create_directories("debug_info/" + std::string(Algorithm_t::id));
    std::ofstream f("debug_info/" + std::string(Algorithm_t::id) + "/errors.dat"s );
  #endif

  for ( auto n_iter : {50,25,20,15,10} ) {
    T dt = Tf/n_iter;
    state_t u_sol = solve_exp(algo,dt,Tf);
    auto e = error( u_exa , u_sol );
    errors.push_back(std::log( e ));
    dts.push_back(std::log(dt));

    #ifdef DEBUG
      f << dt << " " << e << "\n";
    #endif
  }

  #ifdef DEBUG
    f.close();
  #endif

  auto [a,b] = mayor_method(dts,errors);

  return a;
}

template <typename Algorithm_t, typename T=double>
T
long_time_check_order (Algorithm_t & algo)
{
  using state_t = std::valarray<T>;

  T alpha=2./3., beta=4./3., gamma=1., delta=1.;

  // make a problem that can be use also for splitting (in 3 parts) methods
  auto pb = ode::make_problem(
    [=]( T t , state_t const& u ) -> state_t {
      return {
        alpha*u[0] - beta*u[0]*u[1] ,
        0.
      };
    },
    [=]( T t , state_t const& u ) -> state_t {
      return {
        0. ,
        delta*u[0]*u[1]
      };
    },
    [=]( T t , state_t const& u ) -> state_t {
      return {
        0. ,
        - gamma*u[1]
      };
    }
  );

  // invariant calculator
  auto V = [=]( state_t const& u ) -> T {
    return delta*u[0] - gamma*std::log(u[0]) + beta*u[1] - alpha*std::log(u[1]);
  };

  T x0 = 1.9;
  state_t u_ini = {x0, x0};
  T V_ini = V(u_ini);

  std::vector<T> t_span = {0.,1000.};
  std::vector<T> dts = {0.25,0.125,0.1,0.075,0.05}; // find a way to adapt this range to the method
  std::vector<T> relative_errors, ldts;

  for ( auto dt : dts ) {
    state_t u_end = ode::solve(pb, algo, u_ini, t_span, dt,  [](T,state_t const&,T){});
    relative_errors.push_back( std::log10(relative_error(V_ini,V(u_end))) );
    ldts.push_back(std::log10(dt));
  }

  auto [a,b] = mayor_method(ldts,relative_errors);

  return a;
}

template <typename Algorithm_t, typename T=double>
T
check_order (Algorithm_t algo=Algorithm_t())
{
  if constexpr ( ode::butcher::is_embedded<Algorithm_t> ) {
    return Algorithm_t::order;
  } else if constexpr ( Algorithm_t::order >= 8 || ode::splitting::is_splitting_method<Algorithm_t> ) {
    return long_time_check_order(algo);
  } else {
    return short_time_check_order(algo);
  }
}
