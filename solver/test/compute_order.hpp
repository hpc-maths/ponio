#pragma once

#include <vector>
#include <tuple>
#include <numeric>
#include <cmath>
#include <valarray>

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

template <typename Algorithm_t, typename T=double>
auto
solve_brusselator ( T dt )
{
  using state_t = std::valarray<T>;

  auto pb = [](T t,state_t y)->state_t {
      double a = 1. , b=3.;

      return state_t{
        1. - (b+1.)*y[0] + a*y[0]*y[0]*y[1] ,
        b*y[0] - a*y[0]*y[0]*y[1]
      };
    };
  
  state_t y0 = { 1.5 , 3. };
  std::vector<T> t_span = {0.,1.0};

  #ifdef DEBUG
    std::stringstream ss; ss << std::string(Algorithm_t::id) << "/dt_" << dt << ".dat";
    std::string filename = ss.str();
    auto obs = observer::file_observer(filename);
  #else
    auto obs = [](T,state_t,T){};
  #endif
    return ode::solve( pb , Algorithm_t(1e-5) , y0 , t_span , dt , std::forward<decltype(obs)>(obs) );
}


template <typename T=double>
auto
error (std::valarray<T> const& u, std::valarray<T> const& v)
{
  T r = 0.;

  auto it_u = std::begin(u);
  for ( auto it_v=std::begin(v) ; it_v != std::end(v) ; ++it_v,++it_u )
  {
    r += std::abs( *it_u - *it_v );
  }

  return r/u.size();
}

template <typename Algorithm_t, typename T=double>
T
order ()
{
  using state_t = std::valarray<T>;

  std::vector<T> errors, dts;

  state_t u_exa = {1.96873, 1.38722}; // solution of Brusselator at time 1.0

  #ifdef DEBUG
    fs::create_directories(std::string(Algorithm_t::id));
    std::ofstream f( std::string(Algorithm_t::id) + "/errors.dat"s );
  #endif

  for ( auto n_iter : {500,200,100,50,25,10} ) {
    T dt = 2.0/static_cast<double>(n_iter);
    state_t u_sol = solve_brusselator<Algorithm_t>(dt);
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

