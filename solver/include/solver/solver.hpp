#pragma once

#include <vector>
#include <iostream>

#include "method.hpp"

namespace ode {

  /**
   * @brief solve a problem on a specific time range with a specific method
   * 
   * @param pb     problem to solve, it could by any function or functor with a call operator with following parameter `(value_t tn, state_t const& un)`.
   * @param algo   choosen method to solve the problem `pb`
   * @param u0     initial condition \f$u_0 = u(t=0)\f$
   * @param t_span container \f$[t_\text{start} , t_\text{end})\f$ with possible intermediate time value where solver should go
   * @param dt     time step value \f$\Delta t\f$
   * @param obs    observer that do something with current time, solution and time step at each iteration
   * @return returns the last value of solution \f$u^n\f$
   */
  template < typename Problem_t , typename Algorithm_t , typename state_t , typename value_t , typename Observer_t >
  state_t
  solve ( Problem_t & pb , Algorithm_t && algo , state_t const& u0 , std::vector<value_t> const& t_span , value_t dt , Observer_t obs )
  {
    value_t current_time = t_span.front();
    auto it_next_time = t_span.begin();
    ++it_next_time;
    auto it_end_time = t_span.end();

    value_t current_dt = dt;
    bool reset_dt = false;

    state_t un  = u0;
    state_t un1 = u0;

    auto meth = make_method(algo ,un);

    obs( current_time , un , dt );

    while ( it_next_time != it_end_time ) {

      std::tie( current_time , un1 , current_dt ) = meth( pb , current_time , un , current_dt );
      std::swap(un,un1);

      //current_time += current_dt;

      obs( current_time , un , current_dt );

      // prepare next step
      if ( current_time+current_dt > *it_next_time ) {
        current_dt = *it_next_time - current_time;
        reset_dt = true;
        ++it_next_time;
      } else if ( reset_dt ) {
        reset_dt = false;
        current_dt = dt;
      }

    }

    return un;
  }

}
