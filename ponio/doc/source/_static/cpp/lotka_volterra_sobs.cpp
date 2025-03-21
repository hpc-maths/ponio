// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <iostream>
#include <sstream>
#include <valarray>

#include <ponio/observer.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

int
main()
{
    using state_t = std::valarray<double>;

    // parameters
    double alpha = 2. / 3.;
    double beta  = 4. / 3.;
    double gamma = 1.;
    double delta = 1.;

    auto lotka_volterra_pb = [=]( double, state_t const& u ) -> state_t
    {
        double dt_x = alpha * u[0] - beta * u[0] * u[1];
        double dt_y = delta * u[0] * u[1] - gamma * u[1];

        return { dt_x, dt_y };
    };

    std::stringstream buffer;
    auto obs = ponio::observer::stream_observer( buffer );

    ponio::time_span<double> const t_span = { 0., 15. }; // begin and end time
    double const dt                       = 0.1;         // time step

    state_t const u0 = { 1., 1. }; // initial condition

    ponio::solve( lotka_volterra_pb, ponio::runge_kutta::rk_33(), u0, t_span, dt, obs );

    std::cout << buffer.str() << "\n";

    return 0;
}
