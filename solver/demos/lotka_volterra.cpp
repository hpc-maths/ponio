// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>
#include <valarray>

#include <solver/observer.hpp>
#include <solver/problem.hpp>
#include <solver/runge_kutta.hpp>
#include <solver/solver.hpp>
#include <solver/time_span.hpp>

/*
Lotka-Volterra system
---------------------

$$
  \begin{cases}
    \frac{\mathrm{d}x}{\mathrm{d}t} = \alpha x - \beta xy\\
    \frac{\mathrm{d}y}{\mathrm{d}t} = \delta xy - \gamma y\\
  \end{cases}
$$

with $\alpha = \frac{2}{3}$, $\beta = \frac{4}{3}$, and $\delta = \gamma = 1$. Initial condition $(x,y)(t=0) = (x_0,x_0)$ done by user.

This system is solved by RK(11,8) Runge-Kutta method with time step $\Delta t=0.1$ to the time $t\leq 15$.
 */

int
main( int argc, char* argv[] )
{
    // default filename
    std::string const dirname = "lv_data";
    auto filename             = std::filesystem::path( dirname ) / "lv.dat";

    using state_t = std::valarray<double>;

    double x0 = 1.0; // default x0
    if ( argc > 1 )
    {
        filename = argv[1];
        x0       = std::stod( argv[2] );
    }
    observer::file_observer fobs( filename );

    // parameters
    double const alpha = 2. / 3.;
    double const beta  = 4. / 3.;
    double const gamma = 1.;
    double const delta = 1.;

    auto lotka_volterra_pb = ponio::make_simple_problem( // define problem
        [=]( double, state_t const& u ) -> state_t
        {
            return { alpha * u[0] - beta * u[0] * u[1], delta * u[0] * u[1] - gamma * u[1] };
        } );

    ponio::time_span<double> const t_span = { 0., 15. }; // begin and end time

    double const dt  = 0.1;        // time step
    state_t const u0 = { x0, x0 }; // initial condition

    ponio::solve( lotka_volterra_pb, ponio::runge_kutta::rk_118(), u0, t_span, dt, fobs );

    return 0;
}
