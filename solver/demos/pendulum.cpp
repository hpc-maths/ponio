// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <filesystem>
#include <iostream>
#include <numbers>
#include <numeric>
#include <valarray>

#include <solver/butcher_methods.hpp>
#include <solver/observer.hpp>
#include <solver/problem.hpp>
#include <solver/solver.hpp>

/*
solve pendulum problem :
$$
  \ddot{\theta} + b\dot{theta} + c\sin(\theta) = 0
$$
By defining the angular velocity $\omega = \dot{theta}$, we obtain the system:
$$
  \begin{cases}
    \dot{theta} = \omega
    \dot{\omega} = -b\omega - c \sin(\theta)
  \end{cases}
$$

*/

int
main( int, char** )
{
    std::string dirname = "pendulum_data";
    auto filename       = std::filesystem::path( dirname ) / "pendulum.dat";
    observer::file_observer fobs( filename );

    using state_t = std::valarray<double>;

    double dt = 0.1;

    double b = 0.25, c = 5.0;

    auto pendulum_pb = ode::make_simple_problem(
        [=]( double, state_t const& y ) -> state_t
        {
            double theta = y[0], omega = y[1];
            return { omega, -b * omega - c * std::sin( theta ) };
        } );

    state_t yini = { std::numbers::pi - 0.1, 0. };

    ode::solve( pendulum_pb, ode::butcher::rk_44(), yini, { 0., 10.0 }, dt, fobs );

    return 0;
}
