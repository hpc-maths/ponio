// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <valarray>

#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

int
main()
{
    using namespace ponio::observer; // to use _fobs litteral
    using state_t = std::valarray<double>;

    double sigma = 10., rho = 28., beta = 8. / 3.;

    auto lorenz = ponio::make_simple_problem(
        [=]( double /* t */, state_t const& u ) -> state_t
        {
            double dt_u0 = sigma * ( u[1] - u[0] );
            double dt_u1 = rho * u[0] - u[1] - u[0] * u[2];
            double dt_u2 = u[0] * u[1] - beta * u[2];

            return { dt_u0, dt_u1, dt_u2 };
        } );

    state_t const u0 = { 1., 1., 1. };

    ponio::time_span<double> const tspan = { 0., 20. };
    double const dt                      = 0.01;

    ponio::solve( lorenz, ponio::runge_kutta::rk_44(), u0, tspan, dt, "lorenz_rk_pb.txt"_fobs );

    return 0;
}
