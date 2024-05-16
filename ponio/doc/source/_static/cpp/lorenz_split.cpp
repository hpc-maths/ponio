// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <ponio/observer.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/splitting.hpp>
#include <ponio/time_span.hpp>

int
main()
{
    using namespace observer; // to use _fobs litteral
    using state_t = std::valarray<double>;

    double sigma = 10., rho = 28., beta = 8. / 3.;

    auto phi_0 = [=]( double /* t */, state_t& u ) -> state_t
    {
        double dt_u0 = sigma * u[1];
        double dt_u1 = rho * u[0];
        double dt_u2 = u[0] * u[1];

        return { dt_u0, dt_u1, dt_u2 };
    };
    auto phi_1 = [=]( double /* t */, state_t& u ) -> state_t
    {
        double dt_u0 = -sigma * -u[0];
        double dt_u1 = -u[1];
        double dt_u2 = -beta * u[2];

        return { dt_u0, dt_u1, dt_u2 };
    };
    auto phi_2 = [=]( double /* t */, state_t& u ) -> state_t
    {
        double dt_u0 = 0;
        double dt_u1 = -u[0] * u[2];
        double dt_u2 = 0;

        return { dt_u0, dt_u1, dt_u2 };
    };
    auto lorenz = ponio::make_problem( phi_0, phi_1, phi_2 );

    auto strang = ponio::splitting::make_strang_tuple( std::make_pair( ponio::runge_kutta::rk_44_38(), 0.01 ),
        std::make_pair( ponio::runge_kutta::rk_44(), 0.005 ),
        std::make_pair( ponio::runge_kutta::rk_44(), 0.0005 ) );

    state_t const u0 = { 1., 1., 1. };

    ponio::time_span<double> const tspan = { 0., 20. };
    double const dt                      = 0.01;

    ponio::solve( lorenz, strang, u0, tspan, dt, "lorenz_rk44.txt"_fobs );

    return 0;
}
