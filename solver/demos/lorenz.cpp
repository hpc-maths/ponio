// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <filesystem>
#include <iostream>
#include <numeric>
#include <valarray>

#include <solver/butcher_methods.hpp>
#include <solver/observer.hpp>
#include <solver/solver.hpp>
#include <solver/time_span.hpp>

int
main( int, char** )
{
    std::string const dirname = "lorenz_data";
    auto filename             = std::filesystem::path( dirname ) / "lorenz.dat";
    observer::file_observer fobs( filename );

    using state_t = std::valarray<double>;

    double const sigma = 10.;
    double const rho   = 28.;
    double const beta  = 8. / 3.;

    auto lorenz = [=]( double, state_t const& u ) -> state_t
    {
        auto du1 = sigma * ( u[1] - u[0] );
        auto du2 = rho * u[0] - u[1] - u[0] * u[2];
        auto du3 = u[0] * u[1] - beta * u[2];
        return { du1, du2, du3 };
    };

    state_t const u0 = { 1., 1., 1. };

    ponio::time_span<double> const tspan = { 0., 20. };
    double const dt                      = 0.01;

    ode::solve( lorenz, ode::butcher::rk_nssp_53(), u0, tspan, dt, fobs );

    return 0;
}
