// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// NOLINTBEGIN(misc-include-cleaner)

#include <filesystem>
#include <string>
#include <valarray>

#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

// NOLINTEND(misc-include-cleaner)

int
main( int, char** )
{
    std::string const dirname = "lorenz_data";
    auto filename             = std::filesystem::path( dirname ) / "lorenz.dat";
    ponio::observer::file_observer fobs( filename );

    using state_t = std::valarray<double>;

    double const sigma = 10.;
    double const rho   = 28.;
    double const beta  = 8. / 3.;

    auto lorenz = [=]( double, state_t const& u, state_t& du )
    {
        du[0] = sigma * ( u[1] - u[0] );
        du[1] = rho * u[0] - u[1] - u[0] * u[2];
        du[2] = u[0] * u[1] - beta * u[2];
    };

    state_t const u0 = { 1., 1., 1. };

    ponio::time_span<double> const tspan = { 0., 20. };
    double const dt                      = 0.01;

    ponio::solve( lorenz, ponio::runge_kutta::rk_nssp_53(), u0, tspan, dt, fobs );

    return 0;
}
