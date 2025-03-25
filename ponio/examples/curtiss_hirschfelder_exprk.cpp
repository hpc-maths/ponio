// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// NOLINTBEGIN(misc-include-cleaner)

#include <cmath>
#include <filesystem>
#include <string>

#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

// NOLINTEND(misc-include-cleaner)

/*
solve Curtiss and Hirschfelder problem:

$$
    \begin{aligned}
        \dot{y} =  k(\cos(t) - y) \\
        y(0) = y_0
    \end{aligned}
$$

$y_0 = 2$
*/

int
main()
{
    std::string const dirname = "curtiss_hirschfelder_exprk_data";

    using state_t = double;

    double const tf = 2.0;
    double const dt = 0.05;

    double const k = 50;

    auto linear_part    = -k;
    auto nonlinear_part = [=]( double t, double )
    {
        return k * std::cos( t );
    };

    auto pb_curtiss_hirshfelder = ponio::make_lawson_problem( linear_part, nonlinear_part );

    state_t const y_0                    = 2.0;
    ponio::time_span<double> const tspan = { 0., tf };

    { // test with RK(4, 4) method
        auto filename = std::filesystem::path( dirname ) / "rk44.dat";
        auto obs      = ponio::observer::file_observer( filename );

        ponio::solve( pb_curtiss_hirshfelder, ponio::runge_kutta::rk_44(), y_0, tspan, dt, obs );
    }

    { // test with expRK Krogstad method
        auto filename = std::filesystem::path( dirname ) / "krogstad.dat";
        auto obs      = ponio::observer::file_observer( filename );

        ponio::solve( pb_curtiss_hirshfelder, ponio::runge_kutta::krogstad(), y_0, tspan, dt, obs );
    }

    { // test with LRK(4, 4) method
        auto filename = std::filesystem::path( dirname ) / "lrk44.dat";
        auto obs      = ponio::observer::file_observer( filename );

        auto exp = []( double x )
        {
            return std::exp( x );
        };

        ponio::solve( pb_curtiss_hirshfelder, ponio::runge_kutta::lrk_44( exp ), y_0, tspan, dt, obs );
    }

    return 0;
}
