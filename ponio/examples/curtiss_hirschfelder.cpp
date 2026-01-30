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
    std::string const dirname = "ch_data";

    using state_t = double;

    double const tf = 2.0;
    double const dt = 0.05;

    double const k = 50;

    auto ch_pb = ponio::implicit_problem(
        [=]( double t, double y, double& dy )
        {
            dy = k * ( std::cos( t ) - y );
        },
        [=]( double /* t */, double /* y */ )
        {
            return -k;
        } );

    state_t const y_0 = 2.0;

    { // simple use case with explicit method
        auto filename = std::filesystem::path( dirname ) / "sol_rk_33_ralston_cst.dat";
        auto obs      = ponio::observer::file_observer( filename );

        ponio::solve( ch_pb, ponio::runge_kutta::rk_33_ralston(), y_0, { 0., tf }, dt, obs );
    }

    { // simple use case with diagonal implicit method
        auto filename = std::filesystem::path( dirname ) / "sol_dirk23_cst.dat";
        auto obs      = ponio::observer::file_observer( filename );

        ponio::solve( ch_pb, ponio::runge_kutta::dirk23().newton_tol( 1e-5 ).newton_max_iter( 1000 ), y_0, { 0., tf }, dt, obs );
    }

    { // example of time loop with while loop controlled by user
        auto filename = std::filesystem::path( dirname ) / "sol_rk_33_ralston.dat";
        auto obs      = ponio::observer::file_observer( filename );

        auto sol_range = ponio::make_solver_range( ch_pb, ponio::runge_kutta::rk_33_ralston(), y_0, { 0., 0.464, tf }, dt );
        auto it_sol    = sol_range.begin();

        while ( it_sol->time < tf )
        {
            obs( it_sol->time, it_sol->state, it_sol->time_step );

            // pseudo adaptive time-step method
            if ( it_sol->time < 0.5 )
            {
                it_sol->time_step = 0.01;
                ++it_sol;
            }
            else
            {
                it_sol->time_step = 0.05;
                ++it_sol;
            }
        }
        obs( it_sol->time, it_sol->state, tf - it_sol->time ); // to save iteration where it_sol->time == tf
    }

    { // example of time loop with for loop on range
        auto filename = std::filesystem::path( dirname ) / "sol_rk54_6m.dat";
        auto obs      = ponio::observer::file_observer( filename );

        auto sol_range = ponio::make_solver_range( ch_pb, ponio::runge_kutta::rk54_6m(), y_0, { 0., tf }, dt );

        for ( auto ui : sol_range )
        {
            obs( ui.time, ui.state, ui.time_step );
        }
    }

    return 0;
}
