// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <cmath>
#include <iostream>
#include <numbers>
#include <numeric>

#include <solver/butcher_methods.hpp>
#include <solver/observer.hpp>
#include <solver/problem.hpp>
#include <solver/solver.hpp>

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
main( int, char** )
{
    std::string dirname = "ch_data";
    std::string filename;

    using state_t = double;

    double tf = 2.0;
    double dt = 0.01;

    double k = 50;

    auto ch_pb = [=]( double t, double y )
    {
        return k * ( std::cos( t ) - y );
    };

    state_t y_0 = 2.0;

    { // example of time loop with while loop controlled by user
        filename = std::filesystem::path( dirname ) / "sol_rk_33_ralston.dat";
        auto obs = observer::file_observer( filename );

        auto sol_range = ode::make_solver_range( ch_pb, ode::butcher::rk_33_ralston(), y_0, { 0., tf }, dt );
        auto it_sol    = sol_range.begin();

        while ( it_sol->time < tf )
        {
            obs( it_sol->time, it_sol->state, it_sol->time_step );

            if ( it_sol->time < 0.5 )
            {
                ++it_sol;
            }
            else
            {
                it_sol += 0.05;
            }
        }
        obs( it_sol->time, it_sol->state, it_sol->time_step ); // to save iteration where it_sol->time == tf
    }

    { // example of time loop with for loop on range
        filename = std::filesystem::path( dirname ) / "sol_rk54_6m.dat";
        auto obs = observer::file_observer( filename );

        auto sol_range = ode::make_solver_range( ch_pb, ode::butcher::rk54_6m(), y_0, { 0., tf }, dt );

        for ( auto ui : sol_range )
        {
            obs( ui.time, ui.state, ui.time_step );
        }
    }

    return 0;
}
