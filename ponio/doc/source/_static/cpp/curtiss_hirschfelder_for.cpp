// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <ponio/observer.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

int
main()
{
    using namespace ponio::observer;

    double k = 50.;

    auto curtiss_hirschfelder = [=]( double t, double y, double& dy )
    {
        dy = k * ( std::cos( t ) - y );
    };

    ponio::time_span<double> const t_span = { 0., 2. }; // begin and end time
    double const dt                       = 0.01;       // time step

    double y_0 = 2.;

    auto obs = "curtiss_hirschfelder_for.txt"_fobs;

    auto sol_range = ponio::make_solver_range( curtiss_hirschfelder, ponio::runge_kutta::rk_33(), y_0, t_span, dt );

    for ( auto& ui : sol_range )
    {
        obs( ui.time, ui.state, ui.time_step );

        if ( ui.time < 0.5 )
        {
            ui.time_step = 0.01;
        }
        else
        {
            ui.time_step = 0.05;
        }
    }

    return 0;
}
