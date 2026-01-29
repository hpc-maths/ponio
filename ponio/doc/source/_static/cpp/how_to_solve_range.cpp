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

    auto f = []( double /* t */, double y, double& dy )
    {
        dy = -y;
    };

    double const y0 = 1.0;
    double const dt = 0.1;
    auto obs        = "how_to_solve_range.txt"_fobs;

    auto sol_range = ponio::make_solver_range( f, ponio::runge_kutta::euler(), y0, { 0., 2.0 }, dt );
    auto it_sol    = sol_range.begin();

    while ( it_sol->time < 2.0 )
    {
        obs( it_sol->time, it_sol->state, it_sol->time_step );
        ++it_sol;
    }

    return 0;
}
