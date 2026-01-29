// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <filesystem>
#include <string>

#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>

// solve $\dot{y} = -y$ with $y(0) = 1$, and $t\in[0,2]$.

int
main()
{
    using namespace ponio::observer; // to use _fobs litteral

    auto pb = ponio::make_implicit_problem(
        []( double /* t */, double y, double& dy )
        {
            dy = -y;
        },
        []( double /* t */, double /* y */ )
        {
            return -1.;
        } );

    double const y0 = 1.0;
    double const dt = 0.1;

    {
        auto meth = ponio::runge_kutta::rk54_7s().abs_tol( 1e-4 ).rel_tol( 1e-5 );
        ponio::solve( pb, meth, y0, { 0., 2.0 }, dt, "how_to_solve_exp.txt"_fobs );
    }
    {
        auto meth = ponio::runge_kutta::backward_euler().newton_tol( 1e-3 ).newton_max_iter( 1000 );
        ponio::solve( pb, meth, y0, { 0., 2.0 }, dt, "how_to_solve_exp.txt"_fobs );
    }

    return 0;
}
