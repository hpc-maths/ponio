// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <filesystem>
#include <string>

#include <ponio/observer.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>

// solve $\dot{y} = -y$ with $y(0) = 1$, and $t\in[0,2]$.

int
main()
{
    using namespace ponio::observer; // to use _fobs litteral

    auto f = []( double /* t */, double y, double& dy )
    {
        dy = -y;
    };

    double const y0 = 1.0;
    double const dt = 0.1;

    ponio::solve( f, ponio::runge_kutta::euler(), y0, { 0., 2.0 }, dt, "exp.txt"_fobs );

    return 0;
}
