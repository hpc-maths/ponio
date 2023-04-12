// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <filesystem>
#include <iostream>
#include <tuple>

#include <solver/butcher_methods.hpp>
#include <solver/observer.hpp>
#include <solver/solver.hpp>

// solve $\dot{u} = u$ with $u(t=0) = 1$, and $t\in[0,2]$.

int main(int, char**)
{
    std::string dirname = "exp_data";
    auto filename       = std::filesystem::path(dirname) / "exp.dat";
    observer::file_observer fobs(filename);

    auto identity = [](double t, double u)
    {
        return u;
    };
    double x0 = 1.0;
    double dt = 0.1;

    ode::solve(identity, ode::butcher::rk_nssp_21(), x0, {0., 2.0}, dt, fobs);

    return 0;
}
