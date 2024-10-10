// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <filesystem>
#include <string>

#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/splitting.hpp>
#include <ponio/user_defined_method.hpp>

// solve $\dot{u} = u$ with $u(t=0) = 1$, and $t\in[0,2]$.

int
main()
{
    std::string const dirname = "exp_splitting_data";
    auto filename             = std::filesystem::path( dirname ) / "exp_splitting.dat";
    observer::file_observer fobs( filename );

    double lambda = 0.3;

    auto f1 = [=]( double, double u )
    {
        return lambda * u;
    };

    auto f2 = [=]( double, double u )
    {
        return ( 1.0 - lambda ) * u;
    };

    auto pb = ponio::make_problem( f1, f2 );

    auto exact_solver_f1 = [=]( auto /* f */, double tn, double yn, double dt ) -> std::tuple<double, double, double>
    {
        return { tn + dt, std::exp( lambda * dt ) * yn, dt };
    };

    double const y0 = 1.0;
    double const dt = 0.1;

    auto strang = ponio::splitting::make_strang_tuple( std::make_pair( ponio::make_user_defined_method( exact_solver_f1 ), 0.005 ),
        std::make_pair( ponio::runge_kutta::rk_44(), 0.005 ) );

    ponio::solve( pb, strang, y0, { 0., 2.0 }, dt, fobs );

    return 0;
}
