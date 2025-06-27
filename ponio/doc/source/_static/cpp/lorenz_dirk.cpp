// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Eigen/Dense>

#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

#include <ponio/eigen_linear_algebra.hpp>

int
main()
{
    using namespace ponio::observer; // to use _fobs litteral
    using vector_type = Eigen::Vector<double, 3>;
    using matrix_type = Eigen::Matrix<double, 3, 3>;

    using state_t = vector_type;

    double sigma = 10., rho = 28., beta = 8. / 3.;

    auto f = [=]( double /* t */, auto const& u, state_t& du )
    {
        du[0] = sigma * ( u[1] - u[0] );
        du[1] = rho * u[0] - u[1] - u[0] * u[2];
        du[2] = u[0] * u[1] - beta * u[2];
    };
    auto jac_f = [=]( double, state_t const& u ) -> matrix_type
    {
        return matrix_type( {
            {-sigma,      sigma, 0    },
            { rho - u[2], -1,    -u[0]},
            { u[1],       u[0],  -beta}
        } );
    };
    auto lorenz = ponio::make_implicit_problem( f, jac_f );

    state_t const u0 = { 1., 1., 1. };

    ponio::time_span<double> const tspan = { 0., 20. };
    double const dt                      = 0.01;

    ponio::solve( lorenz, ponio::runge_kutta::dirk34(), u0, tspan, dt, "lorenz_dirk.txt"_fobs );

    return 0;
}
