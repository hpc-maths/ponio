// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

int
main()
{
    using namespace ponio::observer; // to use _fobs litteral
    using vector_type = Eigen::Vector<double, 3>;
    using matrix_type = Eigen::Matrix<double, 3, 3>;

    using state_t = vector_type;

    double sigma = 10., rho = 28., beta = 8. / 3.;

    auto L = matrix_type{
        {-sigma, sigma, 0    },
        { rho,   -1,    0    },
        { 0,     0,     -beta}
    };
    auto N = [=]( double, auto&& u, state_t& du )
    {
        du[0] = 0.;
        du[1] = -u[0] * u[2];
        du[2] = u[0] * u[1];
    };
    auto lorenz = ponio::make_lawson_problem( L, N );

    // define matrix exponential
    auto m_exp = []( matrix_type const& A ) -> matrix_type
    {
        return A.exp();
    };

    state_t const u0 = { 1., 1., 1. };

    ponio::time_span<double> const tspan = { 0., 20. };
    double const dt                      = 0.01;

    ponio::solve( lorenz, ponio::runge_kutta::lrk_44( m_exp ), u0, tspan, dt, "lorenz_lrk.txt"_fobs );

    return 0;
}
