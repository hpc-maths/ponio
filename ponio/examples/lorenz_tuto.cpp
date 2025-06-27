// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// NOLINTBEGIN(misc-include-cleaner)

#include <filesystem>
#include <string>
#include <utility>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/splitting.hpp>
#include <ponio/time_span.hpp>

// NOLINTEND(misc-include-cleaner)

// Lorenz
/*
  solve Lorenz equations with multiple methods

  $$
    \begin{cases}
      \dot{x} &= \sigma (y - x) \\
      \dot{y} &= \rho x - y - xz \\
      \dot{z} &= xy - \beta z
    \end{cases}
  $$

  we note $u = (x, y, z)$ and we can rewrite this problem into :

  $$
    \dot{u} = Lu + N(t, u)
  $$

  with :

  $$
    L = \begin{pmatrix}
      -\sigma & \sigma & 0 \\
      \rho    & -1     & 0 \\
      0       & 0      & -\beta
    \end{pmatrix}, \qquad
    N : t,u \mapsto \begin{pmatrix}
      0 \\
      -xz \\
      xy
    \end{pmatrix}
  $$

  we also use that as a splitting with first function defined as $t,u \mapsto Lu$ and the second is $N$.
*/

int
main()
{
    using vector_type = Eigen::Vector<double, 3>;    // NOLINT(misc-include-cleaner)
    using matrix_type = Eigen::Matrix<double, 3, 3>; // NOLINT(misc-include-cleaner)

    std::string const dirname = "lorenz_tuto_data";

    auto filename_1 = std::filesystem::path( dirname ) / "rk44.dat";
    auto filename_2 = std::filesystem::path( dirname ) / "lrk44.dat";
    auto filename_3 = std::filesystem::path( dirname ) / "lie.dat";
    auto filename_4 = std::filesystem::path( dirname ) / "strang.dat";

    ponio::observer::file_observer fobs_1( filename_1 );
    ponio::observer::file_observer fobs_2( filename_2 );
    ponio::observer::file_observer fobs_3( filename_3 );
    ponio::observer::file_observer fobs_4( filename_4 );

    double const sigma = 10.;
    double const rho   = 28.;
    double const beta  = 8. / 3.;

    auto L = matrix_type{
        {-sigma, sigma, 0    },
        { rho,   -1,    0    },
        { 0,     0,     -beta}
    };
    auto N = []( double, vector_type const& u, vector_type& du )
    {
        du = { 0., -u[0] * u[2], u[0] * u[1] };
    };
    auto pb = ponio::make_lawson_problem( L, N );

    // redefine a problem for splitting method
    auto linear_part = [=]( double, vector_type const& u, vector_type& du )
    {
        du = L * u;
    };
    auto nonlinear_part = [=]( double t, vector_type const& u, vector_type& du )
    {
        N( t, u, du );
    };
    auto pb_multi = ponio::make_problem( linear_part, nonlinear_part );

    // define matrix exponential
    auto mexp = []( matrix_type const& A ) -> matrix_type
    {
        return A.exp();
    };

    auto lie    = ponio::splitting::make_lie_tuple( std::make_pair( ponio::runge_kutta::rk_44(), 0.005 ),
        std::make_pair( ponio::runge_kutta::rk_44(), 0.005 ) );
    auto strang = ponio::splitting::make_strang_tuple( std::make_pair( ponio::runge_kutta::rk_44(), 0.005 ),
        std::make_pair( ponio::runge_kutta::rk_44(), 0.005 ) );

    vector_type const u_ini = { 1., 1., 1. };

    ponio::time_span<double> const tspan = ponio::linspace( 0., 20., 2 );
    double const dt                      = 0.01;

    ponio::solve( pb, ponio::runge_kutta::rk_44(), u_ini, tspan, dt, fobs_1 );
    ponio::solve( pb, ponio::runge_kutta::lrk_44( mexp ), u_ini, tspan, dt, fobs_2 );
    ponio::solve( pb_multi, lie, u_ini, tspan, dt, fobs_3 );
    ponio::solve( pb_multi, strang, u_ini, tspan, dt, fobs_4 );

    return 0;
}
