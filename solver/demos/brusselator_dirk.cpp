// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <functional>
#include <iostream>

#include <Eigen/Dense>

#include <solver/butcher_methods.hpp>
#include <solver/eigen_linear_algebra.hpp>
#include <solver/observer.hpp>
#include <solver/problem.hpp>
#include <solver/solver.hpp>
#include <solver/time_span.hpp>

struct lin_alg_2_2
{
    using vector_type = Eigen::Vector<double, 2>;
    using matrix_type = Eigen::Matrix<double, 2, 2>;

    vector_type
    solver( matrix_type const& dfx, vector_type const& fx )
    {
        double det = dfx( 0, 0 ) * dfx( 1, 1 ) - dfx( 1, 0 ) * dfx( 0, 1 );
        return vector_type{ -dfx( 0, 1 ) / det * fx[1] + dfx( 1, 1 ) / det * fx[0], dfx( 0, 0 ) / det * fx[1] - dfx( 1, 0 ) / det * fx[0] };
    }
};

// Brusselator

class brusselator_model
{
    using vector_type = Eigen::Vector<double, 2>;
    using matrix_type = Eigen::Matrix<double, 2, 2>;

    double m_a;
    double m_b;

  public:

    brusselator_model( double a, double b )
        : m_a( a )
        , m_b( b )
    {
    }

    vector_type
    operator()( double t, vector_type const& u )
    {
        double du1 = m_a - ( m_b + 1 ) * u[0] + u[0] * u[0] * u[1];
        double du2 = m_b * u[0] - u[0] * u[0] * u[1];
        return vector_type{ du1, du2 };
    }

    matrix_type
    jacobian( double t, vector_type const& u )
    {
        double y1 = u[0], y2 = u[1];
        return matrix_type{
            {2.0 * y1 * y2 - ( m_b + 1.0 ), y1 * y1 },
            { -2.0 * y1 * y2 + m_b,         -y1 * y1}
        };
    }
};

int
main( int, char** )
{
    using namespace std::placeholders;
    using vector_type = Eigen::Vector<double, 2>;
    using matrix_type = Eigen::Matrix<double, 2, 2>;

    std::string dirname = "brusselator_dirk_data";
    auto filename_1     = std::filesystem::path( dirname ) / "brusselator_dirk23.dat";
    auto filename_2     = std::filesystem::path( dirname ) / "brusselator_dirk23_exact_solver.dat";
    auto filename_3     = std::filesystem::path( dirname ) / "brusselator_rk33.dat";
    observer::file_observer fobs_1( filename_1 );
    observer::file_observer fobs_2( filename_2 );
    observer::file_observer fobs_3( filename_3 );

    auto model          = brusselator_model( 1., 3. );
    auto pb_brusselator = ode::make_implicit_problem( model, std::bind( &brusselator_model::jacobian, &model, _1, _2 ) );

    vector_type uini = { 1.5, 3 };

    ponio::time_span<double> tspan = { 0., 20.0 };

    double dt = 0.25;

    ode::solve( pb_brusselator, ode::butcher::dirk23(), uini, tspan, dt, fobs_1 );
    ode::solve( pb_brusselator, ode::butcher::dirk23<lin_alg_2_2>(), uini, tspan, dt, fobs_2 );
    ode::solve( pb_brusselator, ode::butcher::rk_33(), uini, tspan, dt, fobs_3 );

    return 0;
}
