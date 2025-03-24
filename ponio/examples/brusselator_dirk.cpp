// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// NOLINTBEGIN(misc-include-cleaner)

#include <filesystem>
#include <string>

#include <Eigen/Dense>

#include <ponio/eigen_linear_algebra.hpp>
#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

// NOLINTEND(misc-include-cleaner)

struct lin_alg_2_2
{
    using vector_type = Eigen::Vector<double, 2>;    // NOLINT(misc-include-cleaner)
    using matrix_type = Eigen::Matrix<double, 2, 2>; // NOLINT(misc-include-cleaner)

    static vector_type
    solver( matrix_type const& A, vector_type const& b )
    {
        double const det = A( 0, 0 ) * A( 1, 1 ) - A( 1, 0 ) * A( 0, 1 );

        return vector_type{ -A( 0, 1 ) / det * b[1] + A( 1, 1 ) / det * b[0], A( 0, 0 ) / det * b[1] - A( 1, 0 ) / det * b[0] };
    }
};

// Brusselator

class brusselator_model
{
    using vector_type = Eigen::Vector<double, 2>;    // NOLINT(misc-include-cleaner)
    using matrix_type = Eigen::Matrix<double, 2, 2>; // NOLINT(misc-include-cleaner)

    double m_a;
    double m_b;

  public:

    brusselator_model( double a, double b )
        : m_a( a )
        , m_b( b )
    {
    }

    vector_type
    operator()( double, vector_type const& u ) const
    {
        double const du1 = m_a - ( m_b + 1 ) * u[0] + u[0] * u[0] * u[1];
        double const du2 = m_b * u[0] - u[0] * u[0] * u[1];
        return vector_type{ du1, du2 };
    }

    matrix_type
    jacobian( double, vector_type const& u ) const // NOLINT(readability-convert-member-functions-to-static)
    {
        double const y1 = u[0];
        double const y2 = u[1];

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
    using vector_type = Eigen::Vector<double, 2>; // NOLINT(misc-include-cleaner)

    std::string const dirname = "brusselator_dirk_data";

    auto filename_1 = std::filesystem::path( dirname ) / "brusselator_dirk23.dat";
    auto filename_2 = std::filesystem::path( dirname ) / "brusselator_dirk23_exact_solver.dat";
    auto filename_3 = std::filesystem::path( dirname ) / "brusselator_rk33.dat";

    ponio::observer::file_observer fobs_1( filename_1 );
    ponio::observer::file_observer fobs_2( filename_2 );
    ponio::observer::file_observer fobs_3( filename_3 );

    auto model          = brusselator_model( 1., 3. );
    auto pb_brusselator = ponio::make_implicit_problem( model,
        [&]( double t, vector_type const& u )
        {
            return model.jacobian( t, u );
        } );

    vector_type const uini = { 1.5, 3 };

    ponio::time_span<double> const tspan = { 0., 20.0 };

    double const dt = 0.25;

    ponio::solve( pb_brusselator, ponio::runge_kutta::dirk23(), uini, tspan, dt, fobs_1 );
    ponio::solve( pb_brusselator, ponio::runge_kutta::dirk23<lin_alg_2_2>(), uini, tspan, dt, fobs_2 );
    ponio::solve( pb_brusselator, ponio::runge_kutta::rk_33(), uini, tspan, dt, fobs_3 );

    return 0;
}
