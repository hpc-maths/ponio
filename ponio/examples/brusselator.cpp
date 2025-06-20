// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// NOLINTBEGIN(misc-include-cleaner)

#include <filesystem>
#include <string>
#include <valarray>

#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

// NOLINTEND(misc-include-cleaner)

// Brusselator

class brusselator_model
{
    double m_a;
    double m_b;

  public:

    brusselator_model( double a, double b )
        : m_a( a )
        , m_b( b )
    {
    }

    void
    operator()( double, std::valarray<double>&& u, std::valarray<double>& du ) const
    {
        du[0] = m_a - ( m_b + 1 ) * u[0] + u[0] * u[0] * u[1];
        du[1] = m_b * u[0] - u[0] * u[0] * u[1];
    }
};

int
main( int, char** )
{
    std::string const dirname = "brusselator_data";
    auto filename             = std::filesystem::path( dirname ) / "brusselator.dat";
    ponio::observer::file_observer fobs( filename );

    auto pb_brusselator = brusselator_model( 1., 3. );

    std::valarray<double> const uini = { 1.5, 3 };

    ponio::time_span<double> const tspan = { 0., 20.0 };

    double const dt = 0.01;

    ponio::solve( pb_brusselator, ponio::runge_kutta::rk_86(), uini, tspan, dt, fobs );

    return 0;
}
