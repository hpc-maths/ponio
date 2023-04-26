// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <iostream>
#include <valarray>

#include <solver/butcher_methods.hpp>
#include <solver/observer.hpp>
#include <solver/problem.hpp>
#include <solver/solver.hpp>
#include <solver/time_span.hpp>

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

    std::valarray<double>
    operator()( double t, std::valarray<double> u )
    {
        double du1 = m_a - ( m_b + 1 ) * u[0] + u[0] * u[0] * u[1];
        double du2 = m_b * u[0] - u[0] * u[0] * u[1];
        return std::valarray<double>{ du1, du2 };
    }
};

int
main( int, char** )
{
    std::string dirname = "brusselator_data";
    auto filename       = std::filesystem::path( dirname ) / "brusselator.dat";
    observer::file_observer fobs( filename );

    auto pb_brusselator = ode::make_simple_problem( brusselator_model( 1., 3. ) );

    std::valarray<double> uini = { 1.5, 3 };

    ponio::time_span<double> tspan = { 0., 20.0 };

    double dt = 0.01;

    ode::solve( pb_brusselator, ode::butcher::rk_86(), uini, tspan, dt, fobs );

    return 0;
}
