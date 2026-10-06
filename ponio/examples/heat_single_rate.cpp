// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <cmath>
#include <cstddef>
#include <iostream>
#include <numbers>
#include <valarray>

#include <ponio/observer.hpp>
#include <ponio/runge_kutta/rkc.hpp>
#include <ponio/runge_kutta/rock.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

class heat_model
{
    double m_dx;

  public:

    explicit heat_model( double dx )
        : m_dx( dx )
    {
    }

    void
    operator()( double, std::valarray<double> const& y, std::valarray<double>& ydot ) const
    {
        double const inv_dxdx = 1. / ( m_dx * m_dx );
        std::size_t const nx  = y.size();

        ydot[0] = inv_dxdx * ( -2. * y[0] + y[1] );
        for ( std::size_t i = 1; i < nx - 1; ++i )
        {
            ydot[i] = inv_dxdx * ( y[i - 1] - 2. * y[i] + y[i + 1] );
        }
        ydot[nx - 1] = inv_dxdx * ( y[nx - 2] - 2. * y[nx - 1] );
    }

    static std::valarray<double>
    fundamental_solution( double t, std::valarray<double> const& x )
    {
        double const pi = std::numbers::pi;
        return 1. / ( 2. * std::sqrt( pi * t ) ) * std::exp( -( x * x ) / ( 4. * t ) );
    }
};

double
l1_error( std::valarray<double> const& u, std::valarray<double> const& u_exact )
{
    return std::abs( u - u_exact ).sum() / static_cast<double>( u.size() );
}

int
main()
{
    std::size_t const nx = 1000;
    double const xmin    = -5.;
    double const xmax    = 5.;
    double const dx      = ( xmax - xmin ) / static_cast<double>( nx + 1 );
    double const dt      = 10. * dx * dx;

    std::valarray<double> x( nx );
    for ( std::size_t i = 0; i < nx; ++i )
    {
        x[i] = xmin + static_cast<double>( i + 1 ) * dx;
    }

    double const t_ini = 1.e-3;
    double const t_end = 5.e-1;
    double const rho   = 4. / ( dx * dx );

    auto problem                          = heat_model( dx );
    auto const u_ini                      = heat_model::fundamental_solution( t_ini, x );
    auto const u_exact                    = heat_model::fundamental_solution( t_end, x );
    ponio::time_span<double> const t_span = { t_ini, t_end };

    auto eigmax = [rho]( auto&&, double, auto&, double, auto& )
    {
        return rho;
    };

    auto const u_rkc = ponio::solve( problem, ponio::runge_kutta::chebyshev::rkc1( eigmax ), u_ini, t_span, dt, ponio::observer::null_observer() );

    auto const u_rock = ponio::solve( problem, ponio::runge_kutta::rock::rock2( eigmax ), u_ini, t_span, dt, ponio::observer::null_observer() );

    std::cout << "RKC1:  L1 norm of error = " << l1_error( u_rkc, u_exact ) << '\n';
    std::cout << "ROCK2: L1 norm of error = " << l1_error( u_rock, u_exact ) << '\n';

    return 0;
}
