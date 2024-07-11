// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numbers>
#include <sstream>
#include <string>
#include <valarray>

#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

struct nagumo
{
    double dx;
    double x_min;
    double x_max;
    double x_0;
    double k;
    double d;

    nagumo( double ndx, double nx_min, double nx_max, double nx_0, double nk, double nd )
        : dx( ndx )
        , x_min( nx_min )
        , x_max( nx_max )
        , x_0( nx_0 )
        , k( nk )
        , d( nd )
    {
    }

    std::valarray<double>
    operator()( double, std::valarray<double> const& u ) const
    {
        double const oneoverdxdx = 1. / ( dx * dx );
        std::size_t const nx     = u.size();

        std::valarray<double> udot( nx );

        udot[0] = oneoverdxdx * ( -2. * u[0] + 2. * u[1] );
        for ( std::size_t i = 1; i < nx - 1; ++i )
        {
            udot[i] = oneoverdxdx * ( u[i - 1] - 2. * u[i] + u[i + 1] );
        }
        udot[nx - 1] = oneoverdxdx * ( 2. * u[nx - 2] - 2. * u[nx - 1] );

        udot += k * u * u * ( 1.0 - u );

        return udot;
    }

    std::valarray<double>
    exact_solution( double t, std::valarray<double> const& x ) const
    {
        double const v   = ( 1. / std::numbers::sqrt2 ) * ( std::sqrt( k * d ) );
        double const cst = -( 1. / std::numbers::sqrt2 ) * std::sqrt( k / d );

        return std::exp( cst * ( x - x_0 - v * t ) ) / ( 1.0 + std::exp( cst * ( x - x_0 - v * t ) ) );
    }
};

void
save( std::valarray<double> const& x, std::valarray<double> const& y, std::filesystem::path const& dir, double t )
{
    auto printer_xy = []( double a, double b )
    {
        std::stringstream ss;
        ss << a << " " << b;
        return ss.str();
    };

    std::stringstream filename;
    filename << "u_" << std::setprecision( 3 ) << t << ".dat";

    std::ofstream of( dir / filename.str() );
    std::transform( std::begin( x ), std::end( x ), std::begin( y ), std::ostream_iterator<std::string>( of, "\n" ), printer_xy );
    of.close();
}

int
main()
{
    std::string const dirname = "nagumo_data";
    std::filesystem::create_directories( dirname );

    std::size_t const nx = 501;

    double const x_max = 50.0;
    double const x_min = -50.0;
    double const x_0   = -10.0;
    double const k     = 1.0;
    double const d     = 1.0;

    double const dx = ( x_max - x_min ) / static_cast<double>( nx - 1 );

    std::valarray<double> x( nx );
    std::generate( std::begin( x ),
        std::end( x ),
        [x_min, dx, count = 0]() mutable
        {
            return x_min + dx * ( count++ );
        } );

    auto pb = nagumo( dx, x_min, x_max, x_0, k, d );

    double const t_ini = 0.0;
    double const t_end = 50.0;
    double const dt    = ( t_end - t_ini ) / 100;

    ponio::time_span<double> const t_span = { t_ini, t_end };

    std::valarray<double> const u_n = pb.exact_solution( t_ini, x );

    auto sol_range = ponio::make_solver_range( pb, ponio::runge_kutta::rkc_202(), u_n, t_span, dt );
    auto it_sol    = sol_range.begin();

    std::size_t n_iteration = 0;
    while ( it_sol->time < t_end )
    {
        if ( n_iteration % 10 == 0 )
        {
            save( x, it_sol->state, std::filesystem::path( dirname ), it_sol->time );
        }

        ++n_iteration;
        ++it_sol;
    }
    save( x, it_sol->state, std::filesystem::path( dirname ), it_sol->time );

    return 0;
}
