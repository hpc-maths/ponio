// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <iostream>
#include <valarray>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <numbers>
#include <sstream>

#include <solver/generic_butcher_rk.hpp>
#include <solver/observer.hpp>
#include <solver/solver.hpp>
#include <solver/time_span.hpp>

// Heat
class heat_model
{
    double m_dx;

  public:

    heat_model( double dx )
        : m_dx( dx )
    {
    }

    std::valarray<double>
    operator()( double t, std::valarray<double> const& y )
    {
        double oneoverdxdx = 1. / ( m_dx * m_dx );
        std::size_t nx     = y.size();
        std::valarray<double> ydot( nx );
        ydot[0] = oneoverdxdx * ( -2. * y[0] + y[1] );
        for ( std::size_t i = 1; i < nx - 1; ++i )
        {
            ydot[i] = oneoverdxdx * ( y[i - 1] - 2. * y[i] + y[i + 1] );
        }
        ydot[nx - 1] = oneoverdxdx * ( y[nx - 2] - 2. * y[nx - 1] );
        return ydot;
    }

    std::valarray<double>
    fundamental_sol( double t, std::valarray<double> const& x )
    {
        double pi = std::numbers::pi;
        return ( 1. / ( 2. * std::sqrt( pi * t ) ) ) * ( std::exp( -( x * x ) / ( 4. * t ) ) );
    }
};

int
main( int argc, char** argv )
{
    std::string dirname = "heat_data";
    std::filesystem::create_directories( dirname );

    std::size_t nx = 1000;

    double xmin = -5;
    double xmax = 5;

    double dx = ( xmax - xmin ) / ( nx + 1 );
    double dt = 10 * dx * dx;

    std::valarray<double> x( nx );
    std::generate( std::begin( x ),
        std::end( x ),
        [xmin, dx, count = 0]() mutable
        {
            return xmin + dx * ( ++count );
        } );
    // std::copy(std::begin(x), std::end(x), std::ostream_iterator<double>(std::cout, "  "));
    // std::cout << std::endl;

    auto pb_heat = heat_model( dx );

    double tini                = 0.001;
    double tend                = 0.2;
    std::valarray<double> yini = pb_heat.fundamental_sol( tini, x );
    std::valarray<double> yend( 0., nx );
    ponio::time_span<double> tspan = { tini, tend };

    yend = ode::solve( pb_heat, ode::butcher::chebyshev::explicit_rkc2<28>(), yini, tspan, dt, observer::null_observer() );

    std::valarray<double> yexa = pb_heat.fundamental_sol( tend, x );

    auto err = std::abs( yexa - yend ).sum() / nx;
    std::cout << "L1 norm of error = " << err << std::endl;

    auto printer_xy = []( double x, double y )
    {
        std::stringstream ss;
        ss << x << " " << y;
        return ss.str();
    };

    {
        // save solution
        std::ofstream of( std::filesystem::path( dirname ) / "heat_sol.dat" );
        std::transform( std::begin( x ), std::end( x ), std::begin( yend ), std::ostream_iterator<std::string>( of, "\n" ), printer_xy );
        of.close();
    }

    {
        // save exact solution
        std::ofstream of( std::filesystem::path( dirname ) / "heat_exa.dat" );
        std::transform( std::begin( x ), std::end( x ), std::begin( yexa ), std::ostream_iterator<std::string>( of, "\n" ), printer_xy );
        of.close();
    }

    return 0;
}
