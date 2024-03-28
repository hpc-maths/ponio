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

#include <ponio/observer.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

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
    operator()( double, std::valarray<double> const& y ) const
    {
        double const oneoverdxdx = 1. / ( m_dx * m_dx );
        std::size_t const nx     = y.size();

        std::valarray<double> ydot( nx );

        ydot[0] = oneoverdxdx * ( -2. * y[0] + y[1] );
        for ( std::size_t i = 1; i < nx - 1; ++i )
        {
            ydot[i] = oneoverdxdx * ( y[i - 1] - 2. * y[i] + y[i + 1] );
        }
        ydot[nx - 1] = oneoverdxdx * ( y[nx - 2] - 2. * y[nx - 1] );

        return ydot;
    }

    static std::valarray<double>
    fundamental_sol( double t, std::valarray<double> const& x )
    {
        double const pi = std::numbers::pi;
        return ( 1. / ( 2. * std::sqrt( pi * t ) ) ) * ( std::exp( -( x * x ) / ( 4. * t ) ) );
    }
};

void
save( std::valarray<double> const& x, std::valarray<double> const& y, std::filesystem::path const& filename )
{
    auto printer_xy = []( double a, double b )
    {
        std::stringstream ss;
        ss << a << " " << b;
        return ss.str();
    };

    std::ofstream of( filename );
    std::transform( std::begin( x ), std::end( x ), std::begin( y ), std::ostream_iterator<std::string>( of, "\n" ), printer_xy );
    of.close();
}

int
main()
{
    std::string const dirname = "heat_data";
    std::filesystem::create_directories( dirname );

    std::size_t const nx = 1000;

    double const xmin = -5;
    double const xmax = 5;

    double const dx = ( xmax - xmin ) / static_cast<double>( nx + 1 );
    double const dt = 10 * dx * dx;

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

    double const tini = 0.001;
    double const tend = 0.5;

    std::valarray<double> const yini = heat_model::fundamental_sol( tini, x );
    std::valarray<double> yend( 0., nx );
    ponio::time_span<double> const tspan = { tini, tend };

    save( x, yini, std::filesystem::path( dirname ) / "heat_ini.dat" );

    yend = ponio::solve( pb_heat, ponio::runge_kutta::explicit_rkc2<15>(), yini, tspan, dt, observer::null_observer() );

    std::valarray<double> const yexa = heat_model::fundamental_sol( tend, x );

    auto const err = std::abs( yexa - yend ).sum() / static_cast<double>( nx );
    std::cout << "L1 norm of error = " << err << std::endl;

    save( x, yend, std::filesystem::path( dirname ) / "heat_sol.dat" );
    save( x, yexa, std::filesystem::path( dirname ) / "heat_exa.dat" );

    return 0;
}
