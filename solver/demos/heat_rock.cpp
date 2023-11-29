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

#include <solver/observer.hpp>
#include <solver/runge_kutta.hpp>
#include <solver/solver.hpp>
#include <solver/time_span.hpp>

// Heat
class heat_model
{
    double m_dx;
    double m_xmin;
    double m_xmax;

  public:

    heat_model( double dx, double xmin, double xmax )
        : m_dx( dx )
        , m_xmin( xmin )
        , m_xmax( xmax )
    {
    }

    std::valarray<double>
    operator()( double, std::valarray<double> const& y ) const
    {
        double const oneoverdxdx = 1. / ( m_dx * m_dx );
        std::size_t const nx     = y.size();

        std::valarray<double> ydot( nx );

        ydot[0] = oneoverdxdx * ( -2. * y[0] + 1. * y[1] );
        for ( std::size_t i = 1; i < nx - 1; ++i )
        {
            ydot[i] = oneoverdxdx * ( y[i - 1] - 2. * y[i] + y[i + 1] );
        }
        ydot[nx - 1] = oneoverdxdx * ( 1. * y[nx - 2] - 2. * y[nx - 1] );

        return ydot;
    }

    std::valarray<double>
    fundamental_sol( double t, std::valarray<double> const& x )
    {
        double const xmid = 0.5 * ( m_xmax + m_xmin );
        double const pi   = std::numbers::pi;
        return ( 1. / ( 2. * std::sqrt( pi * t ) ) ) * ( std::exp( -( ( x - xmid ) * ( x - xmid ) ) / ( 4. * t ) ) );
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
    std::string const dirname = "heat_rock_data";
    std::filesystem::create_directories( dirname );

    std::size_t const nx = 101;

    double const xmin = -5.0;
    double const xmax = 5.0;

    double const dx = ( xmax - xmin ) / static_cast<double>( nx - 1 );

    std::valarray<double> x( nx );
    std::generate( std::begin( x ),
        std::end( x ),
        [xmin, dx, count = 0]() mutable
        {
            return xmin + dx * ( count++ );
        } );

    auto pb_heat = heat_model( dx, xmin, xmax );

    double const t_ini = 0.01;
    double const t_end = 0.011; // + 1e-5;

    std::valarray<double> const y_ini = pb_heat.fundamental_sol( t_ini, x );
    std::valarray<double> y_end( 0., nx );
    ponio::time_span<double> const tspan = { t_ini, t_end };

    save( x, y_ini, std::filesystem::path( dirname ) / "heat_ini.dat" );

    std::valarray<double> y_qexa = ponio::solve( pb_heat, ponio::runge_kutta::rkc_202(), y_ini, tspan, 1e-6, observer::null_observer() );
    ;
    // std::ifstream is( std::filesystem::path( dirname ) / "sol_qexa.txt" );
    // std::istream_iterator<double> start( is ), end;
    // std::copy( start, end, std::begin( y_qexa ) );
    // save( x, y_qexa, std::filesystem::path( dirname ) / "heat_qexa.dat" );

    // std::valarray<double> const y_exa = pb_heat.fundamental_sol( t_end, x );
    // save( x, y_exa, std::filesystem::path( dirname ) / "heat_exa.dat" );

    for ( std::size_t N = 1; N < 600; N *= 2 )
    {
        double dt = ( t_end - t_ini ) / static_cast<double>( N );
        // double dt = 1e-5;

        y_end = ponio::solve( pb_heat, ponio::runge_kutta::rock::rock2(), y_ini, tspan, dt, observer::null_observer() );
        // y_end = ponio::solve( pb_heat, ponio::runge_kutta::rkc_202(), y_ini, tspan, dt, observer::null_observer() );

        std::valarray<double> diff = std::abs( y_qexa - y_end );

        auto err = std::sqrt( std::accumulate( std::begin( diff ),
            std::end( diff ),
            0.0,
            [dx]( double partial_sum, double val )
            {
                return partial_sum + val * val * dx;
            } ) );
        std::cout << dt << " " << std::setprecision( 20 ) << err << std::endl;
    }

    save( x, y_end, std::filesystem::path( dirname ) / "heat_sol.dat" );

    return 0;
}
