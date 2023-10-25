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
    operator()( double, std::valarray<double> const& u )
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
    exact_solution( double t, std::valarray<double> const& x )
    {
        double v   = ( 1. / std::sqrt( 2.0 ) ) * ( std::sqrt( k * d ) );
        double cst = -( 1. / std::sqrt( 2.0 ) ) * std::sqrt( k / d );

        return std::exp( cst * ( x - x_0 - v * t ) ) / ( 1.0 + std::exp( cst * ( x - x_0 - v * t ) ) );
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
    std::string const dirname = "nagumo_data";
    std::filesystem::create_directories( dirname );

    std::size_t const nx = 2001;

    double const x_max = 50.0;
    double const x_min = -50.0;
    double const x_0   = -10.0;
    double k           = 1.0;
    double d           = 1.0;

    double const dx = ( x_max - x_min ) / static_cast<double>( nx - 1 );

    std::valarray<double> x( nx );
    std::generate( std::begin( x ),
        std::end( x ),
        [x_min, dx, count = 0]() mutable
        {
            return x_min + dx * ( count++ );
        } );

    auto pb = nagumo( dx, x_min, x_max, x_0, k, d );

    double t_ini                          = 0.0;
    double t_end                          = 10.0;
    ponio::time_span<double> const t_span = { t_ini, t_end };

    std::valarray<double> u_ini = pb.exact_solution( t_ini, x );
    std::valarray<double> u_end( 0., nx );

    save( x, u_ini, std::filesystem::path( dirname ) / "ini.dat" );

    std::valarray<double> u_qexa( nx );
    std::ifstream is( std::filesystem::path( dirname ) / "qexa.txt" );
    std::istream_iterator<double> start( is ), end;
    std::copy( start, end, std::begin( u_qexa ) );
    // save( x, u_qexa, std::filesystem::path( dirname ) / "qexa.dat" );

    for ( std::size_t N = 10; N < 20000; N *= 2 )
    {
        // double dt = 0.01;
        double dt = ( t_end - t_ini ) / static_cast<double>( N );
        u_end     = ponio::solve( pb, ponio::runge_kutta::rock::rock2(), u_ini, t_span, dt, observer::null_observer() );
        // u_end = ponio::solve( pb, ponio::runge_kutta::rkc_202(), u_ini, t_span, dt, observer::null_observer() );

        std::valarray<double> diff = std::abs( u_qexa - u_end );
        // save( x, diff, std::filesystem::path( dirname ) / "diff.dat" );
        double err = std::sqrt( std::accumulate( std::begin( diff ),
            std::end( diff ),
            0.0,
            [dx]( double partial_sum, double val ) -> double
            {
                return partial_sum + val * val * dx;
            } ) );
        std::cout << dt << " " << std::setprecision( 20 ) << err << std::endl;
    }

    save( x, u_end, std::filesystem::path( dirname ) / "end.dat" );

    return 0;
}
