#include <algorithm>
#include <fstream>
#include <numbers>
#include <numeric>
#include <vector>

#include "ascent/Ascent.h"

int
main()
{
    std::size_t n_x = 100;
    double t0       = 0.;
    double tf       = 0.3;
    double dt       = 0.01;

    std::vector<double> x( n_x );
    std::generate( std::begin( x ),
        std::end( x ),
        [i = 0, n_x]() mutable
        {
            return static_cast<double>( i++ ) / static_cast<double>( n_x );
        } );
    double dx = x[1] - x[0];

    // velocity
    double a = 1.0;

    // transport problem with centred difference of order 2
    auto transport = [a, n_x, dx]( asc::state_t const& y, asc::state_t& dy, double t )
    {
        dy[0] = -a * ( y[1] - y[n_x - 1] ) / ( 2. * dx );

        for ( std::size_t i = 1; i < n_x - 1; ++i )
        {
            dy[i] = -a * ( y[i + 1] - y[i - 1] ) / ( 2. * dx );
        }

        dy[n_x - 1] = -a * ( y[0] - y[n_x - 2] ) / ( 2. * dx );
    };

    // initial condition
    double sigma = 0.1;
    asc::state_t y0( n_x );
    for ( std::size_t i = 0; i < n_x; ++i )
    {
        // y0[i] = std::sin( 2. * std::numbers::pi * x[i] );
        // y0[i] = 1. / ( sigma * std::sqrt( 2. * std::numbers::pi ) ) * std::exp( -( x[i] - 0.5 ) * ( x[i] - 0.5 ) / ( sigma * sigma ) );
        y0[i] = 0;
        if ( 0.25 <= x[i] && x[i] < 0.5 )
        {
            y0[i] = x[i] - 0.25;
        }
        if ( 0.5 <= x[i] && x[i] < 0.75 )
        {
            y0[i] = -x[i] + 0.75;
        }
    }

    asc::RK4 integrator;

    double tn       = t0;
    asc::state_t yn = y0;

    std::vector<asc::state_t> sol;
    sol.reserve( static_cast<std::size_t>( ( tf - t0 ) / dt ) );
    sol.push_back( yn );

    while ( tn < tf )
    {
        integrator( transport, yn, tn, dt );
        sol.push_back( yn );
    }

    std::ofstream output( "transport.txt" );
    for ( std::size_t i = 0; i < n_x; ++i )
    {
        output << x[i];
        for ( auto const& y_n : sol )
        {
            output << " " << y_n[i];
        }
        output << "\n";
    }

    return 0;
}
