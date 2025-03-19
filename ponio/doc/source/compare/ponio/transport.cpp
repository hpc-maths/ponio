#include <fstream>
#include <valarray>
#include <vector>

#include <ponio/observer.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

int
main()
{
    using state_t = std::valarray<double>;

    std::size_t n_x = 100;
    double t0       = 0.;
    double tf       = 0.3;
    double dt       = 0.01;

    std::valarray<double> x( n_x );
    std::ranges::generate( x,
        [i = 0, n_x]() mutable
        {
            return static_cast<double>( i++ ) / static_cast<double>( n_x );
        } );
    double dx = x[1] - x[0];

    // velocity
    double a = 1.0;

    // transport problem with centred difference of order 2
    auto transport = [a, n_x, dx]( double t, state_t const& y ) -> state_t
    {
        auto dy = state_t( y.size() );

        dy[0] = -a * ( y[1] - y[n_x - 1] ) / ( 2. * dx );

        for ( std::size_t i = 1; i < n_x - 1; ++i )
        {
            dy[i] = -a * ( y[i + 1] - y[i - 1] ) / ( 2. * dx );
        }

        dy[n_x - 1] = -a * ( y[0] - y[n_x - 2] ) / ( 2. * dx );

        return dy;
    };

    // initial condition
    double sigma = 0.1;
    state_t y0( n_x );
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

    ponio::time_span<double> t_span = { t0, tf };

    std::vector<state_t> sol;
    sol.reserve( static_cast<std::size_t>( ( tf - t0 ) / dt ) );
    auto vec_observer = [&]( double /* t */, state_t const& y, double /* dt */ )
    {
        sol.push_back( y );
    };

    ponio::solve( transport, ponio::runge_kutta::rk_44(), y0, t_span, dt, vec_observer );

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
