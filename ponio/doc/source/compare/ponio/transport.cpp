#include <algorithm>
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

    // space parameter
    std::size_t n_x = 500;
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

    // time parameter
    double t0 = 0.;
    double tf = 0.3;
    double dt = dx / a;

    // initial condition
    state_t y0( n_x );
    for ( std::size_t i = 0; i < n_x; ++i )
    {
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

    auto upwind = [a, n_x, dx]( double t, auto&& y, state_t& dy )
    {
        dy[0] = -( std::max( a, 0. ) * ( y[0] - y[n_x - 1] ) + std::min( a, 0. ) * ( y[1] - y[0] ) ) / dx;

        for ( std::size_t i = 1; i < n_x - 1; ++i )
        {
            dy[i] = -( std::max( a, 0. ) * ( y[i] - y[i - 1] ) + std::min( a, 0. ) * ( y[i + 1] - y[i] ) ) / dx;
        }

        dy[n_x - 1] = -( std::max( a, 0. ) * ( y[n_x - 1] - y[n_x - 2] ) + std::min( a, 0. ) * ( y[0] - y[n_x - 1] ) ) / dx;
    };

    auto vec_observer = ponio::observer::vector_observer<state_t, double>();

    ponio::solve( upwind, ponio::runge_kutta::euler(), y0, { t0, tf }, dt, vec_observer );

    // save solution
    std::ofstream output( "transport.txt" );
    for ( std::size_t i = 0; i < n_x; ++i )
    {
        output << x[i];
        for ( auto const& state_n : vec_observer.solutions )
        {
            output << " " << std::get<1>( state_n )[i];
        }
        output << "\n";
    }

    return 0;
}
