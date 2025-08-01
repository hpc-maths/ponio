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

    double const sigma = 10.0;
    double const rho   = 28.0;
    double const beta  = 8.0 / 3.0;

    auto lorenz = [&]( double t, auto&& y, state_t& dy )
    {
        dy[0] = sigma * ( y[1] - y[0] );
        dy[1] = y[0] * ( rho - y[2] ) - y[1];
        dy[2] = y[0] * y[1] - beta * y[2];
    };

    state_t y0 = {
        {1., 1., 1.}
    };
    ponio::time_span<double> t_span = { 0., 20.0 };
    double dt                       = 0.01;

    ponio::solve( lorenz, ponio::runge_kutta::rk_44(), y0, t_span, dt, ponio::observer::file_observer( "lorenz.txt" ) );

    return 0;
}
