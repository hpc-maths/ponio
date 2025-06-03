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

    double const mu = 0.012277471;

    auto arenstorf = [=]( double t, state_t const& y ) -> state_t
    {
        double const y1 = y[0];
        double const y2 = y[1];
        double const y3 = y[2];
        double const y4 = y[3];

        double const r1 = sqrt( ( y1 + mu ) * ( y1 + mu ) + y2 * y2 );
        double const r2 = sqrt( ( y1 - 1 + mu ) * ( y1 - 1 + mu ) + y2 * y2 );

        double const dy1 = y3;
        double const dy2 = y4;
        double const dy3 = y1 + 2 * y4 - ( 1 - mu ) * ( y1 + mu ) / ( r1 * r1 * r1 ) - mu * ( y1 - 1 + mu ) / ( r2 * r2 * r2 );
        double const dy4 = y2 - 2 * y3 - ( 1 - mu ) * y2 / ( r1 * r1 * r1 ) - mu * y2 / ( r2 * r2 * r2 );

        return { dy1, dy2, dy3, dy4 };
    };

    state_t y0 = {
        {0.994, 0., 0., -2.00158510637908252240537862224}
    };
    ponio::time_span<double> t_span = { 0., 17.0652165601579625588917206249 };
    double dt                       = 1e-5;

    ponio::solve( arenstorf, ponio::runge_kutta::rk54_7m( 1e-5 ), y0, t_span, dt, ponio::observer::file_observer( "arenstorf_rk54_7m.txt" ) );
    ponio::solve( arenstorf, ponio::runge_kutta::rk87_13m( 1e-5 ), y0, t_span, dt, ponio::observer::file_observer( "arenstorf_rk87_13m.txt" ) );

    return 0;
}
