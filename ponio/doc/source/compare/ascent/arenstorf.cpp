#include <fstream>
#include <vector>

#include "ascent/Ascent.h"

int
main()
{
    double const mu = 0.012277471;

    auto arenstorf = [=]( asc::state_t const& y, asc::state_t& dy, double t )
    {
        double const y1 = y[0];
        double const y2 = y[1];
        double const y3 = y[2];
        double const y4 = y[3];

        double const r1 = sqrt( ( y1 + mu ) * ( y1 + mu ) + y2 * y2 );
        double const r2 = sqrt( ( y1 - 1. + mu ) * ( y1 - 1. + mu ) + y2 * y2 );

        dy[0] = y3;
        dy[1] = y4;
        dy[2] = y1 + 2. * y4 - ( 1. - mu ) * ( y1 + mu ) / ( r1 * r1 * r1 ) - mu * ( y1 - 1. + mu ) / ( r2 * r2 * r2 );
        dy[3] = y2 - 2. * y3 - ( 1. - mu ) * y2 / ( r1 * r1 * r1 ) - mu * y2 / ( r2 * r2 * r2 );
    };

    std::vector<std::vector<double>> sol;
    auto vec_observer = [&]( asc::state_t const& y, double const t, double const dt )
    {
        sol.push_back( {
            {t, y[0], y[1], y[2], y[3], dt}
        } );
    };

    asc::state_t y0 = {
        {0.994, 0., 0., -2.00158510637908252240537862224}
    };
    double t0 = 0.;
    double tf = 17.0652165601579625588917206249;
    double dt = 1e-5;

    asc::DOPRI45 integrator;
    asc::AdaptiveT<double> adaptive_settings( 1e-5, 1e-5, 0.9 );

    double tn       = t0;
    asc::state_t yn = y0;

    while ( tn < tf )
    {
        vec_observer( yn, tn, dt );
        integrator( arenstorf, yn, tn, dt, adaptive_settings );
    }

    std::ofstream output( "arenstorf.txt" );
    for ( auto& t_y_dt : sol )
    {
        for ( auto& x : t_y_dt )
        {
            output << x << " ";
        }
        output << "\n";
    }

    return 0;
}
