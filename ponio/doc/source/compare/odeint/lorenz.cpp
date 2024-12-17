#include <fstream>
#include <vector>

#include <boost/numeric/odeint.hpp>

int
main()
{
    using state_t = std::vector<double>;

    double const sigma = 10.0;
    double const rho   = 28.0;
    double const beta  = 8.0 / 3.0;

    auto lorenz = [&]( state_t const& y, state_t& dy, double t )
    {
        dy[0] = sigma * ( y[1] - y[0] );
        dy[1] = y[0] * ( rho - y[2] ) - y[1];
        dy[2] = y[0] * y[1] - beta * y[2];
    };

    std::vector<std::vector<double>> sol;
    auto vec_observer = [&]( state_t const& y, double const t )
    {
        sol.push_back( {
            {t, y[0], y[1], y[2]}
        } );
    };

    state_t y0 = {
        {1., 1., 1.}
    };
    double t0 = 0.;
    double tf = 20.0;
    double dt = 0.05;

    auto rk4 = boost::numeric::odeint::runge_kutta4<state_t>();
    boost::numeric::odeint::integrate_const( rk4, lorenz, y0, t0, tf, dt, vec_observer );

    std::ofstream output( "lorenz.txt" );
    for ( auto& t_y : sol )
    {
        output << t_y[0] << " " << t_y[1] << " " << t_y[2] << " " << t_y[3] << "\n";
    }

    return 0;
}
