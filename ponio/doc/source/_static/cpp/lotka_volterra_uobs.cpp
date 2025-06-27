// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <fstream>
#include <valarray>

#include <ponio/observer.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

struct LV_observer
{
    std::ofstream output;
    double a;
    double b;
    double d;
    double g;

    LV_observer( std::string&& filename, double alpha, double beta, double delta, double gamma )
        : output( filename )
        , a( alpha )
        , b( beta )
        , d( delta )
        , g( gamma )
    {
    }

    ~LV_observer()
    {
        output.close();
    }

    double
    V( std::valarray<double> const& u ) const
    {
        double x = u[0];
        double y = u[1];

        return d * x - g * std::log( x ) + b * y - a * std::log( y );
    }

    void
    operator()( double t, std::valarray<double>& u, double /* dt */ )
    {
        output << t << " " << u[0] << " " << u[1];
        output << " " << V( u ) << "\n";
    }
};

int
main()
{
    using state_t = std::valarray<double>;

    // parameters
    double alpha = 2. / 3.;
    double beta  = 4. / 3.;
    double gamma = 1.;
    double delta = 1.;

    auto lotka_volterra_pb = [=]( double, auto const& u, state_t& du )
    {
        du[0] = alpha * u[0] - beta * u[0] * u[1];
        du[1] = delta * u[0] * u[1] - gamma * u[1];
    };

    auto obs = LV_observer( "lotka_volterra_uobs.txt", alpha, beta, delta, gamma );

    ponio::time_span<double> const t_span = { 0., 15. }; // begin and end time
    double const dt                       = 0.1;         // time step

    state_t const u0 = { 1., 1. }; // initial condition

    ponio::solve( lotka_volterra_pb, ponio::runge_kutta::rk_33(), u0, t_span, dt, obs );

    return 0;
}
