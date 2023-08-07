// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <iostream>
#include <numbers>
#include <numeric>
#include <valarray>

#include <solver/observer.hpp>
#include <solver/problem.hpp>
#include <solver/runge_kutta.hpp>
#include <solver/solver.hpp>

/*
solve Arenstorf orbit problem
*/

struct arenstorf_model
{
    using state_t = std::valarray<double>;

    double mu;

    arenstorf_model( double m )
        : mu( m )
    {
    }

    state_t
    operator()( double, state_t const& y ) const
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
    }
};

int
main( int, char** )
{
    std::string const dirname = "arenstorf_data";
    std::string filename;

    using state_t = std::valarray<double>;

    double const tf = 17.0652165601579625588917206249;
    double const dt = 1e-5;

    double const mu = 0.012277471;

    auto arenstorf_pb = ponio::make_problem( arenstorf_model( mu ) );

    state_t const yini = { 0.994, 0., 0., -2.00158510637908252240537862224 };

    filename = ( std::filesystem::path( dirname ) / "arenstorf_rk546m.dat" ).string();
    ponio::solve( arenstorf_pb, ponio::runge_kutta::rk54_6m( 1e-5 ), yini, { 0., tf }, dt, observer::file_observer( filename ) );

    filename = ( std::filesystem::path( dirname ) / "arenstorf_rk547m.dat" ).string();
    ponio::solve( arenstorf_pb, ponio::runge_kutta::rk54_7m( 1e-5 ), yini, { 0., tf }, dt, observer::file_observer( filename ) );

    filename = ( std::filesystem::path( dirname ) / "arenstorf_rk547s.dat" ).string();
    ponio::solve( arenstorf_pb, ponio::runge_kutta::rk54_7s( 1e-5 ), yini, { 0., tf }, dt, observer::file_observer( filename ) );

    return 0;
}
