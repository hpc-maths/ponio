// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// NOLINTBEGIN(misc-include-cleaner)

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <string>
#include <valarray>

#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/splitting.hpp>

// NOLINTEND(misc-include-cleaner)

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

    void
    phi_1( double, state_t const& y, state_t& dy ) const // NOLINT(readability-convert-member-functions-to-static)
    {
        double const y3 = y[2];
        double const y4 = y[3];

        dy[0] = y3;
        dy[1] = y4;
        dy[2] = 0.;
        dy[3] = 0.;
    }

    void
    phi_2( double, state_t const& y, state_t& dy ) const
    {
        double const y1 = y[0];
        double const y2 = y[1];
        double const y3 = y[2];
        double const y4 = y[3];

        double const r1 = sqrt( ( y1 + mu ) * ( y1 + mu ) + y2 * y2 );
        double const r2 = sqrt( ( y1 - 1 + mu ) * ( y1 - 1 + mu ) + y2 * y2 );

        dy[0] = 0.;
        dy[1] = 0.;
        dy[2] = y1 + 2 * y4 - ( 1 - mu ) * ( y1 + mu ) / ( r1 * r1 * r1 ) - mu * ( y1 - 1 + mu ) / ( r2 * r2 * r2 );
        dy[3] = y2 - 2 * y3 - ( 1 - mu ) * y2 / ( r1 * r1 * r1 ) - mu * y2 / ( r2 * r2 * r2 );
    }

    void
    operator()( double, state_t const& y, state_t& dy ) const
    {
        double const y1 = y[0];
        double const y2 = y[1];
        double const y3 = y[2];
        double const y4 = y[3];

        double const r1 = sqrt( ( y1 + mu ) * ( y1 + mu ) + y2 * y2 );
        double const r2 = sqrt( ( y1 - 1 + mu ) * ( y1 - 1 + mu ) + y2 * y2 );

        dy[0] = y3;
        dy[1] = y4;
        dy[2] = y1 + 2 * y4 - ( 1 - mu ) * ( y1 + mu ) / ( r1 * r1 * r1 ) - mu * ( y1 - 1 + mu ) / ( r2 * r2 * r2 );
        dy[3] = y2 - 2 * y3 - ( 1 - mu ) * y2 / ( r1 * r1 * r1 ) - mu * y2 / ( r2 * r2 * r2 );
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

    auto arenstorf_pb = arenstorf_model( mu );

    state_t const yini = { 0.994, 0., 0., -2.00158510637908252240537862224 };

    filename = ( std::filesystem::path( dirname ) / "arenstorf_rk118.dat" ).string();
    ponio::solve( arenstorf_pb, ponio::runge_kutta::rk_118(), yini, { 0., tf }, dt, ponio::observer::file_observer( filename ) );

    filename = ( std::filesystem::path( dirname ) / "arenstorf_rk546m.dat" ).string();
    ponio::solve( arenstorf_pb, ponio::runge_kutta::rk54_6m( 1e-5 ), yini, { 0., tf }, dt, ponio::observer::file_observer( filename ) );

    filename = ( std::filesystem::path( dirname ) / "arenstorf_rk547m.dat" ).string();
    ponio::solve( arenstorf_pb, ponio::runge_kutta::rk54_7m( 1e-5 ), yini, { 0., tf }, dt, ponio::observer::file_observer( filename ) );

    filename = ( std::filesystem::path( dirname ) / "arenstorf_rk547s.dat" ).string();
    ponio::solve( arenstorf_pb, ponio::runge_kutta::rk54_7s( 1e-5 ), yini, { 0., tf }, dt, ponio::observer::file_observer( filename ) );

    // adaptive Strang method with Lipschitz estimate
    using namespace std::placeholders;
    auto phi_1 = [&]( double t, state_t const& y, state_t& dy )
    {
        arenstorf_pb.phi_1( t, y, dy );
    };
    auto phi_2 = [&]( double t, state_t const& y, state_t& dy )
    {
        arenstorf_pb.phi_2( t, y, dy );
    };

    auto arenstorf_splitting_pb = ponio::make_problem( phi_1, phi_2 );

    auto adaptive_strang = ponio::splitting::make_adaptive_strang_tuple( 0.05,
        1e-5,
        std::make_pair( ponio::runge_kutta::rk_44(), 0.5 * dt ),
        std::make_pair( ponio::runge_kutta::rk_44(), 0.5 * dt ) );

    filename       = ( std::filesystem::path( dirname ) / "arenstorf_adaptive_strang.dat" ).string();
    auto sol_range = ponio::make_solver_range( arenstorf_splitting_pb, adaptive_strang, yini, { 0., tf }, dt );
    auto it_sol    = sol_range.begin();

    auto obs = ponio::observer::file_observer( filename );

    static constexpr std::size_t N_delta = 50;
    std::size_t n_iteration              = 0;

    // constant for safety factor
    static constexpr double beta  = 0.1;
    static constexpr double gamma = 0.95;

    double dt_star = 10;
    double C_delta = 0.;
    double C_0     = 0.;

    while ( it_sol->time < tf )
    {
        obs( it_sol->time, it_sol->state, it_sol->time_step );

        ++it_sol;

        if ( n_iteration % N_delta == 0 || ( it_sol->time_step < beta * dt_star || gamma * dt_star < it_sol->time_step ) )
        {
            auto [omega,
                C0] = it_sol.meth.lipschitz_constant_estimate( arenstorf_splitting_pb, it_sol->time, it_sol->state, it_sol->time_step );
            C_0     = C0;

            dt_star = it_sol.meth.info().error / ( C0 * it_sol->time_step * it_sol->time_step );
            C_delta = it_sol.meth.info().error / ( it_sol->time_step * it_sol->time_step * it_sol.meth.info().delta );
        }

        if ( it_sol->time_step < beta * dt_star || gamma * dt_star < it_sol->time_step )
        {
            dt_star                  = it_sol->time_step;
            it_sol.meth.info().delta = dt_star * C_0 / C_delta;
        }

        if ( it_sol.meth.info().success )
        {
            ++n_iteration;
        }
    }

    return 0;
}
