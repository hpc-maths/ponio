// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Heat equation 1D on a large grid, solved with a fixed number of time steps and no output.
// It is a profiling workload: the run time is set from the command line and the program
// prints a throughput, so that runs with different sizes or methods can be compared.

#include <chrono>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <numbers>
#include <string>
#include <valarray>

#include <CLI/CLI.hpp>

#include <ponio/observer.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

// same right hand side as the heat example
class heat_model
{
    double m_dx;

  public:

    std::size_t n_eval = 0;

    heat_model( double dx )
        : m_dx( dx )
    {
    }

    void
    operator()( double, std::valarray<double> const& y, std::valarray<double>& ydot )
    {
        ++n_eval;

        double const oneoverdxdx = 1. / ( m_dx * m_dx );
        std::size_t const nx     = y.size();

        ydot[0] = oneoverdxdx * ( -2. * y[0] + y[1] );
        for ( std::size_t i = 1; i < nx - 1; ++i )
        {
            ydot[i] = oneoverdxdx * ( y[i - 1] - 2. * y[i] + y[i + 1] );
        }
        ydot[nx - 1] = oneoverdxdx * ( y[nx - 2] - 2. * y[nx - 1] );
    }

    static std::valarray<double>
    fundamental_sol( double t, std::valarray<double> const& x )
    {
        double const pi = std::numbers::pi;
        return ( 1. / ( 2. * std::sqrt( pi * t ) ) ) * ( std::exp( -( x * x ) / ( 4. * t ) ) );
    }
};

// time step of each method, in units of dx^2, chosen inside its stability domain
// (the spectral radius of the discrete laplacian is 4/dx^2)
std::map<std::string, double> const dt_over_dx2 = {
    {"rk_44",    0.5},
    { "rkc_202", 20.},
    { "rkc2_15", 10.},
    { "rock2",   20.},
    { "rock4",   20.},
};

template <typename Algorithm_t>
std::valarray<double>
run( heat_model& pb, Algorithm_t&& algo, std::valarray<double> const& u_ini, ponio::time_span<double> const& t_span, double dt )
{
    return ponio::solve( pb, std::forward<Algorithm_t>( algo ), u_ini, t_span, dt, ponio::observer::null_observer() );
}

int
main( int argc, char** argv )
{
    CLI::App app{ "Heat equation 1D, profiling workload" };

    std::size_t nx     = 100'000;
    std::size_t nsteps = 100;
    std::string method = "rk_44";

    app.add_option( "--nx", nx, "Number of grid points" )->check( CLI::Range( 3ul, 1'000'000'000ul ) );
    app.add_option( "--steps", nsteps, "Number of time steps" )->check( CLI::PositiveNumber );
    app.add_option( "--method", method, "Time integration method" )->check( CLI::IsMember( dt_over_dx2 ) );
    CLI11_PARSE( app, argc, argv );

    double const xmin = -5.;
    double const xmax = 5.;
    double const dx   = ( xmax - xmin ) / static_cast<double>( nx + 1 );
    double const dt   = dt_over_dx2.at( method ) * dx * dx;

    std::valarray<double> x( nx );
    for ( std::size_t i = 0; i < nx; ++i )
    {
        x[i] = xmin + dx * static_cast<double>( i + 1 );
    }

    double const t_ini = 0.001;
    // the half step keeps the number of steps exact whatever the rounding of the accumulated time
    double const t_end = t_ini + ( static_cast<double>( nsteps ) - 0.5 ) * dt;

    ponio::time_span<double> const t_span = { t_ini, t_end };

    std::valarray<double> const u_ini = heat_model::fundamental_sol( t_ini, x );
    auto pb                           = heat_model( dx );

    auto eigmax_computer = [=]( auto&, double, auto&, double, auto& )
    {
        return 4. / ( dx * dx );
    };

    auto const start = std::chrono::steady_clock::now();

    std::valarray<double> u_end;
    if ( method == "rk_44" )
    {
        u_end = run( pb, ponio::runge_kutta::rk_44(), u_ini, t_span, dt );
    }
    else if ( method == "rkc_202" )
    {
        u_end = run( pb, ponio::runge_kutta::rkc_202(), u_ini, t_span, dt );
    }
    else if ( method == "rkc2_15" )
    {
        u_end = run( pb, ponio::runge_kutta::explicit_rkc2<15>(), u_ini, t_span, dt );
    }
    else if ( method == "rock2" )
    {
        u_end = run( pb, ponio::runge_kutta::rock::rock2( eigmax_computer ), u_ini, t_span, dt );
    }
    else if ( method == "rock4" )
    {
        u_end = run( pb, ponio::runge_kutta::rock::rock4( eigmax_computer ), u_ini, t_span, dt );
    }

    std::chrono::duration<double> const elapsed = std::chrono::steady_clock::now() - start;

    // the heat equation with homogeneous Dirichlet conditions dissipates: a growing maximum means an unstable time step
    double const max_ini = u_ini.max();
    double const max_end = std::abs( u_end ).max();
    if ( !std::isfinite( max_end ) || max_end > max_ini )
    {
        std::cerr << "unstable run: max |u| went from " << max_ini << " to " << max_end << "\n";
        return 1;
    }

    double const cell_updates = static_cast<double>( nx ) * static_cast<double>( nsteps );
    double const rhs_updates  = static_cast<double>( nx ) * static_cast<double>( pb.n_eval );

    std::cout << "method " << method << ", nx " << nx << ", " << nsteps << " steps, " << pb.n_eval << " rhs evaluations\n";
    std::cout << "time " << elapsed.count() << " s, " << cell_updates / elapsed.count() * 1e-6 << " Mcell.step/s, "
              << rhs_updates / elapsed.count() * 1e-6 << " Mcell.rhs/s\n";
    std::cout << "checksum " << u_end.sum() * dx << "\n";

    return 0;
}
