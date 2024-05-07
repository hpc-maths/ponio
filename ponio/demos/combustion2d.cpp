// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <iostream>
#include <valarray>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <numbers>
#include <sstream>

#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/runge_kutta/pirock.hpp>
// #include <ponio/samurai_linear_algebra.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

template <typename state_t>
class combustion_2d_model
{
    double m_d;
    double m_alpha;
    double m_delta;
    double m_R;
    std::size_t m_nx;
    std::size_t m_ny;
    double m_dx;
    double m_dy;
    double m_xmin;
    double m_xmax;
    double m_ymin;
    double m_ymax;

  public:

    combustion_2d_model( double d, double alpha, double delta, double R, std::size_t nx, std::size_t ny )
        : m_d( d )
        , m_alpha( alpha )
        , m_delta( delta )
        , m_R( R )
        , m_nx( nx )
        , m_ny( ny )
        , m_dx( 1. / static_cast<double>( nx ) )
        , m_dy( 1. / static_cast<double>( ny ) )
        , m_xmin( 0.0 )
        , m_xmax( 1.0 )
        , m_ymin( 0.0 )
        , m_ymax( 1.0 )
    {
    }

    state_t
    operator()( double, state_t const& y ) const
    {
        double const doverdxdx = m_d / ( m_dx * m_dx );
        double const doverdydy = m_d / ( m_dy * m_dy );

        state_t f( 0., m_nx * m_ny );

        for ( std::size_t j = 0; j < m_ny; ++j )
        {
            for ( std::size_t i = 0; i < m_nx; ++i )
            {
                std::size_t index = i + j * m_nx;

                if ( i == 0 )
                {
                    f[index] += doverdxdx * ( -2 * y[index] + 2 * y[index + 1] );
                }
                else if ( i == m_nx - 1 )
                {
                    f[index] += doverdxdx * ( y[index - 1] - 2 * y[index] + 1 );
                }
                else
                {
                    f[index] += doverdxdx * ( y[index - 1] - 2 * y[index] + y[index + 1] );
                }

                if ( j == 0 )
                {
                    f[index] += doverdydy * ( -2 * y[index] + 2 * y[index + m_nx] );
                }
                else if ( j == m_ny - 1 )
                {
                    f[index] += doverdydy * ( y[index - m_nx] - 2 * y[index] + 1 );
                }
                else
                {
                    f[index] += doverdydy * ( y[index - m_nx] - 2 * y[index] + y[index + m_nx] );
                }
            }
        }

        f += ( m_R / ( m_alpha * m_delta ) ) * ( 1. + m_alpha - y ) * exp( 1. - 1. / y );

        return f;
    }
};

template <typename state_t>
void
save( std::filesystem::path const& path, std::size_t iteration, state_t& u, std::size_t nx, std::size_t ny )
{
    std::stringstream filename;
    filename << "u_" << iteration << ".dat";
    std::ofstream save_file( path / filename.str() );

    for ( std::size_t j = 0; j < ny; ++j )
    {
        for ( std::size_t i = 0; i < nx; ++i )
        {
            std::size_t index = i + j * nx;

            save_file << u[index] << " ";
        }
        save_file << "\n";
    }
    save_file << std::endl;
}

int
main()
{
    using state_t = std::valarray<double>;

    // simulation parameters --------------------------------------------------
    constexpr double d     = 1.;
    constexpr double R     = 5.;
    constexpr double alpha = 1.;
    constexpr double delta = 20.;

    constexpr double t_ini = 0.;
    constexpr double t_end = 0.26;

    constexpr std::size_t nx = 201;
    constexpr std::size_t ny = 201;

    auto pb = combustion_2d_model<state_t>( d, alpha, delta, R, nx, ny );
    state_t u_ini( 1., nx * ny );

    double dx = 1. / static_cast<double>( nx );
    double dy = 1. / static_cast<double>( ny );
    // for ( std::size_t j = 0; j < ny; ++j )
    // {
    //     for ( std::size_t i = 0; i < nx; ++i )
    //     {
    //         std::size_t index = i + j * nx;
    //         double x          = i * dx - 0.5;
    //         double y          = j * dy - 0.5;

    //         u_ini[index] += std::exp( -( x * x + y * y ) / 0.005 );
    //     }
    // }

    auto eigmax_computer = [=]( auto&, double, state_t&, double )
    {
        return 200. * 4. / ( dx * dx );
    };

    // output -----------------------------------------------------------------
    std::string const dirname  = "combustion2d_data";
    std::filesystem::path path = std::filesystem::path( dirname );
    std::filesystem::create_directories( path );

    // time loop  -------------------------------------------------------------
    static constexpr bool is_embedded = false;

    ponio::time_span<double> const t_span = { t_ini, t_end };
    double dt                             = 0.000001; // ( t_end - t_ini ) / 1000;

    // auto sol_range = ponio::make_solver_range( pb, ponio::runge_kutta::rock::rock2<is_embedded>( eigmax_computer ), u_ini, t_span, dt );
    auto sol_range = ponio::make_solver_range( pb, ponio::runge_kutta::rk_44(), u_ini, t_span, dt );
    auto it_sol    = sol_range.begin();

    std::size_t n_save = 0;
    save( path, n_save, it_sol->state, nx, ny );

    while ( it_sol->time < t_end )
    {
        ++it_sol;
        ++n_save;
        std::cout << "tⁿ: " << std::setw( 8 ) << it_sol->time << " (Δt: " << it_sol->time_step << ") " << n_save << "\r";
        save( path, n_save, it_sol->state, nx, ny );
    }
    std::cout << std::endl;

    return 0;
}
