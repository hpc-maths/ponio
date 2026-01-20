// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <iostream>

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
#include <ponio/samurai_linear_algebra.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

#include <samurai/field.hpp>
#include <samurai/io/hdf5.hpp>
#include <samurai/mr/adapt.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/samurai.hpp>
#include <samurai/schemes/fv.hpp>

namespace fs = std::filesystem;

template <class field_t>
void
save( fs::path const& path, std::string const& filename, field_t& u, std::string const& suffix = "" )
{
    auto mesh   = u.mesh();
    auto level_ = samurai::make_scalar_field<std::size_t>( "level", mesh );
    u.name()    = "u";

    if ( !fs::exists( path ) )
    {
        fs::create_directory( path );
    }

    samurai::for_each_cell( mesh,
        [&]( auto& cell )
        {
            level_[cell] = cell.level;
        } );

    samurai::save( path, fmt::format( "{}{}", filename, suffix ), mesh, u, level_ );
}

int
main( int argc, char** argv )
{
    auto& app = samurai::initialize( "Example for the Belousov-Zhabotinsky equation with samurai solved with PIROCK method", argc, argv );
    SAMURAI_PARSE( argc, argv );

    constexpr std::size_t dim = 1; // cppcheck-suppress unreadVariable
    using config_t            = samurai::MRConfig<dim>;
    using box_t               = samurai::Box<double, dim>;
    using point_t             = typename box_t::point_t;

    // simulation parameters --------------------------------------------------
    constexpr double Da = 2.5e-3;
    constexpr double Db = 2.5e-3;
    constexpr double Dc = 1.5e-3;

    constexpr double epsilon = 1e-2;
    constexpr double mu      = 1e-5;
    constexpr double f       = 3.0;
    constexpr double q       = 2e-4;

    constexpr double left_box  = 0.;
    constexpr double right_box = 1.;
    constexpr double t_ini     = 0.;
    constexpr double t_end     = 1.;

    // multiresolution parameters
    std::size_t const min_level = 2;
    std::size_t const max_level = 8;

    // output parameters
    std::string const dirname  = "belousov_zhabotinsky_pirock_data";
    fs::path const path        = std::filesystem::path( dirname );
    std::string const filename = "u";
    fs::create_directories( path );

    // define mesh
    point_t box_corner1;
    point_t box_corner2;
    box_corner1.fill( left_box );
    box_corner2.fill( right_box );
    box_t const box( box_corner1, box_corner2 );
    samurai::MRMesh<config_t> mesh{ box, min_level, max_level };

    // init solution ----------------------------------------------------------
    // auto u_ini = init( mesh );
    auto u_ini = samurai::make_vector_field<double, 3>( "u", mesh );

    double a = 0.;
    double b = 0.;
    double c = 0.;
    u_ini.fill( 0 );
    samurai::for_each_cell( mesh,
        [&]( auto& cell )
        {
            if ( cell.center()[0] < ( right_box - left_box ) / 20. )
            {
                double const y_lim  = 0.05;
                double const x_coor = 0.5;
                double const y_coor = 20.0 * cell.center()[0] - y_lim;

                if ( y_coor >= 0. && y_coor <= 0.3 * x_coor )
                {
                    b = 0.8;
                }
                else
                {
                    b = q * ( f + 1. ) / ( f - 1. );
                }

                if ( y_coor >= 0. )
                {
                    c = q * ( f + 1. ) / ( f - 1. ) + std::atan( y_coor / x_coor ) / ( 8. * std::numbers::pi * f );
                }
                else
                {
                    c = q * ( f + 1. ) / ( f - 1. )
                      + ( std::atan( y_coor / x_coor ) + 2. * std::numbers::pi ) / ( 8. * std::numbers::pi * f );
                }
            }

            a = ( f * c ) / ( q + b );

            // std::cout << cell.center()[0] << "\t a:" << a << " b:" << b << " c:" << c << "\n";

            u_ini[cell]( 0 ) = a;
            u_ini[cell]( 1 ) = b;
            u_ini[cell]( 2 ) = c;
        } );
    samurai::make_bc<samurai::Neumann<1>>( u_ini, 0., 0., 0. );

    // define problem ---------------------------------------------------------

    // diffusion terme
    auto d    = samurai::DiffCoeff<3>( { Da, Db, Dc } );
    auto diff = samurai::make_multi_diffusion_order2<decltype( u_ini )>( d );
    auto fd   = [&]( double /* t */, auto& u, auto& du )
    {
        samurai::make_bc<samurai::Neumann<1>>( u, 0., 0., 0. );
        samurai::update_ghost_mr( u );
        du = -diff( u );
    };

    // reaction terme
    using cfg  = samurai::LocalCellSchemeConfig<samurai::SchemeType::NonLinear, decltype( u_ini ), decltype( u_ini )>;
    auto react = samurai::make_cell_based_scheme<cfg>();
    react.set_name( "Reaction" );
    react.set_scheme_function(
        [&]( auto& scheme_value, auto const& cell, auto const& field )
        {
            auto u  = field[cell];
            auto& a = u[0];
            auto& b = u[1];
            auto& c = u[2];

            scheme_value = { 1. / mu * ( -q * a - a * b + f * c ), 1. / epsilon * ( q * a - a * b + b * ( 1. - b ) ), b - c };
        } );
    // or set option in command line with : -snes_fd -pc_type none
    react.set_jacobian_function(
        [&]( auto& jacobian_matrix, auto const& cell, auto const& field )
        {
            auto u  = field[cell];
            auto& a = u[0];
            auto& b = u[1];
            // auto& c = u[2];

            jacobian_matrix = {
                {( -q - b ) / mu,      -a / mu,                           f / mu},
                { ( q - b ) / epsilon, 1. / epsilon * ( -a - 2 * b + 1 ), 0.    },
                { 0.,                  1.,                                -1.   }
            };
        } );
    auto fr_t = [&]( double /* t */ )
    {
        return react;
    };
    auto fr = [&]( double t, auto&& u, auto& du )
    {
        samurai::make_bc<samurai::Neumann<1>>( u, 0., 0., 0. );
        samurai::update_ghost_mr( u );
        du = fr_t( t )( u );
    };

    ponio::time_span<double> const t_span = { t_ini, t_end };
    double const dt                       = ( t_end - t_ini ) / 2000;

    auto eigmax_computer = [&]( auto&, double, auto&, double, auto& )
    {
        double const dx = mesh.cell_length( mesh.max_level() );
        return 4. / ( dx * dx );
    };

    auto pb = ponio::make_imex_operator_problem( fd, fr, fr_t );

    auto pirock = ponio::runge_kutta::pirock::pirock<1>( ponio::runge_kutta::pirock::beta_0<double>(),
        eigmax_computer,
        ponio::shampine_trick::shampine_trick<decltype( u_ini )>() );

    // time loop  -------------------------------------------------------------
    auto sol_range = ponio::make_solver_range( pb, pirock, u_ini, t_span, dt );
    // auto sol_range = ponio::make_solver_range( pb, ponio::runge_kutta::pirock::pirock_b0( eigmax_computer ), u_ini, t_span, dt );

    auto it_sol = sol_range.begin();

    // preapre MR for solution on iterator
    auto mr_adaptation = samurai::make_MRAdapt( it_sol->state );
    samurai::make_bc<samurai::Neumann<1>>( it_sol->state, 0., 0., 0. );
    auto mra_config = samurai::mra_config().epsilon( 1e-5 ).regularity( 1. );
    mr_adaptation( mra_config );
    samurai::update_ghost_mr( it_sol->state );

    std::size_t n_save = 0;
    save( path, filename, it_sol->state, fmt::format( "_ite_{}", n_save++ ) );

    while ( it_sol->time < t_end )
    {
        samurai::make_bc<samurai::Neumann<1>>( it_sol->state, 0., 0., 0. );
        it_sol.callback_on_stages(
            []( auto& ki )
            {
                ki.resize();
                ki.fill( 0. );
            } );

        ++it_sol;
        std::cout << "tⁿ: " << std::setw( 8 ) << it_sol->time << " (Δt: " << it_sol->time_step << ") " << n_save << "\r";

        mr_adaptation( mra_config );
        samurai::update_ghost_mr( it_sol->state );

        save( path, filename, it_sol->state, fmt::format( "_ite_{}", n_save++ ) );
    }
    std::cout << std::endl;

    samurai::finalize();
    return 0;
}
