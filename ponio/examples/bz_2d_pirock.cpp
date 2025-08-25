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
    PetscInitialize( &argc, &argv, nullptr, nullptr );

    constexpr std::size_t dim = 2; // cppcheck-suppress unreadVariable
    using config_t            = samurai::MRConfig<dim, 3>;
    using box_t               = samurai::Box<double, dim>;
    using point_t             = typename box_t::point_t;

    // simulation parameters --------------------------------------------------
    constexpr double eps = 1e-2;
    constexpr double mu  = 1e-5;
    constexpr double f   = 1.6;
    constexpr double q   = 2e-3;
    constexpr double da  = 2.5e-3;
    constexpr double db  = 2.5e-3;
    constexpr double dc  = 1.5e-3;

    constexpr double left_box  = 0.;
    constexpr double right_box = 1.;
    constexpr double t_ini     = 0.;
    constexpr double t_end     = 1.;

    // multiresolution parameters
    std::size_t const min_level = 2;
    std::size_t const max_level = 6;
    double const mr_epsilon     = 1e-3; // Threshold used by multiresolution
    double const mr_regularity  = 1.;   // Regularity guess for multiresolution

    // output parameters
    std::string const dirname  = "bz_2d_pirock_data";
    fs::path const path        = std::filesystem::path( dirname );
    std::string const filename = "y";
    fs::create_directories( path );

    // define mesh
    point_t box_corner1;
    point_t box_corner2;
    box_corner1.fill( left_box );
    box_corner2.fill( right_box );
    box_t const box( box_corner1, box_corner2 );
    std::array<bool, dim> const periodic = { false, false };
    samurai::MRMesh<config_t> mesh{ box, min_level, max_level, periodic };

    // init solution ----------------------------------------------------------
    auto y_ini = samurai::make_vector_field<double, 3>( "y", mesh );

    y_ini.fill( 0 );
    samurai::for_each_cell( mesh,
        [&]( auto& cell )
        {
            auto& aij = y_ini[cell][0];
            auto& bij = y_ini[cell][1];
            auto& cij = y_ini[cell][2];

            auto xi = cell.center()[0] - 0.5;
            auto yj = cell.center()[1] - 0.5;

            double theta = std::atan( yj / xi );

            double b = 0.;
            double c = 0.;

            if ( xi < 0. )
            {
                theta = std::numbers::pi + theta;
            }
            else if ( yj < 0. )
            {
                theta = 2. * std::numbers::pi + theta;
            }

            b = q * ( ( f + 1. ) / ( f - 1. ) );
            c = b;

            if ( theta < 0.5 )
            {
                bij = 0.8;
            }
            else
            {
                bij = b;
            }
            cij = c + theta / ( 8. * std::numbers::pi * f );
            aij = ( f * cij ) / ( q + bij );
        } );

    // define problem ---------------------------------------------------------

    // diffusion terme
    auto diff_coeff = samurai::DiffCoeff<3>( { da, db, dc } );
    auto diff       = samurai::make_multi_diffusion_order2<decltype( y_ini )>( diff_coeff );
    auto fd         = [&]( double /* t */, auto&& y, auto& dy )
    {
        samurai::update_ghost_mr( y );
        dy = -diff( y );
    };

    // reaction terme
    using cfg  = samurai::LocalCellSchemeConfig<samurai::SchemeType::NonLinear, decltype( y_ini )::n_comp, decltype( y_ini )>;
    auto react = samurai::make_cell_based_scheme<cfg>();
    react.set_name( "Reaction" );
    react.set_scheme_function(
        [&]( auto const& cell, auto const& y ) -> samurai::SchemeValue<cfg>
        {
            auto& aij = y[cell][0];
            auto& bij = y[cell][1];
            auto& cij = y[cell][2];

            // clang-format off
            return {
                (-q*aij - aij*bij + f*cij)/mu,
                (q*aij - aij*bij + bij*(1. - bij))/eps,
                bij - cij
            };
            // clang-format on
        } );
    // or set option in command line with : -snes_fd -pc_type none
    react.set_jacobian_function(
        [&]( auto const& cell, auto const& y ) -> samurai::JacobianMatrix<cfg>
        {
            auto& aij = y[cell][0];
            auto& bij = y[cell][1];
            // auto& cij = y[cell][2];

            // clang-format off
            return {
                { (-q - bij)/mu, -aij/mu                 , f/mu },
                { (q - bij)/eps, (-aij - 2.*bij + 1.)/eps, 0.   },
                { 0.           , 1.                      , -1.  }
            };
            // clang-format on
        } );
    auto fr_t = [&]( double /* t */ )
    {
        return react;
    };
    auto fr = [&]( double t, auto&& uv, auto& dt_uv )
    {
        // samurai::update_ghost_mr( uv );
        dt_uv = fr_t( t )( uv );
    };

    // spectral radius estimator for diffusion term
    auto eigmax_computer = [&]( auto&, double, auto&, double, auto& )
    {
        double const dx = mesh.cell_length( max_level );
        return da * 4. / ( dx * dx );
    };

    auto pb = ponio::make_imex_operator_problem( fd, fr, fr_t );

    // clang-format off
    [[maybe_unused]] auto pirock_b0 = ponio::runge_kutta::pirock::pirock<1>(
        ponio::runge_kutta::pirock::beta_0<double>(),
        eigmax_computer
    );

    [[maybe_unused]] auto pirock_b0_st = ponio::runge_kutta::pirock::pirock<1, true>(
        ponio::runge_kutta::pirock::beta_0<double>(),
        eigmax_computer,
        ponio::shampine_trick::shampine_trick<decltype( y_ini )>()
    );
    // clang-format on

    auto& method = pirock_b0_st;

    // time loop  -------------------------------------------------------------
    ponio::time_span<double> const t_span = { t_ini, t_end };
    double const dt                       = ( t_end - t_ini ) / 1000;

    auto sol_range = ponio::make_solver_range( pb, method, y_ini, t_span, dt );
    auto it_sol    = sol_range.begin();

    // add boundary conditions
    samurai::make_bc<samurai::Neumann<>>( it_sol->state, 0., 0., 0. );

    // preapre MR for solution on iterator
    auto b_field = samurai::make_scalar_field<double>( "b", mesh );
    samurai::for_each_interval( mesh,
        [&]( auto level, auto& i, auto& idx )
        {
            b_field( level, i, idx ) = it_sol->state( 1, level, i, idx );
        } );

    // auto mr_adaptation = samurai::make_MRAdapt( it_sol->state );
    auto mr_adaptation = samurai::make_MRAdapt( b_field );
    mr_adaptation( mr_epsilon, mr_regularity, it_sol->state );
    samurai::update_ghost_mr( it_sol->state );

    std::size_t n_save = 0;
    save( path, filename, it_sol->state, fmt::format( "_ite_{}", n_save++ ) );

    while ( it_sol->time < t_end )
    {
        it_sol.callback_on_stages(
            []( auto& ki )
            {
                ki.resize();
                ki.fill( 0. );
            } );

        ++it_sol;
        std::cout << "tⁿ: " << std::setw( 8 ) << it_sol->time << " (Δt: " << it_sol->time_step << ") " << ++n_save << "\r";

        samurai::for_each_interval( mesh,
            [&]( auto level, auto& i, auto& idx )
            {
                b_field( level, i, idx ) = it_sol->state( 1, level, i, idx );
            } );

        mr_adaptation( mr_epsilon, mr_regularity, it_sol->state );
        samurai::update_ghost_mr( it_sol->state );

        // save( path, filename, it_sol->state, fmt::format( "_ite_{}", n_save ) );
    }
    std::cout << std::endl;
    save( path, filename, it_sol->state, fmt::format( "_final", n_save++ ) );

    PetscFinalize();

    return 0;
}
