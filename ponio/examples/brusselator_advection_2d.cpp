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
#include <samurai/hdf5.hpp>
#include <samurai/mr/adapt.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/samurai.hpp>
#include <samurai/schemes/fv.hpp>

namespace fs = std::filesystem;

namespace samurai
{
    /**
     * Linear convection, discretized by the WENO5 (Jiang & Shu) scheme.
     * @param velocities: constant velocity vectors, one for each component of the field.
     */
    template <class Field>
    auto
    make_multi_convection_weno5( std::array<samurai::VelocityVector<Field::dim>, Field::size> const& velocities )
    {
        static_assert( Field::mesh_t::config::ghost_width >= 3, "WENO5 requires at least 3 ghosts." );

        static constexpr std::size_t dim               = Field::dim;
        static constexpr std::size_t field_size        = Field::size;
        static constexpr std::size_t output_field_size = field_size;
        static constexpr std::size_t stencil_size      = 6;

        using cfg = samurai::FluxConfig<samurai::SchemeType::NonLinear, output_field_size, stencil_size, Field>;

        samurai::FluxDefinition<cfg> weno5;

        samurai::static_for<0, dim>::apply( // for each positive Cartesian direction 'd'
            [&]( auto integral_constant_d )
            {
                static constexpr std::size_t d = decltype( integral_constant_d )::value;

                // Stencil creation:
                //        weno5[0].stencil = {{-2, 0}, {-1, 0}, {0,0}, {1,0}, {2,0}, {3,0}};
                //        weno5[1].stencil = {{ 0,-2}, { 0,-1}, {0,0}, {0,1}, {0,2}, {0,3}};
                weno5[d].stencil = line_stencil<dim, d>( -2, -1, 0, 1, 2, 3 );

                weno5[d].cons_flux_function = [&velocities]( auto& cells, Field const& u )
                {
                    samurai::FluxValue<cfg> flux;
                    for ( std::size_t c = 0; c < Field::size; ++c )
                    {
                        auto& velocity = velocities[c];
                        if ( velocity( d ) >= 0 )
                        {
                            samurai::Array<double, 5> f(
                                { u[cells[0]]( c ), u[cells[1]]( c ), u[cells[2]]( c ), u[cells[3]]( c ), u[cells[4]]( c ) } );
                            f *= velocity( d );
                            flux( c ) = samurai::compute_weno5_flux( f );
                        }
                        else
                        {
                            samurai::Array<double, 5> f(
                                { u[cells[5]]( c ), u[cells[4]]( c ), u[cells[3]]( c ), u[cells[2]]( c ), u[cells[1]]( c ) } );
                            f *= velocity( d );
                            flux( c ) = samurai::compute_weno5_flux( f );
                        }
                    }
                    return flux;
                };
            } );

        return samurai::make_flux_based_scheme( weno5 );
    }
}

template <class field_t>
void
save( fs::path const& path, std::string const& filename, field_t& u, std::string const& suffix = "" )
{
    auto mesh   = u.mesh();
    auto level_ = samurai::make_field<std::size_t, 1>( "level", mesh );
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
    constexpr double A = 1.3;
    constexpr double B = 1.;

    samurai::VelocityVector<dim> const U = { -0.5, 1 };
    samurai::VelocityVector<dim> const V = { 0.4, 0.7 };

    constexpr double nu = 1e-2;
    constexpr double mu = 1.;

    constexpr double left_box  = 0.;
    constexpr double right_box = 1.;
    constexpr double t_ini     = 0.;
    constexpr double t_end     = 1.;

    // multiresolution parameters
    std::size_t const min_level = 6;
    std::size_t const max_level = 6;
    double const mr_epsilon     = 1e-5; // Threshold used by multiresolution
    double const mr_regularity  = 1.;   // Regularity guess for multiresolution

    // output parameters
    std::string const dirname  = "brusselator_advection_2d_data";
    fs::path const path        = std::filesystem::path( dirname );
    std::string const filename = "uv";
    fs::create_directories( path );

    // define mesh
    point_t box_corner1;
    point_t box_corner2;
    box_corner1.fill( left_box );
    box_corner2.fill( right_box );
    box_t const box( box_corner1, box_corner2 );
    std::array<bool, dim> const periodic = { true, true };
    samurai::MRMesh<config_t> mesh{ box, min_level, max_level, periodic };

    // init solution ----------------------------------------------------------
    auto uv_ini = samurai::make_field<2>( "uv", mesh );

    uv_ini.fill( 0 );
    samurai::for_each_cell( mesh,
        [&]( auto& cell )
        {
            auto x = cell.center()[0];
            auto y = cell.center()[1];

            uv_ini[cell]( 0 ) = 22. * y * std::pow( 1 - y, 1.5 );
            uv_ini[cell]( 1 ) = 27. * x * std::pow( 1 - x, 1.5 );
        } );

    // define problem ---------------------------------------------------------

    // diffusion terme
    auto diff = samurai::make_diffusion_order2<decltype( uv_ini )>( nu );
    auto fd   = [&]( double /* t */, auto&& uv, auto& dt_uv )
    {
        samurai::update_ghost_mr( uv );
        dt_uv = -diff( uv );
    };

    // reaction terme
    using cfg  = samurai::LocalCellSchemeConfig<samurai::SchemeType::NonLinear, decltype( uv_ini )::size, decltype( uv_ini )>;
    auto react = samurai::make_cell_based_scheme<cfg>();
    react.set_name( "Reaction" );
    react.set_scheme_function(
        [&]( auto const& cell, auto const& uv ) -> samurai::SchemeValue<cfg>
        {
            auto& u = uv[cell][0];
            auto& v = uv[cell][1];

            return { A + u * u * v - ( B + 1 ) * u, B * u - u * u * v };
        } );
    // or set option in command line with : -snes_fd -pc_type none
    react.set_jacobian_function(
        [&]( auto const& cell, auto const& uv ) -> samurai::JacobianMatrix<cfg>
        {
            auto& u = uv[cell][0];
            auto& v = uv[cell][1];

            return {
                {2. * u * v - ( B + 1 ), u * u },
                { B - 2. * u * v,        -u * u}
            };
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
    auto fr_pb = ponio::make_implicit_operator_problem( fr, fr_t );

    // advection terme
    std::array<samurai::VelocityVector<dim>, 2> const velocities = { mu * U, mu * V };

    auto conv = samurai::make_multi_convection_weno5<decltype( uv_ini )>( velocities );
    auto fa   = [&]( double /* t */, auto&& uv, auto& dt_uv )
    {
        samurai::update_ghost_mr( uv );
        dt_uv = -conv( uv );
    };

    ponio::time_span<double> const t_span = { t_ini, t_end };
    double const dt                       = ( t_end - t_ini ) / 500;

    auto eigmax_computer = [&]( auto&, double, auto&, double, auto& )
    {
        double const dx = mesh.cell_length( mesh.max_level() );
        return nu * 4. / ( dx * dx );
    };

    auto pb = ponio::make_problem( fr_pb, fd, fa );

    [[maybe_unused]] auto pirock_b0    = ponio::runge_kutta::pirock::pirock_RDA<1>( ponio::runge_kutta::pirock::beta_0<double>(),
        eigmax_computer );
    [[maybe_unused]] auto pirock_b0_st = ponio::runge_kutta::pirock::pirock_RDA<1>( ponio::runge_kutta::pirock::beta_0<double>(),
        eigmax_computer,
        ponio::shampine_trick::shampine_trick<decltype( uv_ini )>() );

    auto& method = pirock_b0_st;

    // time loop  -------------------------------------------------------------
    auto sol_range = ponio::make_solver_range( pb, method, uv_ini, t_span, dt );

    auto it_sol = sol_range.begin();

    // preapre MR for solution on iterator
    auto mr_adaptation = samurai::make_MRAdapt( it_sol->state );
    mr_adaptation( mr_epsilon, mr_regularity );
    samurai::update_ghost_mr( it_sol->state );

    std::size_t n_save = 0;
    save( path, filename, it_sol->state, fmt::format( "_ite_{}", n_save++ ) );

    while ( it_sol->time < t_end )
    {
        for ( auto& ki : it_sol.stages() )
        {
            ki.resize();
            ki.fill( 0. );
        }

        ++it_sol;
        std::cout << "tⁿ: " << std::setw( 8 ) << it_sol->time << " (Δt: " << it_sol->time_step << ") " << ++n_save << "\r";

        // mr_adaptation( mr_epsilon, mr_regularity );
        samurai::update_ghost_mr( it_sol->state );

        save( path, filename, it_sol->state, fmt::format( "_ite_{}", n_save ) );
    }
    std::cout << std::endl;
    save( path, filename, it_sol->state, fmt::format( "_final", n_save++ ) );

    PetscFinalize();

    return 0;
}
