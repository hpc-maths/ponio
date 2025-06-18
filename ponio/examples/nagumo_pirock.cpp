// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <iostream>
#include <valarray>

#include <algorithm>
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

#include <filesystem>
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

    constexpr std::size_t dim = 1; // cppcheck-suppress unreadVariable
    using config_t            = samurai::MRConfig<dim>;
    using box_t               = samurai::Box<double, dim>;
    using point_t             = typename box_t::point_t;

    // simulation parameters --------------------------------------------------
    constexpr double d  = .1;
    constexpr double k  = 1. / d;
    constexpr double x0 = -25.;

    constexpr double left_box  = -40;
    constexpr double right_box = 10;
    constexpr double t_ini     = 0.;
    constexpr double t_end     = 35.;

    // multiresolution parameters
    std::size_t const min_level = 7;
    std::size_t const max_level = 7;
    double const mr_epsilon     = 1e-5; // Threshold used by multiresolution
    double const mr_regularity  = 1.;   // Regularity guess for multiresolution

    // output parameters
    std::string const dirname  = "nagumo_pirock_data";
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
    auto u_ini = samurai::make_scalar_field<double>( "u", mesh );

    auto exact_solution = [&]( double x, double t )
    {
        double const v   = ( 1. / std::sqrt( 2. ) ) * std::sqrt( k * d );
        double const cst = -( 1. / std::sqrt( 2. ) ) * std::sqrt( k / d );
        double const e   = std::exp( cst * ( x - x0 - v * t ) );
        return e / ( 1. + e );
    };

    samurai::for_each_cell( mesh,
        [&]( auto& cell )
        {
            u_ini[cell] = exact_solution( cell.center( 0 ), 0 );
        } );
    samurai::make_bc<samurai::Neumann<1>>( u_ini, 0. );

    // define problem ---------------------------------------------------------

    // diffusion terme
    auto diff = samurai::make_diffusion_order2<decltype( u_ini )>( d );
    auto fd   = [&]( double /* t */, auto&& u, auto& du )
    {
        samurai::make_bc<samurai::Neumann<1>>( u, 0. );
        samurai::update_ghost_mr( u );
        du = -diff( u );
    };

    // reaction terme
    using cfg  = samurai::LocalCellSchemeConfig<samurai::SchemeType::NonLinear, decltype( u_ini )::n_comp, decltype( u_ini )>;
    auto react = samurai::make_cell_based_scheme<cfg>();
    react.set_name( "Reaction" );
    react.set_scheme_function(
        [&]( auto const& cell, auto const& field )
        {
            auto u = field[cell];
            return k * u * u * ( 1 - u );
        } );
    react.set_jacobian_function(
        [&]( auto const& cell, auto const& field )
        {
            auto u = field[cell];
            return k * ( 2 * u * ( 1 - u ) - u * u );
        } );
    auto fr_t = [&]( double /* t */ )
    {
        return react;
    };
    auto fr = [&]( double t, auto&& u, auto& du )
    {
        samurai::make_bc<samurai::Neumann<1>>( u, 0. );
        samurai::update_ghost_mr( u );
        du = fr_t( t )( u );
    };

    auto pb = ponio::make_imex_operator_problem( fd, fr, fr_t );

    ponio::time_span<double> const tspan = { t_ini, t_end };
    double const dt                      = ( t_end - t_ini ) / 2000;

    auto eigmax_computer = [&]( auto&, double, auto&, double, auto& )
    {
        double const dx = mesh.cell_length( mesh.max_level() );
        return 4. / ( dx * dx );
    };

    // time loop  -------------------------------------------------------------
    auto pirock = ponio::runge_kutta::pirock::pirock<1>( ponio::runge_kutta::pirock::beta_0<double>(),
        eigmax_computer,
        ponio::shampine_trick::shampine_trick<decltype( u_ini )>() );

    auto sol_range = ponio::make_solver_range( pb, pirock, u_ini, tspan, dt );

    auto it_sol = sol_range.begin();

    // preapre MR for solution on iterator
    auto mr_adaptation = samurai::make_MRAdapt( it_sol->state );
    samurai::make_bc<samurai::Neumann<1>>( it_sol->state, 0. );
    mr_adaptation( mr_epsilon, mr_regularity );
    samurai::update_ghost_mr( it_sol->state );

    std::size_t n_save = 0;
    save( path, filename, it_sol->state, fmt::format( "_ite_{}", n_save++ ) );

    while ( it_sol->time < t_end )
    {
        samurai::make_bc<samurai::Neumann<1>>( it_sol->state, 0. );
        samurai::update_ghost_mr( it_sol->state );

        for ( auto& ki : it_sol.stages() )
        {
            ki.resize();
            ki.fill( 0. );
        }

        ++it_sol;
        std::cout << "tⁿ: " << std::setw( 8 ) << it_sol->time << " (Δt: " << it_sol->time_step << ") " << n_save << "\r";

        mr_adaptation( mr_epsilon, mr_regularity );
        samurai::update_ghost_mr( it_sol->state );

        save( path, filename, it_sol->state, fmt::format( "_ite_{}", n_save++ ) );
    }
    std::cout << std::endl;

    PetscFinalize();

    return 0;
}
