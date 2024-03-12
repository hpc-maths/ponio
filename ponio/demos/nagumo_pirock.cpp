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

template <typename f_explicit_type, typename f_implicit_type, typename f_implicit_t_type>
struct pirock_pb
{
    f_explicit_type f_explicit;
    f_implicit_type f_implicit;
    f_implicit_t_type f_implicit_t;

    pirock_pb( f_explicit_type&& f_explicit_, f_implicit_type&& f_implicit_, f_implicit_t_type&& f_implicit_t_ )
        : f_explicit( f_explicit_ )
        , f_implicit( f_implicit_ )
        , f_implicit_t( f_implicit_t_ )
    {
    }

    // template <typename state_t>
    // auto
    // operator()( double t, state_t&& u )
    // {
    //     return fd( t, u ) + fr( t, u );
    // }
};

template <typename fd_type, typename fr_type, typename fr_t_type>
auto
make_pirock_pb( fd_type&& fd, fr_type&& fr, fr_t_type&& fr_t )
{
    return pirock_pb<fd_type, fr_type, fr_t_type>( std::forward<fd_type>( fd ), std::forward<fr_type>( fr ), std::forward<fr_t_type>( fr_t ) );
}

int
main( int argc, char** argv )
{
    PetscInitialize( &argc, &argv, 0, nullptr );

    constexpr std::size_t dim = 1; // cppcheck-suppress unreadVariable
    using config_t            = samurai::MRConfig<dim>;
    using box_t               = samurai::Box<double, dim>;
    using point_t             = typename box_t::point_t;

    // simulation parameters --------------------------------------------------
    constexpr double d = .1;
    constexpr double k = 1. / d;

    constexpr double left_box  = -40;
    constexpr double right_box = 10;
    constexpr double t_ini     = 0.;
    constexpr double t_end     = 35.;

    // multiresolution parameters
    std::size_t min_level = 0;
    std::size_t max_level = 6;
    double mr_epsilon     = 1e-5; // Threshold used by multiresolution
    double mr_regularity  = 1.;   // Regularity guess for multiresolution

    // output parameters
    std::string const dirname = "nagumo_pirock_data";
    fs::path path             = std::filesystem::path( dirname );
    std::string filename      = "u";
    fs::create_directories( path );

    // define mesh
    point_t box_corner1, box_corner2;
    box_corner1.fill( left_box );
    box_corner2.fill( right_box );
    box_t box( box_corner1, box_corner2 );
    samurai::MRMesh<config_t> mesh{ box, min_level, max_level };

    // init solution ----------------------------------------------------------
    auto u_ini = samurai::make_field<1>( "u", mesh );

    auto exact_solution = [&]( double x, double t )
    {
        double x0  = -25.;
        double v   = ( 1. / std::sqrt( 2. ) ) * std::sqrt( k * d );
        double cst = -( 1. / std::sqrt( 2. ) ) * std::sqrt( k / d );
        double e   = std::exp( cst * ( x - x0 - v * t ) );
        return e / ( 1. + e );
    };

    samurai::for_each_cell( mesh,
        [&]( auto& cell )
        {
            u_ini[cell] = exact_solution( cell.center( 0 ), 0 );
        } );
    samurai::make_bc<samurai::Neumann>( u_ini, 0. );

    // define problem ---------------------------------------------------------

    // diffusion terme
    auto diff = samurai::make_diffusion_order2<decltype( u_ini )>( d );
    auto fd   = [&]( double /* t */, auto&& u )
    {
        samurai::make_bc<samurai::Neumann>( u, 0. );
        samurai::update_ghost_mr( u );
        return -diff( u );
    };

    // reaction terme
    using cfg  = samurai::LocalCellSchemeConfig<samurai::SchemeType::NonLinear, 1, decltype( u_ini )>;
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
    auto fr = [&]( double t, auto&& u )
    {
        samurai::make_bc<samurai::Neumann>( u, 0. );
        samurai::update_ghost_mr( u );
        return fr_t( t )( u );
    };

    auto pb = ponio::make_imex_operator_problem( fd, fr, fr_t );

    ponio::time_span<double> const tspan = { t_ini, t_end };
    double dt                            = ( t_end - t_ini ) / 2000;

    // time loop  -------------------------------------------------------------
    auto sol_range = ponio::make_solver_range( pb, ponio::runge_kutta::pirock::pirock(), u_ini, tspan, dt );

    auto it_sol = sol_range.begin();

    // preapre MR for solution on iterator
    auto mr_adaptation = samurai::make_MRAdapt( it_sol->state );
    samurai::make_bc<samurai::Neumann>( it_sol->state, 0. );
    mr_adaptation( mr_epsilon, mr_regularity );
    samurai::update_ghost_mr( it_sol->state );

    std::size_t n_save = 0;
    save( path, filename, it_sol->state, fmt::format( "_ite_{}", n_save++ ) );

    while ( it_sol->time < t_end )
    {
        // samurai::make_bc<samurai::Neumann>( it_sol->state, 0. );
        //  TODO: add a callback function to make this before each iteration
        for ( auto& ki : it_sol.meth.kis )
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
