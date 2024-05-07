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

int
main( int argc, char** argv )
{
    PetscInitialize( &argc, &argv, 0, nullptr );

    constexpr std::size_t dim = 2; // cppcheck-suppress unreadVariable

    using config_t = samurai::MRConfig<dim, 1>;
    using box_t    = samurai::Box<double, dim>;
    using point_t  = typename box_t::point_t;

    // simulation parameters --------------------------------------------------
    constexpr double d     = 1.;
    constexpr double R     = 5.;
    constexpr double alpha = 1.;
    constexpr double delta = 20.;

    constexpr double t_ini = 0.;
    constexpr double t_end = 0.26;

    // multiresolution parameters
    std::size_t min_level = 6;
    std::size_t max_level = 6;
    double mr_epsilon     = 1e-5; // Threshold used by multiresolution
    double mr_regularity  = 1.;   // Regularity guess for multiresolution

    // output parameters
    std::string const dirname = "combustion_2d_pirock_data";
    fs::path path             = std::filesystem::path( dirname );
    std::string filename      = "u";
    fs::create_directories( path );

    // define mesh
    point_t box_corner_1 = { 0., 0. };
    point_t box_corner_2 = { 1., 1. };
    box_t box( box_corner_1, box_corner_2 );
    samurai::MRMesh<config_t> mesh{ box, min_level, max_level };

    // init solution ----------------------------------------------------------
    auto u_ini = samurai::make_field<1>( "u", mesh );
    u_ini.fill( 1.0 );

    // boundary condition
    samurai::DirectionVector<dim> left   = { -1, 0 };
    samurai::DirectionVector<dim> right  = { 1, 0 };
    samurai::DirectionVector<dim> bottom = { 0, -1 };
    samurai::DirectionVector<dim> top    = { 0, 1 };
    auto set_bc                          = [&]( auto& field )
    {
        samurai::make_bc<samurai::Neumann<1>>( field, 0. )->on( left, bottom );
        samurai::make_bc<samurai::Dirichlet<1>>( field, 1. )->on( right, top );
    };

    set_bc( u_ini );

    // define problem ---------------------------------------------------------

    // diffusion terme
    samurai::DiffCoeff<dim> diff_coeff;
    diff_coeff.fill( d );
    auto diff = samurai::make_diffusion_order2<decltype( u_ini )>( diff_coeff );
    auto fd   = [&]( double /* t */, auto&& u )
    {
        set_bc( u );
        samurai::update_ghost_mr( u );
        // std::cout << xt::amin( u ) << " " << xt::amax( u ) << std::endl;
        return -diff( u );
    };

    // reaction terme
    using cfg  = samurai::LocalCellSchemeConfig<samurai::SchemeType::NonLinear, decltype( u_ini )::size, decltype( u_ini )>;
    auto react = samurai::make_cell_based_scheme<cfg>();
    react.set_name( "Reaction" );
    react.set_scheme_function(
        [&]( auto const& cell, auto const& field ) -> samurai::SchemeValue<cfg>
        {
            auto u = field[cell];

            return ( R / ( alpha * delta ) ) * ( 1 + alpha - u ) * exp( delta * ( 1 - 1 / u ) );
        } );
    // set the jacobian function or set option in command line with : -snes_fd -pc_type none
    react.set_jacobian_function(
        [&]( auto const& cell, auto const& field ) -> samurai::JacobianMatrix<cfg>
        {
            auto u = field[cell];

            return ( R / ( alpha * delta ) ) * exp( delta * ( 1 - 1 / u ) ) * ( -1 + ( 1 + alpha - u ) * ( delta / ( u * u ) ) );
        } );
    auto fr_t = [&]( double /* t */ )
    {
        return react;
    };
    auto fr = [&]( double t, auto&& u )
    {
        set_bc( u );
        samurai::update_ghost_mr( u );
        return fr_t( t )( u );
    };

    auto eigmax_computer = [=]( auto&, double, auto&, double )
    {
        double dx = samurai::cell_length( max_level );
        return 2.01 * 4. * d / ( dx * dx );
    };

    auto pb = ponio::make_imex_operator_problem( fd, fr, fr_t );

    // time loop  -------------------------------------------------------------
    static constexpr bool is_embedded = false;

    ponio::time_span<double> const t_span = { t_ini, t_end };
    double dt                             = 0.00026; // ( t_end - t_ini ) / 1000;

    // auto sol_range = ponio::make_solver_range( pb,
    //     ponio::runge_kutta::pirock::pirock<2, is_embedded>( ponio::runge_kutta::pirock::alpha_fixed<double>( 1.0 ),
    //         eigmax_computer,
    //         ponio::shampine_trick::shampine_trick<decltype( u_ini )>(),
    //         1e-4 ),
    //     u_ini,
    //     t_span,
    //     dt );
    auto sol_range = ponio::make_solver_range( pb, ponio::runge_kutta::rock::rock4<is_embedded>( eigmax_computer ), u_ini, t_span, dt );

    auto it_sol = sol_range.begin();

    // preapre MR for solution on iterator
    auto mr_adaptation = samurai::make_MRAdapt( it_sol->state );

    set_bc( it_sol->state );

    mr_adaptation( mr_epsilon, mr_regularity );
    samurai::update_ghost_mr( it_sol->state );

    std::size_t n_save = 0;
    save( path, filename, it_sol->state, fmt::format( "_ite_{}", n_save++ ) );

    while ( it_sol->time < t_end )
    {
        set_bc( it_sol->state );
        //  TODO: add a callback function to make this before each iteration
        for ( auto& ki : it_sol.meth.kis )
        {
            ki.resize();
            ki.fill( 0. );
            set_bc( ki );
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
