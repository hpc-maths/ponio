// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/runge_kutta/pirock.hpp>
#include <ponio/samurai_linear_algebra.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

#include <samurai/field.hpp>
#include <samurai/mr/adapt.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/samurai.hpp>
#include <samurai/schemes/fv.hpp>

#include "compute_order.hpp"

namespace detail
{
    template <typename value_t, value_t begin, typename lambda_t, value_t... ints>
    void
    static_for_impl( lambda_t&& f, std::integer_sequence<value_t, ints...> )
    {
        ( f( std::integral_constant<value_t, ints>() ), ... );
    }
}

template <typename value_t, value_t begin, value_t end, typename lambda_t>
void
static_for( lambda_t&& f )
{
    detail::static_for_impl<value_t, begin>( std::forward<lambda_t>( f ), std::make_integer_sequence<value_t, end - begin>() );
}

template <typename mesh_t>
auto
init_field( std::string const& name, mesh_t& mesh )
{
    if constexpr ( requires { make_scalar_field<double>( name, mesh ); } )
    {
        return make_scalar_field( name, mesh );
    }
    else
    {
        return make_field<1>( name, mesh );
    }
}

template <typename mesh_t>
auto
min_cell_length( mesh_t& mesh )
{
    if constexpr ( requires { mesh.cell_length( mesh.max_level() ); } )
    {
        return mesh.cell_length( mesh.max_level() );
    }
    else
    {
        return samurai::cell_length( mesh.max_level() );
    }
}

TEST_CASE( "samurai::order::pirock" )
{
    int argc                = 1;
    std::vector<char*> argv = { (char*)"ponio_tests" };
    char** argv_c           = argv.data();

    PetscInitialize( &argc, &argv_c, 0, nullptr );

    constexpr std::size_t dim = 1;
    using config_t            = samurai::MRConfig<dim, 3>;
    using box_t               = samurai::Box<double, dim>;
    using point_t             = typename box_t::point_t;

    // simulation parameters --------------------------------------------------
    constexpr double d = .1;
    constexpr double k = 1. / d;

    constexpr double left_box  = -3;
    constexpr double right_box = 5;
    constexpr double t_ini     = 0.;
    constexpr double t_end     = 3.;

    // multiresolution parameters
    std::size_t min_level = 10;
    std::size_t max_level = 10;

    ponio::time_span<double> const tspan = { t_ini, t_end };

    // define mesh
    point_t box_corner1, box_corner2;
    box_corner1.fill( left_box );
    box_corner2.fill( right_box );
    box_t box( box_corner1, box_corner2 );
    samurai::MRMesh<config_t> mesh{ box, min_level, max_level };

    // init solution ----------------------------------------------------------

    auto u_ini = init_field( "u", mesh );

    auto exact_solution = [&]( double x, double t )
    {
        double x0  = 0.;
        double v   = ( 1. / std::sqrt( 2. ) ) * std::sqrt( k * d );
        double cst = -( 1. / std::sqrt( 2. ) ) * std::sqrt( k / d );
        double e   = std::exp( cst * ( x - x0 - v * t ) );
        return e / ( 1. + e );
    };

    samurai::for_each_cell( mesh,
        [&]( auto& cell )
        {
            u_ini[cell] = exact_solution( cell.center( 0. ), 0. );
        } );
    samurai::make_bc<samurai::Neumann<1>>( u_ini, 0. );

    // define problem ---------------------------------------------------------
    using state_t = decltype( u_ini );

    // diffusion terme
    auto diff = samurai::make_diffusion_order2<decltype( u_ini )>( d );
    auto fd   = [&]( double /* t */, state_t u, auto& du )
    {
        samurai::make_bc<samurai::Neumann<1>>( u, 0. );
        samurai::update_ghost_mr( u );
        du = -diff( u );
    };

    // reaction terme
    using cfg  = samurai::LocalCellSchemeConfig<samurai::SchemeType::NonLinear, 1, decltype( u_ini )>;
    auto react = samurai::make_cell_based_scheme<cfg>();
    react.set_name( "Reaction" );
    react.set_scheme_function(
        [&]( auto const& cell, auto const& field ) -> samurai::SchemeValue<cfg>
        {
            auto u = field[cell];
            return k * u * u * ( 1. - u );
        } );
    react.set_jacobian_function(
        [&]( auto const& cell, auto const& field ) -> samurai::JacobianMatrix<cfg>
        {
            auto u = field[cell];
            return k * ( 2. * u * ( 1. - u ) - u * u );
        } );
    auto fr_t = [&]( double /* t */ )
    {
        return react;
    };
    auto fr = [&]( double t, state_t u, auto& du )
    {
        samurai::make_bc<samurai::Neumann<1>>( u, 0. );
        samurai::update_ghost_mr( u );
        du = fr_t( t )( u );
    };

    auto eigmax_computer = [&mesh]( auto&, double, auto&, double, auto& )
    {
        double dx = min_cell_length( mesh );
        return 4. / ( dx * dx );
    };

    // clang-format off
    auto pirock_methods = std::make_tuple(
        ponio::runge_kutta::pirock::pirock<2>( // PIROCK, alpha=1, l=2, with Shampine trick
            ponio::runge_kutta::pirock::alpha_fixed<double>( 1.0 ),
            eigmax_computer,
            ponio::shampine_trick::shampine_trick<decltype( u_ini )>()
        ),
        ponio::runge_kutta::pirock::pirock<2>( // PIROCK, alpha=1, l=2
            ponio::runge_kutta::pirock::alpha_fixed<double>( 1.0 ),
            eigmax_computer
        ),
        ponio::runge_kutta::pirock::pirock<1>( // PIROCK, beta=0, l=1, with Shampine trick
            ponio::runge_kutta::pirock::beta_0<double>(),
            eigmax_computer,
            ponio::shampine_trick::shampine_trick<decltype( u_ini )>()
        ),
        ponio::runge_kutta::pirock::pirock<1>( // PIROCK, beta=0, l=1
            ponio::runge_kutta::pirock::beta_0<double>(),
            eigmax_computer
        )
    );
    // clang-format on

    static constexpr std::size_t N_methods = std::tuple_size<decltype( pirock_methods )>();

    // clang-format off
    auto pirock_param = std::array<std::string, N_methods>{
        "PIROCK, alpha=1, l=2, with Shampine trick",
        "PIROCK, alpha=1, l=2",
        "PIROCK, beta=0, l=1, with Shampine trick",
        "PIROCK, beta=0, l=1"
    };
    // clang-format on

    samurai::for_each_cell( mesh,
        [&]( auto& cell )
        {
            u_ini[cell] = exact_solution( cell.center( 0. ), 0. );
        } );
    auto pb = ponio::make_imex_operator_problem( fd, fr, fr_t );

    static_for<std::size_t, 0, N_methods>(
        [&]<std::size_t I>( std::integral_constant<std::size_t, I> )
        {
            std::vector<double> errors;
            std::vector<double> time_steps;

            for ( int N_iter = 45; N_iter < 120; N_iter += 10 )
            {
                double dt = ( t_end - t_ini ) / static_cast<double>( N_iter );

                // time loop  ---------------------------------------------------------
                auto sol_range = ponio::make_solver_range( pb, std::get<I>( pirock_methods ), u_ini, tspan, dt );

                auto it_sol = sol_range.begin();

                while ( it_sol->time < t_end )
                {
                    samurai::make_bc<samurai::Neumann<1>>( it_sol->state, 0. );
                    samurai::update_ghost_mr( it_sol->state );
                    it_sol.callback_on_stages(
                        []( auto& ki )
                        {
                            ki.resize();
                            ki.fill( 0. );
                        } );

                    ++it_sol;
                }

                // compute error
                double error = 0.;
                samurai::for_each_cell( mesh,
                    [&]( auto& cell )
                    {
                        error += std::abs( it_sol->state[cell] - exact_solution( cell.center( 0. ), t_end ) ) * cell.length;
                    } );
                errors.push_back( std::log( error ) );
                time_steps.push_back( std::log( dt ) );
            }

            auto [a, b] = mayor_method( time_steps, errors );

            INFO( "test order ", pirock_param[I] );
            CHECK( a >= doctest::Approx( 2 ).epsilon( 0.05 ) );
            WARN( a == doctest::Approx( 2 ).epsilon( 0.05 ) );
        } );
}
