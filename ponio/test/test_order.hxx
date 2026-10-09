// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <array>
#include <cmath>
#include <tuple>

#include <doctest/doctest.h>

#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>

#include "compute_order.hpp"

enum struct class_method
{
    explicit_method,
    diagonal_implicit_method,
    exponential_method,
    additive_method,
    RD_method,
    RDA_method,
    splitting_method
};

template <class_method type, typename exp_t = void>
struct test_order
{
    template <typename rk_t>
    static void
    method_order()
    {
        if constexpr ( type == class_method::explicit_method )
        {
            double computed_order;
            double computed_cst;

            // use std::tie because of a bug in clang++-15
            std::tie( computed_order, computed_cst ) = explicit_method::check_order( rk_t() );

            INFO( "test order of ", rk_t::id );
            INFO( "theoretical order: ", rk_t::order );
            INFO( "computed order   : ", computed_order );
            INFO( "computed error constant: ", computed_cst );

            if ( computed_cst > -8. ) // error is to close than computer error
            {
                CHECK( computed_order >= doctest::Approx( rk_t::order ).epsilon( 0.05 ) );
                WARN( computed_order == doctest::Approx( rk_t::order ).epsilon( 0.05 ) );
            }
            else
            {
                WARN( computed_order == doctest::Approx( rk_t::order ).epsilon( 2. ) );
            }
        }
        else if constexpr ( type == class_method::diagonal_implicit_method )
        {
            using dirk_t = decltype( std::declval<rk_t>()() );

            double computed_order;
            double computed_cst;

            // use std::tie because of a bug in clang++-15
            std::tie( computed_order, computed_cst ) = diagonal_implicit_method::check_order( dirk_t() );

            INFO( "test order of ", dirk_t::id );
            INFO( "theoretical order: ", dirk_t::order );
            INFO( "computed order   : ", computed_order );
            INFO( "computed error constant: ", computed_cst );

            if ( computed_cst > -8. ) // error is to close than computer error
            {
                CHECK( computed_order >= doctest::Approx( dirk_t::order ).epsilon( 0.05 ) );
                WARN( computed_order == doctest::Approx( dirk_t::order ).epsilon( 0.05 ) );
            }
            else
            {
                WARN( computed_order == doctest::Approx( dirk_t::order ).epsilon( 2. ) );
            }
        }
        else if constexpr ( type == class_method::exponential_method )
        {
            using exprk_t = decltype( std::declval<rk_t>()( exp_t() ) );

            // In exponential Runge-Kutta method when coefficient lambda is equal to 1 we get exact solution so we don't test lambda=1.
            for ( auto lambda : { 0.5, 1. / 3., 2. / 3., 0. } )
            {
                double computed_order;
                double computed_cst;

                // use std::tie because of a bug in clang++-15
                std::tie( computed_order, computed_cst ) = exponential_method::check_order( exprk_t( exp_t() ), lambda );

                INFO( "test order of ", exprk_t::id );
                INFO( "lambda: ", lambda );
                INFO( "theoretical order: ", exprk_t::order );
                INFO( "computed order   : ", computed_order );
                INFO( "computed error constant: ", computed_cst );

                if ( computed_cst > -8. ) // error is to close than computer error
                {
                    CHECK( computed_order >= doctest::Approx( exprk_t::order ).epsilon( 0.05 ) );
                    WARN( computed_order == doctest::Approx( exprk_t::order ).epsilon( 0.05 ) );
                }
                else
                {
                    WARN( computed_order == doctest::Approx( exprk_t::order ).epsilon( 2. ) );
                }
            }
        }
        else if constexpr ( type == class_method::additive_method )
        {
            using ark_t = decltype( std::declval<rk_t>()() );

            // In additive Runge-Kutta method, one of method could be higher order than other (so we don't test equality)

            // In additive Runge-Kutta method, one of method could be higher order than other (so we don't test equality)
            for ( auto lambda : { 0.5, 1. / 3., 2. / 3., 0., 1. } )
            {
                double computed_order;
                double computed_cst;

                // use std::tie because of a bug in clang++-15
                std::tie( computed_order, computed_cst ) = additive_method::check_order( ark_t(), lambda );

                INFO( "test order of ", ark_t::id );
                INFO( "lambda: ", lambda );
                INFO( "theoretical order: ", ark_t::order );
                INFO( "computed order   : ", computed_order );
                INFO( "computed error constant: ", computed_cst );

                if ( computed_cst > -8. ) // error is to close than computer error
                {
                    WARN( computed_order >= doctest::Approx( ark_t::order ).epsilon( 0.05 ) );
                }
                else
                {
                    WARN( computed_order >= doctest::Approx( ark_t::order ).epsilon( 2. ) );
                }
            }
        }
        else if constexpr ( type == class_method::RD_method )
        {
            // In additive Runge-Kutta method, one of method could be higher order than other (so we don't test equality)
            for ( auto lambda : { 0.5, 1. / 3., 2. / 3., 0., 1. } )
            {
                double computed_order;
                double computed_cst;

                // use std::tie because of a bug in clang++-15
                std::tie( computed_order, computed_cst ) = RD_method::check_order( rk_t(), lambda );

                INFO( "test order of ", rk_t::id );
                INFO( "lambda: ", lambda );
                INFO( "theoretical order: ", rk_t::order );
                INFO( "computed order   : ", computed_order );
                INFO( "computed error constant: ", computed_cst );

                if ( computed_cst > -8. ) // error is to close than computer error
                {
                    WARN( computed_order >= doctest::Approx( rk_t::order ).epsilon( 0.05 ) );
                }
                else
                {
                    WARN( computed_order >= doctest::Approx( rk_t::order ).epsilon( 2. ) );
                }
            }
        }
        else if constexpr ( type == class_method::RDA_method )
        {
            // In additive Runge-Kutta method, one of method could be higher order than other (so we don't test equality)
            for ( auto lambda : { 0.5, 1. / 3., 2. / 3., 0., 1. } )
            {
                double computed_order;
                double computed_cst;

                // use std::tie because of a bug in clang++-15
                std::tie( computed_order, computed_cst ) = RDA_method::check_order( rk_t(), lambda );

                INFO( "test order of ", rk_t::id );
                INFO( "lambda: ", lambda );
                INFO( "theoretical order: ", rk_t::order );
                INFO( "computed order   : ", computed_order );
                INFO( "computed error constant: ", computed_cst );

                if ( computed_cst > -8. ) // error is to close than computer error
                {
                    WARN( computed_order >= doctest::Approx( rk_t::order ).epsilon( 0.05 ) );
                }
                else
                {
                    WARN( computed_order >= doctest::Approx( rk_t::order ).epsilon( 2. ) );
                }
            }
        }
        else if constexpr ( type == class_method::splitting_method )
        {
            double computed_order;
            double computed_cst;

            // use std::tie because of a bug in clang++-15
            std::tie( computed_order, computed_cst ) = splitting_method::check_order( rk_t() );

            INFO( "test order of ", rk_t::id );
            INFO( "theoretical order: ", rk_t::order );
            INFO( "computed order   : ", computed_order );
            INFO( "computed error constant: ", computed_cst );

            if ( computed_cst > -8. ) // error is to close than computer error
            {
                CHECK( computed_order >= doctest::Approx( rk_t::order ).epsilon( 0.125 ) );
                WARN( computed_order == doctest::Approx( rk_t::order ).epsilon( 0.125 ) );
            }
            else
            {
                WARN( computed_order == doctest::Approx( rk_t::order ).epsilon( 2. ) );
            }
        }
        else
        {
            INFO( "not implemented test for unknown method" );
            WARN( false );
        }
    }

    template <typename rk_tuple, std::size_t... Is>
    static void
    on_impl( std::index_sequence<Is...> )
    {
        ( ( method_order<typename std::tuple_element_t<Is, rk_tuple>>() ), ... );
    }

    template <typename rk_tuple>
    static void
    on()
    {
        on_impl<rk_tuple>( std::make_index_sequence<std::tuple_size_v<rk_tuple>>() );
    }
};

TEST_CASE( "order::explict_runge_kutta" )
{
    test_order<class_method::explicit_method>::on<ponio::runge_kutta::erk_tuple<double>>();
}

TEST_CASE( "order::diagonal_implicit_runge_kutta" )
{
    test_order<class_method::diagonal_implicit_method>::on<ponio::runge_kutta::dirk_tuple<double>>();
}

TEST_CASE( "order::chebychev_runge_kutta" )
{
    // clang-format off
    using rkc_methods = std::tuple<
        decltype( ponio::runge_kutta::chebyshev::explicit_rkc2<10>() ),
        decltype( ponio::runge_kutta::rock::rock2<false>() ),
        decltype( ponio::runge_kutta::rock::rock4<false>() )
    >;
    // clang-format on

    test_order<class_method::explicit_method>::on<rkc_methods>();
}

TEST_CASE( "order::legendre_runge_kutta" )
{
    // clang-format off
    using rkl_methods = std::tuple<
        decltype( ponio::runge_kutta::legendre::explicit_rkl2<10>() ),
        decltype( ponio::runge_kutta::legendre::explicit_rkl2<5>() ),
        decltype( ponio::runge_kutta::legendre::explicit_rkl1<10>() )
    >;
    // clang-format on

    test_order<class_method::explicit_method>::on<rkl_methods>();
}

TEST_CASE( "order::pirock_RD" )
{
    // clang-format off
    using pirock_methods = std::tuple<
        decltype( ponio::runge_kutta::pirock::pirock() ),
        decltype( ponio::runge_kutta::pirock::pirock_a1() ),
        decltype( ponio::runge_kutta::pirock::pirock_b0() )
    >;
    // clang-format on

    test_order<class_method::RD_method>::on<pirock_methods>();
}

TEST_CASE( "order::pirock_DA" )
{
    constexpr double lambda_D = -4.;
    constexpr double lambda_A = 1.;
    constexpr double rho_D    = -lambda_D;
    constexpr double t_end    = 1.;

    auto diffusion = []( double, double const& u, double& du )
    {
        du = lambda_D * u;
    };

    auto advection = []( double, double const& u, double& du )
    {
        du = lambda_A * u;
    };

    auto spectral_radius = [=]( auto&&, double, auto const&, double, auto& )
    {
        return rho_D;
    };

    auto problem = ponio::make_problem( diffusion, advection );

    auto check_method = [&]( auto method )
    {
        auto error = [&]( double dt )
        {
            using method_t = std::decay_t<decltype( method )>;

            std::array<double, method_t::N_storage> work{};

            double t      = 0.;
            double u      = 1.;
            double u_next = 0.;

            auto const n_steps = static_cast<std::size_t>( std::lround( t_end / dt ) );

            for ( std::size_t n = 0; n < n_steps; ++n )
            {
                method( problem, t, u, work, dt, u_next );
                u = u_next;
            }

            return std::abs( u - std::exp( ( lambda_D + lambda_A ) * t_end ) );
        };

        double const err_dt   = error( 1. / 20. );
        double const err_dt_2 = error( 1. / 40. );
        double const order    = std::log2( err_dt / err_dt_2 );

        INFO( "test order of PIROCK-DA" );
        INFO( "computed order: ", order );

        CHECK( order >= doctest::Approx( 2. ).epsilon( 0.1 ) );
    };

    check_method( ponio::runge_kutta::pirock::pirock_DA( spectral_radius ) );
    check_method( ponio::runge_kutta::pirock::pirock_DA_a1( spectral_radius ) );
    check_method( ponio::runge_kutta::pirock::pirock_DA_b0( spectral_radius ) );
}

TEST_CASE( "order::pirock_RDA" )
{
    // clang-format off
    using pirock_methods = std::tuple<
        decltype( ponio::runge_kutta::pirock::pirock_RDA() ),
        decltype( ponio::runge_kutta::pirock::pirock_RDA_a1() ),
        decltype( ponio::runge_kutta::pirock::pirock_RDA_b0() )
    >;
    // clang-format on

    test_order<class_method::RDA_method>::on<pirock_methods>();
}

TEST_CASE( "order::splitting[odd]" )
{
    // for Lie splitting method number of sub-stages can be odd or even (for Strang is always odd)

    // clang-format off
    auto lie_splitting    = ponio::splitting::make_lie_tuple(
        std::make_pair( ponio::runge_kutta::rk_33(), .1 ),
        std::make_pair( ponio::runge_kutta::rk_33(), .1 ),
        std::make_pair( ponio::runge_kutta::rk_33(), .1 )
    );
    auto strang_splitting = ponio::splitting::make_strang_tuple(
        std::make_pair( ponio::runge_kutta::rk_33(), .1 ),
        std::make_pair( ponio::runge_kutta::rk_44(), .1 ),
        std::make_pair( ponio::runge_kutta::rk_44(), .1 )
    );
    // clang-format on

    // Lie splitting method
    double lie_computed_order;
    double lie_computed_cst;
    // use std::tie because of a bug in clang++-15
    std::tie( lie_computed_order, lie_computed_cst ) = splitting_method::check_order( lie_splitting );

    INFO( "test order of ", lie_splitting.id );
    INFO( "theoretical order: ", lie_splitting.order );
    INFO( "computed order   : ", lie_computed_order );
    INFO( "computed error constant: ", lie_computed_cst );

    if ( lie_computed_cst > -8. ) // error is to close than computer error
    {
        CHECK( lie_computed_order >= doctest::Approx( lie_splitting.order ).epsilon( 0.125 ) );
        WARN( lie_computed_order == doctest::Approx( lie_splitting.order ).epsilon( 0.125 ) );
    }
    else
    {
        WARN( lie_computed_order == doctest::Approx( lie_splitting.order ).epsilon( 2. ) );
    }

    // Strang splitting method
    double strang_computed_order;
    double strang_computed_cst;
    // use std::tie because of a bug in clang++-15
    std::tie( strang_computed_order, strang_computed_cst ) = splitting_method::check_order( strang_splitting );

    INFO( "test order of ", strang_splitting.id );
    INFO( "theoretical order: ", strang_splitting.order );
    INFO( "computed order   : ", strang_computed_order );
    INFO( "computed error constant: ", strang_computed_cst );

    if ( strang_computed_cst > -8. ) // error is to close than computer error
    {
        CHECK( strang_computed_order >= doctest::Approx( strang_splitting.order ).epsilon( 0.125 ) );
        WARN( strang_computed_order == doctest::Approx( strang_splitting.order ).epsilon( 0.125 ) );
    }
    else
    {
        WARN( strang_computed_order == doctest::Approx( strang_splitting.order ).epsilon( 2. ) );
    }
}

TEST_CASE( "order::splitting[even]" )
{
    // for Lie splitting method number of sub-stages can be odd or even (for Strang is always odd)

    // clang-format off
    auto lie_splitting    = ponio::splitting::make_lie_tuple(
        std::make_pair( ponio::runge_kutta::rk_44(), .001 ),
        std::make_pair( ponio::runge_kutta::rk_33(), .001 )
    );
    auto strang_splitting = ponio::splitting::make_strang_tuple(
        std::make_pair( ponio::runge_kutta::rk_33(), .1 ),
        std::make_pair( ponio::runge_kutta::rk_44(), .1 )
    );
    // clang-format on

    // Lie splitting method
    double lie_computed_order;
    double lie_computed_cst;

    // use std::tie because of a bug in clang++-15
    std::tie( lie_computed_order, lie_computed_cst ) = splitting_method::check_order( lie_splitting );

    INFO( "test order of ", lie_splitting.id );
    INFO( "theoretical order: ", lie_splitting.order );
    INFO( "computed order   : ", lie_computed_order );
    INFO( "computed error constant: ", lie_computed_cst );

    if ( lie_computed_cst > -8. ) // error is to close than computer error
    {
        CHECK( lie_computed_order >= doctest::Approx( lie_splitting.order ).epsilon( 0.125 ) );
        WARN( lie_computed_order == doctest::Approx( lie_splitting.order ).epsilon( 0.125 ) );
    }
    else
    {
        WARN( lie_computed_order == doctest::Approx( lie_splitting.order ).epsilon( 2. ) );
    }

    // Strang splitting method
    double strang_computed_order;
    double strang_computed_cst;

    // use std::tie because of a bug in clang++-15
    std::tie( strang_computed_order, strang_computed_cst ) = splitting_method::check_order( strang_splitting );

    INFO( "test order of ", strang_splitting.id );
    INFO( "theoretical order: ", strang_splitting.order );
    INFO( "computed order   : ", strang_computed_order );
    INFO( "computed error constant: ", strang_computed_cst );

    if ( strang_computed_cst > -8. ) // error is to close than computer error
    {
        CHECK( strang_computed_order >= doctest::Approx( strang_splitting.order ).epsilon( 0.125 ) );
        WARN( strang_computed_order == doctest::Approx( strang_splitting.order ).epsilon( 0.125 ) );
    }
    else
    {
        WARN( strang_computed_order == doctest::Approx( strang_splitting.order ).epsilon( 2. ) );
    }
}

TEST_CASE( "order::splitting::_split_solve" )
{
    // test the implementation of `ponio::splitting::detail::_split_solve` function

    auto lie_splitting = ponio::splitting::make_lie_tuple( std::make_pair( ponio::runge_kutta::rk_33(), .1 ),
        std::make_pair( ponio::runge_kutta::rk_44(), .1 ),
        std::make_pair( ponio::runge_kutta::rk_22_midpoint(), .1 ) );

    using algos_t = decltype( lie_splitting.algos );

    double y      = 0.;
    auto lie_meth = ponio::make_method<double>( lie_splitting, y );

    static constexpr std::size_t I = 0;
    static constexpr std::size_t J = 1;
    static constexpr std::size_t K = 2;

    double I_computed_order;
    double I_computed_cst;
    std::tie( I_computed_order, I_computed_cst ) = splitting_method::check_order_split_solve<I>( lie_meth );

    double J_computed_order;
    double J_computed_cst;
    std::tie( J_computed_order, J_computed_cst ) = splitting_method::check_order_split_solve<J>( lie_meth );

    double K_computed_order;
    double K_computed_cst;
    std::tie( K_computed_order, K_computed_cst ) = splitting_method::check_order_split_solve<K>( lie_meth );

    WARN( I_computed_order == doctest::Approx( std::tuple_element_t<I, algos_t>::order ).epsilon( 0.125 ) );
    WARN( J_computed_order == doctest::Approx( std::tuple_element_t<J, algos_t>::order ).epsilon( 0.125 ) );
    WARN( K_computed_order == doctest::Approx( std::tuple_element_t<K, algos_t>::order ).epsilon( 0.125 ) );
}

TEST_CASE( "order::additive_runge_kutta " )
{
    test_order<class_method::additive_method>::on<ponio::runge_kutta::ark_tuple<double>>();
}

// TEST_CASE( "order::lawson_runge_kutta" )
// {
//     auto exp = []( double x )
//     {
//         return std::exp( x );
//     };
//     using exp_t = decltype( exp );
//     test_order<class_method::exponential_method, exp_t>::on<ponio::runge_kutta::lrk_tuple<double, exp_t>>();
// }
