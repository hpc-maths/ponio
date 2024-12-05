// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <tuple>

#include <doctest/doctest.h>

#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/splitting.hpp>

#include "compute_order.hpp"

enum struct class_method
{
    explicit_method,
    diagonal_implicit_method,
    exponential_method,
    additive_method,
    RDA_method,
    splitting_method
};

template <class_method type>
struct test_order
{
    template <typename rk_t>
    static void
    method_order()
    {
        if constexpr ( type == class_method::explicit_method )
        {
            INFO( "test order of ", rk_t::id );
            WARN( explicit_method::check_order( rk_t() ) == doctest::Approx( rk_t::order ).epsilon( 0.05 ) );
        }
        else if constexpr ( type == class_method::diagonal_implicit_method )
        {
            INFO( "test order of ", rk_t::id );
            WARN( diagonal_implicit_method::check_order( rk_t() ) == doctest::Approx( rk_t::order ).epsilon( 0.05 ) );
        }
        else if constexpr ( type == class_method::exponential_method )
        {
            INFO( "not implemented test for exponential method" );
            WARN( false );
        }
        else if constexpr ( type == class_method::additive_method )
        {
            // In additive Runge-Kutta method, one of method could be higher order than other (so we don't test equality)
            INFO( "test order of ", rk_t::id );
            WARN( additive_method::check_order( rk_t(), 0.5 ) >= doctest::Approx( rk_t::order ).epsilon( 0.05 ) );
            WARN( additive_method::check_order( rk_t(), 1. / 3. ) >= doctest::Approx( rk_t::order ).epsilon( 0.05 ) );
            WARN( additive_method::check_order( rk_t(), 2. / 3. ) >= doctest::Approx( rk_t::order ).epsilon( 0.05 ) );
            WARN( additive_method::check_order( rk_t(), 1. ) >= doctest::Approx( rk_t::order ).epsilon( 0.05 ) );
            WARN( additive_method::check_order( rk_t(), 0. ) >= doctest::Approx( rk_t::order ).epsilon( 0.05 ) );
        }
        else if constexpr ( type == class_method::RDA_method )
        {
            // In additive Runge-Kutta method, one of method could be higher order than other (so we don't test equality)
            INFO( "test order of ", rk_t::id );
            WARN( RDA_method::check_order( rk_t(), 1. ) >= doctest::Approx( rk_t::order ).epsilon( 0.05 ) );
            WARN( RDA_method::check_order( rk_t(), 2. / 3. ) >= doctest::Approx( rk_t::order ).epsilon( 0.05 ) );
            WARN( RDA_method::check_order( rk_t(), 1. / 3. ) >= doctest::Approx( rk_t::order ).epsilon( 0.05 ) );
            WARN( RDA_method::check_order( rk_t(), 0.5 ) >= doctest::Approx( rk_t::order ).epsilon( 0.05 ) );
            WARN( RDA_method::check_order( rk_t(), 0.25 ) >= doctest::Approx( rk_t::order ).epsilon( 0.05 ) );
            WARN( RDA_method::check_order( rk_t(), 0. ) >= doctest::Approx( rk_t::order ).epsilon( 0.05 ) );
        }
        else if constexpr ( type == class_method::splitting_method )
        {
            INFO( "test order of ", rk_t::id );
            WARN( splitting_method::check_order( rk_t() ) == doctest::Approx( rk_t::order ).epsilon( 0.125 ) );
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
        ( ( method_order<typename std::tuple_element<Is, rk_tuple>::type>() ), ... );
    }

    template <typename rk_tuple>
    static void
    on()
    {
        on_impl<rk_tuple>( std::make_index_sequence<std::tuple_size<rk_tuple>::value>() );
    }
};

TEST_CASE( "order::explict_runge_kutta" )
{
    test_order<class_method::explicit_method>::on<ponio::runge_kutta::erk_tuple<double>>();
}

TEST_CASE( "order::chebychev_runge_kutta" )
{
    auto rkc_methods = std::make_tuple( ponio::runge_kutta::chebyshev::explicit_rkc2<10>(),
        ponio::runge_kutta::rock::rock2<false>(),
        ponio::runge_kutta::rock::rock4<false>() );

    test_order<class_method::explicit_method>::on<decltype( rkc_methods )>();
}

TEST_CASE( "order::pirock" )
{
    auto pirock_methods = std::make_tuple( ponio::runge_kutta::pirock::pirock(),
        ponio::runge_kutta::pirock::pirock_a1(),
        ponio::runge_kutta::pirock::pirock_b0() );

    test_order<class_method::additive_method>::on<decltype( pirock_methods )>();
}

TEST_CASE( "order::pirock_RDA" )
{
    auto pirock_methods = std::make_tuple( ponio::runge_kutta::pirock::pirock_RDA(),
        ponio::runge_kutta::pirock::pirock_RDA_a1(),
        ponio::runge_kutta::pirock::pirock_RDA_b0() );

    test_order<class_method::RDA_method>::on<decltype( pirock_methods )>();
}

TEST_CASE( "order::splitting" )
{
    auto lie_splitting    = ponio::splitting::make_lie_tuple( std::make_pair( ponio::runge_kutta::rk_33(), .1 ),
        std::make_pair( ponio::runge_kutta::rk_33(), .1 ),
        std::make_pair( ponio::runge_kutta::rk_33(), .1 ) );
    auto strang_splitting = ponio::splitting::make_strang_tuple( std::make_pair( ponio::runge_kutta::rk_33(), .1 ),
        std::make_pair( ponio::runge_kutta::rk_44(), .1 ),
        std::make_pair( ponio::runge_kutta::rk_44(), .1 ) );

    WARN( splitting_method::check_order( lie_splitting ) == doctest::Approx( lie_splitting.order ).epsilon( 0.125 ) );
    WARN( splitting_method::check_order( strang_splitting ) == doctest::Approx( strang_splitting.order ).epsilon( 0.125 ) );
}

// TEST_CASE( "order::lawson_runge_kutta" )
// {
//     auto exp = [](double x){ return std::exp(x); };
//     test_order<class_method::exponential_method>::on<ponio::runge_kutta::lrk_tuple<double, decltype(exp)>>();
// }
