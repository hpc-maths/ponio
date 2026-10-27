// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <cstddef>
#include <valarray>

#include <doctest/doctest.h>

#include <ponio/problem.hpp>

namespace problem_test
{
    inline void
    check_state( std::valarray<double> const& u, std::valarray<double> const& expected )
    {
        REQUIRE( u.size() == expected.size() );
        for ( std::size_t i = 0; i < u.size(); ++i )
        {
            CHECK( u[i] == doctest::Approx( expected[i] ) );
        }
    }
} // namespace problem_test

TEST_CASE( "problem::sub_problems_share_the_state" )
{
    // each sub-problem receives the state, even when the state is an rvalue and a sub-problem takes it by value
    auto by_value = []( double, std::valarray<double> u, std::valarray<double>& du )
    {
        du = -1. * u;
    };
    auto by_const_ref = []( double, std::valarray<double> const& u, std::valarray<double>& du )
    {
        du = -2. * u;
    };

    std::valarray<double> const u = { 1., 2., 3. };
    std::valarray<double> du( u.size() );

    SUBCASE( "problem" )
    {
        auto pb = ponio::make_problem( by_value, by_const_ref, by_value );

        pb( 0., std::valarray<double>( u ), du );
        problem_test::check_state( du, -4. * u );

        std::valarray<double> v = u;
        pb( 0., v, du );
        problem_test::check_state( du, -4. * u );
        problem_test::check_state( v, u );
    }

    SUBCASE( "imex_problem" )
    {
        auto pb = ponio::make_imex_jacobian_problem( by_value,
            by_value,
            []( double, std::valarray<double> const& )
            {
                return 0.;
            } );

        pb( 0., std::valarray<double>( u ), du );
        problem_test::check_state( du, -2. * u );
    }

    SUBCASE( "lawson_problem" )
    {
        auto pb = ponio::make_lawson_problem( 3., by_value );

        pb( 0., std::valarray<double>( u ), du );
        problem_test::check_state( du, 2. * u );

        std::valarray<double> v = u;
        pb( 0., v, du );
        problem_test::check_state( du, 2. * u );
        problem_test::check_state( v, u );
    }
}
