// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <array>
#include <cmath>
#include <ponio/detail.hpp>

TEST_CASE( "detail::power" )
{
    CHECK( ponio::detail::power<2>( 16 ) == 256 );
    CHECK( ponio::detail::power<2>( std::sqrt( 2. ) ) == doctest::Approx( 2.0 ).epsilon( 1e-15 ) );
    CHECK( ponio::detail::power<10>( 2 ) == 1024 );
}

TEST_CASE( "detail::init_fill_array" )
{
    SUBCASE( "with a value" )
    {
        std::array<int, 5> const a1 = { 42, 42, 42, 42, 42 };
        auto const a2               = ponio::detail::init_fill_array<5>( 42 );

        REQUIRE( a1.size() == a2.size() );
        for ( auto i = 0u; i < a1.size(); ++i )
        {
            CHECK( a1[i] == a2[i] );
        }
    }
    SUBCASE( "with a function" )
    {
        std::array<int, 5> const a1 = { 0, 1, 4, 9, 16 };
        auto const a2               = ponio::detail::init_fill_array<5>(
            [count = 0]( int i ) mutable
            {
                auto tmp = i + count * count;
                ++count;
                return tmp;
            },
            0 );

        REQUIRE( a1.size() == a2.size() );
        for ( auto i = 0u; i < a1.size(); ++i )
        {
            CHECK( a1[i] == a2[i] );
        }
    }
}

TEST_CASE( "detail::tpl_inner_product" )
{
    constexpr std::size_t N          = 10;
    std::array<int, N + 1> const tab = ponio::detail::init_fill_array<N + 1>(
        [count = 0]( int i ) mutable
        {
            return i + count++;
        },
        0 );

    CHECK( ponio::detail::tpl_inner_product<N + 1>( tab, tab, 0, 1 ) == ( N * ( N + 1 ) * ( 2 * N + 1 ) / 6 ) );

    constexpr std::size_t M = 5;
    CHECK( ponio::detail::tpl_inner_product<M + 1>( tab, tab, 0, 1 ) == ( M * ( M + 1 ) * ( 2 * M + 1 ) / 6 ) );
}
