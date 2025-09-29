// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <array>
#include <valarray>
#include <vector>

#include <ponio/expressions/state.hpp>

/////////////////////////////////////////////////////////////////////
// tests with a `std::array<double, N>`
// expr1: (a+b)^2 == a^2 + b^2 + 2ab
// expr2: (a-b)^2 == a^2 + b^2 - 2ab
// expr3: a^2 - b^2 == (a-b)(a+b)

// ------------------------------------------------------------------
// expr1: (a+b)^2 == a^2 + b^2 + 2ab

TEST_CASE( "expressions::array::expr1.0" )
{
    static constexpr std::size_t N = 5;
    using container_t              = std::array<double, N>;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    double result = 81;

    auto va  = ponio::expression::make_state( a );
    auto vb  = ponio::expression::make_state( b );
    auto two = ponio::expression::make_scalar( 2.0 );

    container_t r1, r2;
    auto vr1 = ponio::expression::make_state( r1 );
    auto vr2 = ponio::expression::make_state( r2 );

    vr1 = ( va + vb ) * ( va + vb );
    vr2 = va * va + vb * vb + two * va * vb;

    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result );
    }
}

TEST_CASE( "expressions::array::expr1.1" )
{
    static constexpr std::size_t N = 5;
    using container_t              = std::array<double, N>;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    double result = 81;

    auto va  = ponio::expression::make_state( a );
    auto vb  = ponio::expression::make_state( b );
    auto two = ponio::expression::make_scalar( 2.0 );

    auto r1 = ( va + vb ) * ( va + vb );
    auto r2 = va * va + vb * vb + two * va * vb;

    CHECK( r1.size() == r2.size() );
    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result );
    }
}

// ------------------------------------------------------------------
// expr2: (a-b)^2 == a^2 + b^2 - 2ab

TEST_CASE( "expressions::array::expr2.0" )
{
    static constexpr std::size_t N = 5;
    using container_t              = std::array<double, N>;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    container_t result = { 81, 49, 25, 9, 1 };

    auto va  = ponio::expression::make_state( a );
    auto vb  = ponio::expression::make_state( b );
    auto two = ponio::expression::make_scalar( 2.0 );

    container_t r1, r2;
    auto vr1 = ponio::expression::make_state( r1 );
    auto vr2 = ponio::expression::make_state( r2 );

    vr1 = ( va - vb ) * ( va - vb );
    vr2 = va * va + vb * vb - two * va * vb;

    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result[i] );
    }
}

TEST_CASE( "expressions::array::expr2.1" )
{
    static constexpr std::size_t N = 5;
    using container_t              = std::array<double, N>;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    container_t result = { 81, 49, 25, 9, 1 };

    auto va  = ponio::expression::make_state( a );
    auto vb  = ponio::expression::make_state( b );
    auto two = ponio::expression::make_scalar( 2.0 );

    auto r1 = ( va - vb ) * ( va - vb );
    auto r2 = va * va + vb * vb - two * va * vb;

    CHECK( r1.size() == r2.size() );
    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result[i] );
    }
}

// ------------------------------------------------------------------
// expr3: a^2 - b^2 == (a-b)(a+b)

TEST_CASE( "expressions::array::expr3.0" )
{
    static constexpr std::size_t N = 5;
    using container_t              = std::array<double, N>;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    container_t result = { -81, -63, -45, -27, -9 };

    auto va = ponio::expression::make_state( a );
    auto vb = ponio::expression::make_state( b );

    container_t r1, r2;
    auto vr1 = ponio::expression::make_state( r1 );
    auto vr2 = ponio::expression::make_state( r2 );

    vr1 = ( va - vb ) * ( va + vb );
    vr2 = va * va - vb * vb;

    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result[i] );
    }
}

TEST_CASE( "expressions::array::expr3.1" )
{
    static constexpr std::size_t N = 5;
    using container_t              = std::array<double, N>;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    container_t result = { -81, -63, -45, -27, -9 };

    auto va = ponio::expression::make_state( a );
    auto vb = ponio::expression::make_state( b );

    auto r1 = ( va - vb ) * ( va + vb );
    auto r2 = va * va - vb * vb;

    CHECK( r1.size() == r2.size() );
    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result[i] );
    }
}

/////////////////////////////////////////////////////////////////////
// tests with a `std::vector<double>`
// expr1: (a+b)^2 == a^2 + b^2 + 2ab
// expr2: (a-b)^2 == a^2 + b^2 - 2ab
// expr3: a^2 - b^2 == (a-b)(a+b)

// ------------------------------------------------------------------
// expr1: (a+b)^2 == a^2 + b^2 + 2ab

TEST_CASE( "expressions::vector::expr1.0" )
{
    static constexpr std::size_t N = 5;
    using container_t              = std::vector<double>;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    double result = 81;

    auto va  = ponio::expression::make_state( a );
    auto vb  = ponio::expression::make_state( b );
    auto two = ponio::expression::make_scalar( 2.0 );

    container_t r1( N ), r2( N );
    auto vr1 = ponio::expression::make_state( r1 );
    auto vr2 = ponio::expression::make_state( r2 );

    vr1 = ( va + vb ) * ( va + vb );
    vr2 = va * va + vb * vb + two * va * vb;

    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result );
    }
}

TEST_CASE( "expressions::vector::expr1.1" )
{
    using container_t = std::vector<double>;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    double result = 81;

    auto va  = ponio::expression::make_state( a );
    auto vb  = ponio::expression::make_state( b );
    auto two = ponio::expression::make_scalar( 2.0 );

    auto r1 = ( va + vb ) * ( va + vb );
    auto r2 = va * va + vb * vb + two * va * vb;

    CHECK( r1.size() == r2.size() );
    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result );
    }
}

// ------------------------------------------------------------------
// expr2: (a-b)^2 == a^2 + b^2 - 2ab

TEST_CASE( "expressions::vector::expr2.0" )
{
    static constexpr std::size_t N = 5;
    using container_t              = std::vector<double>;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    container_t result = { 81, 49, 25, 9, 1 };

    auto va  = ponio::expression::make_state( a );
    auto vb  = ponio::expression::make_state( b );
    auto two = ponio::expression::make_scalar( 2.0 );

    container_t r1( N ), r2( N );
    auto vr1 = ponio::expression::make_state( r1 );
    auto vr2 = ponio::expression::make_state( r2 );

    vr1 = ( va - vb ) * ( va - vb );
    vr2 = va * va + vb * vb - two * va * vb;

    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result[i] );
    }
}

TEST_CASE( "expressions::vector::expr2.1" )
{
    using container_t = std::vector<double>;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    container_t result = { 81, 49, 25, 9, 1 };

    auto va  = ponio::expression::make_state( a );
    auto vb  = ponio::expression::make_state( b );
    auto two = ponio::expression::make_scalar( 2.0 );

    auto r1 = ( va - vb ) * ( va - vb );
    auto r2 = va * va + vb * vb - two * va * vb;

    CHECK( r1.size() == r2.size() );
    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result[i] );
    }
}

// ------------------------------------------------------------------
// expr3: a^2 - b^2 == (a-b)(a+b)

TEST_CASE( "expressions::vector::expr3.0" )
{
    static constexpr std::size_t N = 5;
    using container_t              = std::vector<double>;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    container_t result = { -81, -63, -45, -27, -9 };

    auto va = ponio::expression::make_state( a );
    auto vb = ponio::expression::make_state( b );

    container_t r1( N ), r2( N );
    auto vr1 = ponio::expression::make_state( r1 );
    auto vr2 = ponio::expression::make_state( r2 );

    vr1 = ( va - vb ) * ( va + vb );
    vr2 = va * va - vb * vb;

    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result[i] );
    }
}

TEST_CASE( "expressions::vector::expr3.1" )
{
    using container_t = std::vector<double>;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    container_t result = { -81, -63, -45, -27, -9 };

    auto va = ponio::expression::make_state( a );
    auto vb = ponio::expression::make_state( b );

    auto r1 = ( va - vb ) * ( va + vb );
    auto r2 = va * va - vb * vb;

    CHECK( r1.size() == r2.size() );
    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result[i] );
    }
}
