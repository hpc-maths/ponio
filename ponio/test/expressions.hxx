// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <array>
#include <span>
#include <vector>

#include <ponio/expressions/state.hpp>

/*

Tests on expressions
====================

In this test suite we compute 3 classical expressions with different kind of containers. Expressions are (with :math:`a` and :math:`b` two
containers):

1. :math:`(a+b)^2` and test equality with :math:`a^2 + b^2 + 2ab`;
2. :math:`(a-b)^2` and test equality with :math:`a^b + b^2 - 2ab`;
3. :math:`a^2 - b^2` and test equality with :math:`(a-b)(a+b)`.

For each expression inputs are:

.. math::

    a = (0, 1, 2, 3, 4), \qquad b = (9, 8, 7, 6, 5)

and outputs are:

* For the first expression, the result of :math:`(a+b)^2 = 81`;
* For the second expression, the result of :math:`(a-b)^2 = ( 81, 49, 25, 9, 1)`;
* For the third expression, the result of :math:`a^2 - b^2 = (-81, -63, -45, -27, -9)`.

For each expression we test 3 ways to write expression with ponio expressions:

0. Create a container for result of expressions, compute expression with :cpp:`ponio::expression::state` objects on each container and test
result on output container.
1. Get only the expression object and test result on this object (no storing of results), compute expression with
:cpp:`ponio::expression::state` objects on each container.
2. Create a container for result of expressions, compute expression with only r-value :cpp:`ponio::expression::state` objects and test
result on output container.
3. Get only the expression object and test result on this object (no storing of results), compute expression with only r-value
:cpp:`ponio::expression::state` objects and test result on output container.

.. note::

    That why all tests on expressions are suffixed by ``exprX.Y` this is for expression ``X`` writing with rule ``Y``.

We test this with 3 containers:

* :cpp:`std::array`: container with static size;
* :cpp:`std::vector`: container with dynamic size;
* :cpp:`std::span`: container with dynamic size and without ownership of data.

*/

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

TEST_CASE( "expressions::array::expr1.2" )
{
    static constexpr std::size_t N = 5;
    using container_t              = std::array<double, N>;
    using namespace ponio::expression;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    container_t r1, r2;
    double result = 81;

    make_state( r1 ) = ( make_state( a ) + make_state( b ) ) * ( make_state( a ) + make_state( b ) );
    make_state( r2 ) = make_state( a ) * make_state( a ) + make_state( b ) * make_state( b )
                     + make_scalar( 2.0 ) * make_state( a ) * make_state( b );

    CHECK( r1.size() == r2.size() );
    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result );
    }
}

TEST_CASE( "expressions::array::expr1.3" )
{
    static constexpr std::size_t N = 5;
    using container_t              = std::array<double, N>;
    using namespace ponio::expression;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    double result = 81;

    auto r1 = ( make_state( a ) + make_state( b ) ) * ( make_state( a ) + make_state( b ) );
    auto r2 = make_state( a ) * make_state( a ) + make_state( b ) * make_state( b ) + make_scalar( 2.0 ) * make_state( a ) * make_state( b );

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

TEST_CASE( "expressions::array::expr2.2" )
{
    static constexpr std::size_t N = 5;
    using container_t              = std::array<double, N>;
    using namespace ponio::expression;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    container_t r1, r2;
    container_t result = { 81, 49, 25, 9, 1 };

    make_state( r1 ) = ( make_state( a ) - make_state( b ) ) * ( make_state( a ) - make_state( b ) );
    make_state( r2 ) = make_state( a ) * make_state( a ) + make_state( b ) * make_state( b )
                     - make_scalar( 2.0 ) * make_state( a ) * make_state( b );

    CHECK( r1.size() == r2.size() );
    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result[i] );
    }
}

TEST_CASE( "expressions::array::expr2.3" )
{
    static constexpr std::size_t N = 5;
    using container_t              = std::array<double, N>;
    using namespace ponio::expression;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    container_t result = { 81, 49, 25, 9, 1 };

    auto r1 = ( make_state( a ) - make_state( b ) ) * ( make_state( a ) - make_state( b ) );
    auto r2 = make_state( a ) * make_state( a ) + make_state( b ) * make_state( b ) - make_scalar( 2.0 ) * make_state( a ) * make_state( b );

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

TEST_CASE( "expressions::array::expr3.2" )
{
    static constexpr std::size_t N = 5;
    using container_t              = std::array<double, N>;
    using namespace ponio::expression;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    container_t r1, r2;
    container_t result = { -81, -63, -45, -27, -9 };

    make_state( r1 ) = ( make_state( a ) - make_state( b ) ) * ( make_state( a ) + make_state( b ) );
    make_state( r2 ) = make_state( a ) * make_state( a ) - make_state( b ) * make_state( b );

    CHECK( r1.size() == r2.size() );
    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result[i] );
    }
}

TEST_CASE( "expressions::array::expr3.3" )
{
    static constexpr std::size_t N = 5;
    using container_t              = std::array<double, N>;
    using namespace ponio::expression;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    container_t result = { -81, -63, -45, -27, -9 };

    auto r1 = ( make_state( a ) - make_state( b ) ) * ( make_state( a ) + make_state( b ) );
    auto r2 = make_state( a ) * make_state( a ) - make_state( b ) * make_state( b );

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

TEST_CASE( "expressions::vector::expr1.2" )
{
    using container_t = std::vector<double>;
    using namespace ponio::expression;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    container_t r1, r2;
    double result = 81;

    r1.reserve( a.size() );
    r2.reserve( a.size() );

    make_state( r1 ) = ( make_state( a ) + make_state( b ) ) * ( make_state( a ) + make_state( b ) );
    make_state( r2 ) = make_state( a ) * make_state( a ) + make_state( b ) * make_state( b )
                     + make_scalar( 2.0 ) * make_state( a ) * make_state( b );

    CHECK( r1.size() == r2.size() );
    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result );
    }
}

TEST_CASE( "expressions::vector::expr1.3" )
{
    using container_t = std::vector<double>;
    using namespace ponio::expression;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    double result = 81;

    auto r1 = ( make_state( a ) + make_state( b ) ) * ( make_state( a ) + make_state( b ) );
    auto r2 = make_state( a ) * make_state( a ) + make_state( b ) * make_state( b ) + make_scalar( 2.0 ) * make_state( a ) * make_state( b );

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

TEST_CASE( "expressions::vector::expr2.2" )
{
    using container_t = std::vector<double>;
    using namespace ponio::expression;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    container_t r1, r2;
    container_t result = { 81, 49, 25, 9, 1 };

    r1.reserve( a.size() );
    r2.reserve( a.size() );

    make_state( r1 ) = ( make_state( a ) - make_state( b ) ) * ( make_state( a ) - make_state( b ) );
    make_state( r2 ) = make_state( a ) * make_state( a ) + make_state( b ) * make_state( b )
                     - make_scalar( 2.0 ) * make_state( a ) * make_state( b );

    CHECK( r1.size() == r2.size() );
    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result[i] );
    }
}

TEST_CASE( "expressions::vector::expr2.3" )
{
    using container_t = std::vector<double>;
    using namespace ponio::expression;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    container_t result = { 81, 49, 25, 9, 1 };

    auto r1 = ( make_state( a ) - make_state( b ) ) * ( make_state( a ) - make_state( b ) );
    auto r2 = make_state( a ) * make_state( a ) + make_state( b ) * make_state( b ) - make_scalar( 2.0 ) * make_state( a ) * make_state( b );

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

TEST_CASE( "expressions::vector::expr3.2" )
{
    using container_t = std::vector<double>;
    using namespace ponio::expression;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    container_t r1, r2;
    container_t result = { -81, -63, -45, -27, -9 };

    r1.reserve( a.size() );
    r2.reserve( a.size() );

    make_state( r1 ) = ( make_state( a ) - make_state( b ) ) * ( make_state( a ) + make_state( b ) );
    make_state( r2 ) = make_state( a ) * make_state( a ) - make_state( b ) * make_state( b );

    CHECK( r1.size() == r2.size() );
    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result[i] );
    }
}

TEST_CASE( "expressions::vector::expr3.3" )
{
    using container_t = std::vector<double>;
    using namespace ponio::expression;

    container_t a{ 0., 1., 2., 3., 4. }, b{ 9., 8., 7., 6., 5. };
    container_t result = { -81, -63, -45, -27, -9 };

    auto r1 = ( make_state( a ) - make_state( b ) ) * ( make_state( a ) + make_state( b ) );
    auto r2 = make_state( a ) * make_state( a ) - make_state( b ) * make_state( b );

    CHECK( r1.size() == r2.size() );
    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result[i] );
    }
}

/////////////////////////////////////////////////////////////////////
// tests with a `std::span<double>`
// expr1: (a+b)^2 == a^2 + b^2 + 2ab
// expr2: (a-b)^2 == a^2 + b^2 - 2ab
// expr3: a^2 - b^2 == (a-b)(a+b)

// ------------------------------------------------------------------
// expr1: (a+b)^2 == a^2 + b^2 + 2ab

TEST_CASE( "expressions::span::expr1.0" )
{
    using raw_container_t = std::vector<double>;
    using container_t     = std::span<double>;

    raw_container_t raw_a{ 0., 1., 2., 3., 4. }, raw_b{ 9., 8., 7., 6., 5. };
    container_t a( raw_a.data(), raw_a.size() ), b( raw_b.data(), raw_b.size() );
    double result = 81;

    auto va  = ponio::expression::make_state( a );
    auto vb  = ponio::expression::make_state( b );
    auto two = ponio::expression::make_scalar( 2.0 );

    raw_container_t raw_r1( a.size() ), raw_r2( a.size() );
    container_t r1( raw_r1.data(), raw_r1.size() ), r2( raw_r2.data(), raw_r2.size() );
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

TEST_CASE( "expressions::span::expr1.1" )
{
    using raw_container_t = std::vector<double>;
    using container_t     = std::span<double>;

    raw_container_t raw_a{ 0., 1., 2., 3., 4. }, raw_b{ 9., 8., 7., 6., 5. };
    container_t a( raw_a.data(), raw_a.size() ), b( raw_b.data(), raw_b.size() );
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

TEST_CASE( "expressions::span::expr1.2" )
{
    using raw_container_t = std::vector<double>;
    using container_t     = std::span<double>;
    using namespace ponio::expression;

    raw_container_t raw_a{ 0., 1., 2., 3., 4. }, raw_b{ 9., 8., 7., 6., 5. };
    raw_container_t raw_r1, raw_r2;
    raw_r1.reserve( raw_a.size() );
    raw_r2.reserve( raw_a.size() );

    container_t a( raw_a.data(), raw_a.size() ), b( raw_b.data(), raw_b.size() );
    container_t r1( raw_r1.data(), raw_r1.size() ), r2( raw_r2.data(), raw_r2.size() );
    double result = 81;

    make_state( r1 ) = ( make_state( a ) + make_state( b ) ) * ( make_state( a ) + make_state( b ) );
    make_state( r2 ) = make_state( a ) * make_state( a ) + make_state( b ) * make_state( b )
                     + make_scalar( 2.0 ) * make_state( a ) * make_state( b );

    CHECK( r1.size() == r2.size() );
    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result );
    }
}

TEST_CASE( "expressions::span::expr1.3" )
{
    using raw_container_t = std::vector<double>;
    using container_t     = std::span<double>;
    using namespace ponio::expression;

    raw_container_t raw_a{ 0., 1., 2., 3., 4. }, raw_b{ 9., 8., 7., 6., 5. };

    container_t a( raw_a.data(), raw_a.size() ), b( raw_b.data(), raw_b.size() );
    double result = 81;

    auto r1 = ( make_state( a ) + make_state( b ) ) * ( make_state( a ) + make_state( b ) );
    auto r2 = make_state( a ) * make_state( a ) + make_state( b ) * make_state( b ) + make_scalar( 2.0 ) * make_state( a ) * make_state( b );

    CHECK( r1.size() == r2.size() );
    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result );
    }
}

// ------------------------------------------------------------------
// expr2: (a-b)^2 == a^2 + b^2 - 2ab

TEST_CASE( "expressions::span::expr2.0" )
{
    using raw_container_t = std::vector<double>;
    using container_t     = std::span<double>;

    raw_container_t raw_a{ 0., 1., 2., 3., 4. }, raw_b{ 9., 8., 7., 6., 5. };
    container_t a( raw_a.data(), raw_a.size() ), b( raw_b.data(), raw_b.size() );
    raw_container_t result{ 81, 49, 25, 9, 1 };

    auto va  = ponio::expression::make_state( a );
    auto vb  = ponio::expression::make_state( b );
    auto two = ponio::expression::make_scalar( 2.0 );

    raw_container_t raw_r1( a.size() ), raw_r2( a.size() );
    container_t r1( raw_r1.data(), raw_r1.size() ), r2( raw_r2.data(), raw_r2.size() );
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

TEST_CASE( "expressions::span::expr2.1" )
{
    using raw_container_t = std::vector<double>;
    using container_t     = std::span<double>;

    raw_container_t raw_a{ 0., 1., 2., 3., 4. }, raw_b{ 9., 8., 7., 6., 5. };
    container_t a( raw_a.data(), raw_a.size() ), b( raw_b.data(), raw_b.size() );
    raw_container_t result{ 81, 49, 25, 9, 1 };

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

TEST_CASE( "expressions::span::expr2.2" )
{
    using raw_container_t = std::vector<double>;
    using container_t     = std::span<double>;
    using namespace ponio::expression;

    raw_container_t raw_a{ 0., 1., 2., 3., 4. }, raw_b{ 9., 8., 7., 6., 5. };
    raw_container_t raw_r1, raw_r2;
    raw_r1.reserve( raw_a.size() );
    raw_r2.reserve( raw_a.size() );

    container_t a( raw_a.data(), raw_a.size() ), b( raw_b.data(), raw_b.size() );
    container_t r1( raw_r1.data(), raw_r1.size() ), r2( raw_r2.data(), raw_r2.size() );
    raw_container_t result{ 81, 49, 25, 9, 1 };

    make_state( r1 ) = ( make_state( a ) - make_state( b ) ) * ( make_state( a ) - make_state( b ) );
    make_state( r2 ) = make_state( a ) * make_state( a ) + make_state( b ) * make_state( b )
                     - make_scalar( 2.0 ) * make_state( a ) * make_state( b );

    CHECK( r1.size() == r2.size() );
    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result[i] );
    }
}

TEST_CASE( "expressions::span::expr2.3" )
{
    using raw_container_t = std::vector<double>;
    using container_t     = std::span<double>;
    using namespace ponio::expression;

    raw_container_t raw_a{ 0., 1., 2., 3., 4. }, raw_b{ 9., 8., 7., 6., 5. };

    container_t a( raw_a.data(), raw_a.size() ), b( raw_b.data(), raw_b.size() );
    raw_container_t result{ 81, 49, 25, 9, 1 };

    auto r1 = ( make_state( a ) - make_state( b ) ) * ( make_state( a ) - make_state( b ) );
    auto r2 = make_state( a ) * make_state( a ) + make_state( b ) * make_state( b ) - make_scalar( 2.0 ) * make_state( a ) * make_state( b );

    CHECK( r1.size() == r2.size() );
    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result[i] );
    }
}

// ------------------------------------------------------------------
// expr3: a^2 - b^2 == (a-b)(a+b)

TEST_CASE( "expressions::span::expr3.0" )
{
    using raw_container_t = std::vector<double>;
    using container_t     = std::span<double>;

    raw_container_t raw_a{ 0., 1., 2., 3., 4. }, raw_b{ 9., 8., 7., 6., 5. };
    container_t a( raw_a.data(), raw_a.size() ), b( raw_b.data(), raw_b.size() );
    raw_container_t result{ -81, -63, -45, -27, -9 };

    auto va = ponio::expression::make_state( a );
    auto vb = ponio::expression::make_state( b );

    raw_container_t raw_r1( a.size() ), raw_r2( a.size() );
    container_t r1( raw_r1.data(), raw_r1.size() ), r2( raw_r2.data(), raw_r2.size() );
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

TEST_CASE( "expressions::span::expr3.1" )
{
    using raw_container_t = std::vector<double>;
    using container_t     = std::span<double>;

    raw_container_t raw_a{ 0., 1., 2., 3., 4. }, raw_b{ 9., 8., 7., 6., 5. };
    container_t a( raw_a.data(), raw_a.size() ), b( raw_b.data(), raw_b.size() );
    raw_container_t result{ -81, -63, -45, -27, -9 };

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

TEST_CASE( "expressions::span::expr3.2" )
{
    using raw_container_t = std::vector<double>;
    using container_t     = std::span<double>;
    using namespace ponio::expression;

    raw_container_t raw_a{ 0., 1., 2., 3., 4. }, raw_b{ 9., 8., 7., 6., 5. };
    raw_container_t raw_r1, raw_r2;
    raw_r1.reserve( raw_a.size() );
    raw_r2.reserve( raw_a.size() );

    container_t a( raw_a.data(), raw_a.size() ), b( raw_b.data(), raw_b.size() );
    container_t r1( raw_r1.data(), raw_r1.size() ), r2( raw_r2.data(), raw_r2.size() );
    raw_container_t result{ -81, -63, -45, -27, -9 };

    make_state( r1 ) = ( make_state( a ) - make_state( b ) ) * ( make_state( a ) + make_state( b ) );
    make_state( r2 ) = make_state( a ) * make_state( a ) - make_state( b ) * make_state( b );

    CHECK( r1.size() == r2.size() );
    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result[i] );
    }
}

TEST_CASE( "expressions::span::expr3.3" )
{
    using raw_container_t = std::vector<double>;
    using container_t     = std::span<double>;
    using namespace ponio::expression;

    raw_container_t raw_a{ 0., 1., 2., 3., 4. }, raw_b{ 9., 8., 7., 6., 5. };

    container_t a( raw_a.data(), raw_a.size() ), b( raw_b.data(), raw_b.size() );
    raw_container_t result{ -81, -63, -45, -27, -9 };

    auto r1 = ( make_state( a ) - make_state( b ) ) * ( make_state( a ) + make_state( b ) );
    auto r2 = make_state( a ) * make_state( a ) - make_state( b ) * make_state( b );

    CHECK( r1.size() == r2.size() );
    for ( auto i = 0u; i < r1.size(); ++i )
    {
        CHECK( r1[i] == r2[i] );
        CHECK( r1[i] == result[i] );
    }
}
