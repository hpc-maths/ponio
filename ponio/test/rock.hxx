// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <cmath>
#include <cstddef>
#include <string_view>
#include <utility>

#include <doctest/doctest.h>

#include <ponio/observer.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

namespace rock_test
{
    // error at t = 1 on y' = cos t, y(0) = 0, whose exact solution is sin t
    template <typename Algorithm_t>
    double
    cos_quadrature_error( Algorithm_t&& algo, double dt )
    {
        auto pb = []( double t, double, double& dy )
        {
            dy = std::cos( t );
        };

        double const tf = 1.;
        double const y  = ponio::solve( pb,
            std::forward<Algorithm_t>( algo ),
            0.,
            ponio::time_span<double>{ 0., tf },
            dt,
            ponio::observer::null_observer() );

        return std::abs( y - std::sin( tf ) );
    }
} // namespace rock_test

TEST_CASE_TEMPLATE_DEFINE( "rock::non_autonomous_order", rock_tag, rock_non_autonomous_order )
{
    // y' = cos t does not depend on y, the declared spectral radius only forces many Chebyshev stages,
    // so the error comes from the times at which each stage evaluates the right hand side
    auto eig_computer = []( auto&, double, auto&, double, auto& )
    {
        return 1e3;
    };

    double const err_coarse = rock_test::cos_quadrature_error( rock_tag::make( eig_computer ), 1. / 20. );
    double const err_fine   = rock_test::cos_quadrature_error( rock_tag::make( eig_computer ), 1. / 40. );

    INFO( rock_tag::id );
    CHECK( std::log2( err_coarse / err_fine ) == doctest::Approx( static_cast<double>( rock_tag::order ) ).epsilon( 0.1 ) );
}

struct rock2_tag
{
    static constexpr std::string_view id = "ROCK2";
    static constexpr std::size_t order   = 2;

    template <typename eig_computer_t>
    static auto
    make( eig_computer_t&& eig_computer )
    {
        return ponio::runge_kutta::rock::rock2( std::forward<eig_computer_t>( eig_computer ) );
    }
};

struct rock4_tag
{
    static constexpr std::string_view id = "ROCK4";
    static constexpr std::size_t order   = 4;

    template <typename eig_computer_t>
    static auto
    make( eig_computer_t&& eig_computer )
    {
        return ponio::runge_kutta::rock::rock4( std::forward<eig_computer_t>( eig_computer ) );
    }
};

TEST_CASE_TEMPLATE_INVOKE( rock_non_autonomous_order, rock2_tag, rock4_tag );
