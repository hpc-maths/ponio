// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <vector>

#include <doctest/doctest.h>

#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

namespace solver_test
{
    struct run_record
    {
        std::vector<double> times;
        std::size_t n_eval = 0;
    };

    // times seen by the observer of ponio::solve with explicit Euler method on y' = -y
    inline run_record
    solve_times( ponio::time_span<double> const& t_span, double dt )
    {
        run_record record;

        auto pb = [&]( double, double y, double& dy )
        {
            ++record.n_eval;
            dy = -y;
        };
        auto obs = [&]( double t, double, double )
        {
            record.times.push_back( t );
        };

        ponio::solve( pb, ponio::runge_kutta::euler(), 1., t_span, dt, obs );

        return record;
    }

    // times visited by ponio::make_solver_range with explicit Euler method on y' = -y
    inline run_record
    range_times( ponio::time_span<double> const& t_span, double dt )
    {
        run_record record;

        auto pb = [&]( double, double y, double& dy )
        {
            ++record.n_eval;
            dy = -y;
        };

        for ( auto const& ds : ponio::make_solver_range( pb, ponio::runge_kutta::euler(), 1., t_span, dt ) )
        {
            record.times.push_back( ds.time );
        }

        return record;
    }

    inline void
    check_times( run_record const& record, std::vector<double> const& expected )
    {
        REQUIRE( record.times.size() == expected.size() );
        for ( std::size_t i = 0; i < expected.size(); ++i )
        {
            INFO( "time index ", i );
            CHECK( record.times[i] == doctest::Approx( expected[i] ) );
        }
        // the last time is the end of the time span, not an approximation of it
        CHECK( record.times.back() == expected.back() );
        // explicit Euler method evaluates the right hand side once per step
        CHECK( record.n_eval == expected.size() - 1 );
    }
} // namespace solver_test

TEST_CASE( "solver::time_span" )
{
    SUBCASE( "time step divides the time span" )
    {
        std::vector<double> const expected = { 0., 0.25, 0.5, 0.75, 1. };
        solver_test::check_times( solver_test::solve_times( { 0., 1. }, 0.25 ), expected );
        solver_test::check_times( solver_test::range_times( { 0., 1. }, 0.25 ), expected );
    }

    SUBCASE( "last step is shortened to stop on the final time" )
    {
        std::vector<double> const expected = { 0., 0.3, 0.6, 0.9, 1. };
        solver_test::check_times( solver_test::solve_times( { 0., 1. }, 0.3 ), expected );
        solver_test::check_times( solver_test::range_times( { 0., 1. }, 0.3 ), expected );
    }

    SUBCASE( "intermediate time reached by a regular step" )
    {
        std::vector<double> const expected = { 0., 0.25, 0.5, 0.75, 1. };
        solver_test::check_times( solver_test::solve_times( { 0., 0.5, 1. }, 0.25 ), expected );
        solver_test::check_times( solver_test::range_times( { 0., 0.5, 1. }, 0.25 ), expected );
    }

    SUBCASE( "step shortened to stop on an intermediate time" )
    {
        std::vector<double> const expected = { 0., 0.3, 0.5, 0.8, 1. };
        solver_test::check_times( solver_test::solve_times( { 0., 0.5, 1. }, 0.3 ), expected );
        solver_test::check_times( solver_test::range_times( { 0., 0.5, 1. }, 0.3 ), expected );
    }

    SUBCASE( "time step not exactly representable" )
    {
        // ten steps of 0.1 add up to 0.9999999999999999, the solver must not take an eleventh step of 1e-16
        std::vector<double> const expected = { 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1. };
        solver_test::check_times( solver_test::solve_times( { 0., 1. }, 0.1 ), expected );
        solver_test::check_times( solver_test::range_times( { 0., 1. }, 0.1 ), expected );
    }
}

TEST_CASE( "solver::adaptive_time_step" )
{
    auto pb = []( double, double y, double& dy )
    {
        dy = -y;
    };

    std::vector<double> times;
    std::vector<double> time_steps;
    auto obs = [&]( double t, double, double dt )
    {
        times.push_back( t );
        time_steps.push_back( dt );
    };

    double const y = ponio::solve( pb, ponio::runge_kutta::rk54_6m().abs_tol( 1e-8 ).rel_tol( 1e-8 ), 1., { 0., 0.5, 1. }, 0.1, obs );

    CHECK( y == doctest::Approx( std::exp( -1. ) ).epsilon( 1e-6 ) );
    CHECK( times.back() == 1. );

    // the intermediate time is a time of the solution
    auto const it_middle = std::find( times.begin(), times.end(), 0.5 );
    REQUIRE( it_middle != times.end() );

    // the step after the intermediate time uses the time step chosen by the method before the shortened step,
    // not the initial time step 0.1
    auto const i_middle = static_cast<std::size_t>( std::distance( times.begin(), it_middle ) );
    REQUIRE( i_middle + 1 < times.size() );
    CHECK( times[i_middle + 1] - 0.5 == doctest::Approx( time_steps[i_middle - 1] ) );
}
