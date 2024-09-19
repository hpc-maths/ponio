// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <cmath>

#include <doctest/doctest.h>

#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/splitting.hpp>

TEST_CASE( "number_of_eval::explicit_runge_kutta" )
{
    std::size_t manual_counter = 0;

    double const k            = 50;
    auto curtiss_hirschfelder = [&, k]( double t, double y )
    {
        ++manual_counter;
        return k * ( y + std::cos( t ) );
    };

    double y_0 = 2.0;

    ponio::time_span<double> const t_span = { 0., 2. };
    double const dt                       = 0.05;

    auto sol_range = ponio::make_solver_range( curtiss_hirschfelder, ponio::runge_kutta::rk_33(), y_0, t_span, dt );
    auto it_sol    = sol_range.begin();

    std::size_t cumulative_counter = 0;
    while ( it_sol->time < t_span.back() )
    {
        ++it_sol;
        cumulative_counter += it_sol.info().number_of_eval;
    }

    CHECK( cumulative_counter == manual_counter );
}

TEST_CASE( "number_of_eval::lawson_runge_kutta" )
{
    std::size_t manual_counter = 0;

    double const k = 50;
    auto n         = [&, k]( double t, double )
    {
        ++manual_counter;
        return k * std::cos( t );
    };

    double y_0 = 2.0;

    ponio::time_span<double> const t_span = { 0., 2. };
    double const dt                       = 0.05;

    auto curtiss_hirschfelder = ponio::make_lawson_problem( k, n );

    auto exp = []( double x )
    {
        return std::exp( x );
    };

    auto sol_range = ponio::make_solver_range( curtiss_hirschfelder, ponio::runge_kutta::lrk_33( exp ), y_0, t_span, dt );
    auto it_sol    = sol_range.begin();

    std::size_t cumulative_counter = 0;
    while ( it_sol->time < t_span.back() )
    {
        ++it_sol;
        cumulative_counter += it_sol.info().number_of_eval;
    }

    CHECK( cumulative_counter == manual_counter );
}

TEST_CASE( "number_of_eval::exponential_runge_kutta" )
{
    std::size_t manual_counter = 0;

    double const k = 50;
    auto n         = [&, k]( double t, double )
    {
        ++manual_counter;
        return k * std::cos( t );
    };

    double y_0 = 2.0;

    ponio::time_span<double> const t_span = { 0., 2. };
    double const dt                       = 0.05;

    auto curtiss_hirschfelder = ponio::make_lawson_problem( k, n );
    auto sol_range            = ponio::make_solver_range( curtiss_hirschfelder, ponio::runge_kutta::exprk22(), y_0, t_span, dt );
    auto it_sol               = sol_range.begin();

    std::size_t cumulative_counter = 0;
    while ( it_sol->time < t_span.back() )
    {
        ++it_sol;
        cumulative_counter += it_sol.info().number_of_eval;
    }

    CHECK( cumulative_counter == manual_counter );
}

// TEST_CASE("number_of_eval::rock")
// {

// }

// TEST_CASE("number_of_eval::pirock")
// {

// }

TEST_CASE( "number_of_eval::splitting_lie" )
{
    std::size_t manual_counter_1 = 0;
    std::size_t manual_counter_2 = 0;

    double const k = 50;
    auto f1        = [&]( double, double y )
    {
        ++manual_counter_1;
        return k * y;
    };
    auto f2 = [&, k]( double t, double )
    {
        ++manual_counter_2;
        return k * std::cos( t );
    };

    double y_0 = 2.0;

    ponio::time_span<double> const t_span = { 0., 2. };
    double const dt                       = 0.05;

    auto curtiss_hirschfelder = ponio::make_problem( f1, f2 );
    auto lie                  = ponio::splitting::make_lie_tuple( std::make_pair( ponio::runge_kutta::rk_33(), 0.125 * dt ),
        std::make_pair( ponio::runge_kutta::rk_33(), 0.25 * dt ) );
    auto sol_range            = ponio::make_solver_range( curtiss_hirschfelder, lie, y_0, t_span, dt );
    auto it_sol               = sol_range.begin();

    std::size_t cumulative_counter_1 = 0;
    std::size_t cumulative_counter_2 = 0;
    while ( it_sol->time < t_span.back() )
    {
        ++it_sol;
        cumulative_counter_1 += it_sol.info().number_of_eval[0];
        cumulative_counter_2 += it_sol.info().number_of_eval[1];
    }

    CHECK( cumulative_counter_1 == manual_counter_1 );
    CHECK( cumulative_counter_2 == manual_counter_2 );
}

// TEST_CASE("number_of_eval::splitting_strang")
// {

// }

// TEST_CASE("number_of_eval::splitting_adaptive_strang")
// {

// }
