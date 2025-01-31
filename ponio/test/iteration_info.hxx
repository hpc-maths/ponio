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

/**
 * In this test case we solve the Curtiss and Hirschfelder problem:
 *
 * \f$$
 *  \begin{aligned}
 *    \dot{y} = k(\cos(t) - y) \\
 *    y(0) = 2
 *  \end{aligned}
 * \f$$
 *
 * We solve this equation with a represent of each class of numerical integrator, simply to test if number of evaluation of function is
 * correctly computing.
 *
 */

/**
 * ----------------------------------------------------------------------------
 *
 * # Explicit Runge-Kutta method
 *
 * This class of methods are writing to solve a problem with form as
 *
 * \f$$
 *  \dot{y} = f(t, y)
 * \f$$
 *
 * with simply write the function `curtiss_hirschfelder` for this.
 */
TEST_CASE( "number_of_eval::explicit_runge_kutta" )
{
    std::size_t manual_counter = 0;

    double const k            = 50;
    auto curtiss_hirschfelder = [&, k]( double t, double y )
    {
        ++manual_counter;
        return k * ( std::cos( t ) - y );
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

/**
 * ----------------------------------------------------------------------------
 *
 * # Lawson and exponential Runge-Kutta methods
 *
 * This classes of methods are writing to solve a problem with form as
 *
 * \f$$
 *  \dot{y} = Ly + N(t, y)
 * \f$$
 *
 * with define \f$L=-k\f$ and \f$N:t, y\mapsto k\cos(t)\f$, we build an object `curtiss_hirschfelder` to keep link between two mathematic
 * objects.
 */
TEST_CASE( "number_of_eval::lawson_runge_kutta" )
{
    std::size_t manual_counter = 0;

    double const k = 50;
    double const l = -k;
    auto n         = [&, k]( double t, double )
    {
        ++manual_counter;
        return k * std::cos( t );
    };

    double y_0 = 2.0;

    ponio::time_span<double> const t_span = { 0., 2. };
    double const dt                       = 0.05;

    auto curtiss_hirschfelder = ponio::make_lawson_problem( l, n );

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
    double const l = -k;
    auto n         = [&, k]( double t, double )
    {
        ++manual_counter;
        return k * std::cos( t );
    };

    double y_0 = 2.0;

    ponio::time_span<double> const t_span = { 0., 2. };
    double const dt                       = 0.05;

    auto curtiss_hirschfelder = ponio::make_lawson_problem( l, n );
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

/**
 * ----------------------------------------------------------------------------
 *
 * # Diagonal Implicit Runge-Kutta method
 *
 * This class of methods are writing to solve a problem with form as
 *
 * \f$$
 *  \dot{y} = f(t, y)
 * \f$$
 *
 * we also need the Jacobian function of \f$f\$f (store as `df`) to implicit the problem, we build an object `curtiss_hirschfelder` to keep
 * link between two mathematic objects.
 */
TEST_CASE( "number_of_eval::diagonal_implicit_runge_kutta" )
{
    std::size_t manual_counter = 0;

    double const k = 50;
    auto f         = [&, k]( double t, double y )
    {
        ++manual_counter;
        return k * ( std::cos( t ) - y );
    };

    auto df = [=]( double, double )
    {
        return -k;
    };

    double y_0 = 2.0;

    ponio::time_span<double> const t_span = { 0., 2. };
    double const dt                       = 0.05;

    auto curtiss_hirschfelder = ponio::make_implicit_problem( f, df );
    auto sol_range            = ponio::make_solver_range( curtiss_hirschfelder, ponio::runge_kutta::dirk34(), y_0, t_span, dt );
    auto it_sol               = sol_range.begin();

    std::size_t cumulative_counter = 0;
    while ( it_sol->time < t_span.back() )
    {
        ++it_sol;
        cumulative_counter += it_sol.info().number_of_eval;
    }

    CHECK( cumulative_counter == manual_counter );
}

/**
 * ----------------------------------------------------------------------------
 *
 * # ROCK method
 *
 * This class of methods are writing to solve a problem with form as
 *
 * \f$$
 *  \dot{y} = f(t, y)
 * \f$$
 *
 * with \f$f\f$ an operator with large negative eigenvalues, this kind of methods has an adaptive number of stages, and needs more
 * evaluation of function \f$f\f$ to estimate spectral radius.
 */
TEST_CASE( "number_of_eval::rock2" )
{
    std::size_t manual_counter = 0;

    double const k            = 50;
    auto curtiss_hirschfelder = [&, k]( double t, double y )
    {
        ++manual_counter;
        return k * ( std::cos( t ) - y );
    };

    double y_0 = 2.0;

    ponio::time_span<double> const t_span = { 0., 2. };
    double const dt                       = 0.05;

    auto sol_range = ponio::make_solver_range( curtiss_hirschfelder, ponio::runge_kutta::rock::rock2(), y_0, t_span, dt );
    auto it_sol    = sol_range.begin();

    std::size_t cumulative_counter = 0;
    while ( it_sol->time < t_span.back() )
    {
        ++it_sol;
        cumulative_counter += it_sol.info().number_of_eval;
    }

    CHECK( cumulative_counter == manual_counter );
}

TEST_CASE( "number_of_eval::rock4" )
{
    std::size_t manual_counter = 0;

    double const k            = 50;
    auto curtiss_hirschfelder = [&, k]( double t, double y )
    {
        ++manual_counter;
        return k * ( std::cos( t ) - y );
    };

    double y_0 = 2.0;

    ponio::time_span<double> const t_span = { 0., 2. };
    double const dt                       = 0.05;

    auto sol_range = ponio::make_solver_range( curtiss_hirschfelder, ponio::runge_kutta::rock::rock4(), y_0, t_span, dt );
    auto it_sol    = sol_range.begin();

    std::size_t cumulative_counter = 0;
    while ( it_sol->time < t_span.back() )
    {
        ++it_sol;
        cumulative_counter += it_sol.info().number_of_eval;
    }

    CHECK( cumulative_counter == manual_counter );
}

/**
 * ----------------------------------------------------------------------------
 *
 * # PIROCK method
 *
 * This class of methods are writing to solve a problem with form as
 *
 * \f$$
 *  \dot{y} = f(t, y) + g(t,y)
 * \f$$
 *
 * with \f$f\f$ an operator with large negative eigenvalues, solved with ROCK2 method, and \f$g\f$ an operator which we would like solve
 * with an implicit method (so we need the Jacobian function), we store respectively in `f_ex`, `f_im` (with its Jacobian `df_im`). We
 * choose \f$f:t,y\mapsto k\cos(t)\f$, and \f$g:t,y\mapsto -k y\f$. We build an object `curtiss_hirschfelder` to keep link between three
 * mathematic objects.
 */
TEST_CASE( "number_of_eval::pirock" )
{
    std::size_t manual_counter_im = 0;
    std::size_t manual_counter_ex = 0;

    double const k = 50;
    auto f_im      = [&]( double, double y )
    {
        ++manual_counter_im;
        return -k * y;
    };
    auto df_im = [&]( double, double )
    {
        return -k;
    };
    auto f_ex = [&, k]( double t, double )
    {
        ++manual_counter_ex;
        return k * std::cos( t );
    };

    double y_0 = 2.0;

    ponio::time_span<double> const t_span = { 0., 2. };
    double const dt                       = 0.05;

    auto curtiss_hirschfelder = ponio::make_imex_jacobian_problem( f_ex, f_im, df_im );
    auto sol_range            = ponio::make_solver_range( curtiss_hirschfelder, ponio::runge_kutta::pirock::pirock<1>(), y_0, t_span, dt );
    auto it_sol               = sol_range.begin();

    std::size_t cumulative_counter_ex = 0;
    std::size_t cumulative_counter_im = 0;
    while ( it_sol->time < t_span.back() )
    {
        ++it_sol;
        cumulative_counter_ex += std::get<0>( it_sol.info().number_of_eval );
        cumulative_counter_im += std::get<1>( it_sol.info().number_of_eval );
    }

    CHECK( cumulative_counter_im == manual_counter_im );
    CHECK( cumulative_counter_ex == manual_counter_ex );
}

/**
 * ----------------------------------------------------------------------------
 *
 * # Splitting operator method
 *
 * This class of methods are writing to solve a problem with form as
 *
 * \f$$
 *  \dot{y} = \sum_i f_i(t, y)
 * \f$$
 *
 * We split the problem into some subproblem
 *
 * \f$$
 *  \dot{y} = f_i(t, y)
 * \f$$
 *
 * and solve consecutively each subproblem.
 */
TEST_CASE( "number_of_eval::splitting_lie" )
{
    std::size_t manual_counter_1 = 0;
    std::size_t manual_counter_2 = 0;

    double const k = 50;
    auto f1        = [&]( double, double y )
    {
        ++manual_counter_1;
        return -k * y;
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
        cumulative_counter_1 += std::get<0>( it_sol.info().number_of_eval );
        cumulative_counter_2 += std::get<1>( it_sol.info().number_of_eval );
    }

    CHECK( cumulative_counter_1 == manual_counter_1 );
    CHECK( cumulative_counter_2 == manual_counter_2 );
}

TEST_CASE( "number_of_eval::splitting_strang" )
{
    std::size_t manual_counter_1 = 0;
    std::size_t manual_counter_2 = 0;

    double const k = 50;
    auto f1        = [&]( double, double y )
    {
        ++manual_counter_1;
        return -k * y;
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
    auto strang               = ponio::splitting::make_strang_tuple( std::make_pair( ponio::runge_kutta::rk_33(), 0.125 * dt ),
        std::make_pair( ponio::runge_kutta::rk_33(), 0.25 * dt ) );
    auto sol_range            = ponio::make_solver_range( curtiss_hirschfelder, strang, y_0, t_span, dt );
    auto it_sol               = sol_range.begin();

    std::size_t cumulative_counter_1 = 0;
    std::size_t cumulative_counter_2 = 0;
    while ( it_sol->time < t_span.back() )
    {
        ++it_sol;
        cumulative_counter_1 += std::get<0>( it_sol.info().number_of_eval );
        cumulative_counter_2 += std::get<1>( it_sol.info().number_of_eval );
    }

    CHECK( cumulative_counter_1 == manual_counter_1 );
    CHECK( cumulative_counter_2 == manual_counter_2 );
}

TEST_CASE( "number_of_eval::splitting_adaptive_strang" )
{
    std::size_t manual_counter_1 = 0;
    std::size_t manual_counter_2 = 0;

    double const k = 50;
    auto f1        = [&]( double, double y )
    {
        ++manual_counter_1;
        return -k * y;
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

    double delta         = 5e-3;
    double tol           = 1e-3;
    auto adaptive_strang = ponio::splitting::make_adaptive_strang_tuple( delta,
        tol,
        std::make_pair( ponio::runge_kutta::rk_33(), 0.125 * dt ),
        std::make_pair( ponio::runge_kutta::rk_33(), 0.25 * dt ) );
    auto sol_range       = ponio::make_solver_range( curtiss_hirschfelder, adaptive_strang, y_0, t_span, dt );
    auto it_sol          = sol_range.begin();

    std::size_t cumulative_counter_1 = 0;
    std::size_t cumulative_counter_2 = 0;
    while ( it_sol->time < t_span.back() )
    {
        ++it_sol;
        cumulative_counter_1 += std::get<0>( it_sol.info().number_of_eval );
        cumulative_counter_2 += std::get<1>( it_sol.info().number_of_eval );
    }

    CHECK( cumulative_counter_1 == manual_counter_1 );
    CHECK( cumulative_counter_2 == manual_counter_2 );
}
