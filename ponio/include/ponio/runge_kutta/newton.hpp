// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private

#pragma once

#include <cmath>
#include <concepts>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <utility>

#include "../detail.hpp"
#include "../linear_algebra.hpp"
#include "../ponio_config.hpp"

namespace ponio::runge_kutta::nonlinear_solver
{
    /**
     * @brief Solve a nonlinear system with Newton's method.
     *
     * The Jacobian is evaluated at every nonlinear iteration. The supplied
     * linear solver is therefore responsible for factorizing the current
     * Jacobian and solving the associated correction system.
     *
     * @param f nonlinear residual
     * @param df Jacobian of the nonlinear residual
     * @param x0 initial guess
     * @param solver linear solver
     * @param tol nonlinear tolerance
     * @param max_iter maximum number of nonlinear iterations
     *
     * @return last Newton iterate
     */
    template <typename value_t, typename state_t, typename func_t, typename jacobian_t, typename solver_t>
    state_t
    newton( func_t&& f,
        jacobian_t&& df,
        state_t const& x0,
        solver_t&& solver,
        value_t tol          = ponio::default_config::newton_tolerance,
        std::size_t max_iter = ponio::default_config::newton_max_iterations )
    {
        state_t xk       = x0;
        value_t residual = ::ponio::detail::norm( std::forward<func_t>( f )( xk ) );
        std::size_t iter = 0;

        while ( iter < max_iter && residual > tol )
        {
            auto increment = std::forward<solver_t>( solver )( std::forward<jacobian_t>( df )( xk ), -std::forward<func_t>( f )( xk ) );

            xk       = xk + increment;
            residual = ::ponio::detail::norm( std::forward<func_t>( f )( xk ) );

            iter += 1;
        }

        return xk;
    }

    /**
     * @brief Solve a PIROCK implicit stage with a simplified Newton iteration.
     *
     * The factorization is supplied by the caller. The first corrections reuse
     * the current factorization; after @p refresh_after corrections, the
     * factorization is refreshed at the current iterate before every solve.
     * Convergence is tested on the mean absolute norm of the correction.
     *
     * @param f nonlinear residual
     * @param x0 initial guess
     * @param refresh_factorization callback that rebuilds and factorizes the stage matrix
     * @param tol nonlinear tolerance
     * @param max_iter maximum number of nonlinear iterations
     * @param refresh_after number of corrections performed before refreshing the factorization
     *
     * @return last simplified-Newton iterate
     */
    template <typename value_t, typename state_t, typename matrix_t, typename func_t, typename refresh_t>
        requires std::invocable<refresh_t&, state_t const&>
    state_t
    pirock_newton( func_t&& f,
        state_t const& x0,
        refresh_t&& refresh_factorization,
        value_t tol               = ponio::default_config::newton_tolerance,
        std::size_t max_iter      = ponio::default_config::newton_max_iterations,
        std::size_t refresh_after = 5 )
    {
        state_t xk              = x0;
        state_t residual_vector = std::forward<func_t>( f )( xk );
        state_t residual_rhs;
        state_t increment;

        value_t const n_dof = [&]() -> value_t
        {
            if constexpr ( requires {
                               std::begin( x0 );
                               std::end( x0 );
                           } )
            {
                return static_cast<value_t>( std::distance( std::begin( x0 ), std::end( x0 ) ) );
            }
            else if constexpr ( requires { x0.size(); } )
            {
                return static_cast<value_t>( x0.size() );
            }
            else
            {
                return value_t{ 1 };
            }
        }();

        value_t residual_norm = static_cast<value_t>( ::ponio::detail::mean_abs_norm( residual_vector, n_dof ) );
        if ( !std::isfinite( static_cast<double>( residual_norm ) ) )
        {
            throw std::runtime_error( "PIROCK Newton initial residual is non-finite." );
        }

        std::size_t iter = 0;
        bool necessary   = true;

        while ( necessary )
        {
            if ( iter >= refresh_after )
            {
                refresh_factorization( xk );
            }

            residual_rhs = -residual_vector;
            increment    = ::ponio::linear_algebra::linear_algebra<matrix_t>::solve_factorized( residual_rhs );

            value_t const increment_norm = static_cast<value_t>( ::ponio::detail::mean_abs_norm( increment, n_dof ) );
            if ( !std::isfinite( static_cast<double>( increment_norm ) ) )
            {
                throw std::runtime_error( "PIROCK Newton correction is non-finite." );
            }

            xk += increment;
            ++iter;

            necessary = increment_norm > tol && iter < max_iter;

            // Reevaluate the residual only when another correction is needed.
            if ( necessary )
            {
                residual_vector = std::forward<func_t>( f )( xk );
                residual_norm   = static_cast<value_t>( ::ponio::detail::mean_abs_norm( residual_vector, n_dof ) );

                if ( !std::isfinite( static_cast<double>( residual_norm ) ) )
                {
                    throw std::runtime_error( "PIROCK Newton residual became non-finite." );
                }
            }
        }

        return xk;
    }

    /**
     * @brief Solve a nonlinear system with the chord method (modified Newton).
     *
     * The Jacobian/factorization is frozen by the caller before entering this
     * routine. Each nonlinear iteration therefore reuses the same factorization
     * and only reevaluates the nonlinear residual. This is the generic
     * modified-Newton path used by DIRK methods when factorize() and
     * solve_factorized() are available.
     *
     * Unlike pirock_newton(), convergence is tested on the nonlinear residual
     * and no PIROCK-specific refresh policy is applied.
     *
     * @param f nonlinear residual
     * @param x0 initial guess
     * @param tol nonlinear tolerance
     * @param max_iter maximum number of nonlinear iterations
     *
     * @return last chord-method iterate
     */
    template <typename value_t, typename state_t, typename matrix_t, typename func_t>
    state_t
    chord_method( func_t&& f,
        state_t const& x0,
        value_t tol          = ponio::default_config::newton_tolerance,
        std::size_t max_iter = ponio::default_config::newton_max_iterations )
    {
        state_t xk              = x0;
        state_t residual_vector = std::forward<func_t>( f )( xk );
        state_t residual_rhs;
        state_t increment;

        value_t residual = ::ponio::detail::norm( residual_vector );
        std::size_t iter = 0;

        while ( iter < max_iter && residual > tol )
        {
            residual_rhs = -residual_vector;
            increment    = ::ponio::linear_algebra::linear_algebra<matrix_t>::solve_factorized( residual_rhs );

            xk += increment;

            residual_vector = std::forward<func_t>( f )( xk );
            residual        = ::ponio::detail::norm( residual_vector );

            iter += 1;
        }

        return xk;
    }

} // namespace ponio::runge_kutta::nonlinear_solver
