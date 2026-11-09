// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private

#pragma once

#include <cstddef>
#include <utility>

#include "../detail.hpp"
#include "../ponio_config.hpp"

namespace ponio::runge_kutta::diagonal_implicit_runge_kutta
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

} // namespace ponio::runge_kutta::diagonal_implicit_runge_kutta
