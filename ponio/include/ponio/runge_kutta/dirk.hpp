// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private, include "../runge_kutta.hpp"

#pragma once

#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional> // NOLINT(misc-include-cleaner)
#include <numeric>
#include <ranges>
#include <string_view>
#include <type_traits>

#include "../butcher_tableau.hpp"
#include "../detail.hpp"
#include "../iteration_info.hpp"
#include "../linear_algebra.hpp"
#include "../ponio_config.hpp"
#include "../stage.hpp" // NOLINT(misc-include-cleaner)

namespace ponio::runge_kutta::diagonal_implicit_runge_kutta
{
    template <typename value_t, typename state_t, typename func_t, typename jacobian_t, typename solver_t>
    state_t
    newton( func_t&& f, jacobian_t&& df, state_t const& x0, solver_t&& solver, value_t tol = 1e-10, std::size_t max_iter = 50 )
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

    template <typename tableau_t, typename lin_alg_t = void>
    struct diagonal_implicit_rk_butcher
    {
        tableau_t butcher;
        static constexpr std::size_t N_stages     = tableau_t::N_stages;
        static constexpr bool is_embedded         = butcher::is_embedded_tableau<tableau_t>;
        static constexpr std::size_t order        = tableau_t::order;
        static constexpr std::string_view id      = tableau_t::id;
        static constexpr bool void_linear_algebra = std::is_void_v<lin_alg_t>;
        using linear_algebra_t                    = typename std::conditional_t<void_linear_algebra,
            bool, // just a small valid type
            lin_alg_t>;

        using value_t = typename tableau_t::value_t;

        template <typename>
        diagonal_implicit_rk_butcher( double tol_, std::size_t max_iter_ )
            : butcher()
            , tol( tol_ )
            , max_iter( max_iter_ )
            , linalg( false )
        {
        }

        template <typename... args_t>
        diagonal_implicit_rk_butcher( double tol_, std::size_t max_iter_, args_t... args )
            : butcher()
            , tol( tol_ )
            , max_iter( max_iter_ )
            , linalg( args... )
        {
        }

        template <typename... args_t>
            requires detail::has_newton_method<lin_alg_t>
        diagonal_implicit_rk_butcher( args_t... args )
            : butcher()
            , tol( ponio::default_config::newton_tolerance )
            , max_iter( ponio::default_config::newton_max_iterations )
            , linalg( args... )
        {
        }

        template <typename problem_t, typename state_t, typename array_kj_t, std::size_t I>
            requires detail::problem_operator<problem_t, value_t>
        void
        stage( Stage<I>, problem_t& pb, value_t tn, state_t& un, array_kj_t const& Kj, value_t dt, state_t& ui, state_t& ki )
        {
            if constexpr ( I == 0 )
            {
                _info.reset_eval();
            }

            auto op_i = ::ponio::linear_algebra::operator_algebra<state_t>::identity( un )
                      - dt * butcher.A[I][I] * pb.f_t( tn + butcher.c[I] * dt );
            auto& rhs = ki;
            rhs       = un;
            detail::tpl_inner_product<I>( butcher.A[I], Kj, un, dt, rhs );

            std::size_t n_eval = 0;
            ::ponio::linear_algebra::operator_algebra<state_t>::solve( op_i, ui, rhs, n_eval );

            _info.number_of_eval += n_eval + 1;

            pb.f( tn + butcher.c[I] * dt, ui, ki );
        }

        template <typename problem_t, typename state_t, typename array_kj_t, std::size_t I>
            requires detail::problem_jacobian<problem_t, value_t, state_t>
        void
        stage( Stage<I>, problem_t& pb, value_t tn, state_t& un, array_kj_t const& Kj, value_t dt, state_t& ui, state_t& ki )
        {
            if constexpr ( I == 0 )
            {
                _info.reset_eval();
            }

            using matrix_t = decltype( pb.df( tn, un ) );

            auto identity = [&]( state_t const& u )
            {
                if constexpr ( detail::has_identity_method<lin_alg_t> )
                {
                    return linalg.identity( u );
                }
                else
                {
                    return ::ponio::linear_algebra::linear_algebra<matrix_t>::identity( u );
                }
            }( un );
            // lambda function `g` that equals to :
            // $$
            //   g : k \mapsto k - u^n - \Delta t f( tn+c_i\Delta t, u^n + \Delta t \sum_{j=0}^{i-1} a_{ij}k_j + \Delta t a_{ii}k )
            // $$
            // and compute the lambda function `dg` as $dg = \partial_k g(t,k)$ :
            // $$
            //   dg : k \mapsto I - a_{ii}\Delta t \partial_k f( tn+c_i\Delta t, u^n + \Delta t \sum_{j=0}^{i-1} a_{ij}k_j + \Delta t
            //   a_{ii}k )
            // $$
            auto g = [&]( state_t const& k ) -> state_t
            {
                _info.number_of_eval += 1;
                detail::tpl_inner_product<I>( butcher.A[I], Kj, un, dt, ui );
                ui = ui + dt * butcher.A[I][I] * k;
                pb.f( tn + butcher.c[I] * dt, ui, ki );
                return k - ki;
            };
            auto dg = [&]( state_t const& k ) -> matrix_t
            {
                detail::tpl_inner_product<I>( butcher.A[I], Kj, un, dt, ui );
                ui = ui + dt * butcher.A[I][I] * k;
                return identity - butcher.A[I][I] * dt * pb.df( tn + butcher.c[I] * dt, ui );
            };

            // call newton method
            if constexpr ( detail::has_newton_method<lin_alg_t, decltype( g ), decltype( dg )>
                           || detail::has_newton_method<lin_alg_t, decltype( g ), decltype( dg ), state_t> )
            {
                linalg.newton( g, dg, un );
            }
            else
            {
                auto solver = [&]()
                {
                    if constexpr ( detail::has_solver_method<lin_alg_t, matrix_t, state_t> )
                    {
                        using namespace std::placeholders;
                        return std::bind( &lin_alg_t::solver, linalg, _1, _2 );
                    }
                    else
                    {
                        return &::ponio::linear_algebra::linear_algebra<matrix_t>::solver;
                    }
                }();
                newton<value_t>( g, dg, un, solver, tol, max_iter );
            }
        }

        template <typename problem_t, typename state_t, typename array_kj_t>
        void
        stage( Stage<N_stages>, problem_t&, value_t, state_t& un, array_kj_t const& Kj, value_t dt, state_t&, state_t& ki )
        {
            // last stage is always explicit and just equals to:
            // $$
            //   u^{n+1} = u^n + \Delta t \sum_{i} b_i k_i
            // $$
            detail::tpl_inner_product<N_stages>( butcher.b, Kj, un, dt, ki );
        }

        template <typename problem_t, typename state_t, typename array_kj_t, typename tab_t = tableau_t>
            requires std::same_as<tab_t, tableau_t> && is_embedded
        void
        stage( Stage<N_stages + 1>, problem_t&, value_t, state_t& un, array_kj_t const& Kj, value_t dt, state_t&, state_t& ki )
        {
            detail::tpl_inner_product<N_stages>( butcher.b2, Kj, un, dt, ki );
        }

        /**
         * @brief gets `iteration_info` object
         */
        auto&
        info()
        {
            return _info;
        }

        /**
         * @brief gets `iteration_info` object (constant version)
         */
        auto const&
        info() const
        {
            return _info;
        }

        double tol           = ponio::default_config::newton_tolerance;      // tolerance of Newton method
        std::size_t max_iter = ponio::default_config::newton_max_iterations; // max iterations of Newton method

        linear_algebra_t linalg;
        iteration_info<tableau_t> _info;
    };

    // ---- *helper* ----

    /**
     * factory of DIRK method
     *
     * @tparam tableau_t type of Butcher tableau
     * @tparam lin_alg_t type of linear algebra
     * @tparam args_t    optional and variadic types of arguments
     * @param tol      tolenrence of Newton's method
     * @param max_iter maximum iteration of Newton's method
     * @param args     optional arguments to initialize an instance of type lin_alg_t
     */
    template <typename tableau_t, typename lin_alg_t, typename... args_t>
    auto
    make_dirk( double tol    = ponio::default_config::newton_tolerance,
        std::size_t max_iter = ponio::default_config::newton_max_iterations,
        args_t... args )
    {
        return diagonal_implicit_rk_butcher<tableau_t, lin_alg_t>( tol, max_iter, args... );
    }

} // namespace ponio::runge_kutta::diagonal_implicit_runge_kutta
