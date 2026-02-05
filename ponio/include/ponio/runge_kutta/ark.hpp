// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private

#pragma once

#include <array>
#include <concepts>
#include <cstddef>
#include <string_view> // NOLINT(misc-include-cleaner)

#include "../butcher_tableau.hpp"
#include "../detail.hpp"
#include "../iteration_info.hpp"
#include "../ponio_config.hpp"
#include "../stage.hpp" // NOLINT(misc-include-cleaner)

namespace ponio::runge_kutta::additive_runge_kutta
{

    template <typename tableau_pair_t, typename lin_alg_t = void>
    struct additive_runge_kutta
    {
        using tableau_im_t = typename tableau_pair_t::tableau_im_t;
        using tableau_ex_t = typename tableau_pair_t::tableau_ex_t;
        using value_t      = typename tableau_pair_t::value_t;

        static constexpr std::size_t N_stages     = tableau_pair_t::N_stages;
        static constexpr bool is_embedded         = tableau_pair_t::is_embedded;
        static constexpr bool is_imex_method      = tableau_pair_t::is_imex_method;
        static constexpr std::size_t order        = tableau_pair_t::order;
        static constexpr std::size_t N_operators  = tableau_pair_t::N_operators;
        static constexpr std::string_view id      = detail::join_id_v<tableau_im_t, detail::separator<>, tableau_ex_t>;
        static constexpr bool void_linear_algebra = std::is_void_v<lin_alg_t>;
        using linear_algebra_t                    = typename std::conditional_t<void_linear_algebra,
            bool, // just a small valid type
            lin_alg_t>;

        tableau_im_t butcher_im;
        tableau_ex_t butcher_ex;

        additive_runge_kutta()
            : butcher_im()
            , butcher_ex()
            , _info()
        {
            _info.number_of_eval[0] = N_stages; // explicit evaluation
        }

        template <typename problem_t, typename state_t, typename array_kj_t, std::size_t I>
            requires detail::problem_operator<problem_t, value_t>
        void
        stage( Stage<I>,
            problem_t& pb,
            value_t tn,
            state_t& un,
            array_kj_t const& K_ex_j,
            array_kj_t const& K_im_j,
            value_t dt,
            state_t& ui,
            state_t& u_tmp,
            state_t& k_ex_i,
            state_t& k_im_i )
        {
            if constexpr ( I == 0 )
            {
                _info.reset_eval();
            }

            // u_tmp = un + dt*sum(butcher_ex.A[I]*Kexj) + dt*sum(butcher_im.A[I]*Kimj)
            // ui = u_tmp + dt*butcher_im.A[I+1]*f(ui) <- to solve
            detail::tpl_inner_product<I>( butcher_ex.A[I], K_ex_j, un, dt, ui );
            detail::tpl_inner_product<I>( butcher_im.A[I], K_im_j, ui, dt, u_tmp );

            // solve ui - dt*butcher_im.A[I+1]*f(ui) = u_tmp
            auto op_i = ::ponio::linear_algebra::operator_algebra<state_t>::identity( un )
                      - dt * butcher_im.A[I][I] * pb.f_t( tn + butcher_im.c[I] * dt );
            auto& rhs = u_tmp;

            std::size_t n_eval = 0;
            ::ponio::linear_algebra::operator_algebra<state_t>::solve( op_i, ui, rhs, n_eval );

            _info.number_of_eval += n_eval + 1;

            // call explicit and implicit function on stage ui
            pb.explicit_part( tn + butcher_ex.c[I] * dt, ui, k_ex_i );
            pb.implicit_part( tn + butcher_ex.c[I] * dt, ui, k_im_i );
        }

        template <typename problem_t, typename state_t, typename array_kj_t, std::size_t I>
            requires detail::problem_jacobian<problem_t, value_t, state_t>
        void
        stage( Stage<I>,
            problem_t& pb,
            value_t tn,
            state_t& un,
            array_kj_t const& K_ex_j,
            array_kj_t const& K_im_j,
            value_t dt,
            state_t& ui,
            state_t& u_tmp,
            state_t& k_ex_i,
            state_t& k_im_i )
        {
            if constexpr ( I == 0 )
            {
                _info.reset_eval();
            }

            // u_tmp = un + dt*sum(butcher_ex.A[I]*Kexj) + dt*sum(butcher_im.A[I]*Kimj)
            // ui = u_tmp + dt*butcher_im.A[I+1]*f(ui) <- to solve
            detail::tpl_inner_product<I>( butcher_ex.A[I], K_ex_j, un, dt, ui );
            detail::tpl_inner_product<I>( butcher_im.A[I], K_im_j, ui, dt, u_tmp );

            // solve ui - dt*butcher_im.A[I+1]*f(ui) = u_tmp
            using matrix_t = decltype( pb.implicit_part.df( tn, un ) );

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
            // lambda function `F` that equals to :
            // ..
            //      F(u) = u - dt * ãᵢᵢ * g(tⁿ + cᵢ*dt, u) - u_tmp
            // $$
            // and compute the lambda function `dF` which is the Jacobian function of `F`:
            // ..
            //      dF(u) = I - dt * ãᵢᵢ * dg(tⁿ + cᵢ*dt, u)
            auto F = [&]( state_t const& u ) -> state_t
            {
                double const ti = tn + dt * butcher_im.c[I];
                _info.number_of_eval += 1;
                pb.implicit_part( ti, u, ui );
                return u - dt * butcher_im.A[I][I] * ui - u_tmp;
            };
            auto dF = [&]( state_t const& u ) -> matrix_t
            {
                double const ti = tn + dt * butcher_im.c[I];
                identity - dt* butcher_im.A[I][I] * pb.implicit_part.df( ti, u );
            };

            // call newton method
            if constexpr ( detail::has_newton_method<lin_alg_t, decltype( g ), decltype( dg )>
                           || detail::has_newton_method<lin_alg_t, decltype( g ), decltype( dg ), state_t> )
            {
                ui = linalg.newton( F, dF, un );
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
                ui = newton<value_t>( F, dF, un, solver, tol, max_iter );
            }

            // call explicit and implicit function on stage ui
            pb.explicit_part( tn + butcher_ex.c[I] * dt, ui, k_ex_i );
            pb.implicit_part( tn + butcher_im.c[I] * dt, ui, k_im_i );
        }

        template <typename problem_t, typename state_t, typename array_kj_t>
        void
        stage( Stage<N_stages>,
            problem_t& pb,
            value_t tn,
            state_t& un,
            array_kj_t const& K_ex_j,
            array_kj_t const& K_im_j,
            value_t dt,
            state_t& ui,
            state_t&,
            state_t& unp1,
            state_t& )
        {
            // ui = un + dt*sum( butcher_ex.b[k] * K_ex_j[k] )
            detail::tpl_inner_product<N_stages>( butcher_ex.b, K_ex_j, un, dt, ui );
            // unp1 = ui + dt*sum( butcher_ex.b[k] * K_im_j[k] )
            detail::tpl_inner_product<N_stages>( butcher_im.b, K_im_j, ui, dt, unp1 );
        }

        template <typename problem_t, typename state_t, typename array_kj_t, typename tab_t = tableau_t>
            requires std::same_as<tab_t, tableau_t> && is_embedded
        void
        stage( Stage<N_stages + 1>,
            problem_t& pb,
            value_t tn,
            state_t& un,
            array_kj_t const& K_ex_j,
            array_kj_t const& K_im_j,
            value_t dt,
            state_t& ui,
            state_t&,
            state_t& unp1_bis,
            state_t& )
        {
            // ui = un + dt*sum( butcher_ex.b2[k] * K_ex_j[k] )
            detail::tpl_inner_product<N_stages>( butcher_ex.b2, K_ex_j, un, dt, ui );
            // unp1_bis = ui + dt*sum( butcher_ex.b2[k] * K_im_j[k] )
            detail::tpl_inner_product<N_stages>( butcher_im.b2, K_im_j, ui, dt, unp1_bis );
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

        /**
         * @brief set absolute tolerance in chained config
         *
         * @param tol_ tolerance
         * @return auto& returns this object
         */
        template <typename tab_t = tableau_t>
            requires std::same_as<tab_t, tableau_t> && is_embedded
        auto&
        abs_tol( value_t tol_ )
        {
            info().absolute_tolerance = tol_;
            return *this;
        }

        /**
         * @brief set relative tolerance in chained config
         *
         * @param tol_ tolerance
         * @return auto& returns this object
         */
        template <typename tab_t = tableau_t>
            requires std::same_as<tab_t, tableau_t> && is_embedded
        auto&
        rel_tol( value_t tol_ )
        {
            info().relative_tolerance = tol_;
            return *this;
        }

        /**
         * @brief set tolerance for Newton method (for default Newton method)
         *
         * @param tol_ tolerance
         * @return auto& returns this object
         */
        auto&
        newton_tol( value_t tol_ )
        {
            tol = tol_;
            return *this;
        }

        /**
         * @brief set maximum of iterations for Newton method (for default Newton method)
         *
         * @param max_iter_ maximum of iterations
         * @return auto& returns this object
         */
        auto&
        newton_max_iter( std::size_t max_iter_ )
        {
            max_iter = max_iter_;
            return *this;
        }

        double tol           = ponio::default_config::newton_tolerance;      // tolerance of Newton method
        std::size_t max_iter = ponio::default_config::newton_max_iterations; // max iterations of Newton method

        linear_algebra_t linalg;
        iteration_info<tableau_t> _info;
    };

} // namespace ponio::runge_kutta::additive_runge_kutta
