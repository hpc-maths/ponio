// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private

#pragma once

#include <concepts>
#include <cstddef>
#include <string_view> // NOLINT(misc-include-cleaner)

#include "../butcher_tableau.hpp"
#include "../detail.hpp"
#include "../iteration_info.hpp"
#include "../ponio_config.hpp"
#include "../stage.hpp" // NOLINT(misc-include-cleaner)

namespace ponio::runge_kutta::explicit_runge_kutta
{

    template <typename tableau_t>
    struct explicit_runge_kutta
    {
        tableau_t butcher;
        static constexpr std::size_t N_stages = tableau_t::N_stages;
        static constexpr bool is_embedded     = butcher::is_embedded_tableau<tableau_t>;
        static constexpr std::size_t order    = tableau_t::order;
        static constexpr std::string_view id  = tableau_t::id;

        using value_t = typename tableau_t::value_t;

        explicit_runge_kutta( double tolerance = default_config::tol )
            : butcher()
            , _info( tolerance )
        {
            _info.number_of_eval = N_stages;
        }

        template <typename problem_t, typename state_t, typename array_kj_t, std::size_t I>
        void
        stage( Stage<I>, problem_t& f, value_t tn, state_t& un, array_kj_t const& Kj, value_t dt, state_t& ui, state_t& Ki )
        {
            // ui = un + dt*sum(butcher.A[I]*Kj)
            detail::tpl_inner_product<I>( butcher.A[I], Kj, un, dt, ui );

            // Ki = f(tn + butcher.c[I]*dt, ui)
            f( tn + butcher.c[I] * dt, ui, Ki );
        }

        template <typename problem_t, typename state_t, typename array_kj_t>
        void
        stage( Stage<N_stages>, problem_t&, value_t, state_t& un, array_kj_t const& Kj, value_t dt, state_t&, state_t& Ki )
        {
            // Ki = un + dt*sum(butcher.b*Kj)
            detail::tpl_inner_product<N_stages>( butcher.b, Kj, un, dt, Ki );
        }

        template <typename problem_t, typename state_t, typename array_kj_t, typename tab_t = tableau_t>
            requires std::same_as<tab_t, tableau_t> && is_embedded
        void
        stage( Stage<N_stages + 1>, problem_t&, value_t, state_t& un, array_kj_t const& Kj, value_t dt, state_t&, state_t& Ki )
        {
            // Ki = un + dt*sum(butcher.b2*Kj)
            detail::tpl_inner_product<N_stages>( butcher.b2, Kj, un, dt, Ki );
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
         * @param tol tolerance
         * @return auto& returns this object
         */
        template <typename tab_t = tableau_t>
            requires std::same_as<tab_t, tableau_t> && is_embedded
        auto&
        abs_tol( value_t tol )
        {
            info().absolute_tolerance = tol;
            return *this;
        }

        /**
         * @brief set relative tolerance in chained config
         *
         * @param tol tolerance
         * @return auto& returns this object
         */
        template <typename tab_t = tableau_t>
            requires std::same_as<tab_t, tableau_t> && is_embedded
        auto&
        r_tol( value_t tol )
        {
            info().relative_tolerance = tol;
            return *this;
        }

        iteration_info<tableau_t> _info;
    };

} // namespace ponio::runge_kutta::explicit_runge_kutta
