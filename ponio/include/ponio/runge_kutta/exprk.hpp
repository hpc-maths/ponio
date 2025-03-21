// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private, include "../runge_kutta.hpp"

#pragma once

#include <concepts>
#include <cstddef>
#include <string_view> // NOLINT(misc-include-cleaner)
#include <type_traits> // NOLINT(misc-include-cleaner)
#include <utility>

#include "../butcher_tableau.hpp"
#include "../iteration_info.hpp"
#include "../ponio_config.hpp"
#include "../stage.hpp" // NOLINT(misc-include-cleaner)

namespace ponio::runge_kutta::exponential_runge_kutta
{

    namespace detail
    {
        template <typename func_t, typename linear_t>
            requires std::invocable<func_t, linear_t>
        [[maybe_unused]] auto
        coefficient_eval( func_t&& f, linear_t&& l )
        {
            return f( std::forward<linear_t>( l ) );
        }

        template <typename value_t, typename linear_t>
            requires std::is_arithmetic_v<std::remove_cvref_t<value_t>>
        [[maybe_unused]] auto
        coefficient_eval( value_t&& val, linear_t&& )
        {
            return std::forward<value_t>( val );
        }

        template <std::size_t I, std::size_t J, typename tuple_t>
        auto
        triangular_get( tuple_t& t )
        {
            return std::get<I*( I - 1 ) / 2 + J>( t );
        }

        template <std::size_t I, typename state_t, typename value_t, typename linear_t, typename tuple_t, typename array_t, std::size_t... Is>
        constexpr state_t
        tpl_inner_product_impl( tuple_t const& a,
            array_t const& k,
            state_t const& init,
            [[maybe_unused]] linear_t const& linear_part,
            [[maybe_unused]] value_t mul_coeff,
            std::index_sequence<Is...> )
        {
            return (
                init + ...
                + ( mul_coeff * coefficient_eval( triangular_get<I, Is>( a ), mul_coeff * linear_part ) * ( k[Is] + linear_part * init ) ) );
        }

        template <std::size_t I, typename state_t, typename value_t, typename linear_t, typename tuple_t, typename array_t>
        constexpr state_t
        tpl_inner_product( tuple_t const& a, array_t const& k, state_t const& init, linear_t const& linear_part, value_t mul_coeff )
        {
            return tpl_inner_product_impl<I>( a, k, init, linear_part, mul_coeff, std::make_index_sequence<I>() );
        }

        template <typename state_t, typename value_t, typename linear_t, typename tuple_t, typename array_t, std::size_t... Is>
        constexpr state_t
        tpl_inner_product_b_impl( tuple_t const& b,
            array_t const& k,
            state_t const& init,
            linear_t const& linear_part,
            value_t mul_coeff,
            std::index_sequence<Is...> )
        {
            return ( init + ...
                     + ( mul_coeff * coefficient_eval( std::get<Is>( b ), mul_coeff * linear_part ) * ( k[Is] + linear_part * init ) ) );
        }

        template <std::size_t I, typename state_t, typename value_t, typename linear_t, typename tuple_t, typename array_t>
        constexpr state_t
        tpl_inner_product_b( tuple_t const& b, array_t const& k, state_t const& init, linear_t const& linear_part, value_t mul_coeff )
        {
            return tpl_inner_product_b_impl( b, k, init, linear_part, mul_coeff, std::make_index_sequence<I>() );
        }
    } // namespace detail

    template <typename tableau_t>
    struct explicit_exp_rk_butcher
    {
        tableau_t butcher;
        static constexpr std::size_t N_stages = tableau_t::N_stages;
        static constexpr bool is_embedded     = butcher::is_embedded_tableau<tableau_t>;
        static constexpr std::size_t order    = tableau_t::order;
        static constexpr std::string_view id  = tableau_t::id;

        using value_t = typename tableau_t::value_t;

        iteration_info<tableau_t> _info;

        explicit_exp_rk_butcher( double tolerance = ponio::default_config::tol )
            : butcher()
            , _info( tolerance )
        {
            _info.number_of_eval = N_stages;
        }

        template <typename problem_t, typename state_t, typename value_t, typename array_ki_t, std::size_t i>
        state_t
        stage( Stage<i>, problem_t& pb, value_t tn, state_t& un, array_ki_t const& Ki, value_t dt )
        {
            return pb.n( tn + butcher.c[i] * dt, detail::tpl_inner_product<i>( butcher.a, Ki, un, pb.l, dt ) );
        }

        template <typename problem_t, typename state_t, typename value_t, typename array_ki_t>
        state_t
        stage( Stage<N_stages>, problem_t& pb, value_t, state_t& un, array_ki_t const& Ki, value_t dt )
        {
            return detail::tpl_inner_product_b<N_stages>( butcher.b, Ki, un, pb.l, dt );
        }

        template <typename problem_t, typename state_t, typename value_t, typename array_ki_t, typename tab_t = tableau_t>
            requires std::same_as<tab_t, tableau_t> && is_embedded
        state_t
        stage( Stage<N_stages + 1>, problem_t& pb, value_t, state_t& un, array_ki_t const& Ki, value_t dt )
        {
            return detail::tpl_inner_product_b<N_stages>( butcher.b2, Ki, un, pb.l, dt );
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
    };

} // namespace ponio::runge_kutta::exponential_runge_kutta
