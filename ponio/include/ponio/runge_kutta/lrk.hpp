// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private, include "../runge_kutta.hpp"

#pragma once

#include <concepts>
#include <cstddef>
#include <string_view>

#include "../butcher_tableau.hpp"
#include "../detail.hpp"
#include "../iteration_info.hpp"
#include "../ponio_config.hpp"
#include "../stage.hpp"

namespace ponio::runge_kutta::lawson_runge_kutta
{

    template <typename exp_t>
    struct lawson_base
    {
        exp_t m_exp;

        lawson_base( exp_t exp_ )
            : m_exp( exp_ )
        {
        }
    };

    template <typename tableau_t, typename exp_t>
    struct explicit_runge_kutta : public lawson_base<exp_t>
    {
        using lawson_base<exp_t>::m_exp;

        tableau_t butcher;
        static constexpr std::size_t N_stages = tableau_t::N_stages;
        static constexpr bool is_embedded     = butcher::is_embedded_tableau<tableau_t>;
        static constexpr std::size_t order    = tableau_t::order;
        static constexpr std::string_view id  = tableau_t::id;

        using value_t = typename tableau_t::value_t;

        explicit_runge_kutta( exp_t exp_, double tolerance = ponio::default_config::tol )
            : lawson_base<exp_t>( exp_ )
            , butcher()
            , _info( tolerance )
        {
            _info.number_of_eval = N_stages;
        }

        template <typename problem_t, typename state_t, typename value_t, typename array_ki_t, std::size_t i>
        state_t
        stage( Stage<i>, problem_t& pb, value_t tn, state_t& un, array_ki_t const& Ki, value_t dt )
        {
            return m_exp( -butcher.c[i] * dt * pb.l )
                 * pb.n( tn + butcher.c[i] * dt,
                     m_exp( butcher.c[i] * dt * pb.l ) * detail::tpl_inner_product<i>( butcher.A[i], Ki, un, dt ) );
        }

        template <typename problem_t, typename state_t, typename value_t, typename array_ki_t>
        state_t
        stage( Stage<N_stages>, problem_t& pb, value_t, state_t& un, array_ki_t const& Ki, value_t dt )
        {
            return m_exp( dt * pb.l ) * detail::tpl_inner_product<N_stages>( butcher.b, Ki, un, dt );
        }

        template <typename problem_t, typename state_t, typename value_t, typename array_ki_t, typename tab_t = tableau_t>
            requires std::same_as<tab_t, tableau_t> && is_embedded
        state_t
        stage( Stage<N_stages + 1>, problem_t& pb, value_t, state_t& un, array_ki_t const& Ki, value_t dt )
        {
            return m_exp( dt * pb.l ) * detail::tpl_inner_product<N_stages>( butcher.b2, Ki, un, dt );
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

        iteration_info<tableau_t> _info;
    };

    /**
     * factory of Lawson method to help clang++ which doesn't support C++20
     *
     * @tparam tableau_t type of Butcher tableau
     * @tparam exp_t type of exponential function
     * @param exp_ exponential function
     * @param tol_ tolenrence of method for adaptative time step integrator
     */
    template <typename tableau_t, typename exp_t>
    auto
    make_lawson( exp_t exp_, double tol_ = ponio::default_config::tol )
    {
        return explicit_runge_kutta<tableau_t, exp_t>( exp_, tol_ );
    }

} // namespace ponio::runge_kutta::lawson_runge_kutta
