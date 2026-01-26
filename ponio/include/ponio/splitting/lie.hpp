// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private

#pragma once

#include <array>
#include <cstddef>
#include <string_view> // NOLINT(misc-include-cleaner)
#include <tuple>
#include <utility>

#include "../iteration_info.hpp"
#include "../stage.hpp"
#include "detail.hpp"

namespace ponio::splitting::lie
{

    // ---- class lie -----------------------------------------------

    /** @class lie
     *  Lie splitting method
     *  @tparam value_t   type of time steps
     *  @tparam methods_t list of methods to solve each sub-problem
     */
    template <typename _value_t, typename... methods_t>
    struct lie : detail::splitting_base<_value_t, methods_t...>
    {
        using value_t = _value_t;
        using base_t  = detail::splitting_base<value_t, methods_t...>;

        using base_t::splitting_base;

        using base_t::is_splitting_method;
        using base_t::N_methods;

        using base_t::methods;
        using base_t::time_steps;

        static constexpr std::size_t order   = 1;
        static constexpr std::string_view id = "lie";
        static constexpr std::size_t N_steps = N_methods;

        iteration_info<lie> _info;

        lie( std::tuple<methods_t...> const& meths, std::array<value_t, N_methods> const& dts )
            : base_t( meths, dts )
            , _info( methods )
        {
        }

        // _call_inc can not be outside the class definition due to llvm bug
        // (see https://github.com/llvm/llvm-project/issues/56482)
        template <std::size_t I = 0, typename Problem_t, typename state_t>
            requires( I == N_steps )
        void _call_inc( Problem_t&, value_t, state_t&, value_t, state_t& )
        {
        }

        /**
         * @brief incremental call of each method of each subproblem
         *
         * @param f    problem to solve
         * @param tn   current time \f$\Delta t\f$
         * @param ui   initial solution for step `I`
         * @param dt   time step \f$\Delta t\f$
         * @param uip1 solution \f$\texttt{uip1} = \phi^{[i]}(\texttt{ui})\f$
         */
        template <std::size_t I = 0, typename Problem_t, typename state_t>
            requires( I < N_steps )
        void _call_inc( Problem_t& f, value_t tn, state_t& ui, value_t dt, state_t& uip1 )
        {
            if constexpr ( I == 0 )
            {
                _info.reset_eval();
            }

            detail::_split_solve<I>( f, methods, ui, tn, tn + dt, time_steps[I], uip1, _info );
            _call_inc<I + 1>( f, tn, uip1, dt, ui );
        }

        template <typename Problem_t, typename state_t>
        void
        operator()( Problem_t& f, value_t& tn, state_t& un, value_t& dt, state_t& unp1 );

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
         * @brief gets array of stages of the Ith method
         *
         * @tparam I index of step method
         */
        template <std::size_t I>
        auto&
        stages( sub_method<I> )
        {
            return std::get<I>( methods ).stages();
        }

        /**
         * @brief gets array of stages of the Ith method (constant version)
         *
         * @tparam I index of step method
         */
        template <std::size_t I>
        auto const&
        stages( sub_method<I> ) const
        {
            return std::get<I>( methods ).stages();
        }
    };

    /**
     * call operator to initiate Lie splitting recursion
     * @param f    \ref problem to solve
     * @param tn   current time \f$t^n\f$
     * @param un   current solution \f$u^n \approx u(t^n)\f$
     * @param dt   time step \f$\Delta t\f$
     * @param unp1 solution at time \f$t^{n+1} = t^n + \Delta t\f$
     */
    template <typename value_t, typename... methods_t>
    template <typename Problem_t, typename state_t>
    // auto
    void
    lie<value_t, methods_t...>::operator()( Problem_t& f, value_t& tn, state_t& un, value_t& dt, state_t& unp1 )
    {
        _call_inc( f, tn, un, dt, unp1 );

        if constexpr ( N_methods % 2 == 0 )
        {
            std::swap( un, unp1 );
        }

        tn = tn + dt;
    }

    // ---- *helper* ----

    /**
     * a helper factory for @ref ponio::splitting::detail::splitting_tuple from a tuple of algorithms to build a Lie method
     *
     * @tparam value_t      type of coefficients
     * @tparam algorithms_t variadic list of types of algorithms
     * @param args          variadic list of pairs of algorithm and time step
     * @return a @ref ponio::splitting::detail::splitting_tuple object build from the tuple of methods
     */
    template <typename value_t, typename... algorithms_t>
    auto
    make_lie_tuple( std::pair<algorithms_t, value_t>&&... args )
    {
        return detail::splitting_tuple<lie, value_t, void, algorithms_t...>( std::forward_as_tuple( ( args.first )... ), { args.second... } );
    }

} // namespace ponio::splitting::lie
