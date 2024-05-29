// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <array>
#include <cstddef>
#include <string_view>
#include <tuple>
#include <utility>

#include "detail.hpp"

namespace ponio::splitting::lie
{

    // ---- class lie -----------------------------------------------

    /** @class lie
     *  Lie splitting method
     *  @tparam value_t   type of time steps
     *  @tparam Methods_t list of methods to solve each sub-problem
     */
    template <typename value_t, typename... Methods_t>
    struct lie
    {
        static constexpr std::size_t order        = 1;
        static constexpr bool is_splitting_method = true;
        static constexpr std::string_view id      = "lie";

        std::tuple<Methods_t...> methods;
        std::array<value_t, sizeof...( Methods_t )> time_steps;

        lie( std::tuple<Methods_t...> const& meths, std::array<value_t, sizeof...( Methods_t )> const& dts );

        // _call_inc can not be outside the class definition due to llvm bug
        // (see https://github.com/llvm/llvm-project/issues/56482)
        template <std::size_t I = 0, typename Problem_t, typename state_t>
            requires( I == sizeof...( Methods_t ) )
        inline void _call_inc( Problem_t&, value_t, state_t&, value_t )
        {
        }

        /**
         * incremental call of each method of each subproblem
         * @tparam I solving step
         * @param f          \ref problem to solve
         * @param tn         current time \f$t^n\f$
         * @param[in,out] ui \f$\texttt{ui}=\phi_{\Delta t}^{[f_1]}\circ\cdots\circ\phi_{\Delta t}^{[f_{i-1}]}(t^n,u^n)\f$
         * @param dt         time step \f$\Delta t\f$
         * @details The parameter @p ui is update to \f$\phi_{\Delta t}^{[f_i]}(t^n,\texttt{ui})\f$
         */
        template <std::size_t I = 0, typename Problem_t, typename state_t>
            requires( I < sizeof...( Methods_t ) )
        inline void _call_inc( Problem_t& f, value_t tn, state_t& ui, value_t dt )
        {
            ui = detail::_split_solve<I>( f, methods, ui, tn, tn + dt, time_steps[I] );
            _call_inc<I + 1>( f, tn, ui, dt );
        }

        template <typename Problem_t, typename state_t>
        auto
        operator()( Problem_t& f, value_t tn, state_t const& un, value_t dt );
    };

    /**
     * constructor of \ref lie from a tuple
     */
    template <typename value_t, typename... Methods_t>
    lie<value_t, Methods_t...>::lie( std::tuple<Methods_t...> const& meths, std::array<value_t, sizeof...( Methods_t )> const& dts )
        : methods( meths )
        , time_steps( dts )
    {
    }

    /**
     * call operator to initiate Lie splitting recursion
     * @param f  \ref problem to solve
     * @param tn current time \f$t^n\f$
     * @param un current solution \f$u^n \approx u(t^n)\f$
     * @param dt time step \f$\Delta t\f$
     */
    template <typename value_t, typename... Methods_t>
    template <typename Problem_t, typename state_t>
    auto
    lie<value_t, Methods_t...>::operator()( Problem_t& f, value_t tn, state_t const& un, value_t dt )
    {
        state_t ui = un;
        _call_inc( f, tn, ui, dt );
        return std::make_tuple( tn + dt, ui, dt );
    }

    // ---- *helper* ----

    /**
     * a helper factory for \ref lie functor from a tuple of methods
     *
     * @tparam value_t   type of coefficients
     * @tparam Methods_t variadic list of types of methods
     * @param meths      list of methods
     * @param dts        associated time step foreach method
     * @return a \ref lie object build from the tuple of methods
     */
    template <typename value_t, typename... Methods_t>
    auto
    make_lie_from_tuple( std::tuple<Methods_t...> const& meths, std::array<value_t, sizeof...( Methods_t )> const& dts )
    {
        return lie<value_t, Methods_t...>( meths, dts );
    }

    // ---- class lie_tuple -----------------------------------------

    /** @class lie_tuple
     *  a helper to deduce method for ::ponio::make_method(splitting::lie::lie_tuple<Algorithms_t...> const &, state_t const &)
     *  @tparam value_t      type of time steps
     *  @tparam Algorithms_t variadic template of algorithms to solve each subproblem
     *  @details This is a dummy class to select correct \ref method to solve the problem
     */
    template <typename value_t, typename... Algorithms_t>
    struct lie_tuple
    {
        static constexpr std::size_t order        = 1;
        static constexpr bool is_splitting_method = true;
        static constexpr std::string_view id      = "lie";

        std::tuple<Algorithms_t...> algos;
        std::array<value_t, sizeof...( Algorithms_t )> time_steps;

        lie_tuple( std::tuple<Algorithms_t...>&& algs, std::array<value_t, sizeof...( Algorithms_t )>&& dts );
    };

    /**
     * constructor of \ref lie_tuple from a variadic number of algorithms
     */
    template <typename value_t, typename... Algorithms_t>
    inline lie_tuple<value_t, Algorithms_t...>::lie_tuple( std::tuple<Algorithms_t...>&& algs,
        std::array<value_t, sizeof...( Algorithms_t )>&& dts )
        : algos( std::forward<std::tuple<Algorithms_t...>>( algs ) )
        , time_steps( std::forward<std::array<value_t, sizeof...( Algorithms_t )>>( dts ) )
    {
    }

    // ---- *helper* ----

    /**
     * a helper factory for \ref lie_tuple from a tuple of algorithms
     *
     * @tparam value_t      type of coefficients
     * @tparam Algorithms_t variadic list of types of algorithms
     * @param args          variadic list of pairs of algorithm and time step
     * @return a \ref lie_tuple object build from the tuple of methods
     */
    template <typename value_t, typename... Algorithms_t>
    auto
    make_lie_tuple( std::pair<Algorithms_t, value_t>&&... args )
    {
        return lie_tuple<value_t, Algorithms_t...>( std::forward_as_tuple( ( args.first )... ), { args.second... } );
    }

} // namespace ponio::splitting::lie
