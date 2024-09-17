// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private, include "../splitting.h"

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <tuple>

#include "../ponio_config.hpp"

namespace ponio::splitting::detail
{

    /**
     * tiny solve method only for splitting method, it solves a problem between `ti` and `tf` with the time step `dt` and returns
     * state at final time
     *
     * @tparam I index of subproblem to solve
     * @param pb   problem to solve
     * @param meth numerical method to solve problem `pb`
     * @param ui   initial state
     * @param ti   initial time
     * @param tf   final time
     * @param dt   time step
     */
    template <std::size_t I, typename Problem_t, typename Method_t, typename state_t, typename value_t>
    state_t
    _split_solve( Problem_t& pb, Method_t& meth, state_t& ui, value_t ti, value_t tf, value_t dt )
    {
        value_t current_dt   = std::min( dt, tf - ti );
        value_t current_time = ti;
        while ( current_time != tf )
        {
            std::tie( current_time, ui, current_dt ) = std::get<I>( meth )( std::get<I>( pb.system ), current_time, ui, current_dt );
            if ( current_time + current_dt > tf )
            {
                current_dt = tf - current_time;
            }
        }
        return ui;
    }

    // ---- class splitting_tuple ----------------------------------

    /** @class splitting_tuple
     *  generic class to prevent code duplication between splitting method
     *
     * @tparam _splitting_method_t type of splitting method
     * @tparam value_t             type of coefficient
     * @tparam optional_args_t     type of tuple of optional arguments (void if not needed)
     * @tparam Algorithms_t        type of algorithms to solve each step of splitting
     */
    template <template <typename, typename...> typename _splitting_method_t, typename value_t, typename optional_args_t, typename... Algorithms_t>
    struct splitting_tuple
    {
        using splitting_method_t                  = _splitting_method_t<value_t, Algorithms_t...>;
        static constexpr std::size_t order        = splitting_method_t::order;
        static constexpr bool is_splitting_method = true;
        static constexpr std::string_view id      = splitting_method_t::id;

        static constexpr bool has_optional_args = !std::is_void<optional_args_t>::value;
        using optional_args_container           = typename std::conditional<has_optional_args, optional_args_t, bool>::type;

        std::tuple<Algorithms_t...> algos;
        std::array<value_t, sizeof...( Algorithms_t )> time_steps;
        optional_args_container optional_arguments;

        template <bool _has_optional_args = has_optional_args>
            requires std::same_as<std::bool_constant<has_optional_args>, std::true_type>
        splitting_tuple( std::tuple<Algorithms_t...>&& algs, std::array<value_t, sizeof...( Algorithms_t )>&& dts, optional_args_container args )
            : algos( std::forward<std::tuple<Algorithms_t...>>( algs ) )
            , time_steps( std::forward<std::array<value_t, sizeof...( Algorithms_t )>>( dts ) )
            , optional_arguments( args )
        {
        }

        template <bool _has_optional_args = has_optional_args>
        splitting_tuple( std::tuple<Algorithms_t...>&& algs, std::array<value_t, sizeof...( Algorithms_t )>&& dts )
            : algos( std::forward<std::tuple<Algorithms_t...>>( algs ) )
            , time_steps( std::forward<std::array<value_t, sizeof...( Algorithms_t )>>( dts ) )
        {
        }
    };

    /**
     * @brief factory for generic splitting method (strang, lie)
     *
     * @tparam _splitting_method_t type of splitting method
     * @tparam value_t             type of coefficients
     * @tparam Methods_t           type of methods to solve each step of splitting
     * @param meths tuple of methods
     * @param dts   time step for each method
     */
    template <template <typename, typename...> typename _splitting_method_t, typename value_t, typename... Methods_t>
    auto
    make_splitting_from_tuple( std::tuple<Methods_t...> const& meths, std::array<value_t, sizeof...( Methods_t )> const& dts )
    {
        return _splitting_method_t<value_t, Methods_t...>( meths, dts );
    }

    /**
     * @brief factory for generic splitting method (adaptive_strang)
     *
     * @tparam _splitting_method_t type of splitting method
     * @tparam value_t             type of coefficients
     * @tparam optional_tuple_t    type of tuple of optional arguments
     * @tparam Methods_t           type of methods to solve each step of splitting
     * @param meths         tuple of methods
     * @param dts           time step for each method
     * @param optional_args tuple of optional arguments
     */
    template <template <typename, typename...> typename _splitting_method_t, typename value_t, typename optional_tuple_t, typename... Methods_t>
    auto
    make_splitting_from_tuple( std::tuple<Methods_t...> const& meths,
        std::array<value_t, sizeof...( Methods_t )> const& dts,
        optional_tuple_t optional_args )
    {
        return std::apply(
            [&]<typename... Args_t>( Args_t&&... args )
            {
                return _splitting_method_t<value_t, Methods_t...>( meths, dts, args... );
            },
            optional_args );
    }

    // ---- class splitting_base -----------------------------------

    /**
     * @brief parent class of splitting methods, store list of methods and time steps
     *
     * @tparam _value_t  type of coefficients and time steps
     * @tparam methods_t list of types of methods to solve each sub-problem
     */
    template <typename _value_t, typename... methods_t>
    struct splitting_base
    {
        using value_t = _value_t;
        using tuple_t = std::tuple<methods_t...>;

        static constexpr bool is_splitting_method = true;
        static constexpr std::size_t N_methods    = sizeof...( methods_t );

        tuple_t methods;
        std::array<value_t, N_methods> time_steps;

        splitting_base( std::tuple<methods_t...> const& meths, std::array<value_t, sizeof...( methods_t )> const& dts );
    };

    /**
     * @brief Construct a new splitting base<value t, methods t...>::splitting base object
     *
     * @tparam value_t   type of coefficients and time steps
     * @tparam methods_t list of types of methods to solve each sub-problem
     * @param meths      list of methods to solve each sub-problem
     * @param dts        list of time step to solve each sub-problem (to iterate on each sub-step)
     */
    template <typename value_t, typename... methods_t>
    splitting_base<value_t, methods_t...>::splitting_base( std::tuple<methods_t...> const& meths,
        std::array<value_t, sizeof...( methods_t )> const& dts )
        : methods( meths )
        , time_steps( dts )
    {
    }

} // namespace ponio::splitting::detail
