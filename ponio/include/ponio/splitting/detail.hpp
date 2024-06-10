// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include "../ponio_config.hpp"
#include <algorithm>
#include <cstddef>
#include <tuple>

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

    // ---- class _splitting_tuple ----------------------------------
    template <template <typename, typename...> typename _splitting_method_t, typename value_t, typename... Algorithms_t>
    struct _splitting_tuple
    {
        using splitting_method_t                  = _splitting_method_t<value_t, Algorithms_t...>;
        static constexpr std::size_t order        = splitting_method_t::order;
        static constexpr bool is_splitting_method = true;
        static constexpr std::string_view id      = splitting_method_t::id;

        std::tuple<Algorithms_t...> algos;
        std::array<value_t, sizeof...( Algorithms_t )> time_steps;
        std::array<value_t, 2> optional_arguments;

        _splitting_tuple( std::tuple<Algorithms_t...>&& algs,
            std::array<value_t, sizeof...( Algorithms_t )>&& dts,
            value_t delta     = 0.1,
            value_t tolerance = default_config::tol )
            : algos( std::forward<std::tuple<Algorithms_t...>>( algs ) )
            , time_steps( std::forward<std::array<value_t, sizeof...( Algorithms_t )>>( dts ) )
            , optional_arguments( { delta, tolerance } )
        {
        }
    };

    template <template <typename, typename...> typename _splitting_method_t, typename value_t, typename... Methods_t, typename... Args_t>
    auto
    make_splitting_from_tuple( std::tuple<Methods_t...> const& meths,
        std::array<value_t, sizeof...( Methods_t )> const& dts,
        Args_t... optional_args )
    {
        return _splitting_method_t<value_t, Methods_t...>( meths, dts, optional_args... );
    }

} // namespace ponio::splitting::detail
