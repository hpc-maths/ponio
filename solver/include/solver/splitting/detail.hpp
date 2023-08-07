// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <cstddef>

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

} // namespace ponio::splitting::detail
