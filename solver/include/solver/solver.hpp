// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <cmath>
#include <iterator>
#include <limits>

#include "method.hpp"
#include "time_span.hpp"

namespace ode
{

    /**
     * @brief store \f$(t^n, u^n, \Delta t^n)\f$
     *
     * @tparam value_t type of time
     * @tparam state_t type of solution \f$u^n\f$
     */
    template <typename value_t, typename state_t>
    struct current_solution
    {
        value_t time;
        state_t state;
        value_t time_step;

        current_solution( value_t tn, state_t un, value_t dt )
            : time( tn )
            , state( un )
            , time_step( dt )
        {
        }
    };

    template <typename value_t, typename state_t, typename method_t, typename problem_t>
    struct time_iterator
    {
        using difference_type   = value_t;
        using value_type        = current_solution<value_t, state_t>;
        using pointer           = value_type*;
        using reference         = value_type&;
        using const_pointer     = value_type const*;
        using const_reference   = value_type const&;
        using iterator_category = std::output_iterator_tag;

        value_type sol;
        method_t meth;
        problem_t& pb;
        ponio::time_span<value_t> t_span;
        typename ponio::time_span<value_t>::iterator it_next_time;
        value_t final_time; // should be a ponio::time_span;
        static constexpr value_t sentinel = std::numeric_limits<value_t>::max();

        // time_iterator ( problem_t & pb_, method_t meth_, state_t const& u0, ponio::time_span<value_t> && times, value_t dt )
        // : sol(times.front(), u0, dt)
        // , meth(meth_)
        // , pb(pb_)
        // , t_span(times)
        // , it_next_time(++times.begin())
        // {
        // }

        time_iterator( problem_t& pb_, method_t meth_, state_t const& u0, ponio::time_span<value_t> const& t_span_, value_t dt )
            : sol( ( t_span_.front() == t_span_.back() ) ? sentinel : t_span_.front(), u0, dt )
            , meth( meth_ )
            , pb( pb_ )
            , t_span( t_span_ )
            , it_next_time( std::begin( t_span ) )
            , final_time( t_span.back() )
        {
        }

        time_iterator( time_iterator const& rhs )
            : sol( rhs.sol )
            , meth( rhs.meth )
            , pb( rhs.pb )
            , t_span( rhs.t_span )
            , it_next_time( std::begin( t_span ) + ( rhs.it_next_time - std::begin( rhs.t_span ) ) )
            , final_time( rhs.final_time )
        {
        }

        void
        increment()
        {
            std::tie( sol.time, sol.state, sol.time_step ) = meth( pb, sol.time, sol.state, sol.time_step );
        }

        difference_type
        next_time() const
        {
            return sol.time + sol.time_step;
        }

        time_iterator&
        operator++()
        {
            // if ( next_time() > *it_next_time )
            // {
            //   sol.time_step = final_time - sol.time;
            //   ++it_next_time;
            // }
            if ( sol.time == final_time )
            {
                sol.time = sentinel;
            }

            if ( sol.time == sentinel )
            {
                return *this;
            }

            auto save_dt   = sol.time_step;
            bool change_dt = false;
            if ( next_time() > *it_next_time )
            {
                sol.time_step = *it_next_time - sol.time;
                ++it_next_time;
                change_dt = true;
            }

            increment();
            if ( change_dt )
            {
                sol.time_step = save_dt;
            }
            return *this;
        }

        time_iterator
        operator++( int )
        {
            auto copy = *this;
            ++( *this );
            return copy;
        }

        time_iterator&
        operator+=( difference_type dt )
        {
            sol.time_step = dt;
            increment();
            return *this;
        }

        reference
        operator*()
        {
            return sol;
        }

        const_reference
        operator*() const
        {
            return sol;
        }

        pointer
        operator->()
        {
            return &sol;
        }

        const_pointer
        operator->() const
        {
            return &sol;
        }
    };

    /**
     * @brief equality operator that compares only current time
     */
    template <typename value_t, typename state_t, typename method_t, typename problem_t>
    bool
    operator==( time_iterator<value_t, state_t, method_t, problem_t> const& lhs,
        time_iterator<value_t, state_t, method_t, problem_t> const& rhs )
    {
        return lhs.sol.time == rhs.sol.time;

        //( std::abs( lhs.sol.time - rhs.sol.time ) <= std::numeric_limits<value_t>::epsilon() * std::abs( std::min( lhs.sol.time,
        // rhs.sol.time ) ) * 1 );
    }

    /**
     * @brief three-way comparison operator that compares only current time
     */
    template <typename value_t, typename state_t, typename method_t, typename problem_t>
    auto
    operator<=>( time_iterator<value_t, state_t, method_t, problem_t> const& lhs,
        time_iterator<value_t, state_t, method_t, problem_t> const& rhs )
    {
        return lhs.sol.time <=> rhs.sol.time;
    }

    /**
     * @brief factory of `time_iterator`
     *
     * @tparam value_t   type of time and time step
     * @tparam state_t   type of solution \f$u^n\f$
     * @tparam method_t  type of numerical method to solve problem
     * @tparam problem_t type of problem
     *
     * @param pb           problem to solve, it could be any function of functor with a call operator with following parameter `(value_t tn,
     * state_t const& un)`
     * @param meth         choosen method to solve the problem `pb`
     * @param u0           initial condition \f$u_0 = u(t=0)\f$
     * @param t_init       initial time to start the range
     * @param dt           time step value \f$\Delta t\f$
     * @param t_final      final time to stop the range
     * @return auto
     */
    template <typename value_t, typename state_t, typename method_t, typename problem_t>
    auto
    make_time_iterator( problem_t& pb, method_t meth, state_t const& u0, ponio::time_span<value_t> const& t_span, value_t dt )
    {
        return time_iterator<value_t, state_t, method_t, problem_t>( pb, meth, u0, t_span, dt );
    }

    /**
     * @brief structure that stores begin and end `time_iterator` of a range of solutions
     *
     * @tparam value_t   type of time and time step
     * @tparam state_t   type of solution \f$u^n\f$
     * @tparam method_t  type of numerical method to solve problem
     * @tparam problem_t type of problem
     */
    template <typename value_t, typename state_t, typename method_t, typename problem_t>
    struct solver_range
    {
        using iterator_type = time_iterator<value_t, state_t, method_t, problem_t>;
        iterator_type _begin;
        iterator_type _end;

        solver_range( iterator_type const& first, iterator_type const& last )
            : _begin( first )
            , _end( last )
        {
        }

        inline auto&
        begin()
        {
            return _begin;
        }

        inline auto const&
        cbegin() const
        {
            return _begin;
        }

        inline auto const&
        begin() const
        {
            return cbegin();
        }

        inline auto&
        end()
        {
            return _end;
        }

        inline auto const&
        cend() const
        {
            return _end;
        }

        inline auto const&
        end() const
        {
            return cend();
        }
    };

    /**
     * @brief makes a range of solutions at each time
     *
     * @tparam value_t     type of time and time step
     * @tparam state_t     type of solution \f$u^n\f$
     * @tparam algorithm_t type of numerical method to solve problem
     * @tparam problem_t   type of problem
     *
     * @param pb           problem to solve, it could be any function of functor with a call operator with following parameter `(value_t tn,
     * state_t const& un)`
     * @param algo         choosen method to solve the problem `pb`
     * @param u0           initial condition \f$u_0 = u(t=0)\f$
     * @param t_span       container \f$[t_\text{start} , t_\text{end}]\f$ with possible intermediate time value where solver should go
     * @param dt           time step value \f$\Delta t\f$
     * @return returns the range with all solutions
     */
    template <typename value_t, typename state_t, typename algorithm_t, typename problem_t>
    auto
    make_solver_range( problem_t& pb, algorithm_t&& algo, state_t const& u0, ponio::time_span<value_t> const& t_span, value_t dt )
    {
        auto meth = ode::make_method( algo, u0 );

        auto first = make_time_iterator( pb, meth, u0, t_span, dt );
        auto last  = make_time_iterator( pb, meth, u0, { t_span.back() }, dt );

        return solver_range<value_t, state_t, decltype( meth ), problem_t>( first, last );
    }

    /**
     * @brief solve a problem on a specific time range with a specific method
     *
     * @param pb     problem to solve, it could by any function or functor with a call operator with following parameter `(value_t tn,
     * state_t const& un)`.
     * @param algo   choosen method to solve the problem `pb`
     * @param u0     initial condition \f$u_0 = u(t=0)\f$
     * @param t_span container \f$[t_\text{start} , t_\text{end}]\f$ with possible intermediate time value where solver should go
     * @param dt     time step value \f$\Delta t\f$
     * @param obs    observer that do something with current time, solution and time step at each iteration
     * @return returns the last value of solution \f$u^n\f$
     */
    template <typename Problem_t, typename Algorithm_t, typename state_t, typename value_t, typename Observer_t>
    state_t
    solve( Problem_t& pb, Algorithm_t&& algo, state_t const& u0, ponio::time_span<value_t> const& t_span, value_t dt, Observer_t&& obs )
    {
        value_t current_time = t_span.front();
        auto it_next_time    = t_span.begin() + 1;
        auto it_end_time     = t_span.end();

        value_t current_dt = dt;
        bool reset_dt      = false;

        state_t un  = u0;
        state_t un1 = u0;

        value_t last_time = t_span.back();

        auto meth = make_method( algo, un );

        obs( current_time, un, dt );

        while ( current_time < last_time )
        {
            std::tie( current_time, un1, current_dt ) = meth( pb, current_time, un, current_dt );
            std::swap( un, un1 );

            obs( current_time, un, current_dt );

            // prepare next step
            if ( current_time + current_dt > *it_next_time )
            {
                current_dt = *it_next_time - current_time;
                reset_dt   = true;
                ++it_next_time;
            }
            else if ( reset_dt )
            {
                reset_dt   = false;
                current_dt = dt;
            }
        }

        return un;
    }

}
