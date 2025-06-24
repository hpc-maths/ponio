// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <cmath>
#include <iterator>
#include <limits>
#include <optional>

#include "method.hpp"
#include "time_span.hpp"

namespace ponio
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
            , state( std::move( un ) )
            , time_step( dt )
        {
        }

        current_solution( value_t tn )
            : time( tn )
            , time_step( static_cast<value_t>( 0. ) )
        {
        }
    };

    /**
     * @brief iterator on ponio::solver_range
     *
     * @tparam value_t   type of time value
     * @tparam state_t   type of solution \f$u^n\f$
     * @tparam method_t  type of object used to iterate (an instance of ponio::method)
     * @tparam problem_t type of callable object problem \f$f: t, u \mapsto f(t, u)\f$
     */
    template <typename value_t, typename state_t, typename method_t, typename problem_t>
    struct time_iterator
    {
        /**
         * @typedef difference_type
         * @brief type of difference between two iterators, returns difference of current time of each iterator
         *
         * @typedef value_type
         * @brief type of stored value
         *
         * @typedef pointer
         * @brief type of pointer on stored value
         *
         * @typedef reference
         * @brief type of reference on stored value
         *
         * @typedef const_pointer
         * @brief type of constant pointer on stored value
         *
         * @typedef const_reference
         * @brief type of constant reference on stored value
         *
         * @typedef iterator_category
         * @brief specify the category of iterator, here corresponds to output iterator (see [iterator
         * tags](https://en.cppreference.com/w/cpp/iterator/iterator_tags) for more information)
         */

        using difference_type   = value_t;
        using value_type        = current_solution<value_t, state_t>;
        using pointer           = value_type*;
        using reference         = value_type&;
        using const_pointer     = value_type const*;
        using const_reference   = value_type const&;
        using iterator_category = std::output_iterator_tag;

        value_type sol;
        method_t meth;
        problem_t& pb; // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
        ponio::time_span<value_t> t_span;
        typename ponio::time_span<value_t>::iterator it_next_time;
        std::optional<value_t> dt_reference;
        static constexpr value_t sentinel = std::numeric_limits<value_t>::max();

        /**
         * @brief Construct a new time iterator object
         *
         * @param pb_     problem object that represent \f$f: t, u \mapsto f(t, u)\f$
         * @param meth_   ponio::method object
         * @param u0      initial value of state \f$u\f$
         * @param t_span_ vector of initial time and final time with optional intermediate time
         * @param dt      time step value \f$\Delta t\f$
         *
         * @note In user interface ponio::time_iterator object is build by ponio::solver_range member functions.
         */
        time_iterator( problem_t& pb_, method_t&& meth_, state_t const& u0, ponio::time_span<value_t> const& t_span_, value_t dt )
            : sol( ( t_span_.front() == t_span_.back() ) ? sentinel : t_span_.front(), u0, dt )
            , meth( std::move( meth_ ) )
            , pb( pb_ )
            , t_span( t_span_ )
            , it_next_time( std::next( std::begin( t_span ) ) )
            , dt_reference( std::nullopt )
        {
        }

        time_iterator( problem_t& pb_, value_t t_end )
            : sol( sentinel )
            , pb( pb_ )
            , t_span( { t_end } )
            , it_next_time( std::end( t_span ) )
            , dt_reference( std::nullopt )
        {
        }

        time_iterator() = delete;

        /**
         * @brief Copy constructor of time_iterator
         *
         * @param rhs
         */
        time_iterator( time_iterator const& rhs )
            : sol( rhs.sol )
            , meth( rhs.meth )
            , pb( rhs.pb )
            , t_span( rhs.t_span )
            , it_next_time( std::begin( t_span ) + std::ranges::distance( std::begin( rhs.t_span ), rhs.it_next_time ) )
            , dt_reference( rhs.dt_reference )
        {
        }

        /**
         * @brief Move constructor of time_iterator
         *
         * @param rhs
         */
        time_iterator( time_iterator&& rhs ) noexcept
            : sol( std::move( rhs.sol ) )
            , meth( std::move( rhs.meth ) )
            , pb( std::move( rhs.pb ) )
            , t_span( std::move( rhs.t_span ) )
            , it_next_time( std::move( rhs.it_next_time ) )
            , dt_reference( std::move( rhs.dt_reference ) )
        {
        }

        /**
         * @brief Equality operator of time_iterator
         *
         * @param rhs
         * @return time_iterator&
         */
        time_iterator&
        operator=( time_iterator const& rhs )
        {
            if ( this != &rhs )
            {
                sol          = rhs.sol;
                meth         = rhs.meth;
                pb           = rhs.pb;
                t_span       = rhs.t_span;
                it_next_time = std::begin( t_span ) + std::ranges::distance( std::begin( rhs.t_span ), rhs.it_next_time );
                dt_reference = rhs.dt_reference;
            }

            return *this;
        }

        /**
         * @brief Equality operator of time_iterator
         *
         * @param rhs
         * @return time_iterator&
         */
        time_iterator&
        operator=( time_iterator&& rhs ) noexcept
        {
            if ( this != &rhs )
            {
                sol          = std::move( rhs.sol );
                meth         = std::move( rhs.meth );
                pb           = std::move( rhs.pb );
                t_span       = std::move( rhs.t_span );
                it_next_time = std::move( rhs.it_next_time );
                dt_reference = std::move( rhs.dt_reference );
            }

            return *this;
        }

        ~time_iterator() = default;

        /**
         * @brief increment current state by current time step
         *
         * @details \f$(t^n, u^n, \Delta t^n ) \gets \phi(t^n, u^n, \Delta t^n)\f$ where \f$\phi\f$ represents the method.
         */
        void
        increment()
        {
            std::tie( sol.time, sol.state, sol.time_step ) = meth( pb, sol.time, sol.state, sol.time_step );
        }

        /**
         * @brief get the next time value \f$t^n + \Delta t^n\f$
         *
         * @return difference_type
         */
        difference_type
        next_time() const
        {
            return sol.time + sol.time_step;
        }

        /**
         * @brief increment properly time_iterator in the solver_range (take care of end point and optionally middle points in given t_span)
         *
         */
        time_iterator&
        operator++()
        {
            if ( sol.time == t_span.back() )
            {
                sol.time = sentinel;
            }

            if ( sol.time == sentinel )
            {
                return *this;
            }

            if ( dt_reference.has_value() )
            {
                sol.time_step = dt_reference.value();
                dt_reference  = std::nullopt;
            }

            if ( next_time() > *it_next_time )
            {
                dt_reference  = sol.time_step;
                sol.time_step = *it_next_time - sol.time;
                ++it_next_time;
            }

            increment();
            return *this;
        }

        /**
         * @brief post fix increment time_iterator
         *
         */
        time_iterator // NOLINT(cert-dcl21-cpp): This check is deprecated since it’s no longer part of the CERT standard
        operator++( int )
        {
            auto copy = *this;
            ++( *this );
            return copy;
        }

        /**
         * @brief dereference time_iterator and get current solution data member
         */
        reference
        operator*()
        {
            return sol;
        }

        /**
         * @brief dereference time_iterator and get current solution data member
         */
        const_reference
        operator*() const
        {
            return sol;
        }

        /**
         * @brief accessor to current solution data member
         */
        pointer
        operator->()
        {
            return &sol;
        }

        /**
         * @brief accessor to current solution data member
         */
        const_pointer
        operator->() const
        {
            return &sol;
        }

        /**
         * @brief accessor to informations on iteration with algorithm
         *
         * @return auto&
         */
        auto&
        info() // cppcheck-suppress unusedFunction
        {
            return meth.info();
        }

        /**
         * @brief accessor to informations on iteration with algorithm
         *
         * @return auto const&
         */
        auto const&
        info() const // cppcheck-suppress unusedFunction
        {
            return meth.info();
        }

        /**
         * @brief returns stages if need to access to them to resize or change condition on them
         *
         * @return auto&
         */
        auto&
        stages() // cppcheck-suppress unusedFunction
        {
            return meth.stages();
        }

        /**
         * @brief returns stages if need to access to them
         *
         * @return auto const&
         */
        auto const&
        stages() const // cppcheck-suppress unusedFunction
        {
            return meth.stages();
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
     * @param t_span       vector of initial time and final time with optional intermediate time
     * @param dt           time step value \f$\Delta t\f$
     * @return auto
     */
    template <typename value_t, typename state_t, typename method_t, typename problem_t>
    auto
    make_time_iterator( problem_t& pb, method_t&& meth, state_t const& u0, ponio::time_span<value_t> const& t_span, value_t dt )
    {
        return time_iterator<value_t, state_t, method_t, problem_t>( pb, std::forward<method_t>( meth ), u0, t_span, dt );
    }

    template <typename value_t_>
    struct time_sentinel_iterator
    {
        using value_t         = value_t_;
        using difference_type = value_t;

        value_t time = std::numeric_limits<value_t>::max();
    };

    /**
     * @brief equality operator that compares only current time
     */
    template <typename value_t>
    bool
    operator==( time_sentinel_iterator<value_t> const& lhs, time_sentinel_iterator<value_t> const& rhs )
    {
        return lhs.time == rhs.time;
    }

    /**
     * @brief three-way comparison operator that compares only current time
     */
    template <typename value_t>
    auto
    operator<=>( time_sentinel_iterator<value_t> const& lhs, time_sentinel_iterator<value_t> const& rhs )
    {
        return lhs.time <=> rhs.time;
    }

    /**
     * @brief equality operator that compares only current time
     */
    template <typename value_t, typename state_t, typename method_t, typename problem_t>
    bool
    operator==( time_iterator<value_t, state_t, method_t, problem_t> const& lhs, time_sentinel_iterator<value_t> const& rhs )
    {
        return lhs.sol.time == rhs.time;
    }

    /**
     * @brief three-way comparison operator that compares only current time
     */
    template <typename value_t, typename state_t, typename method_t, typename problem_t>
    auto
    operator<=>( time_iterator<value_t, state_t, method_t, problem_t> const& lhs, time_sentinel_iterator<value_t> const& rhs )
    {
        return lhs.sol.time <=> rhs.time;
    }

    /**
     * @brief equality operator that compares only current time
     */
    template <typename value_t, typename state_t, typename method_t, typename problem_t>
    bool
    operator==( time_sentinel_iterator<value_t> const& lhs, time_iterator<value_t, state_t, method_t, problem_t> const& rhs )
    {
        return lhs.time == rhs.sol.time;
    }

    /**
     * @brief three-way comparison operator that compares only current time
     */
    template <typename value_t, typename state_t, typename method_t, typename problem_t>
    auto
    operator<=>( time_sentinel_iterator<value_t> const& lhs, time_iterator<value_t, state_t, method_t, problem_t> const& rhs )
    {
        return lhs.time <=> rhs.sol.time;
    }

    /**
     * @brief factory of ending `time_sentinel_iterator`
     *
     * @tparam value_t   type of time and time step
     *
     * @return auto
     */
    template <typename value_t>
    auto
    make_sentinel_iterator()
    {
        return time_sentinel_iterator<value_t>();
    }

    /**
     * @brief structure that stores begin and end time_iterator of a range of solutions
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
        using sentinel_type = time_sentinel_iterator<value_t>;
        iterator_type _begin;
        sentinel_type _end;

        /**
         * @brief Construct a new solver range object
         *
         * @param begin initial iterator on solver_range
         * @param end   end iterator on solver_range (sentinel)
         */
        solver_range( iterator_type const& begin, sentinel_type const& end )
            : _begin( begin )
            , _end( end )
        {
        }

        /**
         * @brief returns an iterator to the beginning solver_range
         *
         */
        auto&
        begin()
        {
            return _begin;
        }

        /**
         * @brief returns a constant iterator to the beginning solver_range
         *
         */
        auto const&
        cbegin() const
        {
            return _begin;
        }

        /**
         * @brief returns a constant iterator to the beginning solver_range
         *
         */
        auto const&
        begin() const
        {
            return cbegin();
        }

        /**
         * @brief returns an iterator to the ending solver_range
         *
         */
        auto&
        end()
        {
            return _end;
        }

        /**
         * @brief returns a constant iterator to the ending solver_range
         *
         */
        auto const&
        cend() const
        {
            return _end;
        }

        /**
         * @brief returns a constant iterator to the ending solver_range
         *
         */
        auto const&
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
        auto meth = make_method<value_t>( std::forward<algorithm_t>( algo ), u0 );

        auto begin = make_time_iterator( pb, std::move( meth ), u0, t_span, dt );
        auto end   = make_sentinel_iterator<value_t>();

        return solver_range<value_t, state_t, decltype( meth ), problem_t>( begin, end );
    }

    /**
     * @brief solve a problem on a specific time range with a specific method
     *
     * @param pb     problem to solve, it could by any function or functor with a call operator with following parameter `(value_t tn,
     * state_t& un)`.
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
        // TODO: change this function to use time_iterator
        value_t current_time = t_span.front();
        auto it_next_time    = t_span.begin() + 1;

        value_t current_dt = dt;
        bool reset_dt      = false;

        state_t un  = u0;
        state_t un1 = u0;

        value_t last_time = t_span.back();
        auto last_it      = t_span.end() - 1;

        auto meth = make_method<value_t>( std::forward<Algorithm_t>( algo ), un );

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
                if ( it_next_time != last_it )
                {
                    ++it_next_time;
                }
            }
            else if ( reset_dt )
            {
                reset_dt   = false;
                current_dt = dt;
            }
        }

        return un;
    }

} // namespace ponio
