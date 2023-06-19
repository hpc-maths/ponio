// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <array>
#include <concepts>
#include <tuple>
#include <type_traits>

namespace ode
{
    namespace splitting
    {

        namespace detail
        {
            /**
             * tiny solve method only for splitting method, it solves a problem betwen `ti` and `tf` with the time step `dt` and returns
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
        }

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

            std::tuple<Methods_t...> methods;
            std::array<value_t, sizeof...( Methods_t )> time_steps;

            lie( std::tuple<Methods_t...> const& meths, std::array<value_t, sizeof...( Methods_t )> const& dts );

            template <std::size_t I = 0, typename Problem_t, typename state_t>
                requires( I == sizeof...( Methods_t ) )
            inline void _call_inc( Problem_t& f, value_t tn, state_t& ui, value_t dt );

            template <std::size_t I = 0, typename Problem_t, typename state_t>
                requires( I < sizeof...( Methods_t ) )
            inline void _call_inc( Problem_t& f, value_t tn, state_t& ui, value_t dt );

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

        template <typename value_t, typename... Methods_t>
        template <std::size_t I, typename Problem_t, typename state_t>
            requires( I == sizeof...( Methods_t ) )
        inline void lie<value_t, Methods_t...>::_call_inc( Problem_t&, value_t, state_t&, value_t )
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
        template <typename value_t, typename... Methods_t>
        template <std::size_t I, typename Problem_t, typename state_t>
            requires( I < sizeof...( Methods_t ) )
        inline void lie<value_t, Methods_t...>::_call_inc( Problem_t& f, value_t tn, state_t& ui, value_t dt )
        {
            ui = detail::_split_solve<I>( f, methods, ui, tn, tn + dt, time_steps[I] );
            _call_inc<I + 1>( f, tn, ui, dt );
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

        /**
         * a helper factory for \ref lie functor from a tuple of methods
         * @param t tuple of \ref method
         * @return a \ref lie object build from the tuple of methods
         */
        template <typename value_t, typename... Methods_t>
        auto
        make_lie_from_tuple( std::tuple<Methods_t...> const& meths, std::array<value_t, sizeof...( Methods_t )> const& dts )
        {
            return lie<value_t, Methods_t...>( meths, dts );
        }

        /** @class lie_tuple
         *  a helper to deduce method for ::ode::make_method(splitting::lie_tuple<Algorithms_t...> const &, state_t const &)
         *  @tparam value_t      type of time steps
         *  @tparam Algorithms_t variadic template of algorithms to solve each subproblem
         *  @details This is a dummy class to select correct \ref method to solve the problem
         */
        template <typename value_t, typename... Algorithms_t>
        struct lie_tuple
        {
            static constexpr std::size_t order        = 1;
            static constexpr bool is_splitting_method = true;

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
            : algos( algs )
            , time_steps( dts )
        {
        }

        /**
         * a helper factory for \ref lie_tuple from a tuple of algorithms
         * @param a tuple of \ref method
         * @return a \ref lie_tuple object build from the tuple of methods
         */
        template <typename value_t, typename... Algorithms_t>
        auto
        make_lie_tuple( std::pair<Algorithms_t, value_t>&&... args )
        {
            return lie_tuple<value_t, Algorithms_t...>( std::forward_as_tuple( ( args.first )... ), { args.second... } );
        }

        /** @class strang
         *  Strang splitting method
         *  @tparam Methods_t list of methods to solve each sub-problem
         */
        template <typename value_t, typename... Methods_t>
        struct strang : lie<value_t, Methods_t...>
        {
            using lie<value_t, Methods_t...>::lie;
            using lie<value_t, Methods_t...>::methods;
            using lie<value_t, Methods_t...>::time_steps;

            static constexpr std::size_t order        = 2;
            static constexpr bool is_splitting_method = true;

            template <std::size_t I = 0, typename Problem_t, typename state_t>
                requires( I == sizeof...( Methods_t ) - 1 )
            inline void _call_inc( Problem_t& f, value_t tn, state_t& ui, value_t dt );

            template <std::size_t I = 0, typename Problem_t, typename state_t>
                requires( I < sizeof...( Methods_t ) - 1 )
            inline void _call_inc( Problem_t& f, value_t tn, state_t& ui, value_t dt );

            template <std::size_t I = sizeof...( Methods_t ) - 1, typename Problem_t, typename state_t>
                requires( I == 0 )
            inline void _call_dec( Problem_t& f, value_t tn, state_t& ui, value_t dt );

            template <std::size_t I = sizeof...( Methods_t ) - 1, typename Problem_t, typename state_t>
                requires( I > 0 )
            inline void _call_dec( Problem_t& f, value_t tn, state_t& ui, value_t dt );

            template <typename Problem_t, typename state_t>
            auto
            operator()( Problem_t& f, value_t tn, state_t const& un, value_t dt );
        };

        // end of incremental recursion, start of decremental recursion
        template <typename value_t, typename... Methods_t>
        template <std::size_t I, typename Problem_t, typename state_t>
            requires( I == sizeof...( Methods_t ) - 1 )
        inline void strang<value_t, Methods_t...>::_call_inc( Problem_t& f, value_t tn, state_t& ui, value_t dt )
        {
            ui = detail::_split_solve<I>( f, methods, ui, tn, tn + dt, time_steps[I] );
            _call_dec<I - 1>( f, tn, ui, dt );
        }

        /**
         * incremental call of each method of each subproblem
         * @tparam I solving step
         * @param f          \ref problem to solve
         * @param tn         current time \f$t^n\f$
         * @param[in,out] ui \f$\texttt{ui}=\phi_{^{\Delta t}/_2}^{[f_1]}\circ\cdots\circ\phi_{^{\Delta t}/_2}^{[f_{i-1}]}(t^n,u^n)\f$
         * @param dt         time step \f$\Delta t\f$
         * @details The parameter @p ui is update to \f$\phi_{^{\Delta t}/_2}^{[f_i]}(t^n,\texttt{ui})\f$
         */
        template <typename value_t, typename... Methods_t>
        template <std::size_t I, typename Problem_t, typename state_t>
            requires( I < sizeof...( Methods_t ) - 1 )
        inline void strang<value_t, Methods_t...>::_call_inc( Problem_t& f, value_t tn, state_t& ui, value_t dt )
        {
            ui = detail::_split_solve<I>( f, methods, ui, tn, tn + 0.5 * dt, time_steps[I] );
            _call_inc<I + 1>( f, tn, ui, dt );
        }

        // end of decremental recursion, end of recursion
        template <typename value_t, typename... Methods_t>
        template <std::size_t I, typename Problem_t, typename state_t>
            requires( I == 0 )
        inline void strang<value_t, Methods_t...>::_call_dec( Problem_t& f, value_t tn, state_t& ui, value_t dt )
        {
            ui = detail::_split_solve<I>( f, methods, ui, tn, tn + 0.5 * dt, time_steps[I] );
        }

        /**
         * decremental call of each method of each subproblem
         * @tparam I solving step
         * @param f          \ref problem to solve
         * @param tn         current time \f$t^n\f$
         * @param[in,out] ui \f$\texttt{ui}=\phi_{^{\Delta t}/_2}^{[f_1]}\circ\cdots\circ\phi_{\Delta
         * t}^{[f_{n}]}\circ\cdots\circ\phi_{^{\Delta t}/_2}^{[f_{i+1}]}(t^n,u^n)\f$
         * @param dt         time step \f$\Delta t\f$
         * @details The parameter @p ui is update to \f$\phi_{^{\Delta t}/_2}^{[f_i]}(t^n,\texttt{ui})\f$
         */
        template <typename value_t, typename... Methods_t>
        template <std::size_t I, typename Problem_t, typename state_t>
            requires( I > 0 )
        inline void strang<value_t, Methods_t...>::_call_dec( Problem_t& f, value_t tn, state_t& ui, value_t dt )
        {
            ui = detail::_split_solve<I>( f, methods, ui, tn, tn + 0.5 * dt, time_steps[I] );
            _call_dec<I - 1>( f, tn, ui, dt );
        }

        /**
         * call operator to initiate Strang splitting recursion
         * @param f  \ref problem to solve
         * @param tn current time \f$t^n\f$
         * @param un current solution \f$u^n \approx u(t^n)\f$
         * @param dt time step \f$\Delta t\f$
         */
        template <typename value_t, typename... Methods_t>
        template <typename Problem_t, typename state_t>
        auto
        strang<value_t, Methods_t...>::operator()( Problem_t& f, value_t tn, state_t const& un, value_t dt )
        {
            state_t ui = un;
            _call_inc( f, tn, ui, dt );
            return std::make_tuple( tn + dt, ui, dt );
        }

        /**
         * a helper factory for \ref strang functor from a tuple of methods
         * @param t tuple of \ref method
         * @return a \ref strang object build from the tuple of methods
         */
        template <typename value_t, typename... Methods_t>
        auto
        make_strang_from_tuple( std::tuple<Methods_t...> const& meths, std::array<value_t, sizeof...( Methods_t )> const& dts )
        {
            return strang<value_t, Methods_t...>( meths, dts );
        }

        /** @class strang_tuple
         *  a helper to deduce method for ::ode::make_method(splitting::strang_tuple<Algorithms_t...> const &, state_t const &)
         *  @tparam Algorithms_t variadic template of algorithms to solve each subproblem
         *  @details This is a dummy class to select correct \ref method to solve the problem
         */
        template <typename value_t, typename... Algorithms_t>
        struct strang_tuple
        {
            static constexpr std::size_t order        = 2;
            static constexpr bool is_splitting_method = true;

            std::tuple<Algorithms_t...> algos;
            std::array<value_t, sizeof...( Algorithms_t )> time_steps;

            strang_tuple( std::tuple<Algorithms_t...>&& algs, std::array<value_t, sizeof...( Algorithms_t )>&& dts );
        };

        /**
         * constructor of \ref strang_tuple from a variadic number of algorithms
         */
        template <typename value_t, typename... Algorithms_t>
        inline strang_tuple<value_t, Algorithms_t...>::strang_tuple( std::tuple<Algorithms_t...>&& algs,
            std::array<value_t, sizeof...( Algorithms_t )>&& dts )
            : algos( algs )
            , time_steps( dts )
        {
        }

        /**
         * a helper factory for \ref strang_tuple from a tuple of algorithms
         * @param a tuple of \ref method
         * @return a \ref strang_tuple object build from the tuple of methods
         */
        template <typename value_t, typename... Algorithms_t>
        auto
        make_strang_tuple( std::pair<Algorithms_t, value_t>&&... args )
        {
            return strang_tuple<value_t, Algorithms_t...>( std::forward_as_tuple( ( args.first )... ), { args.second... } );
        }

        template <typename T>
        concept is_splitting_method = requires( T t ) { T::is_splitting_method == true; };

    } // namespace splitting
} // namespace ode
