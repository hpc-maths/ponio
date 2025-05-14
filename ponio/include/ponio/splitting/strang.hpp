// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private, include "../splitting.hpp"

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <string_view> // NOLINT(misc-include-cleaner)
#include <tuple>
#include <type_traits>
#include <utility>

#include "../detail.hpp" // NOLINT(misc-include-cleaner)
#include "../ponio_config.hpp"
#include "detail.hpp"

namespace ponio::splitting::strang
{

    // ---- class strang --------------------------------------------

    /** @class strang
     *  Strang splitting method
     *  @tparam methods_t list of methods to solve each sub-problem
     */
    template <typename _value_t, typename... methods_t>
    struct strang : detail::splitting_base<_value_t, methods_t...>
    {
        using value_t = _value_t;
        using base_t  = detail::splitting_base<value_t, methods_t...>;

        using base_t::splitting_base;

        using base_t::is_splitting_method;
        using base_t::N_methods;

        using base_t::methods;
        using base_t::time_steps;

        static constexpr std::size_t order   = 2;
        static constexpr std::string_view id = "strang";
        static constexpr std::size_t N_steps = 2 * N_methods - 1;

        iteration_info<strang> _info;

        strang( std::tuple<methods_t...> const& meths, std::array<value_t, N_methods> const& dts )
            : base_t( meths, dts )
            , _info( methods )
        {
        }

        // end of incremental recursion, start of decremental recursion
        // _call_inc can not be outside the class definition due to llvm bug
        // (see https://github.com/llvm/llvm-project/issues/56482)
        template <std::size_t I = 0, typename Problem_t, typename state_t>
            requires( I == N_methods - 1 )
        void _call_inc( Problem_t& f, value_t tn, state_t& ui, value_t dt )
        {
            ui = detail::_split_solve<I>( f, methods, ui, tn, tn + dt, time_steps[I], _info );
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
        template <std::size_t I = 0, typename Problem_t, typename state_t>
            requires( I < N_methods - 1 )
        void _call_inc( Problem_t& f, value_t tn, state_t& ui, value_t dt )
        {
            ui = detail::_split_solve<I>( f, methods, ui, tn, tn + 0.5 * dt, time_steps[I], _info );
            _call_inc<I + 1>( f, tn, ui, dt );
        }

        template <std::size_t I = N_methods - 1, typename Problem_t, typename state_t>
            requires( I == 0 )
        void _call_dec( Problem_t& f, value_t tn, state_t& ui, value_t dt );

        template <std::size_t I = N_methods - 1, typename Problem_t, typename state_t>
            requires( I > 0 )
        void _call_dec( Problem_t& f, value_t tn, state_t& ui, value_t dt );

        template <typename Problem_t, typename state_t>
        auto
        operator()( Problem_t& f, value_t tn, state_t const& un, value_t dt );

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
        stages( std::integral_constant<std::size_t, I> )
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
        stages( std::integral_constant<std::size_t, I> ) const
        {
            return std::get<I>( methods ).stages();
        }
    };

    // end of decremental recursion, end of recursion
    template <typename value_t, typename... methods_t>
    template <std::size_t I, typename Problem_t, typename state_t>
        requires( I == 0 )
    inline void strang<value_t, methods_t...>::_call_dec( Problem_t& f, value_t tn, state_t& ui, value_t dt )
    {
        ui = detail::_split_solve<I>( f, methods, ui, tn + 0.5 * dt, tn + dt, time_steps[I], _info );
    }

    /**
     * decremental call of each method of each subproblem
     * @tparam I solving step
     * @param f          \ref problem to solve
     * @param tn         current time \f$t^n\f$
     * @param[in,out] ui current state of the substep I \f$\texttt{ui}=\phi_{^{\Delta t}/_2}^{[f_1]}\circ\cdots\circ\phi_{\Delta
     * t}^{[f_{n}]}\circ\cdots\circ\phi_{^{\Delta t}/_2}^{[f_{i+1}]}(t^n,u^n)\f$
     * @param dt         time step \f$\Delta t\f$
     * @details The parameter @p ui is update to \f$\phi_{^{\Delta t}/_2}^{[f_i]}(t^n,\texttt{ui})\f$
     */
    template <typename value_t, typename... methods_t>
    template <std::size_t I, typename Problem_t, typename state_t>
        requires( I > 0 )
    inline void strang<value_t, methods_t...>::_call_dec( Problem_t& f, value_t tn, state_t& ui, value_t dt )
    {
        ui = detail::_split_solve<I>( f, methods, ui, tn + 0.5 * dt, tn + dt, time_steps[I], _info );
        _call_dec<I - 1>( f, tn, ui, dt );
    }

    /**
     * call operator to initiate Strang splitting recursion
     * @param f  \ref problem to solve
     * @param tn current time \f$t^n\f$
     * @param un current solution \f$u^n \approx u(t^n)\f$
     * @param dt time step \f$\Delta t\f$
     */
    template <typename value_t, typename... methods_t>
    template <typename Problem_t, typename state_t>
    auto
    strang<value_t, methods_t...>::operator()( Problem_t& f, value_t tn, state_t const& un, value_t dt )
    {
        _info.reset_eval();

        state_t ui = un;
        _call_inc( f, tn, ui, dt );
        return std::make_tuple( tn + dt, ui, dt );
    }

    // ---- *helper* ----

    /**
     * a helper factory for @ref ponio::splitting::detail::splitting_tuple from a tuple of algorithms to build a Strang method
     *
     * @tparam value_t      type of coefficients
     * @tparam Algorithms_t variadic list of types of algorithms
     * @param args          variadic list of pairs of algorithm and time step
     * @return a @ref ponio::splitting::detail::splitting_tuple object build from the tuple of methods
     */
    template <typename value_t, typename... Algorithms_t>
    auto
    make_strang_tuple( std::pair<Algorithms_t, value_t>&&... args )
    {
        return detail::splitting_tuple<strang, value_t, void, Algorithms_t...>( std::forward_as_tuple( ( args.first )... ),
            { args.second... } );
    }

    // ---- class adaptive_strang -----------------------------------

    /** @class adaptive_strang
     *  adaptive time step Strang splitting method
     *  @tparam methods_t list of methods to solve each sub-problem
     */
    template <typename _value_t, typename... methods_t>
    struct adaptive_strang : detail::splitting_base<_value_t, methods_t...>
    {
        using value_t = _value_t;
        using base_t  = detail::splitting_base<value_t, methods_t...>;

        using base_t::splitting_base;

        using base_t::is_splitting_method;
        using base_t::N_methods;

        using base_t::methods;
        using base_t::time_steps;

        static constexpr std::size_t order   = 1;
        static constexpr std::string_view id = "adaptive_strang";
        static constexpr std::size_t N_steps = 2 * N_methods - 1;
        static constexpr bool is_embedded    = true;

        iteration_info<adaptive_strang> _info;

        adaptive_strang( std::tuple<methods_t...> const& meths,
            std::array<value_t, N_methods> const& dts,
            value_t delta,
            value_t tol = default_config::tol )
            : base_t( meths, dts )
            , _info( methods, delta, tol )
        {
        }

        template <std::size_t I = 0, typename Problem_t, typename state_t>
            requires( I == N_methods - 1 )
        void _call_inc( Problem_t& f, value_t tn, state_t& ui, value_t dt, value_t shift )
        {
            ui = detail::_split_solve<I>( f, methods, ui, tn, tn + dt, time_steps[I], _info );
            _call_dec<I - 1>( f, tn, ui, dt, shift );
        }

        template <std::size_t I = 0, typename Problem_t, typename state_t>
            requires( 0 < I && I < N_methods - 1 )
        void _call_inc( Problem_t& f, value_t tn, state_t& ui, value_t dt, value_t shift )
        {
            ui = detail::_split_solve<I>( f, methods, ui, tn, tn + 0.5 * dt, time_steps[I], _info );
            _call_inc<I + 1>( f, tn, ui, dt, shift );
        }

        template <std::size_t I = 0, typename Problem_t, typename state_t>
            requires( I == 0 )
        void _call_inc( Problem_t& f, value_t tn, state_t& ui, value_t dt, value_t shift )
        {
            ui = detail::_split_solve<I>( f, methods, ui, tn, tn + ( 0.5 + shift ) * dt, time_steps[I], _info );
            _call_inc<I + 1>( f, tn, ui, dt, shift );
        }

        template <std::size_t I = N_methods - 1, typename Problem_t, typename state_t>
            requires( I > 0 )
        void _call_dec( Problem_t& f, value_t tn, state_t& ui, value_t dt, value_t shift )
        {
            ui = detail::_split_solve<I>( f, methods, ui, tn + 0.5 * dt, tn + dt, time_steps[I], _info );
            _call_dec<I - 1>( f, tn, ui, dt, shift );
        }

        template <std::size_t I = N_methods - 1, typename Problem_t, typename state_t>
            requires( I == 0 )
        void _call_dec( Problem_t& f, value_t tn, state_t& ui, value_t dt, value_t shift )
        {
            ui = detail::_split_solve<I>( f, methods, ui, tn + ( 0.5 + shift ) * dt, tn + dt, time_steps[I], _info );
        }

        template <typename Problem_t, typename state_t>
        auto
        operator()( Problem_t& f, value_t tn, state_t& un, value_t dt )
        {
            _info.reset_eval();

            state_t u_np1_ref   = un;
            state_t u_np1_shift = un;

            // TODO: launch _call_inc in parallel
            // j'ai voulu lancer chacune de ces 2 fonctions dans des threads différents
            // avec std::thread
            // mais j'ai eu une erreur, donc à corriger plus tard
            _call_inc( f, tn, u_np1_ref, dt, 0. );
            _call_inc( f, tn, u_np1_shift, dt, _info.delta );

            _info.error = ::ponio::detail::error_estimate( un, u_np1_ref, u_np1_shift );

            value_t new_dt = 0.9 * std::sqrt( _info.tolerance / _info.error ) * dt;
            new_dt         = std::min( std::max( 0.2 * dt, new_dt ), 5. * dt );

            _info.success = _info.error < _info.tolerance;

            if ( !_info.success )
            {
                return std::make_tuple( tn, un, new_dt );
            }
            return std::make_tuple( tn + dt, u_np1_ref, new_dt );
        }

        /**
         * @brief estimate Lipschitz constant of problem \f$f\f$ at time \f$t^n\f$
         *
         * @param f  callable obect which represents the problem to solve
         * @param tn time \f$t^n\f$ last time where solution is computed
         * @param un computed solution \f$u^n\f$ à time \f$t^n\f$
         * @param dt time step
         *
         * The error of classical splitting method is:
         * \f[
         *      \mathcal{S}^{\Delta t}u_0 - \mathcal{T}^{\Delta t}u_0 = C_0 \Delta t^3
         * \f]
         *
         * where $\mathcal{S} is the time flow of Strang method, and $\mathcal{T}$ is the exact time flow of the equation, with a shifted
         * splitting method, we can compute
         *
         * \f[
         *      \mathcal{S}^{\Delta t}u_0 - \mathcal{S}_\delta^{\Delta t} u_0 = \delta C_\delta \Delta t^2
         * \f]
         *
         * we find the optimal time step \f$\Delta t^\star\f$ as
         *
         * \f[
         *      \| \mathcal{S}^{\Delta t}u_0 - \mathcal{T}^{\Delta t}u_0 \| \leq \| \mathcal{S}^{\Delta t}u_0 - \mathcal{S}_\delta^{\Delta
         * t}u_0 \|,\ \forall \Delta t \leq \Delta t^\star
         * \f]
         *
         * we get an estimate of \f$\Delta t^\star\f$
         *
         * \f[
         *      \Delta t^\star \approx \frac{\delta C_\delta}{C_0}
         * \f]
         *
         * To estimate \f$C_0\f$ we need compute two errors
         *
         * \f[
         *      \begin{aligned}
         *          e_1 &= \| \mathcal{S}^{a_1\Delta t}u_0 \mathcal{S}^{b_1\Delta t}\mathcal{S}^{c_1\Delta t}u_0 \| \\
         *          e_2 &= \| \mathcal{S}^{a_2\Delta t}u_0 \mathcal{S}^{b_2\Delta t}\mathcal{S}^{c_2\Delta t}u_0 \|
         *      \end{aligned}
         * \f]
         *
         * we obtain two inequalities
         *
         * \f[
         *      \begin{aligned}
         *          \|e_1 - (a_1^3 - b_1^3)C_0\Delta t^3\| \leq \omega C_0 c_1^3\Delta t^3 \\
         *          \|e_2 - (a_2^3 - b_2^3)C_0\Delta t^3\| \leq \omega C_0 c_2^3\Delta t^3
         *      \end{aligned}
         * \f]
         *
         * we choose parameters as
         *
         * \f[
         *      a_1 = 1,\quad b_1=\frac{1}{2},\quad c_1=\frac{1}{2},\quad a_2=c_1,\quad b_2=\frac{2}{5},\quad c_2=\frac{1}{10}
         * \f]
         *
         * to obtain estimate of \f$C_0\f$ (error estimate) and \f$\omega\f$ (Lipschitz constant estimate).
         *
         */
        template <typename Problem_t, typename state_t>
        auto
        lipschitz_constant_estimate( Problem_t& f, value_t tn, state_t& un, value_t dt )
        {
            static constexpr value_t a1 = 1.0;
            static constexpr value_t b1 = 0.5;
            static constexpr value_t c1 = 0.5;
            static constexpr value_t a2 = c1;
            static constexpr value_t b2 = 0.4;
            static constexpr value_t c2 = 0.1;

            auto local_error = [&]( value_t a, value_t b, value_t c )
            {
                state_t u_np1_a  = un;
                state_t u_np1_bc = un;

                // compute strang(tn, un, a*dt)
                _call_inc( f, tn, u_np1_a, a * dt, 0. );

                // compute strang(tn+c*dt, strang(tn, un, c*dt), b*dt )
                _call_inc( f, tn, u_np1_bc, c * dt, 0. );
                _call_inc( f, tn + c * dt, u_np1_bc, b * dt, 0. );

                return ::ponio::detail::norm( u_np1_a - u_np1_bc );
            };

            // e1 = || S_{a1*dt}(un) - S_{b1*dt}(S_{c1*dt}(un)) ||
            // e2 = || S_{a2*dt}(un) - S_{b2*dt}(S_{c2*dt}(un)) ||
            auto e1 = local_error( a1, b1, c1 );
            auto e2 = local_error( a2, b2, c2 );

            // we have:
            // || e1 - (a1**3 - b1**3)C0 dt**3 || <= w C0 c1**3 dt**3
            // || e2 - (a2**3 - b2**3)C0 dt**3 || <= w C0 c2**3 dt**3
            // so (if inequalities became equalities):
            // c2**6 * ( e1 - (a1**3 - b1**3)C0 dt**3 )**2 ~= c1**6 * ( e2 - (a2**3 - b2**3)C0 dt**3 )**2
            // then need to solve:
            // alpha C0**2 + beta C0 + gamma = 0
            // with:
            // alpha = c2**6 * dt**6 * (a1**3 - b1**3)**2 - c1**6 * dt**6 * (a2**3 - b2**3)**2
            // beta = -2 * ( c2**6 * dt**3 * (a1**3 - b1**3) * e1  + c1**6 * dt**3 * (a2**3 - b2**3) * e2 )
            // gamma = c2**6 * e2**2 - c1**6 * e2**2
            // we find positive root of this polynomial to compute C0

            using ::ponio::detail::power;

            value_t alpha = power<6>( c2 * dt ) * power<2>( power<3>( a1 ) - power<3>( b1 ) )
                          - power<6>( c1 * dt ) * power<2>( power<3>( a2 ) - power<3>( b2 ) );
            value_t beta = -2 * power<3>( dt )
                         * ( power<6>( c2 ) * ( power<3>( a1 ) - power<3>( b1 ) ) * e1
                             + power<6>( c1 ) * ( power<3>( a2 ) - power<3>( b2 ) ) * e2 );
            value_t gamma = power<6>( c2 ) * power<2>( e1 ) - power<6>( c1 ) * power<2>( e2 );

            value_t discriminent = 2. * beta * beta - alpha * gamma;
            value_t C0           = ( -beta - std::sqrt( discriminent ) ) / ( 2. * alpha );

            value_t omega = std::abs( e1 - ( power<3>( a1 ) - power<3>( b1 ) ) * C0 * power<3>( dt ) ) / ( C0 * power<3>( c1 * dt ) );

            return std::make_pair( omega, C0 );
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

        /**
         * @brief gets array of stages of the Ith method
         *
         * @tparam I index of step method
         */
        template <std::size_t I>
        auto&
        stages( std::integral_constant<std::size_t, I> )
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
        stages( std::integral_constant<std::size_t, I> ) const
        {
            return std::get<I>( methods ).stages();
        }
    };

    // ---- *helper* ----

    /**
     * a helper factory for @ref ponio::splitting::detail::splitting_tuple from a tuple of algorithms to build an adaptive time step Strang
     * method
     *
     * @tparam value_t      type of coefficients
     * @tparam Algorithms_t variadic list of types of algorithms
     * @param delta     shift argument
     * @param tolerance tolerance for adaptive time step algorithm
     * @param args      variadic list of pairs of algorithm and time step
     * @return a @ref ponio::splitting::detail::splitting_tuple object build from the tuple of methods
     */
    template <typename value_t, typename... Algorithms_t>
    auto
    make_adaptive_strang_tuple( value_t delta, value_t tolerance, std::pair<Algorithms_t, value_t>&&... args )
    {
        return detail::splitting_tuple<adaptive_strang, value_t, std::tuple<value_t, value_t>, Algorithms_t...>(
            std::forward_as_tuple( ( args.first )... ),
            { args.second... },
            std::make_tuple( delta, tolerance ) );
    }

} // namespace ponio::splitting::strang
