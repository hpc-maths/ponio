// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private, include "../runge_kutta.hpp"

#pragma once

#include <cstddef>
#include <string_view> // NOLINT(misc-include-cleaner)

#include "../detail.hpp" // NOLINT(misc-include-cleaner)
#include "../iteration_info.hpp"
#include "../stage.hpp" // NOLINT(misc-include-cleaner)

namespace ponio::runge_kutta::chebyshev
{

    /**
     * @brief Chebyshev polynomial of first kind by recursive method: \f$$\begin{aligned}T_0(x) &= 1\\ T_1(x) &= x\\ T_{n+1}(x) &=
     * 2xT_n(x) - T_{n-1}(x)\end{aligned}}\f$$
     *
     * @tparam N degree of polynomial
     * @tparam value_t type of \f$x\f$
     * @param x value where evaluate \f$T_n\f$
     */
    template <std::size_t N, typename value_t>
    value_t
    T( value_t const x )
    {
        if constexpr ( N == 0 )
        {
            return static_cast<value_t>( 1. );
        }
        else if constexpr ( N == 1 )
        {
            return x;
        }
        else
        {
            return 2 * x * T<N - 1>( x ) - T<N - 2>( x );
        }
    }

    /**
     * @brief Chebyshev polynomial of second kind by recursive method: \f$$\begin{aligned}U_0(x) &= 1\\ U_1(x) &= 2x\\ U_{n+1}(x) &=
     * 2xU_n(x) - U_{n-1}(x)\end{aligned}}\f$$
     *
     * @tparam N degree of polynomial
     * @tparam value_t type of \f$x\f$
     * @param x value where evaluate \f$U_n\f$
     */
    template <std::size_t N, typename value_t>
    value_t
    U( value_t const x )
    {
        if constexpr ( N == 0 )
        {
            return static_cast<value_t>( 1. );
        }
        else if constexpr ( N == 1 )
        {
            return 2 * x;
        }
        else
        {
            return 2 * x * U<N - 1>( x ) - U<N - 2>( x );
        }
    }

    /**
     * @brief derivatives of Chebyshev polynomial: \f$\frac{\mathrm{d}T_n}{\mathrm{d}x}(x)=nU_{n-1}(x)\f$
     *
     * @tparam N degree of polynomial
     * @tparam value_t type of \f$x\f$
     * @param x value where evaluate \f$\frac{\mathrm{d}T_n}{\mathrm{d}x}\f$
     */
    template <std::size_t N, typename value_t>
    value_t
    dT( value_t const x )
    {
        if constexpr ( N == 0 )
        {
            return static_cast<value_t>( 0. );
        }
        else
        {
            return N * U<N - 1>( x );
        }
    }

    /**
     * @brief second derivative of Chebyshev polynomial: \f$\frac{\mathrm{d}^2T_n}{\mathrm{d}x^2}(x)=n\frac{nT_n(x) - xU_{n-1}(x)}{x^2
     * -1}\f$
     *
     * @tparam N degree of polynomial
     * @tparam value_t type of \f$x\f$
     * @param x value where evaluate \f$\frac{\mathrm{d}^2T_n}{\mathrm{d}x^2}\f$
     */
    template <std::size_t N, typename value_t>
    value_t
    ddT( value_t const x )
    {
        if constexpr ( N == 0 )
        {
            return static_cast<value_t>( 0. );
        }
        else
        {
            return N * ( N * T<N>( x ) - x * U<N - 1>( x ) ) / ( x * x - static_cast<value_t>( 1. ) );
        }
    }

    /** @class explicit_rkc2
     *  @brief define RKC2 with user defined number of stages
     *
     *  @tparam N_stages_ number of stages
     *  @tparam value_t type of coefficients
     */
    template <std::size_t N_stages_, typename _value_t = double>
    struct explicit_rkc2
    {
        static_assert( N_stages_ > 1, "Number of stages should be at least 2 in eRKC2" );
        static constexpr std::size_t N_stages = N_stages_;
        static constexpr std::size_t order    = 2;
        static constexpr std::string_view id  = "RKC2";
        static constexpr bool is_embedded     = false;
        using value_t                         = _value_t;

        value_t w0;
        value_t w1;
        iteration_info<explicit_rkc2> _info;

        /**
         * @brief computes \f$b_j = \f$ coefficients with following formula: \f$b_0=b_2\f$, \f$b_1=\frac{1}{\omega_0}\f$, \f$b_j =
         * \frac{T_j''(\omega_0)}{(T_j'(\omega_0))^2}\f$
         *
         * @tparam J index \f$j\f$
         * @param x  value of \f$\omega_0\f$
         */
        template <std::size_t J>
        static constexpr value_t
        b( value_t const x )
        {
            if constexpr ( J == 0 || J == 1 )
            {
                return b<2>( x );
            }
            else
            {
                return ddT<J>( x ) / ( detail::power<2>( dT<J>( x ) ) );
            }
        }

        /**
         * @brief Construct a new explicit RKC2 algorithm
         *
         * @param eps value of relaxation
         */
        explicit_rkc2( value_t eps = 2. / 13. )
            : w0( 1. + eps / ( N_stages * N_stages ) )
            , w1( dT<N_stages>( w0 ) / ddT<N_stages>( w0 ) )
            , _info()
        {
            _info.number_of_eval = N_stages;
        }

        /**
         * @brief computes stage \f$j\f$ of RKC2 method
         *
         * @tparam problem_t  type of operator \f$f\f$
         * @tparam state_t    type of current state
         * @tparam array_ki_t type of array of temporary stages
         * @tparam j          integer of stage
         * @param f  operator \f$f\f$
         * @param tn current time
         * @param yn current state
         * @param Yj array of temporary stages
         * @param dt current time step
         * @param yi computed output step
         *
         * @details \f$y_j = (1 - \mu_j - \nu_j)y^n + \mu_j y_{j-1} + \nu_j y_{j-2} + \tilde{\mu}_j\Delta tf(t^n + c_{j-1}\Delta t, y_{j-1})
         * + \tilde{\gamma}_j\Delta t f(t^n, y^n)\f$
         */
        template <typename problem_t, typename state_t, typename array_ki_t, std::size_t j>
        void
        stage( Stage<j>, problem_t& f, value_t tn, state_t& yn, array_ki_t const& Yj, value_t dt, state_t& yi )
        {
            value_t const mj   = 2. * b<j>( w0 ) / b<j - 1>( w0 ) * w0;
            value_t const nj   = -b<j>( w0 ) / b<j - 2>( w0 );
            value_t const mjt  = 2. * b<j>( w0 ) / b<j - 1>( w0 ) * w1;
            value_t const gjt  = -( 1. - b<j - 1>( w0 ) * T<j - 1>( w0 ) ) * mjt;
            value_t const cjm1 = dT<N_stages>( w0 ) / ddT<N_stages>( w0 ) * ddT<j - 1>( w0 ) / dT<j - 1>( w0 );

            f( tn + cjm1 * dt, Yj[j - 1], yi ); // first compute yi = f(t^n + c_{j-1}\Delta t, y_{j-1})
            yi = ( 1. - mj - nj ) * yn + mj * Yj[j - 1] + nj * Yj[j - 2] + mjt * dt * yi + gjt * dt * Yj[0];
        }

        /**
         * @brief computes pseudo first stage of RKC2 method
         *
         * @tparam problem_t  type of operator \f$f\f$
         * @tparam state_t    type of current state
         * @tparam array_ki_t type of array of temporary stages
         * @tparam j          integer of stage
         * @param f  operator \f$f\f$
         * @param tn current time
         * @param yn current state
         * @param yi computed output step
         *
         * @details \f$y_0 = f(t^n, y^n)\f$
         */
        template <typename problem_t, typename state_t, typename array_ki_t>
        void
        stage( Stage<0>, problem_t& f, value_t tn, state_t& yn, array_ki_t const&, value_t, state_t& yi )
        {
            f( tn, yn, yi ); // be careful Yj[0] stores f(tn,yn) not yn!!!
        }

        /**
         * @brief computes first stage of RKC2 method
         *
         * @tparam problem_t  type of operator \f$f\f$
         * @tparam state_t    type of current state
         * @tparam array_ki_t type of array of temporary stages
         * @tparam j          integer of stage
         * @param yn current state
         * @param Yj array of temporary stages
         * @param dt current time step
         * @param yi computed output step
         *
         * @details \f$y_1 = y^n + \tilde{\mu}_1 \Delta t f(t^n, y^n)\f$
         */
        template <typename problem_t, typename state_t, typename array_ki_t>
        void
        stage( Stage<1>, problem_t&, value_t, state_t& yn, array_ki_t const& Yj, value_t dt, state_t& yi )
        {
            value_t const m1t = b<1>( w0 ) * w1;
            yi                = yn + dt * m1t * Yj[0];
        }

        /**
         * @brief computes second stage of RKC2 method
         *
         * @tparam problem_t  type of operator \f$f\f$
         * @tparam state_t    type of current state
         * @tparam array_ki_t type of array of temporary stages
         * @tparam j          integer of stage
         * @param f  operator \f$f\f$
         * @param tn current time
         * @param yn current state
         * @param Yj array of temporary stages
         * @param dt current time step
         * @param yi computed output step
         *
         * @details \f$y_2 = (1-\mu_2 -\nu_2)y^n + \mu_2y_1 + \nu_2y^n + \tilde{\mu}_2\Delta t f(t^n + c_1\Delta t, y_1) + \tilde{\gamma}_2
         * \Delta t f(t^n, y^n)\f$
         */
        template <typename problem_t, typename state_t, typename array_ki_t>
        void
        stage( Stage<2>, problem_t& f, value_t tn, state_t& yn, array_ki_t const& Yj, value_t dt, state_t& yi )
        {
            value_t const m2  = 2. * w0;
            value_t const n2  = -1.;
            value_t const m2t = 2. * w1;
            value_t const c2  = dT<N_stages>( w0 ) / ddT<N_stages>( w0 ) * ddT<2>( w0 ) / dT<2>( w0 );
            value_t const c1  = c2 / dT<2>( w0 );
            value_t const g2t = -( 1. - b<1>( w0 ) * T<1>( w0 ) ) * m2t;

            f( tn + c1 * dt, Yj[1], yi ); // first compute yi = f(t^n + c_1\Delta t, y_1)
            yi = ( 1. - m2 - n2 ) * yn + m2 * Yj[1] + n2 * yn + m2t * dt * yi + g2t * dt * Yj[0];
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
    };

} // namespace ponio::runge_kutta::chebyshev
