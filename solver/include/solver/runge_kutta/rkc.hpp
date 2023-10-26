// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <concepts>
#include <numeric>
#include <ranges>
#include <string_view>
#include <type_traits>

#include "../butcher_tableau.hpp"
#include "../detail.hpp"
#include "../linear_algebra.hpp"
#include "../ponio_config.hpp"
#include "../stage.hpp"

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
    template <std::size_t N_stages_, typename value_t = double>
    struct explicit_rkc2
    {
        static_assert( N_stages_ > 1, "Number of stages should be at least 2 in eRKC2" );
        static constexpr bool is_embedded     = false;
        static constexpr std::size_t N_stages = N_stages_;
        value_t w0;
        value_t w1;

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
                return ddT<J>( x ) / ( ::detail::power<2>( dT<J>( x ) ) );
            }
        }

        explicit_rkc2( value_t eps = 2. / 13. )
            : w0( 1. + eps / ( N_stages * N_stages ) )
            , w1( dT<N_stages>( w0 ) / ddT<N_stages>( w0 ) )
        {
        }

        template <typename problem_t, typename state_t, typename array_ki_t, std::size_t j>
        inline state_t
        stage( Stage<j>, problem_t& f, value_t tn, state_t& un, array_ki_t const& Yi, value_t dt )
        {
            value_t mj   = 2. * b<j>( w0 ) / b<j - 1>( w0 ) * w0;
            value_t nj   = -b<j>( w0 ) / b<j - 2>( w0 );
            value_t mjt  = 2. * b<j>( w0 ) / b<j - 1>( w0 ) * w1;
            value_t gjt  = -( 1. - b<j - 1>( w0 ) * T<j - 1>( w0 ) ) * mjt;
            value_t cjm1 = dT<N_stages>( w0 ) / ddT<N_stages>( w0 ) * ddT<j - 1>( w0 ) / dT<j - 1>( w0 );

            return ( 1. - mj - nj ) * un + mj * Yi[j - 1] + nj * Yi[j - 2] + mjt * dt * f( tn + cjm1 * dt, Yi[j - 1] ) + gjt * dt * Yi[0];
        }

        template <typename problem_t, typename state_t, typename array_ki_t>
        inline state_t
        stage( Stage<0>, problem_t& f, value_t tn, state_t& un, array_ki_t const&, value_t )
        {
            return f( tn, un ); // be careful Yi[0] stores f(tn,un) not un!!!
        }

        template <typename problem_t, typename state_t, typename array_ki_t>
        inline state_t
        stage( Stage<1>, problem_t&, value_t, state_t& un, array_ki_t const& Yi, value_t dt )
        {
            value_t m1t = b<1>( w0 ) * w1;
            return un + dt * m1t * Yi[0];
        }

        template <typename problem_t, typename state_t, typename array_ki_t>
        inline state_t
        stage( Stage<2>, problem_t& f, value_t tn, state_t& un, array_ki_t const& Yi, value_t dt )
        {
            value_t m2  = 2. * w0;
            value_t n2  = -1.;
            value_t m2t = 2. * w1;
            value_t c2  = dT<N_stages>( w0 ) / ddT<N_stages>( w0 ) * ddT<2>( w0 ) / dT<2>( w0 );
            value_t c1  = c2 / dT<2>( w0 );
            value_t g2t = -( 1. - b<1>( w0 ) * T<1>( w0 ) ) * m2t;

            return ( 1. - m2 - n2 ) * un + m2 * Yi[1] + n2 * un + m2t * dt * f( tn + c1 * dt, Yi[1] ) + g2t * dt * Yi[0];
        }
    };

} // namespace ponio::runge_kutta::chebyshev
