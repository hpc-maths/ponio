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

#include "rock_coeff.hpp"

namespace ponio::runge_kutta::rock
{

    /** @class rock2
     *  @brief define ROCK2 with user defined number of stages
     *
     *  @tparam N_stages_ number of stages
     *  @tparam value_t type of coefficients
     */
    template <std::size_t N_stages_, typename value_t = double>
    struct rock2
    {
        static_assert( N_stages_ > 1, "Number of stages should be at least 2 in ROCK2 method" );
        static constexpr bool is_embedded     = false;
        static constexpr std::size_t N_stages = N_stages_;

        std::size_t mz;
        std::size_t mr;
        value_t ci1;
        value_t ci2;
        value_t ci3;
        value_t temp1;
        value_t temp2;
        value_t temp3;

        rock2()
            : mz( 1 )
            , mr( 1 )
            , ci1( 0. )
            , ci2( 0. )
            , ci3( 0. )
            , temp1( 0. )
            , temp2( 0. )
            , temp3( 0. )
        {
        }

        inline value_t
        recf( std::size_t i )
        {
            return rock2_coeff<value_t>::recf[i - 1];
        }

        inline value_t
        fp1( std::size_t i )
        {
            return rock2_coeff<value_t>::fp1[i - 1];
        }

        inline value_t
        fp2( std::size_t i )
        {
            return rock2_coeff<value_t>::fp2[i - 1];
        }

        template <typename problem_t, typename state_t, typename array_ki_t, std::size_t j>
        inline state_t
        stage( Stage<j>, problem_t& f, value_t, state_t const&, array_ki_t const& G, value_t dt )
        {
            auto old_ci1 = ci1;

            ci1 = temp1 + temp2 * ci2 + temp3 * ci3;
            ci2 = old_ci1;
            ci3 = ci2;

            temp1 = dt * recf( mr + 2 * ( j - 2 ) + 1 );
            temp3 = -recf( mr + 2 * ( j - 2 ) + 2 );
            temp2 = 1.0 - temp3;

            return temp1 * f( ci1, G[j - 1] ) + temp2 * G[j - 1] + temp3 * G[j - 2];
        }

        template <typename problem_t, typename state_t, typename array_ki_t>
        inline state_t
        stage( Stage<0>, problem_t& f, value_t tn, state_t const& un, array_ki_t const&, value_t dt )
        {
            // $$ g_1 = un + dt \mu_1 f(un) $$
            temp1 = dt * recf( mr );
            ci1   = tn + temp1;
            ci2   = tn + temp1;
            ci3   = tn;

            return un + temp1 * f( tn, un ); // be careful G[0] doesn't store un!!!
        }

        template <typename problem_t, typename state_t, typename array_ki_t>
        inline state_t
        stage( Stage<1>, problem_t& f, value_t, state_t const& un, array_ki_t const& G, value_t dt )
        {
            // $$ g_2 = dt \mu_2 f(g_1) - \nu_2 g_1 - \kappa_2 un $$
            constexpr std::size_t j = 1;

            temp1 = dt * recf( mr + 2 * ( j - 2 ) + 1 );
            temp3 = -recf( mr + 2 * ( j - 2 ) + 2 );
            temp2 = 1.0 - temp3;

            return temp1 * f( ci1, G[j - 1] ) + temp2 * G[j - 1] + temp3 * un;
        }

        template <typename problem_t, typename state_t, typename array_ki_t>
        inline state_t
        stage( Stage<N_stages - 1>, problem_t& f, value_t, state_t const&, array_ki_t const& G, value_t dt )
        {
            constexpr std::size_t j = N_stages - 1;

            temp1 = dt * fp1( mz );
            temp2 = dt * fp2( mz );

            return G[j - 1] + temp1 * f( ci1, G[j - 2] );
        }

        template <typename problem_t, typename state_t, typename array_ki_t>
        inline state_t
        stage( Stage<N_stages>, problem_t& f, value_t, state_t const&, array_ki_t const& G, value_t )
        {
            // $$ g_s = g_s^\star - dt \sigma (1 - \tau / \sigma^2) ( f(g_{s-1}) - f(g_{s-2}) ) $$
            //
            // we use previous stages to compute $f(g_{s-1}) - f(g_{s-2})$ without evaluation of $f$ function.
            constexpr std::size_t j = N_stages;

            ci1   = ci1 + temp1;
            temp3 = temp2 / temp1;
            return ( 1.0 - temp3 ) * G[j - 1] - temp3 * G[j - 2] + ( temp1 + temp2 ) * f( ci1, G[j - 1] );
        }
    };

} // namespace ponio::runge_kutta::rock
