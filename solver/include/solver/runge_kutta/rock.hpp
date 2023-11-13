// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <algorithm>
#include <cmath>
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

    namespace detail
    {
        template <typename state_t>
        auto
        norm_2( state_t const& u )
        {
            return std::abs( u );
        }

        template <typename state_t>
            requires std::ranges::range<state_t>
        auto
        norm_2( state_t const& u )
        {
            return std::sqrt( std::accumulate( std::ranges::begin( u ),
                std::ranges::end( u ),
                0.,
                []( auto sum, auto x )
                {
                    return sum + std::abs( x ) * std::abs( x );
                } ) );
        }
    } // namespace detail

    /** @class rock2
     *  @brief define ROCK2 with user defined number of stages
     *
     *  @tparam N_stages_ number of stages
     *  @tparam value_t type of coefficients
     */
    template <typename value_t = double>
    struct rock2
    {
        // static_assert( N_stages_ > 1, "Number of stages should be at least 2 in ROCK2 method" );
        static constexpr bool is_embedded     = false;
        static constexpr std::size_t N_stages = stages::dynamic; // N_stages_;

        value_t ci1;
        value_t ci2;
        value_t ci3;
        value_t temp1;
        value_t temp2;
        value_t temp3;

        rock2()
            : ci1( 0. )
            , ci2( 0. )
            , ci3( 0. )
            , temp1( 0. )
            , temp2( 0. )
            , temp3( 0. )
        {
        }

        inline value_t
        recf( std::size_t i ) const
        {
            return rock2_coeff<value_t>::recf[i - 1];
        }

        inline value_t
        fp1( std::size_t i ) const
        {
            return rock2_coeff<value_t>::fp1[i - 1];
        }

        inline value_t
        fp2( std::size_t i ) const
        {
            return rock2_coeff<value_t>::fp2[i - 1];
        }

        inline std::size_t
        ms( std::size_t i ) const
        {
            return rock2_coeff<value_t>::ms[i - 1];
        }

        template <typename problem_t, typename state_t>
        value_t
        rocktrho( problem_t& f, value_t tn, state_t const& un )
        {
            value_t eigmax  = 0.;
            value_t eigmaxo = 0.;

            auto fn = f( tn, un );
            auto fz = fn;

            auto z = f( tn, fz );

            value_t ynor = detail::norm_2( un );
            value_t znor = detail::norm_2( z );

            value_t quot  = 0.;
            value_t dzyn  = 0.;
            value_t dfzfn = 0.;

            value_t const sqrt_eps = std::sqrt( std::numeric_limits<value_t>::epsilon() );

            if ( ynor != 0.0 && znor != 0.0 )
            {
                dzyn = ynor * sqrt_eps;
                quot = dzyn / znor;
                z    = un + quot * z;
            }
            else if ( ynor != 0.0 )
            {
                dzyn = ynor * sqrt_eps;
                z    = ( 1.0 + sqrt_eps ) * un;
            }
            else if ( znor != 0.0 )
            {
                dzyn = sqrt_eps;
                quot = dzyn / znor;
                z    = z * quot;
            }
            else
            {
                dzyn = sqrt_eps;
                z    = dzyn; // TODO: change this line to iterate over z if needed (same trick as norm_2 I think)
            }

            bool necessary                        = true;
            std::size_t iter                      = 0;
            static constexpr std::size_t max_iter = 50;
            static constexpr value_t safe         = 1.2;
            // start power method
            while ( necessary )
            {
                eigmaxo = eigmax;
                fz      = f( tn, z );
                dfzfn   = detail::norm_2( static_cast<state_t>( fz - fn ) );
                eigmax  = safe * dfzfn / dzyn;

                if ( dfzfn != 0.0 )
                {
                    quot = dzyn / dfzfn;
                    z    = un + quot * ( fz - fn );
                }
                else
                {
                    // TODO: make this only on one index (ind = iter % z.size())
                    z = un - ( z - un );
                }

                necessary = ( iter < max_iter ) && !( iter > 2 && std::abs( eigmax - eigmaxo ) <= 0.05 * eigmax );
                ++iter;
            }

            return eigmax;
        }

        std::pair<std::size_t, std::size_t>
        optimal_degree( std::size_t& mdeg )
        {
            std::size_t mz = 1, mr = 1;

            std::size_t iter = 1;
            while ( iter < rock2_coeff<value_t>::ms.size() + 1 && ms( iter ) / mdeg <= 1 )
            {
                mr = mr + ms( iter ) * 2 - 1;
                ++iter;
            }
            mdeg = ms( iter );
            mz   = iter;

            // for ( std::size_t iter = 1; iter < rock2_coeff<value_t>::ms.size() + 1; ++iter )
            // {
            //     if ( ms( iter ) / mdeg > 1 )
            //     {
            //         mdeg = ms( iter );
            //         mz   = iter;
            //         break;
            //     }
            //     mr = mr + ms( iter ) * 2 - 1;
            // }

            return { mz, mr };
        }

        template <typename problem_t, typename state_t>
        std::size_t
        compute_n_stages( problem_t& f, value_t tn, state_t const& un, value_t& dt )
        {
            double eigmax    = rocktrho( f, tn, un );
            std::size_t mdeg = static_cast<std::size_t>( std::ceil( std::sqrt( ( 1.5 + dt * eigmax ) / 0.811 ) ) );
            if ( mdeg > 200 )
            {
                mdeg = 200;
                dt   = 0.8 * ( static_cast<double>( mdeg * mdeg ) * 0.811 - 1.5 ) / eigmax;
            }
            mdeg = std::max( mdeg, static_cast<std::size_t>( 3 ) ) - 2;
            return mdeg;
        }

        template <typename problem_t, typename state_t>
        std::tuple<std::size_t, std::size_t, std::size_t>
        compute_n_stages_optimal_degree( problem_t& f, value_t tn, state_t const& un, value_t& dt )
        {
            std::size_t mdeg = compute_n_stages( f, tn, un, dt );
            auto [mz, mr]    = optimal_degree( mdeg );

            return { mdeg, mz, mr };
        }

        template <typename problem_t, typename state_t, typename array_ki_t>
        inline std::tuple<value_t, state_t, value_t>
        operator()( problem_t& f, value_t& tn, state_t const& un, array_ki_t& G, value_t& dt )
        {
            // std::size_t mdeg = compute_n_stages( f, tn, un, dt );
            // auto [mz, mr]    = optimal_degree( mdeg );
            auto [mdeg, mz, mr] = compute_n_stages_optimal_degree( f, tn, un, dt );
            // mdeg                = 4;
            G.resize( mdeg + 3 );

            temp1 = dt * recf( mr );
            ci1   = tn + temp1;
            ci2   = tn + temp1;
            ci3   = tn;

            G[0] = un + temp1 * f( tn, un );

            // std::cout << "rock2::op() " << mdeg << " " << tn << " " << dt << " " << mr << " . " << mz << "\n";

            if ( mdeg >= 2 )
            {
                temp1 = dt * recf( mr + 1 );
                temp3 = -recf( mr + 2 );
                temp2 = 1.0 - temp3;

                G[1] = temp1 * f( ci1, G[0] ) + temp2 * G[0] + temp3 * un;
                ci1  = temp1 + temp2 * ci2 + temp3 * ci3;
                ci3  = ci2;
                ci2  = ci1;
                for ( std::size_t j = 2; j < mdeg + 1; ++j )
                {
                    temp1 = dt * recf( mr + 2 * ( j - 2 ) + 1 );
                    temp3 = -recf( mr + 2 * ( j - 2 ) + 2 );
                    temp2 = 1.0 - temp3;

                    G[j] = temp1 * f( ci1, G[j - 1] ) + temp2 * G[j - 1] + temp3 * G[j - 2];
                    ci1  = temp1 + temp2 * ci2 + temp3 * ci3;
                    ci3  = ci2;
                    ci2  = ci1;
                }

                temp1 = dt * fp1( mz );
                temp2 = dt * fp2( mz );

                G[mdeg + 1] = G[mdeg] + temp1 * f( ci1, G[mdeg] );

                ci1 = ci1 + temp1;

                // G[mdeg + 2] = G[mdeg + 1] + temp1 * f( ci1, G[mdeg + 1] ) + temp2 * ( f( ci1, G[mdeg + 1] ) - f( ci1 - temp1, G[mdeg] )
                // );

                G[mdeg + 2] = G[mdeg + 1] + ( temp1 + temp2 ) * f( ci1, G[mdeg + 1] ) - temp2 / temp1 * ( G[mdeg + 1] - G[mdeg] );

                return { tn + dt, G[mdeg + 2], dt };
            }
            else
            {
                G.resize( 1 );

                return { tn + dt, G[0], dt };
            }

            // tn += dt;
        }

        // template <typename problem_t, typename state_t, typename array_ki_t>
        // inline state_t
        // stage( std::size_t j, problem_t& f, value_t tn, state_t const& un, array_ki_t const& G, value_t dt )
        // {
        //     mu    = recf( mr + ( j + 1 - 2 ) * 2 + 1 );
        //     kappa = recf( mr + ( j + 1 - 2 ) * 2 + 2 );
        //     nu    = -1 - kappa;

        //     return ( dt * mu ) * f( tn, G[j - 1] ) - nu * G[j - 1] - kappa * G[j - 2];
        // }

        // template <typename problem_t, typename state_t, typename array_ki_t>
        // inline state_t
        // first_stage_0( problem_t& f, value_t tn, state_t const& un, array_ki_t const& G, value_t dt )
        // {
        //     return un + ( dt * recf( mr ) ) * f( tn, un );
        // }

        // template <typename problem_t, typename state_t, typename array_ki_t>
        // inline state_t
        // first_stage_1( problem_t& f, value_t tn, state_t const& un, array_ki_t const& G, value_t dt )
        // {
        //     constexpr std::size_t j = 1;

        //     mu    = recf( mr + ( j + 1 - 2 ) * 2 + 1 );
        //     kappa = recf( mr + ( j + 1 - 2 ) * 2 + 2 );
        //     nu    = -1 - kappa;

        //     return ( dt * mu ) * f( tn, G[j - 1] ) - nu * G[j - 1] - kappa * un;
        // }

        // template <typename problem_t, typename state_t, typename array_ki_t>
        // inline state_t
        // last_stage_0( std::size_t j, problem_t& f, value_t tn, state_t const& un, array_ki_t const& G, value_t dt )
        // {
        //     double delta_t1 = dt * fp1( mz );
        //     double delta_t2 = dt * fp2( mz );

        //     return G[j - 1] + delta_t1 * f( tn, G[j - 1] );
        // }

        // template <typename problem_t, typename state_t, typename array_ki_t>
        // inline state_t
        // last_stage_1( std::size_t j, problem_t& f, value_t tn, state_t const& un, array_ki_t const& G, value_t dt )
        // {
        //     double delta_t1 = dt * fp1( mz );
        //     double delta_t2 = dt * fp2( mz );
        //     // u = -delta_t2*f(tn, G[j-2])
        //     // k = f(tn, G[j-1])

        //     return -delta_t2 * f( tn, G[j - 2] ) + G[j - 1] + ( delta_t1 + delta_t2 ) * f( tn, G[j - 1] );
        // }

        // template <typename problem_t, typename state_t, typename array_ki_t, std::size_t j>
        // inline state_t
        // stage( Stage<j>, problem_t& f, value_t, state_t const&, array_ki_t const& G, value_t dt )
        // {
        //     auto old_ci1 = ci1;

        //     ci1 = temp1 + temp2 * ci2 + temp3 * ci3;
        //     ci2 = old_ci1;
        //     ci3 = ci2;

        //     temp1 = dt * recf( mr + 2 * ( j - 2 ) + 1 );
        //     temp3 = -recf( mr + 2 * ( j - 2 ) + 2 );
        //     temp2 = 1.0 - temp3;

        //     return temp1 * f( ci1, G[j - 1] ) + temp2 * G[j - 1] + temp3 * G[j - 2];
        // }

        // template <typename problem_t, typename state_t, typename array_ki_t>
        // inline state_t
        // stage( Stage<0>, problem_t& f, value_t tn, state_t const& un, array_ki_t const&, value_t dt )
        // {
        //     // $$ g_1 = un + dt \mu_1 f(un) $$
        //     temp1 = dt * recf( mr );
        //     ci1   = tn + temp1;
        //     ci2   = tn + temp1;
        //     ci3   = tn;

        //     return un + temp1 * f( tn, un ); // be careful G[0] doesn't store un!!!
        // }

        // template <typename problem_t, typename state_t, typename array_ki_t>
        // inline state_t
        // stage( Stage<1>, problem_t& f, value_t, state_t const& un, array_ki_t const& G, value_t dt )
        // {
        //     // $$ g_2 = dt \mu_2 f(g_1) - \nu_2 g_1 - \kappa_2 un $$
        //     constexpr std::size_t j = 1;

        //     temp1 = dt * recf( mr + 2 * ( j - 2 ) + 1 );
        //     temp3 = -recf( mr + 2 * ( j - 2 ) + 2 );
        //     temp2 = 1.0 - temp3;

        //     return temp1 * f( ci1, G[j - 1] ) + temp2 * G[j - 1] + temp3 * un;
        // }

        // template <typename problem_t, typename state_t, typename array_ki_t>
        // inline state_t
        // stage( Stage<N_stages - 1>, problem_t& f, value_t, state_t const&, array_ki_t const& G, value_t dt )
        // {
        //     constexpr std::size_t j = N_stages - 1;

        //     temp1 = dt * fp1( mz );
        //     temp2 = dt * fp2( mz );

        //     return G[j - 1] + temp1 * f( ci1, G[j - 2] );
        // }

        // template <typename problem_t, typename state_t, typename array_ki_t>
        // inline state_t
        // stage( Stage<N_stages>, problem_t& f, value_t, state_t const&, array_ki_t const& G, value_t )
        // {
        //     // $$ g_s = g_s^\star - dt \sigma (1 - \tau / \sigma^2) ( f(g_{s-1}) - f(g_{s-2}) ) $$
        //     //
        //     // we use previous stages to compute $f(g_{s-1}) - f(g_{s-2})$ without evaluation of $f$ function.
        //     constexpr std::size_t j = N_stages;

        //     // ci1   = ci1 + temp1;
        //     // temp3 = temp2 / temp1;
        //     // return ( 1.0 - temp3 ) * G[j - 1] - temp3 * G[j - 2] + ( temp1 + temp2 ) * f( ci1, G[j - 1] );

        //     auto new_ci1 = ci1 + temp1;
        //     // return G[j - 1] + temp1 * f( ci1, G[j - 1] ) + temp2 * ( f( new_ci1, G[j - 1] ) - f( ci1, G[j - 2] ) );
        //     return G[j - 1] + temp2 * ( f( new_ci1, G[j - 2] ) - f( ci1, G[j - 3] ) );
        // }
    };

} // namespace ponio::runge_kutta::rock
