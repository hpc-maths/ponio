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

        template <typename value_t, typename rock_coeff_>
        struct degree_computer
        {
            using rock_coeff = rock_coeff_;

            template <typename problem_t, typename state_t>
            static value_t
            rocktrho( problem_t& f, value_t tn, state_t& un )
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

            static std::pair<std::size_t, std::size_t>
            optimal_degree( std::size_t& mdeg )
            {
                std::size_t mz = 1, mr = 1;

                std::size_t iter = 1;
                while ( iter < rock_coeff::ms.size() + 1 && rock_coeff::ms[iter - 1] / mdeg <= 1 )
                {
                    mr = mr + rock_coeff::ms[iter - 1] * 2 - 1;
                    ++iter;
                }
                mdeg = rock_coeff::ms[iter - 1];
                mz   = iter;

                return { mz, mr };
            }

            template <typename problem_t, typename state_t>
            static std::size_t
            compute_n_stages( problem_t& f, value_t tn, state_t& un, value_t& dt )
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
            static std::tuple<std::size_t, std::size_t, std::size_t>
            compute_n_stages_optimal_degree( problem_t& f, value_t tn, state_t& un, value_t& dt )
            {
                std::size_t mdeg = compute_n_stages( f, tn, un, dt );
                auto [mz, mr]    = optimal_degree( mdeg );

                return { mdeg, mz, mr };
            }
        };
    } // namespace detail

    /** @class rock2
     *  @brief define ROCK2
     *
     *  @tparam value_t type of coefficients
     */
    template <typename value_t = double>
    struct rock2
    {
        static constexpr bool is_embedded     = false;
        static constexpr std::size_t N_stages = stages::dynamic;

        using rock_coeff      = rock2_coeff<value_t>;
        using degree_computer = detail::degree_computer<value_t, rock_coeff>;

        rock2()
        {
        }

        template <typename problem_t, typename state_t, typename array_ki_t>
        inline std::tuple<value_t, state_t, value_t>
        operator()( problem_t& f, value_t& tn, state_t& un, array_ki_t& G, value_t& dt )
        {
            auto [mdeg, deg_index, start_index] = degree_computer::compute_n_stages_optimal_degree( f, tn, un, dt );
            G.resize( 3 );

            auto& ujm2 = G[0];
            auto& ujm1 = G[1];
            auto& uj   = G[2];

            uj   = un;
            ujm2 = un;

            value_t const mu1 = rock_coeff::recf[start_index - 1];
            value_t t_jm1     = tn + dt * mu1;
            value_t t_jm2     = tn + dt * mu1;
            value_t t_jm3     = tn;

            ujm1 = un + dt * mu1 * f( tn, un );

            if ( mdeg < 2 )
            {
                uj = ujm1;
            }

            for ( std::size_t j = 2; j < mdeg + 1; ++j )
            {
                value_t const mu    = rock_coeff::recf[start_index + 2 * ( j - 2 ) + 1 - 1];
                value_t const kappa = rock_coeff::recf[start_index + 2 * ( j - 2 ) + 2 - 1];
                value_t const nu    = -1.0 - kappa;

                uj = dt * mu * f( t_jm1, ujm1 ) - nu * ujm1 - kappa * ujm2;

                t_jm1 = dt * mu - nu * t_jm2 - kappa * t_jm3;

                if ( j < mdeg )
                {
                    ujm2 = std::move( ujm1 );
                    ujm1 = std::move( uj );
                }

                t_jm3 = t_jm2;
                t_jm2 = t_jm1;
            }

            value_t const delta_t_1 = dt * rock_coeff::fp1[deg_index - 1]; // equals to $\Delta t \sigma$
            value_t const delta_t_2 = dt * rock_coeff::fp2[deg_index - 1]; // equals to $-\Delta t \sigma(1 - \frac{\tau}{\sigma^2})$

            ujm2 = f( t_jm1, uj );
            ujm1 = uj + delta_t_1 * ujm2;

            t_jm1 = t_jm1 + delta_t_1;

            uj = ujm1 + ( delta_t_1 + delta_t_2 ) * f( t_jm1, ujm1 ) - delta_t_2 * ujm2;

            return { tn + dt, uj, dt };
        }
    };

    /** @class rock4
     *  @brief define ROCK4
     *
     *  @tparam value_t type of coefficients
     */
    template <typename value_t = double>
    struct rock4
    {
        static constexpr bool is_embedded     = false;
        static constexpr std::size_t N_stages = stages::dynamic;

        using rock_coeff      = rock4_coeff<value_t>;
        using degree_computer = detail::degree_computer<value_t, rock_coeff>;

        rock4()
        {
        }

        template <typename problem_t, typename state_t, typename array_ki_t>
        inline std::tuple<value_t, state_t, value_t>
        operator()( problem_t& f, value_t& tn, state_t& un, array_ki_t& G, value_t& dt )
        {
            auto [mdeg, deg_index, start_index] = degree_computer::compute_n_stages_optimal_degree( f, tn, un, dt );
            G.resize( 5 );

            auto& ujm4 = G[0];
            auto& ujm3 = G[1];
            auto& ujm2 = G[2];
            auto& ujm1 = G[3];
            auto& uj   = G[4];

            uj   = un;
            ujm2 = un;

            value_t const mu1 = dt * rock_coeff::recf[start_index - 1];
            value_t t_jm1     = tn + mu1;
            value_t t_jm2     = tn + mu1;
            value_t t_jm3     = tn;

            ujm1 = un + mu1 * f( tn, un );

            // std::cout << "rock2::op() " << mdeg << " " << tn << " " << dt << " " << start_index << " . " << deg_index << "\n";

            if ( mdeg < 2 )
            {
                uj = ujm1;
            }

            for ( std::size_t j = 2; j < mdeg + 1; ++j )
            {
                value_t const mu    = rock_coeff::recf[start_index + 2 * ( j - 2 ) + 1 - 1];
                value_t const kappa = rock_coeff::recf[start_index + 2 * ( j - 2 ) + 2 - 1];
                value_t const nu    = -1.0 - kappa;

                uj = dt * mu * f( t_jm1, ujm1 ) - nu * ujm1 - kappa * ujm2;

                t_jm1 = dt * mu + nu * t_jm2 - kappa * t_jm3;

                if ( j < mdeg )
                {
                    ujm2 = std::move( ujm1 );
                    ujm1 = std::move( uj );
                }

                t_jm3 = t_jm2;
                t_jm2 = t_jm1;
            }

            // the fourth-stage finish procedure

            value_t const a_21 = dt * rock_coeff::fpa[deg_index - 1][0];
            value_t const a_31 = dt * rock_coeff::fpa[deg_index - 1][1];
            value_t const a_32 = dt * rock_coeff::fpa[deg_index - 1][2];
            value_t const a_41 = dt * rock_coeff::fpa[deg_index - 1][3];
            value_t const a_42 = dt * rock_coeff::fpa[deg_index - 1][4];
            value_t const a_43 = dt * rock_coeff::fpa[deg_index - 1][5];
            value_t const b_1  = dt * rock_coeff::fpb[deg_index - 1][0];
            value_t const b_2  = dt * rock_coeff::fpb[deg_index - 1][1];
            value_t const b_3  = dt * rock_coeff::fpb[deg_index - 1][2];
            value_t const b_4  = dt * rock_coeff::fpb[deg_index - 1][3];

            // stage 1.
            ujm1 = f( t_jm1, uj );
            ujm3 = uj + a_21 * ujm1;

            // stage 2.
            t_jm2 = t_jm1 + a_21;
            ujm2  = f( t_jm2, ujm3 );
            ujm4  = uj + a_31 * ujm1 + a_32 * ujm2;

            // stage 3.
            t_jm2 = t_jm1 + a_31 + a_32;
            ujm3  = f( t_jm2, ujm4 );
            ujm4  = uj + a_41 * ujm1 + a_42 * ujm2 + a_43 * ujm3;

            // stage 4.
            t_jm2 = t_jm1 + a_41 + a_42 + a_43;
            uj    = uj + b_1 * ujm1 + b_2 * ujm2 + b_3 * ujm3 + b_4 * f( t_jm2, ujm4 );

            return { tn + dt, uj, dt };
        }
    };

} // namespace ponio::runge_kutta::rock
