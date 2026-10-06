// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private

#pragma once

#include <cstddef>
#include <stdexcept>
#include <type_traits>

#include "../linear_algebra.hpp"
#include "rkc.hpp"

namespace ponio::runge_kutta::chebyshev
{
    /**
     * @brief Dimensioning strategy used by the explicit stabilized multirate methods.
     */
    enum class multirate_dimensioning
    {
        general,
        relaxed
    };

    namespace multirate_detail
    {

        /**
         * @brief Compute the RKC1 coefficient used by the second-order averaged force.
         */
        template <typename value_t>
        value_t
        rkc1_alpha_m( std::size_t m, value_t eps )
        {
            if ( m == 0 )
            {
                throw std::invalid_argument( "mROCK2: RKC1 micro stage number must be positive." );
            }

            value_t const w0 = rkc_detail::omega_0( m, eps );
            value_t const w1 = rkc_detail::omega_1( m, eps );

            value_t Tjm2   = static_cast<value_t>( 1. );
            value_t Tjm1   = w0;
            value_t dTjm2  = static_cast<value_t>( 0. );
            value_t dTjm1  = static_cast<value_t>( 1. );
            value_t ddTjm2 = static_cast<value_t>( 0. );
            value_t ddTjm1 = static_cast<value_t>( 0. );

            if ( m == 1 )
            {
                return static_cast<value_t>( 0. );
            }

            for ( std::size_t j = 2; j <= m; ++j )
            {
                value_t const Tj   = static_cast<value_t>( 2. ) * w0 * Tjm1 - Tjm2;
                value_t const dTj  = static_cast<value_t>( 2. ) * Tjm1 + static_cast<value_t>( 2. ) * w0 * dTjm1 - dTjm2;
                value_t const ddTj = static_cast<value_t>( 4. ) * dTjm1 + static_cast<value_t>( 2. ) * w0 * ddTjm1 - ddTjm2;

                Tjm2   = Tjm1;
                Tjm1   = Tj;
                dTjm2  = dTjm1;
                dTjm1  = dTj;
                ddTjm2 = ddTjm1;
                ddTjm1 = ddTj;
            }

            return w1 * w1 * ddTjm1 / Tjm1;
        }

        /**
         * @brief First-order averaged force used by mRKC.
         *
         * The slow component is frozen at the beginning of the macro stage and
         * the auxiliary fast problem is integrated over the micro interval
         * with RKC1.
         */
        template <typename problem_t, typename value_t, typename state_t, typename array_ki_t>
        void
        averaged_force_order_1( problem_t& pb, value_t t, state_t& y, array_ki_t& U, std::size_t m, value_t eta, value_t eps, state_t& f_bar )
        {
            using fast_op         = std::integral_constant<std::size_t, 0>;
            using slow_op         = std::integral_constant<std::size_t, 1>;
            using state_algebra_t = ::ponio::linear_algebra::state_algebra<state_t>;

            auto& micro_kjm2  = U[4];
            auto& micro_kjm1  = U[5];
            auto& micro_kj    = U[6];
            auto& micro_f_tmp = U[7];
            auto& slow_frozen = U[8];
            auto& u_eta       = U[9];

            pb( slow_op(), t, y, slow_frozen );

            auto auxiliary_rhs = [&]( value_t r, state_t& u, state_t& du )
            {
                pb( fast_op(), r, u, du );
                state_algebra_t::add( du, slow_frozen );
            };

            rkc_detail::rkc1_step( auxiliary_rhs, t, y, eta, m, eps, micro_kjm2, micro_kjm1, micro_kj, micro_f_tmp, u_eta );

            f_bar = ( u_eta - y ) / eta;
        }

        /**
         * @brief Second-order averaged force used by mROCK2.
         *
         * Two RKC1 micro integrations are performed. The second one uses the
         * shifted fast state required by the second-order averaged-force
         * construction.
         */
        template <typename problem_t, typename value_t, typename state_t, typename array_ki_t>
        void
        averaged_force_order_2( problem_t& pb,
            value_t t,
            state_t& y,
            array_ki_t& U,
            std::size_t m,
            value_t eta,
            value_t alpha_m,
            value_t eps,
            state_t& f_bar_2 )
        {
            using fast_op         = std::integral_constant<std::size_t, 0>;
            using slow_op         = std::integral_constant<std::size_t, 1>;
            using state_algebra_t = ::ponio::linear_algebra::state_algebra<state_t>;

            auto& micro_kjm2  = U[4];
            auto& micro_kjm1  = U[5];
            auto& micro_kj    = U[6];
            auto& micro_f_tmp = U[7];
            auto& slow_frozen = U[8];
            auto& u_eta       = U[9];
            auto& f_bar_1     = U[10];
            auto& shifted     = U[11];
            auto& v_eta       = U[12];

            pb( slow_op(), t, y, slow_frozen );

            auto first_auxiliary_rhs = [&]( value_t r, state_t& u, state_t& du )
            {
                pb( fast_op(), r, u, du );
                state_algebra_t::add( du, slow_frozen );
            };

            rkc_detail::rkc1_step( first_auxiliary_rhs, t, y, eta, m, eps, micro_kjm2, micro_kjm1, micro_kj, micro_f_tmp, u_eta );

            f_bar_1 = ( u_eta - y ) / eta;

            value_t const shift_coeff = static_cast<value_t>( 0.5 ) * alpha_m * eta;

            auto second_auxiliary_rhs = [&]( value_t r, state_t& v, state_t& dv )
            {
                shifted = v - shift_coeff * f_bar_1;
                pb( fast_op(), r, shifted, dv );
                state_algebra_t::add( dv, slow_frozen );
            };

            rkc_detail::rkc1_step( second_auxiliary_rhs, t, y, eta, m, eps, micro_kjm2, micro_kjm1, micro_kj, micro_f_tmp, v_eta );

            f_bar_2 = ( v_eta - y ) / eta;
        }
    } // namespace multirate_detail

} // namespace ponio::runge_kutta::chebyshev
