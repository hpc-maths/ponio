// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string_view>
#include <utility>

#include "../detail.hpp"
#include "../iteration_info.hpp"
#include "../stage.hpp"
#include "averaged_forces.hpp"
#include "rkc.hpp"

namespace ponio::runge_kutta::chebyshev
{
    namespace mrkc_detail
    {
        template <typename value_t>
        struct mrkc_dimensioning
        {
            std::size_t s = 1;
            std::size_t m = 2;
            value_t eta   = static_cast<value_t>( 0. );
            value_t beta  = static_cast<value_t>( 0. );
        };

        template <typename value_t>
        struct mrkc_error_coefficients
        {
            value_t a1 = static_cast<value_t>( 0. );
            value_t a2 = static_cast<value_t>( 0. );
            value_t a3 = static_cast<value_t>( 0. );
        };

        /**
         * @brief Coefficients of the RKC1 local error estimator used by the
         *        reference mRKC implementation of Rosilho de Souza.
         */
        template <typename value_t>
        mrkc_error_coefficients<value_t>
        error_coefficients( std::size_t s, value_t eps )
        {
            if ( s < 2 )
            {
                throw std::invalid_argument( "mRKC: the embedded estimator requires at least two macro stages." );
            }

            value_t const w0 = rkc_detail::omega_0( s, eps );
            value_t const w1 = rkc_detail::omega_1( s, eps );

            value_t Tjm2 = static_cast<value_t>( 1. );
            value_t Tjm1 = w0;
            value_t bjm2 = static_cast<value_t>( 1. );
            value_t bjm1 = static_cast<value_t>( 1. ) / w0;

            value_t cjm2 = static_cast<value_t>( 0. );
            value_t cjm1 = w1 / w0;
            value_t djm2 = static_cast<value_t>( 0. );
            value_t djm1 = static_cast<value_t>( 0. );

            value_t csm2 = cjm2;
            value_t csm1 = cjm1;
            value_t dsm2 = djm2;
            value_t dsm1 = djm1;
            value_t ds   = djm1;

            for ( std::size_t j = 2; j <= s; ++j )
            {
                value_t const Tj = static_cast<value_t>( 2. ) * w0 * Tjm1 - Tjm2;
                value_t const bj = static_cast<value_t>( 1. ) / Tj;

                value_t const mu    = static_cast<value_t>( 2. ) * w1 * bj / bjm1;
                value_t const nu    = static_cast<value_t>( 2. ) * w0 * bj / bjm1;
                value_t const kappa = -bj / bjm2;

                value_t const cj = mu + nu * cjm1 + kappa * cjm2;
                value_t const dj = static_cast<value_t>( 2. ) * mu * cjm1 + nu * djm1 + kappa * djm2;

                if ( j == s )
                {
                    csm2 = cjm2;
                    csm1 = cjm1;
                    dsm2 = djm2;
                    dsm1 = djm1;
                    ds   = dj;
                }

                cjm2 = cjm1;
                cjm1 = cj;
                djm2 = djm1;
                djm1 = dj;
                bjm2 = bjm1;
                bjm1 = bj;
                Tjm2 = Tjm1;
                Tjm1 = Tj;
            }

            mrkc_error_coefficients<value_t> coeff;
            if ( s > 2 )
            {
                value_t const mult = ( ds - static_cast<value_t>( 1. ) )
                                   / ( ( csm2 - static_cast<value_t>( 1. ) ) * ( dsm1 - ds )
                                       - ( csm1 - static_cast<value_t>( 1. ) ) * ( dsm2 - ds ) );
                coeff.a1 = ( static_cast<value_t>( 1. ) - csm1 ) * mult;
                coeff.a2 = ( csm2 - static_cast<value_t>( 1. ) ) * mult;
                coeff.a3 = ( csm1 - csm2 ) * mult;
            }
            else
            {
                value_t const mult = ( ds - static_cast<value_t>( 1. ) ) / ( ds * csm1 );
                coeff.a1           = ( static_cast<value_t>( 1. ) - csm1 ) * mult;
                coeff.a2           = -mult;
                coeff.a3           = csm1 * mult;
            }
            return coeff;
        }

        /**
         * @brief Compute mRKC parameters for either the general stability
         *        condition or the relaxed scale-separation/refined-mesh branch.
         */
        template <typename value_t>
        mrkc_dimensioning<value_t>
        get_mrkc_stages( value_t dt,
            value_t rho_fast,
            value_t rho_slow,
            value_t eps,
            multirate_dimensioning mode  = multirate_dimensioning::general,
            bool require_error_estimator = false )
        {
            if ( rho_fast < static_cast<value_t>( 0. ) || rho_slow < static_cast<value_t>( 0. ) )
            {
                throw std::invalid_argument( "mRKC: spectral radii must be non-negative." );
            }

            value_t const beta_value = rkc_detail::beta( eps );
            if ( beta_value <= static_cast<value_t>( 0. ) )
            {
                throw std::runtime_error( "mRKC: damping parameter gives a non-positive beta." );
            }

            value_t const tau = std::abs( dt );
            std::size_t s     = static_cast<std::size_t>( std::floor( std::sqrt( tau * rho_slow / beta_value ) ) ) + 1;
            s                 = std::max<std::size_t>( require_error_estimator ? 2 : 1, s );

            value_t const ss = static_cast<value_t>( s ) * static_cast<value_t>( s );
            std::size_t m    = 1;
            value_t eta      = static_cast<value_t>( 0. );

            if ( mode == multirate_dimensioning::general )
            {
                value_t const lower = static_cast<value_t>( 1. )
                                    + static_cast<value_t>( 6. ) * tau * rho_fast / ( beta_value * beta_value * ss );
                m = std::max<std::size_t>( 2, static_cast<std::size_t>( std::floor( std::sqrt( lower ) ) ) + 1 );

                auto eta_for = [&]( std::size_t micro_stages )
                {
                    value_t const mm = static_cast<value_t>( micro_stages ) * static_cast<value_t>( micro_stages );
                    return static_cast<value_t>( 6. ) * tau / ( beta_value * ss ) * mm / ( mm - static_cast<value_t>( 1. ) );
                };

                eta = eta_for( m );
                while ( eta * rho_fast > beta_value * static_cast<value_t>( m ) * static_cast<value_t>( m ) )
                {
                    ++m;
                    eta = eta_for( m );
                }
            }
            else
            {
                // Relaxed condition used for scale separation / refined meshes:
                // eta = 2 h / (beta s^2), eta rho_F <= beta m^2.
                eta = static_cast<value_t>( 2. ) * tau / ( beta_value * ss );
                m   = std::max<std::size_t>( 1, static_cast<std::size_t>( std::ceil( std::sqrt( eta * rho_fast / beta_value ) ) ) );
            }

            return { s, m, eta, beta_value };
        }
    } // namespace mrkc_detail

    /**
     * @class mrkc_impl
     * @brief First-order multirate Runge--Kutta--Chebyshev method.
     */
    template <bool _is_embedded = false, typename _value_t = double>
    struct mrkc_impl
    {
        static constexpr bool is_embedded        = _is_embedded;
        static constexpr bool is_imex_method     = true;
        static constexpr std::size_t N_operators = 2;
        static constexpr std::size_t N_stages    = stages::dynamic;
        static constexpr std::size_t N_storage   = 10;
        static constexpr std::size_t order       = 1;
        static constexpr std::string_view id     = "mRKC";

        using value_t                        = _value_t;
        using dimensioning_t                 = mrkc_detail::mrkc_dimensioning<value_t>;
        static constexpr value_t default_eps = rkc1_default_eps<value_t>;

        struct runtime_stats
        {
            std::size_t attempts   = 0;
            std::size_t accepted   = 0;
            std::size_t rejected   = 0;
            std::size_t s_min      = std::numeric_limits<std::size_t>::max();
            std::size_t s_max      = 0;
            std::size_t m_min      = std::numeric_limits<std::size_t>::max();
            std::size_t m_max      = 0;
            long double s_sum      = 0.;
            long double m_sum      = 0.;
            value_t eta_min        = std::numeric_limits<value_t>::infinity();
            value_t eta_max        = static_cast<value_t>( 0. );
            long double eta_sum    = 0.;
            value_t dt_min         = std::numeric_limits<value_t>::infinity();
            value_t dt_max         = static_cast<value_t>( 0. );
            long double dt_sum     = 0.;
            std::size_t fast_evals = 0;
            std::size_t slow_evals = 0;

            value_t
            s_mean() const
            {
                return attempts ? static_cast<value_t>( s_sum / attempts ) : 0.;
            }

            value_t
            m_mean() const
            {
                return attempts ? static_cast<value_t>( m_sum / attempts ) : 0.;
            }

            value_t
            eta_mean() const
            {
                return attempts ? static_cast<value_t>( eta_sum / attempts ) : 0.;
            }

            value_t
            dt_mean() const
            {
                return attempts ? static_cast<value_t>( dt_sum / attempts ) : 0.;
            }
        };

        value_t rho_fast;
        value_t rho_slow;
        value_t eps;
        multirate_dimensioning dimensioning_mode;
        iteration_info<mrkc_impl> _info;
        std::shared_ptr<runtime_stats> _stats;

        std::size_t last_s = 1;
        std::size_t last_m = 2;
        value_t last_eta   = static_cast<value_t>( 0. );

        mrkc_impl( value_t spectral_radius_fast,
            value_t spectral_radius_slow,
            multirate_dimensioning mode = multirate_dimensioning::general,
            value_t damping             = default_eps )
            : rho_fast( spectral_radius_fast )
            , rho_slow( spectral_radius_slow )
            , eps( damping )
            , dimensioning_mode( mode )
            , _info()
            , _stats( std::make_shared<runtime_stats>() )
        {
            if ( rho_fast < static_cast<value_t>( 0. ) || rho_slow < static_cast<value_t>( 0. ) )
            {
                throw std::invalid_argument( "mRKC: spectral radii must be non-negative." );
            }
            if ( eps < static_cast<value_t>( 0. ) || rkc_detail::beta( eps ) <= static_cast<value_t>( 0. ) )
            {
                throw std::invalid_argument( "mRKC: invalid damping parameter." );
            }
        }

        dimensioning_t
        parameters( value_t dt ) const
        {
            return mrkc_detail::get_mrkc_stages( dt, rho_fast, rho_slow, eps, dimensioning_mode, is_embedded );
        }

        void
        set_spectral_radii( value_t fast, value_t slow )
        {
            if ( fast < static_cast<value_t>( 0. ) || slow < static_cast<value_t>( 0. ) )
            {
                throw std::invalid_argument( "mRKC: spectral radii must be non-negative." );
            }
            rho_fast = fast;
            rho_slow = slow;
        }

        auto const&
        stats() const
        {
            return *_stats;
        }

        void
        reset_stats()
        {
            *_stats = runtime_stats{};
        }

        std::size_t
        last_macro_stages() const
        {
            return last_s;
        }

        std::size_t
        last_micro_stages() const
        {
            return last_m;
        }

        value_t
        last_micro_step() const
        {
            return last_eta;
        }

        value_t
        stability_beta() const
        {
            return rkc_detail::beta( eps );
        }

        value_t
        damping() const
        {
            return eps;
        }

        multirate_dimensioning
        dimensioning() const
        {
            return dimensioning_mode;
        }

        template <bool embedded = is_embedded>
            requires embedded
        auto&
        abs_tol( value_t tol )
        {
            _info.absolute_tolerance = tol;
            return *this;
        }

        template <bool embedded = is_embedded>
            requires embedded
        auto&
        rel_tol( value_t tol )
        {
            _info.relative_tolerance = tol;
            return *this;
        }

      private:

        void
        update_stats( dimensioning_t const& dim, value_t step_dt, bool success, std::size_t fast_evals, std::size_t slow_evals )
        {
            auto& st = *_stats;
            ++st.attempts;
            success ? ++st.accepted : ++st.rejected;
            st.s_min = std::min( st.s_min, dim.s );
            st.s_max = std::max( st.s_max, dim.s );
            st.s_sum += static_cast<long double>( dim.s );
            st.m_min = std::min( st.m_min, dim.m );
            st.m_max = std::max( st.m_max, dim.m );
            st.m_sum += static_cast<long double>( dim.m );
            st.eta_min = std::min( st.eta_min, dim.eta );
            st.eta_max = std::max( st.eta_max, dim.eta );
            st.eta_sum += dim.eta;
            st.dt_min = std::min( st.dt_min, step_dt );
            st.dt_max = std::max( st.dt_max, step_dt );
            st.dt_sum += step_dt;
            st.fast_evals += fast_evals;
            st.slow_evals += slow_evals;
        }

      public:

        template <typename problem_t, typename state_t, typename array_ki_t>
        void
        operator()( problem_t& pb, value_t& tn, state_t& un, array_ki_t& U, value_t& dt, state_t& unp1 )
        {
            _info.reset_eval();
            value_t const step_dt    = dt;
            dimensioning_t const dim = parameters( step_dt );
            last_s                   = dim.s;
            last_m                   = dim.m;
            last_eta                 = dim.eta;
            _info.number_of_stages   = dim.s;

            auto averaged_rhs = [&]( value_t t, state_t& y, state_t& f_bar )
            {
                multirate_detail::averaged_force_order_1( pb, t, y, U, dim.m, dim.eta, eps, f_bar );
            };

            rkc_detail::rkc1_step( averaged_rhs, tn, un, step_dt, dim.s, eps, U[0], U[1], U[2], U[3], unp1 );

            _info.number_of_eval[0] = dim.s * dim.m;
            _info.number_of_eval[1] = dim.s;

            if constexpr ( is_embedded )
            {
                auto const coeff = mrkc_detail::error_coefficients<value_t>( dim.s, eps );
                U[3]             = coeff.a1 * U[0] + coeff.a2 * U[1] + coeff.a3 * unp1;
                auto const err2  = ::ponio::detail::error_algebra<state_t>::estimate_squared( U[3],
                    unp1,
                    un,
                    _info.absolute_tolerance,
                    _info.relative_tolerance );
                _info.error      = std::sqrt( err2 );
                _info.success    = _info.error <= static_cast<value_t>( 1. );

                value_t fac = static_cast<value_t>( 2. );
                if ( _info.error > static_cast<value_t>( 0. ) )
                {
                    fac = std::min( static_cast<value_t>( 2. ),
                        std::max( static_cast<value_t>( 0.5 ), static_cast<value_t>( 1. ) / _info.error ) );
                }
                value_t const new_dt = static_cast<value_t>( 0.8 ) * fac * step_dt;

                if ( _info.success )
                {
                    tn += step_dt;
                }
                else
                {
                    std::swap( un, unp1 );
                }
                dt = new_dt;
            }
            else
            {
                _info.error   = static_cast<value_t>( 0. );
                _info.success = true;
                tn += step_dt;
            }

            update_stats( dim, step_dt, _info.success, _info.number_of_eval[0], _info.number_of_eval[1] );
        }

        auto&
        info()
        {
            return _info;
        }

        auto const&
        info() const
        {
            return _info;
        }
    };

    template <bool is_embedded = false, typename value_t = double>
    auto
    mrkc( value_t rho_fast,
        value_t rho_slow,
        multirate_dimensioning mode = multirate_dimensioning::general,
        value_t eps                 = mrkc_impl<is_embedded, value_t>::default_eps )
    {
        return mrkc_impl<is_embedded, value_t>( rho_fast, rho_slow, mode, eps );
    }

} // namespace ponio::runge_kutta::chebyshev
