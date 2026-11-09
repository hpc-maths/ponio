// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private

#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>

namespace ponio::runge_kutta::legendre::dynamic
{
    /**
     * @brief Runtime coefficients of the second-order Runge--Kutta--Legendre method.
     *
     * Dynamic counterpart of explicit_rkl2<N_stages> from rkl.hpp.
     */
    template <typename value_t = double>
    struct rkl2_coefficients
    {
        std::size_t stages;
        value_t w1;

        explicit rkl2_coefficients( std::size_t s )
            : stages( s )
            , w1( compute_w1( s ) )
        {
            if ( stages < 2 )
            {
                throw std::invalid_argument( "dynamic RKL2 requires at least two stages" );
            }
        }

        static value_t
        compute_w1( std::size_t s )
        {
            if ( s < 2 )
            {
                throw std::invalid_argument( "dynamic RKL2 requires at least two stages" );
            }

            value_t const sv = static_cast<value_t>( s );
            return static_cast<value_t>( 4 ) / ( sv * sv + sv - static_cast<value_t>( 2 ) );
        }

        static value_t
        b( std::size_t j )
        {
            if ( j <= 1 )
            {
                return static_cast<value_t>( 1 ) / static_cast<value_t>( 3 );
            }

            value_t const jv = static_cast<value_t>( j );
            return ( jv * jv + jv - static_cast<value_t>( 2 ) ) / ( static_cast<value_t>( 2 ) * jv * ( jv + static_cast<value_t>( 1 ) ) );
        }

        static value_t
        a( std::size_t j )
        {
            return static_cast<value_t>( 1 ) - b( j );
        }

        value_t
        mu( std::size_t j ) const
        {
            if ( j < 2 )
            {
                throw std::invalid_argument( "mu_j is defined here for j >= 2" );
            }

            value_t const jv = static_cast<value_t>( j );
            return ( static_cast<value_t>( 2 * j - 1 ) / jv ) * b( j ) / b( j - 1 );
        }

        value_t
        nu( std::size_t j ) const
        {
            if ( j < 2 )
            {
                throw std::invalid_argument( "nu_j is defined here for j >= 2" );
            }

            value_t const jv = static_cast<value_t>( j );
            return -( static_cast<value_t>( j - 1 ) / jv ) * b( j ) / b( j - 2 );
        }

        value_t
        mu_t( std::size_t j ) const
        {
            if ( j == 0 )
            {
                throw std::invalid_argument( "mu_t_j is defined for j >= 1" );
            }
            if ( j == 1 )
            {
                return b( 1 ) * w1;
            }
            return mu( j ) * w1;
        }

        value_t
        gamma_t( std::size_t j ) const
        {
            if ( j < 2 )
            {
                throw std::invalid_argument( "gamma_t_j is defined here for j >= 2" );
            }
            return -a( j - 1 ) * mu_t( j );
        }
    };

    /**
     * @brief Evaluate L_s(x) and L'_s(x) with the three-term Legendre recurrence.
     */
    template <typename value_t = double>
    struct legendre_value_derivative
    {
        value_t value;
        value_t derivative;
    };

    template <typename value_t = double>
    legendre_value_derivative<value_t>
    legendre( std::size_t s, value_t x )
    {
        if ( s == 0 )
        {
            return { static_cast<value_t>( 1 ), static_cast<value_t>( 0 ) };
        }
        if ( s == 1 )
        {
            return { x, static_cast<value_t>( 1 ) };
        }

        value_t p_jm2  = static_cast<value_t>( 1 );
        value_t dp_jm2 = static_cast<value_t>( 0 );
        value_t p_jm1  = x;
        value_t dp_jm1 = static_cast<value_t>( 1 );

        for ( std::size_t j = 2; j <= s; ++j )
        {
            value_t const jv = static_cast<value_t>( j );
            value_t const a  = static_cast<value_t>( 2 * j - 1 ) / jv;
            value_t const b  = static_cast<value_t>( j - 1 ) / jv;

            value_t const p_j  = a * x * p_jm1 - b * p_jm2;
            value_t const dp_j = a * ( p_jm1 + x * dp_jm1 ) - b * dp_jm2;

            p_jm2  = p_jm1;
            dp_jm2 = dp_jm1;
            p_jm1  = p_j;
            dp_jm1 = dp_j;
        }

        return { p_jm1, dp_jm1 };
    }

    /**
     * @brief Parameters of the damped Legendre branch filter used by PIRKL-D.
     *
     * P_br(z) = L_s(w0 + w1 z) / L_s(w0),
     * w0 = 1 + delta/s^2,
     * w1 = L_s(w0)/(2 L'_s(w0)).
     * Hence P'_br(0)=1/2 and beta=0 in the DA coupling.
     */
    template <typename value_t = double>
    struct legendre_filter_parameters
    {
        std::size_t stages;
        value_t delta;
        value_t w0;
        value_t w1;
        value_t Ls_w0;
        value_t dLs_w0;
        value_t stability_interval;
    };

    template <typename value_t = double>
    legendre_filter_parameters<value_t>
    make_legendre_filter_parameters( std::size_t s, value_t delta )
    {
        if ( s == 0 )
        {
            throw std::invalid_argument( "the damped Legendre filter requires at least one stage" );
        }
        if ( delta < static_cast<value_t>( 0 ) )
        {
            throw std::invalid_argument( "the damping parameter delta must be non-negative" );
        }

        value_t const sv = static_cast<value_t>( s );
        value_t const w0 = static_cast<value_t>( 1 ) + delta / ( sv * sv );
        auto const ld    = legendre<value_t>( s, w0 );

        if ( ld.derivative == static_cast<value_t>( 0 ) )
        {
            throw std::runtime_error( "zero derivative while building the damped Legendre filter" );
        }

        value_t const w1 = ld.value / ( static_cast<value_t>( 2 ) * ld.derivative );
        if ( !( w1 > static_cast<value_t>( 0 ) ) )
        {
            throw std::runtime_error( "non-positive w1 while building the damped Legendre filter" );
        }

        return { s, delta, w0, w1, ld.value, ld.derivative, ( w0 + static_cast<value_t>( 1 ) ) / w1 };
    }

    /**
     * @brief Runtime normalized coefficients of one damped Legendre-filter stage.
     */
    template <typename value_t = double>
    struct legendre_filter_stage_coefficients
    {
        value_t mu;
        value_t nu;
        value_t mu_t;
    };

    template <typename value_t = double>
    legendre_filter_stage_coefficients<value_t>
    make_legendre_filter_stage_coefficients( std::size_t j, value_t w0, value_t w1, value_t c_jm2, value_t c_jm1, value_t c_j )
    {
        if ( j < 2 )
        {
            throw std::invalid_argument( "damped Legendre stage coefficients require j >= 2" );
        }
        if ( c_j == static_cast<value_t>( 0 ) )
        {
            throw std::runtime_error( "zero Legendre normalization coefficient" );
        }

        value_t const jv     = static_cast<value_t>( j );
        value_t const factor = static_cast<value_t>( 2 * j - 1 ) / jv;

        value_t const mu   = factor * w0 * c_jm1 / c_j;
        value_t const nu   = -( static_cast<value_t>( j - 1 ) / jv ) * c_jm2 / c_j;
        value_t const mu_t = factor * w1 * c_jm1 / c_j;

        return { mu, nu, mu_t };
    }

    /**
     * @brief Smallest RKL2 stage number whose stability interval covers zmax.
     */
    template <typename value_t = double>
    std::size_t
    required_rkl2_stages( value_t zmax )
    {
        if ( zmax < static_cast<value_t>( 0 ) )
        {
            throw std::invalid_argument( "zmax must be non-negative" );
        }

        value_t const s_real = ( std::sqrt( static_cast<value_t>( 9 ) + static_cast<value_t>( 8 ) * zmax ) - static_cast<value_t>( 1 ) )
                             / static_cast<value_t>( 2 );
        std::size_t const s = static_cast<std::size_t>( std::ceil( s_real ) );
        return s < 2 ? 2 : s;
    }

    /**
     * @brief Portable approximation of the modified Bessel function I_0(x).
     *
     * Piecewise polynomial/asymptotic approximation that avoids depending on
     * std::cyl_bessel_i.
     */
    template <typename value_t = double>
    value_t
    modified_bessel_i0( value_t x )
    {
        value_t const ax = std::abs( x );

        if ( ax < static_cast<value_t>( 3.75 ) )
        {
            value_t const y  = x / static_cast<value_t>( 3.75 );
            value_t const y2 = y * y;

            return static_cast<value_t>( 1.0 )
                 + y2
                       * ( static_cast<value_t>( 3.5156229 )
                           + y2
                                 * ( static_cast<value_t>( 3.0899424 )
                                     + y2
                                           * ( static_cast<value_t>( 1.2067492 )
                                               + y2
                                                     * ( static_cast<value_t>( 0.2659732 )
                                                         + y2
                                                               * ( static_cast<value_t>( 0.0360768 )
                                                                   + y2 * static_cast<value_t>( 0.0045813 ) ) ) ) ) );
        }

        value_t const y = static_cast<value_t>( 3.75 ) / ax;
        value_t const
            poly = static_cast<value_t>( 0.39894228 )
                 + y
                       * ( static_cast<value_t>( 0.01328592 )
                           + y
                                 * ( static_cast<value_t>( 0.00225319 )
                                     + y
                                           * ( static_cast<value_t>( -0.00157565 )
                                               + y
                                                     * ( static_cast<value_t>( 0.00916281 )
                                                         + y
                                                               * ( static_cast<value_t>( -0.02057706 )
                                                                   + y
                                                                         * ( static_cast<value_t>( 0.02635537 )
                                                                             + y
                                                                                   * ( static_cast<value_t>( -0.01647633 )
                                                                                       + y * static_cast<value_t>( 0.00392377 ) ) ) ) ) ) ) );

        return std::exp( ax ) / std::sqrt( ax ) * poly;
    }

    /**
     * @brief Solve I_0(t)=target for t >= 0 by bracketed bisection.
     */
    template <typename value_t = double>
    value_t
    inverse_modified_bessel_i0( value_t target )
    {
        if ( target <= static_cast<value_t>( 1 ) )
        {
            return static_cast<value_t>( 0 );
        }

        auto i0 = []( value_t x ) -> value_t
        {
            return modified_bessel_i0<value_t>( x );
        };

        value_t lo = static_cast<value_t>( 0 );
        value_t hi = static_cast<value_t>( 1 );

        while ( i0( hi ) < target )
        {
            hi *= static_cast<value_t>( 2 );
            if ( !std::isfinite( hi ) || hi > static_cast<value_t>( 1.0e3 ) )
            {
                throw std::runtime_error( "unable to bracket inverse I0" );
            }
        }

        value_t const eps = static_cast<value_t>( 16 ) * std::numeric_limits<value_t>::epsilon();
        for ( std::size_t iter = 0; iter < 128; ++iter )
        {
            value_t const mid = static_cast<value_t>( 0.5 ) * ( lo + hi );
            if ( i0( mid ) < target )
            {
                lo = mid;
            }
            else
            {
                hi = mid;
            }

            if ( hi - lo <= eps * ( static_cast<value_t>( 1 ) + hi ) )
            {
                break;
            }
        }

        return static_cast<value_t>( 0.5 ) * ( lo + hi );
    }

    /**
     * @brief Compute the PIRKL-D damping parameter from I_0(t) >= K zmax, delta = t^2 / 2.
     */
    template <typename value_t = double>
    value_t
    damping_from_zmax( value_t zmax, value_t K = static_cast<value_t>( 10 ) )
    {
        if ( zmax < static_cast<value_t>( 0 ) )
        {
            throw std::invalid_argument( "zmax must be non-negative" );
        }
        if ( !( K > static_cast<value_t>( 0 ) ) )
        {
            throw std::invalid_argument( "K must be strictly positive" );
        }

        value_t const t = inverse_modified_bessel_i0<value_t>( K * zmax );
        return static_cast<value_t>( 0.5 ) * t * t;
    }

    /**
     * @brief Smallest damped-filter stage number whose exact interval covers zmax.
     */
    template <typename value_t = double>
    std::size_t
    required_legendre_filter_stages( value_t zmax, value_t delta, std::size_t s_max = 100000 )
    {
        if ( zmax < static_cast<value_t>( 0 ) )
        {
            throw std::invalid_argument( "zmax must be non-negative" );
        }

        for ( std::size_t s = 1; s <= s_max; ++s )
        {
            if ( make_legendre_filter_parameters<value_t>( s, delta ).stability_interval >= zmax )
            {
                return s;
            }
        }

        throw std::runtime_error( "unable to cover zmax with the requested maximum number of Legendre-filter stages" );
    }

    /**
     * @brief Runtime dimensioning data for PIRKL-D.
     */
    template <typename value_t = double>
    struct pirkl_d_dimensioning
    {
        value_t zmax;
        value_t delta;
        std::size_t branch_stages;
        std::size_t closure_stages;
    };

    template <typename value_t = double>
    pirkl_d_dimensioning<value_t>
    dimension_pirkl_d( value_t dt,
        value_t rho_D,
        value_t safety    = static_cast<value_t>( 1.1 ),
        value_t K         = static_cast<value_t>( 10 ),
        std::size_t s_max = 100000 )
    {
        if ( dt < static_cast<value_t>( 0 ) )
        {
            throw std::invalid_argument( "dt must be non-negative" );
        }
        if ( rho_D < static_cast<value_t>( 0 ) )
        {
            throw std::invalid_argument( "rho_D must be non-negative" );
        }
        if ( !( safety > static_cast<value_t>( 0 ) ) )
        {
            throw std::invalid_argument( "safety must be strictly positive" );
        }

        value_t const zmax          = safety * dt * rho_D;
        value_t const delta         = damping_from_zmax<value_t>( zmax, K );
        std::size_t const s_branch  = required_legendre_filter_stages<value_t>( zmax, delta, s_max );
        std::size_t const s_closure = required_rkl2_stages<value_t>( zmax );

        return { zmax, delta, s_branch, s_closure };
    }

    /**
     * @brief Apply RKL2 with a runtime number of stages.
     *
     * The caller supplies F_D(t_n,u_n) so PIRKL-D can share this evaluation
     * with the damped Legendre branch. The remaining diffusion evaluations are
     * performed through evaluator(t, u, out).
     */
    template <typename evaluator_t, typename state_t, typename value_t = double>
    void
    apply_rkl2( evaluator_t&& evaluator,
        value_t tn,
        state_t const& un,
        value_t dt,
        std::size_t stages,
        state_t const& f_start,
        state_t& y_jm2,
        state_t& y_jm1,
        state_t& y_j,
        state_t& f_tmp,
        state_t& out )
    {
        rkl2_coefficients<value_t> const coeff( stages );

        y_jm2 = un;
        y_jm1 = un + coeff.mu_t( 1 ) * dt * f_start;

        for ( std::size_t j = 2; j <= stages; ++j )
        {
            evaluator( tn, y_jm1, f_tmp );

            value_t const mu      = coeff.mu( j );
            value_t const nu      = coeff.nu( j );
            value_t const mu_t    = coeff.mu_t( j );
            value_t const gamma_t = coeff.gamma_t( j );

            y_j = mu * y_jm1 + nu * y_jm2 + ( static_cast<value_t>( 1 ) - mu - nu ) * un + mu_t * dt * f_tmp + gamma_t * dt * f_start;

            if ( j < stages )
            {
                std::swap( y_jm2, y_jm1 );
                std::swap( y_jm1, y_j );
            }
        }

        out = y_j;
    }

    /**
     * @brief Apply the damped Legendre branch filter of PIRKL-D.
     *
     * The caller supplies F_D(t_n,u_n) so PIRKL-D can share this evaluation
     * with the RKL2 closure. The normalized Legendre coefficients c_j=L_j(w0)
     * are generated on the fly and only three state work arrays are needed.
     */
    template <typename evaluator_t, typename state_t, typename value_t = double>
    void
    apply_legendre_filter( evaluator_t&& evaluator,
        value_t tn,
        state_t const& un,
        value_t dt,
        legendre_filter_parameters<value_t> const& params,
        state_t const& f_start,
        state_t& y_jm2,
        state_t& y_jm1,
        state_t& y_j,
        state_t& f_tmp,
        state_t& out )
    {
        std::size_t const stages = params.stages;

        y_jm2 = un;
        y_jm1 = un + ( params.w1 / params.w0 ) * dt * f_start;

        if ( stages == 1 )
        {
            out = y_jm1;
            return;
        }

        value_t c_jm2 = static_cast<value_t>( 1 );
        value_t c_jm1 = params.w0;

        for ( std::size_t j = 2; j <= stages; ++j )
        {
            value_t const jv  = static_cast<value_t>( j );
            value_t const c_j = ( static_cast<value_t>( 2 * j - 1 ) * params.w0 * c_jm1 - static_cast<value_t>( j - 1 ) * c_jm2 ) / jv;

            auto const coeff = make_legendre_filter_stage_coefficients<value_t>( j, params.w0, params.w1, c_jm2, c_jm1, c_j );

            evaluator( tn, y_jm1, f_tmp );
            y_j = coeff.mu * y_jm1 + coeff.nu * y_jm2 + coeff.mu_t * dt * f_tmp;

            if ( j < stages )
            {
                std::swap( y_jm2, y_jm1 );
                std::swap( y_jm1, y_j );
            }

            c_jm2 = c_jm1;
            c_jm1 = c_j;
        }

        out = y_j;
    }

    /**
     * @brief Approximate diffusion error for adaptive PIRKL-D.
     *
     * @warning This estimator is provisional and is not a dedicated embedded
     * RKL2 estimator.
     */
    template <typename state_t, typename value_t = double>
    void
    approximate_rkl2_diffusion_error( state_t const& f_start, state_t const& f_end, value_t dt, state_t& err_D )
    {
        err_D = static_cast<value_t>( 0.5 ) * dt * ( f_end - f_start );
    }

} // namespace ponio::runge_kutta::legendre::dynamic
