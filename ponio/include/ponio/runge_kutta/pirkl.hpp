// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private

#pragma once

// NOLINTBEGIN(misc-include-cleaner)

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <numbers>
#include <string_view>
#include <type_traits>
#include <utility>

#include "../detail.hpp"
#include "../iteration_info.hpp"
#include "../stage.hpp"
#include "rkl_d.hpp"

// NOLINTEND(misc-include-cleaner)

namespace ponio::runge_kutta::pirkl
{
    /**
     * @class pirkl_DA_impl
     * @brief PIRKL-D method for diffusion-advection problems.
     *
     * The diffusion branch is stabilized with a damped Legendre polynomial,
     * while the second-order diffusion update is computed with a dynamic RKL2
     * recurrence. The diffusion-advection coupling is the PIROCK-DA finishing
     * procedure with beta = 0.
     *
     * @tparam eig_computer_t diffusion spectral-radius provider
     * @tparam _is_embedded   adaptive time-step method when true
     * @tparam _value_t       coefficient type
     */
    template <typename eig_computer_t, bool _is_embedded = false, typename _value_t = double>
    struct pirkl_DA_impl
    {
        static constexpr bool is_embedded = _is_embedded;

        static constexpr bool is_imex_method = true;
        // number_of_eval counts one evaluation per operator: 0 diffusion, 1 advection.
        static constexpr std::size_t N_operators = 2;
        static constexpr std::size_t N_stages    = stages::dynamic;

        // Workspace layout (fixed and adaptive):
        //   0..2  : three-state recurrence workspace (shared by branch and closure)
        //   3     : diffusion temporary
        //   4     : F_D(t_n,u_n), cached across an adaptive retry
        //   5     : branch state K
        //   6     : second-order RKL2 diffusion closure
        //   7..9  : advection evaluations
        //   10..12: DA coupling stages sp3, sp4, sp5
        //   13    : temporary error vector
        static constexpr std::size_t N_storage = 14;

        static constexpr std::size_t order   = 2;
        static constexpr std::string_view id = "PIRKL-D";

        using value_t        = _value_t;
        using dimensioning_t = legendre::dynamic::pirkl_d_dimensioning<value_t>;

        eig_computer_t eig_computer;
        iteration_info<pirkl_DA_impl> _info;

        // Adaptive time-step controller.
        value_t facmax                 = static_cast<value_t>( 5. );
        value_t facd                   = static_cast<value_t>( 0.8 );
        std::size_t attempt_count      = 0;
        std::size_t nrejfac            = 0;
        bool previous_attempt_rejected = false;

        value_t errp = static_cast<value_t>( 0. );
        value_t hp   = static_cast<value_t>( 0. );

        value_t told     = static_cast<value_t>( 0. );
        std::size_t nrej = 0;

        // PIRKL-D dimensioning parameters.
        value_t stability_safety      = static_cast<value_t>( 1.1 );
        value_t damping_K             = static_cast<value_t>( 10. );
        std::size_t max_filter_stages = 100000;

        // Last dimensioning data.
        value_t last_rho_D = static_cast<value_t>( 0. );
        dimensioning_t last_dimensioning{ static_cast<value_t>( 0. ), static_cast<value_t>( 0. ), 0, 0 };

        pirkl_DA_impl() = default;

        explicit pirkl_DA_impl( eig_computer_t&& _eig_computer )
            : eig_computer( std::forward<eig_computer_t>( _eig_computer ) )
            , _info()
        {
        }

        /**
         * @brief One PIRKL-D diffusion-advection step.
         *
         * problem_t::system must contain:
         *   0: diffusion operator F_D
         *   1: advection operator F_A
         */
        template <typename problem_t, typename state_t, typename array_ki_t>
        void
        operator()( problem_t& pb, value_t& tn, state_t& un, array_ki_t& U, value_t& dt, state_t& u_np1 )
        {
            using diffusion_op = std::integral_constant<std::size_t, 0>;
            using advection_op = std::integral_constant<std::size_t, 1>;

            auto& y_j        = U[0];
            auto& y_jm1      = U[1];
            auto& y_jm2      = U[2];
            auto& fd_tmp     = U[3];
            auto& fd_start   = U[4];
            auto& K          = U[5];
            auto& u_diff     = U[6];
            auto& fa_K       = U[7];
            auto& fa_tmp     = U[8];
            auto& fa_tmp_bis = U[9];
            auto& u_sp3      = U[10];
            auto& u_sp4      = U[11];
            auto& u_sp5      = U[12];

            _info.reset_eval();

            auto eval_diffusion = [&]( value_t time, state_t const& state, state_t& output )
            {
                pb( diffusion_op(), time, state, output );
                ++_info.number_of_eval[0];
            };

            auto eval_advection = [&]( value_t time, state_t const& state, state_t& output )
            {
                pb( advection_op(), time, state, output );
                ++_info.number_of_eval[1];
            };

            auto normalized_error_squared = [&]( state_t const& error, state_t const& state_old, state_t const& state_new ) -> value_t
            {
                value_t error_squared = static_cast<value_t>( 0. );
                auto const n          = error.size();

                for ( std::size_t i = 0; i < static_cast<std::size_t>( n ); ++i )
                {
                    value_t const scale = _info.absolute_tolerance
                                        + _info.relative_tolerance * std::max( std::abs( state_old[i] ), std::abs( state_new[i] ) );

                    value_t const normalized_error = error[i] / scale;
                    error_squared += normalized_error * normalized_error;
                }

                return error_squared / static_cast<value_t>( n );
            };

            // The diffusion spectral radius is supplied by the user-provided evaluator.
            last_rho_D = eig_computer( eval_diffusion, tn, un, dt, U );
            last_dimensioning = legendre::dynamic::dimension_pirkl_d<value_t>( dt, last_rho_D, stability_safety, damping_K, max_filter_stages );

            std::size_t const s_branch  = last_dimensioning.branch_stages;
            std::size_t const s_closure = last_dimensioning.closure_stages;

            _info.number_of_stages = s_branch + s_closure + 3;

            value_t const gamma = static_cast<value_t>( 1. ) - static_cast<value_t>( 0.5 ) * std::numbers::sqrt2_v<value_t>;

            // P'_br(0)=1/2 by construction, hence beta=0 exactly.
            value_t const t_diffusion_coupling = tn + static_cast<value_t>( 0.5 ) * dt;
            value_t const t_advection_K        = tn;
            value_t const t_advection_sp4      = tn + dt / static_cast<value_t>( 3. );
            value_t const t_advection_sp5      = tn + static_cast<value_t>( 2. ) * dt / static_cast<value_t>( 3. );

            // Adaptive controller.
            auto raw_fac_from_current_error = [&]() -> value_t
            {
                ++attempt_count;

                if ( attempt_count % 100 == 0 )
                {
                    if ( nrejfac >= 10 )
                    {
                        facd *= static_cast<value_t>( 0.8 );
                    }
                    if ( nrejfac == 0 )
                    {
                        facd *= static_cast<value_t>( 1.02 );
                    }
                    facd    = std::min( static_cast<value_t>( 0.8 ), std::max( static_cast<value_t>( 0.4 ), facd ) );
                    nrejfac = 0;
                }

                // _info.error stores the square of the effective RMS error.
                return std::pow( _info.error, static_cast<value_t>( -0.25 ) );
            };

            auto finalize_fac = [&]( value_t fac ) -> value_t
            {
                if ( previous_attempt_rejected )
                {
                    facmax = static_cast<value_t>( 1. );
                }

                return std::min( facmax, std::max( static_cast<value_t>( 0.1 ), facd * fac ) );
            };

            auto finalize_rejection_dt = [&]( value_t new_dt ) -> value_t
            {
                value_t retry_dt = static_cast<value_t>( 0.8 ) * new_dt;

                if ( told == tn )
                {
                    ++nrej;
                    if ( nrej == 10 )
                    {
                        retry_dt = static_cast<value_t>( 1.0e-5 );
                    }
                }
                told = tn;

                ++nrejfac;
                previous_attempt_rejected = true;

                return retry_dt;
            };

            // F_D(t_n, u_n) is shared by the two Legendre recurrences and can
            // be reused after a rejected step.
            if constexpr ( is_embedded )
            {
                if ( !previous_attempt_rejected )
                {
                    eval_diffusion( tn, un, fd_start );
                }
            }
            else
            {
                eval_diffusion( tn, un, fd_start );
            }

            // -----------------------------------------------------------------
            // 1. Damped Legendre branch filter: K = P_br(dt F_D) u^n.
            // -----------------------------------------------------------------
            auto const branch_parameters = legendre::dynamic::make_legendre_filter_parameters<value_t>( s_branch, last_dimensioning.delta );

            legendre::dynamic::apply_legendre_filter( eval_diffusion, tn, un, dt, branch_parameters, fd_start, y_jm2, y_jm1, y_j, fd_tmp, K );

            // -----------------------------------------------------------------
            // 2. Independent dynamic RKL2 closure: u_diff.
            // -----------------------------------------------------------------
            legendre::dynamic::apply_rkl2( eval_diffusion, tn, un, dt, s_closure, fd_start, y_jm2, y_jm1, y_j, fd_tmp, u_diff );

            value_t err_D_scalar = static_cast<value_t>( 0. );

            if constexpr ( is_embedded )
            {
                // Provisional diffusion estimator. A dedicated embedded RKL2
                // estimator is not available yet.
                eval_diffusion( tn + dt, u_diff, fd_tmp );

                auto& err_D = U[13];
                legendre::dynamic::approximate_rkl2_diffusion_error<state_t, value_t>( fd_start, fd_tmp, dt, err_D );

                err_D_scalar = normalized_error_squared( err_D, un, u_diff );

                // Reject early when the diffusion estimate already exceeds the tolerance.
                if ( err_D_scalar >= static_cast<value_t>( 1. ) )
                {
                    _info.error   = err_D_scalar;
                    _info.success = false;

                    value_t fac = raw_fac_from_current_error();
                    fac         = finalize_fac( fac );

                    value_t const new_dt = dt * fac;
                    u_np1                = un;
                    dt                   = finalize_rejection_dt( new_dt );
                    return;
                }
            }

            // -----------------------------------------------------------------
            // 3. PIROCK diffusion-advection finishing, with beta=0.
            // -----------------------------------------------------------------
            eval_diffusion( t_diffusion_coupling, K, fd_tmp );
            eval_advection( t_advection_K, K, fa_K );

            u_sp3 = K + ( static_cast<value_t>( 1. ) - static_cast<value_t>( 2. ) * gamma ) * dt * fa_K;
            u_sp4 = K + dt / static_cast<value_t>( 3. ) * fa_K;

            eval_advection( t_advection_sp4, u_sp4, fa_tmp );

            // beta=0 exactly: there is no explicit F_D(K) contribution in sp5.
            u_sp5 = K + static_cast<value_t>( 2. ) / static_cast<value_t>( 3. ) * dt * fa_tmp;

            eval_advection( t_advection_sp5, u_sp5, fa_tmp_bis );

            // Final diffusion-advection coupling.
            eval_diffusion( t_diffusion_coupling, u_sp3, y_j );
            y_j = static_cast<state_t>( y_j - fd_tmp );

            u_np1 = u_diff + static_cast<value_t>( 0.25 ) * dt * fa_K
                  + dt / ( static_cast<value_t>( 2. ) - static_cast<value_t>( 4. ) * gamma ) * y_j
                  + static_cast<value_t>( 0.75 ) * dt * fa_tmp_bis;

            if constexpr ( is_embedded )
            {
                // Advection error estimate.
                auto& err_A = u_sp4;
                err_A       = -static_cast<value_t>( 0.15 ) * dt * fa_K + static_cast<value_t>( 0.3 ) * dt * fa_tmp
                      - static_cast<value_t>( 0.15 ) * dt * fa_tmp_bis;

                value_t const err_A_scalar = normalized_error_squared( err_A, un, u_np1 );

                // The advection estimator is third order. In RMS^2 form,
                // ||err_A||^(2/3) becomes err_A_scalar^(2/3).
                value_t const err_A_effective_scalar = std::pow( err_A_scalar, static_cast<value_t>( 2. / 3. ) );

                _info.error   = std::max( err_D_scalar, err_A_effective_scalar );
                _info.success = _info.error < static_cast<value_t>( 1. );

                value_t fac = raw_fac_from_current_error();

                // Controller memory correction when diffusion dominates.
                if ( errp != static_cast<value_t>( 0. ) && !previous_attempt_rejected
                     && err_D_scalar >= static_cast<value_t>( 100. ) * err_A_effective_scalar )
                {
                    value_t const facp = std::pow( errp, static_cast<value_t>( 0.25 ) ) * fac * fac * ( dt / hp );
                    fac                = std::min( fac, facp );
                }

                fac            = finalize_fac( fac );
                value_t new_dt = dt * fac;

                if ( _info.success )
                {
                    tn += dt;

                    errp = _info.error;
                    hp   = dt;

                    if ( previous_attempt_rejected )
                    {
                        new_dt = std::min( new_dt, dt );
                        nrej   = 0;
                    }

                    dt = new_dt;

                    facmax                    = static_cast<value_t>( 2. );
                    previous_attempt_rejected = false;
                }
                else
                {
                    std::swap( un, u_np1 );
                    dt = finalize_rejection_dt( new_dt );
                }
            }
            else
            {
                tn += dt;
            }
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

        /**
         * @brief Last PIRKL-D dimensioning data.
         *
         * Exposes the effective stability range, damping parameter and branch/
         * closure stage counts used during the last attempted step. This is
         * useful for diagnostics, profiling and cost analysis.
         */
        auto const&
        dimensioning() const
        {
            return last_dimensioning;
        }

        /** @brief Last externally supplied diffusion spectral radius. */
        value_t
        spectral_radius() const
        {
            return last_rho_D;
        }

        /** @brief Set the stability safety factor (default: 1.1). */
        auto&
        safety( value_t safety_ )
        {
            stability_safety = safety_;
            return *this;
        }

        /** @brief Set the damping constant K (default: 10). */
        auto&
        damping_constant( value_t K_ )
        {
            damping_K = K_;
            return *this;
        }

        template <bool embedded = is_embedded>
            requires embedded
        auto&
        abs_tol( value_t tol_ )
        {
            info().absolute_tolerance = tol_;
            return *this;
        }

        template <bool embedded = is_embedded>
            requires embedded
        auto&
        rel_tol( value_t tol_ )
        {
            info().relative_tolerance = tol_;
            return *this;
        }
    };

    /**
     * @brief Build a PIRKL-D method for diffusion-advection problems.
     *
     * A diffusion spectral-radius provider must be supplied.
     */
    template <bool is_embedded = false, typename value_t = double, typename eig_computer_t>
    auto
    pirkl_DA( eig_computer_t&& eig_computer )
    {
        return pirkl_DA_impl<eig_computer_t, is_embedded, value_t>( std::forward<eig_computer_t>( eig_computer ) );
    }

} // namespace ponio::runge_kutta::pirkl
