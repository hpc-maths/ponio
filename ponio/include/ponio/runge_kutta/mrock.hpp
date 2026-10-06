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

#include "../iteration_info.hpp"
#include "../stage.hpp"
#include "averaged_forces.hpp"
#include "rkc.hpp"
#include "rock.hpp"

namespace ponio::runge_kutta::mrock
{
    namespace detail
    {
        template <typename value_t>
        struct mrock2_dimensioning
        {
            std::size_t s        = 3;
            std::size_t mdeg     = 1;
            std::size_t m        = 2;
            value_t eta          = static_cast<value_t>( 0. );
            value_t alpha_m      = static_cast<value_t>( 0. );
            value_t beta_micro   = static_cast<value_t>( 0. );
            value_t beta_macro   = static_cast<value_t>( 0.811 );
            value_t rho_macro    = static_cast<value_t>( 0. );
            value_t effective_dt = static_cast<value_t>( 0. );
            bool capped          = false;
        };

        /**
         * @brief Compute mROCK2 parameters from Eq. (4.5) (general) or
         *        Eq. (4.6) (relaxed scale-separation branch).
         *
         * The mROCK2 paper fixes eps=0.05 for the practical 1.35 bound.
         */
        template <typename value_t>
        mrock2_dimensioning<value_t>
        get_mrock2_stages( value_t dt, value_t rho_fast, value_t rho_slow, ::ponio::runge_kutta::chebyshev::multirate_dimensioning mode )
        {
            if ( rho_fast < static_cast<value_t>( 0. ) || rho_slow < static_cast<value_t>( 0. ) )
            {
                throw std::invalid_argument( "mROCK2: spectral radii must be non-negative." );
            }

            value_t constexpr eps        = static_cast<value_t>( 0.05 );
            value_t constexpr beta_macro = static_cast<value_t>( 0.811 );
            value_t const beta_micro     = ::ponio::runge_kutta::chebyshev::rkc_detail::beta( eps );
            value_t const rho_macro      = static_cast<value_t>( 1.35 ) * rho_slow;

            value_t tau             = std::abs( dt );
            bool capped             = false;
            std::size_t preliminary = static_cast<std::size_t>(
                std::sqrt( ( static_cast<value_t>( 1.5 ) + tau * rho_macro ) / beta_macro ) + static_cast<value_t>( 1. ) );
            preliminary = std::max<std::size_t>( preliminary, 3 );

            if ( preliminary > 200 )
            {
                preliminary = 200;
                if ( rho_macro > static_cast<value_t>( 0. ) )
                {
                    tau = static_cast<value_t>( 0.8 )
                        * ( static_cast<value_t>( preliminary * preliminary ) * beta_macro - static_cast<value_t>( 1.5 ) ) / rho_macro;
                }
                capped = true;
            }

            std::size_t mdeg = preliminary - 2;
            using rock_coeff = ::ponio::runge_kutta::rock::rock2_coeff<value_t>;
            ::ponio::runge_kutta::rock::detail::degree_computer<value_t, rock_coeff>::optimal_degree( mdeg );
            std::size_t const s = mdeg + 2;
            value_t const ss    = static_cast<value_t>( s ) * static_cast<value_t>( s );

            std::size_t m = 1;
            value_t eta   = static_cast<value_t>( 0. );
            if ( mode == ::ponio::runge_kutta::chebyshev::multirate_dimensioning::general )
            {
                value_t const lower = static_cast<value_t>( 1. )
                                    + static_cast<value_t>( 6. ) * tau * rho_fast / ( beta_macro * beta_micro * ss );
                m            = std::max<std::size_t>( 2, static_cast<std::size_t>( std::floor( std::sqrt( lower ) ) ) + 1 );
                auto eta_for = [&]( std::size_t micro_stages )
                {
                    value_t const mm = static_cast<value_t>( micro_stages ) * static_cast<value_t>( micro_stages );
                    return static_cast<value_t>( 6. ) * tau / ( beta_macro * ss ) * mm / ( mm - static_cast<value_t>( 1. ) );
                };
                eta = eta_for( m );
                while ( eta * rho_fast > beta_micro * static_cast<value_t>( m ) * static_cast<value_t>( m ) )
                {
                    ++m;
                    eta = eta_for( m );
                }
            }
            else
            {
                // Eq. (4.6): smaller micro interval in the scale-separation regime.
                eta = static_cast<value_t>( 2.8 ) * tau / ( beta_macro * ss );
                m   = std::max<std::size_t>( 1, static_cast<std::size_t>( std::ceil( std::sqrt( eta * rho_fast / beta_micro ) ) ) );
            }

            value_t const alpha_m = ::ponio::runge_kutta::chebyshev::multirate_detail::rkc1_alpha_m( m, eps );
            return { s, mdeg, m, eta, alpha_m, beta_micro, beta_macro, rho_macro, tau, capped };
        }
    } // namespace detail

    /**
     * @brief Second-order explicit stabilized multirate ROCK2 method.
     */
    template <bool _is_embedded = false, typename _value_t = double>
    struct mrock2_impl
    {
        static constexpr bool is_embedded        = _is_embedded;
        static constexpr bool is_imex_method     = true;
        static constexpr std::size_t N_operators = 2;
        static constexpr std::size_t N_stages    = stages::dynamic;
        static constexpr std::size_t N_storage   = 13;
        static constexpr std::size_t order       = 2;
        static constexpr std::string_view id     = "mROCK2";
        static constexpr _value_t epsilon        = static_cast<_value_t>( 0.05 );

        using value_t             = _value_t;
        using dimensioning_t      = detail::mrock2_dimensioning<value_t>;
        using dimensioning_mode_t = ::ponio::runge_kutta::chebyshev::multirate_dimensioning;

        struct runtime_stats
        {
            std::size_t attempts = 0, accepted = 0, rejected = 0, capped_steps = 0;
            std::size_t s_min = std::numeric_limits<std::size_t>::max(), s_max = 0;
            std::size_t mdeg_min = std::numeric_limits<std::size_t>::max(), mdeg_max = 0;
            std::size_t m_min = std::numeric_limits<std::size_t>::max(), m_max = 0;
            long double s_sum = 0., mdeg_sum = 0., m_sum = 0.;
            value_t eta_min = std::numeric_limits<value_t>::infinity(), eta_max = 0.;
            value_t alpha_min = std::numeric_limits<value_t>::infinity(), alpha_max = 0.;
            value_t dt_min = std::numeric_limits<value_t>::infinity(), dt_max = 0.;
            long double eta_sum = 0., alpha_sum = 0., dt_sum = 0.;
            std::size_t fast_evals = 0, slow_evals = 0;

            value_t
            s_mean() const
            {
                return attempts ? static_cast<value_t>( s_sum / attempts ) : 0.;
            }

            value_t
            mdeg_mean() const
            {
                return attempts ? static_cast<value_t>( mdeg_sum / attempts ) : 0.;
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
            alpha_mean() const
            {
                return attempts ? static_cast<value_t>( alpha_sum / attempts ) : 0.;
            }

            value_t
            dt_mean() const
            {
                return attempts ? static_cast<value_t>( dt_sum / attempts ) : 0.;
            }
        };

        struct spectral_radii
        {
            value_t fast = static_cast<value_t>( 0. );
            value_t slow = static_cast<value_t>( 0. );
        };

        dimensioning_mode_t dimensioning_mode;
        iteration_info<mrock2_impl> _info;
        std::shared_ptr<runtime_stats> _stats;
        std::shared_ptr<spectral_radii> _radii;
        std::size_t last_s = 3, last_mdeg = 1, last_m = 2;
        value_t last_eta = 0., last_alpha = 0.;

        mrock2_impl( value_t fast, value_t slow, dimensioning_mode_t mode = dimensioning_mode_t::general )
            : dimensioning_mode( mode )
            , _info()
            , _stats( std::make_shared<runtime_stats>() )
            , _radii( std::make_shared<spectral_radii>( spectral_radii{ fast, slow } ) )
        {
            if ( fast < static_cast<value_t>( 0. ) || slow < static_cast<value_t>( 0. ) )
            {
                throw std::invalid_argument( "mROCK2: spectral radii must be non-negative." );
            }
        }

        /**
         * @brief Update fast/slow spectral-radius estimates.
         *
         * The estimates are shared between copies of the algorithm. Therefore
         * an application can update them after an MR adaptation and the copy
         * stored inside solver_range observes the new values at the next step.
         */
        void
        set_spectral_radii( value_t fast, value_t slow )
        {
            if ( fast < static_cast<value_t>( 0. ) || slow < static_cast<value_t>( 0. ) )
            {
                throw std::invalid_argument( "mROCK2: spectral radii must be non-negative." );
            }

            _radii->fast = fast;
            _radii->slow = slow;
        }

        value_t
        fast_spectral_radius() const
        {
            return _radii->fast;
        }

        value_t
        slow_spectral_radius() const
        {
            return _radii->slow;
        }

        dimensioning_t
        parameters( value_t dt ) const
        {
            return detail::get_mrock2_stages( dt, _radii->fast, _radii->slow, dimensioning_mode );
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
        last_macro_degree() const
        {
            return last_mdeg;
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
        last_alpha_m() const
        {
            return last_alpha;
        }

        dimensioning_mode_t
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
        update_stats( dimensioning_t const& p, value_t step_dt, bool success, std::size_t fast_evals, std::size_t slow_evals )
        {
            auto& st = *_stats;
            ++st.attempts;
            success ? ++st.accepted : ++st.rejected;
            if ( p.capped )
            {
                ++st.capped_steps;
            }

            st.s_min = std::min( st.s_min, p.s );
            st.s_max = std::max( st.s_max, p.s );
            st.s_sum += static_cast<long double>( p.s );
            st.mdeg_min = std::min( st.mdeg_min, p.mdeg );
            st.mdeg_max = std::max( st.mdeg_max, p.mdeg );
            st.mdeg_sum += static_cast<long double>( p.mdeg );
            st.m_min = std::min( st.m_min, p.m );
            st.m_max = std::max( st.m_max, p.m );
            st.m_sum += static_cast<long double>( p.m );
            st.eta_min = std::min( st.eta_min, p.eta );
            st.eta_max = std::max( st.eta_max, p.eta );
            st.eta_sum += p.eta;
            st.alpha_min = std::min( st.alpha_min, p.alpha_m );
            st.alpha_max = std::max( st.alpha_max, p.alpha_m );
            st.alpha_sum += p.alpha_m;
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
            auto const params     = parameters( dt );
            dt                    = params.effective_dt;
            value_t const step_dt = dt;
            last_s                = params.s;
            last_mdeg             = params.mdeg;
            last_m                = params.m;
            last_eta              = params.eta;
            last_alpha            = params.alpha_m;

            auto averaged_rhs = [&]( value_t t, state_t& y, state_t& dy )
            {
                ::ponio::runge_kutta::chebyshev::multirate_detail::averaged_force_order_2( pb,
                    t,
                    y,
                    U,
                    params.m,
                    params.eta,
                    params.alpha_m,
                    epsilon,
                    dy );
            };
            auto macro_eig = [rho = params.rho_macro]( auto&&, value_t, state_t&, value_t, auto& )
            {
                return rho;
            };
            auto macro = ::ponio::runge_kutta::rock::rock2<is_embedded, value_t>( std::move( macro_eig ) );
            if constexpr ( is_embedded )
            {
                macro.abs_tol( _info.absolute_tolerance ).rel_tol( _info.relative_tolerance );
            }
            macro( averaged_rhs, tn, un, U, dt, unp1 );

            _info.number_of_stages           = macro.info().number_of_stages;
            _info.error                      = macro.info().error;
            _info.success                    = macro.info().success;
            std::size_t const averaged_calls = macro.info().number_of_eval;
            _info.number_of_eval[0]          = static_cast<std::size_t>( 2 ) * params.m * averaged_calls;
            _info.number_of_eval[1]          = averaged_calls;
            update_stats( params, step_dt, _info.success, _info.number_of_eval[0], _info.number_of_eval[1] );
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
    mrock2( value_t rho_fast,
        value_t rho_slow,
        ::ponio::runge_kutta::chebyshev::multirate_dimensioning mode = ::ponio::runge_kutta::chebyshev::multirate_dimensioning::general )
    {
        return mrock2_impl<is_embedded, value_t>( rho_fast, rho_slow, mode );
    }

} // namespace ponio::runge_kutta::mrock
