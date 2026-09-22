// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private

#pragma once

// NOLINTBEGIN(misc-include-cleaner)

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <numbers>
#include <numeric>
#include <string_view>
#include <tuple>
#include <type_traits>

#include "../detail.hpp"
#include "../iteration_info.hpp"
#include "../linear_algebra.hpp"
#include "../ponio_config.hpp"
#include "../stage.hpp"
#include "dirk.hpp"
#include "rock.hpp"
#include "rock_coeff.hpp"

// NOLINTEND(misc-include-cleaner)

namespace ponio::runge_kutta::pirock
{
    namespace polynomial
    {
        /**
         * @brief Computes \f$P'_{s-2+\ell}(0)\f$ from number of stages \f$s\f$ and \f$\ell\f$
         *
         * @tparam value_t type of coefficients
         * @param s        number of stages
         * @param l        number of additional stages
         */
        template <typename value_t = double>
        value_t
        Pp_sm2pl_0( std::size_t s, std::size_t l )
        {
            using rock_coeff        = rock::rock2_coeff<value_t>;
            std::size_t const sm2pl = s - 2 + l;

            std::array<value_t, 2> u    = { 0., 0. };
            std::array<value_t, 2> ujm1 = { 0., 0. };
            std::array<value_t, 2> ujm2 = { 0., 0. };

            std::size_t mdeg = s - 2;
            auto [mz, mr]    = rock::detail::degree_computer<value_t, rock_coeff>::optimal_degree( mdeg );

            u[0]    = 1.0;
            ujm2[0] = 1.0;

            auto mu_1 = rock_coeff::recf[mr - 1];

            // in our case f is function that multiply by monome X
            // f only shift degree of its argument
            ujm1[0] = 1.0;  // u_jm1 = un = 1
            ujm1[1] = mu_1; // u_jm1 += mu_1 * f( un )

            if ( mdeg <= 2 )
            {
                u = ujm1;
            }

            for ( std::size_t j = 2; j < sm2pl + 1; ++j )
            {
                value_t mu;
                value_t kappa;

                if ( j <= mdeg )
                {
                    mu    = rock_coeff::recf[mr + 2 * ( j - 2 ) + 1 - 1];
                    kappa = rock_coeff::recf[mr + 2 * ( j - 2 ) + 2 - 1];
                }
                else
                {
                    std::size_t const additional_stage = j - mdeg;
                    std::size_t const recf2_index      = 4 * ( mz - 1 ) + 2 * ( additional_stage - 1 );
                    mu                                 = rock_coeff::recf2[recf2_index];
                    kappa                              = rock_coeff::recf2[recf2_index + 1];
                }

                value_t const nu = -1.0 - kappa;

                // u_j = mu * f(u_jm1) - nu * u_jm1 - kappa * u_jm2
                u[0] = -nu * ujm1[0] - kappa * ujm2[0];
                u[1] = mu * ujm1[0] - nu * ujm1[1] - kappa * ujm2[1];

                ujm2 = ujm1;
                ujm1 = u;
            }

            // we only want the derivative of the polynomial and evaluate it in 0
            return u[1];
        }
    } // namespace polynomial

    /**
     * @class alpha_fixed
     * @brief Computer of \f$\alpha\f$ and \f$\beta\f$ parameters with fixed value of \f$\alpha\f$
     *
     * @tparam value_t Type of coefficients
     */
    template <typename value_t = double>
    struct alpha_fixed
    {
        value_t _alpha;

        alpha_fixed( value_t a = static_cast<value_t>( 1. ) )
            : _alpha( a )
        {
        }

        /**
         * @brief Return \f$\alpha\f$ (given value in constructor)
         */
        value_t
        alpha( std::size_t, std::size_t ) const
        {
            return _alpha;
        }

        /**
         * @brief Return \f$\beta = 1 - 2\alpha P'_{s-2+\ell}(0)\f$
         *
         * @param s number of stages
         * @param l choosen parameter \f$\ell\f$
         */
        value_t
        beta( std::size_t s, std::size_t l ) const
        {
            return 1. - 2. * alpha( s, l ) * polynomial::Pp_sm2pl_0<value_t>( s, l );
        }
    };

    /**
     * @class beta_0
     * @brief Computer of \f$\alpha\f$ and \f$\beta\f$ parameters with fixed value of \f$\beta\f$ to 0
     *
     * @tparam value_t Type of coefficients
     */
    template <typename value_t = double>
    struct beta_0
    {
        /**
         * @brief Return \f$\alpha = \frac{1}{2P'_{s-2+\ell}(0)}\f$ such \f$\beta = 0\f$
         */
        value_t
        alpha( std::size_t s, std::size_t l ) const
        {
            return 1. / ( 2. * polynomial::Pp_sm2pl_0<value_t>( s, l ) );
        }

        /**
         * @brief Return \f$\beta = 0\f$
         */
        value_t
        beta( std::size_t, std::size_t ) const
        {
            return 0.;
        }
    };

    // --- PIROCK REACTION-DIFFUSION ------------------------------------------

    /**
     * @class pirock_impl
     * @brief implementation of PIROCK method
     *
     * @tparam _l                       Number of augmented stages (l = 1, 2)
     * @tparam alpha_beta_computer_t    Choice of computing of \f$\alpha\f$ and \f$\beta\f$ parameters
     * @tparam eig_computer_t           Computing of eigenvalues of explicit part (diffusion)
     * @tparam _shampine_trick_caller_t Computing of Shampine's trick
     * @tparam _is_embedded             Set adaptive time step method (default: false)
     * @tparam value_t                  Type of coefficients
     */
    template <std::size_t _l,
        typename alpha_beta_computer_t,
        typename eig_computer_t,
        typename _shampine_trick_caller_t = void,
        bool _is_embedded                 = false,
        typename _value_t                 = double>
    struct pirock_impl
    {
        static constexpr bool is_embedded           = _is_embedded;
        static constexpr bool shampine_trick_enable = !std::is_void_v<_shampine_trick_caller_t>;

        static constexpr std::size_t l       = _l;
        static constexpr bool is_imex_method = true;
        // number_of_eval counts one evaluation per operator: 0 diffusion, 1 reaction.
        static constexpr std::size_t N_operators = 2;
        static constexpr std::size_t N_stages    = stages::dynamic;
        // clang-format off
        static constexpr std::size_t N_storage   = std::conditional_t<shampine_trick_enable,
                                                    std::integral_constant<std::size_t, 21>, // if shampine's trick
                                                    std::integral_constant<std::size_t, 15>  // else (no error estimation)
                                                >::value;
        // clang-format on
        static constexpr std::array<std::size_t, 1> persistent_storage_indices = { N_storage - 2 };
        static constexpr std::size_t order                                     = 2;
        static constexpr std::string_view id                                   = "PIROCK";

        using value_t                 = _value_t;
        using rock_coeff              = rock::rock2_coeff<value_t>;
        using degree_computer         = rock::detail::degree_computer<value_t, rock_coeff>;
        using shampine_trick_caller_t = typename std::conditional_t<shampine_trick_enable, _shampine_trick_caller_t, bool>;

        alpha_beta_computer_t alpha_beta_computer;
        eig_computer_t eig_computer;
        shampine_trick_caller_t shampine_trick_caller;

        iteration_info<pirock_impl> _info;

        bool compensated_summation_initialized = false;
        value_t time_compensation              = static_cast<value_t>( 0. );

        // Adaptive controller state. These values persist across step attempts,
        // as in the reference Fortran `rockcore` controller.
        value_t facmax                 = static_cast<value_t>( 5. );  // max growth: 5 initially, 1 after rejection, 2 after acceptance
        value_t facd                   = static_cast<value_t>( 0.8 ); // adaptive safety factor, updated every 100 attempts
        std::size_t attempt_count      = 0;
        std::size_t nrejfac            = 0;
        bool previous_attempt_rejected = false;

        // Fortran controller with memory (`iwork(19) = 2`): error and step
        // size of the last accepted step.
        value_t errp = static_cast<value_t>( 0. );
        value_t hp   = static_cast<value_t>( 0. );

        // Consecutive-rejection safeguard from the Fortran controller.
        value_t told     = static_cast<value_t>( 0. );
        std::size_t nrej = 0;

        /**
         * @brief Construct a new pirock impl object
         *
         */
        pirock_impl() = default;

        /**
         * @brief Construct a new pirock impl object
         *
         * @param _alpha_beta_computer computer of parameters alpha and beta object that has two member functions which take number of
         * stages (s) l parameter
         * @param _eig_computer        eigenvalue computer functor (that take as argument the function of explicit part, the current time,
         * the current state and the current time step)
         */
        pirock_impl( alpha_beta_computer_t&& _alpha_beta_computer, eig_computer_t&& _eig_computer )
            : alpha_beta_computer( std::forward<alpha_beta_computer_t>( _alpha_beta_computer ) )
            , eig_computer( std::forward<eig_computer_t>( _eig_computer ) )
            , shampine_trick_caller( false )
            , _info()
        {
        }

        /**
         * @brief Construct a new pirock impl object with Shampine's trick
         *
         * @param _alpha_beta_computer   alpha and beta computer object
         * @param _eig_computer          eigenvalue computer functor
         * @param _shampine_trick_caller Shampine's trick functor
         */
        template <typename _shampine_trick_caller_t_>
            requires std::same_as<_shampine_trick_caller_t_, shampine_trick_caller_t>
                      && std::same_as<std::bool_constant<shampine_trick_enable>, std::true_type>
        pirock_impl( alpha_beta_computer_t&& _alpha_beta_computer,
            eig_computer_t&& _eig_computer,
            _shampine_trick_caller_t_&& _shampine_trick_caller )
            : alpha_beta_computer( std::forward<alpha_beta_computer_t>( _alpha_beta_computer ) )
            , eig_computer( std::forward<eig_computer_t>( _eig_computer ) )
            , shampine_trick_caller( std::forward<_shampine_trick_caller_t_>( _shampine_trick_caller ) )
            , _info()
        {
        }

        /**
         * @brief iteration of PIROCK method
         *
         * @tparam problem_t  type of \f$f\f$
         * @tparam state_t    type of current state
         * @tparam array_ki_t type of temporary stages (only 3 needed for ROCK2)
         * @param pb    problem \f$(F_D, F_R)\f$ and the Jacobian of reaction part \f$\frac{\partial F_R}{\partial u}\f$
         * @param tn    current time
         * @param un    current state
         * @param U     array of temporary stages
         * @param dt    current time step
         * @param u_np1 solution \f$u^{n+1}\f$ at time \f$t^{n+1} = t^n + \Delta t\f$
         */
        template <typename problem_t, typename state_t, typename array_ki_t>
        void
        operator()( problem_t& pb, value_t& tn, state_t& un, array_ki_t& U, value_t& dt, state_t& u_np1 )
        {
            static_assert( detail::problem_operator<decltype( pb.implicit_part ), value_t>
                               || detail::problem_jacobian<decltype( pb.implicit_part ), value_t, state_t>,
                "This kind of problem is not inversible in ponio" );

            // U worker references:
            // | index | variables        | mathematic representation                         |
            // |-------|------------------|---------------------------------------------------|
            // | 0     | u_j              | $u^{(j)}$ current stage in pseudo-ROCK2 method    |
            // | 0     | fe_tmp_bis       | temporary output of explicit part                 |
            // | 1     | u_jm1            | $u^{(j-1)}$ previous stage in pseudo-ROCK2 method |
            // | 1     | fe_tmp_ter       | temporary output of explicit part                 |
            // | 2     | u_jm2            | $u^{(j-2)}$ previous stage in pseudo-ROCK2 method |
            // | 2     | fe_tmp_qua       | temporary output of explicit part                 |
            // | 3     | u_sm2            | $u^{(s-2)}$ stage in PIROCK method                |
            // | 4     | fe_tmp           | temporary output of explicit part (main output)   |
            // | 5     | fi_tmp           | temporary output of implicit part (main output)   |
            // | 6     | f_tmp            | temporary output of implicit part                 |
            // | 7     | us_sm1           | $u^{\star(s-1)}$ stage in PIROCK method           |
            // | 8     | us_s             | $u^{\star(s)}$ stage in PIROCK method             |
            // | 9     | u_sp1            | $u^{(s+1)}$ stage in PIROCK method                |
            // | 10    | u_sp2            | $u^{(s+2)}$ stage in PIROCK method                |
            // | 11    | u_sp3            | $u^{(s+3)}$ stage in PIROCK method                |
            // | 11    | rhs_sp2          | right-hand side to compute $u^{(s+2)}$ stage      |
            // | 12    | shampine_element | Shampine's trick temporary element                |
            // | 13    | f_D_u            | $F_D(u^{(s+3)}) - F_D(u^{(s+1)})$ term             |
            // | 14    | u_tmp            | temporary value for compute Shampine's trick      |
            // | 15    | fd_tn_cache      | cached $F_D(t_n,u_n)$ across rejected retries     |
            // | 16    | rhs_R            | right-hand side to compute $err_R$ error term     |
            // | 17    | err_R            | $err_R$ error on reaction term                    |
            // | last-2| rock_increment    | transported ROCK increment; then $err_D$           |
            // | last-1| compensation      | compensated summation remainder                    |
            // | last  | compensation_tmp  | temporary state for compensated summation          |
            //
            // > if method is called as a constant time step method, only index from 0 to 11 are used
            // > in addition to the three compensated summation work arrays

            _info.reset_eval();

            // PIROCK uses the Fortran ROCK2 stability constants 0.432 for
            // ell=1 and 0.811 for ell=2.
            value_t const rock2_stability_constant = ( l == 1 ) ? static_cast<value_t>( 0.432 ) : static_cast<value_t>( 0.811 );

            auto [mdeg, deg_index, start_index, n_eval] = degree_computer::compute_n_stages_optimal_degree( rock::rock_order::rock_2(),
                eig_computer,
                pb.explicit_part,
                tn,
                un,
                dt,
                U,
                4,
                rock2_stability_constant );

            std::size_t s = mdeg + 2;

            _info.number_of_stages = s + l + 3;
            // Count the spectral-radius evaluations and ROCK2 diffusion stages.
            // On an adaptive retry, F_D(t_n,u_n) is reused from U[15], so one
            // ROCK diffusion evaluation is intentionally skipped.
            _info.number_of_eval[0] = n_eval + s + l;
            if constexpr ( is_embedded )
            {
                if ( previous_attempt_rejected )
                {
                    --_info.number_of_eval[0];
                }
            }

            value_t const alpha = alpha_beta_computer.alpha( s, l );
            value_t const beta  = alpha_beta_computer.beta( s, l );
            value_t const gamma = 1. - 0.5 * std::numbers::sqrt2;

            value_t const t_sp2 = tn + ( 1. - gamma ) * dt;
            // Internal diffusion abscissa c_{s-2+l}; the Fortran code names
            // this quantity `beta`, independently of the PIROCK coefficient beta.
            value_t const t_diffusion_coupling = tn + 0.5 * ( 1. - beta ) * dt;

            auto& u_j    = U[0];
            auto& u_jm1  = U[1];
            auto& u_jm2  = U[2];
            auto& u_sm2  = U[3];
            auto& fe_tmp = U[4];
            auto& fi_tmp = U[5];
            auto& f_tmp  = U[6];

            auto& rock_increment   = U[N_storage - 3];
            auto& compensation     = U[N_storage - 2];
            auto& compensation_tmp = U[N_storage - 1];

            // The compensated sums below update a state in place. state_algebra
            // carries that operation for the state type at hand: compound
            // assignment where the type provides it, as Eigen vectors do, and
            // assignment from an expression where it does not, as for samurai
            // fields. The accumulation itself reads the same in both cases.
            using state_algebra_t = ::ponio::linear_algebra::state_algebra<state_t>;

            if ( !compensated_summation_initialized )
            {
                compensation = un;
                state_algebra_t::scale( compensation, value_t( 0 ) );
                compensated_summation_initialized = true;
            }

            // Compensated summation keeps small stage increments that would
            // otherwise be lost to round-off.
            auto compensated_update = [&]( state_t& state, auto const& increment )
            {
                state_algebra_t::add( compensation, increment );
                compensation_tmp = state;
                state            = state + compensation;
                state_algebra_t::add( compensation, compensation_tmp - state );
            };

            auto compensated_commit = [&]( state_t const& state, state_t& result )
            {
                compensation_tmp = state;
                result           = state + compensation;
                state_algebra_t::add( compensation, compensation_tmp - result );
            };

            auto compensated_time_update = [&]( value_t increment )
            {
                value_t const previous_time = tn;
                time_compensation += increment;
                tn = tn + time_compensation;
                time_compensation += previous_time - tn;
            };

            // Adaptive controller shared by the early-rejection and full-step paths.
            auto raw_fac_from_current_error = [&]() -> value_t
            {
                ++attempt_count;

                // Fortran adapts facd every 100 attempts from the recent rejection rate.
                if ( attempt_count % 100 == 0 )
                {
                    if ( nrejfac >= 10 )
                    {
                        facd *= 0.8;
                    }
                    if ( nrejfac == 0 )
                    {
                        facd *= 1.02;
                    }
                    facd    = std::min( 0.8, std::max( 0.4, facd ) );
                    nrejfac = 0;
                }

                // _info.error stores RMS^2, hence the -1/4 exponent required
                // to reproduce Fortran's fac = RMS^(-1/2).
                return std::pow( _info.error, -0.25 );
            };

            auto finalize_fac = [&]( value_t fac ) -> value_t
            {
                if ( previous_attempt_rejected )
                {
                    facmax = 1.0;
                }

                // Fortran bounds the damped factor facd*fac, not fac alone.
                return std::min( facmax, std::max( 0.1, facd * fac ) );
            };

            auto finalize_rejection_dt = [&]( value_t new_dt ) -> value_t
            {
                // Fortran applies an additional 0.8 reduction after rejection.
                value_t retry_dt = 0.8 * new_dt;

                if ( told == tn )
                {
                    ++nrej;
                    if ( nrej == 10 )
                    {
                        retry_dt = 1.0e-5;
                    }
                }
                told = tn;

                ++nrejfac;
                previous_attempt_rejected = true;

                return retry_dt;
            };

            u_j   = un;
            u_jm2 = un;

            value_t const mu_1 = rock_coeff::recf[start_index - 1];

            value_t t_jm1 = tn + dt * alpha * mu_1;
            value_t t_jm2 = tn + dt * alpha * mu_1;
            value_t t_jm3 = tn;

            // u_1 =u^n + \alpha \mu_1 \Delta F_D( u^n )
            // A rejected adaptive attempt leaves (tn, un) unchanged. Cache
            // F_D(tn, un) and reuse it directly on retries, as Fortran keeps `fn`
            // across a rejected attempt. U[15] is otherwise unused here, so this
            // requires no extra allocation or vector copy.
            if constexpr ( is_embedded )
            {
                auto& fd_tn_cache = U[15];
                if ( !previous_attempt_rejected )
                {
                    pb.explicit_part( tn, un, fd_tn_cache );
                }
                rock_increment = alpha * dt * mu_1 * fd_tn_cache;
            }
            else
            {
                pb.explicit_part( tn, un, fe_tmp );
                rock_increment = alpha * dt * mu_1 * fe_tmp;
            }
            u_jm1 = un;
            compensated_update( u_jm1, rock_increment );

            if ( mdeg < 2 )
            {
                u_j = u_jm1;
            }

            value_t t_sm2 = tn;

            for ( std::size_t j = 2; j < s - 2 + l + 1; ++j )
            {
                value_t mu_j;
                value_t kappa_j;

                if ( j <= mdeg )
                {
                    mu_j    = rock_coeff::recf[start_index + 2 * ( j - 2 ) + 1 - 1];
                    kappa_j = rock_coeff::recf[start_index + 2 * ( j - 2 ) + 2 - 1];
                }
                else
                {
                    std::size_t const additional_stage = j - mdeg;
                    std::size_t const recf2_index      = 4 * ( deg_index - 1 ) + 2 * ( additional_stage - 1 );
                    mu_j                               = rock_coeff::recf2[recf2_index];
                    kappa_j                            = rock_coeff::recf2[recf2_index + 1];
                }

                value_t const nu_j = -1.0 - kappa_j;

                // u_{j} = \alpha \mu_j \Delta t F_D( u_{j-2} ) - \nu_j u_{j-1} - \kappa_j u_{j-2}
                pb.explicit_part( t_jm1, u_jm1, fe_tmp );

                if ( j <= mdeg )
                {
                    state_algebra_t::scale( rock_increment, kappa_j );
                    state_algebra_t::add( rock_increment, alpha * mu_j * dt * fe_tmp );
                    u_j = u_jm1;
                    compensated_update( u_j, rock_increment );
                }
                else
                {
                    // Fortran `rtstep` applies the recf2 stages directly to the
                    // stage values, without compensated summation.
                    u_j = alpha * mu_j * dt * fe_tmp - nu_j * u_jm1 - kappa_j * u_jm2;
                }

                t_jm1 = alpha * dt * mu_j - nu_j * t_jm2 - kappa_j * t_jm3;

                if ( j == s - 2 )
                {
                    u_sm2 = u_j;
                    t_sm2 = t_jm1;
                }
                if ( j < s - 2 + l )
                {
                    std::swap( u_jm2, u_jm1 );
                    std::swap( u_jm1, u_j );
                }

                t_jm3 = t_jm2;
                t_jm2 = t_jm1;
            }

            // if l == 1
            // u_j -> u_{s-2+l} = u_{s-1}
            // u_jm1 -> u_{s-2+l-1} = u_{s-2}
            // u_jm2 -> u_{s-2+l-2} = u_{s-3}

            // if l == 2
            // u_j -> u_{s-2+l} = u_{s}
            // u_jm1 -> u_{s-2+l-1} = u_{s-1}
            // u_jm2 -> u_{s-2+l-2} = u_{s-2}

            value_t sigma   = rock_coeff::fp1[deg_index - 1];
            value_t sigma_a = 0.5 * ( 1.0 - alpha ) + alpha * sigma;
            value_t tau     = sigma * rock_coeff::fp2[deg_index - 1] + sigma * sigma;
            value_t tau_a   = 0.5 * detail::power<2>( alpha - 1. ) + 2. * alpha * ( 1. - alpha ) * sigma + alpha * alpha * tau;

            // The ROCK2 finishing procedure computes err_D once; the embedded
            // controller reuses the same value, as in Fortran `rtstep`.
            value_t err_D_scalar = static_cast<value_t>( 0. );

            // u_{*s-1} = u_{s-2} + \sigma_\alpha \Delta t  F_D( u_{s-2} )
            auto& us_sm1 = U[7];
            pb.explicit_part( t_sm2, u_sm2, fe_tmp );
            rock_increment = fe_tmp;
            us_sm1         = u_sm2;
            compensated_update( us_sm1, sigma_a * dt * fe_tmp );

            // u_{*s} = u_{*s-1} + \sigma_\alpha \Delta t F_D( u_{*s-1} ) - err_D
            // with
            // err_D = \sigma_\alpha(1-\tau_\alpha/\sigma_\alpha^2) \Delta t
            //         [F_D(u_{*s-1}) - F_D(u_{s-2})].
            // Reuse the ROCK increment workspace for err_D. The sign convention
            // matches the Fortran finishing procedure (`temp3 = -err_D`).
            auto& us_s  = U[8];
            auto& err_D = rock_increment;
            pb.explicit_part( t_sm2 + sigma_a * dt, us_sm1, fe_tmp );

            err_D = sigma_a * ( 1. - tau_a / ( sigma_a * sigma_a ) ) * dt * ( fe_tmp - err_D );

            us_s = us_sm1;
            compensated_update( us_s, sigma_a * dt * fe_tmp - err_D );

            if constexpr ( is_embedded )
            {
                // Fortran rejects before the implicit stages when the diffusion
                // estimator alone exceeds one. Normalize err_D with the
                // diffusion candidate us_s, as done in `rtstep`.
                err_D_scalar = detail::error_algebra<state_t>::estimate_squared( err_D,
                    un,
                    us_s,
                    _info.absolute_tolerance,
                    _info.relative_tolerance );

                if ( err_D_scalar >= 1.0 )
                {
                    // Skip the implicit stages when diffusion already rejects the step.

                    _info.error   = err_D_scalar;
                    _info.success = false;

                    value_t new_dt;

                    value_t fac = raw_fac_from_current_error();

                    // facp is unavailable here because err_R is not evaluated.
                    fac = finalize_fac( fac );

                    new_dt = dt * fac;

                    u_np1 = un;
                    dt    = finalize_rejection_dt( new_dt );

                    return;
                }
            }

            // u_{s-2+l} = u_j
            auto& u_sm2pl = u_j;

            auto& u_sp1 = U[9];
            u_sp1       = un;

            auto& u_sp2 = U[10];
            u_sp2       = un;

            if constexpr ( detail::problem_operator<decltype( pb.implicit_part ), value_t> )
            {
                std::size_t n_eval_sp1 = 0;

                auto op_sp1 = ::ponio::linear_algebra::operator_algebra<state_t>::identity( un ) - gamma * dt * pb.implicit_part.f_t( tn );

                auto rhs_sp1 = u_sm2pl;

                // The first implicit stage starts from u_{s-2+l}, the initial
                // guess the nonlinear iteration is given on the matrix path.
                u_sp1 = u_sm2pl;

                ::ponio::linear_algebra::operator_algebra<state_t>::solve( op_sp1, u_sp1, rhs_sp1, n_eval_sp1 );

                std::size_t n_eval_sp2 = 0;

                _info.number_of_eval[0] += 1;
                pb.explicit_part( t_diffusion_coupling, u_sp1, fe_tmp );
                pb.implicit_part( tn, u_sp1, fi_tmp );
                auto& rhs_sp2 = U[11];

                auto op_sp2 = ::ponio::linear_algebra::operator_algebra<state_t>::identity( un ) - gamma * dt * pb.implicit_part.f_t( t_sp2 );

                rhs_sp2 = u_sm2pl + beta * dt * fe_tmp + ( 1. - 2. * gamma ) * dt * fi_tmp;

                // The second implicit stage starts from its own right-hand side
                // c2, as it does on the matrix path.
                u_sp2 = rhs_sp2;

                ::ponio::linear_algebra::operator_algebra<state_t>::solve( op_sp2, u_sp2, rhs_sp2, n_eval_sp2 );

                // Operator-based solves do not retain F_R. Store F_R(u^{s+2})
                // explicitly for the remaining PIROCK formulas.
                pb.implicit_part( t_sp2, u_sp2, f_tmp );

                _info.number_of_eval[1] += n_eval_sp1 + n_eval_sp2 + 2;
            }
            else
            {
                using matrix_t                = std::decay_t<decltype( pb.implicit_part.df( tn, un ) )>;
                using matrix_linear_algebra_t = ::ponio::linear_algebra::linear_algebra<matrix_t>;

                constexpr bool can_reuse_factorization = requires( matrix_t const& matrix, state_t const& rhs ) {
                                                             matrix_linear_algebra_t::identity_minus( matrix, value_t{} );
                                                             matrix_linear_algebra_t::factorize( matrix );
                                                             matrix_linear_algebra_t::solve_factorized( rhs );
                                                         };

                auto g_sp1 = [&]( state_t& u ) -> state_t
                {
                    _info.number_of_eval[1] += 1;

                    pb.implicit_part( tn, u, fi_tmp );

                    return u - gamma * dt * fi_tmp - u_sm2pl;
                };

                if constexpr ( can_reuse_factorization )
                {
                    // Start the simplified Newton process from one frozen reaction
                    // Jacobian. Stage 2 reuses the last factorization left by stage 1,
                    // following the Fortran `ieuler` reuse policy.
                    decltype( auto ) reaction_jacobian = [&]() -> decltype( auto )
                    {
                        return pb.implicit_part.df( tn, u_sm2pl );
                    }();

                    matrix_t frozen_stage_matrix;
                    frozen_stage_matrix = matrix_linear_algebra_t::identity_minus( reaction_jacobian, gamma * dt );

                    matrix_linear_algebra_t::factorize( frozen_stage_matrix );

                    // Fortran `ieuler` reuses the factorization for corrections 1--5
                    // and refreshes it before correction 6 and each later correction.
                    auto refresh_stage1_factorization = [&]( state_t const& u )
                    {
                        decltype( auto ) refreshed_reaction_jacobian = [&]() -> decltype( auto )
                        {
                            return pb.implicit_part.df( tn, u );
                        }();

                        matrix_t refreshed_stage_matrix;
                        refreshed_stage_matrix = matrix_linear_algebra_t::identity_minus( refreshed_reaction_jacobian, gamma * dt );

                        matrix_linear_algebra_t::factorize( refreshed_stage_matrix );
                    };

                    u_sp1 = diagonal_implicit_runge_kutta::pirock_newton<value_t, state_t, matrix_t>( g_sp1,
                        u_sm2pl,
                        refresh_stage1_factorization );
                }
                else
                {
                    auto identity = matrix_linear_algebra_t::identity( un );

                    auto dg_sp1 = [&]( state_t& u ) -> matrix_t
                    {
                        return identity - gamma * dt * pb.implicit_part.df( tn, u );
                    };

                    u_sp1 = diagonal_implicit_runge_kutta::newton<value_t>( g_sp1,
                        dg_sp1,
                        u_sm2pl,
                        matrix_linear_algebra_t::solver,
                        ponio::default_config::newton_tolerance,
                        ponio::default_config::newton_max_iterations );
                }

                _info.number_of_eval[0] += 1;

                pb.explicit_part( t_diffusion_coupling, u_sp1, fe_tmp );

                // Keep the last F_R value produced by stage 1. Fortran `ieuler`
                // returns the same pre-final-correction value in `fnc`.

                // Right-hand side of the second implicit stage:
                // c2 = u_{s-2+l} + beta*dt*F_D(u_{s+1})
                //                    + (1-2*gamma)*dt*F_R(u_{s+1}).
                // U[11] is used as workspace here and overwritten later by u_sp3.
                auto& rhs_sp2 = U[11];
                rhs_sp2       = u_sm2pl + beta * dt * fe_tmp + ( 1. - 2. * gamma ) * dt * fi_tmp;

                auto g_sp2 = [&]( state_t& u ) -> state_t
                {
                    _info.number_of_eval[1] += 1;

                    pb.implicit_part( t_sp2, u, f_tmp );

                    return u - gamma * dt * f_tmp - rhs_sp2;
                };

                if constexpr ( can_reuse_factorization )
                {
                    // Stage 2 starts from the factorization left by stage 1, as in
                    // Fortran (`is_frjac = .false.`), then follows the same refresh rule.
                    auto refresh_stage2_factorization = [&]( state_t const& u )
                    {
                        decltype( auto ) refreshed_reaction_jacobian = [&]() -> decltype( auto )
                        {
                            return pb.implicit_part.df( t_sp2, u );
                        }();

                        matrix_t refreshed_stage_matrix;
                        refreshed_stage_matrix = matrix_linear_algebra_t::identity_minus( refreshed_reaction_jacobian, gamma * dt );

                        matrix_linear_algebra_t::factorize( refreshed_stage_matrix );
                    };

                    u_sp2 = diagonal_implicit_runge_kutta::pirock_newton<value_t, state_t, matrix_t>( g_sp2,
                        rhs_sp2,
                        refresh_stage2_factorization );
                }
                else
                {
                    auto identity = matrix_linear_algebra_t::identity( un );

                    auto dg_sp2 = [&]( state_t& u ) -> matrix_t
                    {
                        return identity - gamma * dt * pb.implicit_part.df( t_sp2, u );
                    };

                    u_sp2 = diagonal_implicit_runge_kutta::newton<value_t>( g_sp2,
                        dg_sp2,
                        rhs_sp2,
                        matrix_linear_algebra_t::solver,
                        ponio::default_config::newton_tolerance,
                        ponio::default_config::newton_max_iterations );
                }
            }

            auto& u_sp3 = U[11];

            // Reuse the stage-1 reaction value left in Fortran `fnc`.
            u_sp3 = u_sm2pl + ( 1. - gamma ) * dt * fi_tmp;

            if constexpr ( shampine_trick_enable )
            {
                auto& shampine_element = U[12];
                auto& f_D_u            = U[13];
                auto& u_tmp            = U[14];

                // Reuse err_D from the ROCK2 finishing procedure; Fortran carries
                // the same estimator from `rtstep` to the final controller.

                // fe_tmp still stores F_D(u^{s+1}); only F_D(u^{s+3}) is new.
                _info.number_of_eval[0] += 1;
                pb.explicit_part( t_diffusion_coupling, u_sp3, f_D_u );

                f_D_u = static_cast<state_t>( f_D_u - fe_tmp );

                if constexpr ( is_embedded )
                {
                    auto& rhs_R = U[16];
                    auto& err_R = U[17];

                    // Fortran `rkstep` builds err_R from the reaction values left
                    // by the two implicit stages.
                    rhs_R = static_cast<state_t>( dt / 6. * ( f_tmp - fi_tmp ) );

                    if constexpr ( detail::problem_operator<decltype( pb.implicit_part ), value_t> )
                    {
                        // Operator-based Shampine's trick.

                        shampine_trick_caller
                            .template operator()<l>( gamma * dt, pb.implicit_part.f_t( tn ), u_sm2pl, f_D_u, u_tmp, shampine_element );

                        // err_R = J_R^{-1} dt/6 (F_R(u^{s+2}) - F_R(u^{s+1}))

                        shampine_trick_caller.template operator()<1>( gamma * dt, pb.implicit_part.f_t( tn ), u_sm2pl, rhs_R, u_tmp, err_R );
                    }
                    else
                    {
                        // Reuse the factorization prepared by the implicit stages.
                        shampine_trick_caller.initialize();

                        shampine_trick_caller.template apply<l>( f_D_u, u_tmp, shampine_element );

                        // err_R = J_R^{-1} dt/6 (F_R(u^{s+2}) - F_R(u^{s+1}))
                        shampine_trick_caller.template apply<1>( rhs_R, u_tmp, err_R );

                        shampine_trick_caller.finalize();
                    }

                    state_algebra_t::add( compensation, 0.5 * dt * fi_tmp );
                    state_algebra_t::add( compensation, 0.5 * dt * f_tmp );
                    state_algebra_t::add( compensation, dt / ( 2. - 4. * gamma ) * shampine_element );
                    compensated_commit( us_s, u_np1 );

                    value_t err_R_scalar = detail::error_algebra<state_t>::estimate_squared( err_R,
                        un,
                        u_np1,
                        _info.absolute_tolerance,
                        _info.relative_tolerance );

                    _info.error   = std::max( err_D_scalar, err_R_scalar );
                    _info.success = _info.error < 1.0;

                    value_t new_dt;

                    value_t fac = raw_fac_from_current_error();

                    // Fortran facp memory correction: use the previous accepted
                    // error/step only when diffusion dominates reaction. Because Ponio
                    // stores RMS^2, the Fortran factor 10 becomes 100 here.
                    if ( errp != 0.0 && !previous_attempt_rejected && err_D_scalar >= 100.0 * err_R_scalar )
                    {
                        // Convert the Fortran RMS formula to Ponio's RMS^2 convention.
                        value_t const facp = std::pow( errp, 0.25 ) * fac * fac * ( dt / hp );
                        fac                = std::min( fac, facp );
                    }

                    fac = finalize_fac( fac );

                    new_dt = dt * fac;

                    // accepted step
                    if ( _info.success )
                    {
                        compensated_time_update( dt );

                        // Store the accepted error and step before updating dt.
                        errp = _info.error;
                        hp   = dt;

                        // After a successful retry, Fortran forbids immediate step
                        // growth and resets the consecutive-rejection counter.
                        if ( previous_attempt_rejected )
                        {
                            new_dt = std::min( new_dt, dt );
                            nrej   = 0;
                        }

                        dt = new_dt;

                        facmax                    = 2.0;
                        previous_attempt_rejected = false;
                    }
                    else
                    {
                        // rejected step: tn is unchanged
                        std::swap( un, u_np1 );
                        dt = finalize_rejection_dt( new_dt );
                    }
                }
                else
                {
                    // Reuse the reaction values produced by the implicit stages.

                    if constexpr ( detail::problem_operator<decltype( pb.implicit_part ), value_t> )
                    {
                        // Operator-based problems keep the historical Shampine interface.

                        shampine_trick_caller
                            .template operator()<l>( gamma * dt, pb.implicit_part.f_t( tn ), u_sm2pl, f_D_u, u_tmp, shampine_element );
                    }
                    else
                    {
                        // Reuse the factorization prepared by the implicit stages.
                        shampine_trick_caller.initialize();

                        shampine_trick_caller.template apply<l>( f_D_u, u_tmp, shampine_element );

                        shampine_trick_caller.finalize();
                    }

                    compensated_time_update( dt );

                    state_algebra_t::add( compensation, 0.5 * dt * fi_tmp );
                    state_algebra_t::add( compensation, 0.5 * dt * f_tmp );
                    state_algebra_t::add( compensation, dt / ( 2. - 4. * gamma ) * shampine_element );
                    compensated_commit( us_s, u_np1 );
                }
            }
            else
            {
                auto& f_D_u = U[0]; // U[0] is free after the ROCK recurrence

                // fe_tmp still stores F_D(u^{s+1}); only F_D(u^{s+3})
                // must be evaluated for the final diffusion correction.
                _info.number_of_eval[0] += 1;
                pb.explicit_part( t_diffusion_coupling, u_sp3, f_D_u );
                f_D_u = static_cast<state_t>( f_D_u - fe_tmp );

                // Reuse the reaction values left by the two implicit stages,
                // matching the Fortran `fnc` semantics.
                compensated_time_update( dt );

                state_algebra_t::add( compensation, 0.5 * dt * fi_tmp );
                state_algebra_t::add( compensation, 0.5 * dt * f_tmp );
                state_algebra_t::add( compensation, dt / ( 2. - 4. * gamma ) * f_D_u );
                compensated_commit( us_s, u_np1 );
            }
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

        /**
         * @brief set absolute tolerance in chained config
         *
         * @param tol_ tolerance
         * @return auto& returns this object
         */
        template <typename rock_t = rock_coeff>
            requires std::same_as<rock_t, rock_coeff> && is_embedded
        auto&
        abs_tol( value_t tol_ )
        {
            info().absolute_tolerance = tol_;
            return *this;
        }

        /**
         * @brief set relative tolerance in chained config
         *
         * @param tol_ tolerance
         * @return auto& returns this object
         */
        template <typename rock_t = rock_coeff>
            requires std::same_as<rock_t, rock_coeff> && is_embedded
        auto&
        rel_tol( value_t tol_ )
        {
            info().relative_tolerance = tol_;
            return *this;
        }
    };

    // cppcheck-suppress-begin unusedFunction

    /**
     * @brief helper function to build PIROCK algorithm
     *
     * @tparam l                       Number of augmented stages (l = 1, 2)
     * @tparam is_embedded             Set adaptive time step method (default: false)
     * @tparam value_t                 Type of coefficients
     * @tparam alpha_beta_computer_t   Choice of computing of \f$\alpha\f$ and \f$\beta\f$ parameters
     * @tparam eig_computer_t          Computing of eigenvalues of explicit part (diffusion)
     * @tparam shampine_trick_caller_t Computing of Shampine's trick
     * @param alpha_beta_computer      \f$\alpha\f$ and \f$\beta\f$ computer object
     * @param eig_computer             Eigenvalue computer of explicit part of the problem
     * @param shampine_trick_caller    Shampine's trick computer
     */
    template <std::size_t l = 1, bool is_embedded = false, typename value_t = double, typename alpha_beta_computer_t, typename eig_computer_t, typename shampine_trick_caller_t>
    auto
    pirock( alpha_beta_computer_t&& alpha_beta_computer, eig_computer_t&& eig_computer, shampine_trick_caller_t&& shampine_trick_caller )
    {
        return pirock_impl<l, alpha_beta_computer_t, eig_computer_t, shampine_trick_caller_t, is_embedded, value_t>(
            std::forward<alpha_beta_computer_t>( alpha_beta_computer ),
            std::forward<eig_computer_t>( eig_computer ),
            std::forward<shampine_trick_caller_t>( shampine_trick_caller ) );
    }

    /**
     * @brief helper function to build PIROCK algorithm
     *
     * @tparam l                     Number of augmented stages (l = 1, 2)
     * @tparam value_t               Type of coefficients
     * @tparam alpha_beta_computer_t Choice of computing of \f$\alpha\f$ and \f$\beta\f$ parameters
     * @tparam eig_computer_t        Computing of eigenvalues of explicit part (diffusion)
     * @param alpha_beta_computer    \f$\alpha\f$ and \f$\beta\f$ computer object
     * @param eig_computer           Eigenvalue computer of explicit part of the problem
     *
     * @note Without a Shampine's trick caller, this method is fixed time step.
     */
    template <std::size_t l = 1, typename value_t = double, typename alpha_beta_computer_t, typename eig_computer_t>
    auto
    pirock( alpha_beta_computer_t&& alpha_beta_computer, eig_computer_t&& eig_computer )
    {
        return pirock_impl<l, alpha_beta_computer_t, eig_computer_t, void, false, value_t>(
            std::forward<alpha_beta_computer_t>( alpha_beta_computer ),
            std::forward<eig_computer_t>( eig_computer ) );
    }

    /**
     * @brief helper function to build PIROCK algorithm
     *
     * @tparam l              Number of augmented stages (l = 1, 2)
     * @tparam value_t        Type of coefficients
     * @tparam eig_computer_t Computing of eigenvalues of explicit part (diffusion)
     * @param eig_computer    Eigenvalue computer of explicit part of the problem
     *
     * @note Without a Shampine's trick caller, this method is fixed time step, and without a \f$\alpha\f$ and \f$\beta\f$ computer,
     * parameters are fixed to \f$\beta = 0\f$.
     */
    template <std::size_t l = 1, typename value_t = double, typename eig_computer_t>
    auto
    pirock( eig_computer_t&& eig_computer )
    {
        return pirock<l, value_t>( beta_0<value_t>(), std::forward<eig_computer_t>( eig_computer ) );
    }

    /**
     * @brief helper function to build PIROCK algorithm
     *
     * @tparam l       Number of augmented stages (l = 1, 2)
     * @tparam value_t Type of coefficients
     *
     * @note Without a Shampine's trick caller, this method is fixed time step, without a \f$\alpha\f$ and \f$\beta\f$ computer, parameters
     * are fixed to \f$\beta = 0\f$ and without a eigenvalues computer this method use power method to estimate spectral radius of explicit
     * part (diffusion operator).
     */
    template <std::size_t l = 1, typename value_t = double>
    auto
    pirock()
    {
        return pirock<l, value_t>( rock::detail::power_method() );
    }

    /**
     * @brief helper function to build PIROCK algorithm
     *
     * @tparam value_t        Type of coefficients
     * @tparam eig_computer_t Computing of eigenvalues of explicit part (diffusion)
     * @param eig_computer    Eigenvalue computer of explicit part of the problem
     *
     * @note This function builds a PIROCK algorithm with parameters \f$\ell=2\f$ and \f$\alpha = 1\f$ and specific spectral radius
     * estimator of explicit part (diffusion operator).
     */
    template <typename value_t = double, typename eig_computer_t>
    auto
    pirock_a1( eig_computer_t&& eig_computer )
    {
        return pirock<2, value_t>( alpha_fixed<value_t>( 1.0 ), std::forward<eig_computer_t>( eig_computer ) );
    }

    /**
     * @brief helper function to build PIROCK algorithm
     *
     * @tparam value_t Type of coefficients
     *
     * @note This function builds a PIROCK algorithm with parameters \f$\ell=2\f$ and \f$\alpha = 1\f$ and power method to estimate spectral
     * radius.
     */
    template <typename value_t = double>
    auto
    pirock_a1()
    {
        return pirock_a1<value_t>( rock::detail::power_method() );
    }

    /**
     * @brief helper function to build PIROCK algorithm
     *
     * @tparam value_t        Type of coefficients
     * @tparam eig_computer_t Computing of eigenvalues of explicit part (diffusion)
     * @param eig_computer    Eigenvalue computer of explicit part of the problem
     *
     * @note This function builds a PIROCK algorithm with parameters \f$\ell=1\f$ and \f$\beta = 0\f$ and specific spectral radius estimator
     * of explicit part (diffusion operator).
     */
    template <typename value_t = double, typename eig_computer_t>
    auto
    pirock_b0( eig_computer_t&& eig_computer )
    {
        return pirock<1, value_t>( beta_0<value_t>(), std::forward<eig_computer_t>( eig_computer ) );
    }

    /**
     * @brief helper function to build PIROCK algorithm
     *
     * @tparam value_t Type of coefficients
     *
     * @note This function builds a PIROCK algorithm with parameters \f$\ell=1\f$ and \f$\beta = 0\f$ and power method to estimate spectral
     * radius.
     */
    template <typename value_t = double>
    auto
    pirock_b0()
    {
        return pirock_b0<value_t>( rock::detail::power_method() );
    }

    // cppcheck-suppress-end unusedFunction

    // --- PIROCK DIFFUSION-ADVECTION ------------------------------------------

    /**
     * @class pirock_DA_impl
     * @brief implementation of PIROCK method for diffusion-advection problems
     *
     * @tparam _l                    Number of augmented stages (l = 1, 2)
     * @tparam alpha_beta_computer_t Choice of computing of \f$\alpha\f$ and \f$\beta\f$ parameters
     * @tparam eig_computer_t        Computing of eigenvalues of diffusion part
     * @tparam _is_embedded          Set adaptive time step method (default: false)
     * @tparam value_t               Type of coefficients
     */
    template <std::size_t _l, typename alpha_beta_computer_t, typename eig_computer_t, bool _is_embedded = false, typename _value_t = double>
    struct pirock_DA_impl
    {
        static constexpr bool is_embedded = _is_embedded;

        static constexpr std::size_t l       = _l;
        static constexpr bool is_imex_method = true;
        // number_of_eval counts one evaluation per operator: 0 diffusion, 1 advection.
        static constexpr std::size_t N_operators                               = 2;
        static constexpr std::size_t N_stages                                  = stages::dynamic;
        static constexpr std::size_t N_storage                                 = is_embedded ? 18 : 17;
        static constexpr std::array<std::size_t, 1> persistent_storage_indices = { N_storage - 2 };
        static constexpr std::size_t order                                     = 2;
        static constexpr std::string_view id                                   = "PIROCK";

        using value_t         = _value_t;
        using rock_coeff      = rock::rock2_coeff<value_t>;
        using degree_computer = rock::detail::degree_computer<value_t, rock_coeff>;

        alpha_beta_computer_t alpha_beta_computer;
        eig_computer_t eig_computer;

        iteration_info<pirock_DA_impl> _info;

        bool compensated_summation_initialized = false;
        value_t time_compensation              = static_cast<value_t>( 0. );

        // Adaptive controller state, shared with the RD and RDA implementations.
        value_t facmax                 = static_cast<value_t>( 5. );
        value_t facd                   = static_cast<value_t>( 0.8 );
        std::size_t attempt_count      = 0;
        std::size_t nrejfac            = 0;
        bool previous_attempt_rejected = false;

        value_t errp = static_cast<value_t>( 0. );
        value_t hp   = static_cast<value_t>( 0. );

        value_t told     = static_cast<value_t>( 0. );
        std::size_t nrej = 0;

        /**
         * @brief Construct a new pirock_DA_impl object
         */
        pirock_DA_impl() = default;

        /**
         * @brief Construct a new pirock_DA_impl object
         *
         * @param _alpha_beta_computer alpha and beta computer object
         * @param _eig_computer        eigenvalue computer functor
         */
        pirock_DA_impl( alpha_beta_computer_t&& _alpha_beta_computer, eig_computer_t&& _eig_computer )
            : alpha_beta_computer( std::forward<alpha_beta_computer_t>( _alpha_beta_computer ) )
            , eig_computer( std::forward<eig_computer_t>( _eig_computer ) )
            , _info()
        {
        }

        /**
         * @brief iteration of PIROCK_DA method
         *
         * @tparam problem_t  type of \f$f\f$
         * @tparam state_t    type of current state
         * @tparam array_ki_t type of temporary stages
         * @param pb    problem with 2 operators: \f$F_D\f$ and \f$F_A\f$
         * @param tn    current time
         * @param un    current state
         * @param U     array of temporary stages
         * @param dt    current time step
         * @param u_np1 solution \f$u^{n+1}\f$ at time \f$t^{n+1}=t^n+\Delta t\f$
         */
        template <typename problem_t, typename state_t, typename array_ki_t>
        void
        operator()( problem_t& pb, value_t& tn, state_t& un, array_ki_t& U, value_t& dt, state_t& u_np1 )
        {
            // In problem_t pb.system (a tuple):
            // 0: diffusion operator
            // 1: advection operator
            using diffusion_op = std::integral_constant<std::size_t, 0>;
            using advection_op = std::integral_constant<std::size_t, 1>;

            // U worker references. Adaptive runs reserve U[6] for the diffusion
            // value at (tn, un); later work arrays are shifted by one position.
            constexpr std::size_t adaptive_shift = is_embedded ? 1 : 0;

            auto& u_j        = U[0];
            auto& u_jm1      = U[1];
            auto& u_jm2      = U[2];
            auto& u_sm2      = U[3];
            auto& fd_tmp     = U[4];
            auto& fd_tmp_bis = U[5];
            auto& fa_K       = U[6 + adaptive_shift];
            auto& fa_tmp     = U[7 + adaptive_shift];
            auto& fa_tmp_bis = U[8 + adaptive_shift];
            auto& us_sm1     = U[9 + adaptive_shift];
            auto& us_s       = U[10 + adaptive_shift];
            auto& u_sp3      = U[11 + adaptive_shift];
            auto& u_sp4      = U[12 + adaptive_shift];
            auto& u_sp5      = U[13 + adaptive_shift];

            auto& rock_increment   = U[N_storage - 3];
            auto& compensation     = U[N_storage - 2];
            auto& compensation_tmp = U[N_storage - 1];

            // The compensated sums below update a state in place. state_algebra
            // carries that operation for the state type at hand: compound
            // assignment where the type provides it, as Eigen vectors do, and
            // assignment from an expression where it does not, as for samurai
            // fields. The accumulation itself reads the same in both cases.
            using state_algebra_t = ::ponio::linear_algebra::state_algebra<state_t>;

            _info.reset_eval();

            auto& diffusion = std::get<diffusion_op::value>( pb.system );

            value_t const rock2_stability_constant = ( l == 1 ) ? static_cast<value_t>( 0.432 ) : static_cast<value_t>( 0.811 );

            auto [mdeg, deg_index, start_index, n_eval] = degree_computer::compute_n_stages_optimal_degree( rock::rock_order::rock_2(),
                eig_computer,
                diffusion,
                tn,
                un,
                dt,
                U,
                4,
                rock2_stability_constant );

            std::size_t const s = mdeg + 2;

            _info.number_of_stages  = s + l + 3;
            _info.number_of_eval[0] = n_eval;

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

            value_t const alpha = alpha_beta_computer.alpha( s, l );
            value_t const beta  = alpha_beta_computer.beta( s, l );
            value_t const gamma = 1. - 0.5 * std::numbers::sqrt2;

            value_t const t_diffusion_coupling = tn + 0.5 * ( 1. - beta ) * dt;
            value_t const t_advection_K        = tn;
            value_t const t_advection_sp4      = tn + dt / 3.;
            value_t const t_advection_sp5      = tn + 2. * dt / 3.;

            if ( !compensated_summation_initialized )
            {
                compensation = un;
                state_algebra_t::scale( compensation, value_t( 0 ) );
                compensated_summation_initialized = true;
            }

            auto compensated_update = [&]( state_t& state, auto const& increment )
            {
                state_algebra_t::add( compensation, increment );
                compensation_tmp = state;
                state            = state + compensation;
                state_algebra_t::add( compensation, compensation_tmp - state );
            };

            auto compensated_commit = [&]( state_t const& state, state_t& result )
            {
                compensation_tmp = state;
                result           = state + compensation;
                state_algebra_t::add( compensation, compensation_tmp - result );
            };

            auto compensated_time_update = [&]( value_t increment )
            {
                value_t const previous_time = tn;
                time_compensation += increment;
                tn = tn + time_compensation;
                time_compensation += previous_time - tn;
            };

            auto raw_fac_from_current_error = [&]() -> value_t
            {
                ++attempt_count;

                if ( attempt_count % 100 == 0 )
                {
                    if ( nrejfac >= 10 )
                    {
                        facd *= 0.8;
                    }
                    if ( nrejfac == 0 )
                    {
                        facd *= 1.02;
                    }
                    facd    = std::min( 0.8, std::max( 0.4, facd ) );
                    nrejfac = 0;
                }

                // _info.error stores the square of the effective RMS error.
                return std::pow( _info.error, -0.25 );
            };

            auto finalize_fac = [&]( value_t fac ) -> value_t
            {
                if ( previous_attempt_rejected )
                {
                    facmax = 1.0;
                }

                return std::min( facmax, std::max( 0.1, facd * fac ) );
            };

            auto finalize_rejection_dt = [&]( value_t new_dt ) -> value_t
            {
                value_t retry_dt = 0.8 * new_dt;

                if ( told == tn )
                {
                    ++nrej;
                    if ( nrej == 10 )
                    {
                        retry_dt = 1.0e-5;
                    }
                }
                told = tn;

                ++nrejfac;
                previous_attempt_rejected = true;

                return retry_dt;
            };

            // Diffusion stabilization procedure.
            u_j   = un;
            u_jm2 = un;

            value_t const mu_1 = rock_coeff::recf[start_index - 1];

            value_t t_jm1 = tn + dt * alpha * mu_1;
            value_t t_jm2 = tn + dt * alpha * mu_1;
            value_t t_jm3 = tn;

            // Reuse F_D(tn, un) after an adaptive rejection. Since (tn, un)
            // is unchanged, this avoids one diffusion evaluation on each retry.
            if constexpr ( is_embedded )
            {
                auto& fd_tn_cache = U[6];
                if ( !previous_attempt_rejected )
                {
                    eval_diffusion( tn, un, fd_tn_cache );
                }
                rock_increment = alpha * dt * mu_1 * fd_tn_cache;
            }
            else
            {
                eval_diffusion( tn, un, fd_tmp );
                rock_increment = alpha * dt * mu_1 * fd_tmp;
            }

            u_jm1 = un;
            compensated_update( u_jm1, rock_increment );

            if ( mdeg < 2 )
            {
                u_j = u_jm1;
            }

            value_t t_sm2 = tn;

            for ( std::size_t j = 2; j < s - 2 + l + 1; ++j )
            {
                value_t mu_j;
                value_t kappa_j;

                if ( j <= mdeg )
                {
                    mu_j    = rock_coeff::recf[start_index + 2 * ( j - 2 ) + 1 - 1];
                    kappa_j = rock_coeff::recf[start_index + 2 * ( j - 2 ) + 2 - 1];
                }
                else
                {
                    std::size_t const additional_stage = j - mdeg;
                    std::size_t const recf2_index      = 4 * ( deg_index - 1 ) + 2 * ( additional_stage - 1 );
                    mu_j                               = rock_coeff::recf2[recf2_index];
                    kappa_j                            = rock_coeff::recf2[recf2_index + 1];
                }

                value_t const nu_j = -1.0 - kappa_j;

                eval_diffusion( t_jm1, u_jm1, fd_tmp );

                if ( j <= mdeg )
                {
                    state_algebra_t::scale( rock_increment, kappa_j );
                    state_algebra_t::add( rock_increment, alpha * mu_j * dt * fd_tmp );
                    u_j = u_jm1;
                    compensated_update( u_j, rock_increment );
                }
                else
                {
                    u_j = alpha * mu_j * dt * fd_tmp - nu_j * u_jm1 - kappa_j * u_jm2;
                }

                t_jm1 = alpha * dt * mu_j - nu_j * t_jm2 - kappa_j * t_jm3;

                if ( j == s - 2 )
                {
                    u_sm2 = u_j;
                    t_sm2 = t_jm1;
                }
                if ( j < s - 2 + l )
                {
                    std::swap( u_jm2, u_jm1 );
                    std::swap( u_jm1, u_j );
                }

                t_jm3 = t_jm2;
                t_jm2 = t_jm1;
            }

            value_t const sigma   = rock_coeff::fp1[deg_index - 1];
            value_t const sigma_a = 0.5 * ( 1.0 - alpha ) + alpha * sigma;
            value_t const tau     = sigma * rock_coeff::fp2[deg_index - 1] + sigma * sigma;
            value_t const tau_a   = 0.5 * detail::power<2>( alpha - 1. ) + 2. * alpha * ( 1. - alpha ) * sigma + alpha * alpha * tau;

            value_t err_D_scalar = static_cast<value_t>( 0. );

            // ROCK2 finishing procedure and diffusion estimator.
            eval_diffusion( t_sm2, u_sm2, fd_tmp );
            rock_increment = fd_tmp;
            us_sm1         = u_sm2;
            compensated_update( us_sm1, sigma_a * dt * fd_tmp );

            auto& err_D = rock_increment;
            eval_diffusion( t_sm2 + sigma_a * dt, us_sm1, fd_tmp );

            err_D = sigma_a * ( 1. - tau_a / ( sigma_a * sigma_a ) ) * dt * ( fd_tmp - err_D );

            us_s = us_sm1;
            compensated_update( us_s, sigma_a * dt * fd_tmp - err_D );

            if constexpr ( is_embedded )
            {
                err_D_scalar = detail::error_algebra<state_t>::estimate_squared( err_D,
                    un,
                    us_s,
                    _info.absolute_tolerance,
                    _info.relative_tolerance );

                // Reject before the advection finishing procedure when diffusion
                // alone already violates the requested tolerance.
                if ( err_D_scalar >= 1.0 )
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

            // Starting value for the advection finishing procedure: K = K_{s-2+ell}.
            auto& K = u_j;

            eval_diffusion( t_diffusion_coupling, K, fd_tmp );
            eval_advection( t_advection_K, K, fa_K );

            // Minimal diffusion-advection finishing stages.
            u_sp3 = K + ( 1. - 2. * gamma ) * dt * fa_K;
            u_sp4 = K + dt / 3. * fa_K;

            eval_advection( t_advection_sp4, u_sp4, fa_tmp );

            u_sp5 = K + 2. / 3. * beta * dt * fd_tmp + 2. / 3. * dt * fa_tmp;

            eval_advection( t_advection_sp5, u_sp5, fa_tmp_bis );

            // Final diffusion-advection coupling.
            eval_diffusion( t_diffusion_coupling, u_sp3, fd_tmp_bis );
            fd_tmp_bis = static_cast<state_t>( fd_tmp_bis - fd_tmp );

            state_algebra_t::add( compensation, 0.25 * dt * fa_K );
            state_algebra_t::add( compensation, dt / ( 2. - 4. * gamma ) * fd_tmp_bis );
            state_algebra_t::add( compensation, 0.75 * dt * fa_tmp_bis );
            compensated_commit( us_s, u_np1 );

            if constexpr ( is_embedded )
            {
                // u_sp4 is no longer needed and can store the advection estimator.
                auto& err_A = u_sp4;
                err_A       = -0.15 * dt * fa_K + 0.3 * dt * fa_tmp - 0.15 * dt * fa_tmp_bis;

                value_t const err_A_scalar = detail::error_algebra<state_t>::estimate_squared( err_A,
                    un,
                    u_np1,
                    _info.absolute_tolerance,
                    _info.relative_tolerance );

                // The advection estimator is third order. In RMS^2 form,
                // ||err_A||^(2/3) becomes err_A_scalar^(2/3).
                value_t const err_A_effective_scalar = std::pow( err_A_scalar, 2. / 3. );

                _info.error   = std::max( err_D_scalar, err_A_effective_scalar );
                _info.success = _info.error < 1.0;

                value_t fac = raw_fac_from_current_error();

                // Apply the controller memory correction only when diffusion
                // clearly dominates the advection error.
                if ( errp != 0.0 && !previous_attempt_rejected && err_D_scalar >= 100.0 * err_A_effective_scalar )
                {
                    value_t const facp = std::pow( errp, 0.25 ) * fac * fac * ( dt / hp );
                    fac                = std::min( fac, facp );
                }

                fac            = finalize_fac( fac );
                value_t new_dt = dt * fac;

                if ( _info.success )
                {
                    compensated_time_update( dt );

                    errp = _info.error;
                    hp   = dt;

                    if ( previous_attempt_rejected )
                    {
                        new_dt = std::min( new_dt, dt );
                        nrej   = 0;
                    }

                    dt = new_dt;

                    facmax                    = 2.0;
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
                compensated_time_update( dt );
            }
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

        /**
         * @brief set absolute tolerance in chained config
         *
         * @param tol_ tolerance
         * @return auto& returns this object
         */
        template <typename rock_t = rock_coeff>
            requires std::same_as<rock_t, rock_coeff> && is_embedded
        auto&
        abs_tol( value_t tol_ )
        {
            info().absolute_tolerance = tol_;
            return *this;
        }

        /**
         * @brief set relative tolerance in chained config
         *
         * @param tol_ tolerance
         * @return auto& returns this object
         */
        template <typename rock_t = rock_coeff>
            requires std::same_as<rock_t, rock_coeff> && is_embedded
        auto&
        rel_tol( value_t tol_ )
        {
            info().relative_tolerance = tol_;
            return *this;
        }
    };

    // cppcheck-suppress-begin unusedFunction

    /**
     * @brief helper function to build PIROCK diffusion-advection algorithm
     *
     * @tparam l                     Number of augmented stages (l = 1, 2)
     * @tparam is_embedded           Set adaptive time step method (default: false)
     * @tparam value_t               Type of coefficients
     * @tparam alpha_beta_computer_t Choice of computing of \f$\alpha\f$ and \f$\beta\f$ parameters
     * @tparam eig_computer_t        Computing of eigenvalues of diffusion part
     * @param alpha_beta_computer    \f$\alpha\f$ and \f$\beta\f$ computer object
     * @param eig_computer           Eigenvalue computer of diffusion part
     */
    template <std::size_t l = 1, bool is_embedded = false, typename value_t = double, typename alpha_beta_computer_t, typename eig_computer_t>
    auto
    pirock_DA( alpha_beta_computer_t&& alpha_beta_computer, eig_computer_t&& eig_computer )
    {
        return pirock_DA_impl<l, alpha_beta_computer_t, eig_computer_t, is_embedded, value_t>(
            std::forward<alpha_beta_computer_t>( alpha_beta_computer ),
            std::forward<eig_computer_t>( eig_computer ) );
    }

    /**
     * @brief helper function to build PIROCK diffusion-advection algorithm
     *
     * @tparam l              Number of augmented stages (l = 1, 2)
     * @tparam value_t        Type of coefficients
     * @tparam eig_computer_t Computing of eigenvalues of diffusion part
     * @param eig_computer    Eigenvalue computer of diffusion part
     *
     * @note Without an alpha/beta computer, parameters are fixed to \f$\beta=0\f$.
     */
    template <std::size_t l = 1, typename value_t = double, typename eig_computer_t>
    auto
    pirock_DA( eig_computer_t&& eig_computer )
    {
        return pirock_DA<l, false, value_t>( beta_0<value_t>(), std::forward<eig_computer_t>( eig_computer ) );
    }

    /**
     * @brief helper function to build PIROCK diffusion-advection algorithm
     *
     * @tparam l       Number of augmented stages (l = 1, 2)
     * @tparam value_t Type of coefficients
     *
     * @note Without an eigenvalue computer, the power method is used for the diffusion part.
     */
    template <std::size_t l = 1, typename value_t = double>
    auto
    pirock_DA()
    {
        return pirock_DA<l, value_t>( rock::detail::power_method() );
    }

    /**
     * @brief helper function to build PIROCK diffusion-advection algorithm
     *
     * @tparam value_t        Type of coefficients
     * @tparam eig_computer_t Computing of eigenvalues of diffusion part
     * @param eig_computer    Eigenvalue computer of diffusion part
     *
     * @note Builds the \f$\ell=2\f$, \f$\alpha=1\f$ fixed-step variant.
     */
    template <typename value_t = double, typename eig_computer_t>
    auto
    pirock_DA_a1( eig_computer_t&& eig_computer )
    {
        return pirock_DA<2, false, value_t>( alpha_fixed<value_t>( 1.0 ), std::forward<eig_computer_t>( eig_computer ) );
    }

    /**
     * @brief helper function to build PIROCK diffusion-advection algorithm
     *
     * @tparam value_t Type of coefficients
     *
     * @note Builds the \f$\ell=2\f$, \f$\alpha=1\f$ fixed-step variant with the power method.
     */
    template <typename value_t = double>
    auto
    pirock_DA_a1()
    {
        return pirock_DA_a1<value_t>( rock::detail::power_method() );
    }

    /**
     * @brief helper function to build PIROCK diffusion-advection algorithm
     *
     * @tparam value_t        Type of coefficients
     * @tparam eig_computer_t Computing of eigenvalues of diffusion part
     * @param eig_computer    Eigenvalue computer of diffusion part
     *
     * @note Builds the \f$\ell=1\f$, \f$\beta=0\f$ fixed-step variant.
     */
    template <typename value_t = double, typename eig_computer_t>
    auto
    pirock_DA_b0( eig_computer_t&& eig_computer )
    {
        return pirock_DA<1, false, value_t>( beta_0<value_t>(), std::forward<eig_computer_t>( eig_computer ) );
    }

    /**
     * @brief helper function to build PIROCK diffusion-advection algorithm
     *
     * @tparam value_t Type of coefficients
     *
     * @note Builds the \f$\ell=1\f$, \f$\beta=0\f$ fixed-step variant with the power method.
     */
    template <typename value_t = double>
    auto
    pirock_DA_b0()
    {
        return pirock_DA_b0<value_t>( rock::detail::power_method() );
    }

    // cppcheck-suppress-end unusedFunction

    // --- PIROCK REACTION-DIFFUSION-ADVECTION --------------------------------

    /**
     * @class pirock_RDA_impl
     * @brief implementation of PIROCK method for reaction-diffusion-advection problem
     *
     * @tparam _l                       Number of augmented stages (l = 1, 2)
     * @tparam alpha_beta_computer_t    Choice of computing of \f$\alpha\f$ and \f$\beta\f$ parameters
     * @tparam eig_computer_t           Computing of eigenvalues of explicit part (diffusion)
     * @tparam _shampine_trick_caller_t Computing of Shampine's trick
     * @tparam _is_embedded             Set adaptive time step method (default: false)
     * @tparam value_t                  Type of coefficients
     */
    template <std::size_t _l,
        typename alpha_beta_computer_t,
        typename eig_computer_t,
        typename _shampine_trick_caller_t = void,
        bool _is_embedded                 = false,
        typename _value_t                 = double>
    struct pirock_RDA_impl
    {
        static constexpr bool is_embedded           = _is_embedded;
        static constexpr bool shampine_trick_enable = !std::is_void_v<_shampine_trick_caller_t>;

        static constexpr std::size_t l       = _l;
        static constexpr bool is_imex_method = true;
        // number_of_eval counts one evaluation per operator: 0 diffusion,
        // 1 reaction, 2 advection. The problem itself lists the reaction first.
        static constexpr std::size_t N_operators = 3;
        static constexpr std::size_t N_stages    = stages::dynamic;
        // clang-format off
        static constexpr std::size_t N_storage   = std::conditional_t<shampine_trick_enable,
                                                    std::integral_constant<std::size_t, 29>,
                                                    std::integral_constant<std::size_t, 23>
                                                >::value;
        // clang-format on
        static constexpr std::array<std::size_t, 1> persistent_storage_indices = { N_storage - 2 };
        static constexpr std::size_t order                                     = 2;
        static constexpr std::string_view id                                   = "PIROCK";

        using value_t                 = _value_t;
        using rock_coeff              = rock::rock2_coeff<value_t>;
        using degree_computer         = rock::detail::degree_computer<value_t, rock_coeff>;
        using shampine_trick_caller_t = typename std::conditional_t<shampine_trick_enable, _shampine_trick_caller_t, bool>;

        alpha_beta_computer_t alpha_beta_computer;
        eig_computer_t eig_computer;
        shampine_trick_caller_t shampine_trick_caller;

        iteration_info<pirock_RDA_impl> _info;

        bool compensated_summation_initialized = false;
        value_t time_compensation              = static_cast<value_t>( 0. );

        // Adaptive controller state. These values persist across step attempts,
        // as in the reference Fortran `rockcore` controller.
        value_t facmax                 = static_cast<value_t>( 5. );
        value_t facd                   = static_cast<value_t>( 0.8 );
        std::size_t attempt_count      = 0;
        std::size_t nrejfac            = 0;
        bool previous_attempt_rejected = false;

        // Error and step size of the last accepted step for the controller
        // with memory.
        value_t errp = static_cast<value_t>( 0. );
        value_t hp   = static_cast<value_t>( 0. );

        // Consecutive-rejection safeguard from the Fortran controller.
        value_t told     = static_cast<value_t>( 0. );
        std::size_t nrej = 0;

        /**
         * @brief Construct a new pirock_RDA_impl object
         *
         */
        pirock_RDA_impl() = default;

        /**
         * @brief Construct a new pirock_RDA_impl object
         *
         * @param _alpha_beta_computer computer of parameters alpha and beta object that has two member functions which take number of
         * stages (s) l parameter
         * @param _eig_computer        eigenvalue computer functor (that take as argument the function of explicit part, the current time,
         * the current state and the current time step)
         */
        pirock_RDA_impl( alpha_beta_computer_t&& _alpha_beta_computer, eig_computer_t&& _eig_computer )
            : alpha_beta_computer( std::forward<alpha_beta_computer_t>( _alpha_beta_computer ) )
            , eig_computer( std::forward<eig_computer_t>( _eig_computer ) )
            , shampine_trick_caller( false )
            , _info()
        {
        }

        /**
         * @brief Construct a new pirock_RDA_impl object with Shampine's trick
         *
         * @param _alpha_beta_computer   alpha and beta computer object
         * @param _eig_computer          eigenvalue computer functor
         * @param _shampine_trick_caller Shampine's trick functor
         */
        template <typename _shampine_trick_caller_t_>
            requires std::same_as<_shampine_trick_caller_t_, shampine_trick_caller_t>
                      && std::same_as<std::bool_constant<shampine_trick_enable>, std::true_type>
        pirock_RDA_impl( alpha_beta_computer_t&& _alpha_beta_computer,
            eig_computer_t&& _eig_computer,
            _shampine_trick_caller_t_&& _shampine_trick_caller )
            : alpha_beta_computer( std::forward<alpha_beta_computer_t>( _alpha_beta_computer ) )
            , eig_computer( std::forward<eig_computer_t>( _eig_computer ) )
            , shampine_trick_caller( std::forward<_shampine_trick_caller_t_>( _shampine_trick_caller ) )
            , _info()
        {
        }

        /**
         * @brief iteration of PIROCK_RDA method
         *
         * @tparam problem_t  type of \f$f\f$
         * @tparam state_t    type of current state
         * @tparam array_ki_t type of temporary stages
         * @param pb    problem with 3 operators: \f$F_R\f$, \f$F_D\f$ and \f$F_A\f$
         * @param tn    current time
         * @param un    current state
         * @param U     array of temporary stages
         * @param dt    current time step
         * @param u_np1 solution \f$u^{n+1}\f$ at time \f$t^{n+1} = t^n + \Delta t\f$
         */
        template <typename problem_t, typename state_t, typename array_ki_t>
        void
        operator()( problem_t& pb, value_t& tn, state_t& un, array_ki_t& U, value_t& dt, state_t& u_np1 )
        {
            // In problem_t pb.system (a tuple):
            // 0: reaction operator
            // 1: diffusion operator
            // 2: advection operator
            using reaction_op  = std::integral_constant<std::size_t, 0>;
            using diffusion_op = std::integral_constant<std::size_t, 1>;
            using advection_op = std::integral_constant<std::size_t, 2>;

            using reaction_problem_t = std::tuple_element_t<reaction_op::value, decltype( pb.system )>;

            static_assert(
                detail::problem_operator<reaction_problem_t, value_t> || detail::problem_jacobian<reaction_problem_t, value_t, state_t>,
                "This kind of problem is not inversible in ponio" );

            // U worker references:
            // | index | variables          | mathematic representation                           |
            // |-------|--------------------|------------------------------------------------------|
            // | 0     | u_j                | current stage in pseudo-ROCK2                         |
            // | 1     | u_jm1              | previous stage in pseudo-ROCK2                        |
            // | 2     | u_jm2              | second previous stage in pseudo-ROCK2                 |
            // | 3     | u_sm2              | \f$u^{(s-2)}\f$                                      |
            // | 4     | fr_tmp             | \f$F_R(u^{(s+1)})\f$                                 |
            // | 5     | fr_tmp_bis         | \f$F_R(u^{(s+2)})\f$                                 |
            // | 6     | fd_tmp             | diffusion temporary / \f$F_D(u^{(s+1)})\f$           |
            // | 7     | fd_tmp_bis         | diffusion temporary                                  |
            // | 8     | fd_tn_cache        | cached \f$F_D(t_n,u_n)\f$ across rejected retries      |
            // | 10    | fa_tmp             | advection temporary                                  |
            // | 11    | fa_tmp_bis         | advection temporary                                  |
            // | 12    | us_sm1             | \f$u^{\star(s-1)}\f$                                |
            // | 13    | us_s               | ROCK2 diffusion solution including \f$-err_D\f$       |
            // | 14    | u_sp1              | \f$u^{(s+1)}\f$                                      |
            // | 15    | u_sp2              | \f$u^{(s+2)}\f$                                      |
            // | 16    | Fa_u_sp1           | \f$F_A(u^{(s+1)})\f$                                 |
            // | 17    | u_sp3              | \f$u^{(s+3)}\f$                                      |
            // | 18    | u_sp4              | \f$u^{(s+4)}\f$                                      |
            // | 19    | u_sp5              | \f$u^{(s+5)}\f$                                      |
            // | 20-25 | Shampine/error work arrays when enabled                          |
            // | last-2| rock_increment     | transported ROCK increment; then \f$err_D\f$         |
            // | last-1| compensation       | compensated summation remainder                       |
            // | last  | compensation_tmp   | temporary state for compensated summation             |

            _info.reset_eval();

            auto& reaction  = std::get<reaction_op::value>( pb.system );
            auto& diffusion = std::get<diffusion_op::value>( pb.system );

            // PIROCK uses the ROCK2 stability constants associated with the two
            // fixed damping choices: 0.432 for ell=1 and 0.811 for ell=2.
            value_t const rock2_stability_constant = ( l == 1 ) ? static_cast<value_t>( 0.432 ) : static_cast<value_t>( 0.811 );

            auto [mdeg, deg_index, start_index, n_eval] = degree_computer::compute_n_stages_optimal_degree( rock::rock_order::rock_2(),
                eig_computer,
                diffusion,
                tn,
                un,
                dt,
                U,
                4,
                rock2_stability_constant );

            std::size_t const s = mdeg + 2;

            _info.number_of_stages  = s + l + 5;
            _info.number_of_eval[0] = n_eval;

            auto eval_diffusion = [&]( value_t time, state_t const& state, state_t& output )
            {
                pb( diffusion_op(), time, state, output );
                ++_info.number_of_eval[0];
            };

            auto eval_reaction = [&]( value_t time, state_t const& state, state_t& output )
            {
                pb( reaction_op(), time, state, output );
                ++_info.number_of_eval[1];
            };

            auto eval_advection = [&]( value_t time, state_t const& state, state_t& output )
            {
                pb( advection_op(), time, state, output );
                ++_info.number_of_eval[2];
            };

            value_t const alpha = alpha_beta_computer.alpha( s, l );
            value_t const beta  = alpha_beta_computer.beta( s, l );
            value_t const gamma = 1. - 0.5 * std::numbers::sqrt2;

            value_t const t_sp2                = tn + ( 1. - gamma ) * dt;
            value_t const t_diffusion_coupling = tn + 0.5 * ( 1. - beta ) * dt;

            // The advection finishing method reduces to the classical three-stage
            // third-order explicit RK method on F_A alone.
            value_t const t_advection_sp1 = tn;
            value_t const t_advection_sp4 = tn + dt / 3.;
            value_t const t_advection_sp5 = tn + 2. * dt / 3.;

            auto& u_j         = U[0];
            auto& u_jm1       = U[1];
            auto& u_jm2       = U[2];
            auto& u_sm2       = U[3];
            auto& fr_tmp      = U[4];
            auto& fr_tmp_bis  = U[5];
            auto& fd_tmp      = U[6];
            auto& fd_tmp_bis  = U[7];
            auto& fd_tn_cache = U[8];
            auto& fa_tmp      = U[10];
            auto& fa_tmp_bis  = U[11];

            auto& rock_increment   = U[N_storage - 3];
            auto& compensation     = U[N_storage - 2];
            auto& compensation_tmp = U[N_storage - 1];

            // The compensated sums below update a state in place. state_algebra
            // carries that operation for the state type at hand: compound
            // assignment where the type provides it, as Eigen vectors do, and
            // assignment from an expression where it does not, as for samurai
            // fields. The accumulation itself reads the same in both cases.
            using state_algebra_t = ::ponio::linear_algebra::state_algebra<state_t>;

            if ( !compensated_summation_initialized )
            {
                compensation = un;
                state_algebra_t::scale( compensation, value_t( 0 ) );
                compensated_summation_initialized = true;
            }

            auto compensated_update = [&]( state_t& state, auto const& increment )
            {
                state_algebra_t::add( compensation, increment );
                compensation_tmp = state;
                state            = state + compensation;
                state_algebra_t::add( compensation, compensation_tmp - state );
            };

            auto compensated_commit = [&]( state_t const& state, state_t& result )
            {
                compensation_tmp = state;
                result           = state + compensation;
                state_algebra_t::add( compensation, compensation_tmp - result );
            };

            auto compensated_time_update = [&]( value_t increment )
            {
                value_t const previous_time = tn;
                time_compensation += increment;
                tn = tn + time_compensation;
                time_compensation += previous_time - tn;
            };

            auto raw_fac_from_current_error = [&]() -> value_t
            {
                ++attempt_count;

                if ( attempt_count % 100 == 0 )
                {
                    if ( nrejfac >= 10 )
                    {
                        facd *= 0.8;
                    }
                    if ( nrejfac == 0 )
                    {
                        facd *= 1.02;
                    }
                    facd    = std::min( 0.8, std::max( 0.4, facd ) );
                    nrejfac = 0;
                }

                // _info.error stores the square of the effective RMS error.
                return std::pow( _info.error, -0.25 );
            };

            auto finalize_fac = [&]( value_t fac ) -> value_t
            {
                if ( previous_attempt_rejected )
                {
                    facmax = 1.0;
                }

                return std::min( facmax, std::max( 0.1, facd * fac ) );
            };

            auto finalize_rejection_dt = [&]( value_t new_dt ) -> value_t
            {
                value_t retry_dt = 0.8 * new_dt;

                if ( told == tn )
                {
                    ++nrej;
                    if ( nrej == 10 )
                    {
                        retry_dt = 1.0e-5;
                    }
                }
                told = tn;

                ++nrejfac;
                previous_attempt_rejected = true;

                return retry_dt;
            };

            // Diffusion stabilization procedure.
            u_j   = un;
            u_jm2 = un;

            value_t const mu_1 = rock_coeff::recf[start_index - 1];

            value_t t_jm1 = tn + dt * alpha * mu_1;
            value_t t_jm2 = tn + dt * alpha * mu_1;
            value_t t_jm3 = tn;

            // A rejected adaptive attempt leaves (tn, un) unchanged. Cache
            // F_D(tn, un) and reuse it directly on retries, as Fortran keeps `fn`
            // across a rejected attempt. U[8] is otherwise unused here, so this
            // requires no extra allocation or vector copy.
            if constexpr ( is_embedded )
            {
                if ( !previous_attempt_rejected )
                {
                    eval_diffusion( tn, un, fd_tn_cache );
                }
                rock_increment = alpha * dt * mu_1 * fd_tn_cache;
            }
            else
            {
                eval_diffusion( tn, un, fd_tmp );
                rock_increment = alpha * dt * mu_1 * fd_tmp;
            }
            u_jm1 = un;
            compensated_update( u_jm1, rock_increment );

            if ( mdeg < 2 )
            {
                u_j = u_jm1;
            }

            value_t t_sm2 = tn;

            for ( std::size_t j = 2; j < s - 2 + l + 1; ++j )
            {
                value_t mu_j;
                value_t kappa_j;

                if ( j <= mdeg )
                {
                    mu_j    = rock_coeff::recf[start_index + 2 * ( j - 2 ) + 1 - 1];
                    kappa_j = rock_coeff::recf[start_index + 2 * ( j - 2 ) + 2 - 1];
                }
                else
                {
                    std::size_t const additional_stage = j - mdeg;
                    std::size_t const recf2_index      = 4 * ( deg_index - 1 ) + 2 * ( additional_stage - 1 );
                    mu_j                               = rock_coeff::recf2[recf2_index];
                    kappa_j                            = rock_coeff::recf2[recf2_index + 1];
                }

                value_t const nu_j = -1.0 - kappa_j;

                eval_diffusion( t_jm1, u_jm1, fd_tmp );

                if ( j <= mdeg )
                {
                    state_algebra_t::scale( rock_increment, kappa_j );
                    state_algebra_t::add( rock_increment, alpha * mu_j * dt * fd_tmp );
                    u_j = u_jm1;
                    compensated_update( u_j, rock_increment );
                }
                else
                {
                    // The supplementary recf2 stages are applied directly to the
                    // stage values, as in the ROCK2 reference implementation.
                    u_j = alpha * mu_j * dt * fd_tmp - nu_j * u_jm1 - kappa_j * u_jm2;
                }

                t_jm1 = alpha * dt * mu_j - nu_j * t_jm2 - kappa_j * t_jm3;

                if ( j == s - 2 )
                {
                    u_sm2 = u_j;
                    t_sm2 = t_jm1;
                }
                if ( j < s - 2 + l )
                {
                    std::swap( u_jm2, u_jm1 );
                    std::swap( u_jm1, u_j );
                }

                t_jm3 = t_jm2;
                t_jm2 = t_jm1;
            }

            value_t const sigma   = rock_coeff::fp1[deg_index - 1];
            value_t const sigma_a = 0.5 * ( 1.0 - alpha ) + alpha * sigma;
            value_t const tau     = sigma * rock_coeff::fp2[deg_index - 1] + sigma * sigma;
            value_t const tau_a   = 0.5 * detail::power<2>( alpha - 1. ) + 2. * alpha * ( 1. - alpha ) * sigma + alpha * alpha * tau;

            value_t err_D_scalar = static_cast<value_t>( 0. );

            // ROCK2 finishing procedure and diffusion estimator.
            auto& us_sm1 = U[12];
            eval_diffusion( t_sm2, u_sm2, fd_tmp );
            rock_increment = fd_tmp;
            us_sm1         = u_sm2;
            compensated_update( us_sm1, sigma_a * dt * fd_tmp );

            auto& us_s  = U[13];
            auto& err_D = rock_increment;
            eval_diffusion( t_sm2 + sigma_a * dt, us_sm1, fd_tmp );

            err_D = sigma_a * ( 1. - tau_a / ( sigma_a * sigma_a ) ) * dt * ( fd_tmp - err_D );

            us_s = us_sm1;
            compensated_update( us_s, sigma_a * dt * fd_tmp - err_D );

            if constexpr ( is_embedded )
            {
                err_D_scalar = detail::error_algebra<state_t>::estimate_squared( err_D,
                    un,
                    us_s,
                    _info.absolute_tolerance,
                    _info.relative_tolerance );

                if ( err_D_scalar >= 1.0 )
                {
                    _info.error   = err_D_scalar;
                    _info.success = false;

                    value_t fac = raw_fac_from_current_error();
                    fac         = finalize_fac( fac );

                    value_t const new_dt = dt * fac;

                    u_np1 = un;
                    dt    = finalize_rejection_dt( new_dt );

                    return;
                }
            }

            auto& u_sm2pl  = u_j;
            auto& u_sp1    = U[14];
            auto& u_sp2    = U[15];
            auto& Fa_u_sp1 = U[16];

            u_sp1 = un;
            u_sp2 = un;

            // Two diagonally implicit reaction stages. The sparse matrix path
            // mirrors the validated RD implementation and reuses the latest
            // factorization between stages.
            if constexpr ( detail::problem_operator<reaction_problem_t, value_t> )
            {
                std::size_t n_eval_sp1 = 0;

                auto op_sp1 = ::ponio::linear_algebra::operator_algebra<state_t>::identity( un ) - gamma * dt * reaction.f_t( tn );

                auto rhs_sp1 = u_sm2pl;

                // The first implicit stage starts from u_{s-2+l}, the initial
                // guess the nonlinear iteration is given on the matrix path.
                u_sp1 = u_sm2pl;

                ::ponio::linear_algebra::operator_algebra<state_t>::solve( op_sp1, u_sp1, rhs_sp1, n_eval_sp1 );
                _info.number_of_eval[1] += n_eval_sp1;

                eval_reaction( tn, u_sp1, fr_tmp );
                eval_diffusion( t_diffusion_coupling, u_sp1, fd_tmp );
                eval_advection( t_advection_sp1, u_sp1, Fa_u_sp1 );

                auto& rhs_sp2 = U[17];
                rhs_sp2       = u_sm2pl + beta * dt * fd_tmp + dt * Fa_u_sp1 + ( 1. - 2. * gamma ) * dt * fr_tmp;

                std::size_t n_eval_sp2 = 0;

                auto op_sp2 = ::ponio::linear_algebra::operator_algebra<state_t>::identity( un ) - gamma * dt * reaction.f_t( t_sp2 );

                // The second implicit stage starts from its own right-hand side
                // c2, as it does on the matrix path.
                u_sp2 = rhs_sp2;

                ::ponio::linear_algebra::operator_algebra<state_t>::solve( op_sp2, u_sp2, rhs_sp2, n_eval_sp2 );
                _info.number_of_eval[1] += n_eval_sp2;

                eval_reaction( t_sp2, u_sp2, fr_tmp_bis );
            }
            else
            {
                using matrix_t                = std::decay_t<decltype( reaction.df( tn, un ) )>;
                using matrix_linear_algebra_t = ::ponio::linear_algebra::linear_algebra<matrix_t>;

                constexpr bool can_reuse_factorization = requires( matrix_t const& matrix, state_t const& rhs ) {
                                                             matrix_linear_algebra_t::identity_minus( matrix, value_t{} );
                                                             matrix_linear_algebra_t::factorize( matrix );
                                                             matrix_linear_algebra_t::solve_factorized( rhs );
                                                         };

                auto g_sp1 = [&]( state_t& u ) -> state_t
                {
                    eval_reaction( tn, u, fr_tmp );
                    return u - gamma * dt * fr_tmp - u_sm2pl;
                };

                if constexpr ( can_reuse_factorization )
                {
                    decltype( auto ) reaction_jacobian = reaction.df( tn, u_sm2pl );

                    matrix_t frozen_stage_matrix;
                    frozen_stage_matrix = matrix_linear_algebra_t::identity_minus( reaction_jacobian, gamma * dt );

                    matrix_linear_algebra_t::factorize( frozen_stage_matrix );

                    auto refresh_stage1_factorization = [&]( state_t const& u )
                    {
                        decltype( auto ) refreshed_reaction_jacobian = reaction.df( tn, u );

                        matrix_t refreshed_stage_matrix;
                        refreshed_stage_matrix = matrix_linear_algebra_t::identity_minus( refreshed_reaction_jacobian, gamma * dt );

                        matrix_linear_algebra_t::factorize( refreshed_stage_matrix );
                    };

                    u_sp1 = diagonal_implicit_runge_kutta::pirock_newton<value_t, state_t, matrix_t>( g_sp1,
                        u_sm2pl,
                        refresh_stage1_factorization );
                }
                else
                {
                    auto identity = matrix_linear_algebra_t::identity( un );

                    auto dg_sp1 = [&]( state_t& u ) -> matrix_t
                    {
                        return identity - gamma * dt * reaction.df( tn, u );
                    };

                    u_sp1 = diagonal_implicit_runge_kutta::newton<value_t>( g_sp1,
                        dg_sp1,
                        u_sm2pl,
                        matrix_linear_algebra_t::solver,
                        ponio::default_config::newton_tolerance,
                        ponio::default_config::newton_max_iterations );
                }

                // Keep the last F_R value produced by stage 1, matching the RD
                // path and the reference ieuler semantics.
                eval_diffusion( t_diffusion_coupling, u_sp1, fd_tmp );
                eval_advection( t_advection_sp1, u_sp1, Fa_u_sp1 );

                auto& rhs_sp2 = U[17];
                rhs_sp2       = u_sm2pl + beta * dt * fd_tmp + dt * Fa_u_sp1 + ( 1. - 2. * gamma ) * dt * fr_tmp;

                auto g_sp2 = [&]( state_t& u ) -> state_t
                {
                    eval_reaction( t_sp2, u, fr_tmp_bis );
                    return u - gamma * dt * fr_tmp_bis - rhs_sp2;
                };

                if constexpr ( can_reuse_factorization )
                {
                    auto refresh_stage2_factorization = [&]( state_t const& u )
                    {
                        decltype( auto ) refreshed_reaction_jacobian = reaction.df( t_sp2, u );

                        matrix_t refreshed_stage_matrix;
                        refreshed_stage_matrix = matrix_linear_algebra_t::identity_minus( refreshed_reaction_jacobian, gamma * dt );

                        matrix_linear_algebra_t::factorize( refreshed_stage_matrix );
                    };

                    u_sp2 = diagonal_implicit_runge_kutta::pirock_newton<value_t, state_t, matrix_t>( g_sp2,
                        rhs_sp2,
                        refresh_stage2_factorization );
                }
                else
                {
                    auto identity = matrix_linear_algebra_t::identity( un );

                    auto dg_sp2 = [&]( state_t& u ) -> matrix_t
                    {
                        return identity - gamma * dt * reaction.df( t_sp2, u );
                    };

                    u_sp2 = diagonal_implicit_runge_kutta::newton<value_t>( g_sp2,
                        dg_sp2,
                        rhs_sp2,
                        matrix_linear_algebra_t::solver,
                        ponio::default_config::newton_tolerance,
                        ponio::default_config::newton_max_iterations );
                }
            }

            // Remaining advection-reaction finishing stages from Eq. (31).
            auto& u_sp3 = U[17];
            auto& u_sp4 = U[18];
            auto& u_sp5 = U[19];

            u_sp3 = u_sm2pl + ( 1. - 2. * gamma ) * dt * Fa_u_sp1 + ( 1. - gamma ) * dt * fr_tmp;

            u_sp4 = u_sm2pl + dt / 3. * Fa_u_sp1;

            eval_advection( t_advection_sp4, u_sp4, fa_tmp );

            if constexpr ( shampine_trick_enable )
            {
                auto& filtered_A = U[20];
                auto& u_tmp      = U[22];

                if constexpr ( detail::problem_operator<reaction_problem_t, value_t> )
                {
                    shampine_trick_caller.template operator()<1>( gamma * dt, reaction.f_t( tn ), u_sm2pl, fa_tmp, u_tmp, filtered_A );
                }
                else
                {
                    // Reuse the factorization left by the implicit reaction stages.
                    shampine_trick_caller.initialize();
                    shampine_trick_caller.template apply<1>( fa_tmp, u_tmp, filtered_A );
                    shampine_trick_caller.finalize();
                }

                u_sp5 = u_sm2pl + 2. / 3. * beta * dt * fd_tmp + 2. / 3. * dt * filtered_A + ( 2. / 3. - gamma ) * dt * fr_tmp
                      + 2. / 3. * gamma * dt * fr_tmp_bis;
            }
            else
            {
                u_sp5 = u_sm2pl + 2. / 3. * beta * dt * fd_tmp + 2. / 3. * dt * fa_tmp + ( 2. / 3. - gamma ) * dt * fr_tmp
                      + 2. / 3. * gamma * dt * fr_tmp_bis;
            }

            // Third and last advection evaluation, reused by the final formula
            // and the advection error estimator.
            eval_advection( t_advection_sp5, u_sp5, fa_tmp_bis );

            // Final diffusion coupling F_D(u^{s+3}) - F_D(u^{s+1}).
            eval_diffusion( t_diffusion_coupling, u_sp3, fd_tmp_bis );
            fd_tmp_bis = static_cast<state_t>( fd_tmp_bis - fd_tmp );

            if constexpr ( shampine_trick_enable )
            {
                auto& filtered_D = U[21];
                auto& u_tmp      = U[22];

                if constexpr ( is_embedded )
                {
                    auto& rhs_R = U[23];
                    auto& err_R = U[24];
                    auto& err_A = U[25];

                    rhs_R = static_cast<state_t>( dt / 6. * ( fr_tmp_bis - fr_tmp ) );

                    if constexpr ( detail::problem_operator<reaction_problem_t, value_t> )
                    {
                        shampine_trick_caller.template operator()<l>( gamma * dt, reaction.f_t( tn ), u_sm2pl, fd_tmp_bis, u_tmp, filtered_D );

                        shampine_trick_caller.template operator()<1>( gamma * dt, reaction.f_t( tn ), u_sm2pl, rhs_R, u_tmp, err_R );
                    }
                    else
                    {
                        shampine_trick_caller.initialize();
                        shampine_trick_caller.template apply<l>( fd_tmp_bis, u_tmp, filtered_D );
                        shampine_trick_caller.template apply<1>( rhs_R, u_tmp, err_R );
                        shampine_trick_caller.finalize();
                    }

                    err_A = -0.15 * dt * Fa_u_sp1 + 0.3 * dt * fa_tmp - 0.15 * dt * fa_tmp_bis;

                    state_algebra_t::add( compensation, 0.5 * dt * fr_tmp );
                    state_algebra_t::add( compensation, 0.25 * dt * Fa_u_sp1 );
                    state_algebra_t::add( compensation, 0.5 * dt * fr_tmp_bis );
                    state_algebra_t::add( compensation, dt / ( 2. - 4. * gamma ) * filtered_D );
                    state_algebra_t::add( compensation, 0.75 * dt * fa_tmp_bis );
                    compensated_commit( us_s, u_np1 );

                    value_t const err_R_scalar = detail::error_algebra<state_t>::estimate_squared( err_R,
                        un,
                        u_np1,
                        _info.absolute_tolerance,
                        _info.relative_tolerance );

                    value_t const err_A_scalar = detail::error_algebra<state_t>::estimate_squared( err_A,
                        un,
                        u_np1,
                        _info.absolute_tolerance,
                        _info.relative_tolerance );

                    // The advection estimator has order three. In RMS^2 form,
                    // ||err_A||^(2/3) becomes err_A_scalar^(2/3).
                    value_t const err_A_effective_scalar = std::pow( err_A_scalar, 2. / 3. );

                    _info.error   = std::max( { err_D_scalar, err_A_effective_scalar, err_R_scalar } );
                    _info.success = _info.error < 1.0;

                    value_t fac = raw_fac_from_current_error();

                    // Apply the controller memory correction only when diffusion
                    // dominates the other error components. The factor 100 is the
                    // squared-RMS counterpart of the Fortran factor 10.
                    if ( errp != 0.0 && !previous_attempt_rejected && err_D_scalar >= 100.0 * std::max( err_A_effective_scalar, err_R_scalar ) )
                    {
                        value_t const facp = std::pow( errp, 0.25 ) * fac * fac * ( dt / hp );
                        fac                = std::min( fac, facp );
                    }

                    fac            = finalize_fac( fac );
                    value_t new_dt = dt * fac;

                    if ( _info.success )
                    {
                        compensated_time_update( dt );

                        errp = _info.error;
                        hp   = dt;

                        if ( previous_attempt_rejected )
                        {
                            new_dt = std::min( new_dt, dt );
                            nrej   = 0;
                        }

                        dt = new_dt;

                        facmax                    = 2.0;
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
                    if constexpr ( detail::problem_operator<reaction_problem_t, value_t> )
                    {
                        shampine_trick_caller.template operator()<l>( gamma * dt, reaction.f_t( tn ), u_sm2pl, fd_tmp_bis, u_tmp, filtered_D );
                    }
                    else
                    {
                        shampine_trick_caller.initialize();
                        shampine_trick_caller.template apply<l>( fd_tmp_bis, u_tmp, filtered_D );
                        shampine_trick_caller.finalize();
                    }

                    compensated_time_update( dt );

                    state_algebra_t::add( compensation, 0.5 * dt * fr_tmp );
                    state_algebra_t::add( compensation, 0.25 * dt * Fa_u_sp1 );
                    state_algebra_t::add( compensation, 0.5 * dt * fr_tmp_bis );
                    state_algebra_t::add( compensation, dt / ( 2. - 4. * gamma ) * filtered_D );
                    state_algebra_t::add( compensation, 0.75 * dt * fa_tmp_bis );
                    compensated_commit( us_s, u_np1 );
                }
            }
            else
            {
                // Fixed-step variant without Shampine filtering.
                compensated_time_update( dt );

                state_algebra_t::add( compensation, 0.5 * dt * fr_tmp );
                state_algebra_t::add( compensation, 0.25 * dt * Fa_u_sp1 );
                state_algebra_t::add( compensation, 0.5 * dt * fr_tmp_bis );
                state_algebra_t::add( compensation, dt / ( 2. - 4. * gamma ) * fd_tmp_bis );
                state_algebra_t::add( compensation, 0.75 * dt * fa_tmp_bis );
                compensated_commit( us_s, u_np1 );
            }
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

        /**
         * @brief set absolute tolerance in chained config
         *
         * @param tol_ tolerance
         * @return auto& returns this object
         */
        template <typename rock_t = rock_coeff>
            requires std::same_as<rock_t, rock_coeff> && is_embedded
        auto&
        abs_tol( value_t tol_ )
        {
            info().absolute_tolerance = tol_;
            return *this;
        }

        /**
         * @brief set relative tolerance in chained config
         *
         * @param tol_ tolerance
         * @return auto& returns this object
         */
        template <typename rock_t = rock_coeff>
            requires std::same_as<rock_t, rock_coeff> && is_embedded
        auto&
        rel_tol( value_t tol_ )
        {
            info().relative_tolerance = tol_;
            return *this;
        }
    };

    // cppcheck-suppress-begin unusedFunction

    /**
     * @brief helper function to build PIROCK algorithm
     *
     * @tparam l                       Number of augmented stages (l = 1, 2)
     * @tparam is_embedded             Set adaptive time step method (default: false)
     * @tparam value_t                 Type of coefficients
     * @tparam alpha_beta_computer_t   Choice of computing of \f$\alpha\f$ and \f$\beta\f$ parameters
     * @tparam eig_computer_t          Computing of eigenvalues of explicit part (diffusion)
     * @tparam shampine_trick_caller_t Computing of Shampine's trick
     * @param alpha_beta_computer      \f$\alpha\f$ and \f$\beta\f$ computer object
     * @param eig_computer             Eigenvalue computer of explicit part of the problem
     * @param shampine_trick_caller    Shampine's trick computer
     */
    template <std::size_t l = 1, bool is_embedded = false, typename value_t = double, typename alpha_beta_computer_t, typename eig_computer_t, typename shampine_trick_caller_t>
    auto
    pirock_RDA( alpha_beta_computer_t&& alpha_beta_computer, eig_computer_t&& eig_computer, shampine_trick_caller_t&& shampine_trick_caller )
    {
        return pirock_RDA_impl<l, alpha_beta_computer_t, eig_computer_t, shampine_trick_caller_t, is_embedded, value_t>(
            std::forward<alpha_beta_computer_t>( alpha_beta_computer ),
            std::forward<eig_computer_t>( eig_computer ),
            std::forward<shampine_trick_caller_t>( shampine_trick_caller ) );
    }

    /**
     * @brief helper function to build PIROCK algorithm
     *
     * @tparam l                     Number of augmented stages (l = 1, 2)
     * @tparam value_t               Type of coefficients
     * @tparam alpha_beta_computer_t Choice of computing of \f$\alpha\f$ and \f$\beta\f$ parameters
     * @tparam eig_computer_t        Computing of eigenvalues of explicit part (diffusion)
     * @param alpha_beta_computer    \f$\alpha\f$ and \f$\beta\f$ computer object
     * @param eig_computer           Eigenvalue computer of explicit part of the problem
     *
     * @note Without a Shampine's trick caller, this method is fixed time step.
     */
    template <std::size_t l = 1, typename value_t = double, typename alpha_beta_computer_t, typename eig_computer_t>
    auto
    pirock_RDA( alpha_beta_computer_t&& alpha_beta_computer, eig_computer_t&& eig_computer )
    {
        return pirock_RDA_impl<l, alpha_beta_computer_t, eig_computer_t, void, false, value_t>(
            std::forward<alpha_beta_computer_t>( alpha_beta_computer ),
            std::forward<eig_computer_t>( eig_computer ) );
    }

    /**
     * @brief helper function to build PIROCK algorithm
     *
     * @tparam l              Number of augmented stages (l = 1, 2)
     * @tparam value_t        Type of coefficients
     * @tparam eig_computer_t Computing of eigenvalues of explicit part (diffusion)
     * @param eig_computer    Eigenvalue computer of explicit part of the problem
     *
     * @note Without a Shampine's trick caller, this method is fixed time step, and without a \f$\alpha\f$ and \f$\beta\f$ computer,
     * parameters are fixed to \f$\beta = 0\f$.
     */
    template <std::size_t l = 1, typename value_t = double, typename eig_computer_t>
    auto
    pirock_RDA( eig_computer_t&& eig_computer )
    {
        return pirock_RDA<l, value_t>( beta_0<value_t>(), std::forward<eig_computer_t>( eig_computer ) );
    }

    /**
     * @brief helper function to build PIROCK algorithm
     *
     * @tparam l       Number of augmented stages (l = 1, 2)
     * @tparam value_t Type of coefficients
     *
     * @note Without a Shampine's trick caller, this method is fixed time step, without a \f$\alpha\f$ and \f$\beta\f$ computer, parameters
     * are fixed to \f$\beta = 0\f$ and without a eigenvalues computer this method use power method to estimate spectral radius of explicit
     * part (diffusion operator).
     */
    template <std::size_t l = 1, typename value_t = double>
    auto
    pirock_RDA()
    {
        return pirock_RDA<l, value_t>( rock::detail::power_method() );
    }

    /**
     * @brief helper function to build PIROCK algorithm
     *
     * @tparam value_t        Type of coefficients
     * @tparam eig_computer_t Computing of eigenvalues of explicit part (diffusion)
     * @param eig_computer    Eigenvalue computer of explicit part of the problem
     *
     * @note This function builds a PIROCK algorithm with parameters \f$\ell=2\f$ and \f$\alpha = 1\f$ and specific spectral radius
     * estimator of explicit part (diffusion operator).
     */
    template <typename value_t = double, typename eig_computer_t>
    auto
    pirock_RDA_a1( eig_computer_t&& eig_computer )
    {
        return pirock_RDA<2, value_t>( alpha_fixed<value_t>( 1.0 ), std::forward<eig_computer_t>( eig_computer ) );
    }

    /**
     * @brief helper function to build PIROCK algorithm
     *
     * @tparam value_t Type of coefficients
     *
     * @note This function builds a PIROCK algorithm with parameters \f$\ell=2\f$ and \f$\alpha = 1\f$ and power method to estimate spectral
     * radius.
     */
    template <typename value_t = double>
    auto
    pirock_RDA_a1()
    {
        return pirock_RDA_a1<value_t>( rock::detail::power_method() );
    }

    /**
     * @brief helper function to build PIROCK algorithm
     *
     * @tparam value_t        Type of coefficients
     * @tparam eig_computer_t Computing of eigenvalues of explicit part (diffusion)
     * @param eig_computer    Eigenvalue computer of explicit part of the problem
     *
     * @note This function builds a PIROCK algorithm with parameters \f$\ell=1\f$ and \f$\beta = 0\f$ and specific spectral radius estimator
     * of explicit part (diffusion operator).
     */
    template <typename value_t = double, typename eig_computer_t>
    auto
    pirock_RDA_b0( eig_computer_t&& eig_computer )
    {
        return pirock_RDA<1, value_t>( beta_0<value_t>(), std::forward<eig_computer_t>( eig_computer ) );
    }

    /**
     * @brief helper function to build PIROCK algorithm
     *
     * @tparam value_t Type of coefficients
     *
     * @note This function builds a PIROCK algorithm with parameters \f$\ell=1\f$ and \f$\beta = 0\f$ and power method to estimate spectral
     * radius.
     */
    template <typename value_t = double>
    auto
    pirock_RDA_b0()
    {
        return pirock_RDA_b0<value_t>( rock::detail::power_method() );
    }

    // cppcheck-suppress-end unusedFunction

} // namespace ponio::runge_kutta::pirock
