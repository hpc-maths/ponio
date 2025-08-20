// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private, include "../runge_kutta.hpp"

#pragma once

// NOLINTBEGIN(misc-include-cleaner)

#include <array>
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
                double const mu    = rock_coeff::recf[mr + 2 * ( j - 2 ) + 1 - 1];
                double const kappa = rock_coeff::recf[mr + 2 * ( j - 2 ) + 2 - 1];
                double const nu    = -1.0 - kappa;

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

        static constexpr std::size_t l           = _l;
        static constexpr bool is_imex_method     = true;
        static constexpr std::size_t N_operators = 2;
        static constexpr std::size_t N_stages    = stages::dynamic;
        // clang-format off
        static constexpr std::size_t N_storage   = std::conditional_t<shampine_trick_enable,
                                                    std::integral_constant<std::size_t, 18>, // if shampine's trick
                                                    std::integral_constant<std::size_t, 12>  // else (no error estimation)
                                                >::value;
        // clang-format on
        static constexpr std::size_t order   = 2;
        static constexpr std::string_view id = "PIROCK";

        using value_t                 = _value_t;
        using rock_coeff              = rock::rock2_coeff<value_t>;
        using degree_computer         = rock::detail::degree_computer<value_t, rock_coeff>;
        using shampine_trick_caller_t = typename std::conditional_t<shampine_trick_enable, _shampine_trick_caller_t, bool>;

        alpha_beta_computer_t alpha_beta_computer;
        eig_computer_t eig_computer;
        shampine_trick_caller_t shampine_trick_caller;

        iteration_info<pirock_impl> _info;

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
         * @param a_tol                  absolute tolerance
         * @param r_tol                  relative tolerance
         */
        template <typename _shampine_trick_caller_t_>
            requires std::same_as<_shampine_trick_caller_t_, shampine_trick_caller_t>
                      && std::same_as<std::bool_constant<shampine_trick_enable>, std::true_type>
        pirock_impl( alpha_beta_computer_t&& _alpha_beta_computer,
            eig_computer_t&& _eig_computer,
            _shampine_trick_caller_t_&& _shampine_trick_caller,
            value_t a_tol = default_config::tol,
            value_t r_tol = default_config::tol )
            : alpha_beta_computer( std::forward<alpha_beta_computer_t>( _alpha_beta_computer ) )
            , eig_computer( std::forward<eig_computer_t>( _eig_computer ) )
            , shampine_trick_caller( std::forward<_shampine_trick_caller_t_>( _shampine_trick_caller ) )
            , _info( a_tol, r_tol )
        {
        }

        /**
         * @brief iteration of PIROCK method
         *
         * @tparam problem_t  type of \f$f\f$
         * @tparam state_t    type of current state
         * @tparam array_ki_t type of temporary stages (only 3 needed for ROCK2)
         * @param pb    problem \f$(F_D, F_R)\f$ and the Jacibian of reaction part \f$\frac{\partial F_R}{\partial u}\f$
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
            // | 13    | f_D_u            | $F_D(u^{\star(s-1)}) - F_D(u^{(s-2)})$ term       |
            // | 14    | u_tmp            | temporary value for compute Shampine's trick      |
            // | 15    | fe_tmp_bis       | temporary output of explicit part                 |
            // | 15    | err_D            | $err_D$ error on diffusion term                   |
            // | 16    | rhs_R            | right-hand side to compute $err_R$ error term     |
            // | 17    | err_R            | $err_R$ error on reaction term                    |
            //
            // > if method is called as a constant time step method, only index from 0 to 12 are used

            _info.reset_eval();

            auto [mdeg, deg_index, start_index, n_eval] = degree_computer::compute_n_stages_optimal_degree( rock::rock_order::rock_2(),
                eig_computer,
                pb.explicit_part,
                tn,
                un,
                dt,
                U,
                4 );

            std::size_t s = mdeg + 2;

            _info.number_of_stages  = s + l + 3;
            _info.number_of_eval[0] = n_eval + s + l + 4; // explicit evaluation

            value_t const alpha = alpha_beta_computer.alpha( s, l );
            value_t const beta  = alpha_beta_computer.beta( s, l );
            value_t const gamma = 1. - 0.5 * std::numbers::sqrt2;

            auto& u_j    = U[0];
            auto& u_jm1  = U[1];
            auto& u_jm2  = U[2];
            auto& u_sm2  = U[3];
            auto& fe_tmp = U[4];
            auto& fi_tmp = U[5];
            auto& f_tmp  = U[6];

            u_j   = un;
            u_jm2 = un;

            value_t const mu_1 = rock_coeff::recf[start_index - 1];

            value_t t_jm1 = tn + dt * alpha * mu_1;
            value_t t_jm2 = tn + dt * alpha * mu_1;
            value_t t_jm3 = tn;

            // u_1 =u^n + \alpha \mu_1 \Delta F_D( u^n )
            pb.explicit_part( tn, un, fe_tmp );
            u_jm1 = un + alpha * dt * mu_1 * fe_tmp;

            if ( mdeg < 2 )
            {
                u_j = u_jm1;
            }

            for ( std::size_t j = 2; j < s - 2 + l + 1; ++j )
            {
                value_t const mu_j    = rock_coeff::recf[start_index + 2 * ( j - 2 ) + 1 - 1];
                value_t const kappa_j = rock_coeff::recf[start_index + 2 * ( j - 2 ) + 2 - 1];
                value_t const nu_j    = -1.0 - kappa_j;

                // u_{j} = \alpha \mu_j \Delta t F_D( u_{j-2} ) - \nu_j u_{j-1} - \kappa_j u_{j-2}
                pb.explicit_part( t_jm1, u_jm1, fe_tmp );
                u_j = alpha * mu_j * dt * fe_tmp - nu_j * u_jm1 - kappa_j * u_jm2;

                t_jm1 = alpha * dt * mu_j - nu_j * t_jm2 - kappa_j * t_jm3;

                if ( j == s - 2 )
                {
                    u_sm2 = u_j;
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

            // u_{*s-1} = u_{s-2} + \sigma_\alpha \Delta t  F_D( u_{s-2} )
            auto& us_sm1 = U[7];
            pb.explicit_part( t_jm1, u_sm2, fe_tmp );
            us_sm1 = u_sm2 + sigma_a * dt * fe_tmp;

            // u_{*s} = u_{*s-1} + \sigma_\alpha \Delta t  F_D( u_{*s-1} )
            auto& us_s = U[8];
            pb.explicit_part( t_jm1, us_sm1, fe_tmp );
            us_s = us_sm1 + sigma_a * dt * fe_tmp;

            // u_{s-2+l} = u_j
            auto& u_sm2pl = u_j;

            auto& u_sp1 = U[9];
            u_sp1       = un;

            auto& u_sp2 = U[10];
            u_sp2       = un;

            if constexpr ( detail::problem_operator<decltype( pb.implicit_part ), value_t> )
            {
                std::size_t n_eval_sp1 = 0;

                auto op_sp1  = ::ponio::linear_algebra::operator_algebra<state_t>::identity( un ) - gamma * dt * pb.implicit_part.f_t( tn );
                auto rhs_sp1 = u_sm2pl;
                ::ponio::linear_algebra::operator_algebra<state_t>::solve( op_sp1, u_sp1, rhs_sp1, n_eval_sp1 );

                std::size_t n_eval_sp2 = 0;

                pb.explicit_part( tn, u_sp1, fe_tmp );
                pb.implicit_part( tn, u_sp1, fi_tmp );
                auto& rhs_sp2 = U[11]; // temporary use of U[10] before u_sp3

                auto op_sp2 = ::ponio::linear_algebra::operator_algebra<state_t>::identity( un ) - gamma * dt * pb.implicit_part.f_t( tn );
                rhs_sp2     = u_sm2pl + beta * dt * fe_tmp + ( 1. - 2. * gamma ) * dt * fi_tmp;
                ::ponio::linear_algebra::operator_algebra<state_t>::solve( op_sp2, u_sp2, rhs_sp2, n_eval_sp2 );

                _info.number_of_eval[1] += n_eval_sp1 + n_eval_sp2 + 1;
            }
            else
            {
                using matrix_t = decltype( pb.implicit_part.df( tn, un ) );

                auto identity = ::ponio::linear_algebra::linear_algebra<matrix_t>::identity( un );
                auto g_sp1    = [&]( state_t& u ) -> state_t
                {
                    _info.number_of_eval[1] += 1;
                    pb.implicit_part( tn, u, fi_tmp );
                    return u - gamma * dt * fi_tmp - u_sm2pl;
                };
                auto dg = [&]( state_t const& u ) -> matrix_t
                {
                    return identity - gamma * dt * pb.implicit_part.df( tn, u );
                };
                u_sp1 = diagonal_implicit_runge_kutta::newton<value_t>( g_sp1,
                    dg,
                    u_sm2pl,
                    ::ponio::linear_algebra::linear_algebra<matrix_t>::solver,
                    ponio::default_config::newton_tolerance,
                    ponio::default_config::newton_max_iterations );

                _info.number_of_eval[0] += 1;
                _info.number_of_eval[1] += 1;
                pb.explicit_part( tn, u_sp1, fe_tmp );
                pb.implicit_part( tn, u_sp1, fi_tmp );

                auto g_sp2 = [&]( state_t& u ) -> state_t
                {
                    _info.number_of_eval[1] += 1;

                    pb.implicit_part( tn, u, f_tmp );
                    return u - gamma * dt * f_tmp - ( u_sm2pl + beta * dt * fe_tmp + ( 1. - 2. * gamma ) * dt * fi_tmp );
                };
                u_sp2 = diagonal_implicit_runge_kutta::newton<value_t>( g_sp2,
                    dg,
                    u_sm2pl,
                    ::ponio::linear_algebra::linear_algebra<matrix_t>::solver,
                    ponio::default_config::newton_tolerance,
                    ponio::default_config::newton_max_iterations );
            }

            _info.number_of_eval[1] += 3;

            auto& u_sp3 = U[11];
            pb.implicit_part( tn, u_sp1, fi_tmp );
            u_sp3 = u_sm2pl + ( 1. - gamma ) * dt * fi_tmp;

            value_t tau   = sigma * rock_coeff::fp2[deg_index - 1] + sigma * sigma;
            value_t tau_a = 0.5 * detail::power<2>( alpha - 1. ) + 2. * alpha * ( 1. - alpha ) * sigma + alpha * alpha * tau;

            // auto& u_np1 = U[12];

            if constexpr ( shampine_trick_enable && detail::problem_operator<decltype( pb.implicit_part ), value_t> )
            {
                auto& shampine_element = U[12];
                auto& f_D_u            = U[13];
                auto& u_tmp            = U[14];
                auto& fe_tmp_bis       = U[14]; // temporary use of U[15] before u_tmp

                // for embedded method
                auto& err_D = U[15];

                // $err_D = \sigma_\alpha(1-\tau_a/\sigma_a^2)\Delta t (F_D(u^{*(s-1)}) - F_D(u^{(s-2)}))$
                pb.explicit_part( tn, us_sm1, fe_tmp );
                pb.explicit_part( tn, u_sm2, fe_tmp_bis );
                err_D = sigma_a * ( 1. - tau_a / ( sigma_a * sigma_a ) ) * dt * ( fe_tmp - fe_tmp_bis );

                pb.explicit_part( tn, u_sp3, fe_tmp );
                pb.explicit_part( tn, u_sp1, fe_tmp_bis );
                f_D_u = static_cast<state_t>( fe_tmp - fe_tmp_bis );

                shampine_trick_caller.template operator()<l>( gamma * dt, pb.implicit_part.f_t( tn ), u_sm2pl, f_D_u, u_tmp, shampine_element );

                if constexpr ( is_embedded )
                {
                    auto& rhs_R = U[16];
                    auto& err_R = U[17];

                    _info.number_of_eval[1] += 2;

                    pb.implicit_part( tn, u_sp1, fi_tmp );
                    pb.implicit_part( tn, u_sp2, f_tmp );
                    rhs_R = static_cast<state_t>( dt / 6. * ( fi_tmp - f_tmp ) );

                    // $err_R = J_R^{-1} \Delta t/6 (F_R(u^{s+1}) - F_R(u^{s+2}))$
                    // to compute it, get $rhs_R = \Delta t/6 (F_R(u^{s+1}) - F_R(u^{s+2}))$
                    // then solve $J_R err_R = rhs_R$ (that what Shampine's trick does, it build $J_R$ and solve it)
                    shampine_trick_caller.template operator()<1>( gamma * dt, pb.implicit_part.f_t( tn ), u_sm2pl, rhs_R, u_tmp, err_R );

                    u_np1 = us_s - err_D + 0.5 * dt * fi_tmp + 0.5 * dt * f_tmp + dt / ( 2. - 4. * gamma ) * shampine_element;

                    // auto accumulator_error_gen = []( auto const& yn, auto const& ynp1, value_t a_tol, value_t r_tol )
                    // {
                    //     return [=, it_yn = yn.begin(), it_ynp1 = ynp1.begin()]( value_t const& acc, value_t const err_i ) mutable
                    //     {
                    //         return acc
                    //              + detail::power<2>( err_i / ( a_tol + r_tol * std::max( std::abs( *it_yn++ ), std::abs( *it_ynp1++ ) ) )
                    //              );
                    //     };
                    // };

                    // TODO: this couple of lines works only with samurai (because of err_D.array())
                    // value_t err_R_scalar = std::accumulate( err_R.array().begin(),
                    //     err_R.array().end(),
                    //     static_cast<value_t>( 0. ),
                    //     accumulator_error_gen( un.array(), u_np1.array(), _info.absolute_tolerance, _info.relative_tolerance ) );
                    // value_t err_D_scalar = std::accumulate( err_D.array().begin(),
                    //     err_D.array().end(),
                    //     static_cast<value_t>( 0. ),
                    //     accumulator_error_gen( un.array(), u_np1.array(), _info.absolute_tolerance, _info.relative_tolerance ) );

                    value_t err_R_scalar = shampine_trick_caller.norm_error( err_R,
                        un,
                        u_np1,
                        _info.absolute_tolerance,
                        _info.relative_tolerance );
                    value_t err_D_scalar = shampine_trick_caller.norm_error( err_D,
                        un,
                        u_np1,
                        _info.absolute_tolerance,
                        _info.relative_tolerance );

                    _info.error   = std::max( err_D_scalar, err_R_scalar );
                    _info.success = _info.error < 1.0;

                    // std::cout << "tn " << tn << " dt " << dt << " mdeg " << mdeg << "\n";
                    // std::cout << "err_D " << err_D_scalar << " err_R " << err_R_scalar << "\n";
                    // std::cout << "a_tol " << _info.absolute_tolerance << " r_tol " << _info.relative_tolerance << "\n";

                    value_t fac    = std::min( 2.0, std::max( 0.5, std::sqrt( 1.0 / _info.error ) ) );
                    value_t new_dt = 0.8 * fac * dt;

                    // accepted step
                    if ( _info.success )
                    {
                        tn = tn + dt;
                        dt = new_dt;

                        // return { tn + dt, u_np1, new_dt };
                    }
                    else
                    {
                        // tn = tn
                        std::swap( un, u_np1 );
                        dt = new_dt;

                        // return { tn, un, new_dt };
                    }
                }
                else
                {
                    pb.implicit_part( tn, u_sp1, fi_tmp );
                    pb.implicit_part( tn, u_sp2, f_tmp );

                    tn    = tn + dt;
                    u_np1 = us_s - err_D + 0.5 * dt * fi_tmp + 0.5 * dt * f_tmp + dt / ( 2. - 4. * gamma ) * shampine_element;
                }
            }
            else
            {
                auto& fe_tmp_bis = U[0]; // temporary use of U[0] after u_j
                auto& fe_tmp_ter = U[1]; // temporary use of U[1] after u_jm1
                auto& fe_tmp_qua = U[2]; // temporary use of U[2] after u_jm2

                pb.explicit_part( tn, us_sm1, fe_tmp );
                pb.explicit_part( tn, u_sm2, fe_tmp_bis );
                pb.explicit_part( tn, u_sp3, fe_tmp_ter );
                pb.explicit_part( tn, u_sp1, fe_tmp_qua );

                pb.implicit_part( tn, u_sp1, fi_tmp );
                pb.implicit_part( tn, u_sp2, f_tmp );

                tn    = tn + dt;
                u_np1 = us_s - sigma_a * ( 1. - tau_a / ( sigma_a * sigma_a ) ) * dt * ( fe_tmp - fe_tmp_bis ) + 0.5 * dt * fi_tmp
                      + 0.5 * dt * f_tmp + dt / ( 2. - 4. * gamma ) * ( fe_tmp_ter - fe_tmp_qua );
            }

            // return { tn + dt, u_np1, dt };
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
     * @param absolute_tolerance       Absolute tolerance for adaptive time step method (default: ponio::default_config::tol)
     * @param relative_tolerance       Relative tolerance for adaptive time step method (default: ponio::default_config::tol)
     */
    template <std::size_t l = 1, bool is_embedded = false, typename value_t = double, typename alpha_beta_computer_t, typename eig_computer_t, typename shampine_trick_caller_t>
    auto
    pirock( alpha_beta_computer_t&& alpha_beta_computer,
        eig_computer_t&& eig_computer,
        shampine_trick_caller_t&& shampine_trick_caller,
        value_t absolute_tolerance = default_config::tol,
        value_t relative_tolerance = default_config::tol )
    {
        return pirock_impl<l, alpha_beta_computer_t, eig_computer_t, shampine_trick_caller_t, is_embedded, value_t>(
            std::forward<alpha_beta_computer_t>( alpha_beta_computer ),
            std::forward<eig_computer_t>( eig_computer ),
            std::forward<shampine_trick_caller_t>( shampine_trick_caller ),
            absolute_tolerance,
            relative_tolerance );
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

        static constexpr std::size_t l           = _l;
        static constexpr bool is_imex_method     = true;
        static constexpr std::size_t N_operators = 3;
        static constexpr std::size_t N_stages    = stages::dynamic;
        static constexpr std::size_t N_storage   = std::
            conditional_t<shampine_trick_enable, std::integral_constant<std::size_t, 27>, std::integral_constant<std::size_t, 20>>::value;
        static constexpr std::size_t order   = 2;
        static constexpr std::string_view id = "PIROCK";

        using value_t                 = _value_t;
        using rock_coeff              = rock::rock2_coeff<value_t>;
        using degree_computer         = rock::detail::degree_computer<value_t, rock_coeff>;
        using shampine_trick_caller_t = typename std::conditional_t<shampine_trick_enable, _shampine_trick_caller_t, bool>;

        alpha_beta_computer_t alpha_beta_computer;
        eig_computer_t eig_computer;
        shampine_trick_caller_t shampine_trick_caller;

        iteration_info<pirock_RDA_impl> _info;

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
         * @param a_tol                  absolute tolerance
         * @param r_tol                  relative tolerance
         */
        template <typename _shampine_trick_caller_t_>
            requires std::same_as<_shampine_trick_caller_t_, shampine_trick_caller_t>
                      && std::same_as<std::bool_constant<shampine_trick_enable>, std::true_type>
        pirock_RDA_impl( alpha_beta_computer_t&& _alpha_beta_computer,
            eig_computer_t&& _eig_computer,
            _shampine_trick_caller_t_&& _shampine_trick_caller,
            value_t a_tol = default_config::tol,
            value_t r_tol = default_config::tol )
            : alpha_beta_computer( std::forward<alpha_beta_computer_t>( _alpha_beta_computer ) )
            , eig_computer( std::forward<eig_computer_t>( _eig_computer ) )
            , shampine_trick_caller( std::forward<_shampine_trick_caller_t_>( _shampine_trick_caller ) )
            , _info( a_tol, r_tol )
        {
        }

        /**
         * @brief iteration of PIROCK_RDA method
         *
         * @tparam problem_t  type of \f$f\f$
         * @tparam state_t    type of current state
         * @tparam array_ki_t type of temporary stages (13 or 18 stages needed)
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
            // In problem_t pb.system (a tuple) :
            // 0: reaction operator
            // 1: diffusion operator
            // 2: advection operator
            using reaction_op  = std::integral_constant<std::size_t, 0>;
            using diffusion_op = std::integral_constant<std::size_t, 1>;
            using advection_op = std::integral_constant<std::size_t, 2>;

            static_assert( detail::problem_operator<std::tuple_element_t<reaction_op::value, decltype( pb.system )>, value_t>
                               || detail::problem_jacobian<std::tuple_element_t<reaction_op::value, decltype( pb.system )>, value_t, state_t>,
                "This kind of problem is not inversible in ponio" );

            _info.reset_eval();

            auto [mdeg, deg_index, start_index, n_eval] = degree_computer::compute_n_stages_optimal_degree( rock::rock_order::rock_2(),
                eig_computer,
                std::get<diffusion_op::value>( pb.system ),
                tn,
                un,
                dt,
                U,
                4 );
            std::size_t s                               = mdeg + 2;

            _info.number_of_stages  = s + l + 5;
            _info.number_of_eval[0] = n_eval + s + l + 4; // explicit evaluation

            value_t const alpha = alpha_beta_computer.alpha( s, l );
            value_t const beta  = alpha_beta_computer.beta( s, l );
            value_t const gamma = 1. - 0.5 * std::numbers::sqrt2;

            auto& u_j   = U[0];
            auto& u_jm1 = U[1];
            auto& u_jm2 = U[2];
            auto& u_sm2 = U[3];

            auto& fr_tmp     = U[4];
            auto& fr_tmp_bis = U[5];
            auto& fd_tmp     = U[6];
            auto& fd_tmp_bis = U[7];
            auto& fa_tmp     = U[10];
            auto& fa_tmp_bis = U[11];

            u_j   = un;
            u_jm2 = un;

            value_t const mu_1 = rock_coeff::recf[start_index - 1];

            value_t t_jm1 = tn + dt * alpha * mu_1;
            value_t t_jm2 = tn + dt * alpha * mu_1;
            value_t t_jm3 = tn;

            // u_1 =u^n + \alpha \mu_1 \Delta F_D( u^n )
            pb( diffusion_op(), tn, un, fd_tmp );
            u_jm1 = un + alpha * dt * mu_1 * fd_tmp;

            if ( mdeg < 2 )
            {
                u_j = u_jm1;
            }

            for ( std::size_t j = 2; j < s - 2 + l + 1; ++j )
            {
                value_t const mu_j    = rock_coeff::recf[start_index + 2 * ( j - 2 ) + 1 - 1];
                value_t const kappa_j = rock_coeff::recf[start_index + 2 * ( j - 2 ) + 2 - 1];
                value_t const nu_j    = -1.0 - kappa_j;

                // u_{j} = \alpha \mu_j \Delta t F_D( u_{j-2} ) - \nu_j u_{j-1} - \kappa_j u_{j-2}
                pb( diffusion_op(), t_jm1, u_jm1, fd_tmp );
                u_j = alpha * mu_j * dt * fd_tmp - nu_j * u_jm1 - kappa_j * u_jm2;

                t_jm1 = alpha * dt * mu_j - nu_j * t_jm2 - kappa_j * t_jm3;

                if ( j == s - 2 )
                {
                    u_sm2 = u_j;
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

            // u_{*s-1} = u_{s-2} + \sigma_\alpha \Delta t  F_D( u_{s-2} )
            auto& us_sm1 = U[12];
            pb( diffusion_op(), t_jm1, u_sm2, fd_tmp );
            us_sm1 = u_sm2 + sigma_a * dt * fd_tmp;

            // u_{*s} = u_{*s-1} + \sigma_\alpha \Delta t  F_D( u_{*s-1} )
            auto& us_s = U[13];
            pb( diffusion_op(), t_jm1, us_sm1, fd_tmp );
            us_s = us_sm1 + sigma_a * dt * fd_tmp;

            // u_{s-2+l} = u_j
            auto& u_sm2pl = u_j;

            // u^{(s+1)} = u^{(s-2+\ell)} + \gamma \Delta t F_R(u^{(s+1)})
            auto& u_sp1 = U[14];
            u_sp1       = un;

            // u^{(s+2)} = u^{(s-2+\ell)} + \beta \Delta t F_D(u^{(s+1)}) + \Delta t F_A(u^{(s+1)}) + (1-2\gamma)\Delta t F_R(u^{(s+1)}) +
            // \gamma \Delta t F_R(u^{(s+2)})
            auto& u_sp2 = U[15];
            u_sp2       = un;

            // store F_A(u^{(s+1)})
            auto& Fa_u_sp1 = U[16];

            if constexpr ( detail::problem_operator<std::tuple_element_t<reaction_op::value, decltype( pb.system )>, value_t> )
            {
                std::size_t n_eval_sp1 = 0;

                // I - \gamma \Delta t F_R
                auto op_sp1 = ::ponio::linear_algebra::operator_algebra<state_t>::identity( un )
                            - gamma * dt * std::get<reaction_op::value>( pb.system ).f_t( tn );
                // u^{(s-2+\ell)}
                auto rhs_sp1 = u_sm2pl;

                ::ponio::linear_algebra::operator_algebra<state_t>::solve( op_sp1, u_sp1, rhs_sp1, n_eval_sp1 );

                pb( advection_op(), tn, u_sp1, Fa_u_sp1 );

                std::size_t n_eval_sp2 = 0;

                pb( diffusion_op(), tn, u_sp1, fd_tmp );
                pb( reaction_op(), tn, u_sp1, fr_tmp );
                // I - \gamma \Delta t F_R
                auto op_sp2 = ::ponio::linear_algebra::operator_algebra<state_t>::identity( un )
                            - gamma * dt * std::get<reaction_op::value>( pb.system ).f_t( tn );
                // u^{(s-2+\ell)} + \beta \Delta t F_D(u^{(s+1)}) + \Delta t F_A(u^{(s+1)}) + (1-2\gamma)\Delta t F_R(u^{(s+1)})
                auto rhs_sp2 = static_cast<state_t>( u_sm2pl + beta * dt * fd_tmp + dt * Fa_u_sp1 + ( 1. - 2. * gamma ) * dt * fr_tmp );

                ::ponio::linear_algebra::operator_algebra<state_t>::solve( op_sp2, u_sp2, rhs_sp2, n_eval_sp2 );

                _info.number_of_eval[1] += n_eval_sp1 + n_eval_sp2 + 1;
            }
            else
            {
                using matrix_t = decltype( std::get<reaction_op::value>( pb.system ).df( tn, un ) );

                auto identity = ::ponio::linear_algebra::linear_algebra<matrix_t>::identity( un );
                auto g_sp1    = [&]( state_t& u ) -> state_t
                {
                    _info.number_of_eval[1] += 1;
                    pb( reaction_op(), tn, u, fr_tmp );

                    // u - \gamma \Delta t F_R(u) - u^{(s-2+\ell)}
                    return u - gamma * dt * fr_tmp - u_sm2pl;
                };
                auto dg = [&]( state_t& u ) -> matrix_t
                {
                    // Jacobian of g
                    // I - \gamma \Delta t \partial_u F_R(u)
                    return identity - gamma * dt * std::get<reaction_op::value>( pb.system ).df( tn, u );
                };
                u_sp1 = diagonal_implicit_runge_kutta::newton<value_t>( g_sp1,
                    dg,
                    u_sm2pl,
                    ::ponio::linear_algebra::linear_algebra<matrix_t>::solver,
                    ponio::default_config::newton_tolerance,
                    ponio::default_config::newton_max_iterations );

                pb( advection_op(), tn, u_sp1, Fa_u_sp1 );
                pb( diffusion_op(), tn, u_sp1, fd_tmp );
                pb( reaction_op(), tn, u_sp1, fr_tmp_bis );

                auto g_sp2 = [&]( state_t& u ) -> state_t
                {
                    _info.number_of_eval[0] += 1;
                    _info.number_of_eval[1] += 2;

                    pb( reaction_op(), tn, u, fr_tmp );

                    // u - \gamma \Delta t F_R(u) - (u^{(s-2+\ell)} + \beta \Delta t F_D(u^{(s+1)}) + \Delta t F_A(u^{(s+1)}) +
                    // (1-2\gamma)\Delta t F_R(u^{(s+1)}))
                    return u - gamma * dt * fr_tmp - ( u_sm2pl + beta * dt * fd_tmp + dt * Fa_u_sp1 + ( 1. - 2. * gamma ) * dt * fr_tmp_bis );
                };
                u_sp2 = diagonal_implicit_runge_kutta::newton<value_t>( g_sp2,
                    dg,
                    u_sm2pl,
                    ::ponio::linear_algebra::linear_algebra<matrix_t>::solver,
                    ponio::default_config::newton_tolerance,
                    ponio::default_config::newton_max_iterations );
            }

            _info.number_of_eval[1] += 3;

            // u^{(s+3)} = u^{(s-2+\ell)} + (1-2\gamma)\Delta t F_A(u^{(s+1)}) + (1-\gamma)\Delta t F_R(u^{(s+1)})
            pb( reaction_op(), tn, u_sp1, fr_tmp );
            auto& u_sp3 = U[17];
            u_sp3       = u_sm2pl + ( 1. - 2. * gamma ) * dt * Fa_u_sp1 + ( 1. - gamma ) * dt * fr_tmp;

            // u^{(s+4)} = u^{(s+2-\ell)} + 1/3 \Delta t F_A(u^{(s+1)})
            auto& u_sp4 = U[18];
            u_sp4       = u_sm2pl + 1. / 3. * dt * Fa_u_sp1;

            // u^{(s+5)} = u^{(s+2-\ell)} + 2/3 \beta \Delta t F_D(u^{(s+1)}) + 2/3 \Delta t J_R^{-1} F_A(u^{(s+4)}) + (2/3-\gamma) \Delta t
            // F_R(u^{(s+1)}) + 2/3 \gamma \Delta t F_R(u^{(s+2)}) with Shampine's trick J_R = I - \gamma \Delta t \partial_u
            // F_R(u^{(s-2+\ell)}) else J_R = I

            auto& u_sp5 = U[19];

            if constexpr ( shampine_trick_enable
                           && detail::problem_operator<std::tuple_element_t<reaction_op::value, decltype( pb.system )>, value_t> )
            {
                auto& shampine_element = U[21]; // solution of (I - \gamma \Delta t \partial_u F_R ) X = F_A(u^{(s+4)})
                auto& f_A_u            = U[22];
                auto& u_tmp            = U[23];

                pb( advection_op(), tn, u_sp4, f_A_u );
                f_A_u = dt * f_A_u;

                shampine_trick_caller.template operator()<1>( gamma * dt,
                    std::get<reaction_op::value>( pb.system ).f_t( tn ),
                    u_sm2pl,
                    f_A_u,
                    u_tmp,
                    shampine_element );

                pb( diffusion_op(), tn, u_sp1, fd_tmp );
                pb( reaction_op(), tn, u_sp1, fr_tmp );
                pb( reaction_op(), u_sp2, fr_tmp_bis );
                u_sp5 = u_sm2pl + 2. / 3. * beta * dt * fd_tmp + 2. / 3. * shampine_element + ( 2. / 3. - gamma ) * dt * fr_tmp
                      + 2. / 3. * gamma * dt * fr_tmp_bis;
            }
            else
            {
                pb( diffusion_op(), tn, u_sp1, fd_tmp );
                pb( advection_op(), tn, u_sp4, fa_tmp );
                pb( reaction_op(), tn, u_sp1, fr_tmp );
                pb( reaction_op(), tn, u_sp2, fr_tmp_bis );
                u_sp5 = u_sm2pl + 2. / 3. * beta * dt * fd_tmp + 2. / 3. * dt * fa_tmp + ( 2. / 3. - gamma ) * dt * fr_tmp
                      + 2. / 3. * gamma * dt * fr_tmp_bis;
            }

            value_t tau   = sigma * rock_coeff::fp2[deg_index - 1] + sigma * sigma;
            value_t tau_a = 0.5 * detail::power<2>( alpha - 1. ) + 2. * alpha * ( 1. - alpha ) * sigma + alpha * alpha * tau;

            if constexpr ( shampine_trick_enable
                           && detail::problem_operator<std::tuple_element_t<reaction_op::value, decltype( pb.system )>, value_t> )
            {
                auto& shampine_element = U[20];
                auto& f_D_u            = U[21];
                auto& u_tmp            = U[22];

                // for embedded method
                auto& err_D = U[23];
                pb( diffusion_op(), tn, us_sm1, fd_tmp );
                pb( diffusion_op(), tn, u_sm2, fd_tmp_bis );
                // $err_D = \sigma_\alpha(1-\tau_a/\sigma_a^2)\Delta t (F_D(u^{*(s-1)}) - F_D(u^{(s-2)}))$
                err_D = sigma_a * ( 1. - tau_a / ( sigma_a * sigma_a ) ) * dt * ( fd_tmp - fd_tmp_bis );

                pb( diffusion_op(), tn, u_sp3, fd_tmp );
                pb( diffusion_op(), tn, u_sp1, fd_tmp_bis );
                f_D_u = static_cast<state_t>( fd_tmp - fd_tmp_bis );

                shampine_trick_caller.template operator()<l>( gamma * dt,
                    std::get<reaction_op::value>( pb.system ).f_t( tn ),
                    u_sm2pl,
                    f_D_u,
                    u_tmp,
                    shampine_element );

                if constexpr ( is_embedded )
                {
                    auto& rhs_R = U[24];
                    auto& err_R = U[25];
                    auto& err_A = U[26];

                    _info.number_of_eval[1] += 2;

                    pb( reaction_op(), tn, u_sp1, fr_tmp );
                    pb( reaction_op(), tn, u_sp2, fr_tmp_bis );
                    rhs_R = dt * ( fr_tmp - fr_tmp_bis ) / 6.;

                    // $err_R = J_R^{-1} \Delta t/6 (F_R(u^{s+1}) - F_R(u^{s+2}))$
                    // to compute it, get $rhs_R = \Delta t/6 (F_R(u^{s+1}) - F_R(u^{s+2}))$
                    // then solve $J_R err_R = rhs_R$ (that what Shampine's trick does, it build $J_R$ and solve it)
                    shampine_trick_caller.template operator()<1>( gamma * dt,
                        std::get<reaction_op::value>( pb.system ).f_t( tn ),
                        u_sm2pl,
                        rhs_R,
                        u_tmp,
                        err_R );

                    pb( advection_op(), tn, u_sp4, fa_tmp );
                    pb( advection_op(), tn, u_sp5, fa_tmp_bis );
                    // err_A = -3/20 \Delta t F_A(u^{(s+1)}) + 3/10 \Delta t F_A(u^{(s+4)}) - 3/20 \Delta t F_A(u^{(s+5)})
                    err_A = -0.15 * dt * Fa_u_sp1 + 0.3 * dt * fa_tmp - 0.15 * dt * fa_tmp_bis;

                    u_np1 = us_s - err_D + 0.5 * dt * fr_tmp + 0.25 * dt * Fa_u_sp1 + 0.75 * dt * fa_tmp_bis + 0.5 * dt * fr_tmp
                          + 0.5 * dt * fr_tmp_bis + 1.0 / ( 2. - 4. * gamma ) * shampine_element;

                    auto accumulator_error_gen = []( auto const& yn, auto const& ynp1, value_t a_tol, value_t r_tol )
                    {
                        return [=, it_yn = yn.begin(), it_ynp1 = ynp1.begin()]( value_t const& acc, value_t const err_i ) mutable
                        {
                            return acc
                                 + detail::power<2>( err_i / ( a_tol + r_tol * std::max( std::abs( *it_yn++ ), std::abs( *it_ynp1++ ) ) ) );
                        };
                    };

                    // TODO: this couple of lines works only with samurai (because of err_D.array())
                    value_t err_R_scalar = norm_error( err_R, un, u_np1, _info.absolute_tolerance, _info.relative_tolerance );
                    value_t err_D_scalar = norm_error( err_D, un, u_np1, _info.absolute_tolerance, _info.relative_tolerance );
                    value_t err_A_scalar = norm_error( err_A, un, u_np1, _info.absolute_tolerance, _info.relative_tolerance );

                    _info.error   = std::max( err_D_scalar, err_R_scalar, std::pow( err_A_scalar, 2. / 3. ) );
                    _info.success = _info.error < 1.0;

                    value_t fac    = std::min( 2.0, std::max( 0.5, std::sqrt( 1.0 / _info.error ) ) );
                    value_t new_dt = 0.8 * fac * dt;

                    // accepted step
                    if ( _info.success )
                    {
                        tn = tn + dt;
                        dt = new_dt;

                        // return { tn + dt, u_np1, new_dt };
                    }
                    else
                    {
                        // tn = tn;
                        std::swap( un, u_np1 );
                        dt = new_dt;

                        // return { tn, un, new_dt };
                    }
                }

                pb( reaction_op(), tn, u_sp1, fr_tmp );
                pb( advection_op(), tn, u_sp5, fa_tmp );
                pb( reaction_op(), tn, u_sp2, fr_tmp_bis );

                tn    = tn + dt;
                u_np1 = us_s - err_D + 0.5 * dt * fr_tmp + 0.25 * dt * Fa_u_sp1 + 0.75 * dt * fa_tmp + 0.5 * dt * fr_tmp
                      + 0.5 * dt * fr_tmp_bis + dt / ( 2. - 4. * gamma ) * shampine_element;
            }
            else
            {
                auto& fd_tmp_ter = U[8];
                auto& fd_tmp_qua = U[9];

                pb( diffusion_op(), tn, us_sm1, fd_tmp );
                pb( diffusion_op(), tn, u_sm2, fd_tmp_bis );
                pb( advection_op(), tn, u_sp5, fa_tmp );
                pb( reaction_op(), tn, u_sp1, fr_tmp );
                pb( reaction_op(), tn, u_sp2, fr_tmp_bis );
                pb( diffusion_op(), tn, u_sp3, fd_tmp_ter );
                pb( diffusion_op(), tn, u_sp1, fd_tmp_qua );

                tn    = tn + dt;
                u_np1 = us_s - sigma_a * ( 1. - tau_a / ( sigma_a * sigma_a ) ) * dt * ( fd_tmp - fd_tmp_bis ) + 0.25 * dt * Fa_u_sp1
                      + 0.75 * dt * fa_tmp + 0.5 * dt * fr_tmp + 0.5 * dt * fr_tmp_bis
                      + dt / ( 2. - 4. * gamma ) * ( fd_tmp_ter - fd_tmp_qua );
            }

            // return { tn + dt, u_np1, dt };
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
     * @param absolute_tolerance       Absolute tolerance for adaptive time step method (default: ponio::default_config::tol)
     * @param relatove_tolerance       Relative tolerance for adaptive time step method (default: ponio::default_config::tol)
     */
    template <std::size_t l = 1, bool is_embedded = false, typename value_t = double, typename alpha_beta_computer_t, typename eig_computer_t, typename shampine_trick_caller_t>
    auto
    pirock_RDA( alpha_beta_computer_t&& alpha_beta_computer,
        eig_computer_t&& eig_computer,
        shampine_trick_caller_t&& shampine_trick_caller,
        value_t absolute_tolerance = default_config::tol,
        value_t relative_tolerance = default_config::tol )
    {
        return pirock_RDA_impl<l, alpha_beta_computer_t, eig_computer_t, shampine_trick_caller_t, is_embedded, value_t>(
            std::forward<alpha_beta_computer_t>( alpha_beta_computer ),
            std::forward<eig_computer_t>( eig_computer ),
            std::forward<shampine_trick_caller_t>( shampine_trick_caller ),
            absolute_tolerance,
            relative_tolerance );
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
