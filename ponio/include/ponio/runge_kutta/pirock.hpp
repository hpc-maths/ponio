// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once
#include <array>
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
#include "dirk.hpp"
#include "rock.hpp"

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
     * @brief Computer of \f$\alpha\f$ and \f$\beta\f$ parameters with fixed value of \f$\alpha\f$
     *
     * @tparam value_t
     */
    template <typename value_t = double>
    struct alpha_fixed
    {
        value_t _alpha;

        alpha_fixed( value_t a = static_cast<value_t>( 1. ) )
            : _alpha( a )
        {
        }

        inline value_t
        alpha( std::size_t, std::size_t ) const
        {
            return _alpha;
        }

        inline value_t
        beta( std::size_t s, std::size_t l ) const
        {
            return 1. - 2. * alpha( s, l ) * polynomial::Pp_sm2pl_0<value_t>( s, l );
        }
    };

    /**
     * @brief Computer of \f$\alpha\f$ and \f$\beta\f$ parameters with fixed value of \f$\beta\f$ to 0
     *
     * @tparam value_t
     */
    template <typename value_t = double>
    struct beta_0
    {
        inline value_t
        alpha( std::size_t s, std::size_t l ) const
        {
            return 1. / ( 2. * polynomial::Pp_sm2pl_0<value_t>( s, l ) );
        }

        inline value_t
        beta( std::size_t, std::size_t ) const
        {
            return 0.;
        }
    };

    /**
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
        typename value_t                  = double>
    struct pirock_impl
    {
        static constexpr bool is_embedded           = _is_embedded;
        static constexpr bool shampine_trick_enable = !std::is_void<_shampine_trick_caller_t>::value;

        static constexpr std::size_t l         = _l;
        static constexpr std::size_t N_stages  = stages::dynamic;
        static constexpr std::size_t N_storage = std::
            conditional_t<shampine_trick_enable, std::integral_constant<std::size_t, 16>, std::integral_constant<std::size_t, 10>>::value;
        static constexpr std::size_t order   = 2;
        static constexpr std::string_view id = "PIROCK";

        using rock_coeff              = rock::rock2_coeff<value_t>;
        using degree_computer         = rock::detail::degree_computer<value_t, rock_coeff>;
        using shampine_trick_caller_t = typename std::conditional<shampine_trick_enable, _shampine_trick_caller_t, bool>::type;

        alpha_beta_computer_t alpha_beta_computer;
        eig_computer_t eig_computer;
        shampine_trick_caller_t shampine_trick_caller;
        value_t tolerance;

        /**
         * @brief Construct a new pirock impl object
         *
         */
        pirock_impl()
        {
        }

        /**
         * @brief Construct a new pirock impl object
         *
         * @param _alpha_beta_computer computer of parameters alpha and beta object that has two member functions which take number of
         * stages (s) l parameter
         * @param _eig_computer        eigenvalue computer functor (that take as argument the function of explicit part, the current time,
         * the current state and the current time step)
         */
        pirock_impl( alpha_beta_computer_t&& _alpha_beta_computer, eig_computer_t&& _eig_computer )
            : alpha_beta_computer( _alpha_beta_computer )
            , eig_computer( _eig_computer )
        {
        }

        /**
         * @brief Construct a new pirock impl object with Shampine's trick
         *
         * @param _alpha_beta_computer   alpha and beta computer object
         * @param _eig_computer          eigenvalue computer functor
         * @param _shampine_trick_caller Shampine's trick functor
         * @param tol                    tolerance
         */
        template <typename _shampine_trick_caller_t_>
            requires std::same_as<_shampine_trick_caller_t_, shampine_trick_caller_t>
                      && std::same_as<std::bool_constant<shampine_trick_enable>, std::true_type>
        pirock_impl( alpha_beta_computer_t&& _alpha_beta_computer,
            eig_computer_t&& _eig_computer,
            _shampine_trick_caller_t_&& _shampine_trick_caller,
            value_t tol = default_config::tol )
            : alpha_beta_computer( _alpha_beta_computer )
            , eig_computer( _eig_computer )
            , shampine_trick_caller( _shampine_trick_caller )
            , tolerance( tol )
        {
        }

        template <typename problem_t, typename state_t, typename array_ki_t>
        inline std::tuple<value_t, state_t, value_t>
        operator()( problem_t& pb, value_t& tn, state_t& un, array_ki_t& U, value_t& dt )
        {
            static_assert( detail::problem_operator<decltype( pb.implicit_part ), value_t>
                               || detail::problem_jacobian<decltype( pb.implicit_part ), value_t, state_t>,
                "This kind of problem is not inversible in ponio" );

            auto [mdeg, deg_index, start_index] = degree_computer::compute_n_stages_optimal_degree( eig_computer, pb, tn, un, dt, 4 );
            std::size_t s                       = mdeg + 2;

            value_t const alpha = alpha_beta_computer.alpha( s, l );
            value_t const beta  = alpha_beta_computer.beta( s, l );
            value_t const gamma = 1. - 0.5 * std::sqrt( 2. );

            auto& u_j   = U[0];
            auto& u_jm1 = U[1];
            auto& u_jm2 = U[2];
            auto& u_sm2 = U[3];

            u_j   = un;
            u_jm2 = un;

            value_t const mu_1 = rock_coeff::recf[start_index - 1];

            value_t t_jm1 = tn + dt * alpha * mu_1;
            value_t t_jm2 = tn + dt * alpha * mu_1;
            value_t t_jm3 = tn;

            // u_1 =u^n + \alpha \mu_1 \Delta F_D( u^n )
            u_jm1 = un + alpha * dt * mu_1 * pb.explicit_part( tn, un );

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
                u_j = alpha * mu_j * dt * pb.explicit_part( t_jm1, u_jm1 ) - nu_j * u_jm1 - kappa_j * u_jm2;

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
            auto& us_sm1 = U[4];
            us_sm1       = u_sm2 + sigma_a * dt * pb.explicit_part( t_jm1, u_sm2 );

            // u_{*s} = u_{*s-1} + \sigma_\alpha \Delta t  F_D( u_{*s-1} )
            auto& us_s = U[5];
            us_s       = us_sm1 + sigma_a * dt * pb.explicit_part( t_jm1, us_sm1 );

            // u_{s-2+l} = u_j
            auto& u_sm2pl = u_j;

            auto& u_sp1 = U[6];
            u_sp1       = un;

            auto& u_sp2 = U[7];
            u_sp2       = un;

            if constexpr ( detail::problem_operator<decltype( pb.implicit_part ), value_t> )
            {
                auto op_sp1  = ::ponio::linear_algebra::operator_algebra<state_t>::identity( un ) - gamma * dt * pb.implicit_part.f_t( tn );
                auto rhs_sp1 = u_sm2pl;
                ::ponio::linear_algebra::operator_algebra<state_t>::solve( op_sp1, u_sp1, rhs_sp1 );

                auto op_sp2  = ::ponio::linear_algebra::operator_algebra<state_t>::identity( un ) - gamma * dt * pb.implicit_part.f_t( tn );
                auto rhs_sp2 = static_cast<state_t>(
                    u_sm2pl + beta * dt * pb.explicit_part( tn, u_sp1 ) + ( 1. - 2. * gamma ) * dt * pb.implicit_part( tn, u_sp1 ) );
                ::ponio::linear_algebra::operator_algebra<state_t>::solve( op_sp2, u_sp2, rhs_sp2 );
            }
            else
            {
                using matrix_t = decltype( pb.implicit_part.df( tn, un ) );

                auto identity = ::ponio::linear_algebra::linear_algebra<matrix_t>::identity( un );
                auto g_sp1    = [&]( state_t const& u ) -> state_t
                {
                    return u - gamma * dt * pb.implicit_part.f( tn, u ) - u_sm2pl;
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

                auto g_sp2 = [&]( state_t const& u ) -> state_t
                {
                    return u - gamma * dt * pb.implicit_part.f( tn, u )
                         - ( u_sm2pl + beta * dt * pb.explicit_part( tn, u_sp1 ) + ( 1. - 2. * gamma ) * dt * pb.implicit_part( tn, u_sp1 ) );
                };
                u_sp2 = diagonal_implicit_runge_kutta::newton<value_t>( g_sp2,
                    dg,
                    u_sm2pl,
                    ::ponio::linear_algebra::linear_algebra<matrix_t>::solver,
                    ponio::default_config::newton_tolerance,
                    ponio::default_config::newton_max_iterations );
            }

            auto& u_sp3 = U[8];
            u_sp3       = u_sm2pl + ( 1. - gamma ) * dt * pb.implicit_part( tn, u_sp1 );

            value_t tau   = sigma * rock_coeff::fp2[deg_index - 1] + sigma * sigma;
            value_t tau_a = 0.5 * detail::power<2>( alpha - 1. ) + 2. * alpha * ( 1. - alpha ) * sigma + alpha * alpha * tau;

            auto& u_np1 = U[9];

            if constexpr ( shampine_trick_enable && detail::problem_operator<decltype( pb.implicit_part ), value_t> )
            {
                auto& shampine_element = U[10];
                auto& f_D_u            = U[11];
                auto& u_tmp            = U[12];

                // for embedded method
                auto& err_D = U[13];

                // $err_D = \sigma_\alpha(1-\tau_a/\sigma_a^2)\Delta t (F_D(u^{*(s-1)}) - F_D(u^{(s-2)}))$
                err_D = sigma_a * ( 1. - tau_a / ( sigma_a * sigma_a ) ) * dt
                      * ( pb.explicit_part( tn, us_sm1 ) - pb.explicit_part( tn, u_sm2 ) );

                f_D_u = static_cast<state_t>( pb.explicit_part( tn, u_sp3 ) - pb.explicit_part( tn, u_sp1 ) );

                // shampine_trick_caller.template operator()<l>( gamma * dt, pb.implicit_part.f_t( tn ), u_sm2pl, f_D_u, u_tmp,
                // shampine_element );
                shampine_element = f_D_u;

                // if constexpr ( is_embedded )
                // {
                //     auto& rhs_R = U[14];
                //     auto& err_R = U[15];

                //     rhs_R = dt * ( pb.implicit_part( tn, u_sp1 ) - pb.implicit_part( tn, u_sp2 ) ) / 6.;

                //     // $err_R = J_R^{-1} \Delta t/6 (F_R(u^{s+1}) - F_R(u^{s+2}))$
                //     // to compute it, get $rhs_R = \Delta t/6 (F_R(u^{s+1}) - F_R(u^{s+2}))$
                //     // then solve $J_R err_R = rhs_R$ (that what Shampine's trick does, it build $J_R$ and solve it)
                //     shampine_trick_caller.template operator()<1>( gamma * dt, pb.implicit_part.f_t( tn ), u_sm2pl, rhs_R, u_tmp, err_R );

                //     // TODO: this couple of lines works only with samurai (because of err_D.array())
                //     auto err = std::max( std::accumulate( err_D.array().begin(),
                //                              err_D.array().end(),
                //                              static_cast<value_t>( 0. ),
                //                              []( value_t const& acc, value_t const xi )
                //                              {
                //                                  return acc + std::abs( xi );
                //                              } ),
                //         std::accumulate( err_R.array().begin(),
                //             err_R.array().end(),
                //             static_cast<value_t>( 0. ),
                //             []( value_t const& acc, value_t const xi )
                //             {
                //                 return acc + std::abs( xi );
                //             } ) );

                //     value_t new_dt = dt;
                //     if ( err > 0. )
                //     {
                //         new_dt = 0.8 * std::sqrt( tolerance / err ) * dt;
                //     }

                //     if ( err > tolerance )
                //     {
                //         return { tn, un, new_dt };
                //     }

                //     u_np1 = us_s - err_D + 0.5 * dt * pb.implicit_part( tn, u_sp1 ) + 0.5 * dt * pb.implicit_part( tn, u_sp2 )
                //           + 1.0 / ( 2. - 4. * gamma ) * shampine_element;
                //     return { tn + dt, u_np1, new_dt };
                // }

                u_np1 = us_s - err_D + 0.5 * dt * pb.implicit_part( tn, u_sp1 ) + 0.5 * dt * pb.implicit_part( tn, u_sp2 )
                      + dt / ( 2. - 4. * gamma ) * shampine_element;
            }
            else
            {
                u_np1 = us_s
                      - sigma_a * ( 1. - tau_a / ( sigma_a * sigma_a ) ) * dt
                            * ( pb.explicit_part( tn, us_sm1 ) - pb.explicit_part( tn, u_sm2 ) )
                      + 0.5 * dt * pb.implicit_part( tn, u_sp1 ) + 0.5 * dt * pb.implicit_part( tn, u_sp2 )
                      + dt / ( 2. - 4. * gamma ) * ( pb.explicit_part( tn, u_sp3 ) - pb.explicit_part( tn, u_sp1 ) );
            }

            return { tn + dt, u_np1, dt };
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
     * @param tolerance                Tolerance for adaptive time step method (default: ponio::default_config::tol)
     */
    template <std::size_t l = 1, bool is_embedded = false, typename value_t = double, typename alpha_beta_computer_t, typename eig_computer_t, typename shampine_trick_caller_t>
    auto
    pirock( alpha_beta_computer_t&& alpha_beta_computer,
        eig_computer_t&& eig_computer,
        shampine_trick_caller_t&& shampine_trick_caller,
        value_t tolerance = default_config::tol )
    {
        return pirock_impl<l, alpha_beta_computer_t, eig_computer_t, shampine_trick_caller_t, is_embedded, value_t>(
            std::forward<alpha_beta_computer_t>( alpha_beta_computer ),
            std::forward<eig_computer_t>( eig_computer ),
            std::forward<shampine_trick_caller_t>( shampine_trick_caller ),
            tolerance );
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

} // namespace ponio::runge_kutta::pirock
