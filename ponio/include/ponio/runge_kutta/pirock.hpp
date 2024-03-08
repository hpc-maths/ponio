// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once
#include <concepts>
#include <numeric>
#include <ranges>
#include <string_view>
#include <type_traits>
#include <valarray>

#include "../butcher_tableau.hpp"
#include "../detail.hpp"
#include "../linear_algebra.hpp"
#include "../ponio_config.hpp"
#include "../stage.hpp"
#include "dirk.hpp"
#include "rock.hpp"

namespace ponio::runge_kutta::pirock
{
    /*
    template <typename value_t = double>
    struct piheun
    {
        static constexpr bool is_embedded     = false;
        static constexpr std::size_t N_stages = 5;

        value_t delta;
        value_t gamma;
        value_t beta;

        piheun()
            : delta( 1.0 )
            , gamma( 1.0 )
            , beta( 1.0 )
        {
        }

        template <typename problem_t, typename state_t, typename array_ki_t, std::size_t j>
        inline state_t
        stage( Stage<0>, problem_t& pb, value_t tn, state_t& un, array_ki_t const& Ki, value_t dt )
        {
            return un + dt * pb.F_d( un );
        }

        template <typename problem_t, typename state_t, typename array_ki_t, std::size_t j>
        inline state_t
        stage( Stage<1>, problem_t& pb, value_t tn, state_t& un, array_ki_t const& Ki, value_t dt )
        {
            return un + dt * delta * pb.F_d( un );
        }

        template <typename problem_t, typename state_t, typename array_ki_t, std::size_t j>
        inline state_t
        stage( Stage<2>, problem_t& pb, value_t tn, state_t& un, array_ki_t const& Ki, value_t dt )
        {
            // TODO: Newton method
            auto g = [&]( k3 )
            {
                return k3 - ( Ki[1] + dt * gamma * pb.F_r( k3 ) );
            };

            auto dg = [&]( k3 )
            {
                return Identity - dt * gamma * pb.JF_r( k3 );
            };

            // ....

            return Ki[1] + dt * gamma * pb.F_r( Ki[2] );
        }

        template <typename problem_t, typename state_t, typename array_ki_t, std::size_t j>
        inline state_t
        stage( Stage<3>, problem_t& pb, value_t tn, state_t& un, array_ki_t const& Ki, value_t dt )
        {
            // TODO: Newton method
            auto g = [&]( k4 )
            {
                return k4 - ( Ki[1] + dt * beta * pb.F_d( Ki[2] ) + ( 1 - 2 * gamma ) * dt * pb.F_r( Ki[2] ) + gamma * dt * pb.F_r( k4 ) );
            };

            auto Jg = [&]( k4 )
            {
                return Identity - gamma * dt * pb.JF_r( k4 );
            };

            // ....

            return Ki[1] + dt * beta * pb.F_d( Ki[2] ) + ( 1 - 2 * gamma ) * dt * pb.F_r( Ki[2] ) + gamma * dt * pb.F_r( Ki[3] );
        }

        template <typename problem_t, typename state_t, typename array_ki_t, std::size_t j>
        inline state_t
        stage( Stage<4>, problem_t& pb, value_t tn, state_t& un, array_ki_t const& Ki, value_t dt )
        {
            return Ki[1] + dt * ( 1.0 - gamma ) * pb.F_r( Ki[2] );
        }

        template <typename problem_t, typename state_t, typename array_ki_t, std::size_t j>
        inline state_t
        stage( Stage<5>, problem_t& pb, value_t tn, state_t& un, array_ki_t const& Ki, value_t dt )
        {
            return un + 0.5 * dt * pb.F_d( Ki[0] ) + 0.5 * dt * pb.F_d( Ki[1] )
                 + dt / ( 2.0 - 4.0 * gamma ) * ( pb.F_d( Ki[4] ) - pb.F_d( Ki[2] ) ) + 0.5 * dt * pb.F_r( Ki[2] )
                 + 0.5 * dt * pb.F_r( Ki[3] );
        }
    };
    */

    namespace polynomial
    {
        template <std::size_t N, typename value_t = double>
        struct polynomial
        {
            static constexpr std::size_t degree = N;

            std::array<value_t, degree + 1> coeff;

            constexpr double
            operator[]( std::size_t idx ) const
            {
                if ( idx < coeff.size() )
                {
                    return coeff[idx];
                }
                else
                {
                    return static_cast<value_t>( 0. );
                }
            }

            template <std::size_t M>
            polynomial<N, value_t>&
            operator=( polynomial<M, value_t> const& rhs )
            {
                for ( std::size_t i = 0; i < coeff.size(); ++i )
                {
                    coeff[i] = rhs[i];
                }

                return *this;
            }

            template <std::size_t... Is>
            constexpr value_t
            eval_impl( value_t&& x, std::index_sequence<Is...> ) const
            {
                return ( 0. + ... + ( coeff[Is] * detail::power<Is>( std::forward<value_t>( x ) ) ) );
            }

            constexpr value_t
            operator()( value_t&& x ) const
            {
                return eval_impl( std::forward<value_t>( x ), std::make_index_sequence<degree + 1>{} );
            }
        };

        template <std::size_t N, std::size_t M, typename value_t>
        constexpr polynomial<std::max( N, M ), value_t> const&
        max_degree( polynomial<N, value_t> const& P, polynomial<M, value_t> const& Q )
        {
            if constexpr ( N >= M )
            {
                return P;
            }
            else
            {
                return Q;
            }
        }

        template <std::size_t N, std::size_t M, typename value_t>
        constexpr polynomial<std::min( N, M ), value_t> const&
        min_degree( polynomial<N, value_t> const& P, polynomial<M, value_t> const& Q )
        {
            if constexpr ( N >= M )
            {
                return Q;
            }
            else
            {
                return P;
            }
        }

        /* operator + (polynomial, polynomial) */
        template <std::size_t N, std::size_t M, typename value_t, std::size_t... Is, std::size_t... Js>
        constexpr polynomial<N, value_t>
        add_poly_poly( polynomial<N, value_t> const& P, polynomial<M, value_t> const& Q, std::index_sequence<Is...>, std::index_sequence<Js...> )
        {
            return {
                {( P[Is] + Q[Is] )..., P[Js + M]...}
            };
        }

        template <std::size_t N, std::size_t M, typename value_t>
        polynomial<std::max( N, M ), value_t>
        operator+( polynomial<N, value_t> const& P, polynomial<M, value_t> const& Q )
        {
            return add_poly_poly( max_degree( P, Q ),
                min_degree( P, Q ),
                std::make_index_sequence<std::min( N, M ) + 1>{},
                std::make_index_sequence<std::max( N, M ) - std::min( N, M )>{} );
        }

        /* operator - (polynomial, polynomial) */
        template <std::size_t N, std::size_t M, typename value_t>
        polynomial<std::max( N, M ), value_t>
        operator-( polynomial<N, value_t> const& P, polynomial<M, value_t> const& Q )
        {
            return P + ( -1. ) * Q;
        }

        /* operator * (value_t, polynomial) */
        template <std::size_t N, typename value_t, std::size_t... Is>
        constexpr polynomial<N, value_t>
        mul_scalar_poly( value_t a, polynomial<N, value_t> const& P, std::index_sequence<Is...> )
        {
            return polynomial<N, value_t>( { { ( a * P[Is] )... } } );
        }

        template <std::size_t N, typename value_t>
        polynomial<N, value_t>
        operator*( value_t a, polynomial<N, value_t> const& P )
        {
            return mul_scalar_poly( a, P, std::make_index_sequence<N + 1>{} );
        }

        /* operator * (polynomial, polynomial) */
        template <std::size_t K, std::size_t N, std::size_t M, typename value_t, std::size_t... Is>
        constexpr value_t
        mul_poly_poly_coeff( polynomial<N, value_t> const& P, polynomial<M, value_t> const& Q, std::index_sequence<Is...> )
        {
            return ( static_cast<value_t>( 0. ) + ... + ( P[Is] * Q[K - Is] ) );
        }

        template <std::size_t N, std::size_t M, typename value_t, std::size_t... Is>
        constexpr polynomial<N + M, value_t>
        mul_poly_poly( polynomial<N, value_t> const& P, polynomial<M, value_t> const& Q, std::index_sequence<Is...> )
        {
            return { { mul_poly_poly_coeff<Is>( P, Q, std::make_index_sequence<N + M>{} )... } };
        }

        template <std::size_t N, std::size_t M, typename value_t>
        constexpr polynomial<N + M, value_t>
        operator*( polynomial<N, value_t> const& P, polynomial<M, value_t> const& Q )
        {
            return mul_poly_poly( max_degree( P, Q ), min_degree( P, Q ), std::make_index_sequence<N + M + 1>{} );
        }

        /* derivation */
        template <std::size_t N, typename value_t, std::size_t... Is>
        constexpr polynomial<std::max( 0ul, N - 1 ), value_t>
        deriv_impl( polynomial<N, value_t> const& P, std::index_sequence<Is...> )
        {
            return { { ( ( Is + 1 ) * P[Is + 1] )... } };
        }

        template <std::size_t N, typename value_t>
        constexpr polynomial<std::max( 0ul, N - 1 ), value_t>
        deriv( polynomial<N, value_t> const& P )
        {
            return deriv_impl( P, std::make_index_sequence<N>{} );
        }

        /* display */
        template <std::size_t N, typename value_t>
        std::ostream&
        operator<<( std::ostream& os, polynomial<N, value_t> const& P )
        {
            os << P[0];
            if constexpr ( N > 0 )
            {
                for ( std::size_t i = 0; i < P.coeff.size(); ++i )
                {
                    os << " + " << P[i] << "*X";
                    if ( i > 1 )
                    {
                        os << "**" << i;
                    }
                }
            }

            return os;
        }

        template <std::size_t s, typename rock_coeff>
        std::pair<std::size_t, std::size_t>
        optimal_degree()
        {
            std::size_t mdeg = s - 2;
            std::size_t mz = 1, mr = 1;

            if ( mdeg < 2 )
            {
                return { mz, mr };
            }

            std::size_t i = 1;
            for ( auto ms_i : rock_coeff::ms )
            {
                if ( ms_i / mdeg >= 1 )
                {
                    // mdeg = rock_coeff::ms[i - 1];
                    mz = i;
                    break;
                }
                mr = mr + rock_coeff::ms[i - 1] * 2 - 1;

                ++i;
            }

            return { mz, mr };
        }

        template <std::size_t s, std::size_t l, typename value_t>
        polynomial<s - 2 + l>
        P_sm2pl()
        {
            using rock_coeff = rock::rock2_coeff<value_t>;

            polynomial<1, value_t> z   = { 0., 1. };
            polynomial<0, value_t> one = { 1. };

            auto f = [=]( auto u )
            {
                return z * u;
            };

            std::size_t mdeg = s - 2;
            auto [mz, mr]    = rock::detail::degree_computer<value_t, rock_coeff>::optimal_degree( mdeg );

            polynomial<s - 2 + l, value_t> u;
            polynomial<s - 2 + l, value_t> ujm1;
            polynomial<s - 2 + l, value_t> ujm2;

            u    = one;
            ujm2 = one;

            auto mu_1 = rock_coeff::recf[mr - 1];
            ujm1      = u + mu_1 * f( u );

            if constexpr ( s <= 2 )
            {
                u = ujm1;
            }

            for ( std::size_t j = 2; j < s - 2 + l + 1; ++j )
            {
                double const mu    = rock_coeff::recf[mr + 2 * ( j - 2 ) + 1 - 1];
                double const kappa = rock_coeff::recf[mr + 2 * ( j - 2 ) + 2 - 1];
                double const nu    = -1.0 - kappa;

                u = mu * f( ujm1 ) - nu * ujm1 - kappa * ujm2;

                ujm2 = ujm1;
                ujm1 = u;
            }

            return u;
        }

        template <std::size_t s, std::size_t l, typename value_t = double>
        value_t
        Pp_sm2pl_0()
        {
            return deriv( P_sm2pl<s, l, value_t>() )( 0. );
        }

        template <typename value_t = double>
        value_t
        Pp_sm2pl_0( std::size_t s, std::size_t l )
        {
            using rock_coeff        = rock::rock2_coeff<value_t>;
            std::size_t const sm2pl = s - 2 + l;

            std::valarray<value_t> u( 0., sm2pl + 1 );
            std::valarray<value_t> ujm1( 0., sm2pl + 1 );
            std::valarray<value_t> ujm2( 0., sm2pl + 1 );

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

                u = -nu * ujm1 - kappa * ujm2; // u_j = -nu * u_jm1 - kappa * u_jm2

                std::valarray<value_t> slice = ujm1[std::slice( 0, sm2pl, 1 )];
                u[std::slice( 1, sm2pl + 1, 1 )] += mu * slice; // u_j += mu * f( u_jm1 )

                ujm2 = ujm1;
                ujm1 = u;
            }

            return u[1];
        }

    } // namespace polynomial

    /** @class pirock_impl
     *  @brief define PIROCK method
     *
     *  @warning the implementation is only with l=1 and beta=0 with a fixed number of stages
     */
    template <std::size_t _s, typename value_t = double>
    struct pirock
    {
        static constexpr bool is_embedded      = false;
        static constexpr std::size_t N_stages  = stages::dynamic;
        static constexpr std::size_t N_storage = 10;
        static constexpr std::size_t order     = 2;
        static constexpr std::string_view id   = "PIROCK";

        static constexpr std::size_t s = _s;
        static constexpr std::size_t l = 1;

        using rock_coeff      = rock::rock2_coeff<value_t>;
        using degree_computer = rock::detail::degree_computer<value_t, rock_coeff>;

        // eig_computer_t eig_computer_t;
        value_t alpha;

        pirock()
            : alpha( 1. / ( 2. * polynomial::Pp_sm2pl_0<s, l, value_t>() ) )
        {
        }

        pirock( value_t alpha_ )
            : alpha( alpha_ )
        {
        }

        template <typename problem_t, typename state_t, typename array_ki_t>
        inline std::tuple<value_t, state_t, value_t>
        operator()( problem_t& pb, value_t& tn, state_t& un, array_ki_t& U, value_t& dt )
        {
            // TODO: change here s and mdeg computation to compute this from pb.explicit_part
            std::size_t mdeg              = s - 2;
            auto [deg_index, start_index] = polynomial::optimal_degree<s, rock_coeff>();

            // TODO: change this polynomial into a runtime computation
            // value_t const beta  = 1. - 2. * alpha * polynomial::Pp_sm2pl_0<s, l, value_t>();
            value_t const beta  = 1. - 2. * alpha * polynomial::Pp_sm2pl_0<value_t>( s, l );
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

                u_j = alpha * dt * mu_j * pb.explicit_part( t_jm1, u_jm1 ) - nu_j * u_jm1 - kappa_j * u_jm2;

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
            // u_j -> u_{s-2+l} = u_{s-1}
            // u_jm1 -> u_{s-2+l-1} = u_{s-2}
            // u_jm2 -> u_{s-2+l-2} = u_{s-3}

            value_t sigma   = rock_coeff::fp1[deg_index - 1];
            value_t sigma_a = 0.5 * ( 1.0 - alpha ) + alpha * sigma;
            auto& us_sm1    = U[4];
            us_sm1          = u_sm2 + sigma_a * dt * pb.explicit_part( t_jm1, u_sm2 );

            auto& us_s = U[5];
            us_s       = us_sm1 + sigma_a * dt * pb.explicit_part( t_jm1, us_sm1 );

            auto& u_sm2pl = u_j;

            auto& u_sp1 = U[6];
            u_sp1       = u_sm2pl;

            auto& u_sp2 = U[7];
            u_sp2       = u_sm2pl;

            if constexpr ( detail::problem_operator<decltype( pb.implicit_part ), value_t> )
            {
                auto op_sp1  = ::ponio::linear_algebra::operator_algebra<state_t>::identity( un ) - gamma * dt * pb.implicit_part.f_t( tn );
                auto rhs_sp1 = u_sm2pl;
                ::ponio::linear_algebra::operator_algebra<state_t>::solve( op_sp1, u_sp1, rhs_sp1 );

                auto op_sp2  = ::ponio::linear_algebra::operator_algebra<state_t>::identity( un ) - gamma * dt * pb.implicit_part.f_t( tn );
                auto rhs_sp2 = static_cast<state_t>(
                    u_sm2pl + beta * pb.explicit_part( tn, u_sp1 ) + ( 1. - 2. * gamma ) * dt * pb.implicit_part( tn, u_sp1 ) );
                ::ponio::linear_algebra::operator_algebra<state_t>::solve( op_sp2, u_sp2, rhs_sp2 );
            }
            else if constexpr ( detail::problem_jacobian<decltype( pb.implicit_part ), value_t, state_t> )
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
                         - ( u_sm2pl + beta * pb.explicit_part( tn, u_sp1 ) + ( 1. - 2. * gamma ) * dt * pb.implicit_part( tn, u_sp1 ) );
                };
                u_sp2 = diagonal_implicit_runge_kutta::newton<value_t>( g_sp2,
                    dg,
                    u_sm2pl,
                    ::ponio::linear_algebra::linear_algebra<matrix_t>::solver,
                    ponio::default_config::newton_tolerance,
                    ponio::default_config::newton_max_iterations );
            }
            else
            {
                static_assert( detail::problem_operator<decltype( pb.implicit_part ), value_t>
                                   || detail::problem_jacobian<decltype( pb.implicit_part ), value_t, state_t>,
                    "This kind of problem is not inversible in ponio" );
            }

            auto& u_sp3 = U[8];
            u_sp3       = u_sm2pl + ( 1. - gamma ) * dt * pb.implicit_part( tn, u_sp1 );

            value_t tau   = sigma * rock_coeff::fp2[deg_index - 1] + sigma * sigma;
            value_t tau_a = 0.5 * detail::power<2>( alpha - 1. ) + 2. * alpha * ( 1. - alpha ) * sigma + alpha * alpha * tau;

            auto& u_np1 = U[9];
            u_np1       = us_s
                  - sigma_a * ( 1. - tau_a / ( sigma_a * sigma_a ) ) * dt * ( pb.explicit_part( tn, us_sm1 ) - pb.explicit_part( tn, u_sm2 ) )
                  + 0.5 * dt * pb.implicit_part( tn, u_sp1 ) + 0.5 * dt * pb.implicit_part( tn, u_sp2 )
                  + 1. / ( 2. - 4. * gamma ) * dt * ( pb.explicit_part( tn, u_sp3 ) - pb.explicit_part( tn, u_sp1 ) );

            return { tn + dt, u_np1, dt };
        }
    };

} // namespace ponio::runge_kutta::pirock
