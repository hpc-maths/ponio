// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <numeric>
#include <ranges>
#include <string_view>
#include <tuple>
#include <utility>

#include "../detail.hpp"
#include "../iteration_info.hpp"
#include "../ponio_config.hpp"
#include "../stage.hpp"

#include "rock_coeff.hpp"

namespace ponio::runge_kutta::rock
{
    /**
     * @brief selectors with types
     *
     */
    namespace rock_order
    {
        struct rock_2
        {
        };

        struct rock_4
        {
        };
    } // namespace rock_order

    namespace detail
    {
        template <typename state_t>
        auto
        norm_2( state_t const& u )
        {
            using namespace std;
            return abs( u );
        }

        /**
         * @brief compute norm 2 of a range
         *
         * @tparam state_t type of range
         * @param u        range where compute \f$||u||_2 = \left(\sum_j u_j^2\right)^{\frac{1}{2}}\f$
         */
        template <typename state_t>
            requires std::ranges::range<state_t>
        auto
        norm_2( state_t const& u )
        {
            using namespace std;
            return sqrt( std::accumulate( std::ranges::begin( u ),
                std::ranges::end( u ),
                0.,
                []( auto sum, auto x )
                {
                    return sum + ::ponio::detail::power<2>( abs( x ) );
                } ) );
        }

        /**
         * @brief functor of power method to estimate spectral radius of an operator \f$f\f$
         *
         */
        struct power_method
        {
            /**
             * @brief implementation of power method
             *
             * @tparam problem_t type of operator \f$f\f$
             * @tparam value_t   type of coefficients
             * @tparam state_t   type of state
             * @param f          operator \f$f\f$
             * @param tn         current time where estimate spectral radius
             * @param un         current state where estimate spectral radius
             * @param dt         current time step (unused in power method)
             * @return value_t   estimation of spectral radius
             */
            template <typename problem_t, typename value_t, typename state_t, typename array_work_t>
            value_t
            operator()( problem_t&& f, value_t tn, state_t& un, [[maybe_unused]] value_t dt, array_work_t& du_work )
            {
                value_t eigmax  = 0.;
                value_t eigmaxo = 0.;

                auto& fn = du_work[0];
                auto& fz = du_work[1];
                auto& z  = du_work[2];

                std::forward<problem_t>( f )( tn, un, fn );
                fz = fn;

                std::forward<problem_t>( f )( tn, fz, z );

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
                    z    = quot * z;
                }
                else
                {
                    dzyn = sqrt_eps;
                    z    = dzyn; // TODO: change this line to iterate over z if needed (same trick as norm_2 I think)
                }

                bool necessary   = true;
                std::size_t iter = 0;

                static constexpr std::size_t max_iter  = 50;
                static constexpr value_t safety_factor = 1.2;

                // start power method
                while ( necessary )
                {
                    eigmaxo = eigmax;
                    std::forward<problem_t>( f )( tn, z, fz );
                    dfzfn  = detail::norm_2( static_cast<state_t>( fz - fn ) );
                    eigmax = safety_factor * dfzfn / dzyn;

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

                    using namespace std;
                    necessary = ( iter < max_iter ) && !( iter > 2 && abs( eigmax - eigmaxo ) <= 0.05 * eigmax );
                    ++iter;
                }

                return eigmax;
            }
        };

        /**
         * @brief computes degree of ROCK polynomial
         *
         * @tparam value_t     type of coefficient
         * @tparam rock_coeff_ structure of ROCK coefficients (ROCK2 or ROCK4)
         */
        template <typename value_t, typename rock_coeff_>
        struct degree_computer
        {
            using rock_coeff = rock_coeff_;

            /**
             * @brief computes optimal degree of ROCK polynomial
             *
             * @param mdeg                                 number of stages estimates by spectral radius
             * @return std::pair<std::size_t, std::size_t> shift indexes to read ROCK coefficients
             */
            static std::pair<std::size_t, std::size_t>
            optimal_degree( std::size_t& mdeg )
            {
                std::size_t mz = 1;
                std::size_t mr = 1;

                // TODO : verifier si ce test est utile
                // if ( mdeg < 2 )
                // {
                //     return { mz, mr };
                // }

                std::size_t i = 1;
                for ( auto ms_i : rock_coeff::ms )
                {
                    if ( ms_i / mdeg >= 1 )
                    {
                        mdeg = rock_coeff::ms[i - 1];
                        mz   = i;
                        break;
                    }
                    mr = mr + rock_coeff::ms[i - 1] * 2 - 1;

                    ++i;
                }

                return { mz, mr };
            }

            /**
             * @brief computes number of stages needed to stabilized spectral radius of \f$f\f$ and returns also number of evaluation of
             * function \f$f\f$
             *
             * @tparam rock_method    type selector for ROCK2 or ROCK4 method (with `rock_order::rock_2` or `rock_order::rock_4`)
             * @tparam eig_computer_t type of computer of spectral radius
             * @tparam problem_t      type of operator \f$f\f$
             * @tparam state_t        type of state
             * @param eig_computer computer of spectral radius
             * @param f            operator \f$f\f$
             * @param tn           current time
             * @param un           current state
             * @param dt           current time step
             * @param s_min        minimal number of stages (3 for ROCK2, 5 for ROCK4)
             */
            template <typename rock_method, typename eig_computer_t, typename problem_t, typename state_t, typename array_work_t>
            static std::tuple<std::size_t, std::size_t>
            compute_n_stages( rock_method,
                eig_computer_t&& eig_computer,
                problem_t& f,
                value_t tn,
                state_t& un,
                value_t& dt,
                array_work_t& du_work,
                std::size_t s_min )
            {
                std::size_t n_eval = 0;
                auto f_counter     = [&n_eval, &f]( value_t t, state_t& u, state_t& du )
                {
                    ++n_eval;
                    f( t, u, du );
                };

                value_t const eigmax = std::forward<eig_computer_t>( eig_computer )( f_counter, tn, un, dt, du_work );
                auto mdeg            = s_min;
                if constexpr ( std::same_as<rock_method, rock_order::rock_2> )
                {
                    mdeg = static_cast<std::size_t>( std::ceil( std::sqrt( ( 1.5 + dt * eigmax ) / 0.811 ) ) );
                    if ( mdeg > 200 )
                    {
                        mdeg = 200;
                        dt   = 0.8 * ( static_cast<double>( mdeg * mdeg ) * 0.811 - 1.5 ) / eigmax;
                    }

                    mdeg = std::max( mdeg, s_min ) - 2;
                }
                else if constexpr ( std::same_as<rock_method, rock_order::rock_4> )
                {
                    mdeg = static_cast<std::size_t>( ( std::sqrt( ( 3.0 + dt * eigmax ) / 0.353 ) ) ) + 1;
                    if ( mdeg > 152 )
                    {
                        mdeg = 152;
                        dt   = 0.8 * ( static_cast<double>( mdeg * mdeg ) * 0.353 - 3.0 ) / eigmax;
                    }

                    mdeg = std::max( mdeg, s_min ) - 4;
                }
                else
                {
                    static_assert( std::same_as<rock_method, rock_order::rock_2>, "Unknow ROCK method" );
                }

                return { mdeg, n_eval };
            }

            /**
             * @brief complete procedure to compute number of stages of ROCK method and shift indexes to read ROCK coefficients
             *
             * @tparam rock_method    type selector for ROCK2 or ROCK4 method (with `rock_order::rock_2` or `rock_order::rock_4`)
             * @tparam eig_computer_t type of computer of spectral radius
             * @tparam problem_t      type of operator \f$f\f$
             * @tparam state_t        type of state
             * @param eig_computer computer of spectral radius
             * @param f            operator \f$f\f$
             * @param tn           current time
             * @param un           current state
             * @param dt           current time step
             * @param s_min        minimal number of stages (3 for ROCK2, 5 for ROCK4)
             * @return std::tuple<std::size_t, std::size_t, std::size_t> tuple with number of stages of ROCK method, shift index for last
             * stages, shift index of ROCK stages
             */
            template <typename rock_method, typename eig_computer_t, typename problem_t, typename state_t, typename array_work_t>
            static std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>
            compute_n_stages_optimal_degree( rock_method,
                eig_computer_t&& eig_computer,
                problem_t& f,
                value_t tn,
                state_t& un,
                value_t& dt,
                array_work_t& du_work,
                std::size_t s_min = 3 )
            {
                auto [mdeg,
                    n_eval] = compute_n_stages( rock_method(), std::forward<eig_computer_t>( eig_computer ), f, tn, un, dt, du_work, s_min );
                auto [mz, mr] = optimal_degree( mdeg );

                return { mdeg, mz, mr, n_eval };
            }
        };
    } // namespace detail

    /** @class rock2_impl
     *  @brief define ROCK2 method
     *
     *  @tparam eig_computer_t type of computer of maximal eigenvalue
     *  @tparam _is_embedded   define if method is used as adaptive or constant time step method [default is false]
     *  @tparam _value_t       type of coefficients
     */
    template <typename eig_computer_t, bool _is_embedded = false, typename _value_t = double>
    struct rock2_impl
    {
        static constexpr bool is_embedded      = _is_embedded;
        static constexpr std::size_t N_stages  = stages::dynamic;
        static constexpr std::size_t N_storage = 4;
        static constexpr std::size_t order     = 2;
        static constexpr std::string_view id   = "ROCK2";

        using value_t         = _value_t;
        using rock_coeff      = rock2_coeff<value_t>;
        using degree_computer = detail::degree_computer<value_t, rock_coeff>;

        iteration_info<rock2_impl> _info;

        eig_computer_t eig_computer;

        rock2_impl()
            : _info( default_config::tol, default_config::tol )
            , eig_computer( detail::power_method() )
        {
        }

        /**
         * @brief Construct a new ROCK2 algorithm object
         *
         * @param _eig_computer estimator of spectral radius
         */
        rock2_impl( eig_computer_t&& _eig_computer )
            : _info( default_config::tol, default_config::tol )
            , eig_computer( std::forward<eig_computer_t>( _eig_computer ) )
        {
        }

        /**
         * @brief computes an error for adaptive time step method
         *
         * @tparam state_t type of current state
         * @param unp1     estimation of solution at time \f$t^{n+1} = t^n+\Delta t\f$
         * @param un       solution at time \f$t^n\f$
         * @param tmp      \f$-\Delta t \sigma(1 - \frac{\tau}{\sigma^2})(f(u_{s-1}) - f(u_{s-2}))\f$
         * @return auto    estimation of error to compare to 1
         */
        template <typename state_t>
        auto
        error( state_t&& unp1, state_t&& un, state_t&& tmp )
        {
            using namespace std;
            return abs( std::forward<state_t>( tmp )
                        / ( info().absolute_tolerance
                            + info().relative_tolerance * max( abs( std::forward<state_t>( unp1 ) ), abs( std::forward<state_t>( un ) ) ) ) );
        }

        // same with ranges
        template <typename state_t>
            requires std::ranges::range<state_t>
        auto
        error( state_t&& unp1, state_t&& un, state_t&& tmp )
        {
            auto it_un  = std::ranges::begin( std::forward<state_t>( un ) );
            auto it_tmp = std::ranges::begin( std::forward<state_t>( tmp ) );

            using namespace std;
            return sqrt( std::accumulate( std::ranges::begin( std::forward<state_t>( unp1 ) ),
                             std::ranges::end( std::forward<state_t>( unp1 ) ),
                             0.,
                             [&]( auto sum, auto unp1_i )
                             {
                                 return sum + ::ponio::detail::power<2>( error( unp1_i, *it_un++, *it_tmp++ ) );
                             } )
                         / static_cast<value_t>( std::size( unp1 ) ) );
        }

        // same with something which contains a range
        template <typename state_t>
            requires ::ponio::detail::has_array_range<state_t>
        auto
        error( state_t&& unp1, state_t&& un, state_t&& tmp )
        {
            return error( std::forward<state_t>( unp1 ).array(), std::forward<state_t>( un ).array(), std::forward<state_t>( tmp ).array() );
        }

        /**
         * @brief iteration of ROCK2 method
         *
         * @tparam problem_t  type of \f$f\f$
         * @tparam state_t    type of current state
         * @tparam array_ki_t type of temporary stages (only 3 needed for ROCK2)
         * @param f  operator \f$f\f$
         * @param tn current time
         * @param un current state
         * @param G  array of temporary stages
         * @param dt current time step
         */
        template <typename problem_t, typename state_t, typename array_ki_t>
        // std::tuple<value_t, state_t, value_t>
        void
        operator()( problem_t& f, value_t& tn, state_t& un, array_ki_t& G, value_t& dt, state_t& unp1 )
        {
            _info.reset_eval();

            auto [mdeg,
                deg_index,
                start_index,
                n_eval] = degree_computer::compute_n_stages_optimal_degree( rock_order::rock_2(), eig_computer, f, tn, un, dt, G );

            _info.number_of_stages = mdeg + 2;
            _info.number_of_eval   = n_eval + mdeg + 2;

            auto& uj    = G[0];
            auto& ujm1  = G[1];
            auto& ujm2  = G[2];
            auto& f_tmp = G[3];

            uj   = un;
            ujm2 = un;

            value_t const mu1 = rock_coeff::recf[start_index - 1];

            value_t t_jm1 = tn + dt * mu1;
            value_t t_jm2 = tn + dt * mu1;
            value_t t_jm3 = tn;

            f( tn, un, f_tmp );
            ujm1 = un + dt * mu1 * f_tmp;

            if ( mdeg < 2 )
            {
                uj = ujm1;
            }
            // std::cout << "\nmdeg = " << mdeg << "\n";
            for ( std::size_t j = 2; j < mdeg + 1; ++j )
            {
                value_t const mu    = rock_coeff::recf[start_index + 2 * ( j - 2 ) + 1 - 1];
                value_t const kappa = rock_coeff::recf[start_index + 2 * ( j - 2 ) + 2 - 1];
                value_t const nu    = -1.0 - kappa;

                f( t_jm1, ujm1, f_tmp );
                uj = dt * mu * f_tmp - nu * ujm1 - kappa * ujm2;

                t_jm1 = dt * mu - nu * t_jm2 - kappa * t_jm3;

                if ( j < mdeg )
                {
                    std::swap( ujm2, ujm1 );
                    std::swap( ujm1, uj );
                }

                t_jm3 = t_jm2;
                t_jm2 = t_jm1;
            }

            // the two-stages finish procedure

            value_t const delta_t_1 = dt * rock_coeff::fp1[deg_index - 1]; // equals to $\Delta t \sigma$
            value_t const delta_t_2 = dt * rock_coeff::fp2[deg_index - 1]; // equals to $-\Delta t \sigma(1 - \frac{\tau}{\sigma^2})$

            f( t_jm1, uj, ujm2 );
            ujm1 = uj + delta_t_1 * ujm2;

            t_jm1 = t_jm1 + delta_t_1;

            if constexpr ( is_embedded )
            {
                f( t_jm1, ujm1, uj );
                f_tmp = delta_t_2 * ( uj - ujm2 );

                uj = ujm1 + delta_t_1 * uj + f_tmp;

                _info.error   = error( uj, un, f_tmp );
                _info.success = _info.error < 1.0;
                _info.number_of_eval += 1;

                value_t fac    = std::min( 2.0, std::max( 0.5, std::sqrt( 1.0 / _info.error ) ) );
                value_t new_dt = 0.8 * fac * dt;

                // accepted step
                if ( _info.success )
                {
                    // return { tn + dt, uj, new_dt };

                    tn = tn + dt;
                    std::swap( uj, unp1 );
                    dt = new_dt;
                }
                else
                {
                    // return { tn, un, new_dt };

                    // tn = tn;
                    std::swap( un, unp1 );
                    dt = new_dt;
                }
            }
            else
            {
                f( t_jm1, ujm1, f_tmp );
                // uj = ujm1 + ( delta_t_1 + delta_t_2 ) * f_tmp - delta_t_2 * ujm2;

                tn   = tn + dt;
                unp1 = ujm1 + ( delta_t_1 + delta_t_2 ) * f_tmp - delta_t_2 * ujm2;

                // return { tn + dt, uj, dt };
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

    /**
     * @brief helper to build a `rock2_impl` object
     *
     * @tparam is_embedded    define if method is used as adaptive or constant time step method [default is false]
     * @tparam value_t        type of coefficients
     * @tparam eig_computer_t type of computer of maximal eigenvalue
     * @param eig_computer    computer of maximal eigenvalues
     */
    template <bool is_embedded = false, typename value_t = double, typename eig_computer_t>
    auto
    rock2( eig_computer_t&& eig_computer )
    {
        return rock2_impl<eig_computer_t, is_embedded, value_t>( std::forward<eig_computer_t>( eig_computer ) );
    }

    /**
     * @brief helper to build a `rock2_impl` object with power method to compute maximal eigen value
     *
     * @tparam is_embedded    define if method is used as adaptive or constant time step method [default is false]
     * @tparam value_t        type of coefficients
     */
    template <bool is_embedded = false, typename value_t = double>
    auto
    rock2()
    {
        return rock2<is_embedded, value_t>( detail::power_method() );
    }

    /** @class rock4_impl
     *  @brief define ROCK4 method
     *
     *  @tparam eig_computer_t type of computer of maximal eigenvalue
     *  @tparam _is_embedded   define if method is used as adaptive or constant time step method [default is false]
     *  @tparam _value_t       type of coefficients
     */
    template <typename eig_computer_t, bool _is_embedded = false, typename _value_t = double>
    struct rock4_impl
    {
        static constexpr bool is_embedded      = _is_embedded;
        static constexpr std::size_t N_stages  = stages::dynamic;
        static constexpr std::size_t N_storage = 7;
        static constexpr std::size_t order     = 4;
        static constexpr std::string_view id   = "ROCK4";

        using value_t         = _value_t;
        using rock_coeff      = rock4_coeff<value_t>;
        using degree_computer = detail::degree_computer<value_t, rock_coeff>;

        iteration_info<rock4_impl> _info;

        eig_computer_t eig_computer;

        rock4_impl()
            : _info( default_config::tol, default_config::tol )
            , eig_computer( detail::power_method() )
        {
        }

        /**
         * @brief Construct a new ROCK4 algorithm object
         *
         * @param _eig_computer estimator of spectral radius
         */
        rock4_impl( eig_computer_t&& _eig_computer )
            : _info( default_config::tol, default_config::tol )
            , eig_computer( std::forward<eig_computer_t>( _eig_computer ) )
        {
        }

        /**
         * @brief computes an error for adaptive time step method
         *
         * @tparam state_t type of current state
         * @param unp1     estimation of solution at time \f$t^{n+1} = t^n+\Delta t\f$
         * @param tmp      other estimation of solution at time \f$t^{n+1} = t^n+\Delta t\f$
         * @return auto    estimation of error to compare to 1
         */
        template <typename state_t>
        auto
        error( state_t const& unp1, state_t const& tmp )
        {
            using namespace std;
            return abs( tmp / ( info().absolute_tolerance + info().relative_tolerance * abs( unp1 ) ) );
        }

        // same with ranges
        template <typename state_t>
            requires std::ranges::range<state_t>
        auto
        error( state_t const& unp1, state_t const& tmp )
        {
            auto it_tmp = std::ranges::begin( tmp );

            using namespace std;
            return sqrt( std::accumulate( std::ranges::begin( unp1 ),
                             std::ranges::end( unp1 ),
                             0.,
                             [&]( auto sum, auto unp1_i )
                             {
                                 return sum + ::ponio::detail::power<2>( error( unp1_i, *it_tmp++ ) );
                             } )
                         / static_cast<value_t>( std::size( unp1 ) ) );
        }

        // same with something which contains a range
        template <typename state_t>
            requires ::ponio::detail::has_array_range<state_t>
        auto
        error( state_t const& unp1, state_t const& tmp )
        {
            return error( unp1.array(), tmp.array() );
        }

        /**
         * @brief iteration of ROCK4 method
         *
         * @tparam problem_t  type of \f$f\f$
         * @tparam state_t    type of current state
         * @tparam array_ki_t type of temporary stages (only 6 needed for ROCK4)
         * @param f     operator \f$f\f$
         * @param tn    current time
         * @param un    current state
         * @param G     array of temporary stages
         * @param dt    current time step
         * @param u_np1 solution \f$u^{n+1}\f$ at time \f$t^{n+1} = t^n + \Delta t\f$
         */
        template <typename problem_t, typename state_t, typename array_ki_t>
        void
        operator()( problem_t& f, value_t& tn, state_t& un, array_ki_t& G, value_t& dt, state_t& unp1 )
        {
            _info.reset_eval();

            auto [mdeg,
                deg_index,
                start_index,
                n_eval] = degree_computer::compute_n_stages_optimal_degree( rock_order::rock_4(), eig_computer, f, tn, un, dt, G, 5 );

            _info.number_of_stages = mdeg + 4;
            _info.number_of_eval   = n_eval + mdeg + 4;

            auto& uj    = G[0];
            auto& ujm1  = G[1];
            auto& ujm2  = G[2];
            auto& ujm3  = G[3];
            auto& ujm4  = G[4];
            auto& f_tmp = G[5];

            uj   = un;
            ujm2 = un;

            value_t const mu1 = dt * rock_coeff::recf[start_index - 1];
            value_t t_jm1     = tn + mu1;
            value_t t_jm2     = tn + mu1;
            value_t t_jm3     = tn;

            f( tn, un, f_tmp );
            ujm1 = un + mu1 * f_tmp;

            if ( mdeg < 2 )
            {
                uj = ujm1;
            }

            for ( std::size_t j = 2; j < mdeg + 1; ++j )
            {
                value_t const mu    = rock_coeff::recf[start_index + 2 * ( j - 2 ) + 1 - 1];
                value_t const kappa = rock_coeff::recf[start_index + 2 * ( j - 2 ) + 2 - 1];
                value_t const nu    = -1.0 - kappa;

                f( t_jm1, ujm1, f_tmp );
                uj = dt * mu * f_tmp - nu * ujm1 - kappa * ujm2;

                t_jm1 = dt * mu + nu * t_jm2 - kappa * t_jm3;

                if ( j < mdeg )
                {
                    std::swap( ujm2, ujm1 );
                    std::swap( ujm1, uj );
                }

                t_jm3 = t_jm2;
                t_jm2 = t_jm1;
            }

            // the fourth-stages finish procedure

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
            f( t_jm1, uj, ujm1 );
            ujm3 = uj + a_21 * ujm1;

            // stage 2.
            t_jm2 = t_jm1 + a_21;
            f( t_jm2, ujm3, ujm2 );
            ujm4 = uj + a_31 * ujm1 + a_32 * ujm2;

            // stage 3.
            t_jm2 = t_jm1 + a_31 + a_32;
            f( t_jm2, ujm4, ujm3 );
            ujm4 = uj + a_41 * ujm1 + a_42 * ujm2 + a_43 * ujm3;

            // stage 4.
            t_jm2 = t_jm1 + a_41 + a_42 + a_43;

            if constexpr ( is_embedded )
            {
                auto& tmp = G[6];

                // for embedded method for error estimation
                value_t const bh_1 = dt * ( rock_coeff::fpbe[deg_index - 1][0] - rock_coeff::fpb[deg_index - 1][0] );
                value_t const bh_2 = dt * ( rock_coeff::fpbe[deg_index - 1][1] - rock_coeff::fpb[deg_index - 1][1] );
                value_t const bh_3 = dt * ( rock_coeff::fpbe[deg_index - 1][2] - rock_coeff::fpb[deg_index - 1][2] );
                value_t const bh_4 = dt * ( rock_coeff::fpbe[deg_index - 1][3] - rock_coeff::fpb[deg_index - 1][3] );
                value_t const bh_5 = dt * rock_coeff::fpbe[deg_index - 1][4];

                f( t_jm2, ujm4, tmp );
                uj = uj + b_1 * ujm1 + b_2 * ujm2 + b_3 * ujm3 + b_4 * tmp;

                f( t_jm2, uj, f_tmp );
                tmp = bh_1 * ujm1 + bh_2 * ujm2 + bh_3 * ujm3 + bh_4 * tmp + bh_5 * f_tmp;

                _info.error   = error( uj, tmp );
                _info.success = _info.error < 1.0;
                _info.number_of_eval += 1; // one of two evaluations already count

                value_t fac = std::pow( ( 1. / _info.error ), 0.25 );
                fac         = std::min( 5., std::max( 0.1, 0.8 * fac ) );

                value_t new_dt = fac * dt;

                // accepted step
                if ( _info.success )
                {
                    tn = tn + dt;
                    std::swap( uj, unp1 );
                    dt = new_dt;
                }
                else
                {
                    // tn = tn;
                    std::swap( un, unp1 );
                    dt = new_dt;
                }
            }
            else
            {
                f( t_jm2, ujm4, f_tmp );

                tn   = tn + dt;
                unp1 = uj + b_1 * ujm1 + b_2 * ujm2 + b_3 * ujm3 + b_4 * f_tmp;
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

    /**
     * @brief helper to build a `rock4_impl` object
     *
     * @tparam is_embedded    define if method is used as adaptive or constant time step method [default is false]
     * @tparam value_t        type of coefficients
     * @tparam eig_computer_t type of computer of maximal eigenvalue
     * @param eig_computer    computer of maximal eigenvalues
     */
    template <bool is_embedded = false, typename value_t = double, typename eig_computer_t>
    auto
    rock4( eig_computer_t&& eig_computer )
    {
        return rock4_impl<eig_computer_t, is_embedded, value_t>( std::forward<eig_computer_t>( eig_computer ) );
    }

    /**
     * @brief helper to build a `rock4_impl` object with power method to compute maximal eigen value
     *
     * @tparam is_embedded    define if method is used as adaptive or constant time step method [default is false]
     * @tparam value_t        type of coefficients
     */
    template <bool is_embedded = false, typename value_t = double>
    auto
    rock4()
    {
        return rock4<is_embedded, value_t>( detail::power_method() );
    }

} // namespace ponio::runge_kutta::rock
