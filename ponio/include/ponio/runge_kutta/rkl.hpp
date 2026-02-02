// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private

#pragma once

#include <cstddef>
#include <string_view> // NOLINT(misc-include-cleaner)

#include "../iteration_info.hpp"
#include "../stage.hpp" // NOLINT(misc-include-cleaner)

namespace ponio::runge_kutta::legendre
{

    /** @class explicit_rkl1
     * @brief define RKL1 (Runge-Kutta-Legendre method of ordre 1) with user defined number of stages
     *
     * @tparam N_stages_ number of stages
     * @tparam _value_t type of coefficients
     *
     * @warning the method is only presented with an autonomous problem, ie \f$\dot{y} = f(y)\f$
     */
    template <std::size_t N_stages_, typename _value_t = double>
    struct explicit_rkl1
    {
        static_assert( N_stages_ > 0, "Number of stages should be at least 1 in eRKL1" );
        static constexpr std::size_t N_stages = N_stages_;
        static constexpr std::size_t order    = 1;
        static constexpr std::string_view id  = "RKL1";
        static constexpr bool is_embedded     = false;
        using value_t                         = _value_t;

        iteration_info<explicit_rkl1> _info;

        /**
         * @brief compute \f$\mu_j = \frac{2j-1}{j}\f$
         *
         * @tparam j index \f$j\f$
         */
        template <std::size_t j>
        static constexpr value_t
        mu()
        {
            return static_cast<value_t>( 2 * j - 1 ) / static_cast<value_t>( j );
        }

        /**
         * @brief compute \f$\nu_j = \frac{1-j}{j}\f$
         *
         * @tparam j index \f$j\f$
         */
        template <std::size_t j>
        static constexpr value_t
        nu()
        {
            return ( 1 - static_cast<value_t>( j ) ) / static_cast<value_t>( j );
        }

        /**
         * @brief compute \f$\tilde{\mu}_j = \frac{2j-1}{j}\frac{2}{s^2 + s}\f$ with \f$s\f$ the number of stages
         *
         * @tparam j index \f$j\f$
         */
        template <std::size_t j>
        static constexpr value_t
        mu_t()
        {
            return static_cast<value_t>( 2 * j - 1 ) / static_cast<value_t>( j ) * 2.
                 / static_cast<value_t>( N_stages * N_stages + N_stages );
        }

        explicit_rkl1()
        {
            _info.number_of_eval = N_stages;
        }

        /**
         * @brief compute stage \f$j\f$ of RKL1 method
         *
         * @tparam problem_t  type of operator \f$f\f$
         * @tparam state_t    type of current state
         * @tparam array_ki_t type of array of temporary stages
         * @tparam j          integer of stage
         * @param f  operator \f$f\f$
         * @param tn current time
         * @param yn current state
         * @param Yj array of temporary stages
         * @param dt current time step
         * @param ui temporary step
         * @param yi computed output step
         *
         * @details \f$y^{(j)} = \mu_j y^{(j-1)} + \nu_j y^{(j-2)} + \tilde{\mu}_j \Delta t f(t^n, y^{(j-1)})\f$
         */
        template <typename problem_t, typename state_t, typename array_ki_t, std::size_t j>
        void
        stage( Stage<j>, problem_t& f, value_t tn, state_t&, array_ki_t const& Yj, value_t dt, state_t& ui, state_t& yi )
        {
            f( tn, Yj[j - 1], ui );
            yi = mu<j>() * Yj[j - 1] + nu<j>() * Yj[j - 2] + mu_t<j>() * dt * ui; // be careful Yj[j] is y^{(j+1)}
        }

        /**
         * @brief compute stage \f$0\f$ of RKL1 method
         *
         * @tparam problem_t  type of operator \f$f\f$
         * @tparam state_t    type of current state
         * @tparam array_ki_t type of array of temporary stages
         * @param yn current state
         * @param yi computed output step
         *
         * @details \f$y^{(0)} = y^n\f$
         */
        template <typename problem_t, typename state_t, typename array_ki_t>
        void
        stage( Stage<0>, problem_t&, value_t, state_t& yn, array_ki_t const&, value_t, state_t&, state_t& yi )
        {
            yi = yn;
        }

        /**
         * @brief compute stage \f$1\f$ of RKL1 method
         *
         * @tparam problem_t  type of operator \f$f\f$
         * @tparam state_t    type of current state
         * @tparam array_ki_t type of array of temporary stages
         * @param f  operator \f$f\f$
         * @param tn current time
         * @param yn current state
         * @param dt current time step
         * @param ui temporary step
         * @param yi computed output step
         *
         * @details \f$y^{(1)} = y^n + \tilde{\mu}_1 \Delta t f(t^n, y^n)\f$
         */
        template <typename problem_t, typename state_t, typename array_ki_t>
        void
        stage( Stage<1>, problem_t& f, value_t tn, state_t& yn, array_ki_t const&, value_t dt, state_t& ui, state_t& yi )
        {
            f( tn, yn, ui );
            yi = yn + mu_t<1>() * dt * ui;
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

    namespace details
    {
        /**
         * @brief compute \f$b_j\f$ for \f$j\geq 2\f$
         *
         * @tparam value_t type of coefficient
         * @tparam j       index \f$j\f$
         *
         * @details \f$b_j = \frac{j^2 + j - 2}{2j(j+1)}\f$
         */
        template <typename value_t, std::size_t j>
        struct b
        {
            static constexpr value_t value = static_cast<value_t>( j * j + j - 2 ) / static_cast<value_t>( 2 * j * ( j + 1 ) );
        };

        /**
         * @brief specialization of \f$b_j\f$ for \f$j=0\f$
         *
         * @tparam value_t type of coefficient
         *
         * @details \f$b_0 = \frac{1}{3}\f$
         */
        template <typename value_t>
        struct b<value_t, 0>
        {
            static constexpr value_t value = static_cast<value_t>( 1. / 3. );
        };

        /**
         * @brief specialization of \f$b_j\f$ for \f$j=1\f$
         *
         * @tparam value_t type of coefficient
         *
         * @details \f$b_1 = \frac{1}{3}\f$
         */
        template <typename value_t>
        struct b<value_t, 1>
        {
            static constexpr value_t value = static_cast<value_t>( 1. / 3. );
        };

        template <typename value_t, std::size_t j>
        constexpr value_t b_v = b<value_t, j>::value;

        /**
         * @brief compute \f$a_j\f$ coefficient
         *
         * @tparam value_t type of coefficient
         * @tparam j       index \f$j\f$
         *
         * @details \f$a_j = 1 - b_j\f$
         */
        template <typename value_t, std::size_t j>
        struct a
        {
            static constexpr value_t value = static_cast<value_t>( 1. ) - b_v<value_t, j>;
        };

        template <typename value_t, std::size_t j>
        constexpr value_t a_v = a<value_t, j>::value;
    } // namespace details

    /** @class explicit_rkl2
     * @brief define RKL2 (Runge-Kutta-Legendre method of ordre 2) with user defined number of stages
     *
     * @tparam N_stages_ number of stages
     * @tparam _value_t type of coefficients
     *
     * @warning the method is only presented with an autonomous problem, ie \f$\dot{y} = f(y)\f$
     */
    template <std::size_t N_stages_, typename _value_t = double>
    struct explicit_rkl2
    {
        static_assert( N_stages_ > 0, "Number of stages should be at least 2 in eRKL2" );
        static constexpr std::size_t N_stages = N_stages_;
        static constexpr std::size_t order    = 2;
        static constexpr std::string_view id  = "RKL2";
        static constexpr bool is_embedded     = false;
        using value_t                         = _value_t;

        iteration_info<explicit_rkl2> _info;

        /**
         * @brief compute \f$w_1\f$ coefficient
         *
         *
         * @details \f$w_1 = \frac{4}{s^2 + s -2}\f$ with \f$s\f$ the number of stages of the method
         */
        static constexpr value_t
        w1()
        {
            return static_cast<value_t>( 4. ) / static_cast<value_t>( N_stages * N_stages + N_stages - 2 );
        }

        /**
         * @brief compute \f$\mu_j\f$ coefficient
         *
         * @tparam j       index \f$j\f$
         *
         * @details \f$b_1 = \frac{1}{3}\f$
         */
        template <std::size_t j>
        static constexpr value_t
        mu()
        {
            return static_cast<value_t>( 2 * j - 1 ) * details::b_v<value_t, j> / static_cast<value_t>( j * details::b_v<value_t, j - 1> );
        }

        /**
         * @brief compute \f$\mu_j\f$ coefficient
         *
         * @tparam j       index \f$j\f$
         *
         * @details \f$mu_j = \frac{2j-1}{j}\frac{b_j}{b_{j-1}}\f$
         */
        template <std::size_t j>
        static constexpr value_t
        nu()
        {
            return -static_cast<value_t>( ( j - 1 ) * details::b_v<value_t, j> ) / static_cast<value_t>( j * details::b_v<value_t, j - 2> );
        }

        /**
         * @brief compute \f$\tilde{\mu}_j\f$ coefficient
         *
         * @tparam j       index \f$j\f$
         *
         * @details \f$\tilde{\mu}_j = \mu_j w_1\f$ for \f$1<j\f$, \f$\tilde{\mu}_1 - b_1w_1\f$
         */
        template <std::size_t j>
        static constexpr value_t
        mu_t()
        {
            if constexpr ( j == 1 )
            {
                return details::b_v<value_t, 1> * w1();
            }
            else
            {
                return mu<j>() * w1();
            }
        }

        /**
         * @brief compute \f$\gamma_j\f$ coefficient
         *
         * @tparam j       index \f$j\f$
         *
         * @details \f$gamma_j = -a_{j-1}\tilde{\mu}_j\f$
         */
        template <std::size_t j>
        static constexpr value_t
        gamma_t()
        {
            return -details::a_v<value_t, j - 1> * mu_t<j>();
        }

        explicit_rkl2()
        {
            _info.number_of_eval = N_stages;
        }

        /**
         * @brief compute stage \f$j\f$ of RKL2 method
         *
         * @tparam problem_t  type of operator \f$f\f$
         * @tparam state_t    type of current state
         * @tparam array_ki_t type of array of temporary stages
         * @tparam j          integer of stage
         * @param f  operator \f$f\f$
         * @param tn current time
         * @param yn current state
         * @param Yj array of temporary stages
         * @param dt current time step
         * @param ui temporary step
         * @param yi computed output step
         *
         * @details \f$y^{(j)} = \mu_j y^{(j-1)} + \nu_j y^{(j-2)} + (1-\mu_j-\nu_j)y^{(0)} + \tilde{\mu}_j \Delta t f(t^n, y^{(j-1)}) +
         * \gamma_j\Delta t f(t^n, y^{(0)})\f$
         */
        template <typename problem_t, typename state_t, typename array_ki_t, std::size_t j>
        void
        stage( Stage<j>, problem_t& f, value_t tn, state_t& yn, array_ki_t const& Yj, value_t dt, state_t& ui, state_t& yi )
        {
            f( tn, Yj[j - 1], ui );
            yi = mu<j>() * Yj[j - 1] + nu<j>() * Yj[j - 2] + ( 1. - mu<j>() - nu<j>() ) * yn + mu_t<j>() * dt * ui + gamma_t<j>() * Yj[0];
        }

        /**
         * @brief compute \f$\Delta t f(t^n, y^n)\f$ term (needs for all other stages)
         *
         * @tparam problem_t  type of operator \f$f\f$
         * @tparam state_t    type of current state
         * @tparam array_ki_t type of array of temporary stages
         * @param f  operator \f$f\f$
         * @param tn current time
         * @param yn current state
         * @param dt current time step
         * @param ui temporary step
         * @param yi computed output step
         *
         * @warning first stage of RKL2 method is \f$y^n\f$ but here we compute \f$\Delta t f(t^n, y^n)\f$ term
         */
        template <typename problem_t, typename state_t, typename array_ki_t>
        void
        stage( Stage<0>, problem_t& f, value_t tn, state_t& yn, array_ki_t const&, value_t dt, state_t& ui, state_t& yi )
        {
            f( tn, yn, ui );
            yi = dt * ui; // be careful Yj[0] is dt*f(tn, yn)
        }

        /**
         * @brief compute stage \f$1\f$ of RKL2 method
         *
         * @tparam problem_t  type of operator \f$f\f$
         * @tparam state_t    type of current state
         * @tparam array_ki_t type of array of temporary stages
         * @param yn current state
         * @param Yj array of temporary stages
         * @param yi computed output step
         *
         * @details \f$y^{(1)} = y^n + \tilde{\mu}_1 \Delta t f(t^n, y^n)\f$
         */
        template <typename problem_t, typename state_t, typename array_ki_t>
        void
        stage( Stage<1>, problem_t&, value_t, state_t& yn, array_ki_t const& Yj, value_t, state_t&, state_t& yi )
        {
            yi = yn + mu_t<1>() * Yj[0];
        }

        /**
         * @brief compute stage \f$2\f$ of RKL2 method
         *
         * @tparam problem_t  type of operator \f$f\f$
         * @tparam state_t    type of current state
         * @tparam array_ki_t type of array of temporary stages
         * @param f  operator \f$f\f$
         * @param tn current time
         * @param yn current state
         * @param Yj array of temporary stages
         * @param dt current time step
         * @param ui temporary step
         * @param yi computed output step
         *
         * @details \f$y^{(2)} = \mu_2 y^{(1)} + \nu_2y^n + (1-\mu_2-\nu_2)y^n + \tilde{\mu}_2\Delta tf(t^n, y^{(1)}) + \gamma_2\Delta t
         * f(t^n, y^n)\f$
         */
        template <typename problem_t, typename state_t, typename array_ki_t>
        void
        stage( Stage<2>, problem_t& f, value_t tn, state_t& yn, array_ki_t const& Yj, value_t dt, state_t& ui, state_t& yi )
        {
            f( tn, Yj[1], ui );
            yi = mu<2>() * Yj[1] + nu<2>() * yn + ( 1. - mu<2>() - nu<2>() ) * yn + mu_t<2>() * dt * ui + gamma_t<2>() * Yj[0];
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

} // namespace ponio::runge_kutta::legendre
