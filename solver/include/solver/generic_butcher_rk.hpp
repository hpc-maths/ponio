// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <concepts>

#include "butcher_tableau.hpp"
#include "detail.hpp"
#include "ponio_config.hpp"
#include "stage.hpp"

namespace ode::butcher
{

    namespace runge_kutta
    {

        template <typename Tableau>
        struct explicit_rk_butcher
        {
            Tableau butcher;
            static constexpr std::size_t N_stages = Tableau::N_stages;
            static constexpr bool is_embedded     = is_embedded_tableau<Tableau>;
            static constexpr std::size_t order    = Tableau::order;
            static constexpr char const* id       = Tableau::id;

            explicit_rk_butcher(double tol_ = ponio::default_config::tol)
                : butcher()
                , tol(tol_)
            {
            }

            template <typename Problem_t, typename state_t, typename value_t, typename ArrayKi_t, std::size_t I>
            inline state_t stage(Stage<I>, Problem_t& f, value_t tn, state_t const& un, ArrayKi_t const& Ki, value_t dt)
            {
                return f(tn + butcher.c[I] * dt, ::detail::tpl_inner_product<I>(butcher.A[I], Ki, un, dt));
            }

            template <typename Problem_t, typename state_t, typename value_t, typename ArrayKi_t>
            inline state_t stage(Stage<N_stages>, Problem_t& f, value_t tn, state_t const& un, ArrayKi_t const& Ki, value_t dt)
            {
                return ::detail::tpl_inner_product<N_stages>(butcher.b, Ki, un, dt);
            }

            template <typename Problem_t, typename state_t, typename value_t, typename ArrayKi_t, typename Tab = Tableau>
                requires std::same_as<Tab, Tableau> && is_embedded
            inline state_t stage(Stage<N_stages + 1>, Problem_t& f, value_t tn, state_t const& un, ArrayKi_t const& Ki, value_t dt)
            {
                return ::detail::tpl_inner_product<N_stages>(butcher.b2, Ki, un, dt);
            }

            double tol;
        };

    } // namespace runge_kutta

    namespace lawson
    {

        template <typename Exp_t>
        struct lawson_base
        {
            Exp_t m_exp;

            lawson_base(Exp_t exp_)
                : m_exp(exp_)
            {
            }
        };

        template <typename Tableau, typename Exp_t>
        struct explicit_rk_butcher : public lawson_base<Exp_t>
        {
            using lawson_base<Exp_t>::m_exp;

            Tableau butcher;
            static constexpr std::size_t N_stages = Tableau::N_stages;
            static constexpr bool is_embedded     = is_embedded_tableau<Tableau>;
            static constexpr std::size_t order    = Tableau::order;
            static constexpr char const* id       = Tableau::id;

            explicit_rk_butcher(Exp_t exp_, double tol_ = ponio::default_config::tol)
                : lawson_base<Exp_t>(exp_)
                , butcher()
                , tol(tol_)
            {
            }

            template <typename Problem_t, typename state_t, typename value_t, typename ArrayKi_t, std::size_t i>
            inline state_t stage(Stage<i>, Problem_t& pb, value_t tn, state_t const& un, ArrayKi_t const& Ki, value_t dt)
            {
                return m_exp(-butcher.c[i] * dt * pb.l)
                     * pb.n(tn + butcher.c[i] * dt,
                         m_exp(butcher.c[i] * dt * pb.l) * ::detail::tpl_inner_product<i>(butcher.A[i], Ki, un, dt));
            }

            template <typename Problem_t, typename state_t, typename value_t, typename ArrayKi_t>
            inline state_t stage(Stage<N_stages>, Problem_t& pb, value_t tn, state_t const& un, ArrayKi_t const& Ki, value_t dt)
            {
                return m_exp(dt * pb.l) * ::detail::tpl_inner_product<N_stages>(butcher.b, Ki, un, dt);
            }

            template <typename Problem_t, typename state_t, typename value_t, typename ArrayKi_t, typename Tab = Tableau>
                requires std::same_as<Tab, Tableau> && is_embedded
            inline state_t stage(Stage<N_stages + 1>, Problem_t& pb, value_t tn, state_t const& un, ArrayKi_t const& Ki, value_t dt)
            {
                return m_exp(dt * pb.l) * ::detail::tpl_inner_product<N_stages>(butcher.b2, Ki, un, dt);
            }

            double tol;
        };

        /**
         * factory of Lawson method to help clang++ which doesn't support C++20
         *
         * @tparam Tableau type of Butcher tableau
         * @tparam Exp_t type of exponential function
         * @param exp_ exponential function
         * @param tol_ tolenrence of method for adaptative time step integrator
         */
        template <typename Tableau, typename Exp_t>
        auto make_lawson(Exp_t exp_, double tol_ = ponio::default_config::tol)
        {
            return explicit_rk_butcher<Tableau, Exp_t>(exp_, tol_);
        }

    } // namespace lawson

    namespace exp_runge_kutta
    {

        namespace detail
        {
            template <typename func_t, typename linear_t>
                requires std::invocable<func_t, linear_t>
            auto coefficient_eval(func_t&& f, linear_t&& l)
            {
                return f(std::move(l));
            }

            template <typename value_t, typename linear_t>
                requires std::is_arithmetic_v<std::remove_cvref_t<value_t>>
            auto coefficient_eval(value_t&& val, linear_t&&)
            {
                return val;
            }

            template <std::size_t I, std::size_t J, typename tuple_t>
            auto triangular_get(tuple_t& t)
            {
                return std::get<I*(I - 1) / 2 + J>(t);
            }

            template <std::size_t I, typename state_t, typename value_t, typename linear_t, typename tuple_t, typename array_t, std::size_t... Is>
            constexpr state_t tpl_inner_product_impl(tuple_t const& a,
                array_t const& k,
                state_t const& init,
                linear_t const& l,
                value_t mul_coeff,
                std::index_sequence<Is...>)
            {
                return (init + ... + (mul_coeff * coefficient_eval(triangular_get<I, Is>(a), mul_coeff * l) * (k[Is] + l * init)));
            }

            template <std::size_t I, typename state_t, typename value_t, typename linear_t, typename tuple_t, typename array_t>
            constexpr state_t tpl_inner_product(tuple_t const& a, array_t const& k, state_t const& init, linear_t const& l, value_t mul_coeff)
            {
                return tpl_inner_product_impl<I>(a, k, init, l, mul_coeff, std::make_index_sequence<I>());
            }

            template <typename state_t, typename value_t, typename linear_t, typename tuple_t, typename array_t, std::size_t... Is>
            constexpr state_t tpl_inner_product_b_impl(tuple_t const& b,
                array_t const& k,
                state_t const& init,
                linear_t const& l,
                value_t mul_coeff,
                std::index_sequence<Is...>)
            {
                return (init + ... + (mul_coeff * coefficient_eval(std::get<Is>(b), mul_coeff * l) * (k[Is] + l * init)));
            }

            template <std::size_t I, typename state_t, typename value_t, typename linear_t, typename tuple_t, typename array_t>
            constexpr state_t
            tpl_inner_product_b(tuple_t const& b, array_t const& k, state_t const& init, linear_t const& l, value_t mul_coeff)
            {
                return tpl_inner_product_b_impl(b, k, init, l, mul_coeff, std::make_index_sequence<I>());
            }
        } // namespace detail

        template <typename Tableau>
        struct explicit_exp_rk_butcher
        {
            Tableau butcher;
            static constexpr std::size_t N_stages = Tableau::N_stages;
            static constexpr bool is_embedded     = is_embedded_tableau<Tableau>;
            static constexpr std::size_t order    = Tableau::order;
            static constexpr char const* id       = Tableau::id;

            explicit_exp_rk_butcher(double tol_ = ponio::default_config::tol)
                : butcher()
                , tol(tol_)
            {
            }

            template <typename problem_t, typename state_t, typename value_t, typename arrayKi_t, std::size_t i>
            inline state_t stage(Stage<i>, problem_t& pb, value_t tn, state_t const& un, arrayKi_t const& Ki, value_t dt)
            {
                return pb.n(tn + butcher.c[i] * dt, detail::tpl_inner_product<i>(butcher.a, Ki, un, pb.l, dt));
            }

            template <typename problem_t, typename state_t, typename value_t, typename arrayKi_t>
            inline state_t stage(Stage<N_stages>, problem_t& pb, value_t tn, state_t const& un, arrayKi_t const& Ki, value_t dt)
            {
                return detail::tpl_inner_product_b<N_stages>(butcher.b, Ki, un, pb.l, dt);
            }

            template <typename problem_t, typename state_t, typename value_t, typename arrayKi_t, typename tab_t = Tableau>
                requires std::same_as<tab_t, Tableau> && is_embedded
            inline state_t stage(Stage<N_stages + 1>, problem_t& pb, value_t tn, state_t const& un, arrayKi_t const& Ki, value_t dt)
            {
                return detail::tpl_inner_product_b<N_stages>(butcher.b2, Ki, un, pb.l, dt);
            }

            double tol;
        };

    } // namespace exp_runge_kutta

    namespace chebyshev
    {

        /**
         * @brief Chebyshev polynomial of first kind by recursive method: \f$$\begin{aligned}T_0(x) &= 1\\ T_1(x) &= x\\ T_{n+1}(x) &=
         * 2xT_n(x) - T_{n-1}(x)\end{aligned}}\f$$
         *
         * @tparam N degree of polynomial
         * @tparam value_t type of \f$x\f$
         * @param x value where evaluate \f$T_n\f$
         */
        template <std::size_t N, typename value_t>
        value_t T(const value_t x)
        {
            if constexpr (N == 0)
            {
                return static_cast<value_t>(1.);
            }
            else if constexpr (N == 1)
            {
                return x;
            }
            else
            {
                return 2 * x * T<N - 1>(x) - T<N - 2>(x);
            }
        }

        /**
         * @brief Chebyshev polynomial of second kind by recursive method: \f$$\begin{aligned}U_0(x) &= 1\\ U_1(x) &= 2x\\ U_{n+1}(x) &=
         * 2xU_n(x) - U_{n-1}(x)\end{aligned}}\f$$
         *
         * @tparam N degree of polynomial
         * @tparam value_t type of \f$x\f$
         * @param x value where evaluate \f$U_n\f$
         */
        template <std::size_t N, typename value_t>
        value_t U(const value_t x)
        {
            if constexpr (N == 0)
            {
                return static_cast<value_t>(1.);
            }
            else if constexpr (N == 1)
            {
                return 2 * x;
            }
            else
            {
                return 2 * x * U<N - 1>(x) - U<N - 2>(x);
            }
        }

        /**
         * @brief derivatives of Chebyshev polynomial: \f$\frac{\mathrm{d}T_n}{\mathrm{d}x}(x)=nU_{n-1}(x)\f$
         *
         * @tparam N degree of polynomial
         * @tparam value_t type of \f$x\f$
         * @param x value where evaluate \f$\frac{\mathrm{d}T_n}{\mathrm{d}x}\f$
         */
        template <std::size_t N, typename value_t>
        value_t dT(const value_t x)
        {
            if constexpr (N == 0)
            {
                return static_cast<value_t>(0.);
            }
            else
            {
                return N * U<N - 1>(x);
            }
        }

        /**
         * @brief second derivative of Chebyshev polynomial: \f$\frac{\mathrm{d}^2T_n}{\mathrm{d}x^2}(x)=n\frac{nT_n(x) - xU_{n-1}(x)}{x^2
         * -1}\f$
         *
         * @tparam N degree of polynomial
         * @tparam value_t type of \f$x\f$
         * @param x value where evaluate \f$\frac{\mathrm{d}^2T_n}{\mathrm{d}x^2}\f$
         */
        template <std::size_t N, typename value_t>
        value_t ddT(const value_t x)
        {
            if constexpr (N == 0)
            {
                return static_cast<value_t>(0.);
            }
            else
            {
                return N * (N * T<N>(x) - x * U<N - 1>(x)) / (x * x - static_cast<value_t>(1.));
            }
        }

        /** @class explicit_rkc2
         *  @brief define RKC2 with user defined number of stages
         *
         *  @tparam _Nstages number of stages
         *  @tparam value_t type of coefficients
         */
        template <std::size_t _N_stages, typename value_t = double>
        struct explicit_rkc2
        {
            static_assert(_N_stages > 1, "Number of stages should be at least 2 in eRKC2");
            static constexpr bool is_embedded     = false;
            static constexpr std::size_t N_stages = _N_stages;
            value_t w0;
            value_t w1;

            template <std::size_t J>
            static constexpr value_t b(const value_t x)
            {
                if constexpr (J == 0 || J == 1)
                {
                    return b<2>(x);
                }
                else
                {
                    return ddT<J>(x) / (::detail::power<2>(dT<J>(x)));
                }
            }

            explicit_rkc2(value_t eps = 2. / 13.)
                : w0(1. + eps / (N_stages * N_stages))
            {
                w1 = dT<N_stages>(w0) / ddT<N_stages>(w0);
            }

            template <typename Problem_t, typename state_t, typename ArrayKi_t, std::size_t j>
            inline state_t stage(Stage<j>, Problem_t& f, value_t tn, state_t const& un, ArrayKi_t const& Yi, value_t dt)
            {
                value_t mj   = 2. * b<j>(w0) / b<j - 1>(w0) * w0;
                value_t nj   = -b<j>(w0) / b<j - 2>(w0);
                value_t mjt  = 2. * b<j>(w0) / b<j - 1>(w0) * w1;
                value_t gjt  = -(1. - b<j - 1>(w0) * T<j - 1>(w0)) * mjt;
                value_t cjm1 = dT<N_stages>(w0) / ddT<N_stages>(w0) * ddT<j - 1>(w0) / dT<j - 1>(w0);

                return (1. - mj - nj) * un + mj * Yi[j - 1] + nj * Yi[j - 2] + mjt * dt * f(tn + cjm1 * dt, Yi[j - 1]) + gjt * dt * Yi[0];
            }

            template <typename Problem_t, typename state_t, typename ArrayKi_t>
            inline state_t stage(Stage<0>, Problem_t& f, value_t tn, state_t const& un, ArrayKi_t const& Yi, value_t dt)
            {
                return f(tn, un); // be careful Yi[0] stores f(tn,un) not un!!!
            }

            template <typename Problem_t, typename state_t, typename ArrayKi_t>
            inline state_t stage(Stage<1>, Problem_t& f, value_t tn, state_t const& un, ArrayKi_t const& Yi, value_t dt)
            {
                value_t m1t = b<1>(w0) * w1;
                return un + dt * m1t * Yi[0];
            }

            template <typename Problem_t, typename state_t, typename ArrayKi_t>
            inline state_t stage(Stage<2>, Problem_t& f, value_t tn, state_t const& un, ArrayKi_t const& Yi, value_t dt)
            {
                value_t m2  = 2. * w0;
                value_t n2  = -1.;
                value_t m2t = 2. * w1;
                value_t c2  = dT<N_stages>(w0) / ddT<N_stages>(w0) * ddT<2>(w0) / dT<2>(w0);
                value_t c1  = c2 / dT<2>(w0);
                value_t g2t = -(1. - b<1>(w0) * T<1>(w0)) * m2t;

                return (1. - m2 - n2) * un + m2 * Yi[1] + n2 * un + m2t * dt * f(tn + c1 * dt, Yi[1]) + g2t * dt * Yi[0];
            }
        };

    } // namespace chebyshev

} // namespace ode::butcher
