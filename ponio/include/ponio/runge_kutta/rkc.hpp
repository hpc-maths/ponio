
// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string_view> // NOLINT(misc-include-cleaner)
#include <utility>

#include "../detail.hpp" // NOLINT(misc-include-cleaner)
#include "../iteration_info.hpp"
#include "../stage.hpp" // NOLINT(misc-include-cleaner)

namespace ponio::runge_kutta::chebyshev
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
    value_t
    T( value_t const x )
    {
        if constexpr ( N == 0 )
        {
            return static_cast<value_t>( 1. );
        }
        else if constexpr ( N == 1 )
        {
            return x;
        }
        else
        {
            return 2 * x * T<N - 1>( x ) - T<N - 2>( x );
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
    value_t
    U( value_t const x )
    {
        if constexpr ( N == 0 )
        {
            return static_cast<value_t>( 1. );
        }
        else if constexpr ( N == 1 )
        {
            return 2 * x;
        }
        else
        {
            return 2 * x * U<N - 1>( x ) - U<N - 2>( x );
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
    value_t
    dT( value_t const x )
    {
        if constexpr ( N == 0 )
        {
            return static_cast<value_t>( 0. );
        }
        else
        {
            return N * U<N - 1>( x );
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
    value_t
    ddT( value_t const x )
    {
        if constexpr ( N == 0 )
        {
            return static_cast<value_t>( 0. );
        }
        else
        {
            return N * ( N * T<N>( x ) - x * U<N - 1>( x ) ) / ( x * x - static_cast<value_t>( 1. ) );
        }
    }

    /** @class explicit_rkc2
     *  @brief define RKC2 with user defined number of stages
     *
     *  @tparam N_stages_ number of stages
     *  @tparam value_t type of coefficients
     */
    template <std::size_t N_stages_, typename _value_t = double>
    struct explicit_rkc2
    {
        static_assert( N_stages_ > 1, "Number of stages should be at least 2 in eRKC2" );
        static constexpr std::size_t N_stages = N_stages_;
        static constexpr std::size_t order    = 2;
        static constexpr std::string_view id  = "RKC2";
        static constexpr bool is_embedded     = false;
        using value_t                         = _value_t;

        value_t w0;
        value_t w1;
        iteration_info<explicit_rkc2> _info;

        /**
         * @brief computes \f$b_j = \f$ coefficients with following formula: \f$b_0=b_2\f$, \f$b_1=\frac{1}{\omega_0}\f$, \f$b_j =
         * \frac{T_j''(\omega_0)}{(T_j'(\omega_0))^2}\f$
         *
         * @tparam J index \f$j\f$
         * @param x  value of \f$\omega_0\f$
         */
        template <std::size_t J>
        static constexpr value_t
        b( value_t const x )
        {
            if constexpr ( J == 0 || J == 1 )
            {
                return b<2>( x );
            }
            else
            {
                return ddT<J>( x ) / ( detail::power<2>( dT<J>( x ) ) );
            }
        }

        /**
         * @brief Construct a new explicit RKC2 algorithm
         *
         * @param eps value of relaxation
         */
        explicit_rkc2( value_t eps = 2. / 13. )
            : w0( 1. + eps / ( N_stages * N_stages ) )
            , w1( dT<N_stages>( w0 ) / ddT<N_stages>( w0 ) )
            , _info()
        {
            _info.number_of_eval = N_stages;
        }

        /**
         * @brief computes stage \f$j\f$ of RKC2 method
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
         * @details \f$y_j = (1 - \mu_j - \nu_j)y^n + \mu_j y_{j-1} + \nu_j y_{j-2} + \tilde{\mu}_j\Delta tf(t^n + c_{j-1}\Delta t, y_{j-1})
         * + \tilde{\gamma}_j\Delta t f(t^n, y^n)\f$
         */
        template <typename problem_t, typename state_t, typename array_ki_t, std::size_t j>
        void
        stage( Stage<j>, problem_t& f, value_t tn, state_t& yn, array_ki_t const& Yj, value_t dt, state_t& ui, state_t& yi )
        {
            value_t const mj   = 2. * b<j>( w0 ) / b<j - 1>( w0 ) * w0;
            value_t const nj   = -b<j>( w0 ) / b<j - 2>( w0 );
            value_t const mjt  = 2. * b<j>( w0 ) / b<j - 1>( w0 ) * w1;
            value_t const gjt  = -( 1. - b<j - 1>( w0 ) * T<j - 1>( w0 ) ) * mjt;
            value_t const cjm1 = dT<N_stages>( w0 ) / ddT<N_stages>( w0 ) * ddT<j - 1>( w0 ) / dT<j - 1>( w0 );

            f( tn + cjm1 * dt, Yj[j - 1], ui ); // first compute ui = f(t^n + c_{j-1}\Delta t, y_{j-1})
            yi = ( 1. - mj - nj ) * yn + mj * Yj[j - 1] + nj * Yj[j - 2] + mjt * dt * ui + gjt * dt * Yj[0];
        }

        /**
         * @brief computes pseudo first stage of RKC2 method
         *
         * @tparam problem_t  type of operator \f$f\f$
         * @tparam state_t    type of current state
         * @tparam array_ki_t type of array of temporary stages
         * @tparam j          integer of stage
         * @param f  operator \f$f\f$
         * @param tn current time
         * @param yn current state
         * @param yi computed output step
         *
         * @details \f$y_0 = f(t^n, y^n)\f$
         */
        template <typename problem_t, typename state_t, typename array_ki_t>
        void
        stage( Stage<0>, problem_t& f, value_t tn, state_t& yn, array_ki_t const&, value_t, state_t&, state_t& yi )
        {
            f( tn, yn, yi ); // be careful Yj[0] stores f(tn,yn) not yn!!!
        }

        /**
         * @brief computes first stage of RKC2 method
         *
         * @tparam problem_t  type of operator \f$f\f$
         * @tparam state_t    type of current state
         * @tparam array_ki_t type of array of temporary stages
         * @tparam j          integer of stage
         * @param yn current state
         * @param Yj array of temporary stages
         * @param dt current time step
         * @param yi computed output step
         *
         * @details \f$y_1 = y^n + \tilde{\mu}_1 \Delta t f(t^n, y^n)\f$
         */
        template <typename problem_t, typename state_t, typename array_ki_t>
        void
        stage( Stage<1>, problem_t&, value_t, state_t& yn, array_ki_t const& Yj, value_t dt, state_t&, state_t& yi )
        {
            value_t const m1t = b<1>( w0 ) * w1;
            yi                = yn + dt * m1t * Yj[0];
        }

        /**
         * @brief computes second stage of RKC2 method
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
         * @details \f$y_2 = (1-\mu_2 -\nu_2)y^n + \mu_2y_1 + \nu_2y^n + \tilde{\mu}_2\Delta t f(t^n + c_1\Delta t, y_1) + \tilde{\gamma}_2
         * \Delta t f(t^n, y^n)\f$
         */
        template <typename problem_t, typename state_t, typename array_ki_t>
        void
        stage( Stage<2>, problem_t& f, value_t tn, state_t& yn, array_ki_t const& Yj, value_t dt, state_t& ui, state_t& yi )
        {
            value_t const m2  = 2. * w0;
            value_t const n2  = -1.;
            value_t const m2t = 2. * w1;
            value_t const c2  = dT<N_stages>( w0 ) / ddT<N_stages>( w0 ) * ddT<2>( w0 ) / dT<2>( w0 );
            value_t const c1  = c2 / dT<2>( w0 );
            value_t const g2t = -( 1. - b<1>( w0 ) * T<1>( w0 ) ) * m2t;

            f( tn + c1 * dt, Yj[1], ui ); // first compute yi = f(t^n + c_1\Delta t, y_1)
            yi = ( 1. - m2 - n2 ) * yn + m2 * Yj[1] + n2 * yn + m2t * dt * ui + g2t * dt * Yj[0];
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

    namespace rkc_detail
    {
        /**
         * @brief RKC1 damping/stability constant.
         *
         * For the damped first-order RKC polynomial, beta = 2 - 4 eps / 3
         * gives the sufficient real-axis stability interval beta * m^2.
         */
        template <typename value_t>
        constexpr value_t
        beta( value_t eps )
        {
            return static_cast<value_t>( 2. ) - static_cast<value_t>( 4. ) * eps / static_cast<value_t>( 3. );
        }

        /**
         * @brief Compute T_m(x) and T'_m(x) by forward recurrence.
         */
        template <typename value_t>
        void
        chebyshev_T_dT( std::size_t m, value_t x, value_t& Tm, value_t& dTm )
        {
            if ( m == 0 )
            {
                Tm  = static_cast<value_t>( 1. );
                dTm = static_cast<value_t>( 0. );
                return;
            }

            if ( m == 1 )
            {
                Tm  = x;
                dTm = static_cast<value_t>( 1. );
                return;
            }

            value_t Tjm2  = static_cast<value_t>( 1. );
            value_t Tjm1  = x;
            value_t dTjm2 = static_cast<value_t>( 0. );
            value_t dTjm1 = static_cast<value_t>( 1. );

            for ( std::size_t j = 2; j <= m; ++j )
            {
                value_t const Tj  = static_cast<value_t>( 2. ) * x * Tjm1 - Tjm2;
                value_t const dTj = static_cast<value_t>( 2. ) * Tjm1 + static_cast<value_t>( 2. ) * x * dTjm1 - dTjm2;

                Tjm2  = Tjm1;
                Tjm1  = Tj;
                dTjm2 = dTjm1;
                dTjm1 = dTj;
            }

            Tm  = Tjm1;
            dTm = dTjm1;
        }

        template <typename value_t>
        value_t
        omega_0( std::size_t m, value_t eps )
        {
            value_t const mm = static_cast<value_t>( m ) * static_cast<value_t>( m );
            return static_cast<value_t>( 1. ) + eps / mm;
        }

        template <typename value_t>
        value_t
        omega_1( std::size_t m, value_t eps )
        {
            value_t const w0 = omega_0( m, eps );

            value_t Tm  = static_cast<value_t>( 0. );
            value_t dTm = static_cast<value_t>( 0. );
            chebyshev_T_dT( m, w0, Tm, dTm );

            return Tm / dTm;
        }

        template <typename value_t>
        value_t
        stability_interval( std::size_t m, value_t eps )
        {
            value_t const w0 = omega_0( m, eps );
            value_t const w1 = omega_1( m, eps );
            return static_cast<value_t>( 2. ) * w0 / w1;
        }

        /**
         * @brief Smallest stage number for standalone dynamic RKC1.
         *
         * The beta*m^2 bound gives an initial estimate; the exact real-axis
         * interval 2*omega_0/omega_1 is then used to return the smallest
         * admissible stage number.
         */
        template <typename value_t>
        std::size_t
        dynamic_rkc1_stages( value_t dt, value_t rho, value_t eps )
        {
            value_t const z = std::abs( dt ) * rho;

            if ( z == static_cast<value_t>( 0. ) )
            {
                return 1;
            }

            value_t const beta_value = beta( eps );
            if ( beta_value <= static_cast<value_t>( 0. ) )
            {
                throw std::runtime_error( "RKC1: damping parameter gives a non-positive beta." );
            }

            std::size_t m = std::max<std::size_t>( 1, static_cast<std::size_t>( std::ceil( std::sqrt( z / beta_value ) ) ) );

            while ( m > 1 && z <= stability_interval( m - 1, eps ) )
            {
                --m;
            }

            while ( z > stability_interval( m, eps ) )
            {
                ++m;
            }

            return m;
        }

        /**
         * @brief One RKC1 step with a prescribed number of stages.
         *
         * This is the common RKC1 kernel used by standalone dynamic RKC1,
         * by the mRKC micro-solver, and by the mRKC macro-solver.
         *
         * Stage abscissae are propagated with the same RKC recurrence, so the
         * kernel also supports non-autonomous right-hand sides.
         */
        template <typename problem_t, typename value_t, typename state_t>
        void
        rkc1_step( problem_t&& f,
            value_t t,
            state_t& y,
            value_t dt,
            std::size_t n_stages,
            value_t eps,
            state_t& kjm2,
            state_t& kjm1,
            state_t& kj,
            state_t& f_tmp,
            state_t& y_out )
        {
            if ( n_stages == 0 )
            {
                throw std::invalid_argument( "RKC1: number of stages must be at least one." );
            }

            value_t const w0 = omega_0( n_stages, eps );
            value_t const w1 = omega_1( n_stages, eps );

            value_t Tjm2 = static_cast<value_t>( 1. );
            value_t Tjm1 = w0;

            kjm2 = y;

            value_t const mu1 = w1 / w0;

            value_t cjm2 = static_cast<value_t>( 0. );
            value_t cjm1 = mu1;

            std::forward<problem_t>( f )( t, kjm2, f_tmp );
            kjm1 = kjm2 + mu1 * dt * f_tmp;

            if ( n_stages == 1 )
            {
                y_out = kjm1;
                return;
            }

            for ( std::size_t j = 2; j <= n_stages; ++j )
            {
                value_t const Tj = static_cast<value_t>( 2. ) * w0 * Tjm1 - Tjm2;

                value_t const bjm2 = static_cast<value_t>( 1. ) / Tjm2;
                value_t const bjm1 = static_cast<value_t>( 1. ) / Tjm1;
                value_t const bj   = static_cast<value_t>( 1. ) / Tj;

                value_t const mu_j    = static_cast<value_t>( 2. ) * w1 * bj / bjm1;
                value_t const nu_j    = static_cast<value_t>( 2. ) * w0 * bj / bjm1;
                value_t const kappa_j = -bj / bjm2;

                std::forward<problem_t>( f )( t + cjm1 * dt, kjm1, f_tmp );
                kj = nu_j * kjm1 + kappa_j * kjm2 + mu_j * dt * f_tmp;

                value_t const cj = nu_j * cjm1 + kappa_j * cjm2 + mu_j;

                if ( j < n_stages )
                {
                    std::swap( kjm2, kjm1 );
                    std::swap( kjm1, kj );

                    cjm2 = cjm1;
                    cjm1 = cj;
                }

                Tjm2 = Tjm1;
                Tjm1 = Tj;
            }

            y_out = kj;
        }
    } // namespace rkc_detail

    template <typename value_t>
    inline constexpr value_t rkc1_default_eps = static_cast<value_t>( 0.05 );

    /**
     * @brief First-order explicit Runge--Kutta--Chebyshev method (RKC1)
     *        with a number of stages selected dynamically at run time.
     *
     * The spectral radius is supplied by the caller through a callable, using
     * the same interface as ROCK2. No power iteration is performed internally
     * unless the caller explicitly provides such a spectral-radius estimator.
     */
    template <typename eig_computer_t, typename _value_t = double>
    struct explicit_rkc1
    {
        static constexpr bool is_embedded      = false;
        static constexpr std::size_t N_stages  = stages::dynamic;
        static constexpr std::size_t N_storage = 4;
        static constexpr std::size_t order     = 1;
        static constexpr std::string_view id   = "RKC1";

        using value_t = _value_t;

        static constexpr value_t default_eps = rkc1_default_eps<value_t>;

        iteration_info<explicit_rkc1> _info;

        eig_computer_t eig_computer;
        value_t eps;

        explicit explicit_rkc1( eig_computer_t&& _eig_computer, value_t damping = default_eps )
            : _info()
            , eig_computer( std::forward<eig_computer_t>( _eig_computer ) )
            , eps( damping )
        {
            if ( eps < static_cast<value_t>( 0. ) || rkc_detail::beta( eps ) <= static_cast<value_t>( 0. ) )
            {
                throw std::invalid_argument( "RKC1: damping parameter eps must be non-negative and satisfy 2 - 4 eps / 3 > 0." );
            }
        }

        std::size_t
        required_stages( value_t dt, value_t spectral_radius ) const
        {
            if ( spectral_radius < static_cast<value_t>( 0. ) )
            {
                throw std::invalid_argument( "RKC1: spectral radius must be non-negative." );
            }
            return rkc_detail::dynamic_rkc1_stages( dt, spectral_radius, eps );
        }

        template <typename problem_t, typename state_t, typename array_ki_t>
        void
        operator()( problem_t& f, value_t& tn, state_t& un, array_ki_t& K, value_t& dt, state_t& unp1 )
        {
            _info.reset_eval();

            std::size_t n_eval = 0;
            auto f_counter     = [&n_eval, &f]( value_t t, state_t& u, state_t& du )
            {
                ++n_eval;
                f( t, u, du );
            };

            value_t const spectral_radius = eig_computer( f_counter, tn, un, dt, K );
            std::size_t const m           = required_stages( dt, spectral_radius );

            _info.number_of_stages = m;
            _info.number_of_eval   = n_eval + m;

            rkc_detail::rkc1_step( f, tn, un, dt, m, eps, K[0], K[1], K[2], K[3], unp1 );

            tn += dt;
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

    template <typename value_t = double, typename eig_computer_t>
    auto
    rkc1( eig_computer_t&& eig_computer, value_t eps = rkc1_default_eps<value_t> )
    {
        return explicit_rkc1<eig_computer_t, value_t>( std::forward<eig_computer_t>( eig_computer ), eps );
    }

} // namespace ponio::runge_kutta::chebyshev
