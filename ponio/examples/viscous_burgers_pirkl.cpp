// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// NOLINTBEGIN(misc-include-cleaner)

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <numbers>
#include <valarray>
#include <vector>

#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta/pirkl.hpp>
#include <ponio/solver.hpp>

// NOLINTEND(misc-include-cleaner)

/*
Solve the one-dimensional viscous Burgers equation

    u_t + (u^2 / 2)_x = nu u_xx

with periodic boundary conditions. Diffusion is discretized with centered
finite differences and advection with a fifth-order WENO flux.
*/

namespace
{
    using state_t = std::valarray<double>;

    double
    burgers_flux( double u )
    {
        return 0.5 * u * u;
    }

    double
    weno5_flux( std::array<double, 6> const& u )
    {
        constexpr double eps = 1e-14;

        constexpr double ceno[3][4] = {
            {1. / 3.,   -1. / 6., 1. / 3.,  11. / 6.},
            { -7. / 6., 5. / 6.,  5. / 6.,  -7. / 6.},
            { 11. / 6., 1. / 3.,  -1. / 6., 1. / 3. }
        };
        constexpr double ckrp[3] = { 1. / 10., 6. / 10., 3. / 10. };
        constexpr double ckrm[3] = { 3. / 10., 6. / 10., 1. / 10. };

        double f[6];
        for ( std::size_t i = 0; i < 6; ++i )
        {
            f[i] = burgers_flux( u[i] );
        }

        double g[3][6] = {};
        for ( int i = -2; i <= 3; ++i )
        {
            g[0][i + 2] = f[i + 2];
        }
        for ( int k = 1; k <= 2; ++k )
        {
            for ( int i = -2; i <= 3 - k; ++i )
            {
                g[k][i + 2] = g[k - 1][i + 3] - g[k - 1][i + 2];
            }
        }

        double smoothness[3] = {};
        double alpha[3]      = {};
        double omega[3]      = {};

        bool const positive_velocity = 0.5 * ( u[2] + u[3] ) > 0.;
        auto const& weights          = positive_velocity ? ckrp : ckrm;

        for ( int stencil = 0; stencil < 3; ++stencil )
        {
            for ( int k = 1; k <= 2; ++k )
            {
                for ( int j = 1; j <= 3 - k; ++j )
                {
                    int const index = positive_velocity ? stencil + j - 3 : stencil + j - 2;
                    double const v  = g[k][index + 2];
                    smoothness[stencil] += v * v / static_cast<double>( 3 - k );
                }
            }
            alpha[stencil] = weights[stencil] / std::pow( eps + smoothness[stencil], 2 );
        }

        double const alpha_sum = alpha[0] + alpha[1] + alpha[2];
        for ( int stencil = 0; stencil < 3; ++stencil )
        {
            omega[stencil] = alpha[stencil] / alpha_sum;
        }

        double numerical_flux = 0.;
        for ( int stencil = 0; stencil < 3; ++stencil )
        {
            for ( int k = 0; k < 3; ++k )
            {
                int const index  = positive_velocity ? -2 + k + stencil : -1 + k + stencil;
                int const column = positive_velocity ? stencil : stencil + 1;
                numerical_flux += omega[stencil] * ceno[k][column] * g[0][index + 2];
            }
        }

        return numerical_flux;
    }

    class burgers_diffusion
    {
        double m_nu;
        double m_inv_dx2;

      public:

        burgers_diffusion( double nu, double dx )
            : m_nu( nu )
            , m_inv_dx2( 1. / ( dx * dx ) )
        {
        }

        void
        operator()( double, state_t const& u, state_t& du ) const
        {
            std::size_t const n = u.size();

            du[0] = m_nu * m_inv_dx2 * ( u[n - 1] - 2. * u[0] + u[1] );
            for ( std::size_t i = 1; i + 1 < n; ++i )
            {
                du[i] = m_nu * m_inv_dx2 * ( u[i - 1] - 2. * u[i] + u[i + 1] );
            }
            du[n - 1] = m_nu * m_inv_dx2 * ( u[n - 2] - 2. * u[n - 1] + u[0] );
        }

        double
        spectral_radius() const
        {
            return 4. * m_nu * m_inv_dx2;
        }
    };

    class burgers_advection
    {
        double m_inv_dx;

      public:

        explicit burgers_advection( double dx )
            : m_inv_dx( 1. / dx )
        {
        }

        void
        operator()( double, state_t const& u, state_t& du ) const
        {
            constexpr std::size_t ng = 3;
            std::size_t const n      = u.size();

            std::vector<double> ghosted( n + 2 * ng );
            for ( std::size_t i = 0; i < n; ++i )
            {
                ghosted[ng + i] = u[i];
            }
            for ( std::size_t i = 0; i < ng; ++i )
            {
                ghosted[i]          = u[n - ng + i];
                ghosted[ng + n + i] = u[i];
            }

            std::vector<double> flux( n + 1 );
            for ( std::size_t interface = 0; interface <= n; ++interface )
            {
                std::array<double, 6> stencil{};
                for ( int offset = -2; offset <= 3; ++offset )
                {
                    auto const index = static_cast<std::size_t>(
                        static_cast<long long>( ng ) - 1 + static_cast<long long>( interface ) + offset );
                    stencil[static_cast<std::size_t>( offset + 2 )] = ghosted[index];
                }
                flux[interface] = weno5_flux( stencil );
            }

            for ( std::size_t i = 0; i < n; ++i )
            {
                du[i] = -m_inv_dx * ( flux[i + 1] - flux[i] );
            }
        }
    };
}

int
main()
{
    constexpr std::size_t nx = 200;
    constexpr double xmin    = -1.;
    constexpr double xmax    = 1.;
    constexpr double nu      = 2e-2;
    constexpr double t_ini   = 0.;
    constexpr double t_end   = 5e-2;

    double const dx = ( xmax - xmin ) / static_cast<double>( nx );
    double const dt = 2.5e-3;

    state_t u_ini( nx );
    for ( std::size_t i = 0; i < nx; ++i )
    {
        double const x = xmin + ( static_cast<double>( i ) + 0.5 ) * dx;
        u_ini[i]       = std::sin( std::numbers::pi * x );
    }

    burgers_diffusion diffusion( nu, dx );
    burgers_advection advection( dx );
    auto problem = ponio::make_problem( diffusion, advection );

    auto spectral_radius = [&diffusion]( auto&, double, auto const&, double, auto& )
    {
        return diffusion.spectral_radius();
    };

    auto pirkl          = ponio::runge_kutta::pirkl::pirkl_DA( spectral_radius );
    state_t const u_end = ponio::solve( problem, pirkl, u_ini, { t_ini, t_end }, dt, ponio::observer::null_observer() );

    std::cout << "u(0, " << t_end << ") = " << u_end[nx / 2] << '\n';

    return 0;
}
