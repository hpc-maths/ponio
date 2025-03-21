// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// NOLINTBEGIN(misc-include-cleaner)

#include <cstddef>
#include <filesystem>
#include <initializer_list>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include <ponio/eigen_linear_algebra.hpp>
#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

// NOLINTEND(misc-include-cleaner)

namespace detail
{
    template <class tuple_t, class func_t, std::size_t... I>
    constexpr func_t
    for_each_impl( tuple_t const& t, func_t&& f, std::index_sequence<I...> )
    {
        return (void)std::initializer_list<int>{ ( std::forward<func_t>( f )( std::get<I>( t ) ), 0 )... }, f;
    }
}

template <class tuple_t, class func_t>
constexpr func_t
for_each( tuple_t const& t, func_t&& f )
{
    return detail::for_each_impl( t,
        std::forward<func_t>( f ),
        std::make_index_sequence<std::tuple_size<std::remove_reference_t<tuple_t>>::value>{} );
}

int
main( int, char** )
{
    std::string const dirname = "lorenz_all_data";

    using vector_type = Eigen::Vector<double, 3>;
    using matrix_type = Eigen::Matrix<double, 3, 3>;
    using state_t     = vector_type;

    double const sigma = 10.;
    double const rho   = 28.;
    double const beta  = 8. / 3.;

    state_t const u0 = { 1., 1., 1. };

    ponio::time_span<double> const t_span = { 0., 15. };
    double const dt                       = 0.01;

    // explicit methods
    {
        // define problem
        auto lorenz = [=]( double, state_t const& u ) -> state_t
        {
            auto du1 = sigma * ( u[1] - u[0] );
            auto du2 = rho * u[0] - u[1] - u[0] * u[2];
            auto du3 = u[0] * u[1] - beta * u[2];
            return { du1, du2, du3 };
        };

        auto explicit_solve_lorenz = [&]( auto& algo )
        {
            std::cout << algo.id << "...";

            std::stringstream basename;
            basename << algo.id << ".dat";
            auto filename = std::filesystem::path( dirname ) / basename.str();

            ponio::solve( lorenz, algo, u0, t_span, dt, ponio::observer::file_observer( filename ) );

            std::cout << "\n";
        };

        auto erk_algo = ponio::runge_kutta::erk_tuple<double>();
        for_each( erk_algo, explicit_solve_lorenz );
    }

    // lawson methods
    {
        // define problem
        auto L = matrix_type{
            {-sigma, sigma, 0    },
            { rho,   -1,    0    },
            { 0,     0,     -beta}
        };
        auto N = [=]( double, vector_type const& u ) -> vector_type
        {
            return { 0., -u[0] * u[2], u[0] * u[1] };
        };
        auto lorenz = ponio::make_lawson_problem( L, N );

        // define matrix exponential
        auto mexp = []( matrix_type const& A ) -> matrix_type
        {
            return A.exp();
        };

        auto lawson_solve_lorenz = [&]( auto& algo )
        {
            using algo_t = decltype( algo( mexp ) );
            std::cout << "l" << algo_t::id << "...";

            std::stringstream basename;
            basename << "l" << algo_t::id << ".dat";
            auto filename = std::filesystem::path( dirname ) / basename.str();

            ponio::solve( lorenz, algo( mexp ), u0, t_span, dt, ponio::observer::file_observer( filename ) );

            std::cout << "\n";
        };

        auto lrk_algo = ponio::runge_kutta::lrk_tuple<double, decltype( mexp )>();
        for_each( lrk_algo, lawson_solve_lorenz );
    }

    // dirk methods
    {
        // define problem
        auto f = [=]( double, state_t const& u ) -> state_t
        {
            auto du1 = sigma * ( u[1] - u[0] );
            auto du2 = rho * u[0] - u[1] - u[0] * u[2];
            auto du3 = u[0] * u[1] - beta * u[2];
            return { du1, du2, du3 };
        };

        auto jac_f = [=]( double, state_t const& u ) -> matrix_type
        {
            return matrix_type( {
                {-sigma,      sigma, 0    },
                { rho - u[2], -1,    -u[0]},
                { u[1],       u[0],  -beta}
            } );
        };

        auto lorenz = ponio::make_implicit_problem( f, jac_f );

        auto dirk_solve_lorenz = [&]( auto& algo )
        {
            std::cout << algo.id << "...";

            std::stringstream basename;
            basename << algo.id << ".dat";
            auto filename = std::filesystem::path( dirname ) / basename.str();

            ponio::solve( lorenz, algo, u0, t_span, dt, ponio::observer::file_observer( filename ) );

            std::cout << "\n";
        };

        auto dirk_algo = std::make_tuple( ponio::runge_kutta::backward_euler(),
            ponio::runge_kutta::crancknicolson(),
            ponio::runge_kutta::crancknicolson_2(),
            ponio::runge_kutta::dirk23(),
            ponio::runge_kutta::dirk23_crouzeix(),
            ponio::runge_kutta::dirk34(),
            ponio::runge_kutta::dirk_qin_zhang(),
            ponio::runge_kutta::implicit_midpoint(),
            ponio::runge_kutta::lsdirk22qz(),
            ponio::runge_kutta::lsdirk33(),
            ponio::runge_kutta::lsdirk43(),
            ponio::runge_kutta::sspirk33() );
        for_each( dirk_algo, dirk_solve_lorenz );
    }

    // splitting methods
    {
        auto linear_part = [=]( double, vector_type const& u ) -> vector_type
        {
            return { -sigma * u[0] + sigma * u[1], rho * u[0] - u[1], -beta * u[2] };
        };
        auto nonlinear_part = [=]( double, vector_type const& u ) -> vector_type
        {
            return { 0., -u[0] * u[2], u[0] * u[1] };
        };
        auto lorenz = ponio::make_problem( linear_part, nonlinear_part );

        auto splitting_solve_lorenz = [&]( auto& algo )
        {
            std::cout << algo.id << "...";

            std::stringstream basename;
            basename << algo.id << ".dat";
            auto filename = std::filesystem::path( dirname ) / basename.str();

            ponio::solve( lorenz, algo, u0, t_span, dt, ponio::observer::file_observer( filename ) );

            std::cout << "\n";
        };

        auto splitting_algo = std::make_tuple( ponio::splitting::make_lie_tuple( std::make_pair( ponio::runge_kutta::rk_44(), 0.005 ),
                                                   std::make_pair( ponio::runge_kutta::rk_44(), 0.005 ) ),
            ponio::splitting::make_strang_tuple( std::make_pair( ponio::runge_kutta::rk_44(), 0.005 ),
                std::make_pair( ponio::runge_kutta::rk_44(), 0.005 ) ),
            ponio::splitting::make_adaptive_strang_tuple( 5e-3,
                1e-3,
                std::make_pair( ponio::runge_kutta::rk_44(), 0.005 ),
                std::make_pair( ponio::runge_kutta::rk_44(), 0.005 ) ) );
        for_each( splitting_algo, splitting_solve_lorenz );
    }

    return 0;
}
