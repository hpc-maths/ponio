// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <cmath>
#include <numeric>
#include <tuple>
#include <valarray>
#include <vector>

#define DEBUG

#ifdef DEBUG

#include <filesystem>
#include <sstream>
namespace fs = std::filesystem;
#include <ponio/observer.hpp>
#include <string>
#include <string_view>
using namespace std::string_literals;

#endif

#include <ponio/detail.hpp>
#include <ponio/problem.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

template <typename Container>
auto
mayor_method( Container const& x, Container const& y )
{
    using value_t = typename Container::value_type;

    auto x_mid = std::begin( x ) + ( std::end( x ) - std::begin( x ) ) / 2;
    auto y_mid = std::begin( y ) + ( std::end( y ) - std::begin( y ) ) / 2;

    value_t x1_avg = std::accumulate( x.begin(), x_mid, value_t{ 0.0 } ) / static_cast<value_t>( x_mid - x.begin() );
    value_t x2_avg = std::accumulate( x_mid, x.end(), value_t{ 0.0 } ) / static_cast<value_t>( x.end() - x_mid );
    value_t y1_avg = std::accumulate( y.begin(), y_mid, value_t{ 0.0 } ) / static_cast<value_t>( y_mid - y.begin() );
    value_t y2_avg = std::accumulate( y_mid, y.end(), value_t{ 0.0 } ) / static_cast<value_t>( y.end() - y_mid );

    value_t a = ( y2_avg - y1_avg ) / ( x2_avg - x1_avg );
    value_t b = y1_avg - a * x1_avg;

    return std::make_tuple( a, b );
}

template <typename T = double>
auto
error( T u, T v )
{
    return std::abs( u - v );
}

template <typename T = double>
auto
relative_error( T u, T v )
{
    return std::abs( ( u - v ) / u );
}

namespace explicit_method
{
    template <typename Algorithm_t, typename T = double>
    auto
    solve_exp( Algorithm_t& algo, T dt, T Tf )
    {
        using state_t = T;

        auto pb = []( T, state_t y, state_t& dy )
        {
            dy = y;
        };

        state_t y0                 = 1.0;
        ponio::time_span<T> t_span = { 0., Tf };

#ifdef DEBUG
        std::stringstream ss;
        ss << "debug_info/" << std::string( Algorithm_t::id ) << "/dt_" << dt << ".dat";
        std::string const filename = ss.str();
        auto obs                   = ponio::observer::file_observer( filename );
#else
        auto obs = []( T, state_t, T ) {};
#endif
        return ponio::solve( pb, algo, y0, t_span, dt, obs );
    }

    template <typename Algorithm_t, typename T = double>
    T
    short_time_check_order( Algorithm_t algo = Algorithm_t() )
    {
        using state_t = T;

        std::vector<T> errors;
        std::vector<T> dts;

        T Tf = 1.0;

        state_t u_exa = std::exp( Tf );

#ifdef DEBUG
        fs::create_directories( "debug_info/" + std::string( Algorithm_t::id ) );
        std::ofstream f( "debug_info/" + std::string( Algorithm_t::id ) + "/errors.dat"s );
#endif

        for ( auto n_iter : { 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30 } )
        {
            T dt          = Tf / n_iter;
            state_t u_sol = solve_exp( algo, dt, Tf );
            auto e        = error( u_exa, u_sol );
            errors.push_back( std::log( e ) );
            dts.push_back( std::log( dt ) );

#ifdef DEBUG
            f << dt << " " << e << "\n";
#endif
        }

#ifdef DEBUG
        f.close();
#endif

        auto [a, b] = mayor_method( dts, errors );

        return a;
    }

    template <typename Algorithm_t, typename T = double>
    T
    long_time_check_order( Algorithm_t& algo )
    {
        using state_t = std::valarray<T>;

        T alpha = 2. / 3.;
        T beta  = 4. / 3.;
        T gamma = 1.;
        T delta = 1.;

        auto pb = [=]( T, auto&& u, state_t& du )
        {
            du[0] = alpha * u[0] - beta * u[0] * u[1];
            du[1] = delta * u[0] * u[1] - gamma * u[1];
        };

        // invariant calculator
        auto V = [=]( state_t const& u ) -> T
        {
            return delta * u[0] - gamma * std::log( u[0] ) + beta * u[1] - alpha * std::log( u[1] );
        };

        T x0          = 1.9;
        state_t u_ini = { x0, x0 };
        T V_ini       = V( u_ini );

        ponio::time_span<T> t_span = { 0., 1000. };
        std::vector<T> dts         = { 0.25, 0.125, 0.1, 0.075, 0.05 }; // find a way to adapt this range to the method
        std::vector<T> relative_errors;
        std::vector<T> log_dts;

        for ( auto dt : dts )
        {
            state_t u_end = ::ponio::solve( pb, algo, u_ini, t_span, dt, []( T, state_t const&, T ) {} );
            relative_errors.push_back( std::log10( relative_error( V_ini, V( u_end ) ) ) );
            log_dts.push_back( std::log10( dt ) );
        }

        auto [a, b] = mayor_method( log_dts, relative_errors );

        return a;
    }

    template <typename Algorithm_t, typename T = double>
    T
    check_order( Algorithm_t algo = Algorithm_t() )
    {
        if constexpr ( ponio::runge_kutta::butcher::is_embedded<Algorithm_t> )
        {
            return Algorithm_t::order;
        }
        else if constexpr ( Algorithm_t::order >= 8 || ponio::splitting::is_splitting_method<Algorithm_t> )
        {
            return long_time_check_order( algo );
        }
        else
        {
            return short_time_check_order( algo );
        }
    }

} // namespace explicit_method

namespace additive_method
{

    template <typename Algorithm_t, typename T = double>
    auto
    solve_exp( Algorithm_t& algo, T dt, T Tf, T lambda )
    {
        using state_t = T;

        auto pb = ponio::make_imex_jacobian_problem(
            [=]( T, state_t y, state_t& dy )
            {
                dy = lambda * y;
            },
            [=]( T, state_t y, state_t& dy )
            {
                dy = ( 1. - lambda ) * y;
            },
            [=]( T, state_t )
            {
                return 1. - lambda;
            } );

        state_t y0                 = 1.0;
        ponio::time_span<T> t_span = { 0., Tf };

        auto obs = []( T, state_t, T ) {};
        return ::ponio::solve( pb, algo, y0, t_span, dt, obs );
    }

    template <typename Algorithm_t, typename T = double>
    T
    short_time_check_order( Algorithm_t algo = Algorithm_t(), T lambda = 0.33 )
    {
        using state_t = T;

        std::vector<T> errors;
        std::vector<T> dts;

        T Tf = 1.0;

        state_t u_exa = std::exp( Tf );

        for ( auto n_iter : { 50, 25, 20, 15, 10 } )
        {
            T dt          = Tf / n_iter;
            state_t u_sol = solve_exp( algo, dt, Tf, lambda );
            auto e        = error( u_exa, u_sol );
            errors.push_back( std::log( e ) );
            dts.push_back( std::log( dt ) );
        }

        auto [a, b] = mayor_method( dts, errors );

        return a;
    }

    template <typename Algorithm_t, typename T = double>
    T
    check_order( Algorithm_t algo = Algorithm_t(), T lambda = 0.33 )
    {
        return short_time_check_order( algo, lambda );
    }

} // namespace additive_method

namespace RDA_method
{
    template <typename Algorithm_t, typename T = double>
    auto
    solve_exp( Algorithm_t& algo, T dt, T Tf, T lambda )
    {
        using state_t = T;

        auto pb = ponio::make_problem(
            ponio::make_implicit_problem(
                [=]( T, state_t y, state_t& dy )
                {
                    dy = 0.125 * lambda * y;
                },
                [=]( T, state_t ) -> state_t
                {
                    return 0.125 * lambda;
                } ),
            [=]( T, state_t y, state_t& dy )
            {
                dy = 0.375 * lambda * y;
            },
            [=]( T, state_t y, state_t& dy )
            {
                dy = ( 1. - 0.5 * lambda ) * y;
            } );

        state_t y0                 = 1.0;
        ponio::time_span<T> t_span = { 0., Tf };

        auto obs = []( T, state_t, T ) {};
        return ::ponio::solve( pb, algo, y0, t_span, dt, obs );
    }

    template <typename Algorithm_t, typename T = double>
    T
    short_time_check_order( Algorithm_t algo = Algorithm_t(), T lambda = 0.33 )
    {
        using state_t = T;

        std::vector<T> errors;
        std::vector<T> dts;

        T Tf = 1.0;

        state_t u_exa = std::exp( Tf );

        for ( auto n_iter : { 50, 25, 20, 15, 10 } )
        {
            T dt          = Tf / n_iter;
            state_t u_sol = solve_exp( algo, dt, Tf, lambda );
            auto e        = error( u_exa, u_sol );
            errors.push_back( std::log( e ) );
            dts.push_back( std::log( dt ) );
        }

        auto [a, b] = mayor_method( dts, errors );

        return a;
    }

    template <typename Algorithm_t, typename T = double>
    T
    check_order( Algorithm_t algo = Algorithm_t(), T lambda = 0.33 )
    {
        return short_time_check_order( algo, lambda );
    }

} // namespace RDA_method

namespace splitting_method
{

    namespace detail
    {

        template <typename func_t, std::size_t... Is>
        constexpr auto
        tuple_init_impl( func_t const& f, std::index_sequence<Is...> )
        {
            return std::make_tuple( ( static_cast<void>( Is ), f )... );
        }

        template <std::size_t N, typename func_t>
        constexpr auto
        tuple_init( func_t const& f )
        {
            return tuple_init_impl( f, std::make_index_sequence<N>() );
        }

        template <typename problem_tuple_t, size_t... Is>
        auto
        make_variadic_problem_impl( problem_tuple_t& pb_tuple, std::index_sequence<Is...> )
        {
            return ponio::make_problem( std::get<Is>( pb_tuple )... );
        }

        template <typename problem_tuple_t>
        auto
        make_variadic_problem( problem_tuple_t& pb_tuple )
        {
            static constexpr std::size_t size = std::tuple_size_v<problem_tuple_t>;
            return make_variadic_problem_impl( pb_tuple, std::make_index_sequence<size>{} );
        }
    }

    template <typename Algorithm_t, typename T = double>
    auto
    long_time_check_order_odd( Algorithm_t& algo )
    {
        using state_t = std::valarray<T>;

        T alpha = 2. / 3.;
        T beta  = 4. / 3.;
        T gamma = 1.;
        T delta = 1.;

        auto pb = ponio::make_problem(
            [=]( T, auto&& u, state_t& du )
            {
                du[0] = alpha * u[0] - beta * u[0] * u[1];
                du[1] = 0.;
            },
            [=]( T, auto&& u, state_t& du )
            {
                du[0] = 0.;
                du[1] = delta * u[0] * u[1];
            },
            [=]( T, auto&& u, state_t& du )
            {
                du[0] = 0.;
                du[1] = -gamma * u[1];
            } );

        // invariant calculator
        auto V = [=]( state_t const& u ) -> T
        {
            return delta * u[0] - gamma * std::log( u[0] ) + beta * u[1] - alpha * std::log( u[1] );
        };

        T x0          = 1.9;
        state_t u_ini = { x0, x0 };
        T V_ini       = V( u_ini );

        ponio::time_span<T> t_span = { 0., 1000. };
        std::vector<T> dts         = { 0.25, 0.125, 0.1, 0.075, 0.05 }; // find a way to adapt this range to the method
        std::vector<T> relative_errors;
        std::vector<T> log_dts;

        for ( auto dt : dts )
        {
            state_t u_end = ::ponio::solve( pb, algo, u_ini, t_span, dt, []( T, state_t const&, T ) {} );
            relative_errors.push_back( std::log10( relative_error( V_ini, V( u_end ) ) ) );
            log_dts.push_back( std::log10( dt ) );
        }

        auto [a, b] = mayor_method( log_dts, relative_errors );

        return a;
    }

    template <typename Algorithm_t, typename T = double>
    auto
    long_time_check_order_even( Algorithm_t& algo )
    {
        using state_t = std::valarray<T>;

        T alpha = 2. / 3.;
        T beta  = 4. / 3.;
        T gamma = 1.;
        T delta = 1.;

        auto pb = ponio::make_problem(
            [=]( T, auto&& u, state_t& du )
            {
                du[0] = alpha * u[0] - beta * u[0] * u[1];
                du[1] = -gamma * u[1];
            },
            [=]( T, auto&& u, state_t& du )
            {
                du[0] = 0.;
                du[1] = delta * u[0] * u[1];
            } );

        // invariant calculator
        auto V = [=]( state_t const& u ) -> T
        {
            return delta * u[0] - gamma * std::log( u[0] ) + beta * u[1] - alpha * std::log( u[1] );
        };

        T x0          = 1.9;
        state_t u_ini = { x0, x0 };
        T V_ini       = V( u_ini );

        ponio::time_span<T> t_span = { 0., 1000. };
        std::vector<T> dts         = { 0.5, 0.3333, 0.25, 0.125, 0.1 }; // find a way to adapt this range to the method
        std::vector<T> relative_errors;
        std::vector<T> log_dts;

        for ( auto dt : dts )
        {
            state_t u_end = ::ponio::solve( pb, algo, u_ini, t_span, dt, []( T, state_t const&, T ) {} );
            relative_errors.push_back( std::log10( relative_error( V_ini, V( u_end ) ) ) );
            log_dts.push_back( std::log10( dt ) );
        }

        auto [a, b] = mayor_method( log_dts, relative_errors );

        return a;
    }

    template <typename Algorithm_t, typename T = double>
    T
    check_order( Algorithm_t algo = Algorithm_t() )
    {
        using algos_t = decltype( algo.algos );

        if constexpr ( std::tuple_size_v<algos_t> % 2 == 0 )
        {
            return long_time_check_order_even( algo );
        }
        else
        {
            return long_time_check_order_odd( algo );
        }
    }

    template <std::size_t I, typename splitting_tuple_t, typename T = double>
    T
    check_order_split_solve( splitting_tuple_t split_tuple = splitting_tuple_t() )
    {
        using state_t                             = T;
        static constexpr std::size_t split_length = std::tuple_size_v<decltype( split_tuple.methods )>;

        std::vector<T> errors;
        std::vector<T> dts;

        T t0 = 0.0;
        T t1 = 1.0;

        // make split_length problems for _split_solve function
        auto pb_tuple = detail::tuple_init<split_length>(
            []( T, state_t y, state_t& dy )
            {
                dy = y;
            } );
        auto pb = detail::make_variadic_problem( pb_tuple );

        state_t y_exa = std::exp( t1 );

        for ( auto n_iter : { 50, 25, 20, 15, 10 } )
        {
            T dt = ( t1 - t0 ) / n_iter;
            T y0 = 1.0;
            T y1 = 0.0;

            ponio::splitting::detail::_split_solve<I>( pb, split_tuple.methods, y0, t0, t1, dt, y1, split_tuple.info() );

            auto e = error( y_exa, y1 );
            errors.push_back( std::log( e ) );
            dts.push_back( std::log( dt ) );
        }

        auto [a, b] = mayor_method( dts, errors );

        return a;
    }

} // splitting_method

namespace diagonal_implicit_method
{

    template <typename Algorithm_t, typename T = double>
    auto
    solve_exp( Algorithm_t& algo, T dt, T Tf )
    {
        using state_t = T;

        auto pb = ponio::make_implicit_problem(
            [=]( T, state_t y, state_t& dy )
            {
                dy = y;
            },
            [=]( T, state_t ) -> state_t
            {
                return 1.;
            } );

        state_t y0                 = 1.0;
        ponio::time_span<T> t_span = { 0., Tf };

        auto obs = []( T, state_t, T ) {};
        return ::ponio::solve( pb, algo, y0, t_span, dt, obs );
    }

    template <typename Algorithm_t, typename T = double>
    T
    short_time_check_order( Algorithm_t algo = Algorithm_t() )
    {
        using state_t = T;

        std::vector<T> errors;
        std::vector<T> dts;

        T Tf = 1.0;

        state_t u_exa = std::exp( Tf );

        for ( auto n_iter : { 50, 25, 20, 15, 10 } )
        {
            T dt          = Tf / n_iter;
            state_t u_sol = solve_exp( algo, dt, Tf );
            auto e        = error( u_exa, u_sol );
            errors.push_back( std::log( e ) );
            dts.push_back( std::log( dt ) );
        }

        auto [a, b] = mayor_method( dts, errors );

        return a;
    }

    template <typename Algorithm_t, typename T = double>
    T
    check_order( Algorithm_t algo = Algorithm_t() )
    {
        return short_time_check_order( algo );
    }

} // namespace diagonal_implicit_method
