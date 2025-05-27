// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <concepts>
#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>

namespace ponio
{

    // --- SIMPLE_PROBLEM ----------------------------------------------------------
    /** @class simple_problem
     *  define a problem with a unique function
     *  @tparam Callable_t type of callable object (or function) stored in problem
     *
     *  This class represent a problem of the form \f( \dot{y}=f(t,y) \f)
     */
    template <typename Callable_t>
    struct simple_problem
    {
        Callable_t f;

        simple_problem( Callable_t&& f_ );

        template <typename state_t, typename value_t>
            requires std::invocable<Callable_t, value_t, state_t, state_t&>
        void
        operator()( value_t t, state_t&& y, state_t& dy );

        template <typename state_t, typename value_t>
            requires std::invocable<Callable_t, value_t, state_t, state_t&>
        void
        operator()( value_t t, state_t& y, state_t& dy );

        template <typename state_t, typename value_t>
            requires std::invocable<Callable_t, value_t, state_t, state_t&>
        void
        operator()( value_t t, state_t const& y, state_t& dy );

        template <typename state_t, typename value_t>
            requires std::invocable<Callable_t, value_t, state_t>
        void
        operator()( value_t t, state_t&& y, state_t& dy );

        template <typename state_t, typename value_t>
            requires std::invocable<Callable_t, value_t, state_t>
        void
        operator()( value_t t, state_t& y, state_t& dy );

        template <typename state_t, typename value_t>
            requires std::invocable<Callable_t, value_t, state_t>
        void
        operator()( value_t t, state_t const& y, state_t& dy );
    };

    /**
     * constructor of \ref simple_problem from a callable and hints
     * @param f_       callable object
     */
    template <typename Callable_t>
    inline simple_problem<Callable_t>::simple_problem( Callable_t&& f_ )
        : f( std::forward<Callable_t>( f_ ) )
    {
    }

    /**
     * call operator to evaluate \f$f(t,y)\f$
     * @param t  time value \f$t\f$
     * @param y  variable of the problem value \f$y\f$
     * @param dy \f$\dot{y}\f$, returns value of \f$f(t, y)\f$ of the problem
     */
    template <typename Callable_t>
    template <typename state_t, typename value_t>
        requires std::invocable<Callable_t, value_t, state_t, state_t&>
    inline void
    simple_problem<Callable_t>::operator()( value_t t, state_t&& y, state_t& dy )
    {
        f( t, std::forward<state_t>( y ), dy );
    }

    template <typename Callable_t>
    template <typename state_t, typename value_t>
        requires std::invocable<Callable_t, value_t, state_t, state_t&>
    inline void
    simple_problem<Callable_t>::operator()( value_t t, state_t& y, state_t& dy )
    {
        f( t, y, dy );
    }

    template <typename Callable_t>
    template <typename state_t, typename value_t>
        requires std::invocable<Callable_t, value_t, state_t, state_t&>
    inline void
    simple_problem<Callable_t>::operator()( value_t t, state_t const& y, state_t& dy )
    {
        f( t, y, dy );
    }

    template <typename Callable_t>
    template <typename state_t, typename value_t>
        requires std::invocable<Callable_t, value_t, state_t>
    inline void
    simple_problem<Callable_t>::operator()( value_t t, state_t&& y, state_t& dy )
    {
        dy = f( t, std::forward<state_t>( y ) );
    }

    template <typename Callable_t>
    template <typename state_t, typename value_t>
        requires std::invocable<Callable_t, value_t, state_t>
    inline void
    simple_problem<Callable_t>::operator()( value_t t, state_t& y, state_t& dy )
    {
        dy = f( t, y );
    }

    template <typename Callable_t>
    template <typename state_t, typename value_t>
        requires std::invocable<Callable_t, value_t, state_t>
    inline void
    simple_problem<Callable_t>::operator()( value_t t, state_t const& y, state_t& dy )
    {
        dy = f( t, y );
    }

    /**
     * factory of \ref simple_problem
     * @param c        callable object (function or functor) which represent the function of the problem
     */
    template <typename Callable_t>
    simple_problem<Callable_t>
    make_simple_problem( Callable_t&& c )
    {
        return simple_problem<Callable_t>( std::forward<Callable_t>( c ) );
    }

    // --- IMPLICIT_PROBLEM --------------------------------------------------------
    /** @class implicit_problem
     *  define a problem with its jacobian to use implicit Runge-Kutta method
     *  @tparam Callable_t type of callable object (or function) stored in problem
     *  @tparam Jacobian_t type of callable object (or function) that represents the jacobian
     */
    template <typename Callable_t, typename Jacobian_t>
    struct implicit_problem : public simple_problem<Callable_t>
    {
        using simple_problem<Callable_t>::simple_problem;

        Jacobian_t df;

        implicit_problem( Callable_t&& f_, Jacobian_t&& df_ );
    };

    /**
     * constructor of \ref implicit_problem from a callable and hints
     * @param f_       callable object that represents problem
     * @param df_      callable object that represents Jacobian of problem
     */
    template <typename Callable_t, typename Jacobian_t>
    inline implicit_problem<Callable_t, Jacobian_t>::implicit_problem( Callable_t&& f_, Jacobian_t&& df_ )
        : simple_problem<Callable_t>( std::forward<Callable_t>( f_ ) )
        , df( std::forward<Jacobian_t>( df_ ) )
    {
    }

    /**
     * factory of \ref implicit_problem
     * @param f        callable object (function or functor) which represent the function of the problem
     * @param df       jacobian of \f$f\f$ function
     */
    template <typename Callable_t, typename Jacobian_t>
    implicit_problem<Callable_t, Jacobian_t>
    make_implicit_problem( Callable_t&& f, Jacobian_t&& df )
    {
        return implicit_problem<Callable_t, Jacobian_t>( std::forward<Callable_t>( f ), std::forward<Jacobian_t>( df ) );
    }

    // --- IMPLICIT_OPERATOR_PROBLEM -----------------------------------------------
    /** @class implicit_operator_problem
     *  define a problem with its operator to use implicit Runge-Kutta method with Samurai
     *  @tparam Callable1_t type of callable object (or function) stored in problem
     *  @tparam Callable2_t type of callable object (or function) that represents the \f$f:\mapsto f(t,\cdot)\f$
     */
    template <typename Callable1_t, typename Callable2_t>
    struct implicit_operator_problem : public simple_problem<Callable1_t>
    {
        using simple_problem<Callable1_t>::simple_problem;

        Callable2_t f_t;

        implicit_operator_problem( Callable1_t&& f_, Callable2_t&& f_t_ );
    };

    /**
     * @brief Construct a new implicit operator problem<Callable1 t, Callable2 t>::implicit operator problem object
     *
     * @tparam Callable1_t
     * @tparam Callable2_t
     * @param f_   callable object that represents problem, \f$f:t,u\mapsto f(t,u)\f$
     * @param f_t_ callable object that returns operator, \f$f_t:t\mapsto f(t,\cdot)\f$
     */
    template <typename Callable1_t, typename Callable2_t>
    inline implicit_operator_problem<Callable1_t, Callable2_t>::implicit_operator_problem( Callable1_t&& f_, Callable2_t&& f_t_ )
        : simple_problem<Callable1_t>( std::forward<Callable1_t>( f_ ) )
        , f_t( std::forward<Callable2_t>( f_t_ ) )
    {
    }

    /**
     * factory of \ref implicit_operator_problem
     * @param f        callable object (function or functor) which represent the function of the problem
     * @param f_t      callable of \f$f:t\mapsto f(t,\cdot)\f$ function
     */
    template <typename Callable1_t, typename Callable2_t>
    implicit_operator_problem<Callable1_t, Callable2_t>
    make_implicit_operator_problem( Callable1_t&& f, Callable2_t&& f_t )
    {
        return implicit_operator_problem<Callable1_t, Callable2_t>( std::forward<Callable1_t>( f ), std::forward<Callable2_t>( f_t ) );
    }

    // --- IMEX_PROBLEM -----------------------------------------------------------
    /** @class imex_problem
     *  define a problem with its explicit and implicit part
     *  @tparam Callable_explicit_t type of callable object (or function) that represents the explicit part of the problem
     *  @tparam Implicit_problem_t  type of implicit problem (by operator or jacobian) that represents the implict part of the problem
     */
    template <typename Callable_explicit_t, typename Implicit_problem_t>
    struct imex_problem
    {
        Callable_explicit_t explicit_part;
        Implicit_problem_t implicit_part;

        imex_problem( Callable_explicit_t&& f_explicit, Implicit_problem_t&& pb_implicit );

        template <typename state_t, typename value_t>
        void
        operator()( value_t t, state_t&& y, state_t& dy )
        {
            state_t dy_exp = dy;
            state_t dy_imp = dy;
            explicit_part( t, std::forward<state_t>( y ), dy_exp );
            implicit_part( t, std::forward<state_t>( y ), dy_imp );
            dy = dy_exp + dy_imp;
        }

        template <typename state_t, typename value_t>
        void
        operator()( value_t t, state_t& y, state_t& dy )
        {
            state_t dy_exp = dy;
            state_t dy_imp = dy;
            explicit_part( t, y, dy_exp );
            implicit_part( t, y, dy_imp );
            dy = dy_exp + dy_imp;
        }
    };

    /**
     * @brief Construct a new imex_problem<Callable explicit t, Implicit problem t>::imex_problem object
     *
     * @tparam Callable_explicit_t
     * @tparam Implicit_problem_t
     * @param f_explicit  explicit part of problem
     * @param pb_implicit implicit part of problem (a implicit_operator_problem or a implicit_problem)
     */
    template <typename Callable_explicit_t, typename Implicit_problem_t>
    inline imex_problem<Callable_explicit_t, Implicit_problem_t>::imex_problem( Callable_explicit_t&& f_explicit,
        Implicit_problem_t&& pb_implicit )
        : explicit_part( std::forward<Callable_explicit_t>( f_explicit ) )
        , implicit_part( std::forward<Implicit_problem_t>( pb_implicit ) )
    {
    }

    // cppcheck-suppress-begin unusedFunction

    /**
     * @brief factory of imex_problem from an explicit part and a implicit part which is a implicit_operator_problem
     *
     * @tparam Callable_explicit_t
     * @tparam Callable_implicit_t
     * @tparam Callable_implicit_op_t
     * @param f explicit part
     * @param g implicit part
     * @param g_t operator on the implicit part
     */
    template <typename Callable_explicit_t, typename Callable_implicit_t, typename Callable_implicit_op_t>
    auto
    make_imex_operator_problem( Callable_explicit_t&& f, Callable_implicit_t&& g, Callable_implicit_op_t&& g_t )
    {
        return imex_problem<Callable_explicit_t, implicit_operator_problem<Callable_implicit_t, Callable_implicit_op_t>>(
            std::forward<Callable_explicit_t>( f ),
            make_implicit_operator_problem( std::forward<Callable_implicit_t>( g ), std::forward<Callable_implicit_op_t>( g_t ) ) );
    }

    template <typename Callable_explicit_t, typename Callable_implicit_t, typename Callable_implicit_jac_t>
    auto
    make_imex_jacobian_problem( Callable_explicit_t&& f, Callable_implicit_t&& g, Callable_implicit_jac_t&& dg )
    {
        return imex_problem<Callable_explicit_t, implicit_problem<Callable_implicit_t, Callable_implicit_jac_t>>(
            std::forward<Callable_explicit_t>( f ),
            make_implicit_problem( std::forward<Callable_implicit_t>( g ), std::forward<Callable_implicit_jac_t>( dg ) ) );
    }

    // cppcheck-suppress-end unusedFunction

    // --- LAWSON_PROBLEM ----------------------------------------------------------
    /** @class lawson_problem
     *  define a problem with a linear part and non-linear part
     *  @tparam Linear_t     type of the linear part (a matrix)
     *  @tparam Nonlinear_t  type of callable of the non-lineart part (function or functor)
     *
     *  This class represent a problem of the form \f( \dot{u}=Lu + N(t,u) \f)
     */
    template <typename Linear_t, typename Nonlinear_t>
    struct lawson_problem
    {
        Linear_t l;
        Nonlinear_t n;

        lawson_problem( Linear_t&& l_, Nonlinear_t&& n_ );

        template <typename state_t, typename value_t>
        void
        operator()( value_t t, state_t&& y, state_t& dy );
    };

    /**
     * constructor of \ref lawson_problem
     * @param l_       linerar part \f$L\f$ of Lawson problem
     * @param n_       nonlinerar part \f$N(t,u)\f$ of Lawson problem
     */
    template <typename Linear_t, typename Nonlinear_t>
    lawson_problem<Linear_t, Nonlinear_t>::lawson_problem( Linear_t&& l_, Nonlinear_t&& n_ )
        : l( std::forward<Linear_t>( l_ ) )
        , n( std::forward<Nonlinear_t>( n_ ) )
    {
    }

    /**
     * call operator to evaluate \f$f(t,u)\f$
     * @param t time value \f$t\f$
     * @param u variable of the problem value \f$u\f$
     * @return Returns \f$Lu + N(t,u)\f$
     */
    template <typename Linear_t, typename Nonlinear_t>
    template <typename state_t, typename value_t>
    void
    lawson_problem<Linear_t, Nonlinear_t>::operator()( value_t t, state_t&& y, state_t& dy )
    {
        n( t, std::forward<state_t>( y ), dy );
        dy = l * y + dy;
    }

    /**
     * factory of \ref lawson_problem
     * @param l        linerar part \f$L\f$ of Lawson problem
     * @param n        nonlinerar part \f$N(t,u)\f$ of Lawson problem
     */
    template <typename Linear_t, typename Nonlinear_t>
    lawson_problem<Linear_t, Nonlinear_t>
    make_lawson_problem( Linear_t&& l, Nonlinear_t&& n )
    {
        return lawson_problem<Linear_t, Nonlinear_t>( std::forward<Linear_t>( l ), std::forward<Nonlinear_t>( n ) );
    }

    // --- PROBLEM -----------------------------------------------------------------
    /** @class problem
     *  main problem that could contains multiple sub-problem
     *  @details Main goal of this class is to represent a problem of the form
     *  \f( \dot{u} = \sum_i f_i(t,u) \f)
     */
    template <typename... Callables_t>
    struct problem
    {
        static std::size_t const size = sizeof...( Callables_t );
        std::tuple<Callables_t...> system;

        problem( Callables_t... args );

        template <typename value_t, typename state_t, std::size_t... Is>
        state_t
        _sum_components_impl( value_t t, state_t&& y, state_t& dy, std::index_sequence<Is...> );

        template <typename value_t, typename state_t, std::size_t... Is>
        state_t
        _sum_components_impl( value_t t, state_t& y, state_t& dy, std::index_sequence<Is...> );

        template <typename value_t, typename state_t>
        void
        operator()( value_t t, state_t&& y, state_t& dy );

        template <typename value_t, typename state_t>
        void
        operator()( value_t t, state_t& y, state_t& dy );

        template <std::size_t I, typename value_t, typename state_t>
        void
        operator()( std::integral_constant<std::size_t, I>, value_t t, state_t&& y, state_t& dy );

        template <std::size_t I, typename value_t, typename state_t>
        void
        operator()( std::integral_constant<std::size_t, I>, value_t t, state_t& y, state_t& dy );

        template <std::size_t I, typename value_t, typename state_t>
        void
        operator()( std::integral_constant<std::size_t, I>, value_t t, state_t const& y, state_t& dy );

        template <std::size_t Index, typename value_t, typename state_t>
        state_t
        call( value_t t, state_t&& y, state_t& dy );

        template <std::size_t Index, typename value_t, typename state_t>
        state_t
        call( value_t t, state_t& y, state_t& dy );

        template <std::size_t Index, typename value_t, typename state_t>
        state_t
        call( value_t t, state_t const& y, state_t& dy );
    };

    /**
     * constructor of \ref problem from variadic number of sub-problem
     * @param args list of sub-problems
     */
    template <typename... Callables_t>
    inline problem<Callables_t...>::problem( Callables_t... args )
        : system( args... )
    {
    }

    /**
     * sum all call of each sub-problem
     * @tparam Is index sequence to iterate over tuple of sub-problems
     * @param t time \f$t\f$
     * @param u time \f$u(t)\f$
     * @return returns \f$\sum_i f_i(t,u)\f$
     */
    template <typename... Callables_t>
    template <typename value_t, typename state_t, std::size_t... Is>
    inline state_t
    problem<Callables_t...>::_sum_components_impl( value_t t, state_t&& y, state_t& dy, std::index_sequence<Is...> )
    {
        return ( call<Is>( t, std::forward<state_t>( y ), dy ) + ... );
    }

    template <typename... Callables_t>
    template <typename value_t, typename state_t, std::size_t... Is>
    inline state_t
    problem<Callables_t...>::_sum_components_impl( value_t t, state_t& y, state_t& dy, std::index_sequence<Is...> )
    {
        return ( call<Is>( t, y, dy ) + ... );
    }

    /**
     * call operator
     * @param t time \f$t\f$
     * @param u time \f$u(t)\f$
     * @return returns \f$\sum_i f_i(t,u)\f$
     */
    template <typename... Callables_t>
    template <typename value_t, typename state_t>
    inline void
    problem<Callables_t...>::operator()( value_t t, state_t&& y, state_t& dy )
    {
        dy = _sum_components_impl( t, std::forward<state_t>( y ), dy, std::make_index_sequence<size>{} );
    }

    template <typename... Callables_t>
    template <typename value_t, typename state_t>
    inline void
    problem<Callables_t...>::operator()( value_t t, state_t& y, state_t& dy )
    {
        dy = _sum_components_impl( t, y, dy, std::make_index_sequence<size>{} );
    }

    /**
     * call operator for the I operator in Callables_t
     *
     * @tparam I  Index index of tuple of sub-problems to call
     * @param t   time \f$t\f$
     * @param u   time \f$u(t)\f$
     * @return    returns \f$f_i(t,u)\f$
     */
    template <typename... Callables_t>
    template <std::size_t I, typename value_t, typename state_t>
    inline void
    problem<Callables_t...>::operator()( std::integral_constant<std::size_t, I>, value_t t, state_t&& y, state_t& dy )
    {
        call<I>( t, std::forward<state_t>( y ), dy );
    }

    template <typename... Callables_t>
    template <std::size_t I, typename value_t, typename state_t>
    inline void
    problem<Callables_t...>::operator()( std::integral_constant<std::size_t, I>, value_t t, state_t& y, state_t& dy )
    {
        call<I>( t, y, dy );
    }

    template <typename... Callables_t>
    template <std::size_t I, typename value_t, typename state_t>
    inline void
    problem<Callables_t...>::operator()( std::integral_constant<std::size_t, I>, value_t t, state_t const& y, state_t& dy )
    {
        call<I>( t, y, dy );
    }

    /**
     * call the `Index` sub-problem
     * @tparam Index index of tuple of sub-problems to call
     * @param t time \f$t\f$
     * @param u time \f$u(t)\f$
     * @return returns \f$f_i(t,u)\f$
     */
    template <typename... Callables_t>
    template <std::size_t Index, typename value_t, typename state_t>
    inline state_t
    problem<Callables_t...>::call( value_t t, state_t&& y, state_t& dy )
    {
        std::get<Index>( system )( t, std::forward<state_t>( y ), dy );
        return dy;
    }

    template <typename... Callables_t>
    template <std::size_t Index, typename value_t, typename state_t>
    inline state_t
    problem<Callables_t...>::call( value_t t, state_t& y, state_t& dy )
    {
        std::get<Index>( system )( t, y, dy );
        return dy;
    }

    template <typename... Callables_t>
    template <std::size_t Index, typename value_t, typename state_t>
    inline state_t
    problem<Callables_t...>::call( value_t t, state_t const& y, state_t& dy )
    {
        std::get<Index>( system )( t, y, dy );
        return dy;
    }

    /**
     * factory of \ref problem
     * @param f list of sub-problems
     */
    template <typename... Callables_t>
    problem<Callables_t...>
    make_problem( Callables_t... f )
    {
        return problem<Callables_t...>( f... );
    }

} // namespace ponio
