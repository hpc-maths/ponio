// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <bitset>
#include <tuple>

namespace ponio
{

    // --- PARENT_PROBLEM ----------------------------------------------------------
    /** @class parent_problem
     *  parent class of each sub-problem, simple propose of store hint about problen
     *  @details The main goal of this class is propose uniform hints on each problem
     *  to select the best method. In particular it propose the minimum order of
     *  time integrator (if solving method is choosen automaticaly) and indicator
     *  if problem is sitff. All hints are stored in a `std::bitset` for the state
     *  of computer science.
     */
    struct parent_problem
    {
        // std::bitset<3 + 1> property;
        // 3 bits for minimal order and 1 bit for stiffness

        // parent_problem( std::size_t order = 3, bool is_stiff = false );

        // void
        // set_order( std::size_t o );

        // unsigned long
        // order();

        // void
        // set_stiffness( bool s );

        // bool
        // stiffness();
    };

    /**
     * constructor of \ref parent_problem from hint of minimal order and stiffness
     * of problem
     * @param order    minimum order for solving the probleme
     * @param is_stiff hint on stiffness of considering problem
     */
    // inline parent_problem::parent_problem( std::size_t order, bool is_stiff )
    //     : property( ( order << 1 ) | is_stiff )
    // {
    // }

    /**
     * change the hint of order
     * @param o new minimum order of solving method
     */
    // inline void
    // parent_problem::set_order( std::size_t o )
    // {
    //     property = ( o << 1 ) | property[0];
    // }

    /**
     * get the hint on order
     * @return minimum order of time integrator choosen for solving the method
     */
    // inline unsigned long
    // parent_problem::order()
    // {
    //     return ( property >> 1 ).to_ulong();
    // }

    /**
     * set the stiffness hint on a specific value
     * @param s hint is problem is stiff
     */
    // inline void
    // parent_problem::set_stiffness( bool s )
    // {
    //     property[0] = s;
    // }

    /**
     * get the hint on stifness
     * @return Returns `true` if problem is indicate as stiff
     */
    // inline bool
    // parent_problem::stiffness()
    // {
    //     return static_cast<bool>( property[0] );
    // }

    // --- SIMPLE_PROBLEM ----------------------------------------------------------
    /** @class simple_problem
     *  define a problem with a unique function
     *  @tparam Callable_t type of callable object (or function) stored in problem
     *
     *  This class represent a problem of the form \f( \dot{u}=f(t,u) \f)
     */
    template <typename Callable_t>
    struct simple_problem : public parent_problem
    {
        // using parent_problem::parent_problem;
        Callable_t f;

        simple_problem( Callable_t& f_ );

        template <typename state_t, typename value_t>
        state_t
        operator()( value_t t, state_t const& u );
    };

    /**
     * constructor of \ref simple_problem from a callable and hints
     * @param f_       callable object
     */
    template <typename Callable_t>
    inline simple_problem<Callable_t>::simple_problem( Callable_t& f_ )
        : f( f_ )
    {
    }

    /**
     * call operator to evaluate \f$f(t,u)\f$
     * @param t time value \f$t\f$
     * @param u variable of the problem value \f$u\f$
     * @return Returns \f$f(t,u)\f$
     */
    template <typename Callable_t>
    template <typename state_t, typename value_t>
    inline state_t
    simple_problem<Callable_t>::operator()( value_t t, state_t const& u )
    {
        return f( t, u );
    }

    /**
     * factory of \ref simple_problem
     * @param c        callable object (function or functor) which represent the function of the problem
     */
    template <typename Callable_t>
    simple_problem<Callable_t>
    make_simple_problem( Callable_t c )
    {
        return simple_problem<Callable_t>( c );
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

        implicit_problem( Callable_t& f_, Jacobian_t& df_ );
    };

    /**
     * constructor of \ref implicit_problem from a callable and hints
     * @param f_       callable object
     */
    template <typename Callable_t, typename Jacobian_t>
    inline implicit_problem<Callable_t, Jacobian_t>::implicit_problem( Callable_t& f_, Jacobian_t& df_ )
        : simple_problem<Callable_t>( f_ )
        , df( df_ )
    {
    }

    /**
     * factory of \ref implicit_problem
     * @param f        callable object (function or functor) which represent the function of the problem
     * @param df       jacobian of \f$f\f$ function
     */
    template <typename Callable_t, typename Jacobian_t>
    implicit_problem<Callable_t, Jacobian_t>
    make_implicit_problem( Callable_t f, Jacobian_t df )
    {
        return implicit_problem<Callable_t, Jacobian_t>( f, df );
    }

    // --- LAWSON_PROBLEM ----------------------------------------------------------
    /** @class lawson_problem
     *  define a problem with a linear part and non-linear part
     *  @tparam Linear_t     type of the linear part (a matrix)
     *  @tparam Nonlinear_t  type of callable of the non-lineart part (function or functor)
     *
     *  This class represent a problem of the form \f( \dot{u}=Lu + N(t,u) \f)
     */
    template <typename Linear_t, typename Nonlinear_t>
    struct lawson_problem : public parent_problem
    {
        // using parent_problem::parent_problem;
        Linear_t l;
        Nonlinear_t n;

        lawson_problem( Linear_t& l_, Nonlinear_t& n_ );

        template <typename state_t, typename value_t>
        state_t
        operator()( value_t t, state_t const& u );
    };

    /**
     * constructor of \ref lawson_problem
     * @param l_       linerar part \f$L\f$ of Lawson problem
     * @param n_       nonlinerar part \f$N(t,u)\f$ of Lawson problem
     */
    template <typename Linear_t, typename Nonlinear_t>
    lawson_problem<Linear_t, Nonlinear_t>::lawson_problem( Linear_t& l_, Nonlinear_t& n_ )
        : l( l_ )
        , n( n_ )
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
    state_t
    lawson_problem<Linear_t, Nonlinear_t>::operator()( value_t t, state_t const& u )
    {
        return l * u + n( t, u );
    }

    /**
     * factory of \ref lawson_problem
     * @param l        linerar part \f$L\f$ of Lawson problem
     * @param n        nonlinerar part \f$N(t,u)\f$ of Lawson problem
     */
    template <typename Linear_t, typename Nonlinear_t>
    lawson_problem<Linear_t, Nonlinear_t>
    make_lawson_problem( Linear_t l, Nonlinear_t n )
    {
        return lawson_problem<Linear_t, Nonlinear_t>( l, n );
    }

    // --- IMEX_PROBLEM ------------------------------------------------------------
    /** @class imex_problem
     *  define a problem with a easy to implicit part and part which will be solve explicitly
     *  @tparam Implicit_t  type of the easy to implicit part
     *  @tparam Explicit_t  type of callable of the explicit part (function or functor)
     *
     *  This class represent a problem of the form \f( \dot{u}=I(t,u) + E(t,u) \f)
     */
    template <typename Implicit_t, typename Explicit_t>
    struct imex_problem : public parent_problem
    {
        using parent_problem::parent_problem;
        Implicit_t i;
        Explicit_t e;

        imex_problem( Implicit_t& i_, Explicit_t& e_ );

        template <typename state_t, typename value_t>
        state_t
        operator()( value_t t, state_t const& u );
    };

    /**
     * constructor of \ref imex_problem
     * @param i_       easy to implicit part \f$I(t,u)\f$ of IMEX problem
     * @param e_       explicit part \f$E(t,u)\f$ of IMEX problem
     */
    template <typename Implicit_t, typename Explicit_t>
    imex_problem<Implicit_t, Explicit_t>::imex_problem( Implicit_t& i_, Explicit_t& e_ )
        : i( i_ )
        , e( e_ )
    {
    }

    /**
     * call operator to evaluate \f$f(t,u)\f$
     * @param t time value \f$t\f$
     * @param u variable of the problem value \f$u\f$
     * @return Returns \f$I(t,u) + E(t,u)\f$
     */
    template <typename Implicit_t, typename Explicit_t>
    template <typename state_t, typename value_t>
    state_t
    imex_problem<Implicit_t, Explicit_t>::operator()( value_t t, state_t const& u )
    {
        return i * u + e( t, u );
    }

    /**
     * factory of \ref imex_problem
     * @param i        easy to implicit part \f$I(t,u)\f$ of IMEX problem
     * @param e        explicit part \f$E(t,u)\f$ of IMEX problem
     */
    template <typename Implicit_t, typename Explicit_t>
    imex_problem<Implicit_t, Explicit_t>
    make_imex_problem( Implicit_t i, Explicit_t e )
    {
        return imex_problem<Implicit_t, Explicit_t>( i, e );
    }

    // --- PROBLEM -----------------------------------------------------------------
    /** @class problem
     *  main problem that could contains multiple sub-problem
     *  @todo this class has no inheristence from \ref parent_problem so it has no
     *  hints (because of variadic template... can't have default value parameter
     *  also)
     *  @details Main goal of this class is to represent a problem of the form
     *  \f( \dot{u} = \sum_i f_i(t,u) \f)
     */
    template <typename... Callables_t>
    struct problem
    {
        static const std::size_t size = sizeof...( Callables_t );
        std::tuple<Callables_t...> system;

        problem( Callables_t... args );

        template <typename value_t, typename state_t, std::size_t... Is>
        state_t
        _sum_components_impl( value_t t, state_t const& u, std::index_sequence<Is...> );

        template <typename value_t, typename state_t>
        state_t
        operator()( value_t t, state_t const& u );

        template <std::size_t Index, typename value_t, typename state_t>
        state_t
        call( value_t t, state_t const& u );
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
    problem<Callables_t...>::_sum_components_impl( value_t t, state_t const& u, std::index_sequence<Is...> )
    {
        return ( call<Is>( t, u ) + ... );
    }

    /**
     * call operator
     * @param t time \f$t\f$
     * @param u time \f$u(t)\f$
     * @return returns \f$\sum_i f_i(t,u)\f$
     */
    template <typename... Callables_t>
    template <typename value_t, typename state_t>
    inline state_t
    problem<Callables_t...>::operator()( value_t t, state_t const& u )
    {
        return _sum_components_impl( t, u, std::make_index_sequence<size>{} );
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
    problem<Callables_t...>::call( value_t t, state_t const& u )
    {
        return std::get<Index>( system )( t, u );
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
