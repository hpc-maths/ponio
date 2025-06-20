// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <ranges>
#include <tuple>
#include <type_traits>

#include "detail.hpp"
#include "splitting.hpp" // NOLINT(misc-include-cleaner)
#include "stage.hpp"
#include "user_defined_method.hpp" // NOLINT(misc-include-cleaner)

namespace ponio
{
    template <typename Algorithm_t, typename state_t>
    struct method;

    template <typename Algorithm_t>
    concept is_implemented_method = stages::has_static_number_of_stages<Algorithm_t> || stages::has_dynamic_number_of_stages<Algorithm_t>;

    ///////////////////////////////////////////////////////////////////////////
    // method with static number of stages

    /** @class method
     *  @brief define a time method
     *
     *  A method is define by its algorithm (how compute \f$u^{n+1}\f$ from
     *  \f$(t^n,u^n,\Delta t)\f$) and store possible substeps. Actually this class
     *  avoid the user to specify a template parameter, while making it easy to
     *  add new methods.
     *
     *  @tparam Algorithm_t type of the algorithm which define properly the method
     *  @tparam state_t type of \f$u^n\f$
     */
    template <typename Algorithm_t, typename state_t>
        requires stages::has_static_number_of_stages<Algorithm_t>
    struct method<Algorithm_t, state_t>
    {
        static constexpr bool is_embedded = Algorithm_t::is_embedded;
        using step_storage_t              = typename std::
            conditional_t<is_embedded, std::array<state_t, Algorithm_t::N_stages + 2>, std::array<state_t, Algorithm_t::N_stages + 1>>;

        Algorithm_t alg;
        step_storage_t kis;

        /**
         * constructor of \ref method from its stages and a \f$u_0\f$ (only for preallocation)
         * of temporary substeps
         * @param alg_         a `Algorithm_t` objet with predifined stages of the method
         * @param shadow_of_u0 an object with the same size of computed value for allocation
         */
        method( Algorithm_t const& alg_, state_t const& shadow_of_u0 )
            : alg( alg_ )
            , kis( ::ponio::detail::init_fill_array<std::tuple_size<step_storage_t>::value>( shadow_of_u0 ) )
        {
        }

        /**
         * call operator which process all stages of underlying algorithm
         * @param f  callable obect which represents the problem to solve
         * @param tn time \f$t^n\f$ last time where solution is computed
         * @param un computed solution \f$u^n\f$ à time \f$t^n\f$
         * @param dt time step
         * @return tuple of result of iteration \f$(t^{n+1},u^{n+1},\Delta t_{opt})\f$ if
         * iteration is accepted, \f$(t^{n},u^{n},\Delta t_{opt})\f$ otherwise. If `Algorithm_t`
         * is a constant time step method, so \f$\Delta t_{opt}\f$ is alway equal to the same
         * initial value.
         */
        template <typename Problem_t, typename value_t>
        std::tuple<value_t, state_t, value_t>
        operator()( Problem_t& f, value_t tn, state_t& un, value_t dt )
        {
            _call_stage( f, tn, un, dt );

            return _return( tn, un, dt );
        }

        // NOLINTBEGIN(modernize-type-traits,modernize-use-constraints)
        // TODO: change to get expression (I == N_stages+1) into a requires expression

        template <std::size_t I = 0, typename Problem_t, typename value_t, typename Algo_t = Algorithm_t>
            requires std::same_as<Algo_t, Algorithm_t> && Algorithm_t::is_embedded
        typename std::enable_if<( I == Algorithm_t::N_stages + 1 ), void>::type
        _call_stage( Problem_t& f, value_t tn, state_t& un, value_t dt )
        {
            alg.stage( Stage<I>{}, f, tn, un, kis, dt, kis[I] );
        }

        template <std::size_t I = 0, typename Problem_t, typename value_t, typename Algo_t = Algorithm_t>
            requires std::same_as<Algo_t, Algorithm_t>
        typename std::enable_if<( I == Algorithm_t::N_stages + 1 ), void>::type
        _call_stage( Problem_t&, value_t, state_t&, value_t )
        {
        }

        /**
         * unroll all stages of `Algorithm_t` with templated recursion
         * @tparam I stage of `Algorithm_t` to compute
         * @param f  problem whose solution must be computed
         * @param tn time \f$t^n\f$ last time where solution is computed
         * @param un computed solution \f$u^n\f$ à time \f$t^n\f$
         * @param dt time step
         * @return this function store its result in specific attribut of \ref method
         */
        template <std::size_t I = 0, typename Problem_t, typename value_t, typename Algo_t = Algorithm_t>
            requires std::same_as<Algo_t, Algorithm_t>
        typename std::enable_if<( I < Algorithm_t::N_stages + 1 ), void>::type
        _call_stage( Problem_t& f, value_t tn, state_t& un, value_t dt )
        {
            alg.stage( Stage<I>{}, f, tn, un, kis, dt, kis[I] );
            _call_stage<I + 1>( f, tn, un, dt );
        }

        // NOLINTEND(modernize-type-traits,modernize-use-constraints)

        /**
         * return values \f$(t^n,u^n,\Delta t)\f$ after call of all stages
         * @param tn time at the begining of the step
         * @param un state at the begining of the step
         * @param dt time step of the step
         * @return return \f$(t^{n+1},u^{n+1},\Delta t_{opt})\f$ if iteration is accepted,
         * and return \f$(t^{n},u^{n},\Delta t_{opt})\f$ otherwise.
         * @details This member function differs if the algorithm is adaptive time stepping or not.
         */
        template <typename value_t, typename Algo_t = Algorithm_t>
        std::tuple<value_t, state_t, value_t>
        _return( value_t tn, [[maybe_unused]] state_t const& un, value_t dt )
        {
            return std::forward_as_tuple( tn + dt, kis.back(), dt );
        }

        template <typename value_t, typename Algo_t = Algorithm_t>
            requires std::same_as<Algo_t, Algorithm_t> && Algorithm_t::is_embedded
        std::tuple<value_t, state_t, value_t>
        _return( value_t tn, state_t const& un, value_t dt )
        {
            alg.info().error = ::ponio::detail::error_estimate( un, kis[Algorithm_t::N_stages], kis[Algorithm_t::N_stages + 1] );
            // std::cout << "alg.info().error = " << alg.info().error << std::endl;

            value_t new_dt = 0.9 * std::pow( alg.info().tolerance / alg.info().error, 1. / static_cast<value_t>( Algorithm_t::order ) ) * dt;
            new_dt = std::min( std::max( 0.2 * dt, new_dt ), 5. * dt );

            if ( alg.info().error > alg.info().tolerance )
            {
                alg.info().success = false;
                return std::forward_as_tuple( tn, un, new_dt );
            }

            alg.info().success = true;
            return std::forward_as_tuple( tn + dt, kis[Algorithm_t::N_stages], new_dt );
        }

        /**
         * @brief returns iteration_info object on algorithm
         */
        auto&
        info()
        {
            return alg.info();
        }

        /**
         * @brief returns iteration_info object on algorithm
         */
        auto const&
        info() const
        {
            return alg.info();
        }

        /**
         * @brief returns array of stages
         *
         * @return auto&
         */
        auto&
        stages()
        {
            return kis;
        }

        /**
         * @brief returns array of stages
         *
         * @return auto const&
         */
        auto const&
        stages() const
        {
            return kis;
        }
    };

    ///////////////////////////////////////////////////////////////////////////
    // method with dynamic number of stages

    template <typename Algorithm_t, typename state_t>
        requires stages::has_dynamic_number_of_stages<Algorithm_t>
    struct method<Algorithm_t, state_t>
    {
        static constexpr bool is_embedded = Algorithm_t::is_embedded;
        using step_storage_t              = std::array<state_t, Algorithm_t::N_storage>;

        Algorithm_t alg;
        step_storage_t kis;

        method( Algorithm_t const& alg_, state_t const& shadow_of_u0 )
            : alg( alg_ )
            , kis( ::ponio::detail::init_fill_array<std::tuple_size<step_storage_t>::value>( shadow_of_u0 ) )
        {
        }

        template <typename Problem_t, typename value_t>
        std::tuple<value_t, state_t, value_t>
        operator()( Problem_t& f, value_t tn, state_t& un, value_t dt )
        {
            return alg( f, tn, un, kis, dt );
        }

        /**
         * @brief returns iteration_info object on algorithm
         */
        auto&
        info()
        {
            return alg.info();
        }

        /**
         * @brief returns iteration_info object on algorithm
         */
        auto const&
        info() const
        {
            return alg.info();
        }

        /**
         * @brief returns array of stages
         *
         * @return auto&
         */
        auto&
        stages()
        {
            return kis;
        }

        /**
         * @brief returns array of stages
         *
         * @return auto const&
         */
        auto const&
        stages() const
        {
            return kis;
        }
    };

    ///////////////////////////////////////////////////////////////////////////
    // method defined by user

    template <typename user_defined_algorithm_t>
    concept is_user_method = is_user_defined_method<user_defined_algorithm_t>;

    template <typename user_defined_algorithm_t, typename state_t>
        requires is_user_method<user_defined_algorithm_t>
    struct method<user_defined_algorithm_t, state_t>
    {
        static constexpr bool is_embedded = false;

        user_defined_algorithm_t alg;

        method( user_defined_algorithm_t const& user_defined_algorithm, state_t const& )
            : alg( user_defined_algorithm )
        {
        }

        template <typename Problem_t, typename value_t>
        std::tuple<value_t, state_t, value_t>
        operator()( Problem_t& f, value_t tn, state_t& un, value_t dt )
        {
            return alg( f, tn, un, dt );
        }

        /**
         * @brief returns iteration_info object on algorithm
         */
        auto&
        info()
        {
            return alg.info();
        }

        /**
         * @brief returns iteration_info object on algorithm
         */
        auto const&
        info() const
        {
            return alg.info();
        }

        /**
         * @brief returns array of stages
         *
         * @return auto&
         */
        auto&
        stages()
        {
            return std::array<state_t, 0>{};
        }

        /**
         * @brief returns array of stages
         *
         * @return auto const&
         */
        auto const&
        stages() const
        {
            return std::array<state_t, 0>{};
        }
    };

    ///////////////////////////////////////////////////////////////////////////

    /**
     *  generic factory to build a method from an algoritm, it only reuses `method`
     *  constructor
     *  @tparam value_t type of coefficients
     *  @param algo         a `Algorithm_t` objet with predifined stages of the method
     *  @param shadow_of_u0 an object with the same size of computed value for allocation
     */
    template <typename value_t, typename Algorithm_t, typename state_t>
        requires is_implemented_method<Algorithm_t>
    auto
    make_method( Algorithm_t const& algo, state_t const& shadow_of_u0 )
    {
        return method<Algorithm_t, state_t>( algo, shadow_of_u0 );
    }

    /**
     * factory of tuple of methods from a tuple of `Algorithm_t`
     * @param algos        a tuple of `Algorithm_t` objets with predifined stages
     * @param shadow_of_u0 an object with the same size of computed value for allocation
     * @details this factory is to prevent duplucation of code in factory of methods for
     * splitting methods (Lie or Strang method).
     */
    template <typename value_t, typename state_t, typename... Algorithms_t>
    auto
    make_tuple_methods( std::tuple<Algorithms_t...> const& algos, state_t const& shadow_of_u0 )
    {
        return std::apply(
            [&]( auto const&... args )
            {
                auto maker = [&]( auto const& arg )
                {
                    return make_method<value_t>( arg, shadow_of_u0 ); // maybe should use std::bind
                };
                return std::make_tuple( maker( args )... );
            },
            algos );
    }

    /**
     * @brief generic factory for splitting methods to build a method from a tuple of algorithms and a state
     *
     * @tparam _splitting_method_t splitting method (Lie or Strang splitting)
     * @tparam value_t             type of coefficients and time step
     * @tparam state_t             type of state
     * @tparam optional_args_t     type of tuple of optional arguments to build _splitting_method_t object (void if not needed)
     * @tparam Algorithms_t        types of algorithms to solve each step of splitting
     * @param algos        tuple of algorithms and splitting method
     * @param shadow_of_u0 an object with the same sixe of computed value for allocation
     */
    template <typename value_t, template <typename, typename...> typename _splitting_method_t, typename state_t, typename optional_args_t, typename... Algorithms_t>
    auto
    make_method( splitting::detail::splitting_tuple<_splitting_method_t, value_t, optional_args_t, Algorithms_t...> const& algos,
        state_t const& shadow_of_u0 )
    {
        using splitting_tuple = splitting::detail::splitting_tuple<_splitting_method_t, value_t, optional_args_t, Algorithms_t...>;

        auto methods = make_tuple_methods<value_t>( algos.algos, shadow_of_u0 );

        if constexpr ( splitting_tuple::has_optional_args )
        {
            return splitting::detail::make_splitting_from_tuple<_splitting_method_t>( methods, algos.time_steps, algos.optional_arguments );
        }
        else
        {
            return splitting::detail::make_splitting_from_tuple<_splitting_method_t>( methods, algos.time_steps );
        }
    }

    /**
     * @brief helper function to build a method from a `user_defined_method`
     *
     * @tparam value_t               type of coefficients
     * @tparam user_defined_method_t type of user defined method
     * @tparam state_t               type of current state
     * @param u_meth       user defined method with the underlying user function
     * @param shadow_of_u0 an object with the same sixe of computed value for allocation
     */
    template <typename value_t, typename user_defined_method_t, typename state_t>
        requires is_user_defined_method<user_defined_method_t>
    auto
    make_method( user_defined_method_t const& u_meth, state_t const& shadow_of_u0 )
    {
        auto algo = make_user_defined_algorithm<value_t>( u_meth );
        return method<decltype( algo ), state_t>( algo, shadow_of_u0 );
    }

} // namespace ponio
