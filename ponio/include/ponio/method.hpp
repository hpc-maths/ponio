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
        static constexpr std::size_t
            step_storage_size = detail::conditional_v<is_embedded, std::size_t, Algorithm_t::N_stages + 2, Algorithm_t::N_stages + 1>;
        using step_storage_t  = std::array<state_t, step_storage_size>;

        Algorithm_t alg;
        step_storage_t kis;
        state_t ui;

        /**
         * constructor of \ref method from its stages and a \f$u_0\f$ (only for preallocation)
         * of temporary substeps
         * @param alg_         a `Algorithm_t` objet with predifined stages of the method
         * @param shadow_of_u0 an object with the same size of computed value for allocation
         */
        method( Algorithm_t alg_, state_t const& shadow_of_u0 )
            : alg( std::move( alg_ ) )
            , kis( ::ponio::detail::init_fill_array<std::tuple_size_v<step_storage_t>>( shadow_of_u0 ) )
            , ui( shadow_of_u0 )
        {
        }

        method() = default;

        /**
         * call operator which process all stages of underlying algorithm
         * @param f    callable obect which represents the problem to solve
         * @param tn   time \f$t^n\f$ last time where solution is computed
         * @param un   computed solution \f$u^n\f$ à time \f$t^n\f$
         * @param dt   time step
         * @param unp1 computed solution \f$u^{n+1}\f$ à time \f$t^{n+1}\f$
         * @details Current time `tn`, time step `dt` and state `unp1` are updated. If this is an adaptive time step method, and the
         * iteration failed with time step `dt`, `unp1` is step to initial solution `un` and current time `tn` isn't updated.
         */
        template <typename Problem_t, typename value_t>
        void
        operator()( Problem_t& f, value_t& tn, state_t& un, value_t& dt, state_t& unp1 )
        {
            _call_stage( f, tn, un, dt );

            _return( tn, un, dt, unp1 );
        }

        // NOLINTBEGIN(modernize-type-traits,modernize-use-constraints)
        // TODO: change to get expression (I == N_stages+1) into a requires expression

        template <std::size_t I = 0, typename Problem_t, typename value_t, typename Algo_t = Algorithm_t>
            requires std::same_as<Algo_t, Algorithm_t> && Algorithm_t::is_embedded
        typename std::enable_if<( I == Algorithm_t::N_stages + 1 ), void>::type
        _call_stage( Problem_t& f, value_t tn, state_t& un, value_t dt )
        {
            alg.stage( Stage<I>{}, f, tn, un, kis, dt, ui, kis[I] );
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
            alg.stage( Stage<I>{}, f, tn, un, kis, dt, ui, kis[I] );
            _call_stage<I + 1>( f, tn, un, dt );
        }

        // NOLINTEND(modernize-type-traits,modernize-use-constraints)

        /**
         * return values \f$(t^n,u^n,\Delta t)\f$ after call of all stages
         * @param tn   time at the begining of the step
         * @param un   state at the begining of the step
         * @param dt   time step of the step
         * @param unp1 state at the begining of the step
         * @details This member function differs if the algorithm is adaptive time stepping or not.
         */
        template <typename value_t, typename Algo_t = Algorithm_t>
        void
        _return( value_t& tn, [[maybe_unused]] state_t& un, value_t& dt, state_t& unp1 )
        {
            tn = tn + dt;
            std::swap( kis.back(), unp1 );
        }

        template <typename value_t, typename Algo_t = Algorithm_t>
            requires std::same_as<Algo_t, Algorithm_t> && Algorithm_t::is_embedded
        void
        _return( value_t& tn, state_t& un, value_t& dt, state_t& unp1 )
        {
            alg.info().error = ::ponio::detail::error_estimate( un,
                kis[Algorithm_t::N_stages],
                kis[Algorithm_t::N_stages + 1],
                info().absolute_tolerance,
                info().relative_tolerance );
            // std::cout << "alg.info().error = " << alg.info().error << std::endl;

            value_t new_dt = 0.9 * std::pow( alg.info().tolerance / alg.info().error, 1. / static_cast<value_t>( Algorithm_t::order ) ) * dt;
            new_dt = std::min( std::max( 0.2 * dt, new_dt ), 5. * dt );

            if ( alg.info().error > static_cast<value_t>( 1.0 ) )
            {
                alg.info().success = false;

                // tn = tn;
                std::swap( un, unp1 );
                dt = new_dt;
            }
            else
            {
                alg.info().success = true;

                tn = tn + dt;
                std::swap( kis[Algorithm_t::N_stages], unp1 );
                dt = new_dt;
            }
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
            , kis( ::ponio::detail::init_fill_array<std::tuple_size_v<step_storage_t>>( shadow_of_u0 ) )
        {
        }

        /**
         * call operator which process the underlying algorithm
         * @param f    callable obect which represents the problem to solve
         * @param tn   time \f$t^n\f$ last time where solution is computed
         * @param un   computed solution \f$u^n\f$ à time \f$t^n\f$
         * @param dt   time step
         * @param unp1 computed solution \f$u^{n+1}\f$ à time \f$t^{n+1}\f$
         */
        template <typename Problem_t, typename value_t>
        void
        operator()( Problem_t& f, value_t& tn, state_t& un, value_t& dt, state_t& unp1 )
        {
            return alg( f, tn, un, kis, dt, unp1 );
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
    // method for additive Runge-Kutta

    template <typename Algorithm_t, typename state_t>
        requires stages::has_static_number_of_stages<Algorithm_t> && Algorithm_t::is_imex_method
    struct method<Algorithm_t, state_t>
    {
        static constexpr bool is_embedded = Algorithm_t::is_embedded;
        static constexpr std::size_t
            step_storage_size = detail::conditional_v<is_embedded, std::size_t, Algorithm_t::N_stages + 2, Algorithm_t::N_stages + 1>;
        using step_storage_t  = std::array<state_t, step_storage_size>;

        Algorithm_t alg;
        step_storage_t k_ex_is;
        step_storage_t k_im_is;
        state_t ui;
        state_t u_tmp;

        method( Algorithm_t const& alg_, state_t const& shadow_of_u0 )
            : alg( alg_ )
            , k_ex_is( ::ponio::detail::init_fill_array<std::tuple_size_v<step_storage_t>>( shadow_of_u0 ) )
            , k_im_is( ::ponio::detail::init_fill_array<std::tuple_size_v<step_storage_t>>( shadow_of_u0 ) )
            , ui( shadow_of_u0 )
            , u_tmp( shadow_of_u0 )
        {
        }

        method() = default;

        /**
         * call operator which process all stages of underlying algorithm
         * @param f    callable obect which represents the problem to solve
         * @param tn   time \f$t^n\f$ last time where solution is computed
         * @param un   computed solution \f$u^n\f$ à time \f$t^n\f$
         * @param dt   time step
         * @param unp1 computed solution \f$u^{n+1}\f$ à time \f$t^{n+1}\f$
         * @details Current time `tn`, time step `dt` and state `unp1` are updated. If this is an adaptive time step method, and the
         * iteration failed with time step `dt`, `unp1` is step to initial solution `un` and current time `tn` isn't updated.
         */
        template <typename Problem_t, typename value_t>
        void
        operator()( Problem_t& f, value_t& tn, state_t& un, value_t& dt, state_t& unp1 )
        {
            _call_stage( f, tn, un, dt );

            _return( tn, un, dt, unp1 );
        }

        // NOLINTBEGIN(modernize-type-traits,modernize-use-constraints)
        // TODO: change to get expression (I == N_stages+1) into a requires expression

        template <std::size_t I = 0, typename Problem_t, typename value_t, typename Algo_t = Algorithm_t>
            requires std::same_as<Algo_t, Algorithm_t> && Algorithm_t::is_embedded
        typename std::enable_if<( I == Algorithm_t::N_stages + 1 ), void>::type
        _call_stage( Problem_t& f, value_t tn, state_t& un, value_t dt )
        {
            alg.stage( Stage<I>{}, f, tn, un, k_ex_is, k_im_is, dt, ui, u_tmp, k_ex_is[I], k_im_is[I] );
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
            alg.stage( Stage<I>{}, f, tn, un, k_ex_is, k_im_is, dt, ui, u_tmp, k_ex_is[I], k_im_is[I] );
            _call_stage<I + 1>( f, tn, un, dt );
        }

        // NOLINTEND(modernize-type-traits,modernize-use-constraints)

        /**
         * return values \f$(t^n,u^n,\Delta t)\f$ after call of all stages
         * @param tn   time at the begining of the step
         * @param un   state at the begining of the step
         * @param dt   time step of the step
         * @param unp1 state at the begining of the step
         * @details This member function differs if the algorithm is adaptive time stepping or not.
         */
        template <typename value_t, typename Algo_t = Algorithm_t>
        void
        _return( value_t& tn, [[maybe_unused]] state_t& un, value_t& dt, state_t& unp1 )
        {
            tn = tn + dt;
            std::swap( k_ex_is.back(), unp1 );
        }

        template <typename value_t, typename Algo_t = Algorithm_t>
            requires std::same_as<Algo_t, Algorithm_t> && Algorithm_t::is_embedded
        void
        _return( value_t& tn, state_t& un, value_t& dt, state_t& unp1 )
        {
            alg.info().error = ::ponio::detail::error_estimate( un,
                k_ex_is[Algorithm_t::N_stages],
                k_ex_is[Algorithm_t::N_stages + 1],
                info().absolute_tolerance,
                info().relative_tolerance );
            // std::cout << "alg.info().error = " << alg.info().error << std::endl;

            value_t new_dt = 0.9 * std::pow( alg.info().tolerance / alg.info().error, 1. / static_cast<value_t>( Algorithm_t::order ) ) * dt;
            new_dt = std::min( std::max( 0.2 * dt, new_dt ), 5. * dt );

            if ( alg.info().error > static_cast<value_t>( 1.0 ) )
            {
                alg.info().success = false;

                // tn = tn;
                std::swap( un, unp1 );
                dt = new_dt;
            }
            else
            {
                alg.info().success = true;

                tn = tn + dt;
                std::swap( k_ex_is[Algorithm_t::N_stages], unp1 );
                dt = new_dt;
            }
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

        // /**
        //  * @brief returns array of stages
        //  *
        //  * @return auto&
        //  */
        // auto&
        // stages_ex()
        // {
        //     return k_ex_is;
        // }

        // /**
        //  * @brief returns array of stages
        //  *
        //  * @return auto const&
        //  */
        // auto const&
        // stages_ex() const
        // {
        //     return k_ex_is;
        // }

        // /**
        //  * @brief returns array of stages
        //  *
        //  * @return auto&
        //  */
        // auto&
        // stages_im()
        // {
        //     return k_im_is;
        // }

        // /**
        //  * @brief returns array of stages
        //  *
        //  * @return auto const&
        //  */
        // auto const&
        // stages_im() const
        // {
        //     return k_im_is;
        // }
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
        void
        operator()( Problem_t& f, value_t& tn, state_t& un, value_t& dt, state_t& unp1 )
        {
            std::tie( tn, unp1, dt ) = alg( f, tn, un, dt );
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
