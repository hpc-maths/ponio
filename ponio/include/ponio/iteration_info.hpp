// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <array>
#include <cstddef>
#include <ranges>
#include <tuple> // NOLINT(misc-include-cleaner)

#include "detail.hpp"
#include "stage.hpp"

namespace ponio
{
    /**
     * @brief stores information for and of current iteration in algorithm
     *
     * @tparam tableau_t type of Butcher tableau
     */
    template <typename tableau_t>
    struct iteration_info
    {
        /**
         * @brief type of coefficient of Butcher tableau, same type to store error and tolerance
         *
         */
        using value_t = typename tableau_t::value_t;

        value_t error;                /**< error makes on time iteration for adaptive time step method */
        bool success = true;          /**< sets as true only for success iteration */
        bool is_step = false;         /**< sets as true only if iterator is on a step given in solver */
        std::size_t number_of_stages; /**< number of stages of method */
        std::size_t number_of_eval;   /**< number of evaluation of function */
        value_t tolerance;            /**< tolerance for the method (for adaptive time step method) */

        /**
         * @brief Construct a new iteration info object
         *
         * @param tol tolerance for adaptive time step method
         */
        iteration_info( value_t tol = static_cast<value_t>( 0 ) )
            : error( static_cast<value_t>( 0 ) )
            , number_of_stages( 0 )
            , number_of_eval( 0 )
            , tolerance( tol )
        {
        }

        iteration_info( value_t tol = static_cast<value_t>( 0 ) )
            requires stages::has_static_number_of_stages<tableau_t>
            : error( static_cast<value_t>( 0 ) )
            , number_of_stages( tableau_t::N_stages )
            , number_of_eval( 0 )
            , tolerance( tol )
        {
        }

        /**
         * @brief reset number of evaluations to zero
         *
         */
        void
        reset_eval()
        {
            number_of_eval = 0;
        }
    };

    /**
     * @brief template specialization of iteration_info for IMEX methods
     *
     * @tparam tableaus_t type of IMEX method
     */
    template <typename tableaus_t>
        requires tableaus_t::is_imex_method
    struct iteration_info<tableaus_t>
    {
        /**
         * @brief type of coefficient of Butcher tableau, same type to store error and tolerance
         *
         */
        using value_t = typename tableaus_t::value_t;

        value_t error;                                                   /**< error makes on time iteration for adaptive time step method */
        bool success = true;                                             /**< sets as true only for success iteration */
        bool is_step = false;                                            /**< sets as true only if iterator is on a step given in solver */
        std::size_t number_of_stages;                                    /**< number of stages of method */
        std::array<std::size_t, tableaus_t::N_operators> number_of_eval; /**< number of evaluation of function */
        value_t tolerance;                                               /**< tolerance for the method (for adaptive time step method) */

        /**
         * @brief Construct a new iteration info object
         *
         * @param tol tolerance for adaptive time step method
         */
        iteration_info( value_t tol = static_cast<value_t>( 0 ) )
            : error( static_cast<value_t>( 0 ) )
            , number_of_stages( 0 )
            , number_of_eval( detail::init_fill_array<tableaus_t::N_operators, std::size_t>( 0 ) )
            , tolerance( tol )
        {
        }

        iteration_info( value_t tol = static_cast<value_t>( 0 ) )
            requires stages::has_static_number_of_stages<tableaus_t>
            : error( static_cast<value_t>( 0 ) )
            , number_of_stages( tableaus_t::N_stages )
            , number_of_eval( detail::init_fill_array<tableaus_t::N_operators, std::size_t>( 0 ) )
            , tolerance( tol )
        {
        }

        /**
         * @brief reset number of evaluations to zero
         */
        void
        reset_eval()
        {
            number_of_eval.fill( 0 );
        }
    };

    namespace details
    {
        template <typename tuple_t, std::size_t... I>
        constexpr auto
        tuple_of_number_of_eval_impl( std::index_sequence<I...> )
        {
            // for each type elements of a tuple, get a `iteration_info` on this type, get a `number_of_element` instance to get its type
            return std::make_tuple( decltype( std::declval<std::tuple_element_t<I, tuple_t>>().info().number_of_eval )()... );
        }

        template <typename tuple_t>
        struct tuple_of_number_of_eval
        {
            using type = decltype( tuple_of_number_of_eval_impl<tuple_t>( std::make_index_sequence<std::tuple_size<tuple_t>{}>{} ) );
        };

        template <typename tuple_t>
        using tuple_of_number_of_eval_t = typename tuple_of_number_of_eval<tuple_t>::type;

        template <typename T>
        void
        set_to_zero( T& x )
        {
            x = 0;
        }

        template <typename T>
            requires std::ranges::range<T>
        void
        set_to_zero( T& x )
        {
            x.fill( 0 );
        }

        template <typename T>
        void
        increment( T& x, T const& y )
        {
            x += y;
        }

        template <typename T>
            requires std::ranges::range<T>
        void
        increment( T& x, T const& y )
        {
            for ( std::size_t i = 0; i < x.size(); ++i )
            {
                x[i] += y[i];
            }
        }

    } // namespace details

    /**
     * @brief template specialization of iteration_info for splitting methods
     *
     * @tparam splitting_t type of splitting method
     */
    template <typename splitting_t>
        requires splitting_t::is_splitting_method
    struct iteration_info<splitting_t>
    {
        /**
         * @brief type of coefficient of Butcher tableau, same type to store error and tolerance
         */
        using value_t = typename splitting_t::value_t;
        using tuple_t = typename splitting_t::tuple_t;

        value_t delta;        /**< parameter of shifting for adaptive time step method */
        value_t error;        /**< error makes on time iteration for adaptive time step method */
        bool success = true;  /**< sets as true only for success iteration */
        bool is_step = false; /**< sets as true only if iterator is on a step given in solver */

        std::size_t number_of_steps;                                /**< number of stages of method */
        details::tuple_of_number_of_eval_t<tuple_t> number_of_eval; /**< number of evaluation of function */

        value_t tolerance; /**< tolerance for the method (for adaptive time step method) */

        tuple_t* ptr_methods; /**< pointer to tuple of methods to access to iteration_info of each substep */

        iteration_info( tuple_t& methods, value_t delta_ = static_cast<value_t>( 0 ), value_t tol = static_cast<value_t>( 0 ) )
            : delta( delta_ )
            , error( static_cast<value_t>( 0 ) )
            , number_of_steps( splitting_t::N_steps )
            , number_of_eval()
            , tolerance( tol )
            , ptr_methods( &methods )
        {
            reset_eval();
        }

        /**
         * @brief get information on substep algorithm for splitting method
         *
         * @tparam I    rank of substep
         * @return auto
         */
        template <std::size_t I>
        auto
        get( std::integral_constant<std::size_t, I> )
        {
            return std::get<I>( *ptr_methods ).info();
        }

        template <std::size_t... Is>
        void
        reset_eval_impl( std::index_sequence<Is...> )
        {
            [[maybe_unused]] auto l = { ( details::set_to_zero( std::get<Is>( number_of_eval ) ), 0 )... };
        }

        /**
         * @brief reset number of evaluations to zero
         */
        void
        reset_eval()
        {
            reset_eval_impl( std::make_index_sequence<std::tuple_size<tuple_t>{}>{} );
        }

        /**
         * @brief increment number of evaluations of operator I
         */
        template <std::size_t I, typename evaluations_t>
        void
        increment( evaluations_t const& evals )
        {
            details::increment( std::get<I>( number_of_eval ), evals );
        }
    };

    /**
     * @brief template specialization of iteration_info for external method provide by user
     *
     * @tparam user_defined_method_t type of user method
     */
    template <typename user_defined_method_t>
        requires user_defined_method_t::is_user_defined_method
    struct iteration_info<user_defined_method_t>
    {
        /**
         * @brief type of coefficient of Butcher tableau, same type to store error and tolerance
         */
        using value_t = typename user_defined_method_t::value_t;

        value_t error;                        /**< error makes on time iteration for adaptive time step method */
        bool success                 = true;  /**< sets as true only for success iteration */
        bool is_step                 = false; /**< sets as true only if iterator is on a step given in solver */
        std::size_t number_of_stages = 0;     /**< number of stages of method */
        std::size_t number_of_eval   = 1;     /**< number of evaluation of function */
        value_t tolerance;                    /**< tolerance for the method (for adaptive time step method) */

        /**
         * @brief Construct a new iteration info object
         *
         * @param tol tolerance for adaptive time step method
         */
        iteration_info( value_t tol = static_cast<value_t>( 0 ) )
            : error( static_cast<value_t>( 0 ) )
            , tolerance( tol )
        {
        }

        /**
         * @brief reset number of evaluations to zero
         */
        void
        reset_eval()
        {
            number_of_eval = 0;
        }
    };

} // namespace ponio
