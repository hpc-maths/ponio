// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <array>
#include <cstddef>
#include <tuple>

#include "butcher_tableau.hpp"
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
        bool success;                 /**< sets as true only for success iteration */
        bool is_step;                 /**< sets as true only if iterator is on a step given in solver */
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
            , success( true )
            , is_step( false )
            , number_of_stages( 0 )
            , number_of_eval( 0 )
            , tolerance( tol )
        {
        }

        iteration_info( value_t tol = static_cast<value_t>( 0 ) )
            requires stages::has_static_number_of_stages<tableau_t>
            : error( static_cast<value_t>( 0 ) )
            , success( true )
            , is_step( false )
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
        bool success;                                                    /**< sets as true only for success iteration */
        bool is_step;                                                    /**< sets as true only if iterator is on a step given in solver */
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
            , success( true )
            , is_step( false )
            , number_of_stages( 0 )
            , number_of_eval( ::detail::init_fill_array<tableaus_t::N_operators, std::size_t>( 0 ) )
            , tolerance( tol )
        {
        }

        iteration_info( value_t tol = static_cast<value_t>( 0 ) )
            requires stages::has_static_number_of_stages<tableaus_t>
            : error( static_cast<value_t>( 0 ) )
            , success( true )
            , is_step( false )
            , number_of_stages( tableaus_t::N_stages )
            , number_of_eval( ::detail::init_fill_array<tableaus_t::N_operators, std::size_t>( 0 ) )
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
            number_of_eval.fill( 0 );
        }
    };

    template <typename splitting_t>
        requires splitting_t::is_splitting_method
    struct iteration_info<splitting_t>
    {
        /**
         * @brief type of coefficient of Butcher tableau, same type to store error and tolerance
         *
         */
        using value_t = typename splitting_t::value_t;
        using tuple_t = typename splitting_t::tuple_t;

        value_t error; /**< error makes on time iteration for adaptive time step method */
        bool success;  /**< sets as true only for success iteration */
        bool is_step;  /**< sets as true only if iterator is on a step given in solver */

        std::size_t number_of_steps;                                    /**< number of stages of method */
        std::array<std::size_t, splitting_t::N_methods> number_of_eval; /**< number of evaluation of function */

        value_t tolerance; /**< tolerance for the method (for adaptive time step method) */

        tuple_t* ptr_methods; /**< pointer to tuple of methods to access to iteration_info of each substep */

        iteration_info( tuple_t& methods, value_t tol = static_cast<value_t>( 0 ) )
            : error( static_cast<value_t>( 0 ) )
            , success( true )
            , is_step( false )
            , number_of_steps( splitting_t::N_steps )
            , number_of_eval( ::detail::init_fill_array<splitting_t::N_methods, std::size_t>( 0 ) )
            , tolerance( tol )
            , ptr_methods( &methods )
        {
        }

        /**
         * @brief get information on substep algorithm for splitting method
         *
         * @tparam I    rank of substep
         * @return auto
         */
        template <std::size_t I>
        auto
        get()
        {
            return std::get<I>( *ptr_methods ).info();
        }

        /**
         * @brief reset number of evaluations to zero
         *
         */
        void
        reset_eval()
        {
            number_of_eval.fill( 0 );
        }
    };

    template <typename user_defined_method_t>
        requires user_defined_method_t::is_user_defined_method
    struct iteration_info<user_defined_method_t>
    {
        /**
         * @brief type of coefficient of Butcher tableau, same type to store error and tolerance
         *
         */
        using value_t = typename user_defined_method_t::value_t;

        value_t error;                /**< error makes on time iteration for adaptive time step method */
        bool success;                 /**< sets as true only for success iteration */
        bool is_step;                 /**< sets as true only if iterator is on a step given in solver */
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
            , success( true )
            , is_step( false )
            , number_of_stages( 0 )
            , number_of_eval( 1 )
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

} // namespace ponio
