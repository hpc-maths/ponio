// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <cstddef>

#include "butcher_tableau.hpp"
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
        bool success;                 /**< set as true only for success iteration */
        bool is_step;                 /**< set as true only if iterator is on a step given in solver */
        std::size_t number_of_stages; /**< number of stages of method */
        std::size_t number_of_eval;   /**< number of evaluation of function */
        value_t tolerance;            /** tolerance for the method (for adaptive time step method) */

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
    };

} // namespace ponio
