// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <cstddef>

#include "butcher_tableau.hpp"
#include "stage.hpp"

namespace ponio
{

    template <typename tableau_t>
    struct iteration_info
    {
        using value_t = typename tableau_t::value_t;

        value_t error;
        bool success;
        bool is_step;
        std::size_t number_of_stages;
        std::size_t number_of_eval;
        value_t tolerance;

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
