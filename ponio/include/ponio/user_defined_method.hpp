// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <concepts>
#include <limits>
#include <type_traits>

#include "iteration_info.hpp"

namespace ponio
{

    template <typename u_method_t>
    concept is_user_defined_method = requires( u_method_t u_method ) {
                                         {
                                             std::bool_constant<u_method_t::is_user_defined_method>()
                                             } -> std::same_as<std::true_type>;
                                     };

    /** @class user_defined_method
     * @brief helper class to provide user having its solver
     *
     * @tparam _user_function_t type of user function
     */
    template <typename _user_function_t>
    struct user_defined_method
    {
        using user_function_t = _user_function_t;
        user_function_t* user_function;

        user_defined_method( user_function_t& func )
            : user_function( &func )
        {
        }
    };

    template <typename _value_t, typename _user_function_t>
    struct user_defined_algorithm
    {
        static constexpr bool is_user_defined_method = true;
        static constexpr std::size_t N_stages        = 0;
        static constexpr bool is_embedded            = false;
        static constexpr std::size_t order           = std::numeric_limits<std::size_t>::infinity();
        static constexpr std::string_view id         = "user_defined";

        using value_t         = _value_t;
        using user_function_t = _user_function_t;

        iteration_info<user_defined_algorithm> info;

        user_function_t* user_function;

        user_defined_algorithm() = default;

        user_defined_algorithm( user_defined_method<user_function_t>&& u_method )
            : user_function( u_method.user_function )
        {
        }

        template <typename problem_t, typename state_t, typename array_k0_t>
        std::tuple<value_t, state_t, value_t>
        operator()( problem_t& f, value_t& tn, state_t& un, array_k0_t&, value_t& dt )
        {
            return ( *user_function )( f, tn, un, dt );
        }
    };

    template <typename value_t, typename user_defined_method_t>
    auto
    make_user_defined_algorithm( user_defined_method_t&& u_meth )
    {
        return user_defined_algorithm<value_t, typename user_defined_method_t::user_function_t>( u_meth );
    }

} // namespace ponio
