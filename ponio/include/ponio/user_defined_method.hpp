// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <concepts> // NOLINT(misc-include-cleaner)
#include <cstddef>
#include <limits>
#include <string_view>
#include <tuple>
#include <type_traits> // NOLINT(misc-include-cleaner)

#include "iteration_info.hpp"

namespace ponio
{

    /**
     * @brief test if method or algorithm is defined by user by testing `is_user_defined_method` static attribut
     */
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
        static constexpr bool is_user_defined_method = true;

        using user_function_t = _user_function_t;
        user_function_t* user_function;

        user_defined_method( user_function_t& func )
            : user_function( &func )
        {
        }
    };

    /**
     * @brief helper function to build a user_defined_method
     *
     * @tparam user_function_t type of function (lambda, function or functor)
     * @param algo     function that takes $(f, t^n, u^n, \Delta t)$ (could be by reference) and returns a tuple \f$(t^{n+1}, u^{n+1},
     * \Delta t)\f$
     */
    template <typename user_function_t>
    auto
    make_user_defined_method( user_function_t& algo )
    {
        return user_defined_method<user_function_t>( algo );
    }

    /** @class user_defined_algorithm
     * @brief algorithm in term of ponio with the function provides by user
     *
     * @tparam _value_t         type of coefficients and values in provides info
     * @tparam _user_function_t type of user_defined_method with the underlying function provides by user
     */
    template <typename _value_t, typename _user_function_t>
    struct user_defined_algorithm
    {
        static constexpr bool is_user_defined_method = true;
        static constexpr bool is_embedded            = false;
        static constexpr std::size_t order           = std::numeric_limits<std::size_t>::infinity();
        static constexpr std::string_view id         = "user_defined";

        using value_t         = _value_t;
        using user_function_t = _user_function_t;

        iteration_info<user_defined_algorithm> _info;

        user_function_t* user_function = nullptr;

        user_defined_algorithm() = default;

        user_defined_algorithm( user_defined_method<user_function_t> const& u_method )
            : user_function( u_method.user_function )
        {
        }

        template <typename problem_t, typename state_t>
        std::tuple<value_t, state_t, value_t>
        operator()( problem_t& f, value_t& tn, state_t& un, value_t& dt )
        {
            return ( *user_function )( f, tn, un, dt );
        }

        /**
         * @brief gets `iteration_info` object
         */
        auto&
        info()
        {
            return _info;
        }

        /**
         * @brief gets `iteration_info` object (constant version)
         */
        auto const&
        info() const
        {
            return _info;
        }
    };

    /**
     * @brief helper function to build a `user_defined_algorithm`
     *
     * @tparam value_t               type of coefficients
     * @tparam user_defined_method_t type of `user_defined_method` with the underlying user function
     * @param u_meth `user_defined_method` with the underlying user function
     */
    template <typename value_t, typename user_defined_method_t>
    auto
    make_user_defined_algorithm( user_defined_method_t const& u_meth )
    {
        return user_defined_algorithm<value_t, typename user_defined_method_t::user_function_t>( u_meth );
    }

} // namespace ponio
