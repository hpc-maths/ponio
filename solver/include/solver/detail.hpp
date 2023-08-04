// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <array>
#include <concepts>
#include <functional>
#include <ranges>
#include <utility>

namespace detail
{

    /* tpl_inner_product */
    template <typename state_t, typename value_t, typename ArrayA_t, typename ArrayB_t, std::size_t... Is>
    constexpr state_t
    tpl_inner_product_impl( ArrayA_t const& a,
        ArrayB_t const& b,
        state_t const& init,
        [[maybe_unused]] value_t mul_coeff,
        std::index_sequence<Is...> )
    {
        return ( init + ... + ( mul_coeff * a[Is] * b[Is] ) );
    }

    /**
     * @brief inner product between two array from 0 to N
     *
     * @tparam N        number of elements to compute
     * @tparam state_t  type of computed value
     * @tparam value_t  type of coefficents
     * @tparam ArrayA_t type of first array
     * @tparam ArrayB_t type of second array
     *
     * @param a         first array
     * @param b         second array
     * @param init      starting value to add other values to
     * @param mul_coeff coefficient to multiply each multiplication of inner product
     *
     * @details This function compute \f$\texttt{init} + \sum_{i=0}^N \texttt{mul_coeff}a_ib_i\f$ without loop thanks to template.
     */
    template <std::size_t N, typename state_t, typename value_t, typename ArrayA_t, typename ArrayB_t>
    constexpr state_t
    tpl_inner_product( ArrayA_t const& a, ArrayB_t const& b, state_t const& init, value_t mul_coeff = value_t{ 1.0 } )
    {
        return tpl_inner_product_impl( a, b, init, mul_coeff, std::make_index_sequence<N>() );
    }

    /* init_fill_array */
    // first version with a value
    template <typename T, std::size_t... Is>
    constexpr std::array<std::remove_cvref_t<T>, sizeof...( Is )>
    init_fill_array_impl( T&& value, std::index_sequence<Is...> )
    {
        return { { ( static_cast<void>( Is ), value )... } };
    }

    /**
     * @brief fill an uninitialize array
     *
     * @tparam N size of array to fill
     * @tparam T type of stored value in array
     *
     * @param value value to fill in uninitialize array
     *
     * @code{,cpp}
     *   int i = 42;
     *   const std::array<int,8> arr = detail::init_fill_array<8>( i ); // all values of `arr` are `42`
     * @endcode
     *
     * @details `value` can also be an invokable object that take an unsigned integer and return type stored in array `T`.
     *
     * @code{,cpp}
     *   const std::array<int,8> arr = detail::init_fill_array<8>([](int i){ return i*i; }); // get {0,1,4,9,16,25,36,49}
     * @endcode
     *
     */
    template <std::size_t N, typename T>
    constexpr std::array<std::remove_cvref_t<T>, N>
    init_fill_array( T&& value )
    {
        return init_fill_array_impl( std::forward<T>( value ), std::make_index_sequence<N>() );
    }

    /* applyable_impl */
    // unpack tuple to test if `function_t` is invocable with types in tuple
    template <typename function_t, typename tuple_args_t, std::size_t... Is>
    constexpr bool
    applyable_impl( std::index_sequence<Is...> )
    {
        return std::invocable<function_t, std::tuple_element_t<Is, tuple_args_t>...>;
    }

    /**
     * @brief concept that specifies that a callable type `function_t` can be called with the unpacked `tuple_args_t` tuple
     *
     * @tparam function_t   type of callable
     * @tparam tuple_args_t type of tuple to unpack
     */
    template <typename function_t, typename tuple_args_t>
    concept applyable = applyable_impl<function_t, tuple_args_t>( std::make_index_sequence<std::tuple_size<tuple_args_t>::value>{} );

#ifndef IN_DOXYGEN
    // second version with a invocable parameter and its arguments
    template <typename Return_t, typename Function_t, typename TupleArgs_t, std::size_t... Is>
        requires applyable<Function_t, TupleArgs_t>
    constexpr std::array<Return_t, sizeof...( Is )>
    init_fill_array_impl( Function_t&& f, TupleArgs_t&& args, std::index_sequence<Is...> )
    {
        return { { ( static_cast<void>( Is ), std::apply( f, args ) )... } };
    }

    /**
     * @brief fill an uninitialize array with a generator
     *
     * @tparam N size of array to fill
     * @tparam Function_t type of invocable object
     * @tparam Args       types of arguments to call Function_t
     *
     * @param f invocable object that need to get a `Args...` and return the value type of output array
     * @param args arguments of function `f`
     *
     * @code{,cpp}
     *   const std::array<int,8> arr = detail::init_fill_array<3>([](int i){ return i*i; }, 2); // get {4,4,4}
     * @endcode
     */
    template <std::size_t N, typename Function_t, typename... Args>
        requires std::invocable<Function_t, Args...>
    constexpr std::array<typename std::invoke_result<Function_t, Args...>::type, N>
    init_fill_array( Function_t&& f, Args&&... args )
    {
        return init_fill_array_impl<typename std::invoke_result<Function_t, Args...>::type>( std::forward<Function_t>( f ),
            std::forward_as_tuple( args... ),
            std::make_index_sequence<N>() );
    }
#endif

    template <typename Arithmetic, std::size_t... Is>
    constexpr Arithmetic
    power_impl( Arithmetic&& value, std::index_sequence<Is...> )
    {
        return ( static_cast<Arithmetic>( 1.0 ) * ... * ( static_cast<void>( Is ), value ) );
    }

    /**
     * @brief return `Iexp` multiplication of `value` to make an efficient power function for integers
     *
     * @tparam Iexp exponent of expression
     * @tparam Arithmetic_t arithmetic type (need multiplication and `1` can be converted into this type)
     * @param value value to powered
     * @return constexpr Arithmetic_t \f$\texttt{value}^{\texttt{Iexp}}=1\times\underbrace{\texttt{value}\times\cdots\times
     * \texttt{value}}_{\texttt{Iexp}\ \text{times}}\f$ with \f$\texttt{Iexp}\in\mathbb{N}\f$ (could be equal to zero)
     *
     * @note if we remove condition \p Iexp could be equal to zero, we could remove the requirement `1` can be converted into
     * `Arithmetic_t`.
     */
    template <std::size_t Iexp, typename Arithmetic_t>
    constexpr Arithmetic_t
    power( Arithmetic_t&& value )
    {
        return power_impl( std::forward<Arithmetic_t>( value ), std::make_index_sequence<Iexp>() );
    }

} // namespace detail
