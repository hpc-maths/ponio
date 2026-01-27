// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <cstddef>
#include <limits>
#include <type_traits> // NOLINT(misc-include-cleaner)

namespace ponio // NOLINT(modernize-concat-nested-namespaces)
{

    /** @class Stage
     *  @brief indicator for overloading each stages of multistages method
     */
    template <std::size_t s>
    using Stage = std::integral_constant<std::size_t, s>;

    namespace stages
    {
        /**
         * @brief constant to set dynamic number of stages in an algorithm
         * @warning don't check if `N_stages` is lower than `ponio::stages::dynamic` because `std::numeric_limits<std::size_t>::infinity()
         == 0`

         */
        static constexpr std::size_t dynamic = std::numeric_limits<std::size_t>::infinity();

        /**
         * @brief test is algorithm has a static number of stages
         *
         * @tparam Algorithm_t algorithm (Runge-Kutta method) to check
         */
        template <typename Algorithm_t>
        concept has_static_number_of_stages = requires( Algorithm_t alg ) {
                                                  // clang-format off
            // NOLINTBEGIN(clang-diagnostic-error)
            { Algorithm_t::N_stages } -> std::convertible_to<std::size_t>;
            { std::bool_constant<Algorithm_t::N_stages != dynamic>() } -> std::same_as<std::true_type>;
            // NOLINTEND(clang-diagnostic-error)
                                                  // clang-format on
                                              };

        /**
         * @brief test is algorithm has a dynamic number of stages
         *
         * @tparam Algorithm_t algorithm (Runge-Kutta method) to check
         */
        template <typename Algorithm_t>
        concept has_dynamic_number_of_stages = requires( Algorithm_t alg ) {
                                                   // clang-format off
            // NOLINTBEGIN(clang-diagnostic-error)
            { Algorithm_t::N_storage } -> std::convertible_to<std::size_t>;
            { std::bool_constant<!has_static_number_of_stages<Algorithm_t>>() } -> std::same_as<std::true_type>;
            // NOLINTEND(clang-diagnostic-error)
                                                   // clang-format on
                                               };

    } // namespace stages

    /** @class sub_method
     *  @brief indicator for sub-method in splitting method
     */
    template <std::size_t I>
    using sub_method = std::integral_constant<std::size_t, I>;

} // namespace ponio
