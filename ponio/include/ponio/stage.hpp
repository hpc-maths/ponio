// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <cstddef>
#include <limits>
#include <type_traits>

namespace ponio
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
        concept has_static_number_of_stages = Algorithm_t::N_stages !=
        dynamic;

        /**
         * @brief test is algorithm has a dynamic number of stages
         *
         * @tparam Algorithm_t algorithm (Runge-Kutta method) to check
         */
        template <typename Algorithm_t>
        concept has_dynamic_number_of_stages = !
        has_static_number_of_stages<Algorithm_t>;

    }

}
