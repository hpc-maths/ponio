// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <concepts>

namespace ponio::linear_algebra
{

    template <typename matrix_t>
    struct linear_algebra
    {
        template <typename>
        static constexpr bool dependent_false = false;

        static_assert( dependent_false<matrix_t>, "non implemented linear_algebra structure for this type" );
    };

    template <typename scalar_t>
        requires std::floating_point<scalar_t>
    struct linear_algebra<scalar_t>
    {
        using matrix_type = scalar_t;
        using vector_type = scalar_t;

        static constexpr matrix_type
        identity( vector_type const& )
        {
            return static_cast<matrix_type>( 1.0 );
        }

        static vector_type
        solver( matrix_type const& dfx, vector_type const& fx )
        {
            return fx / dfx;
        }
    };

} // namespace ponio::linear_algebra
