// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <cmath>
#include <concepts>
#include <cstddef>

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

        static constexpr matrix_type
        identity_minus( matrix_type const& matrix, scalar_t alpha )
        {
            return static_cast<matrix_type>( 1.0 ) - alpha * matrix;
        }

        static vector_type
        solver( matrix_type const& dfx, vector_type const& fx )
        {
            return fx / dfx;
        }
    };

    template <typename state_t>
    struct operator_algebra
    {
        template <typename>
        static constexpr bool dependent_false = false;

        static_assert( dependent_false<state_t>, "non implemented operator_algebra structure for this type" );
    };

    template <typename scalar_t>
        requires std::floating_point<scalar_t>
    struct operator_algebra<scalar_t>
    {
        static constexpr auto
        identity( scalar_t const& )
        {
            return static_cast<scalar_t>( 1.0 );
        }

        template <typename operator_t, typename rhs_t>
        static void
        solve( operator_t& op, scalar_t& xn, rhs_t const& rhs )
        {
            static constexpr std::size_t max_iter = 100;
            static constexpr scalar_t tol         = 1e-5;

            auto xnm1    = xn * ( 1 - 1e-4 );
            auto op_xn   = op( xn );
            auto op_xnm1 = op( xnm1 );

            scalar_t residual = abs( op_xn - rhs );

            std::size_t iter = 0;
            while ( iter < max_iter && residual > tol )
            {
                auto h = xn - xnm1;

                auto xnp1 = xn - ( op_xn - rhs ) / ( ( op_xn - op_xnm1 ) / h );

                op_xnm1 = op_xn;
                xnm1    = xn;
                xn      = xnp1;
                op_xn   = op( xn );

                residual = abs( op_xn - rhs );
                iter += 1;
            }
        }
    };

    /**
     * @brief in-place update of a state, as needed by compensated summation
     *
     * The default implementation uses compound assignment, which is provided by
     * arithmetic types and by Eigen vectors. A state type that offers addition
     * and assignment but no compound assignment specializes this structure.
     *
     * @tparam state_t type of state
     */
    template <typename state_t>
    struct state_algebra
    {
        /**
         * @brief computes \f$y \leftarrow \alpha y\f$
         *
         * @param y     state to scale
         * @param alpha scaling factor
         */
        template <typename value_t>
        static void
        scale( state_t& y, value_t alpha )
        {
            y *= alpha;
        }

        /**
         * @brief computes \f$y \leftarrow y + x\f$
         *
         * @param y state to update
         * @param x increment, possibly an unevaluated expression
         */
        template <typename increment_t>
        static void
        add( state_t& y, increment_t const& x )
        {
            y += x;
        }
    };

} // namespace ponio::linear_algebra

namespace ponio::shampine_trick
{
    template <typename state_t>
    struct shampine_trick
    {
        template <typename>
        static constexpr bool dependent_false = false;

        static_assert( dependent_false<state_t>, "non implemented Shampine's trick structure for this type" );
    };
} // namespace ponio::shampine_trick
