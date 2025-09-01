// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <algorithm>
#include <concepts>
#include <limits>
#include <type_traits>

namespace ponio::expression
{
    namespace detail
    {
        /**
         * @brief test if template element is a member of a ponio expression
         *
         * @tparam expression_t
         */
        template <typename expression_t>
        concept is_ponio_expression = std::same_as<std::true_type, std::integral_constant<bool, expression_t::is_ponio_expression>>;
    }

    /**
     * @brief one of leaf of ponio expression that store a container
     *
     * @tparam container_t
     */
    template <typename container_t>
    struct state
    {
        static constexpr bool is_ponio_expression = true;

        using container_type = container_t;
        container_type& _data; // store only a reference on data

        /**
         * @brief Construct a new state object from the reference on a container
         *
         * @param dat container to provide expressions
         */
        state( container_type& dat )
            : _data( dat )
        {
        }

        /**
         * @brief accessor operator (const and non-const version)
         *
         * @param i index to access
         */
        auto const&
        operator[]( std::size_t i ) const
        {
            return _data[i];
        }

        auto&
        operator[]( std::size_t i )
        {
            return _data[i];
        }

        /**
         * @brief Compute expression only here in the loop
         *
         * @tparam expression_t type of expression tree
         * @param expr          expression
         */
        template <typename expression_t>
            requires detail::is_ponio_expression<expression_t>
        state&
        operator=( expression_t const& expr )
        {
            for ( std::size_t i = 0ul; i < expr.size(); ++i )
            {
                _data[i] = expr[i];
            }

            return *this;
        }

        /**
         * @brief returns the size of container
         */
        std::size_t
        size() const
        {
            return _data.size();
        }

        /**
         * @brief raw data accessor (const and non-const version)
         */
        container_type const&
        data() const
        {
            return _data;
        }

        container_type&
        data()
        {
            return _data;
        }

        void
        array_reader()
        {
        }

        void
        array_restore()
        {
        }
    };

    /**
     * @brief helper function to construct a `state` object from a container
     *
     * @tparam container_t type of container
     * @param lhs          container on which we provide expressions
     */
    template <typename container_t>
    auto
    make_state( container_t& lhs )
    {
        return state<container_t>( lhs );
    }

    /**
     * @brief one leaf of ponio expression that store a scalar
     *
     * @tparam value_t type of scalar
     */
    template <typename value_t>
        requires( std::floating_point<value_t> || std::integral<value_t> )
    struct state<value_t>
    {
        static constexpr bool is_ponio_expression = true;

        using container_type = value_t;
        container_type _value;

        state( container_type const& dat )
            : _value( dat )
        {
        }

        auto const&
        operator[]( std::size_t ) const
        {
            return _value;
        }

        std::size_t
        size() const
        {
            // un hack pour itérer à l'infini sur cette valeur
            return std::numeric_limits<std::size_t>::max();
        }

        container_type const&
        data() const
        {
            return _value;
        }

        container_type&
        data()
        {
            return _value;
        }

        void
        array_reader()
        {
        }

        void
        array_restore()
        {
        }
    };

    template <typename value_t>
        requires( std::floating_point<value_t> || std::integral<value_t> )
    auto make_scalar( value_t const& lhs )
    {
        return state<value_t>( lhs );
    }

    /////////////////////////////////////////////////////////////////
    // UNARY OPERATOR
    /////////////////////////////////////////////////////////////////

    enum struct unary_operation
    {
        plus,
        minus
    };

    /**
     * @brief concept to test equality between element of `unary_operation` enumeration
     *
     * @tparam I
     * @tparam J
     */
    template <unary_operation I, unary_operation J>
    concept equals_unary_op = std::same_as<std::true_type, std::integral_constant<bool, I == J>>;

    /**
     * @brief base function for unary operators
     *
     * @tparam op    type of unary operation
     * @tparam lhs_t type of operand
     */
    template <unary_operation op, typename lhs_t>
    auto
    operation( lhs_t const&, std::size_t )
    {
    }

    /**
     * @brief function to compute \f$+\cdot\f$
     *
     * @tparam op    type of unary operation (here \f$+\f$)
     * @tparam lhs_t type of operand
     * @param lhs    operand
     * @param i      index to compute this expression
     */
    template <unary_operation op, typename lhs_t>
        requires equals_unary_op<op, unary_operation::plus>
    auto
    operation( lhs_t const& lhs, std::size_t i )
    {
        return +lhs[i];
    }

    /**
     * @brief function to compute \f$-\cdot\f$
     *
     * @tparam op    type of unary operation (here \f$-\f$)
     * @tparam lhs_t type of operand
     * @param lhs    operand
     * @param i      index to compute this expression
     */
    template <unary_operation op, typename lhs_t>
        requires equals_unary_op<op, unary_operation::minus>
    auto
    operation( lhs_t const& lhs, std::size_t i )
    {
        return -lhs[i];
    }

    /**
     * @brief base class to compute unary operation
     *
     * @tparam op    type of unary operation
     * @tparam lhs_t type of operand
     */
    template <unary_operation op, typename lhs_t>
        requires detail::is_ponio_expression<lhs_t>
    struct unary_op
    {
        static constexpr bool is_ponio_expression = true;

        lhs_t lhs; // should copy the part of expression tree or keep it by reference on temporary object?

        unary_op( lhs_t const& l )
            : lhs( l )
        {
        }

        auto
        operator[]( std::size_t i ) const
        {
            return operation<op>( lhs, i );
        }

        std::size_t
        size() const
        {
            return lhs.size();
        }

        void
        array_reader()
        {
            lhs.array_reader();
        }

        void
        array_restore()
        {
            lhs.array_restore();
        }
    };

    /**
     * @brief helper function to construct a new object `unary_op` for plus
     *
     * @tparam lhs_t type of operand
     * @param lhs    operand
     */
    template <typename lhs_t>
        requires detail::is_ponio_expression<lhs_t>
    unary_op<unary_operation::plus, lhs_t>
    operator+( lhs_t const& lhs )
    {
        return unary_op<unary_operation::plus, lhs_t>( lhs );
    }

    /**
     * @brief helper function to construct a new object `unary_op` for minus
     *
     * @tparam lhs_t type of operand
     * @param lhs    operand
     */
    template <typename lhs_t>
        requires detail::is_ponio_expression<lhs_t>
    unary_op<unary_operation::minus, lhs_t>
    operator-( lhs_t const& lhs )
    {
        return unary_op<unary_operation::minus, lhs_t>( lhs );
    }

    /////////////////////////////////////////////////////////////////
    // BINARY OPERATOR
    /////////////////////////////////////////////////////////////////

    enum struct binary_operation
    {
        add,
        sub,
        mul,
        div
    };

    /**
     * @brief concept to test equality between element of `binary_operation` enumeration
     *
     * @tparam I
     * @tparam J
     */
    template <binary_operation I, binary_operation J>
    concept equals_binary_op = std::same_as<std::true_type, std::integral_constant<bool, I == J>>;

    /**
     * @brief base function for binary operators
     *
     * @tparam op    type of binary operation
     * @tparam lhs_t type of operand
     */
    template <binary_operation op, typename lhs_t, typename rhs_t>
    auto
    operation( lhs_t const&, rhs_t const&, std::size_t )
    {
    }

    /**
     * @brief function to compute \f$\cdot + \cdot\f$
     *
     * @tparam op    type of binary operation (here \f$+\f$)
     * @tparam lhs_t type of left-hand-side
     * @tparam rhs_t type of right-hand-side
     * @param lhs    lhs operand
     * @param rhs    rhs operand
     * @param i      index to compute expression
     */
    template <binary_operation op, typename lhs_t, typename rhs_t>
        requires equals_binary_op<op, binary_operation::add>
    auto
    operation( lhs_t const& lhs, rhs_t const& rhs, std::size_t i )
    {
        return lhs[i] + rhs[i];
    }

    /**
     * @brief function to compute \f$\cdot - \cdot\f$
     *
     * @tparam op    type of binary operation (here \f$-\f$)
     * @tparam lhs_t type of left-hand-side
     * @tparam rhs_t type of right-hand-side
     * @param lhs    lhs operand
     * @param rhs    rhs operand
     * @param i      index to compute expression
     */
    template <binary_operation op, typename lhs_t, typename rhs_t>
        requires equals_binary_op<op, binary_operation::sub>
    auto
    operation( lhs_t const& lhs, rhs_t const& rhs, std::size_t i )
    {
        return lhs[i] - rhs[i];
    }

    /**
     * @brief function to compute \f$\cdot \times \cdot\f$
     *
     * @tparam op    type of binary operation (here \f$\times\f$)
     * @tparam lhs_t type of left-hand-side
     * @tparam rhs_t type of right-hand-side
     * @param lhs    lhs operand
     * @param rhs    rhs operand
     * @param i      index to compute expression
     */
    template <binary_operation op, typename lhs_t, typename rhs_t>
        requires equals_binary_op<op, binary_operation::mul>
    auto
    operation( lhs_t const& lhs, rhs_t const& rhs, std::size_t i )
    {
        return lhs[i] * rhs[i];
    }

    /**
     * @brief function to compute \f$\cdot / \cdot\f$
     *
     * @tparam op    type of binary operation (here \f$/\f$)
     * @tparam lhs_t type of left-hand-side
     * @tparam rhs_t type of right-hand-side
     * @param lhs    lhs operand
     * @param rhs    rhs operand
     * @param i      index to compute expression
     */
    template <binary_operation op, typename lhs_t, typename rhs_t>
        requires equals_binary_op<op, binary_operation::div>
    auto
    operation( lhs_t const& lhs, rhs_t const& rhs, std::size_t i )
    {
        return lhs[i] / rhs[i];
    }

    /**
     * @brief base class to compute binary operation
     *
     * @tparam op    type of operation
     * @tparam lhs_t type of lhs
     * @tparam rhs_t type of rhs
     */
    template <binary_operation op, typename lhs_t, typename rhs_t>
        requires( detail::is_ponio_expression<lhs_t> && detail::is_ponio_expression<rhs_t> )
    struct binary_op
    {
        static constexpr bool is_ponio_expression = true;

        lhs_t lhs; // should copy the part of expression tree or keep it by reference on temporary object?
        rhs_t rhs; // should copy the part of expression tree or keep it by reference on temporary object?

        binary_op( lhs_t const& l, rhs_t const& r )
            : lhs( l )
            , rhs( r )
        {
        }

        auto
        operator[]( std::size_t i ) const
        {
            return operation<op>( lhs, rhs, i );
        }

        std::size_t
        size() const
        {
            return std::min( lhs.size(), rhs.size() );
        }

        void
        array_reader()
        {
            lhs.array_reader();
            rhs.array_reader();
        }

        void
        array_restore()
        {
            lhs.array_restore();
            rhs.array_restore();
        }
    };

    /**
     * @brief helper function to construct a new object `binary_op` for addition
     *
     * @tparam lhs_t lhs type
     * @tparam rhs_t rhs type
     */
    template <typename lhs_t, typename rhs_t>
        requires( detail::is_ponio_expression<lhs_t> && detail::is_ponio_expression<rhs_t> )
    binary_op<binary_operation::add, lhs_t, rhs_t>
    operator+( lhs_t const& lhs, rhs_t const& rhs )
    {
        return binary_op<binary_operation::add, lhs_t, rhs_t>( lhs, rhs );
    }

    /**
     * @brief helper function to construct a new object `binary_op` for substraction
     *
     * @tparam lhs_t lhs type
     * @tparam rhs_t rhs type
     */
    template <typename lhs_t, typename rhs_t>
        requires( detail::is_ponio_expression<lhs_t> && detail::is_ponio_expression<rhs_t> )
    binary_op<binary_operation::sub, lhs_t, rhs_t>
    operator-( lhs_t const& lhs, rhs_t const& rhs )
    {
        return binary_op<binary_operation::sub, lhs_t, rhs_t>( lhs, rhs );
    }

    /**
     * @brief helper function to construct a new object `binary_op` for multiplication
     *
     * @tparam lhs_t lhs type
     * @tparam rhs_t rhs type
     */
    template <typename lhs_t, typename rhs_t>
        requires( detail::is_ponio_expression<lhs_t> && detail::is_ponio_expression<rhs_t> )
    binary_op<binary_operation::mul, lhs_t, rhs_t>
    operator*( lhs_t const& lhs, rhs_t const& rhs )
    {
        return binary_op<binary_operation::mul, lhs_t, rhs_t>( lhs, rhs );
    }

    /**
     * @brief helper function to construct a new object `binary_op` for division
     *
     * @tparam lhs_t lhs type
     * @tparam rhs_t rhs type
     */
    template <typename lhs_t, typename rhs_t>
        requires( detail::is_ponio_expression<lhs_t> && detail::is_ponio_expression<rhs_t> )
    binary_op<binary_operation::div, lhs_t, rhs_t>
    operator/( lhs_t const& lhs, rhs_t const& rhs )
    {
        return binary_op<binary_operation::div, lhs_t, rhs_t>( lhs, rhs );
    }
} // namespace ponio::expression
