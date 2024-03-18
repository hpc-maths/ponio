// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <array>
#include <concepts>
#include <cstddef>
#include <utility>

namespace ponio::runge_kutta::butcher
{

    template <std::size_t N, typename _value_t = double>
    struct butcher_tableau
    {
        static constexpr std::size_t N_stages = N;

        using value_t  = _value_t;
        using matrix_t = std::array<std::array<value_t, N_stages>, N_stages>;
        using vector_t = std::array<value_t, N_stages>;

        constexpr butcher_tableau( matrix_t&& A_, vector_t&& b_, vector_t&& c_ )
            : A( std::forward<matrix_t>( A_ ) )
            , b( std::forward<vector_t>( b_ ) )
            , c( std::forward<vector_t>( c_ ) )
        {
        }

        matrix_t A;
        vector_t b;
        vector_t c;
    };

    template <typename Tableau>
    concept is_butcher_tableau = std::derived_from<butcher_tableau<Tableau::N_stages, typename Tableau::value_t>, Tableau>;

    template <std::size_t N, typename _value_t = double>
    struct adaptive_butcher_tableau : public butcher_tableau<N, _value_t>
    {
        using base_t   = butcher_tableau<N, _value_t>;
        using value_t  = typename base_t::value_t;
        using matrix_t = typename base_t::matrix_t;
        using vector_t = typename base_t::vector_t;

        using base_t::N_stages;

        constexpr adaptive_butcher_tableau( matrix_t&& A_, vector_t&& b1_, vector_t&& b2_, vector_t&& c_ )
            : base_t( std::forward<matrix_t>( A_ ), std::forward<vector_t>( b1_ ), std::forward<vector_t>( c_ ) )
            , b2( std::forward<vector_t>( b2_ ) )
        {
        }

        vector_t b2;
    };

    template <typename Tableau>
    concept is_embedded_tableau = requires( Tableau t ) { t.b2; };

    template <typename Algorithm_t>
    concept is_embedded = requires { Algorithm_t::is_embedded; };

} // namespace ponio::butcher
