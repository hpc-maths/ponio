// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <array>
#include <concepts>
#include <cstddef>
#include <type_traits> // NOLINT(misc-include-cleaner)
#include <utility>

#include "detail.hpp" // NOLINT(misc-include-cleaner)

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
            : A( std::move( A_ ) )
            , b( std::move( b_ ) )
            , c( std::move( c_ ) )
        {
        }

        matrix_t A;
        vector_t b;
        vector_t c;
    };

    template <typename Tableau>
    concept is_butcher_tableau = std::derived_from<Tableau, butcher_tableau<Tableau::N_stages, typename Tableau::value_t>>;

    template <std::size_t N, typename _value_t = double>
    struct adaptive_butcher_tableau : public butcher_tableau<N, _value_t>
    {
        using base_t   = butcher_tableau<N, _value_t>;
        using value_t  = typename base_t::value_t;
        using matrix_t = typename base_t::matrix_t;
        using vector_t = typename base_t::vector_t;

        using base_t::N_stages;

        constexpr adaptive_butcher_tableau( matrix_t&& A_, vector_t&& b1_, vector_t&& b2_, vector_t&& c_ )
            : base_t( std::move( A_ ), std::move( b1_ ), std::move( c_ ) )
            , b2( std::move( b2_ ) )
        {
        }

        vector_t b2;
    };

    template <typename Tableau>
    concept is_embedded_tableau = requires( Tableau t ) { t.b2; };

    template <typename Algorithm_t>
    concept is_embedded = requires( Algorithm_t algo ) {
                              {
                                  std::bool_constant<Algorithm_t::is_embedded>()
                                  } -> std::same_as<std::true_type>;
                          };

    template <typename Tableau>
    concept is_exprk_tableau = requires( Tableau t ) {
                                   typename Tableau::value_t;
                                   typename Tableau::linear_t;
                                   typename Tableau::func_t;
                                   {
                                       t.a
                                   };
                                   {
                                       t.b
                                   };
                                   {
                                       t.c
                                   };
                               };

    template <typename Tableau>
    constexpr bool
    is_explicit_impl()
    {
        if constexpr ( is_exprk_tableau<Tableau> )
        {
            return true;
        }
        else
        {
            using value_t       = typename Tableau::value_t;
            value_t partial_sum = 0;

            for ( std::size_t i = 0ul; i < Tableau::N_stages; ++i )
            {
                for ( std::size_t j = 0ul; j < Tableau::N_stages; ++j )
                {
                    partial_sum += ponio::detail::power<2>( Tableau::A[i][j] );
                }
            }

            return partial_sum == static_cast<value_t>( 0. );
        }
    }

    template <typename Tableau>
    concept is_explicit = requires( Tableau t ) {
                              {
                                  std::bool_constant<is_explicit_impl<Tableau>()>()
                                  } -> std::same_as<std::true_type>;
                          };

    template <typename Tableau>
    constexpr bool
    is_diagonal_implicit_impl()
    {
        if constexpr ( is_exprk_tableau<Tableau> )
        {
            return true;
        }
        else
        {
            using value_t       = typename Tableau::value_t;
            value_t partial_sum = 0;

            for ( std::size_t i = 0ul; i < Tableau::N_stages; ++i )
            {
                for ( std::size_t j = 0ul; j < Tableau::N_stages - 1; ++j )
                {
                    partial_sum += ponio::detail::power<2>( Tableau::A[i][j] );
                }
            }

            return partial_sum == static_cast<value_t>( 0. );
        }
    }

    template <typename Tableau>
    concept is_diagonal_implicit = requires( Tableau t ) {
                                       {
                                           std::bool_constant<is_diagonal_implicit_impl<Tableau>()>()
                                           } -> std::same_as<std::true_type>;
                                   };

    template <typename tableau_im_t, typename tableau_ex_t>
    concept is_compatible_tableaus = requires( tableau_im_t t_im, tableau_ex_t t_ex ) {
                                         // check if arguments are Butcher tableaus
                                         requires is_butcher_tableau<tableau_im_t>;
                                         requires is_butcher_tableau<tableau_ex_t>;

                                         // check size of each tableau
                                         {
                                             std::bool_constant<tableau_im_t::N_stages == tableau_ex_t::N_stages>()
                                             } -> std::same_as<std::true_type>;

                                         // check type stored in tableau
                                         requires std::same_as<typename tableau_im_t::value_t, typename tableau_ex_t::value_t>;

                                         // check compatibility in ponio
                                         requires is_diagonal_implicit<tableau_ex_t>;
                                         requires is_explicit<tableau_ex_t>;
                                     };

    template <typename _tableau_im_t, typename _tableau_ex_t>
        requires is_compatible_tableaus<_tableau_im_t, _tableau_ex_t>
    struct pair_butcher_tableau
    {
        using tableau_im_t = _tableau_im_t;
        using tableau_ex_t = _tableau_ex_t;

        static constexpr std::size_t N_stages    = tableau_im_t::N_stages;
        static constexpr bool is_embedded        = is_embedded_tableau<tableau_im_t> && is_embedded_tableau<tableau_ex_t>;
        static constexpr bool is_imex_method     = true;
        static constexpr std::size_t N_operators = 2;

        using value_t  = tableau_im_t::value_t;
        using matrix_t = tableau_im_t::matrix_t;
        using vector_t = tableau_im_t::vector_t;

        pair_butcher_tableau() = default;

        tableau_im_t tableau_im;
        tableau_ex_t tableau_ex;
    };

} // namespace ponio::runge_kutta::butcher
