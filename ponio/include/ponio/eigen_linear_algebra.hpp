// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

// NOLINTBEGIN(misc-include-cleaner)

#if __has_include( <eigen3/Eigen/Dense>)
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#elif __has_include( <Eigen/Dense>)
#include <Eigen/Dense>
#include <Eigen/Sparse>
#else
#error "Eigen should be included"
#endif

// NOLINTEND(misc-include-cleaner)

#include "linear_algebra.hpp"

namespace ponio::linear_algebra
{

    template <typename scalar_t, int size, int options, int maxrows, int maxcols>
    struct linear_algebra<Eigen::Matrix<scalar_t, size, size, options, maxrows, maxcols>> // NOLINT(misc-include-cleaner)
    {
        using matrix_type = Eigen::Matrix<scalar_t, size, size>; // NOLINT(misc-include-cleaner)
        using vector_type = Eigen::Vector<scalar_t, size>;       // NOLINT(misc-include-cleaner)

        static matrix_type
        identity( vector_type const& )
        {
            return matrix_type::Identity();
        }

        static vector_type
        solver( matrix_type const& dfx, vector_type const& fx )
        {
            return dfx.householderQr().solve( fx );
        }
    };

    template <typename scalar_t>
    struct linear_algebra<Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>> // NOLINT(misc-include-cleaner)
    {
        using matrix_type = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>; // NOLINT(misc-include-cleaner)
        using vector_type = Eigen::Vector<scalar_t, Eigen::Dynamic>;                 // NOLINT(misc-include-cleaner)

        static matrix_type
        identity( vector_type const& u )
        {
            return matrix_type::Identity( u.rows(), u.rows() );
        }

        static vector_type
        solver( matrix_type const& dfx, vector_type const& fx )
        {
            return dfx.colPivHouseholderQr().solve( fx );
        }
    };

    template <typename scalar_t>
    struct linear_algebra<Eigen::SparseMatrix<scalar_t>> // NOLINT(misc-include-cleaner)
    {
        using matrix_type = Eigen::SparseMatrix<scalar_t>;           // NOLINT(misc-include-cleaner)
        using vector_type = Eigen::Vector<scalar_t, Eigen::Dynamic>; // NOLINT(misc-include-cleaner)
        using solver_type = Eigen::SimplicialCholesky<matrix_type>;  // NOLINT(misc-include-cleaner)

      private:

        inline static matrix_type I; // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)

      public:

        static matrix_type const&
        identity( vector_type const& u )
        {
            I.resize( u.size(), u.size() );
            I.setIdentity();

            return I;
        }

        static vector_type
        solver( matrix_type const& dfx, vector_type const& fx )
        {
            solver_type chol( dfx );
            return chol.solve( fx );
        }
    };

} // namespace ponio::linear_algebra

namespace ponio::shampine_trick
{
    template <typename scalar_t, int size, int options, int maxrows>
    struct shampine_trick<Eigen::Matrix<scalar_t, size, 1, options, maxrows, 1>>
    {
        using value_t = scalar_t;

        template <std::size_t ell, typename jacobian_reac_t, typename state_t>
        void
        operator()( value_t alpha, jacobian_reac_t&& jacobian_reac, state_t& initial_guess, state_t& rhs, state_t& u_tmp, state_t& shampine_result )
        {
            auto J_reac = jacobian_reac( initial_guess );

            using matrix_t = std::decay_t<decltype( J_reac )>;

            matrix_t J_R = ::ponio::linear_algebra::linear_algebra<matrix_t>::identity( initial_guess ) - alpha * J_reac;
            J_R.makeCompressed();

            auto solve_linear_system = []( auto const& A, auto const& b )
            {
                using matrix_t = std::decay_t<decltype( A )>;
                using vector_t = std::decay_t<decltype( b )>;

                if constexpr ( std::same_as<matrix_t, Eigen::SparseMatrix<value_t>> )
                {
                    Eigen::SparseLU<matrix_t> solver;

                    solver.analyzePattern( A );
                    solver.factorize( A );

                    if ( solver.info() != Eigen::Success )
                    {
                        throw std::runtime_error( "SparseLU factorization failed in Eigen Shampine trick." );
                    }

                    vector_t x = solver.solve( b );

                    if ( solver.info() != Eigen::Success )
                    {
                        throw std::runtime_error( "SparseLU solve failed in Eigen Shampine trick." );
                    }

                    return x;
                }
                else
                {
                    vector_t x = A.colPivHouseholderQr().solve( b );
                    return x;
                }
            };

            if constexpr ( ell == 1 )
            {
                shampine_result = solve_linear_system( J_R, rhs );
            }
            else if constexpr ( ell == 2 )
            {
                u_tmp = solve_linear_system( J_R, rhs );

                shampine_result = solve_linear_system( J_R, u_tmp );
            }
        }
    };
} // namespace ponio::shampine_trick
