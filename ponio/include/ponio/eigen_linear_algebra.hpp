// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

// NOLINTBEGIN(misc-include-cleaner)

#if __has_include( <eigen3/Eigen/Dense> )
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#elif __has_include( <Eigen/Dense> )
#include <Eigen/Dense>
#include <Eigen/Sparse>
#else
#error "Eigen should be included"
#endif

// NOLINTEND(misc-include-cleaner)

#include <stdexcept>
#include <type_traits>
#include <utility>

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

        static matrix_type
        identity_minus( matrix_type const& matrix, scalar_t alpha )
        {
            return matrix_type::Identity() - alpha * matrix;
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

        static matrix_type
        identity_minus( matrix_type const& matrix, scalar_t alpha )
        {
            return matrix_type::Identity( matrix.rows(), matrix.cols() ) - alpha * matrix;
        }

        static vector_type
        solver( matrix_type const& dfx, vector_type const& fx )
        {
            return dfx.colPivHouseholderQr().solve( fx );
        }
    };

    /**
     * @brief Linear algebra specialization for Eigen sparse matrices.
     *
     * SparseLU is used instead of a Cholesky factorization because the linear
     * systems arising from implicit methods are not assumed to be symmetric or
     * positive definite. This makes the sparse Eigen backend applicable to a
     * broader class of problems.
     *
     * The symbolic analysis and the numerical factorization are deliberately
     * separated. For a fixed sparse pattern, the symbolic analysis is performed
     * only once, while the numerical factorization may be updated whenever the
     * matrix coefficients change.
     *
     * The factorization is stored in a static solver so that several solves can
     * reuse it. In particular, this is required by simplified Newton iterations and
     * by Shampine's trick, where the same factorized matrix may be applied more
     * than once.
     */
    template <typename scalar_t>
    struct linear_algebra<Eigen::SparseMatrix<scalar_t>> // NOLINT(misc-include-cleaner)
    {
        using matrix_type = Eigen::SparseMatrix<scalar_t>;           // NOLINT(misc-include-cleaner)
        using vector_type = Eigen::Vector<scalar_t, Eigen::Dynamic>; // NOLINT(misc-include-cleaner)
        using solver_type = Eigen::SparseLU<matrix_type>;            // NOLINT(misc-include-cleaner)

      private:

        // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
        inline static matrix_type I;

        /**
         * The sparse solver is shared by all calls using this specialization.
         * Its numerical factorization therefore remains available after
         * factorize() returns and can be reused by solve_factorized().
         */
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
        inline static solver_type sparse_solver;

        /**
         * The current implementation assumes a fixed sparse pattern. The symbolic
         * analysis can therefore be reused after the first successful factorization.
         */
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
        inline static bool pattern_analyzed = false;

      public:

        static matrix_type const&
        identity( vector_type const& u )
        {
            I.resize( u.size(), u.size() );
            I.setIdentity();

            return I;
        }

        /**
         * @brief Build the shifted matrix I - alpha A.
         *
         * Sparse storage details are kept inside the Eigen backend so that
         * algorithms using linear_algebra do not depend on Eigen-specific
         * operations such as coeffRef() or makeCompressed().
         */
        static matrix_type
        identity_minus( matrix_type const& matrix, scalar_t alpha )
        {
            matrix_type shifted_matrix = matrix;
            shifted_matrix *= -alpha;

            for ( Eigen::Index i = 0; i < shifted_matrix.rows(); ++i )
            {
                shifted_matrix.coeffRef( i, i ) += scalar_t( 1 );
            }

            shifted_matrix.makeCompressed();

            return shifted_matrix;
        }

        /**
         * @brief Analyze and factorize a sparse matrix.
         *
         * The symbolic analysis is performed only on the first call. Subsequent
         * calls reuse the same symbolic structure and update only the numerical
         * factorization.
         *
         * @param matrix Matrix to factorize.
         */
        static void
        factorize( matrix_type const& matrix )
        {
            if ( !pattern_analyzed )
            {
                sparse_solver.analyzePattern( matrix );

                // SparseLU::info() must not be queried here because the numerical factorization has not been initialized yet.
            }

            sparse_solver.factorize( matrix );

            if ( sparse_solver.info() != Eigen::Success )
            {
                throw std::runtime_error( "Eigen SparseLU factorization failed in "
                                          "ponio::linear_algebra." );
            }

            // The symbolic pattern is considered valid only after a successful numerical factorization.
            pattern_analyzed = true;
        }

        /**
         * @brief Solve a system using the factorization stored by factorize().
         *
         * This function intentionally does not perform a new factorization. It
         * is used when several right-hand sides, Newton corrections, or
         * Shampine inverse applications share the same linear system.
         *
         * @param rhs Right-hand side of the linear system.
         *
         * @return Solution of the factorized system.
         */
        static vector_type
        solve_factorized( vector_type const& rhs )
        {
            vector_type solution;

            solution = sparse_solver.solve( rhs );

            if ( sparse_solver.info() != Eigen::Success )
            {
                throw std::runtime_error( "Eigen SparseLU solve failed in "
                                          "ponio::linear_algebra." );
            }

            return solution;
        }

        /**
         * @brief Factorize a matrix and solve the associated linear system.
         *
         * This function preserves the original linear_algebra interface. Code
         * that needs to reuse the same factorization should call factorize()
         * once and solve_factorized() for each right-hand side.
         */
        static vector_type
        solver( matrix_type const& matrix, vector_type const& rhs )
        {
            factorize( matrix );

            return solve_factorized( rhs );
        }
    };

} // namespace ponio::linear_algebra

namespace ponio::shampine_trick
{

    template <typename state_t>
    struct shampine_trick;

    /**
     * @brief Sparse Eigen implementation of Shampine's trick.
     *
     * Shampine's trick applies the inverse of
     *
     *     I - alpha J
     *
     * one or two times, depending on the value of ell. The initialization,
     * application and finalization steps are kept explicit so that PIROCK only
     * depends on a small Shampine interface and not on Eigen-specific solver
     * details.
     *
     * The object can either prepare its own SparseLU factorization or reuse the
     * factorization already computed by the sparse linear algebra backend. In
     * both cases, every call to apply() uses the same prepared factorization.
     */
    template <typename scalar_t, int size, int options, int maxrows>
    struct shampine_trick<Eigen::Matrix<scalar_t, size, 1, options, maxrows, 1>> // NOLINT(misc-include-cleaner)
    {
        using value_t         = scalar_t;
        using sparse_matrix_t = Eigen::SparseMatrix<value_t>;     // NOLINT(misc-include-cleaner)
        using sparse_solver_t = Eigen::SparseLU<sparse_matrix_t>; // NOLINT(misc-include-cleaner)

      private:

        sparse_matrix_t J_reac_;
        sparse_matrix_t J_R_;
        sparse_solver_t sparse_solver_;

        /**
         * The sparse pattern is fixed because no mesh adaptation occurs in the
         * Eigen path. The symbolic analysis can therefore be reused throughout
         * the lifetime of the Shampine object.
         */
        bool sparse_pattern_analyzed_ = false;

        /**
         * When true, apply() uses the factorization previously prepared by
         * ponio::linear_algebra. Otherwise, the factorization owned by this
         * Shampine object is used.
         */
        bool use_external_factorization_ = false;

        void
        analyze_pattern_if_needed()
        {
            if ( !sparse_pattern_analyzed_ )
            {
                sparse_solver_.analyzePattern( J_R_ );

                // SparseLU::info() must not be queried here because the numerical factorization has not been initialized yet.
            }
        }

        /**
         * @brief Solve a system with the factorization selected by initialize().
         */
        template <typename state_t>
        state_t
        solve( state_t const& rhs )
        {
            if ( use_external_factorization_ )
            {
                return ::ponio::linear_algebra::linear_algebra<sparse_matrix_t>::solve_factorized( rhs );
            }

            state_t solution = sparse_solver_.solve( rhs );

            if ( sparse_solver_.info() != Eigen::Success )
            {
                throw std::runtime_error( "Eigen SparseLU solve failed in Shampine's trick." );
            }

            return solution;
        }

      public:

        /**
         * @brief Initialize Shampine's trick with a locally factorized system.
         */
        template <typename jacobian_reac_t, typename state_t>
        void
        initialize( value_t alpha, jacobian_reac_t&& jacobian_reac, state_t const& initial_guess )
        {
            using returned_matrix_t = std::decay_t<decltype( jacobian_reac( initial_guess ) )>;

            static_assert( std::is_same_v<returned_matrix_t, sparse_matrix_t>,
                "The Eigen implementation of Shampine's trick expects "
                "Eigen::SparseMatrix<value_t>." );

            J_reac_ = jacobian_reac( initial_guess );
            J_R_    = ::ponio::linear_algebra::linear_algebra<sparse_matrix_t>::identity_minus( J_reac_, alpha );

            analyze_pattern_if_needed();
            sparse_solver_.factorize( J_R_ );

            if ( sparse_solver_.info() != Eigen::Success )
            {
                throw std::runtime_error( "Eigen SparseLU factorization failed in "
                                          "Shampine's trick." );
            }

            // The symbolic pattern is considered reusable only after a successful numerical factorization.
            sparse_pattern_analyzed_    = true;
            use_external_factorization_ = false;
        }

        /**
         * @brief Initialize Shampine's trick from an existing factorization.
         *
         * PIROCK uses this overload after its implicit stages have already
         * factorized the same matrix I - gamma dt J_R. No symbolic analysis,
         * matrix assembly or numerical factorization is repeated here.
         */
        void
        initialize()
        {
            use_external_factorization_ = true;
        }

        /**
         * @brief Apply (I - alpha J)^(-ell) to a right-hand side.
         */
        template <std::size_t ell, typename state_t>
        void
        apply( state_t const& rhs, state_t& u_tmp, state_t& result )
        {
            if constexpr ( ell == 1 )
            {
                result = solve( rhs );
            }
            else if constexpr ( ell == 2 )
            {
                u_tmp  = solve( rhs );
                result = solve( u_tmp );
            }
            else
            {
                static_assert( ell == 1 || ell == 2, "Eigen Shampine's trick only supports ell = 1 or ell = 2." );
            }
        }

        /**
         * @brief Finalize the current Shampine application sequence.
         */
        void
        finalize()
        {
            use_external_factorization_ = false;
        }

        /**
         * @brief Assemble, factorize and apply Shampine's trick.
         *
         * This operator preserves the original direct-call interface. New code
         * can use initialize(), apply() and finalize() explicitly when several
         * inverse applications share the same factorization.
         */
        template <std::size_t ell, typename jacobian_reac_t, typename state_t>
        void
        operator()( value_t alpha, jacobian_reac_t&& jacobian_reac, state_t& initial_guess, state_t& rhs, state_t& u_tmp, state_t& shampine_result )
        {
            initialize( alpha, std::forward<jacobian_reac_t>( jacobian_reac ), initial_guess );
            apply<ell>( rhs, u_tmp, shampine_result );
            finalize();
        }
    };

} // namespace ponio::shampine_trick
