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
            return dfx.colPivHouseholderQr().solve( fx );
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
        using matrix_type = Eigen::SparseMatrix<scalar_t>;          // NOLINT(misc-include-cleaner)
        using vector_type = Eigen::SparseVector<scalar_t>;          // NOLINT(misc-include-cleaner)
        using solver_type = Eigen::SimplicialCholesky<matrix_type>; // NOLINT(misc-include-cleaner)

      private:

        static matrix_type I; // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)

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
