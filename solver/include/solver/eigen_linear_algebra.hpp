// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#if __has_include( <eigen3/Eigen/Dense>)
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#elif __has_include( <Eigen/Dense>)
#include <Eigen/Dense>
#include <Eigen/Sparse>
#else
#error "Eigen should be included"
#endif

#include "linear_algebra.hpp"

namespace ponio::linear_algebra
{

    template <typename scalar_t, int size, int options, int maxrows, int maxcols>
    struct linear_algebra<Eigen::Matrix<scalar_t, size, size, options, maxrows, maxcols>>
    {
        using matrix_type = Eigen::Matrix<scalar_t, size, size>;
        using vector_type = Eigen::Vector<scalar_t, size>;

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
    struct linear_algebra<Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>>
    {
        using matrix_type = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;
        using vector_type = Eigen::Vector<scalar_t, Eigen::Dynamic>;

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
    struct linear_algebra<Eigen::SparseMatrix<scalar_t>>
    {
        using matrix_type = Eigen::SparseMatrix<scalar_t>;
        using vector_type = Eigen::SparseVector<scalar_t>;
        using solver_type = Eigen::SimplicialCholesky<matrix_type>;

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
