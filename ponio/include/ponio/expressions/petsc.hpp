// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include "state.hpp"
#include <petsc.h>

namespace ponio::expression
{

    template <>
    struct state<Vec>
    {
        static constexpr bool is_ponio_expression = true;

        using container_type = Vec;
        container_type& _data; // store only a reference on data
        double const* _reader_data = nullptr;

        state( container_type& dat )
            : _data( dat )
        {
        }

        auto
        operator[]( std::size_t i ) const
        {
            if ( _reader_data == nullptr )
            {
                int ni = 1;
                int ix = static_cast<int>( i );

                double data_i;

                VecGetValues( _data, ni, &ix, &data_i );

                return data_i;
            }

            return _reader_data[i];
        }

        template <typename expression_t>
            requires detail::is_ponio_expression<expression_t>
        state&
        operator=( expression_t& expr )
        {
            double* _ptr_data;
            VecGetArray( _data, &_ptr_data );

            expr.raw_data_open();

            for ( std::size_t i = 0ul; i < expr.size(); ++i )
            {
                _ptr_data[i] = expr[i];
            }

            expr.raw_data_close();

            VecRestoreArray( _data, &_ptr_data );

            return *this;
        }

        template <typename expression_t>
            requires detail::is_ponio_expression<expression_t>
        state&
        operator=( expression_t&& expr )
        {
            double* _ptr_data;
            VecGetArray( _data, &_ptr_data );

            expr.raw_data_open();

            for ( std::size_t i = 0ul; i < expr.size(); ++i )
            {
                _ptr_data[i] = expr[i];
            }

            expr.raw_data_close();

            VecRestoreArray( _data, &_ptr_data );

            return *this;
        }

        std::size_t
        size() const
        {
            int ssize;
            VecGetSize( _data, &ssize );
            return static_cast<std::size_t>( ssize );
        }

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
        raw_data_open()
        {
            VecGetArrayRead( _data, &_reader_data );
        }

        void
        raw_data_close()
        {
            VecRestoreArrayRead( _data, &_reader_data );
        }
    };

    template <>
    struct state<Mat>
    {
        static constexpr bool is_ponio_expression = true;

        using container_type = Mat;
        container_type& _data; // store only a reference on data
        double const* _reader_data = nullptr;

        state( container_type& dat )
            : _data( dat )
        {
        }

        auto
        operator[]( std::size_t idx ) const
        {
            if ( _reader_data == nullptr )
            {
                int nb_elm = 1;

                int m, n;
                MatGetSize( _data, &m, &n );

                // only for row-oriented dense matrix
                int idxm = static_cast<int>( idx / m );
                int idxn = static_cast<int>( idx % n );

                double data_i;

                MatGetValues( _data, nb_elm, &idxm, nb_elm, &idxn, &data_i );

                return data_i;
            }

            return _reader_data[idx];
        }

        template <typename expression_t>
            requires detail::is_ponio_expression<expression_t>
        state&
        operator=( expression_t& expr )
        {
            double* _ptr_data;
            MatDenseGetArray( _data, &_ptr_data );

            expr.raw_data_open();

            for ( std::size_t i = 0ul; i < expr.size(); ++i )
            {
                _ptr_data[i] = expr[i];
            }

            expr.raw_data_close();

            MatDenseRestoreArray( _data, &_ptr_data );

            return *this;
        }

        template <typename expression_t>
            requires detail::is_ponio_expression<expression_t>
        state&
        operator=( expression_t&& expr )
        {
            double* _ptr_data;
            MatDenseGetArray( _data, &_ptr_data );

            expr.raw_data_open();

            for ( std::size_t i = 0ul; i < expr.size(); ++i )
            {
                _ptr_data[i] = expr[i];
            }

            expr.raw_data_close();

            MatDenseRestoreArray( _data, &_ptr_data );

            return *this;
        }

        std::size_t
        size() const
        {
            int m, n;
            MatGetSize( _data, &m, &n );
            return static_cast<std::size_t>( m * n );
        }

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
        raw_data_open()
        {
            MatDenseGetArrayRead( _data, &_reader_data );
        }

        void
        raw_data_close()
        {
            MatDenseRestoreArrayRead( _data, &_reader_data );
        }
    };
} // namespace ponio::expression
