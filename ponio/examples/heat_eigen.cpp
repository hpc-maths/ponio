// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// NOLINTBEGIN(misc-include-cleaner)

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numbers>
#include <numeric>
#include <sstream>
#include <string>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <ponio/eigen_linear_algebra.hpp>
#include <ponio/observer.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

// NOLINTEND(misc-include-cleaner)

// Heat
struct heat_model
{
    using vector_type = Eigen::VectorX<double>;
    using matrix_type = Eigen::SparseMatrix<double>;

    int m_nx;
    double m_xmin;
    double m_xmax;
    double m_dx;

    vector_type x;
    matrix_type laplacian;

    heat_model( int nx, double xmin, double xmax )
        : m_nx( nx )
        , m_xmin( xmin )
        , m_xmax( xmax )
        , m_dx( ( xmax - xmin ) / static_cast<double>( nx ) )
        , x( nx )
        , laplacian( nx, nx )
    {
        std::vector<Eigen::Triplet<double>> tripletList( static_cast<std::size_t>( m_nx ) );

        tripletList.push_back( Eigen::Triplet<double>( 0, 0, -2 ) );
        tripletList.push_back( Eigen::Triplet<double>( 0, 1, 1 ) );

        x[0] = m_xmin;

        for ( int i = 1; i < m_nx - 1; ++i )
        {
            tripletList.push_back( Eigen::Triplet<double>( i, i - 1, 1 ) );
            tripletList.push_back( Eigen::Triplet<double>( i, i, -2 ) );
            tripletList.push_back( Eigen::Triplet<double>( i, i + 1, 1 ) );

            x[i] = m_xmin + i * m_dx;
        }
        tripletList.push_back( Eigen::Triplet<double>( m_nx - 1, m_nx - 2, 1 ) );
        tripletList.push_back( Eigen::Triplet<double>( m_nx - 1, m_nx - 1, -2 ) );

        x[m_nx - 1] = m_xmax;

        laplacian.setFromTriplets( tripletList.begin(), tripletList.end() );
    }

    void
    operator()( double t, vector_type const& y, vector_type& ydot ) const
    {
        f( t, y, ydot );
    }

    void
    f( double, vector_type const& y, vector_type& ydot ) const
    {
        double const oneoverdxdx = 1. / ( m_dx * m_dx );

        ydot[0] = oneoverdxdx * ( -2. * y[0] + 1. * y[1] );
        for ( int i = 1; i < m_nx - 1; ++i )
        {
            ydot[i] = oneoverdxdx * ( y[i - 1] - 2. * y[i] + y[i + 1] );
        }
        ydot[m_nx - 1] = oneoverdxdx * ( 1. * y[m_nx - 2] - 2. * y[m_nx - 1] );
    }

    matrix_type
    df( double, vector_type const& ) const
    {
        return laplacian;
    }

    vector_type
    fundamental_sol( double t ) const
    {
        double const xmid = 0.5 * ( m_xmax + m_xmin );
        double const pi   = std::numbers::pi;

        vector_type y( m_nx );
        for ( int i = 0; i < y.size(); ++i )
        {
            y[i] = ( 1. / ( 2. * std::sqrt( pi * t ) ) ) * ( std::exp( -( ( x[i] - xmid ) * ( x[i] - xmid ) ) / ( 4. * t ) ) );
        }

        return y;
    }

    double
    dx() const
    {
        return m_dx;
    }
};

template <typename vector_type>
void
save( vector_type const& x, vector_type const& y, std::filesystem::path const& filename )
{
    static auto printer_xy = []( double a, double b )
    {
        std::stringstream ss;
        ss << a << " " << b;
        return ss.str();
    };

    std::ofstream of( filename );
    std::transform( std::begin( x ), std::end( x ), std::begin( y ), std::ostream_iterator<std::string>( of, "\n" ), printer_xy );
    of.close();
}

template <typename vector_type>
double
error_l2( vector_type const& a, vector_type const& b, double dx )
{
    auto it_a = std::begin( a );

    double partial_sum = 0.;

    for ( auto it_b = std::begin( b ); it_b < std::end( b ); ++it_a, ++it_b )
    {
        partial_sum += ( ( *it_a ) - ( *it_b ) ) * ( ( *it_a ) - ( *it_b ) ) * dx;
    }

    return std::sqrt( partial_sum );
}

int
main()
{
    std::string const dirname = "heat_eigen_data";
    std::filesystem::create_directories( dirname );

    std::size_t const nx = 101;

    double const xmin = -5.0;
    double const xmax = 5.0;

    auto pb_heat         = heat_model( nx, xmin, xmax );
    auto eigmax_computer = [=]( auto&, double, heat_model::vector_type&, double, auto& )
    {
        return 4. / ( pb_heat.dx() * pb_heat.dx() );
    };

    double const t_ini = 0.1;
    double const t_end = 0.2;

    heat_model::vector_type const y_ini = pb_heat.fundamental_sol( t_ini );
    heat_model::vector_type y2_end( nx );
    heat_model::vector_type y4_end( nx );
    ponio::time_span<double> const tspan = { t_ini, t_end };

    save( pb_heat.x, y_ini, std::filesystem::path( dirname ) / "heat_ini.dat" );

    // make quasi-exact solution from SDIRK(3,4) with a small time step
    heat_model::vector_type const y_qexa = ponio::solve( pb_heat,
        ponio::runge_kutta::sdirk_34(),
        y_ini,
        tspan,
        1e-6,
        ponio::observer::null_observer() );
    ;
    save( pb_heat.x, y_qexa, std::filesystem::path( dirname ) / "heat_qexa.dat" );

    std::ofstream errors_file( std::filesystem::path( dirname ) / "errors.dat" );

    for ( std::size_t N = 1; N < 513; N *= 2 )
    {
        double const dt = ( t_end - t_ini ) / static_cast<double>( N );

        y2_end = ponio::solve( pb_heat, ponio::runge_kutta::rock::rock2( eigmax_computer ), y_ini, tspan, dt, ponio::observer::null_observer() );
        y4_end = ponio::solve( pb_heat, ponio::runge_kutta::rock::rock4( eigmax_computer ), y_ini, tspan, dt, ponio::observer::null_observer() );

        errors_file << dt << " " << std::setprecision( 20 ) << error_l2( y_qexa, y2_end, pb_heat.dx() ) << " "
                    << error_l2( y_qexa, y4_end, pb_heat.dx() ) << "\n";
    }

    errors_file.close();

    save( pb_heat.x, y2_end, std::filesystem::path( dirname ) / "heat_sol_rock2.dat" );
    save( pb_heat.x, y4_end, std::filesystem::path( dirname ) / "heat_sol_rock4.dat" );

    return 0;
}
