// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <numbers>
#include <set>
#include <stdexcept>
#include <string_view>
#include <type_traits>

#include <ponio/problem.hpp>
#include <ponio/runge_kutta/mrkc.hpp>
#include <ponio/runge_kutta/mrock.hpp>
#include <ponio/samurai_linear_algebra.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

#include <samurai/field.hpp>
#include <samurai/interface.hpp>
#include <samurai/mesh_config.hpp>
#include <samurai/mr/adapt.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/samurai.hpp>

struct partition_t
{
    std::set<std::size_t> active_levels;
    std::size_t min_level   = 0;
    std::size_t max_level   = 0;
    std::size_t split_level = 0;
    double rho_fast         = 0.;
    double rho_slow         = 0.;
};

enum class part_t
{
    fast,
    slow
};

template <std::size_t dim>
double
exact_solution( xt::xtensor_fixed<double, xt::xshape<dim>> const& x, double t, double nu )
{
    double result   = 1.;
    double const pi = std::numbers::pi;

    for ( std::size_t d = 0; d < dim; ++d )
    {
        result *= 1. / ( 2. * std::sqrt( pi * nu * t ) ) * std::exp( -x( d ) * x( d ) / ( 4. * nu * t ) );
    }

    return result;
}

template <class Field>
void
apply_bc( Field& u )
{
    u.get_bc().clear();
    samurai::make_bc<samurai::Neumann<1>>( u, 0. );
}

template <class Mesh>
partition_t
make_partition( Mesh& mesh, double nu )
{
    partition_t partition;

    samurai::for_each_cell( mesh,
        [&]( auto& cell )
        {
            partition.active_levels.insert( cell.level );
        } );

    if ( partition.active_levels.empty() )
    {
        throw std::runtime_error( "Cannot build a multirate partition on an empty mesh." );
    }

    partition.min_level = *partition.active_levels.begin();
    partition.max_level = *partition.active_levels.rbegin();

    if ( partition.min_level == partition.max_level )
    {
        partition.split_level = partition.max_level + 1;
        double const h        = mesh.cell_length( partition.max_level );
        partition.rho_slow    = 4. * nu / ( h * h );
        return partition;
    }

    std::size_t const delta = partition.max_level - partition.min_level;
    partition.split_level   = partition.min_level + ( 2 * delta + 2 ) / 3;

    std::size_t slow_level = partition.min_level;
    std::size_t fast_level = partition.max_level;

    for ( auto level : partition.active_levels )
    {
        if ( level < partition.split_level )
        {
            slow_level = level;
        }
        else
        {
            fast_level = level;
        }
    }

    double const h_slow = mesh.cell_length( slow_level );
    double const h_fast = mesh.cell_length( fast_level );

    partition.rho_slow = 4. * nu / ( h_slow * h_slow );
    partition.rho_fast = 4. * nu / ( h_fast * h_fast );

    return partition;
}

template <class Field>
class restricted_diffusion
{
    using field_t                    = Field;
    static constexpr std::size_t dim = field_t::dim;

    double m_nu;
    partition_t const* m_partition;

    bool
    keep_level( std::size_t level, part_t part ) const
    {
        bool const fast = level >= m_partition->split_level;
        return part == part_t::fast ? fast : !fast;
    }

    void
    add_face( field_t& u, field_t& du, auto const& interface_cells, auto const& comput_cells ) const
    {
        auto& mesh                   = u.mesh();
        std::size_t const face_level = std::max( interface_cells[0].level, interface_cells[1].level );
        double const h_face          = mesh.cell_length( face_level );
        double const h_left          = mesh.cell_length( interface_cells[0].level );
        double const h_right         = mesh.cell_length( interface_cells[1].level );

        double const flux = m_nu / h_face * ( static_cast<double>( u[comput_cells[0]] ) - static_cast<double>( u[comput_cells[1]] ) );

        du[interface_cells[0]] -= flux / h_left;
        du[interface_cells[1]] += flux / h_right;
    }

  public:

    restricted_diffusion( double nu, partition_t const& partition )
        : m_nu( nu )
        , m_partition( &partition )
    {
        static_assert( field_t::is_scalar );
    }

    void
    operator()( field_t& u, field_t& du, part_t part ) const
    {
        du.fill( 0. );
        auto& mesh = u.mesh();

        for ( std::size_t d = 0; d < dim; ++d )
        {
            samurai::DirectionVector<dim> direction;
            direction.fill( 0 );
            direction[d] = 1;

            samurai::Stencil<2, dim> stencil = samurai::in_out_stencil<dim>( direction );
            auto analyzer                    = samurai::make_stencil_analyzer( stencil );

            for ( auto level : m_partition->active_levels )
            {
                if ( !keep_level( level, part ) )
                {
                    continue;
                }

                samurai::for_each_interior_interface__same_level( mesh,
                    level,
                    direction,
                    analyzer,
                    [&]( auto& interface_cells, auto& comput_cells )
                    {
                        add_face( u, du, interface_cells, comput_cells );
                    } );
            }

            for ( std::size_t level = m_partition->min_level; level < m_partition->max_level; ++level )
            {
                std::size_t const fine_level = level + 1;
                if ( !keep_level( fine_level, part ) )
                {
                    continue;
                }

                samurai::for_each_interior_interface__level_jump_direction( mesh,
                    level,
                    direction,
                    analyzer,
                    [&]( auto& interface_cells, auto& comput_cells )
                    {
                        add_face( u, du, interface_cells, comput_cells );
                    } );

                samurai::for_each_interior_interface__level_jump_opposite_direction( mesh,
                    level,
                    direction,
                    analyzer,
                    [&]( auto& interface_cells, auto& comput_cells )
                    {
                        add_face( u, du, interface_cells, comput_cells );
                    } );
            }
        }
    }
};

template <class Field>
double
l1_error( Field const& u, double t, double nu )
{
    double error  = 0.;
    double volume = 0.;

    samurai::for_each_cell( u.mesh(),
        [&]( auto& cell )
        {
            double const dx    = cell.length;
            double const exact = exact_solution<Field::dim>( cell.center(), t, nu );
            error += std::abs( static_cast<double>( u[cell] ) - exact ) * dx;
            volume += dx;
        } );

    return error / volume;
}

template <class MethodFactory>
void
run_multirate( std::string_view name, MethodFactory&& make_method )
{
    constexpr std::size_t dim = 1;
    double const nu           = 1.;
    double const t_ini        = 1.e-2;
    double const t_end        = 2.e-1;

    samurai::Box<double, dim> const box( { -20. }, { 20. } );
    auto config = samurai::mesh_config<dim>().min_level( 8 ).max_level( 9 ).max_stencil_size( 2 ).disable_minimal_ghost_width();

    auto mesh  = samurai::mra::make_mesh( box, config );
    auto u_ini = samurai::make_scalar_field<double>( "u", mesh );

    samurai::for_each_cell( mesh,
        [&]( auto& cell )
        {
            u_ini[cell] = exact_solution<dim>( cell.center(), t_ini, nu );
        } );

    apply_bc( u_ini );

    auto mra_config         = samurai::mra_config().epsilon( 1.e-5 ).regularity( 1. ).relative_detail( true );
    auto initial_adaptation = samurai::make_MRAdapt( u_ini );
    initial_adaptation( mra_config );
    apply_bc( u_ini );
    samurai::update_ghost_mr( u_ini );

    partition_t partition = make_partition( mesh, nu );
    restricted_diffusion<decltype( u_ini )> diffusion( nu, partition );

    auto fast = [&]( double, auto const& u_const, auto& du )
    {
        using field_t = std::remove_const_t<std::remove_reference_t<decltype( u_const )>>;
        auto& u       = const_cast<field_t&>( u_const );
        apply_bc( u );
        samurai::update_ghost_mr( u );
        diffusion( u, du, part_t::fast );
    };

    auto slow = [&]( double, auto const& u_const, auto& du )
    {
        using field_t = std::remove_const_t<std::remove_reference_t<decltype( u_const )>>;
        auto& u       = const_cast<field_t&>( u_const );
        apply_bc( u );
        samurai::update_ghost_mr( u );
        diffusion( u, du, part_t::slow );
    };

    auto problem = ponio::make_problem( fast, slow );
    auto method  = make_method( partition.rho_fast, partition.rho_slow );

    double const h_slow                   = std::sqrt( 4. * nu / partition.rho_slow );
    double const dt                       = 0.5 * h_slow * h_slow / nu;
    ponio::time_span<double> const t_span = { t_ini, t_end };

    auto sol_range     = ponio::make_solver_range( problem, method, u_ini, t_span, dt );
    auto it_sol        = sol_range.begin();
    auto mr_adaptation = samurai::make_MRAdapt( it_sol->state );

    while ( it_sol->time < t_end )
    {
        it_sol.callback_on_transient_stages(
            []( auto& stage )
            {
                stage.resize();
                stage.fill( 0. );
            } );

        ++it_sol;

        mr_adaptation( mra_config );
        apply_bc( it_sol->state );
        samurai::update_ghost_mr( it_sol->state );

        partition = make_partition( it_sol->state.mesh(), nu );
        it_sol.meth.alg.set_spectral_radii( partition.rho_fast, partition.rho_slow );
    }

    std::cout << name << ": L1 norm of error = " << l1_error( it_sol->state, it_sol->time, nu ) << '\n';
}

int
main( int argc, char** argv )
{
    samurai::initialize( "Heat equation with multirate methods", argc, argv );
    SAMURAI_PARSE( argc, argv );

    run_multirate( "mRKC",
        []( double rho_fast, double rho_slow )
        {
            return ponio::runge_kutta::chebyshev::mrkc( rho_fast, rho_slow );
        } );

    run_multirate( "mROCK2",
        []( double rho_fast, double rho_slow )
        {
            return ponio::runge_kutta::mrock::mrock2( rho_fast, rho_slow );
        } );

    samurai::finalize();
    return 0;
}
