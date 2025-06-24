// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <samurai/field.hpp>
#include <samurai/io/hdf5.hpp>
#include <samurai/mr/adapt.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/schemes/fv.hpp>

#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/samurai_linear_algebra.hpp>
#include <ponio/solver.hpp>

#include <filesystem>

template <class Field>
void
save( std::filesystem::path const& path, std::string const& filename, Field& u, std::string const& suffix = "" )
{
    auto mesh   = u.mesh();
    auto level_ = samurai::make_scalar_field<std::size_t>( "level", mesh );
    u.name()    = "u";

    if ( !std::filesystem::exists( path ) )
    {
        std::filesystem::create_directory( path );
    }

    samurai::for_each_cell( mesh,
        [&]( auto& cell )
        {
            level_[cell] = cell.level;
        } );

    samurai::save( path, fmt::format( "{}{}", filename, suffix ), mesh, u, level_ );
}

template <class Mesh>
auto
init( Mesh& mesh )
{
    auto u = samurai::make_scalar_field<double>( "u", mesh );
    u.fill( 0. );

    samurai::for_each_cell( mesh,
        [&]( auto& cell )
        {
            auto center         = cell.center();
            double const radius = .2;

            double const x_center = 0;
            if ( std::abs( center[0] - x_center ) <= radius )
            {
                u[cell] = 1;
            }
        } );

    return u;
}

int
main( int argc, char** argv )
{
    PetscInitialize( &argc, &argv, nullptr, nullptr );

    constexpr std::size_t dim = 1; // cppcheck-suppress unreadVariable
    using Config              = samurai::MRConfig<dim>;

    // Simulation parameters
    double const left_box  = -5;
    double const right_box = 5;
    double const Tf        = 0.5;
    double const cfl       = 0.5;

    // Multiresolution parameters
    std::size_t const min_level = 2;
    std::size_t const max_level = 5;
    double const mr_epsilon     = 2.e-4; // Threshold used by multiresolution
    double const mr_regularity  = 1.;    // Regularity guess for multiresolution

    // Output parameters
    std::string const dirname  = "heat_samurai_data";
    std::filesystem::path path = std::filesystem::path( dirname );
    std::string const filename = "sol_1d";

    // Define mesh
    samurai::Box<double, dim> const box( { left_box }, { right_box } );
    samurai::MRMesh<Config> mesh( box, min_level, max_level );

    // Initial condition
    auto un_ini = init( mesh );

    // Define problem (diffusion), solve : d_t u = - d_xx u
    samurai::DiffCoeff<dim> diff_coeff;
    diff_coeff.fill( 1.0 );
    auto diff = samurai::make_diffusion_order2<decltype( un_ini )>( diff_coeff );

    auto f_t = [&]( [[maybe_unused]] double t )
    {
        return -diff;
    };
    auto f = [&]( [[maybe_unused]] double t, auto&& u, auto& du )
    {
        samurai::make_bc<samurai::Neumann<1>>( u, 0. );
        samurai::update_ghost_mr( u );
        du = f_t( t )( u );
    };
    auto pb = ponio::make_implicit_operator_problem( f, f_t );

    // Time step
    double const dx_min = mesh.cell_length( mesh.max_level() );
    double const dt     = cfl * dx_min * dx_min;

    auto eigmax_computer = [=]( auto&, double, auto&, double, auto& )
    {
        return 4. / ( dx_min * dx_min );
    };

    // Range to iterate over solution
    // auto sol_range = ponio::make_solver_range( pb, ponio::runge_kutta::rkc_202(), un_ini, { 0., Tf }, dt );
    // auto sol_range = ponio::make_solver_range( pb, ponio::runge_kutta::dirk23(), un_ini, { 0., Tf }, dt );
    auto sol_range = ponio::make_solver_range( pb, ponio::runge_kutta::rock::rock4<true>( eigmax_computer ), un_ini, { 0., Tf }, dt );

    auto it_sol = sol_range.begin();
    samurai::make_bc<samurai::Neumann<1>>( it_sol->state, 0. );

    // Prepare MR for solution on iterator
    auto MRadaptation = samurai::make_MRAdapt( it_sol->state );
    MRadaptation( mr_epsilon, mr_regularity );
    samurai::update_ghost_mr( it_sol->state );

    std::size_t n_save = 0;
    save( path, filename, it_sol->state, fmt::format( "_ite_{}", n_save++ ) );

    std::cerr << "> time loop" << std::endl;
    while ( it_sol->time < Tf )
    {
        // samurai::make_bc<samurai::Neumann<1>>( it_sol->state, 0. );

        for ( auto& ki : it_sol.stages() )
        {
            ki.resize();
            ki.fill( 0. );
        }

        ++it_sol;

        std::cout << "tⁿ: " << std::setw( 8 ) << it_sol->time << " (Δt: " << it_sol->time_step << ")  "
                  << "N stages:" << it_sol.info().number_of_stages << "   \r";

        MRadaptation( mr_epsilon, mr_regularity );
        samurai::update_ghost_mr( it_sol->state );

        save( path, filename, it_sol->state, fmt::format( "_ite_{}", n_save++ ) );
    }
    std::cout << std::endl;

    PetscFinalize();

    return 0;
}
