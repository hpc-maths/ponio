// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <samurai/field.hpp>
#include <samurai/hdf5.hpp>
#include <samurai/mr/adapt.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/schemes/fv.hpp>

#include <solver/observer.hpp>
#include <solver/problem.hpp>
#include <solver/runge_kutta.hpp>
#include <solver/samurai_linear_algebra.hpp>
#include <solver/solver.hpp>

#include <filesystem>
namespace fs = std::filesystem;

template <class Field>
void
save( fs::path const& path, std::string const& filename, Field& u, std::string const& suffix = "" )
{
    auto mesh   = u.mesh();
    auto level_ = samurai::make_field<std::size_t, 1>( "level", mesh );
    u.name()    = "u";

    if ( !fs::exists( path ) )
    {
        fs::create_directory( path );
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
    auto u = samurai::make_field<double, 1>( "u", mesh );
    u.fill( 0. );

    samurai::for_each_cell( mesh,
        [&]( auto& cell )
        {
            auto center         = cell.center();
            const double radius = .2;

            const double x_center = 0;
            if ( std::abs( center[0] - x_center ) <= radius )
            {
                u[cell] = 1;
            }
        } );

    return u;
}

int
main()
{
    constexpr std::size_t dim = 1; // cppcheck-suppress unreadVariable
    using Config              = samurai::MRConfig<dim>;

    // Simulation parameters
    double left_box  = -5;
    double right_box = 5;
    bool is_periodic = false;
    double Tf        = 0.5;
    double cfl       = 0.5;

    // Multiresolution parameters
    std::size_t min_level = 2;
    std::size_t max_level = 5;
    double mr_epsilon     = 2.e-4; // Threshold used by multiresolution
    double mr_regularity  = 1.;    // Regularity guess for multiresolution

    // Output parameters
    std::string const dirname = "heat_samurai_data";
    fs::path path             = std::filesystem::path( dirname );
    std::string filename      = "sol_1d";

    // Define mesh
    samurai::Box<double, dim> const box( { left_box }, { right_box } );
    samurai::MRMesh<Config> mesh( box, min_level, max_level, { is_periodic } );

    // Initial condition
    auto un_ini = init( mesh );

    // Define problem (diffusion), solve : d_t u = - d_xx u
    samurai::DiffCoeff<dim> diff_coeff;
    diff_coeff.fill( 1.0 );
    auto diff = samurai::make_diffusion<decltype( un_ini )>( diff_coeff );

    auto f_t = [&]( [[maybe_unused]] double t )
    {
        return -diff;
    };

    auto f = [&]( [[maybe_unused]] double t, auto&& u )
    {
        samurai::make_bc<samurai::Neumann>( u, 0. );
        samurai::update_ghost_mr( u );
        return f_t( t )( u );
    };
    auto pb = ponio::make_implicit_operator_problem( f, f_t );

    // Time step
    double dt = cfl * samurai::cell_length( max_level ) * samurai::cell_length( max_level );

    // Range to iterate over solution
    // auto sol_range = ponio::make_solver_range( pb, ponio::runge_kutta::rkc_202(), un_ini, { 0., Tf }, dt );
    auto sol_range = ponio::make_solver_range( pb, ponio::runge_kutta::dirk23(), un_ini, { 0., Tf }, dt );
    auto it_sol    = sol_range.begin();

    // Prepare MR for solution on iterator
    auto MRadaptation = samurai::make_MRAdapt( it_sol->state );
    samurai::make_bc<samurai::Neumann>( it_sol->state, 0. );
    MRadaptation( mr_epsilon, mr_regularity );
    samurai::update_ghost_mr( it_sol->state );

    std::size_t n_save = 0;
    save( path, filename, it_sol->state, fmt::format( "_ite_{}", n_save++ ) );

    while ( it_sol->time < Tf )
    {
        // TODO: add a callback function to make this before each iteration
        for ( auto& ki : it_sol.meth.kis )
        {
            ki.resize();
            ki.fill( 0. );
        }

        ++it_sol;

        std::cout << "tⁿ: " << std::setw( 8 ) << it_sol->time << " (Δt: " << it_sol->time_step << ")\r";

        MRadaptation( mr_epsilon, mr_regularity );
        samurai::update_ghost_mr( it_sol->state );

        save( path, filename, it_sol->state, fmt::format( "_ite_{}", n_save++ ) );
    }
    std::cout << std::endl;

    return 0;
}
