// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <xtensor/xfixed.hpp>

#include <samurai/algorithm.hpp>
#include <samurai/bc.hpp>
#include <samurai/field.hpp>
#include <samurai/hdf5.hpp>
#include <samurai/mr/adapt.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/petsc.hpp>
#include <samurai/stencil_field.hpp>
#include <samurai/subset/subset_op.hpp>

#include <solver/observer.hpp>
#include <solver/problem.hpp>
#include <solver/runge_kutta.hpp>
#include <solver/solver.hpp>

#include <filesystem>
namespace fs = std::filesystem;

template <class Field>
void
save( fs::path const& path, std::string const& filename, Field const& u, std::string const& suffix = "" )
{
    auto mesh   = u.mesh();
    auto level_ = samurai::make_field<std::size_t, 1>( "level", mesh );

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
    double left_box  = -2;
    double right_box = 2;
    bool is_periodic = false;
    double a         = 1.;
    double Tf        = 1.;
    double cfl       = 0.95;

    // Multiresolution parameters
    std::size_t min_level = 4;
    std::size_t max_level = 10;
    double mr_epsilon     = 2.e-4; // Threshold used by multiresolution
    double mr_regularity  = 1.;    // Regularity guess for multiresolution
    bool correction       = false;

    // Output parameters
    fs::path path        = fs::current_path();
    std::string filename = "FV_advection_1d";
    std::size_t nfiles   = 1;

    samurai::Box<double, dim> const box( { left_box }, { right_box } );
    samurai::MRMesh<Config> mesh( box, min_level, max_level, { is_periodic } );

    double dt            = cfl * samurai::cell_length( max_level );
    double const dt_save = Tf / static_cast<double>( nfiles );
    double t             = 0.;

    auto un = init( mesh );

    auto MRadaptation = samurai::make_MRAdapt( un );
    samurai::make_bc<samurai::Neumann>( un, 0. );
    MRadaptation( mr_epsilon, mr_regularity );
    samurai::update_ghost_mr( un );

    std::size_t n_save = 0, nt = 0;
    save( path, filename, un, fmt::format( "_ite_{}", n_save++ ) );

    samurai::DiffCoeff<dim> K;
    K.fill( 1.0 );
    auto diff = samurai::make_diffusion<decltype( un )>( K );

    auto pb = [&]( double t, auto ui )
    {
        return diff( ui );
    };

    auto sol_range = ponio::make_solver_range( pb, ponio::runge_kutta::rkc_202(), un, { 0., Tf }, dt );

    for ( auto ui : sol_range )
    {
        MRadaptation( mr_epsilon, mr_regularity );
        samurai::update_ghost_mr( un );
    }

    save( path, filename, un, fmt::format( "_ite_{}", n_save++ ) );

    return 0;
}
