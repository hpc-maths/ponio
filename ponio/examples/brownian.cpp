// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// NOLINTBEGIN(misc-include-cleaner)

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <valarray>

#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>

#include <CLI/CLI.hpp>

// NOLINTEND(misc-include-cleaner)

struct C_random_device
{
    using result_type = std::uint_fast32_t;

    static constexpr result_type
    min()
    {
        return static_cast<result_type>( 0 );
    }

    static constexpr result_type
    max()
    {
        return static_cast<result_type>( RAND_MAX );
    }

    static double
    entropy() noexcept // cppcheck-suppress unusedFunction
    {
        return 42.0;
    }

    result_type
    operator()()
    {
        return static_cast<result_type>( std::rand() ); // NOLINT: GitHub Action doesn't allow `std::random_device`
    }
};

int
main( int argc, char* argv[] )
{
    CLI::App app{ "Launch N brownian motion process solved by a RK(3,3) method" }; // NOLINT(misc-include-cleaner)

    std::string const dirname = "brownian_data";

    std::size_t n = 10;
    app.add_option( "N", n, "Number of brownian motion" );
    CLI11_PARSE( app, argc, argv ); // NOLINT(misc-include-cleaner)

    using state_t = std::valarray<double>;

    // std::random_device rd; // can't be use in GitHub Action...
    C_random_device rd;
    std::mt19937 gen( rd() );

    std::normal_distribution<> d{ 0., 2 };

    double const dt = 1e-3;

    auto brownian_pb = ponio::make_simple_problem(
        [&]( double, state_t const& ) -> state_t
        {
            return { d( gen ), d( gen ) };
        } );

    state_t const yini = { 0., 0. };

    for ( unsigned int i = 0; i < n; ++i )
    {
        std::stringstream ssfilename;
        ssfilename << "brownian_" << i << ".dat";
        auto filename = std::filesystem::path( dirname ) / ssfilename.str();
        ponio::observer::file_observer fobs( filename );
        ponio::solve( brownian_pb, ponio::runge_kutta::rk_33(), yini, { 0., 10. }, dt, fobs );
    }

    return 0;
}
