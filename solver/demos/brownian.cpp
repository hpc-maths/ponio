// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <valarray>

#include <solver/butcher_methods.hpp>
#include <solver/observer.hpp>
#include <solver/problem.hpp>
#include <solver/solver.hpp>

struct C_random_device
{
    using result_type = unsigned int;

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
    entropy() noexcept
    {
        return 42.0;
    }

    result_type
    operator()()
    {
        return static_cast<result_type>( std::rand() );
    }
};

int
main( int argc, char* argv[] )
{
    std::string const dirname = "brownian_data";

    std::size_t n = 10;
    if ( argc > 1 )
    {
        n = std::stoul( argv[1] );
    }

    using state_t = std::valarray<double>;

    // std::random_device rd;
    C_random_device rd;
    std::mt19937 gen( rd() );

    std::normal_distribution<> d{ 0., 2 };

    double const dt = 1e-3;

    auto brownian_pb = ode::make_simple_problem(
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
        observer::file_observer fobs( filename );
        ode::solve( brownian_pb, ode::butcher::rk_33(), yini, { 0., 10. }, dt, fobs );
    }

    return 0;
}
