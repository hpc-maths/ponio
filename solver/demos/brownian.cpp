// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
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
        return 0u;
    }

    static constexpr result_type
    max()
    {
        return RAND_MAX;
    }

    double
    entropy() const noexcept
    {
        return 42.0;
    }

    result_type
    operator()()
    {
        return std::rand();
    }
};

int
main( int argc, char** argv )
{
    std::string dirname = "brownian_data";

    std::size_t n = 10;
    if ( argc > 1 )
    {
        n = std::atoi( argv[1] );
    }

    using state_t = std::valarray<double>;

    // std::random_device rd;
    C_random_device rd;
    std::mt19937 gen( rd() );

    std::normal_distribution<> d{ 0., 2 };

    double dt = 1e-3;

    auto brownian_pb = ode::make_simple_problem(
        [&]( double t, state_t const& y ) -> state_t
        {
            return { d( gen ), d( gen ) };
        } );

    state_t yini = { 0., 0. };

    for ( auto i = 0u; i < n; ++i )
    {
        std::stringstream ssfilename;
        ssfilename << "brownian_" << i << ".dat";
        auto filename = std::filesystem::path( dirname ) / ssfilename.str();
        observer::file_observer fobs( filename );
        ode::solve( brownian_pb, ode::butcher::rk_33(), yini, { 0., 10. }, dt, fobs );
    }

    return 0;
}
