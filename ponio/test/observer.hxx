// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <ponio/observer.hpp>

TEST_CASE( "observer::file_observer_currentpath" )
{
    auto pwd = std::filesystem::current_path();

    auto test_path = pwd / "test.txt";
    int tn = 0, dt = 1;
    double un = 0.5;

    { // create an observer only in this scope
        observer::file_observer obs( "test.txt" );
        obs( tn, un, dt );
    }

    auto ls = std::filesystem::directory_iterator( pwd );

    auto it = std::find_if( std::filesystem::begin( ls ),
        std::filesystem::end( ls ),
        [&]( auto const& p )
        {
            return p.path() == test_path;
        } );

    // check creation of file
    CHECK( it != std::filesystem::end( ls ) );
    CHECK( *it == test_path );

    // check values in file
    std::ifstream input( test_path );
    int i, j;
    double x;
    input >> i >> x >> j;
    input.close();

    CHECK( i == tn );
    CHECK( x == un );
    CHECK( j == dt );

    std::filesystem::remove( test_path );
}

TEST_CASE( "observer::file_observer_newdir" )
{
    std::filesystem::path test_path = "my_new_unique_dir/test.txt";
    int tn = 0, dt = 1;
    double un = 0.5;

    { // create an observer only in this scope
        observer::file_observer obs( test_path );
        obs( tn, un, dt );
    }

    auto ls = std::filesystem::directory_iterator( test_path.parent_path() );

    auto it = std::find_if( std::filesystem::begin( ls ),
        std::filesystem::end( ls ),
        [&]( auto const& p )
        {
            return p.path() == test_path;
        } );

    // check creation of file
    CHECK( it != std::filesystem::end( ls ) );
    CHECK( *it == test_path );

    // check values in file
    std::ifstream input( test_path );
    int i, j;
    double x;
    input >> i >> x >> j;
    input.close();

    CHECK( i == tn );
    CHECK( x == un );
    CHECK( j == dt );

    std::filesystem::remove_all( test_path.parent_path() );
}
