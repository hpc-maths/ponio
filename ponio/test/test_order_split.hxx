// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include <doctest/doctest.h>

#include <ponio/runge_kutta.hpp>
#include <ponio/splitting.hpp>

#include "compute_order.hpp"

TEST_CASE( "detail::power" )
{
    auto lie_splitting    = ponio::splitting::make_lie_tuple( std::make_pair( ponio::runge_kutta::rk_33(), .1 ),
        std::make_pair( ponio::runge_kutta::rk_33(), .1 ),
        std::make_pair( ponio::runge_kutta::rk_33(), .1 ) );
    auto strang_splitting = ponio::splitting::make_strang_tuple( std::make_pair( ponio::runge_kutta::rk_44(), .1 ),
        std::make_pair( ponio::runge_kutta::rk_44(), .1 ),
        std::make_pair( ponio::runge_kutta::rk_44(), .1 ) );

    WARN( check_order( lie_splitting ) == doctest::Approx( lie_splitting.order ).epsilon( 0.125 ) );
    WARN( check_order( strang_splitting ) == doctest::Approx( strang_splitting.order ).epsilon( 0.125 ) );
}
