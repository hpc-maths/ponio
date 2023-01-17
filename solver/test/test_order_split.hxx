// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

#include <doctest/doctest.h>

#include <solver/splitting.hpp>
#include <solver/butcher_methods.hpp>

#include "compute_order.hpp"

TEST_CASE("detail::power")
{
  auto lie_splitting = ode::splitting::make_lie_tuple(
    std::make_pair(ode::butcher::rk_33(), .1),
    std::make_pair(ode::butcher::rk_33(), .1),
    std::make_pair(ode::butcher::rk_33(), .1)
  );
  auto strang_splitting = ode::splitting::make_strang_tuple(
    std::make_pair(ode::butcher::rk_44(), .1),
    std::make_pair(ode::butcher::rk_44(), .1),
    std::make_pair(ode::butcher::rk_44(), .1)
  );

  WARN( check_order(lie_splitting) == doctest::Approx(lie_splitting.order).epsilon(0.125) );
  WARN( check_order(strang_splitting) == doctest::Approx(strang_splitting.order).epsilon(0.125) );
}
