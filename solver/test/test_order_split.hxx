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
    ode::butcher::rk_33<>(),
    ode::butcher::rk_33<>(),
    ode::butcher::rk_33<>()
  );
  auto strang_splitting = ode::splitting::make_strang_tuple(
    ode::butcher::rk_33<>(),
    ode::butcher::rk_33<>(),
    ode::butcher::rk_33<>()
  );

  WARN( check_order(lie_splitting) == doctest::Approx(lie_splitting.order).epsilon(0.05) );
  WARN( check_order(strang_splitting) == doctest::Approx(strang_splitting.order).epsilon(0.05) );
}
