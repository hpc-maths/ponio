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
  WARN( order_lv(ode::splitting::make_lie_tuple(
      ode::butcher::rk_33<>(),
      ode::butcher::rk_33<>(),
      ode::butcher::rk_33<>()
    )) == doctest::Approx(1.0).epsilon(0.05) );
  WARN( order_lv(ode::splitting::make_strang_tuple(
      ode::butcher::rk_33<>(),
      ode::butcher::rk_33<>(),
      ode::butcher::rk_33<>()
    )) == doctest::Approx(2.0).epsilon(0.05) );
}
