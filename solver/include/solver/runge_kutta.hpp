// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include "runge_kutta/dirk.hpp"
#include "runge_kutta/erk.hpp"
#include "runge_kutta/exprk.hpp"
#include "runge_kutta/lrk.hpp"
#include "runge_kutta/rkc.hpp"

namespace ponio::runge_kutta
{
    // using diagonal_implicit_runge_kutta::make_dirk;

    // using explicit_runge_kutta::explicit_runge_kutta;

    // using exponential_runge_kutta::explicit_exp_rk_butcher;

    // using lawson_runge_kutta::make_lawson;

    using chebyshev::explicit_rkc2; // NOLINT(misc-unused-using-decls): using to improve interface

} // namespace ponio::runge_kutta

#include "runge_kutta/butcher_methods.hpp"
