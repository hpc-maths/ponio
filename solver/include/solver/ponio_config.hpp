// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <cstddef> // for std::size_t

namespace ponio::default_config
{
    static constexpr double tol                        = 1e-4;
    static constexpr double newton_tolerance           = 1e-10;
    static constexpr std::size_t newton_max_iterations = 50;
}
