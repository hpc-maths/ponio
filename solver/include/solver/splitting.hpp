// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <concepts>

#include "splitting/lie.hpp"
#include "splitting/strang.hpp"

namespace ponio::splitting
{
    using lie::make_lie_tuple;
    using strang::make_strang_tuple;

    template <typename T>
    concept is_splitting_method = requires( T t ) { static_cast<bool>( T::is_splitting_method ); };

} // namespace ponio::splitting
